## Helen Kang
## 230520
## Script to convert cNMF output in yaml files according to IGVF CPN Focus Group requirement
## See: https://github.com/mortazavilab/Topyfic/blob/main/tutorials/topic_modeling_model.md

## Libraries

suppressPackageStartupMessages(library(conflicted))
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("first", "dplyr")
conflict_prefer("rename", "dplyr")

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(reshape2)
    ## library(ggplot2)
    ## library(cowplot)
    ## library(ggpubr) ## ggarrange
    ## library(gplots) ## heatmap.2
    ## library(scales) ## geom_tile gradient rescale
    ## library(ggrepel)
    library(stringr)
    library(stringi)
    library(svglite)
    ## library(Seurat)
    ## library(SeuratObject)
    library(xlsx)
    library(yaml)
})


##########################################################################################
## Constants and Directories
option.list <- list(
    make_option("--sampleName", type="character", default="WeissmanK562gwps", help="Name of the sample"),
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/", help="Output directory"),
    make_option("--K.val", type="numeric", default=90, help="K value to analyze"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    ## make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    ## make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--perturbSeq", type="logical", default=TRUE, help="Whether this is a Perturb-seq experiment")
)
opt <- parse_args(OptionParser(option_list=option.list))



SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
OUTDIR=opt$outdir
## TMDIR=opt$topic.model.result.dir
## SEP=opt$sep
k <- opt$K.val
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
## FIGDIR=opt$figdir
## FIGDIRSAMPLE=paste0(FIGDIR, "/", SAMPLE, "/K",k,"/")
## FIGDIRTOP=paste0(FIGDIRSAMPLE,"/",SAMPLE,"_K",k,"_dt_", DENSITY.THRESHOLD,"_")
OUTDIRSAMPLE=paste0(OUTDIR, "/", SAMPLE, "/K",k,"/threshold_", DENSITY.THRESHOLD, "/")
OUTDIRSAMPLEIGVF = paste0(OUTDIRSAMPLE, "IGVF_format/")

check.dir <- c(OUTDIRSAMPLEIGVF)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))


## subscript for files
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
## SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## cNMF direct output file (GEP)
cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
## cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.k_", k, ".dt_", density.threshold, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
} else {
    print(paste0("file ", cNMF.result.file, " not found"))
}

db <- ifelse(grepl("mouse", SAMPLE), "org.Mm.eg.db", "org.Hs.eg.db")
library(!!db) ## load the appropriate database
topic.gene.names <- rownames(theta.zscore)
topic.gene.name.type <- ifelse(grepl("^ENSG", topic.gene.names) %>% as.numeric %>% sum == length(topic.gene.names), "ENSGID", "SYMBOL")
if(topic.gene.name.type == "ENSGID") {
    ENSGID.gene.names <- topic.gene.names
    SYMBOL.gene.names <- mapIds(get(db), keys=topic.gene.names, keytype = "ENSEMBL", column = "SYMBOL")
    SYMBOL.gene.names[is.na(SYMBOL.gene.names)] <- ENSGID.gene.names[is.na(SYMBOL.gene.names)]
} else {
    SYMBOL.gene.names <- topic.gene.names
    ENSGID.gene.names <- mapIds(get(db), keys=topic.gene.names, keytype = "SYMBOL", column = "ENSEMBL")
    ENSGID.gene.names[is.na(ENSGID.gene.names)] <- SYMBOL.gene.names[is.na(ENSGID.gene.names)]
}

## 1. Model YAML file: Capture all the information about the dataset that you used and which method and how you run it along with Topic_ID which points to Topic YAML files and Cell-Topic participation ID for pointing out the Cell-Topic participation h5ad file.
out <- list("Assay" = NULL,
            "Cell-Topic participation ID" = NULL,
            "Experiment ID" = SAMPLE,
            "Name of method" = "cNMF",
            "Number of topics" = k,
            "Technology" = "10x",
            "cNMF spectra threshold" = opt$density.thr,
            "Topic IDs" = paste0(SAMPLE, "_K", k, "_", 1:k),
            "level" = "cell line",
            "cell type" = "K562")
write_yaml(out, paste0(OUTDIRSAMPLEIGVF, SAMPLE, ".", SUBSCRIPT.SHORT, ".modelYAML.yaml"))

## 2. Topics YAML files: Capture all the information about Topics including Topic_ID, gene_weight, gene_id and gene_name and any other information that suits your data
theta.zscore.long <- theta.zscore %>%
    as.data.frame %>%
    `colnames<-`(paste0(SAMPLE, "_K", k, "_", colnames(.))) %>%
    mutate(Gene = SYMBOL.gene.names,
           ENSGID = ENSGID.gene.names) %>%
    melt(id.vars=c("Gene", "ENSGID"), value.name="Gene weights", variable.name="Topic ID") %>%
    rename("gene_id" = "ENSGID")

for( t in 1:k ) {
    theta.zscore.long.here <- theta.zscore.long %>%
        subset(grepl(paste0("_", t, "$"), `Topic ID`))
    duplicated.index <- duplicated(theta.zscore.long.here$Gene)
    theta.zscore.long.here$Gene[duplicated.index] <- paste0(theta.zscore.long.here$Gene[duplicated.index], "_", theta.zscore.long.here$ENSGID[duplicated.index])

    ## create output list
    out <- list("gene_id" = theta.zscore.long.here %>%
                    `rownames<-`(.$Gene) %>%
                    select(gene_id) %>% t %>% as.data.frame,
                "Gene weights" = theta.zscore.long.here %>%
                    `rownames<-`(.$Gene) %>%
                    select(`Gene weights`) %>% t %>% as.data.frame,
                "Topic ID" = paste0(SAMPLE, "_K", k, "_", t))
    
    write_yaml(out, paste0(OUTDIRSAMPLEIGVF, SAMPLE, "_", SUBSCRIPT.SHORT, "_program", t, "_topicYAML.yaml"))
}
