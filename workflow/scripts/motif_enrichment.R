## Helen Kang
## Motif Enrichment
## 220310


library(conflicted)
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")
packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
              "tidyr", "grid", "gtable", "gridExtra","ggrepel","ramify",
              "ggpubr","gridExtra",
              "org.Hs.eg.db","limma","fgsea",
              "cluster","textshape","readxl", 
              "ggdist", "gghalves", "writexl") #              "GGally","RNOmni","usedist","GSEA","clusterProfiler","IsoplotR","wesanderson",
## packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
##               "CountClust", "Hmisc", "tidyr", "grid", "gtable", "gridExtra","ggrepel","ramify",
##               "GGally","RNOmni","usedist","ggpubr","gridExtra","GSEA",
##               "org.Hs.eg.db","limma","clusterProfiler","fgsea", "conflicted",
##               "cluster","textshape","readxl", "IsoplotR", "wesanderson", 
##               "ggdist", "gghalves", "Seurat", "writexl")
                                        # library(Seurat)
xfun::pkg_attach(packages)
conflict_prefer("slice", "dplyr")


## source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")
source("workflow/scripts/motif_enrichment_functions.R")

option.list <- list(
    make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220307_prioritized_topic_motif_enrichment/figures/", help="Figure directory"),
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220307_prioritized_topic_motif_enrichment/outputs/", help="Output directory"),
                                        # make_option("--olddatadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/", help="Input 10x data directory"),
    make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/", help="Input 10x data directory"),
                                        # make_option("--topic.model.result.dir", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/210707_snakemake_maxParallel/all_genes_acrossK/2kG.library/", help="Topic model results directory"),
    make_option("--sampleName", type="character", default="FT010_fresh_4min", help="Name of Samples to be processed, separated by commas"),
                                        # make_option("--sep", type="logical", default=F, help="Whether to separate replicates or samples"),
    make_option("--K.list", type="character", default="2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25", help="K values available for analysis"),
    make_option("--K.val", type="numeric", default=20, help="K value to analyze"),
    make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--ABCdir",type="character", default="/oak/stanford/groups/engreitz/Projects/ABC/200220_CAD/ABC_out/TeloHAEC_Ctrl/Neighborhoods/", help="Path to ABC enhancer directory"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
                                        # make_option("--raw.mtx.dir",type="character",default="stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/data/no_IL1B_filtered.normalized.ptb.by.gene.mtx.filtered.txt", help="input matrix to cNMF pipeline"),
                                        # make_option("--raw.mtx.RDS.dir",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.RDS", help="input matrix to cNMF pipeline"), # the first lane: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.expandedMultiTargetGuide.RDS"
                                        # make_option("--subsample.type", type="character", default="", help="Type of cells to keep. Currently only support ctrl"),
                                        # make_option("--barcode.names", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/barcodes.tsv", help="barcodes.tsv for all cells"),
    make_option("--reference.table", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/210702_2kglib_adding_more_brief_ca0713.xlsx"),
    
    ## fisher motif enrichment
    ## make_option("--outputTable", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.k_14.df_0_2.txt", help="Output directory"),
    ## make_option("--outputTableBinary", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210607_snakemake_output/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.binary.k_14.df_0_2.txt", help="Output directory"),
    ## make_option("--outputEnrichment", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210607_snakemake_output/outputs/no_IL1B/topic.top.100.zscore.gene.motif.fisher.enrichment.k_14.df_0_2.txt", help="Output directory"),
    make_option("--motif.promoter.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv", help="All promoter's motif matches"),
    make_option("--motif.enhancer.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv", help="All enhancer's motif matches specific to {no,plus}_IL1B"),
    make_option("--enhancer.fimo.threshold", type="character", default="1.0E-4", help="Enhancer fimo motif match threshold"),
    make_option("--ep.type", type="character", default="enhancer", help="motif enrichment for enhancer or promoter, specify 'enhancer' or 'promoter'"),

                                        #summary plot parameters
    make_option("--test.type", type="character", default="per.guide.wilcoxon", help="Significance test to threshold perturbation results"),
    make_option("--adj.p.value.thr", type="numeric", default=0.1, help="adjusted p-value threshold"),
    make_option("--recompute", type="logical", default=F, help="T for recomputing statistical tests and F for not recompute"),
    make_option("--motif.match.thr.str", type="character", default="pval0.0001", help="threshold for subsetting motif matches"),

    ## Organism flag
    make_option("--organism", type="character", default="human", help="Organism type, accept org.Hs.eg.db. Only support human and mouse.")

    
)
opt <- parse_args(OptionParser(option_list=option.list))

## ## ## all genes directories (for sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/"
## opt$K.val <- 60
## opt$sampleName <- "2kG.library"


# ## ## all genes directories (for sdev)
# opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/"
# opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/"
# opt$K.val <- 60
# opt$sampleName <- "2kG.library"



## ## debug ctrl
## opt$topic.model.result.dir <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/210810_snakemake_ctrls/all_genes_acrossK/2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation"
## opt$sampleName <- "2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210810_snakemake_ctrls/figures/2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210810_snakemake_ctrls/analysis/2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation/all_genes/"
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210806_curate_ctrl_mtx/outputs/2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation.barcodes.tsv"
## opt$K.val <- 60

## ## K562 gwps sdev
## opt$sampleName <- "WeissmanK562gwps"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/"
## opt$K.val <- 100
## opt$ep.type <- "enhancer"
## opt$motif.match.thr.str <- "qval0.1"
## opt$motif.enhancer.background <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/fimo/fimo_out/fimo.formatted.tsv"
## opt$motif.promoter.background <- "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv"

## ## ENCODE mouse heart
## opt$sampleName <- "mouse_ENCODE_heart"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/figures/top2000VariableGenes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/analysis/top2000VariableGenes"
## opt$K.val <- 5
## opt$ep.type <- "enhancer"
## opt$motif.match.thr.str <- "pval1e-4"
## opt$motif.enhancer.background <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/analysis/top2000VariableGenes/mouse_ENCODE_heart/fimo/fimo_out/fimo.txt"
## opt$motif.promoter.background <- "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv"

## ## no_IL1B 200 gene library
## opt$sampleName <- "no_IL1B"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/figures/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/analysis/all_genes"
## opt$K.val <- 14
## opt$ep.type <- "promoter"
## opt$motif.match.thr.str <- "qval0.1"
## opt$motif.enhancer.background <- "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.tsv"
## opt$motif.promoter.background <- "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv"

## ## IGVF b01_LeftCortex sdev
## opt$sampleName <- "IGVF_b01_LeftCortex"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230706_snakemake_igvf_b01_LeftCortex/figures/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230706_snakemake_igvf_b01_LeftCortex/analysis/all_genes"
## opt$K.val <- 20
## opt$ep.type <- "enhancer"
## opt$organism <- "mouse"
## opt$motif.match.thr.str <- "pval1e-6"
## opt$motif.enhancer.background <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230706_snakemake_igvf_b01_LeftCortex/analysis/all_genes/IGVF_b01_LeftCortex/fimo/fimo_out/fimo.txt"
## opt$motif.promoter.background <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230706_snakemake_igvf_b01_LeftCortex/analysis/all_genes/IGVF_b01_LeftCortex/fimo/fimo_out/fimo.txt"

## ## PD DA neuron
## opt$sampleName <- "dan"
## opt$figdir <- "/broad/macosko/ferris/NMF/240606_snakemake_nurr_da/figures/top2000VariableGenes/"
## opt$outdir <- "/broad/macosko/ferris/NMF/240606_snakemake_nurr_da/analysis/top2000VariableGenes/"
## opt$K.val <- 5
## opt$density.thr <- 0.2
## opt$ep.type <- "promoter"
## opt$organism <- "human"
## opt$motif.match.thr.str <- "pval1e-6"
## opt$motif.enhancer.background <- "/broad/macosko/ferris/NMF/240606_snakemake_nurr_da/analysis/top2000VariableGenes/dan/fimo/fimo_out/fimo.tsv"
## opt$motif.promoter.background <- "/broad/macosko/ferris/NMF/240606_snakemake_nurr_da/analysis/top2000VariableGenes/dan/fimo/fimo_out/fimo.tsv"



mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))

SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
DATADIR=opt$olddatadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
                                        # TMDIR=opt$topic.model.result.dir
                                        # SEP=opt$sep
                                        # K.list <- strsplit(opt$K.list,",") %>% unlist() %>% as.numeric()
k <- opt$K.val
num.top.genes <- 300 ## number of topic defining genes
ep.type <- opt$ep.type
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
FIGDIR=opt$figdir
FIGDIRSAMPLE=paste0(FIGDIR, "/", SAMPLE, "/K",k,"/")
FIGDIRTOP=paste0(FIGDIRSAMPLE,"/",SAMPLE,"_K",k,"_dt_", DENSITY.THRESHOLD,"_")
OUTDIRSAMPLE=paste0(OUTDIR, "/", SAMPLE, "/K",k,"/threshold_", DENSITY.THRESHOLD, "/")
FGSEADIR=paste0(OUTDIRSAMPLE,"/fgsea/")
FGSEAFIG=paste0(FIGDIRSAMPLE,"/fgsea/")

## subscript for files
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
## SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr
motif.match.thr.str <- opt$motif.match.thr.str


## create dir if not already
check.dir <- c(OUTDIR, FIGDIR, paste0(FIGDIR,SAMPLE,"/"), paste0(FIGDIR,SAMPLE,"/K",k,"/"), paste0(OUTDIR,SAMPLE,"/"), OUTDIRSAMPLE, FIGDIRSAMPLE, FGSEADIR, FGSEAFIG)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

## palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)


######################################################################
## Load topic model results

cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
}


## not used
# ## load hg38 promoter region file
# promoter.region.hg38.original <- read.table(file=paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS300bp.hg38.bed"), header = F, stringsAsFactors = F)  %>%
#     `colnames<-`(c("chr","start","end", "gene","cell.type","strand")) %>%
#     mutate(sequence_name = paste0(chr, ":", start, "-", end, "(", strand, ")"))


# ## keep only the expressed gene's motifs for background
# expressed.genes <- theta %>% rownames()
# promoter.region.hg38 <- promoter.region.hg38.original %>%
#     mutate(expressed = (gene %in% expressed.genes) ) %>%
#     filter(expressed)
###########


## load FIMO matched motifs ## ifelse on promoter vs enhancer
if (ep.type == "promoter") {
    ## load promoter motif matches
    motif.background <- read.delim(file=paste0(ifelse(opt$motif.promoter.background!="", opt$motif.promoter.background, "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv")), header=F, stringsAsFactors=F) ## %>% filter(!grepl("#", motif_id))  # 30 seconds
    if(ncol(motif.background) > 9) {
        motif.background <- motif.background %>%
            `colnames<-`(c("motif_id", "motif_alt_id", "enhancer_region", "enhancer_type", "gene_region","sequence_name","start","stop","motif.matched.strand","score","p.value","q.value","matched_sequence")) %>% filter(!grepl("#|motif_id", motif_id))  # more than 30 seconds, minutes?
        motif.background <- motif.background %>% filter(grepl("promoter", enhancer_type))
        
    } else {
        motif.background <- motif.background %>%
            `colnames<-`(c("motif_id", "sequence_name", "start", "stop", "motif.matched.strand", "score", "p.value", "q.value", "matched_sequence")) 
    }
    ## old
    ## motif.background <- read.delim(file=paste0(ifelse(opt$motif.promoter.background!="", opt$motif.promoter.background, "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv")), header=T, stringsAsFactors=F) ## %>% filter(!grepl("#", motif_id))  # 30 seconds
    ## ## colnames(motif.background) <- c("motif_id", "motif_alt_id", "enhancer_region", "enhancer_type", "gene_region","sequence_name","start","stop","motif.matched.strand","score","p.value","q.value","matched_sequence")
    ## colnames(motif.background)[colnames(motif.background) == "strand"] <- "motif.matched.strand"
    ## motif.background <- motif.background %>%
    ##     mutate(motif.short = strsplit(motif_id, split="_") %>% sapply("[[", 1) %>% as.character)
    ## end of old

} else {
    ## load enhancer motif matches
    print(opt$motif.enhancer.background)
    motif.background <- read.delim(file=paste0(ifelse(opt$motif.enhancer.background!="", opt$motif.enhancer.background, "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv")), header=F, stringsAsFactors=F)
        if(ncol(motif.background) > 9) {
        motif.background <- motif.background %>%
            `colnames<-`(c("motif_id", "motif_alt_id", "enhancer_region", "enhancer_type", "gene_region","sequence_name","start","stop","motif.matched.strand","score","p.value","q.value","matched_sequence")) %>% filter(!grepl("#|motif_id", motif_id))  # more than 30 seconds, minutes?
        motif.background <- motif.background %>% filter(!grepl("promoter", enhancer_type))
        
    } else {
        motif.background <- motif.background %>%
            `colnames<-`(c("motif_id", "sequence_name", "start", "stop", "motif.matched.strand", "score", "p.value", "q.value", "matched_sequence")) 
    }
}
motif.background <- motif.background %>% mutate(motif.short = strsplit(motif_id, split="_") %>% sapply("[[", 1) %>% as.character)
message("finished loading motif input")


## subset to q.value < 0.1 
if (grepl("qval", motif.match.thr.str)) {
    motif.background <- motif.background %>%
        subset(q.value < 0.1)
} else {
    print(paste0("subset to ", motif.match.thr.str))
    threshold <- gsub("pval", "", motif.match.thr.str) %>% as.numeric
    motif.background <- motif.background %>%
        subset(p.value < threshold)
}


if(ep.type == "enhancer") {
    if(ncol(motif.background) == 9 | sum(as.numeric(grepl("|", motif.background$sequence_name))) == nrow(motif.background)) {
        motif.background <- motif.background %>%
            filter(!grepl("promoter", sequence_name) & !grepl("start", start)) %>%
            separate(col="sequence_name", into=c("enhancer_region", "enhancer_type", "gene_region", "gene_name_sequence_region"), sep="[|]") %>%
            separate(col="gene_name_sequence_region", into=c("sequence_name", "to_remove"), sep="::") %>%
            select(-to_remove)
    }
    colnames(motif.background)[colnames(motif.background) == "strand"] <- "motif.matched.strand"
    motif.background <- motif.background %>%
        mutate(motif.short = strsplit(motif_id, split="_") %>% sapply("[[", 1) %>% as.character)
}


expressed.genes <- rownames(theta.zscore)
## todo: convert expressed genes to symbol if they are not in symbol
db <- ifelse(grepl("mouse|org.Mm.eg.db", opt$organism), "org.Mm.eg.db", "org.Hs.eg.db")
gene.type <- ifelse(length(expressed.genes) == sum(as.numeric(grepl("^ENS", expressed.genes))), "ENSGID", "Gene")
if(gene.type == "ENSGID") expressed.genes = mapIds(get(db), keys=expressed.genes, keytype="ENSEMBL", column="SYMBOL")
motif.background <- motif.background %>%
    subset(sequence_name %in% expressed.genes)





####################################################################################################
## get list of topic defining genes
theta.rank.list <- vector("list", ncol(theta.zscore))## initialize storage list
for(i in 1:ncol(theta.zscore)) {
    topic <- paste0("topic_", colnames(theta.zscore)[i])
    theta.rank.list[[i]] <- theta.zscore[,i] %>%
        as.data.frame %>%
        `colnames<-`("topic.zscore") %>%
        mutate(Gene = rownames(.)) %>%
        arrange(desc(topic.zscore), .before="topic.zscore") %>%
        mutate(zscore.specificity.rank = 1:n()) %>% ## add rank column
        mutate(Topic = topic) ## add topic column
}
theta.rank.df <- do.call(rbind, theta.rank.list) ## combine list to df
topic.defining.gene.df <- theta.rank.df %>%
    subset(zscore.specificity.rank <= num.top.genes) ## select top 300 genes for each topic
gene.type <- ifelse(nrow(topic.defining.gene.df) == sum(as.numeric(grepl("^ENS", topic.defining.gene.df$Gene))), "ENSGID", "Gene")
if(gene.type=="ENSGID") topic.defining.gene.df <- topic.defining.gene.df %>% mutate(ENSGID = Gene) %>% mutate(Gene = mapIds(get(db), keys=.$ENSGID, keytype="ENSEMBL", column="SYMBOL"))

## topic.defining.gene.df <- topic.defining.gene.df %>% mutate(Gene = toupper(Gene))

topic.motif.match.df <- merge(topic.defining.gene.df, motif.background %>% 
                                                      select(motif_id, motif.short, sequence_name, score, p.value, q.value, motif.matched.strand), by.x="Gene", by.y="sequence_name", all.y=T) ## filtered motif.background to genes expressed in this data set, so keep all


motif.id.type <- "motif.short" ## or "motif_id"
topic.motif.match.df.long <- topic.motif.match.df %>%
    group_by(Gene, get(motif.id.type), Topic) %>%
    summarize(count = n()) %>%
    `colnames<-`(c("gene", motif.id.type, "Topic", "count")) %>%
    as.data.frame
 
wide <- topic.motif.match.df.long %>%
    spread(key=motif.id.type, value="count", fill=0) %>% ## get each motif's count in each top promoter
    mutate(topic.value=1) %>%
    spread(key=Topic, value=topic.value, fill=0) ## %>% select(-topic_NA) ## get presence/absense of promoter in topic ## 5 seconds


print(paste0(ep.type, " motif (", motif.match.thr.str, ") count t-test"))
ttest.df <- ttest.on.motifs(wide)


## significant.motifs.df <- ttest.df %>%
##     group_by(topic) %>%
##     subset(p.adjust < 0.1 &
##            enrichment.log2fc > 0) %>%
##     arrange(p.adjust) %>%
##     summarize(motifs = paste0(motif %>% sort, collapse=",")) %>%
##     as.data.frame


####################################################################################################
## save results
all.ttest.df.path <- paste0(OUTDIRSAMPLE,"/", ep.type, ".topic.top.", num.top.genes, ".zscore.gene_motif.count.ttest.enrichment_motif.thr.", motif.match.thr.str, "_", SUBSCRIPT.SHORT,".txt")
print(paste0("saving to ", all.ttest.df.path))
write.table(ttest.df, all.ttest.df.path, sep="\t", quote=F, row.names=F)
# write.table(significant.motifs.df, "./outputs/significant.motifs.df.csv", sep=",", quote=F, row.names=F)
# write.table(significant.motifs.df, "./outputs/significant.motifs.df.txt", sep="\t", quote=F, row.names=F)



## venn diagram of this result



