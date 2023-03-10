## Helen Kang
## run MAST with batch correction
## 220217




library(conflicted)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("first", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("melt", "reshape2")
conflict_prefer("filter", "dplyr")

packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
              "tidyr") #, "grid", "gtable", "gridExtra","ggrepel",#"ramify",
              ## "ggpubr","gridExtra", "parallel", "future",
              ## "org.Hs.eg.db","limma","conflicted", #"fgsea", 
              ## "cluster","textshape","readxl", 
              ## "ggdist", "gghalves", "Seurat", "writexl", "SingleCellExperiment", "MAST") #              "GGally","RNOmni","usedist","GSEA","clusterProfiler","IsoplotR","wesanderson",
xfun::pkg_attach(packages)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")


##########################################################################################
## Constants and Directories
option.list <- list(
    make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
    make_option("--sampleName", type="character", default="2kG.library", help="sample name"),
    make_option("--barcode.names", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/2kG.library.barcodes.tsv", help="barcodes.tsv for all cells"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--outdirsample", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/", help="path to cNMF analysis results"), ## or for 2n1.99x: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/Perturb_2kG_dup4/K60/threshold_0_2/"
    make_option("--num.genes.per.MAST.runGroup", type="numeric", default=494, help="Number of MAST parallel processes to create"), 
    make_option("--scatteroutput", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/MAST/", help="path to gene breakdown table output"),

    ## script dir
    make_option("--scriptdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts/", help="location for this script and functions script")

)
opt <- parse_args(OptionParser(option_list=option.list))


## ## sdev debug K562 gwps
## opt$sampleName <- "WeissmanK562gwps"
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/data/K562_gwps_raw_singlecell_01_metadata.txt"
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K35/threshold_0_2/"
## opt$scatteroutput <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/MAST/"
## opt$scriptdir <- "/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts"
## opt$K.val <- 35


k <- opt$K.val 
SAMPLE <- opt$sampleName
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
## OUTDIRSAMPLE <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K31/threshold_0_2/"
## OUTDIRSAMPLE <- opt$outdirsample
SCATTEROUTDIR <- opt$scatteroutput

## INPUTDIR <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/inputs/"
## OUTDIR <- OUTDIRSAMPLE
## check.dir <- c(INPUTDIR, OUTDIR, SCATTEROUTDIR)
check.dir <- c(SCATTEROUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

## mytheme <- theme_classic() + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, face = "bold"))
## palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)
## p.adjust.thr <- 0.1



## ##########################################################################################
## ## load data
## ## load known CAD genes
## CAD.genes <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/known_CAD_gene_set.txt", header=F, stringsAsFactors=F) %>% as.matrix %>% as.character

## ## load topic model results
## cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
## print(cNMF.result.file)
## if(file.exists(cNMF.result.file)) {
##     print("loading cNMF result file")
##     load(cNMF.result.file)
## }

## Organize metadata
print("organize meta data")
## get metadata from log2.X.full rownames
# if(SAMPLE %in% c("2kG.library", "Perturb_2kG_dup4")) {
if( grepl("2kG.library|Perturb_2kG_dup4", SAMPLE) ) {
    barcode.names <- read.table(opt$barcode.names, header=F, stringsAsFactors=F) %>% `colnames<-`("long.CBC")
    ## rownames(omega) <- barcode.names %>% pull(long.CBC) %>% gsub("CSNK2B-and-CSNK2B", "CSNK2B",.)
    barcode.names <- barcode.names %>%
        rownames %>%
        as.data.frame %>%
        `colnames<-`("long.CBC") %>%
        separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
        ## separate(col="CBC", into=c("CBC", "sample"), sep="-", remove=F) %>%
        separate(col="CBC", into=c("CBC", "sample"), sep="-scRNAseq_2kG_", remove=F) %>%
        mutate(Gene = gsub("-TSS2", "", Gene.full.name),
               CBC = gsub("RHOA-and-", "", CBC)) %>%
        as.data.frame
    sample.to.10X.lane <- data.frame(sample_num = 1:20,
                                     sample = meta_data$sample %>% unique %>% sort)

    ## meta_data <- merge(meta_data, sample.to.10X.lane, by="sample") %>%
    ##     mutate(CBC_10x = paste0(CBC, "-", sample_num)) %>%
    ##     merge(guidePerCell.df %>% select(CBC_10x, guides_per_cbc, max_umi_ct), by="CBC_10x") %>%
    ##     merge(umiPerCell.df, by="long.CBC") %>%
    ##     merge(geneDetectedPerCell.df, by.x="long.CBC", by.y=0)
} else {
    barcode.names <- read.table(opt$barcode.names, header=T, stringsAsFactors=F) ## %>% `colnames<-`("long.CBC")
}    

cat("finished loading barcode names\n")
## print(paste0("omega dimensions: ", dim(omega)))
cat(paste0("barcode names dimensions: ", dim(barcode.names)[1], " x ", dim(barcode.names)[2], "\n"))
num_genes <- barcode.names %>% pull(Gene) %>% unique %>% length
cat(paste0("total number of perturbations: ", num_genes, "\n"))
numGenesPerRun <- opt$num.genes.per.MAST.runGroup
num_MAST_runs <- floor(num_genes / numGenesPerRun) + 1
cat(paste0("total number of MAST runs: ", num_MAST_runs, "\n"))
gene.ary <- barcode.names %>% pull(Gene) %>% unique %>% sort
## separate into numGenesPerRun perturbations per group for MAST
for (i in 1:(num_MAST_runs-1)) {
    gene.ary.i <- gene.ary[((i-1)*numGenesPerRun+1) : (numGenesPerRun*i)]
    write.table(gene.ary.i, paste0(SCATTEROUTDIR, "/GeneNames_Group", i, ".txt"), sep="\n", quote=F, row.names=F, col.names=F)
}   
i <- i + 1 # the last one
gene.ary.i <- gene.ary[((i-1)*numGenesPerRun+1) : length(gene.ary)]
write.table(gene.ary.i, paste0(SCATTEROUTDIR, "/GeneNames_Group", i, ".txt"), sep="\n", quote=F, row.names=F, col.names=F)
