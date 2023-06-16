## Helen Kang
## Batch and topic correlation (Perturb-seq only)
## 211216

library(conflicted)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")


packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
              "tidyr", "grid", "gtable", "gridExtra","ggrepel","ramify",
              "ggpubr","gridExtra",
              "org.Hs.eg.db","limma","fgsea", "conflicted",
              "cluster","textshape","readxl", 
              "ggdist", "gghalves", "Seurat", "writexl", "purrr") #              "GGally","RNOmni","usedist","GSEA","clusterProfiler","IsoplotR","wesanderson",
xfun::pkg_attach(packages)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")



## source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")

option.list <- list(
  make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/figures/all_genes/", help="Figure directory"),
  make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/all_genes/", help="Output directory"),
  # make_option("--olddatadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/", help="Input 10x data directory"),
  make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/", help="Input 10x data directory"),
  make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
  make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
  make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
  make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
  make_option("--ABCdir",type="character", default="/oak/stanford/groups/engreitz/Projects/ABC/200220_CAD/ABC_out/TeloHAEC_Ctrl/Neighborhoods/", help="Path to ABC enhancer directory"),
  make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
  make_option("--reference.table", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/210702_2kglib_adding_more_brief_ca0713.xlsx"),
  make_option("--barcode.names", type="character", default="", help="metadata CBC and sample information data table"),
  
  
  #summary plot parameters
  make_option("--test.type", type="character", default="per.guide.wilcoxon", help="Significance test to threshold perturbation results"),
  make_option("--adj.p.value.thr", type="numeric", default=0.1, help="adjusted p-value threshold"),
  make_option("--recompute", type="logical", default=F, help="T for recomputing statistical tests and F for not recompute")
  
)
opt <- parse_args(OptionParser(option_list=option.list))

## ## all genes directories (for sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/"
## opt$K.val <- 60

# ## control genes directories (for sdev)
# opt$sampleName <- "2kG.library.ctrl.only"
# opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/figures/all_genes/"
# opt$topic.model.result.dir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/all_genes_acrossK/2kG.library.ctrl.only/"
# opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/analysis/all_genes/"
# opt$K.val <- 60

## ## overdispersed gene directories (for sdev)
## opt$sampleName <- "2kG.library_overdispersedGenes"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220716_snakemake_overdispersedGenes/figures/top2000VariableGenes"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220716_snakemake_overdispersedGenes/analysis/top2000VariableGenes"
## opt$K.val <- 120

## ## K562 gwps sdev
## opt$sampleName <- "WeissmanK562gwps"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/"
## opt$K.val <- 90
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/data/K562_gwps_raw_singlecell_01_metadata.txt"

## ## ENCODE Mouse Heart data
## opt$sampleName <- "mouse_ENCODE_heart"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/analysis/top2000VariableGenes"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/figures/top2000VariableGenes"
## opt$K.val <- 55
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/collab_data/IGVF/mouse_ENCODE_heart/auxiliary_data/snrna/heart_Parse_10x_integrated_metadata.csv"

## ## teloHAEC no_IL1B 200 gene library
## opt$sampleName <- "no_IL1B"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/analysis/all_genes/"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/figures/all_genes/"
## opt$K.val <- 20
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/data/no_IL1B.barcodes.txt"

 
SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
## DATADIR=opt$olddatadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
k <- opt$K.val
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
FIGDIR=opt$figdir
FIGDIRSAMPLE=paste0(FIGDIR, "/", SAMPLE, "/K",k,"/")
FIGDIRTOP=paste0(FIGDIRSAMPLE,"/",SAMPLE,"_K",k,"_dt_", DENSITY.THRESHOLD,"_")
OUTDIRSAMPLE=paste0(OUTDIR, "/", SAMPLE, "/K",k,"/threshold_", DENSITY.THRESHOLD, "/")


## subscript for files
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr

# create dir if not already
check.dir <- c(OUTDIR, FIGDIR, OUTDIRSAMPLE, FIGDIRSAMPLE)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))


## graphing constants and helpers
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)


##################################################
## load data
cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
if(file.exists(cNMF.result.file)) {
    message(paste0("loading cNMF result file: \n", cNMF.result.file))
    load(cNMF.result.file)
} else {
	print(paste0(cNMF.result.file, " does not exist"))
}


## annotate omega to get Gene, Guide, Sample, and CBC
if(grepl("2kG.library", SAMPLE)) {
    ann.omega <- cbind(omega, barcode.names)  ## %>%
} else {
    if(grepl("[.]csv", opt$barcode.names)) barcode.names <- read.delim(opt$barcode.names, stringsAsFactors=F, sep=",") else barcode.names <- read.delim(opt$barcode.names, stringsAsFactors=F)
    ann.omega <- merge(omega, barcode.names %>% select(CBC, sample), by.x=0, by.y="CBC", all.x=T)
}


## Batch Effect QC: correlate batch binary labels with topic expression
## Batch binary label
if(grepl("2kG.library", SAMPLE)) {
    ann.omega.batch <- ann.omega %>% mutate(sample.short = gsub("scRNAseq_2kG_", "", sample) %>% gsub("_.*$","", .))
    ann.omega.batch.binary <- ann.omega.batch %>% mutate(tmp.value = 1) %>% spread(key="sample.short", fill=0, value="tmp.value")
    ann.omega.batch.binary.mtx <- ann.omega.batch.binary %>% select(-long.CBC,-Gene.full.name,-Guide,-CBC,-sample,-Gene) %>% as.matrix()
    m <- cor(ann.omega.batch.binary.mtx, method="pearson") %>% as.matrix()
    batch.correlation.mtx <- m[1:k,(k+1):(dim(m)[2])]
} 
ann.omega.sample.batch.binary <- ann.omega %>% mutate(tmp.value = 1) %>% spread(key="sample", fill=0, value="tmp.value")
ann.omega.sample.batch.binary.mtx <- ann.omega.sample.batch.binary %>% select(-Row.names) %>% as.matrix()
m <- cor(ann.omega.sample.batch.binary.mtx, method="pearson") %>% as.matrix()
sample.batch.correlation.mtx <- m[1:k, (k+1):(dim(m)[2])]


## calculate percent of topics with correlation past a threshold (0.1, 0.2, 0.4, 0.6)
correlation.threshold.list <- c(0.1, 0.2, 0.4, 0.6)
batch.passed.threshold.df <- do.call(rbind, lapply(correlation.threshold.list, function(threshold) {
    df <- sample.batch.correlation.mtx %>% apply(1, function(x) (x > threshold) %>% as.numeric %>% sum) %>% as.data.frame %>% `colnames<-`("num.batch.correlated") %>% mutate(batch.thr = threshold) %>% mutate(ProgramID = rownames(.), K = k)
}))
batch.percent.df <- batch.passed.threshold.df %>% group_by(batch.thr) %>% summarize(percent.correlated = ((num.batch.correlated > 0) %>% as.numeric %>% sum) / k) %>% mutate(K = k)

## max batch correlation per topic
max.batch.correlation.df <- sample.batch.correlation.mtx %>%
    apply(1, function(x) {
        out <- max(abs(x))
    }) %>%
    as.data.frame %>%
    `colnames<-`("maxPearsonCorrelation") %>%
    mutate(ProgramID = row.names(.)) %>%
    as.data.frame


## store batch and sample correlation matrix
if(grepl("2kG.library", SAMPLE)) write.table(batch.correlation.mtx, file=paste0(OUTDIRSAMPLE, "/batch.correction.mtx.txt"), sep="\t", quote=F)
write.table(sample.batch.correlation.mtx, file=paste0(OUTDIRSAMPLE, "/sample.batch.correction.mtx.txt"), sep="\t", quote=F)
write.table(batch.passed.threshold.df, file=paste0(OUTDIRSAMPLE, "/batch.passed.thr.df.txt"), sep="\t", quote=F, row.names=F)
write.table(batch.percent.df, file=paste0(OUTDIRSAMPLE, "/batch.percent.df.txt"), sep="\t", quote=F, row.names=F)
write.table(max.batch.correlation.df, file=paste0(OUTDIRSAMPLE, "/max.batch.correlation.df.txt"), sep="\t", quote=F, row.names=F)
if(grepl("2kG.library", SAMPLE)) {
    save(batch.correlation.mtx, sample.batch.correlation.mtx, batch.passed.threshold.df, batch.percent.df, max.batch.correlation.df,
     file=paste0(OUTDIRSAMPLE, "/batch.correlation.RDS"))
} else {
    save(sample.batch.correlation.mtx, batch.passed.threshold.df, batch.percent.df, max.batch.correlation.df,
         file=paste0(OUTDIRSAMPLE, "/batch.correlation.RDS"))
}


## Batch correlation heatmap
plotHeatmap <- function(mtx, title){
    heatmap.2(
        mtx, 
        Rowv=T, 
        Colv=T,
        trace='none',
        key=T,
        col=palette,
        labCol=colnames(mtx),
        ## margins=c(15,5), 
        cex.main=0.1, 
        cexCol=1/(nrow(mtx)^(1/7)), cexRow=1/(ncol(mtx)^(1/7)),
        main=title
    )
}


pdf(paste0(FIGDIRTOP, "batch.correlation.heatmap.pdf"),width=0.15*ncol(sample.batch.correlation.mtx)+5, height=0.1*nrow(sample.batch.correlation.mtx)+5)
if(grepl("2kG.library", SAMPLE)) plotHeatmap(batch.correlation.mtx, title=paste0(SAMPLE, ", K=", k, ", topic batch correlation"))
plotHeatmap(sample.batch.correlation.mtx, title=paste0(SAMPLE, ", K=", k, ", topic sample correlation"))
dev.off()

## automate batch topic selection
Pearson.correlation.threshold <- 0.1 ## is this a good threshold for all values of K? ## check CDF of average correlation?
batch.topic <- apply(sample.batch.correlation.mtx, 1, 
                    function(x) sum(as.numeric(abs(x) > Pearson.correlation.threshold))) %>% 
                keep(function(x) x > 0) %>% 
                names
write.table(batch.topic, file=paste0(OUTDIRSAMPLE, "batch.topics.txt"), quote=F, sep="\t", row.names=F, col.names=F)
