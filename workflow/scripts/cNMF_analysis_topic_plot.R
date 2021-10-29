## Helen Kang
## Topic Model Analysis 
## 210503

## .libPaths("/home/groups/engreitz/Software/R_3.6.1")


# packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
#               "CountClust", "Hmisc", "tidyr", "grid", "gtable", "gridExtra","ggrepel","ramify",
#               "GGally","RNOmni","usedist","ggpubr","gridExtra","GSEA",
#               "org.Hs.eg.db","limma","clusterProfiler","fgsea", "conflicted",
#               "cluster","textshape","readxl", "IsoplotR", "wesanderson", 
#               "ggdist", "patchwork", "gghalves", "Seurat", "writexl")
# library(textshape)
# library(readxl)
packages <- c("optparse","dplyr", "ggplot2", "reshape2", "ggrepel", "conflicted")
## library(Seurat)
xfun::pkg_attach(packages)
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")


## source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")

option.list <- list(
  make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/figures/all_genes/", help="Figure directory"),
  make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/", help="Output directory"),
  # make_option("--olddatadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/", help="Input 10x data directory"),
  make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/", help="Input 10x data directory"),
  # make_option("--topic.model.result.dir", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/210625_snakemake_output/top3000VariableGenes_acrossK/2kG.library/", help="Topic model results directory"),
  make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
  # make_option("--sep", type="logical", default=F, help="Whether to separate replicates or samples"),
  # make_option("--K.list", type="character", default="2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25", help="K values available for analysis"),
  make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
  # make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
  # make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
  # make_option("--ABCdir",type="character", default="/oak/stanford/groups/engreitz/Projects/ABC/200220_CAD/ABC_out/TeloHAEC_Ctrl/Neighborhoods/", help="Path to ABC enhancer directory"),
  make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
  # make_option("--raw.mtx.dir",type="character",default="stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/data/no_IL1B_filtered.normalized.ptb.by.gene.mtx.filtered.txt", help="input matrix to cNMF pipeline"),
  # make_option("--raw.mtx.RDS.dir",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.RDS", help="input matrix to cNMF pipeline"), # the first lane: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.expandedMultiTargetGuide.RDS"
  # make_option("--subsample.type", type="character", default="", help="Type of cells to keep. Currently only support ctrl"),
  # make_option("--barcode.names", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/barcodes.tsv", help="barcodes.tsv for all cells"),
  # make_option("--reference.table", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/210702_2kglib_adding_more_brief_ca0713.xlsx"),
  
  ## fisher motif enrichment
  ## make_option("--outputTable", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.k_14.df_0_2.txt", help="Output directory"),
  ## make_option("--outputTableBinary", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210607_snakemake_output/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.binary.k_14.df_0_2.txt", help="Output directory"),
  ## make_option("--outputEnrichment", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210607_snakemake_output/outputs/no_IL1B/topic.top.100.zscore.gene.motif.fisher.enrichment.k_14.df_0_2.txt", help="Output directory"),
  # make_option("--motif.promoter.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv", help="All promoter's motif matches"),
  # make_option("--motif.enhancer.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv", help="All enhancer's motif matches specific to {no,plus}_IL1B"),
  # make_option("--enhancer.fimo.threshold", type="character", default="1.0E-4", help="Enhancer fimo motif match threshold"),

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

mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))

SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
# STATIC.SAMPLE=c("Telo_no_IL1B_T200_1", "Telo_no_IL1B_T200_2", "Telo_plus_IL1B_T200_1", "Telo_plus_IL1B_T200_2", "no_IL1B", "plus_IL1B",  "pooled")
# DATADIR=opt$olddatadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
## TMDIR=opt$topic.model.result.dir
## SEP=opt$sep
# K.list <- strsplit(opt$K.list,",") %>% unlist() %>% as.numeric()
k <- opt$K.val
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
FIGDIR=opt$figdir
FIGDIRSAMPLE=paste0(FIGDIR, "/", SAMPLE, "/K",k,"/")
FIGDIRTOP=paste0(FIGDIRSAMPLE,"/",SAMPLE,"_K",k,"_dt_", DENSITY.THRESHOLD,"_")
OUTDIRSAMPLE=paste0(OUTDIR, "/", SAMPLE, "/K",k,"/threshold_", DENSITY.THRESHOLD, "/")
FGSEADIR=paste0(OUTDIRSAMPLE,"/fgsea/")
FGSEAFIG=paste0(FIGDIRSAMPLE,"/fgsea/")

message(FIGDIRTOP)

## subscript for files
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
# SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr

 
# create dir if not already
check.dir <- c(OUTDIR, FIGDIR, paste0(FIGDIR,SAMPLE,"/"), paste0(FIGDIR,SAMPLE,"/K",k,"/"), paste0(OUTDIR,SAMPLE,"/"), OUTDIRSAMPLE, FIGDIRSAMPLE, FGSEADIR, FGSEAFIG)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)




######################################################################
## Process topic model results
cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
} else {
    warning(paste0(cNMF.result.file, " does not exist"))
}


## End of data loading



##########################################################################
## Plots


##########################################################################
## topic gene z-score list
pdf(file=paste0(FIGDIRTOP,"top50GeneInTopics.zscore.pdf"), width=4, height=6)
topFeatures.raw.weight <- theta.zscore %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene", variable.name="topic", value.name="scores") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:50)
for ( t in 1:dim(theta)[2] ) {
    toPlot <- data.frame(Gene=topFeatures.raw.weight %>% subset(topic == t) %>% pull(Gene),
                         Score=topFeatures.raw.weight %>% subset(topic == t) %>% pull(scores))
    p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col() + theme_minimal()
    p <- p + coord_flip() + xlab("Top 50 Genes") + ylab("z-score (Specificity)") + ggtitle(paste(SAMPLE, ", Topic ", t, sep="")) + mytheme
    print(p)
}
dev.off()


##########################################################################
## Topic's top gene list, ranked by raw weight
pdf(file=paste0(FIGDIRTOP,"top50GeneInTopics.rawWeight.pdf"), width=4, height=6)
topFeatures <- theta %>% as.data.frame() %>% mutate(genes=rownames(.)) %>% melt(id.vars="genes",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:50)
for ( t in 1:dim(theta)[2] ) {
    toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
                         Score=topFeatures %>% subset(topic == t) %>% pull(scores))
    p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col() + theme_minimal()
    p <- p + coord_flip() + xlab("Top 50 Genes") + ylab("Raw Score (gene's weight in topic)") + ggtitle(paste(SAMPLE, ", Topic ", t, sep="")) + mytheme
    print(p)
}
dev.off()


##########################################################################
## raw program TPM list with annotataion
pdf(file=paste0(FIGDIRTOP,"top10GeneInTopics.rawWeight.pdf"), width=4.5, height=5)
topFeatures <- theta %>% as.data.frame() %>% mutate(genes=rownames(.)) %>% melt(id.vars="genes",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:10)
for ( t in 1:dim(theta)[2] ) {
    toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
                         Score=topFeatures %>% subset(topic == t) %>% pull(scores)) # %>%
        ## merge(., gene.def.pathways, by="Gene", all.x=T)
    ## toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
    p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(width=0.5, fill="#38b4f7") + theme_minimal()
    p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("Raw Score (gene's weight in topic)") +
        mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ", K = ", k, ", Topic ", t))
    print(p4)
}
dev.off()



##########################################################################
## raw program zscore list (top 10)  (can potentially include annotation)
pdf(file=paste0(FIGDIRTOP,"top10GeneInTopics.zscore.pdf"), width=4.5, height=5)
topFeatures <- theta.zscore %>% as.data.frame() %>% mutate(genes=rownames(.)) %>% melt(id.vars="genes",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:10)
for ( t in 1:dim(theta)[2] ) {
    toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
                         Score=topFeatures %>% subset(topic == t) %>% pull(scores)) # %>%
        ## merge(., gene.def.pathways, by="Gene", all.x=T)
    ## toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
    p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(width=0.5, fill="#38b4f7") + theme_minimal()
    p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("z-score") +
        mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ", K = ", k, ", Topic ", t))
    print(p4)
}
dev.off()

