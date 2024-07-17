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
library(conflicted)
conflict_prefer("Position", "base")
packages <- c("optparse","dplyr", "ggplot2", "reshape2", "ggrepel", "conflicted", "gplots", "org.Hs.eg.db")
## library(Seurat)
loadPackages <- function(package) suppressPackageStartupMessages(require(package, character.only=T))
loadedPackages.boolean <- lapply(packages, loadPackages)
if(loadedPackages.boolean %>% as.numeric %>% sum != length(loadedPackages.boolean)) warning(paste0(packages[!(loadedPackages.boolean %>% unlist())] %>% paste0(collapse = ", "), " are not loaded"))
## xfun::pkg_attach(packages)
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
  make_option("--recompute", type="logical", default=F, help="T for recomputing statistical tests and F for not recompute"),
  make_option("--organism", type="character", default="human", help="Identify organism to load database for ENSEMBL and Gene Symbol conversion")
  
)
opt <- parse_args(OptionParser(option_list=option.list))

## ## all genes directories (for sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/"
## opt$K.val <- 60

## ## mouse ENCODE adrenal data sdev
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230117_snakemake_mouse_ENCODE_adrenal/figures/top2000VariableGenes/"
## opt$sampleName <- "mouse_ENCODE_adrenal"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230117_snakemake_mouse_ENCODE_adrenal/analysis/top2000VariableGenes"
## opt$K.val <- 60

## ## mouse ENCODE heart data sdev
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/figures/top2000VariableGenes/"
## opt$sampleName <- "mouse_ENCODE_heart"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/analysis/top2000VariableGenes"
## opt$K.val <- 15

## ## K562 gwps sdev
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/"
## opt$sampleName <- "WeissmanK562gwps"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/"
## opt$K.val <- 20

## ## brain cell data MB
## opt$sampleName <- "DopEx_Sox6_Aldh1a1"
## opt$outdir <- "/broad/macosko/helenk/brain_cell_data/240417_snakemake_MB/analysis/all_genes/"
## opt$figdir <- "/broad/macosko/helenk/brain_cell_data/240417_snakemake_MB/figures/all_genes/"
## opt$K.val <- 60
## opt$density.thr <- 0.2

## ## PD DA neurons
## opt$sampleName <- "dan"
## opt$outdir <- "/broad/macosko/ferris/NMF/240606_snakemake_nurr_da/analysis/top2000VariableGenes/"
## opt$figdir <- "/broad/macosko/ferris/NMF/240606_snakemake_nurr_da/figures/top2000VariableGenes/"
## opt$K.val <- 50
## opt$density.thr <- 0.2



mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))
mytheme <- theme_classic() + theme(axis.text = element_text(size = 7),
                                   axis.title = element_text(size = 8),
                                   plot.title = element_text(hjust = 0.5, face = "bold", size=10),
                                   axis.line = element_line(color = "black", size = 0.25),
                                   axis.ticks = element_line(color = "black", size = 0.25),
                                   legend.key.size = unit(10, units="pt"),
                                   legend.text = element_text(size=7),
                                   legend.title = element_text(size=8)
                                   )


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
## FGSEADIR=paste0(OUTDIRSAMPLE,"/fgsea/")
## FGSEAFIG=paste0(FIGDIRSAMPLE,"/fgsea/")

message(FIGDIRTOP)

## subscript for files
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
# SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr

 
# create dir if not already
check.dir <- c(OUTDIR, FIGDIR, paste0(FIGDIR,SAMPLE,"/"), paste0(FIGDIR,SAMPLE,"/K",k,"/"), paste0(OUTDIR,SAMPLE,"/"), OUTDIRSAMPLE, FIGDIRSAMPLE)
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


## (patch up cNMF_analysis.R ENSGID to Gene conversion)
## gene mapping function
map.ENSGID.SYMBOL <- function(topFeatures) {
    gene.type <- ifelse(sum(as.numeric(colnames(topFeatures) %in% "Gene")) > 0,
                        ##(median.spectra.zscore.df) == sum(as.numeric(grepl("^ENS", median.spectra.zscore %>% rownames))),
                        "Gene",
                        "ENSGID")
    db <- ifelse(grepl("mouse", SAMPLE) | grepl("mouse", opt$organism), "org.Mm.eg.db", "org.Hs.eg.db")
    library(db, character.only = T)
    if(gene.type == "Gene") {
        if (nrow(topFeatures) == sum(as.numeric(grepl("^ENS", topFeatures$Gene)))) {
            ## put median spectra zscore into ENSGID format for PoPS
            mapped.genes <- mapIds(get(db), keys=topFeatures$Gene, keytype = "ENSEMBL", column = "SYMBOL")
            topFeatures <- topFeatures %>%
                mutate(ENSGID = Gene) %>%
                as.data.frame %>%
                mutate(Gene = mapped.genes)
            na.index <- which(is.na(topFeatures$Gene))
            if(length(na.index) > 0) topFeatures$Gene[na.index] <- topFeatures$ENSGID[na.index]
        }

    } else {
        mapped.genes <- mapIds(get(db), keys=topFeatures$ENSGID, keytype = "ENSEMBL", column = "SYMBOL")
        topFeatures <- topFeatures %>%
            as.data.frame %>%
            mutate(Gene = mapped.genes)
        
    }
    return(topFeatures)
}

## end of ENSGID to Gene conversion



##########################################################################
## Plots


##########################################################################
## topic gene z-score list
pdf(file=paste0(FIGDIRTOP,"top50GeneInTopics.zscore.pdf"), width=2.5, height=4.5)
topFeatures.raw.weight <- theta.zscore %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene", variable.name="topic", value.name="scores") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:50) %>% map.ENSGID.SYMBOL
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
pdf(file=paste0(FIGDIRTOP,"top50GeneInTopics.rawWeight.pdf"), width=2.5, height=4.5)
topFeatures <- theta %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:50) %>% map.ENSGID.SYMBOL
for ( t in 1:dim(theta)[2] ) {
    toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(Gene),
                         Score=topFeatures %>% subset(topic == t) %>% pull(scores))
    p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col() + theme_minimal()
    p <- p + coord_flip() + xlab("Top 50 Genes") + ylab("Raw Score (gene's weight in topic)") + ggtitle(paste(SAMPLE, ", Topic ", t, sep="")) + mytheme
    print(p)
}
dev.off()


##########################################################################
## raw program TPM list with annotataion
pdf(file=paste0(FIGDIRTOP,"top10GeneInTopics.rawWeight.pdf"), width=2.5, height=3)
topFeatures <- theta %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:10) %>% map.ENSGID.SYMBOL
for ( t in 1:dim(theta)[2] ) {
    toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(Gene),
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
pdf(file=paste0(FIGDIRTOP,"top10GeneInTopics.zscore.pdf"), width=2.5, height=3)
topFeatures <- theta.zscore %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:10) %>% map.ENSGID.SYMBOL
for ( t in 1:dim(theta)[2] ) {
    toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(Gene),
                         Score=topFeatures %>% subset(topic == t) %>% pull(scores)) # %>%
        ## merge(., gene.def.pathways, by="Gene", all.x=T)
    ## toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
    p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(width=0.5, fill="#38b4f7") + theme_minimal()
    p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("z-score") +
        mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, "\nK = ", k, ", Topic ", t))
    print(p4)
}
dev.off()



##########################################################################
## median spectra list (top 10)  (can potentially include annotation)
pdf(file=paste0(FIGDIRTOP,"top10GeneInTopics.median_spectra.zscore.pdf"), width=2.5, height=3)
topFeatures <- median.spectra.zscore.df %>% as.data.frame() %>% group_by(ProgramID) %>% arrange(desc(median.spectra.zscore)) %>% slice(1:10)  %>% map.ENSGID.SYMBOL
for ( t in 1:dim(theta)[2] ) {
    ProgramID.here <- paste0("K", k, "_", t)
    toPlot <- data.frame(Gene=topFeatures %>% subset(ProgramID == ProgramID.here) %>% pull(Gene),
                         Score=topFeatures %>% subset(ProgramID == ProgramID.here) %>% pull(median.spectra.zscore)) # %>%
        ## merge(., gene.def.pathways, by="Gene", all.x=T)
    ## toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
    p7 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(width=0.5, fill="#38b4f7") + theme_minimal()
    p7 <- p7 + coord_flip() + xlab("Top 10 Genes") + ylab("z-score") +
        mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ",\nK = ", k, ", Program ", t, "\nMedian Spectra"))
    print(p7)
}
dev.off()



##########################################################################
## median spectra raw list (top 10)  (can potentially include annotation)
pdf(file=paste0(FIGDIRTOP,"top10GeneInTopics.median_spectra.raw.pdf"), width=2.5, height=3)
topFeatures <- median.spectra %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene",value.name="Score", variable.name="ProgramID") %>% group_by(ProgramID) %>% arrange(desc(Score)) %>% slice(1:10) %>% mutate(ProgramID = paste0("K", k, "_", ProgramID)) %>% map.ENSGID.SYMBOL
for ( t in 1:dim(theta)[2] ) {
    ProgramID.here <- paste0("K", k, "_", t)
    toPlot <- data.frame(Gene=topFeatures %>% subset(ProgramID == ProgramID.here) %>% pull(Gene),
                         Score=topFeatures %>% subset(ProgramID == ProgramID.here) %>% pull(Score)) # %>%
        ## merge(., gene.def.pathways, by="Gene", all.x=T)
    ## toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
    p8 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(width=0.5, fill="#38b4f7") + theme_minimal()
    p8 <- p8 + coord_flip() + xlab("Top 10 Genes") + ylab("Weight in Program") +
        mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ",\nK = ", k, ", Program ", t, "\nMedian Spectra"))
    print(p8)
}
dev.off()



##########################################################################
## Topic's top gene list, ranked by median spectra zscore
pdf(file=paste0(FIGDIRTOP,"top50GeneInProgram.median_spectra.zscore.pdf"), width=2.5, height=4.5)
topFeatures <- median.spectra.zscore.df %>% as.data.frame() %>% group_by(ProgramID) %>% arrange(desc(median.spectra.zscore)) %>% slice(1:50) %>% map.ENSGID.SYMBOL
for ( t in 1:dim(theta)[2] ) {
    ProgramID.here <- paste0("K", k, "_", t)
    toPlot <- data.frame(Gene=topFeatures %>% subset(ProgramID == ProgramID.here) %>% pull(Gene),
                         Score=topFeatures %>% subset(ProgramID == ProgramID.here) %>% pull(median.spectra.zscore)) # %>%
        ## merge(., gene.def.pathways, by="Gene", all.x=T)
    ## toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
    p9 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(width=0.5, fill="gray30") + theme_minimal()
    p9 <- p9 + coord_flip() + xlab("Top 50 Genes") + ylab("z-score") +
        mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ",\nK = ", k, ", Program ", t, "\nMedian Spectra"))
    print(p9)
}
dev.off()


##########################################################################
## median spectra raw list (top 50)  (can potentially include annotation)
pdf(file=paste0(FIGDIRTOP,"top50GeneInProgram.median_spectra.raw.pdf"), width=2.5, height=4.5)
topFeatures <- median.spectra %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene",value.name="Score", variable.name="ProgramID") %>% group_by(ProgramID) %>% arrange(desc(Score)) %>% slice(1:50) %>% mutate(ProgramID = paste0("K", k, "_", ProgramID)) %>% map.ENSGID.SYMBOL
for ( t in 1:dim(theta)[2] ) {
    ProgramID.here <- paste0("K", k, "_", t)
    toPlot <- data.frame(Gene=topFeatures %>% subset(ProgramID == ProgramID.here) %>% pull(Gene),
                         Score=topFeatures %>% subset(ProgramID == ProgramID.here) %>% pull(Score)) # %>%
        ## merge(., gene.def.pathways, by="Gene", all.x=T)
    ## toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
    p10 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(width=0.5, fill="gray30") + theme_minimal()
    p10 <- p10 + coord_flip() + xlab("Top 50 Genes") + ylab("Weight in Program") +
        mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ",\nK = ", k, ", Program ", t, "\nMedian Spectra"))
    print(p10)
}
dev.off()


##########################################################################
## topic Pearson correlation heatmap
## remove NA from theta.zscore
tokeep <- (!is.na(theta.zscore)) %>% apply(1, sum) == k ## remove genes with NA ## why is there NA?
d <- cor(theta.zscore[tokeep,], method="pearson")
m <- as.matrix(d)

## Function for plotting heatmap  # new version (adjusted font size)
plotHeatmap <- function(mtx, labCol, title, margins=c(12,6), ...) { #original
  heatmap.2(
    mtx %>% t(), 
    Rowv=T, 
    Colv=T,
    trace='none',
    key=T,
    col=palette,
    labCol=labCol,
    margins=margins, 
    cex.main=0.8, 
    cexCol=4.8/sqrt(ncol(mtx)), cexRow=4.8/sqrt(ncol(mtx)), #4.8/sqrt(nrow(mtx))
    ## cexCol=1/(ncol(mtx)^(1/3)), cexRow=1/(ncol(mtx)^(1/3)), #4.8/sqrt(nrow(mtx))
    main=title,
    ...
  )
}

pdf(file=paste0(FIGDIRTOP, "topic.Pearson.correlation.pdf"))
plotHeatmap(m, labCol=rownames(m), margins=c(3,3), title=paste0("cNMF, topic zscore clustering by Pearson Correlation"))
dev.off()


tokeep <- (!is.na(median.spectra)) %>% apply(1, sum) == k ## remove genes with NA ## why is there NA?
d <- cor(median.spectra[tokeep,], method="pearson")
m <- as.matrix(d)

pdf(file=paste0(FIGDIRTOP, "topic.median_spectra.Pearson.correlation.pdf"))
plotHeatmap(m, labCol=rownames(m), margins=c(3,3), title=paste0("cNMF, Program median spectra clustering by Pearson Correlation"))
dev.off()
