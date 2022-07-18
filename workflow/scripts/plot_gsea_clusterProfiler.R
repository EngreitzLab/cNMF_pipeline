## Helen Kang
## Plot enriched GO term for each component
## 220718

 
suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(reshape2)
    library(ggplot2)
    library(cowplot)
    library(ggpubr) ## ggarrange
    library(gplots) ## heatmap.2
    library(scales) ## geom_tile gradient rescale
    library(ggrepel)
    library(stringr)
    library(svglite)
    library(ggseqlogo)
    library(universalmotif)
})



##########################################################################################
## Constants and Directories
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211205_sig_topic_TPM/figures/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/"
option.list <- list(
    make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/", help="Figure directory"),
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/", help="Output directory"),
    make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
    make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--raw.mtx.RDS.dir",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.RDS", help="input matrix to cNMF pipeline"),
    make_option("--adj.p.value.thr", type="numeric", default=0.1, help="adjusted p-value threshold"),
    
    ## GSEA parameters
    make_option("--ranking.type", type="character", default="zscore", help="{zscore, raw} ranking for the top program genes"),
    make_option("--GSEA.type", type="character", default="GOEnrichment", help="{GOEnrichment, ByWeightGSEA, GSEA}")
)
opt <- parse_args(OptionParser(option_list=option.list))


OUTDIR <- opt$outdir
FIGDIR <- opt$figdir
SAMPLE <- opt$sampleName
k <- opt$K.val
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
FIGDIRSAMPLE=paste0(FIGDIR, "/", SAMPLE, "/K",k,"/")
FIGDIRTOP=paste0(FIGDIRSAMPLE,"/",SAMPLE,"_K",k,"_dt_", DENSITY.THRESHOLD,"_")
OUTDIRSAMPLE=paste0(OUTDIR, "/", SAMPLE, "/K",k,"/threshold_", DENSITY.THRESHOLD, "/")
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
# SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## GSEA specific parameters
ranking.type.here <- opt$ranking.type
GSEA.type <- opt$GSEA.type


message(FIGDIRTOP)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr

# create dir if not already
check.dir <- c(OUTDIR, FIGDIR, paste0(FIGDIR,SAMPLE,"/"), paste0(FIGDIR,SAMPLE,"/K",k,"/"), paste0(OUTDIR,SAMPLE,"/"), OUTDIRSAMPLE, FIGDIRSAMPLE)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

mytheme <- theme_classic() + theme(axis.text = element_text(size = 7),
                                   axis.title = element_text(size = 8),
                                   plot.title = element_text(hjust = 0.5, face = "bold", size=10),
                                   axis.line = element_line(color = "black", size = 0.25),
                                   axis.ticks = element_line(color = "black", size = 0.25),
                                   legend.key.size = unit(10, units="pt"),
                                   legend.text = element_text(size=7),
                                   legend.title = element_text(size=8)
                                   )

file.name <- paste0(OUTDIRSAMPLE, "/clusterProfiler_GeneRankingType", ranking.type.here, "_EnrichmentType", GSEA.type,".txt")
gsea.df <- read.delim(file.name, stringsAsFactors=F)
toplot <- gsea.df %>%
    subset(fdr.across.ont < fdr.thr) %>%
    group_by(ProgramID) %>%
    arrange(fdr.across.ont) %>%
    unique %>%
    slice(1:10) %>%
    mutate(TruncatedDescription = str_trunc(paste0(ID, "; ", Description), width=50, side="right"),
           t = gsub("K60_", "", ProgramID) %>% as.numeric) %>%
    arrange(t, fdr.across.ont) %>%
    as.data.frame

plot.title <- paste0(ifelse(grepl("GO", GSEA.type), "GO Term Enrichment", "MSigDB Pathway Enrichment"),
                     "\non ",
                     ifelse(grepl("zscore", ranking.type.here), "Program Gene Specificity", "Raw Weight"),
                     "\nby ",
                     ifelse(grepl("ByWeight", GSEA.type), "All Gene Weight", "Top 300 Gene Set"))


pdf(file=paste0(FIGDIRTOP,"top10EnrichedPathways_GeneRankingType", ranking.type.here, "_EnrichmentType", GSEA.type, ".pdf"), width=4, height=4)
for(program in (toplot$ProgramID %>% unique)) {
    t <- strsplit(program, split="_") %>% sapply(`[[`,2)
    toplot.here <- toplot %>%
        subset(ProgramID %in% program) %>%
        arrange(fdr.across.ont)
    labels <- toplot.here$TruncatedDescription %>% unique %>% rev
    toplot.here <- toplot.here %>%
        mutate(TruncatedDescription = factor(TruncatedDescription, levels = labels))
    p <- toplot.here %>% ggplot(aes(x=TruncatedDescription, y=-log10(fdr.across.ont))) + geom_col(fill="gray") + coord_flip() + mytheme +
        xlab(ifelse(grepl("GO", GSEA.type), "GO Terms", "Pathways")) + ylab("FDR (-log10)") +
        ggtitle(paste0(plot.title, "\nProgram ", t))    
    print(p)
}
dev.off()
