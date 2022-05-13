## Helen Kang
## Motif enrichment plots
## 210503

packages <- c("optparse","dplyr", "ggplot2", "reshape2", "ggrepel", "conflicted")
xfun::pkg_attach(packages)
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")



option.list <- list(
  make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/figures/all_genes/", help="Figure directory"),
  make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/", help="Output directory"),
  make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/", help="Input 10x data directory"),
  make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
  make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
  make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
  
  #summary plot parameters
  make_option("--test.type", type="character", default="per.guide.wilcoxon", help="Significance test to threshold perturbation results"),
  make_option("--ep.type", type="character", default="enhancer", help="motif enrichment for enhancer or promoter, specify 'enhancer' or 'promoter'"),
  make_option("--adj.p.value.thr", type="numeric", default=0.05, help="adjusted p-value threshold"),
  make_option("--recompute", type="logical", default=F, help="T for recomputing statistical tests and F for not recompute"),
  make_option("--motif.match.thr.str", type="character", default="pval1e-6", help="threshold for subsetting motif matches")
  
)
opt <- parse_args(OptionParser(option_list=option.list))

## ## all genes directories (for sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/"
## opt$K.val <- 60

mytheme <- theme_classic() + theme(axis.text = element_text(size = 7),
                                   axis.title = element_text(size = 8),
                                   plot.title = element_text(hjust = 0.5, face = "bold", size=8))

SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
# STATIC.SAMPLE=c("Telo_no_IL1B_T200_1", "Telo_no_IL1B_T200_2", "Telo_plus_IL1B_T200_1", "Telo_plus_IL1B_T200_2", "no_IL1B", "plus_IL1B",  "pooled")
# DATADIR=opt$olddatadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
## TMDIR=opt$topic.model.result.dir
## SEP=opt$sep
# K.list <- strsplit(opt$K.list,",") %>% unlist() %>% as.numeric()
k <- opt$K.val
num.top.genes <- 300 ## number of top topic defining genes
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
# SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr
motif.match.thr.str <- opt$motif.match.thr.str

## ## directories for factor motif enrichment
## FILENAME=opt$filename


## ## modify motif.enhancer.background input directory ##HERE: perhaps do a for loop for all the desired thresholds (use strsplit on enhancer.fimo.threshold)
## opt$motif.enhancer.background <- paste0(opt$motif.enhancer.background, opt$enhancer.fimo.threshold, "/fimo.formatted.tsv")

 
# create dir if not already
check.dir <- c(OUTDIR, FIGDIR, paste0(FIGDIR,SAMPLE,"/"), paste0(FIGDIR,SAMPLE,"/K",k,"/"), paste0(OUTDIR,SAMPLE,"/"), OUTDIRSAMPLE, FIGDIRSAMPLE, FGSEADIR, FGSEAFIG)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)



######################################################################
## load data
cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
} else {
    warning(paste0(cNMF.result.file, " does not exist"))
}

## load motif enrichment results
all.ttest.df.path <- paste0(OUTDIRSAMPLE,"/", ep.type, ".topic.top.", num.top.genes, ".zscore.gene_motif.count.ttest.enrichment_motif.thr.", motif.match.thr.str, "_", SUBSCRIPT.SHORT,".txt")
ttest.df <- read.delim(all.ttest.df.path, stringsAsFactors=F)



# file.name <- paste0(OUTDIRSAMPLE,"/cNMFAnalysis.factorMotifEnrichment.",SUBSCRIPT.SHORT,".RData")
# print(file.name)
# if(file.exists((file.name))) { 
#     load(file.name)
#     print(paste0("loading ", file.name))
# }
# motif.enrichment.variables <- c("all.enhancer.fisher.df", "all.promoter.fisher.df", 
#                                 "promoter.wide", "enhancer.wide", "promoter.wide.binary", "enhancer.wide.binary",
#                                 "enhancer.wide.10en6", "enhancer.wide.binary.10en6", "all.enhancer.fisher.df.10en6",
#                                 "promoter.wide.10en6", "promoter.wide.binary.10en6", "all.promoter.fisher.df.10en6",
#                                 "all.promoter.ttest.df", "all.promoter.ttest.df.10en6", "all.enhancer.ttest.df", "all.enhancer.ttest.df.10en6")
# motif.enrichment.variables.missing <- (!(motif.enrichment.variables %in% ls())) %>% as.numeric %>% sum 
# if ( motif.enrichment.variables.missing > 0 ) {
#     warning(paste0(motif.enrichment.variables[!(motif.enrichment.variables %in% ls())], " not available"))
# }


## End of data loading



##########################################################################
## Plots


## volcano plots
volcano.plot <- function(toplot, ep.type, ranking.type, label.type="") {
    if( label.type == "pos") {
        label <- toplot %>% subset(two.sided.p.adjust < fdr.thr & enrichment.log2fc > 0) %>% mutate(motif.toshow = gsub("HUMAN.H11MO.", "", motif))
    } else {
        label <- toplot %>% subset(two.sided.p.adjust < fdr.thr) %>% mutate(motif.toshow = gsub("HUMAN.H11MO.", "", motif))
    }
    t <- gsub("topic_", "", toplot$topic[1])
    p <- toplot %>% ggplot(aes(x=enrichment.log2fc, y=-log10(two.sided.p.adjust))) + geom_point(size=0.5) + mytheme +
        ggtitle(paste0(SAMPLE[1], " Topic ", t, " Top ", num.top.genes, " ", ranking.type,"\n", ifelse(ep.type=="promoter", "Promoter", "Enhancer"), " Motif Enrichment")) + xlab("Motif Enrichment (log2FC)") + ylab("-log10(adjusted p-value)") +
        xlim(0,max(toplot$enrichment.log2fc)) +
        geom_hline(yintercept=-log10(fdr.thr), linetype="dashed", color="gray") +
        geom_text_repel(data=label, box.padding = 0.25,
                        aes(label=motif.toshow), size=2.5,
                        max.overlaps = 15,
                        color="black")# + theme(text=element_text(size=16), axis.title=element_text(size=16), axis.text=element_text(size=16), plot.title=element_text(size=14))
    print(p)
    p <- toplot %>% ggplot(aes(x=enrichment.log2fc, y=-log10(p.value))) + geom_point(size=0.25) + mytheme +
        ggtitle(paste0(SAMPLE[1], " Topic ", t, " Top ", num.top.genes," ", ranking.type,"\n", ifelse(ep.type=="promoter", "Promoter", "Enhancer"), " Motif Enrichment")) + xlab("Motif Enrichment (log2FC)") + ylab("-log10(p-value)") +
        xlim(0,max(toplot$enrichment.log2fc)) +
        geom_hline(yintercept=-log10(fdr.thr), linetype="dashed", color="gray") +
        geom_text_repel(data=label, box.padding = 0.25,
                        aes(label=motif.toshow), size=2.5,
                        max.overlaps = 15,
                        color="black") #+ theme(text=element_text(size=16), axis.title=element_text(size=16), axis.text=element_text(size=16), plot.title=element_text(size=14))
    return(p)
}

## function for all volcano plots
all.volcano.plots <- function(all.fisher.df, ep.type, ranking.type, label.type="") {
    for ( t in 1:k ){
        toplot <- all.fisher.df %>% subset(topic==paste0("topic_",t))
        volcano.plot(toplot, ep.type, ranking.type, label.type) %>% print()
    }
}


##########################################################################
## motif enrichment plot
pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.count.ttest.enrichment_motif.thr.", motif.match.thr.str, ".pdf"), width=3, height=3)
all.volcano.plots(get(paste0("ttest.df")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score")
dev.off()
    

# ep.names <- c("enhancer", "promoter")
# for (ep.type in ep.names) {
#     pdf(file=paste0(FIGDIRTOP,"zscore.",ep.type,".motif.enrichment.pdf"))
#     all.volcano.plots(get(paste0("all.",ep.type,".fisher.df")), ep.type, ranking.type="z-score")
#     dev.off()
#     pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment_motif.thr.10e-6.pdf"))
#     all.volcano.plots(get(paste0("all.",ep.type,".fisher.df.10en6")), ep.type, ranking.type="z-score")
#     dev.off()
#     pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest.pdf"))
#     all.volcano.plots(get(paste0("all.",ep.type,".ttest.df")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score")
#     dev.off() 
#     pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest_motif.thr.10e-6.pdf"))
#     all.volcano.plots(get(paste0("all.",ep.type,".ttest.df.10en6")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score")
#     dev.off() 
#     pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest.labelPos.pdf"), width=6, height=6)
#     all.volcano.plots(get(paste0("all.",ep.type,".ttest.df")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score", label.type="pos")
#     dev.off() 
#     pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest_motif.thr.10e-6.labelPos.pdf"))
#     all.volcano.plots(get(paste0("all.",ep.type,".ttest.df.10en6")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score", label.type="pos")
#     dev.off() 
# }

