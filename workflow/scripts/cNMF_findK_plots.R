## Find K plots
## Helen Kang
## 210715

## .libPaths("/home/groups/engreitz/Software/R_3.6.1")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(gplots))

## source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")

option.list <- list(
  make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/figures/all_genes/2kG.library/acrossK/", help="Figure directory"),
  make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/acrossK/", help="Output directory"),
  make_option("--reference.table", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/210702_2kglib_adding_more_brief_ca0713.xlsx"),
  make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
  make_option("--p.adj.threshold", type="numeric", default=0.1, help="Threshold for fdr and adjusted p-value"),
  make_option("--aggregated.data", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/acrossK/aggregated.outputs.findK.RData")
)
opt <- parse_args(OptionParser(option_list=option.list))


## ## for all genes (in sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/2kG.library/acrossK/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/acrossK/"
## opt$aggregated.data <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/acrossK//aggregated.outputs.findK.RData"


## ## ## for testing cNMF_ pipeline
## ## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/all_genes/2kG.library/acrossK/"
## opt$aggregated.data <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/all_genes/FT010_fresh_2min/acrossK/aggregated.outputs.findK.RData"


## ## for testing cNMF_pipeline with FT010_fresh_4min
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/figures/all_genes/FT010_fresh_4min/acrossK/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/all_genes/FT010_fresh_4min/acrossK/"
## opt$sampleName <- "FT010_fresh_4min"
## opt$aggregated.data <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/all_genes/FT010_fresh_4min/acrossK/aggregated.outputs.findK.RData"


## ## for testing findK_plots for scRNAseq_2kG_11AMDox_1
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/figures/all_genes/scRNAseq_2kG_11AMDox_1/acrossK/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/acrossK/"
## opt$sampleName <- "scRNAseq_2kG_11AMDox_1"
## opt$aggregated.data <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/acrossK/aggregated.outputs.findK.RData"


## Directories and Constants
SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
DATADIR=opt$datadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
FIGDIR=opt$figdir
check.dir <- c(OUTDIR, FIGDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))


## load("/Volumes/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2105_findK/analysis/no_IL1B/aggregated.outputs.findK.RData")
## load("/Volumes/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/differential_expression/210526_SeuratDE/outputs/no_IL1B/de.markers.RData")
mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)

## Load data
load(opt$aggregated.data) ## all.fdr.df, all.test.df, enhancer.fisher.df, count.by.GWAS, fgsea.results.df, promoter.fisher.df
## ref.table <- read_xlsx(opt$reference.table, sheet="2000_gene_library_annotated") 


## ## Update the files for 2kG library
## ## Intersecting EdgeR and Topic Model Results
## threshold <- opt$p.adj.threshold
## # % of perturbations with a large number of DE genes in the per-gene analysis that ALSO have significantly DE topics  (p-value 0.005, logFC 1.2)
## strict.DE.all <- read_excel(path="/Users/helenkang/Documents/EngreitzLab/Pertube-seq_Analysis/200_gene_lib_EdgeR_DE.xlsx", sheet="EdgeR_DE_p.001_lfc1.2")
## lenient.DE.all <- read_excel(path="/Users/helenkang/Documents/EngreitzLab/Pertube-seq_Analysis/200_gene_lib_EdgeR_DE.xlsx", sheet="EdgeR_DE_p.01_lfc1.15")

## strict.DE <- strict.DE.all %>% select(Target...13, `NoI_#DE_genes`, `PlusI_#DE_genes`) %>% `colnames<-`(c("Gene", "NoI.num.genes", "PlusI.num.genes")) %>% mutate(cutoff.type = "strict", cutoff.detail = "p < 0.001, logFC > 1.2")
## lenient.DE <- lenient.DE.all %>% select(Target...13, NoI_DE...14, PlusI_DE...15) %>% `colnames<-`(c("Gene", "NoI.num.genes", "PlusI.num.genes")) %>% mutate(cutoff.type = "lenient", cutoff.detail = "p < 0.01, logFC > 1.15")
## edgeR.DE <- rbind(strict.DE, lenient.DE)

## perturbation.test.stat.df <- all.test.df %>% subset(adjusted.p.value < threshold) %>% select(K, Gene, test.type, Topic) %>% group_by(K, Gene, test.type) %>% summarize(topic.count = n())
## edgeR.test.df <- merge(edgeR.DE, perturbation.test.stat.df, by="Gene")


## MSigDB Pathway Enrichment 
## Number of GO enrichment per model


## Notes 210518
## also output pdf
## plot top enriched MSigDB pathways for each topic and K
## also consider eps file
## make sure that p.adjust is caluclated over all topics for a model


## plot number of GO enrichment on raw score ranking per K
par(mar = c(4, 4, .1, .1))
threshold = opt$p.adj.threshold
pdf(file=paste0(FIGDIR, "/GO.enrichment.on.raw.score.ranking.threshold", as.character(threshold), ".pdf"))
fgsea.summary.df <- fgsea.results.df %>% subset(padj.over.topics < threshold & type == "raw.score") %>% group_by(database,K) %>% summarize(total.per.K.count = n()) %>% mutate(average.topic.count.per.K = total.per.K.count / K)
K <- fgsea.summary.df %>% pull(K) %>% unique() # K for x tick labels
# GO only plot for sum for sum of significant pathways
p <- fgsea.summary.df %>% subset(database == "msigdb.c3") %>% ggplot(aes(x=K, y = total.per.K.count)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of GO Pathways Under p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +
  ggtitle("FGSEA on raw weights of each gene in the topics")
print(p)

# GO only plot for average number of significant pathways per topic
p <- fgsea.summary.df %>% subset(database == "msigdb.c3") %>% ggplot(aes(x=K, y = average.topic.count.per.K)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of GO Pathways Per Topic \nUnder p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +
  ggtitle("FGSEA on raw weights of each gene in the topics")
print(p)

# Number of selected MSigDB pathways per model
p <- fgsea.summary.df %>% ggplot(aes(x=K, y = total.per.K.count, color = database)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of Pathways Under p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +   
  ggtitle("FGSEA on raw weights of each gene in the topics")
print(p)

# selected MSigDB pathway plot for average number of significant pathways per topic
p <- fgsea.summary.df %>% ggplot(aes(x=K, y = average.topic.count.per.K, color = database)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Pathways Under p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +   
  ggtitle("FGSEA on raw weights of each gene in the topics")
print(p)
dev.off()


pdf(file=paste0(FIGDIR, "/GO.enrichment.on.z-score.ranking.threshold", as.character(threshold), ".pdf"))
fgsea.summary.df <- fgsea.results.df %>% subset(padj < threshold & type == "z.score") %>% group_by(database,K) %>% summarize(total.per.K.count = n()) %>% mutate(average.topic.count.per.K = total.per.K.count / K)
K <- fgsea.summary.df %>% pull(K) %>% unique() # K for x tick labels
# GO only plot for sum for sum of significant pathways
p <- fgsea.summary.df %>% subset(database == "msigdb.c3") %>% ggplot(aes(x=K, y = total.per.K.count)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of GO Pathways Under p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +   
  ggtitle("FGSEA on z-score of each gene in the topics")
print(p)

# GO only plot for average number of significant pathways per topic
p <- fgsea.summary.df %>% subset(database == "msigdb.c3") %>% ggplot(aes(x=K, y = average.topic.count.per.K)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of GO Pathways Per Topic \nUnder p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +
  ggtitle("FGSEA on z-score of each gene in the topics")
print(p)

# Number of selected MSigDB pathways per model
p <- fgsea.summary.df %>% ggplot(aes(x=K, y = total.per.K.count, color = database)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of Pathways Under p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +   
  ggtitle("FGSEA on z-score of each gene in the topics")
print(p)

# selected MSigDB pathway plot for average number of significant pathways per topic
p <- fgsea.summary.df %>% ggplot(aes(x=K, y = average.topic.count.per.K, color = database)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Pathways Under p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +   
  ggtitle("FGSEA on z-score of each gene in the topics")
print(p)
dev.off()


## plot number of unique GO enrichment per K
threshold <- opt$p.adj.threshold
fgsea.results.raw.df <- fgsea.results.df %>% subset(padj < threshold & type == "raw.weight")
fgsea.results.zscore.df <- fgsea.results.df %>% subset(padj < threshold & type == "z.score")
## calculate unique pathways per K
fgsea.summary.unique.raw.df <- fgsea.results.raw.df %>% select(K, pathway, database) %>% group_by(K) %>% unique() %>% group_by(database,K) %>% summarize(total.per.K.count = n()) %>% mutate(average.topic.count.per.K = total.per.K.count / K)
fgsea.summary.unique.zscore.df <- fgsea.results.zscore.df %>% select(K, pathway, database) %>% group_by(K) %>% unique() %>% group_by(database,K) %>% summarize(total.per.K.count = n()) %>% mutate(average.topic.count.per.K = total.per.K.count / K)

pdf(file=paste0(FIGDIR,"/GO.unique.p.adj.", threshold, ".pdf"))
# GO only plot for sum of significant pathways
K <- fgsea.summary.unique.raw.df$K %>% unique()
p <- fgsea.summary.unique.raw.df %>% subset(database == "msigdb.c3") %>% ggplot(aes(x=K, y = total.per.K.count)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of Unique GO Pathways Under p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +   
  ggtitle("FGSEA on raw score of each gene in the factors")
print(p)
K <- fgsea.summary.unique.zscore.df$K %>% unique()
p <- fgsea.summary.unique.zscore.df %>% subset(database == "msigdb.c3") %>% ggplot(aes(x=K, y = total.per.K.count)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of Unique GO Pathways Under p.adj of ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K)  +   
  ggtitle("FGSEA on z-score of each gene in the factors")
print(p)
dev.off()


 
## plot number of unique MSigDB pathways per K
# selected MSigDB pathway plot for sum of significant pathways
pdf(file=paste0(FIGDIR, "/MSigDB.unique.p.adj.", threshold, ".pdf"))
p <- fgsea.summary.unique.raw.df %>% ggplot(aes(x=K, y = total.per.K.count, color = database)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of Unique Pathways Under p.adj of ", threshold, ", raw weight")) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K)
print(p)
p <- fgsea.summary.unique.zscore.df %>% ggplot(aes(x=K, y = total.per.K.count, color = database)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of Unique Pathways Under p.adj of ", threshold, ", z-score")) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K)
print(p)
dev.off()


## plot average number of unique GO enrichment per topic vs K
# GO only plot for sum for average number of significant pathways per topic
pdf(file=paste0(FIGDIR, "/GO.avg.unique.over.K_p.adj.", threshold, ".pdf"))
K <- fgsea.summary.unique.raw.df$K %>% unique()
p <- fgsea.summary.unique.raw.df %>% subset(database == "msigdb.c3") %>% ggplot(aes(x=K, y = average.topic.count.per.K)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Unique GO Pathways Per Topic \nUnder p.adj of ", threshold, ", raw weight")) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K)
print(p)
K <- fgsea.summary.unique.zscore.df$K %>% unique()
p <- fgsea.summary.unique.zscore.df %>% subset(database == "msigdb.c3") %>% ggplot(aes(x=K, y = average.topic.count.per.K)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Unique GO Pathways Per Topic \nUnder p.adj of ", threshold, ", z-score")) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K)
print(p)
dev.off()


## plot average number of unique MSigDB pathways per topic vs K
# selected MSigDB pathway plot for average number of significant pathways per topic
pdf(file=paste0(FIGDIR, "/MsigDB.avg.unique.over.K_p.adj.", threshold, ".pdf"))
K <- fgsea.summary.unique.raw.df$K %>% unique()
p <- fgsea.summary.unique.raw.df %>% ggplot(aes(x=K, y = average.topic.count.per.K, color = database)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Unique Pathways Under p.adj of ", threshold, ", raw weight")) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K)
print(p)
K <- fgsea.summary.unique.zscore.df$K %>% unique()
p <- fgsea.summary.unique.zscore.df %>% ggplot(aes(x=K, y = average.topic.count.per.K, color = database)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Unique Pathways Under p.adj of ", threshold, ", z-score")) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K)
print(p)
dev.off()


## ## plot number of significant DE topics per K
## threshold = opt$p.adj.threshold
## pdf(file=paste0(FIGDIR, "/Topics.DE.sig.threshold", as.character(threshold), ".pdf"), width=8, height=6)
## all.test.summary.df <- all.test.df %>% subset(adjusted.p.value < threshold & grepl("wilcoxon", test.type)) %>% group_by(K, test.type) %>% summarize(total=n())
## K <- all.test.summary.df %>% pull(K) %>% unique() # K for x tick labels
## p <- all.test.summary.df %>% ggplot(aes(x=K, y=total, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Number of DE (topic x gene) pairs (adjusted p-value < ", threshold, ")")) + mytheme +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## dev.off()


## ## Number of perturbations with at least one DE topic per K ## need to update DE file
## pdf(file=paste0(FIGDIR, "/Perturbation.DE.sig.threshold", as.character(threshold), ".pdf"), width=8, height=6)
## num.ptb <- all.test.df %>% pull(Gene) %>% unique() %>% length()
## all.test.summary.unique.df <- all.test.df %>% subset(adjusted.p.value < threshold) %>% select(K, test.type, Gene) %>% group_by(K, test.type) %>% unique() %>% summarize(total=n()) %>% mutate(fraction.perturbation = total / num.ptb)
## K <- all.test.summary.unique.df %>% pull(K) %>% unique() # K for x tick labels
## p <- all.test.summary.unique.df %>% ggplot(aes(x=K, y=total, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Number of Perturbations with at least one DE topic \n(adjusted p-value < ", threshold, ")")) + mytheme +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)

## p <- all.test.summary.unique.df %>% ggplot(aes(x=K, y=fraction.perturbation, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Fraction of Perturbations with at least one DE topic \n(adjusted p-value < ", threshold, ")")) + mytheme +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## dev.off()


## ## Fraction of significant topics by adjusted p-value versus K
## pdf(file=paste0(FIGDIR, "/Topics.DE.fraction.sig.threshold", as.character(threshold), ".pdf"), width=8, height=6)
## all.test.fraction.df <- all.test.df %>% subset(adjusted.p.value < threshold) %>% select(K, Topic, test.type) %>% unique() %>% group_by(K, test.type) %>% summarise(topic.count = n()) %>% mutate(topic.fraction = topic.count / K)
## p <- all.test.fraction.df %>% ggplot(aes(x=K, y=topic.fraction, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Fraction of topics with adjusted p-value < ", threshold)) + mytheme +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## dev.off()


## ## Plots number of topics per model based on emprical FDR
## threshold=opt$p.adj.threshold
## pdf(file=paste0(FIGDIR, "/Number.of.DE.topics.fdr.sig.threshold", as.character(threshold), ".pdf"), width=8, height=6)
## all.fdr.summary.df <- all.fdr.df %>% subset(fdr < threshold) %>% group_by(K, test.type) %>% summarize(total=n())
## K <- all.fdr.summary.df %>% pull(K) %>% unique() # K for x tick labels
## p <- all.fdr.summary.df %>% ggplot(aes(x=K, y=total, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Number of DE topics (FDR < ", threshold, ")")) + mytheme +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## dev.off()


## ## Number of perturbations with FDR < 0.1 versus K ## need updated EdgeR DE results
## pdf(file=paste0(FIGDIR, "/Perturbation.DE.fdr.sig.threshold.", as.character(threshold), ".pdf"), width=8, height=6)
## num.ptb <- all.fdr.df %>% pull(Gene) %>% unique() %>% length()
## all.fdr.summary.unique.df <- all.fdr.df %>% subset(fdr < threshold) %>% select(K, test.type, Gene) %>% group_by(K, test.type) %>% unique() %>% summarize(total=n()) %>% mutate(fraction.perturbation = total / num.ptb)
## K <- all.fdr.summary.unique.df %>% pull(K) %>% unique() # K for x tick labels
## p <- all.fdr.summary.unique.df %>% ggplot(aes(x=K, y=total, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Number of Perturbations with at least one DE topic \n(FDR < ", threshold, ")")) + mytheme +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## dev.off()


## ## Fraction of perturbations with FDR < 0.1 versus K
## pdf(file=paste0(FIGDIR, "/Perturbation.DE.fraction.fdr.sig.threshold.", as.character(threshold), ".pdf"), width=8, height=6)
## p <- all.fdr.summary.unique.df %>% ggplot(aes(x=K, y=fraction.perturbation, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Fraction of Perturbations with at least one DE topic \n(FDR < ", threshold, ")")) + mytheme +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## dev.off()


## ## fraction of FDR significant topics versus K
## threshold <- opt$p.adj.threshold
## pdf(file=paste0(FIGDIR, "/Topics.DE.fdr.sig.threshold.", threshold, ".pdf"))
## all.fdr.fraction.df <- all.fdr.df %>% subset(fdr < threshold) %>% select(K, Topic, test.type) %>% unique() %>% group_by(K, test.type) %>% summarise(topic.count = n()) %>% mutate(topic.fraction = topic.count / K)
## p <- all.fdr.fraction.df %>% ggplot(aes(x=K, y=topic.fraction, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Fraction of topics with FDR < ", threshold)) + mytheme +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## dev.off()


## Plots for Transcription Factor Motif Enrichment in ABC Linked Enhancers
threshold <- opt$p.adj.threshold
enrichment.thr <- 0 ### all.enhancer.ttest.df is empty
pdf(file=paste0(FIGDIR, "/TF.motif.enrichment.ABC.linked.enhancers.sig.threshold", as.character(threshold), ".pdf"), width=8, height=6)
## ABC linked enhancers
enhancer.summary.df <- all.enhancer.ttest.df %>% subset(p.adjust < threshold & enrichment > enrichment.thr) %>% group_by(K) %>% summarise(total=n()) %>% mutate(average.per.topic = total / K)
K <- enhancer.summary.df %>% pull(K) %>% unique() # K for x tick labels
# total number of significant TF plots
p1 <- enhancer.summary.df %>% ggplot(aes(x=K, y=total)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Total Number of Transcription Factors\n with adjusted p-value < ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K)
print(p1)
p2 <- enhancer.summary.df %>% ggplot(aes(x=K, y=average.per.topic)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Transcription Factors per Topic \n with adjusted p-value < ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +
  ggtitle(paste0("Transcription Factors Enriched in \nABC Enhancers of the Top 100 Genes of Each Topic "))
print(p2)


## plot number of unique TF enriched per K for ABC linked enhancers
enhancer.summary.unique.df <- all.enhancer.ttest.df %>% subset(p.adjust < threshold & enrichment > enrichment.thr) %>% select(K, motif) %>% unique() %>% group_by(K) %>% summarise(total=n()) %>% mutate(average.per.topic = total / K)
K <- enhancer.summary.unique.df %>% pull(K) %>% unique() # K for x tick labels
# total number of significant TF plots
p3 <- enhancer.summary.unique.df %>% ggplot(aes(x=K, y=total)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of Unique Transcription Factors\n with adjusted p-value < ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +
  ggtitle(paste0("Unique Transcription Factors Enriched in \nABC Enhancers of the Top 100 Genes of Each Topic "))
print(p3)
p4 <- enhancer.summary.unique.df %>% ggplot(aes(x=K, y=average.per.topic)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Unique Transcription Factors per Topic \n with adjusted p-value < ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +
  ggtitle(paste0("Unique Transcription Factors Enriched in \nABC Enhancers of the Top 100 Genes of Each Topic "))
print(p4)

p <- ggarrange(p1 + ylab("Total TF count"), p2 + ggtitle("") + ylab("Total TF per factor"),
               p3 + ggtitle("") + ylab("Unique TF count"), p4 + ggtitle("") + ylab("Unique TF per factor"), nrow=4)
print(p)
dev.off()



## Notes:
## why is there a sharp change from K=19 to K=21?
## Fraction of topics that has at least one significant (metric) versus K

## plot number of TF enriched per K for promoters}
pdf(file=paste0(FIGDIR, "/TF.motif.enrichment.promoters.sig.threshold", as.character(threshold), ".pdf"), width=8, height=6)
## promoters
promoter.summary.df <- all.promoter.ttest.df %>% subset(p.adjust < threshold & enrichment > enrichment.thr) %>% group_by(K) %>% summarise(total=n()) %>% mutate(average.per.topic = total / K)
K <- promoter.summary.df %>% pull(K) %>% unique() # K for x tick labels
# total number of significant TF plots
p1 <- promoter.summary.df %>% ggplot(aes(x=K, y=total)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of Transcription Factors\n with adjusted p-value < ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +
  ggtitle(paste0("Transcription Factors Enriched in \nPromoters of the Top 100 Genes of Each Topic "))
print(p1)
p2 <- promoter.summary.df %>% ggplot(aes(x=K, y=average.per.topic)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Transcription Factors per Topic \n with adjusted p-value < ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +
  ggtitle(paste0("Transcription Factors Enriched in \nPromoter of the Top 100 Genes of Each Topic "))
print(p2)


## plot number of unique TF enriched per K for promoters
promoter.summary.unique.df <- all.promoter.ttest.df %>% subset(p.adjust < threshold & enrichment > enrichment.thr) %>% select(K, motif) %>% unique() %>% group_by(K) %>% summarise(total=n()) %>% mutate(average.per.topic = total / K)
K <- promoter.summary.unique.df %>% pull(K) %>% unique() # K for x tick labels
# total number of significant TF plots
p3 <- promoter.summary.unique.df %>% ggplot(aes(x=K, y=total)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Number of Transcription Factors with adjusted p-value < ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K)
print(p3)
p4 <- promoter.summary.unique.df %>% ggplot(aes(x=K, y=average.per.topic)) + geom_line() + geom_point() +
  xlab("K") + ylab(paste0("Average Number of Unique Transcription Factors per Topic \n with adjusted p-value < ", threshold)) + mytheme +
  scale_x_continuous("K", labels = as.character(K), breaks = K) +
  ggtitle(paste0("Unique Transcription Factors Enriched in \nPromoters of the Top 100 Genes of Each Topic "))
print(p4)

p <- ggarrange(p1 + ggtitle("") + ylab("Total TF count"), p2 + ggtitle("") + ylab("Total TF per factor"),
               p3 + ggtitle("") + ylab("Unique TF count"), p4 + ggtitle("") + ylab("Unique TF per factor"), nrow=4) +
    ggtitle(paste0("Transcription Factors Enriched in \nPromoters of the Top 100 Genes of Each Topic "))
## p <- ggarrange(p1 + ggtitle("") + ylab("Total TF count"), p2 + ggtitle("") + ylab("Total TF per factor"),
##                p3 + ggtitle("") + ylab("Unique TF count"), p4 + ggtitle("") + ylab("Unique TF per factor"), nrow=4) +
##     ggtitle(paste0("Transcription Factors Enriched in \nPromoters of the Top 100 Genes of Each Topic "))
print(p)
dev.off()


## ## cluster theta.zscore across topics
## ## old 211115
## theta.zscore.df.wide <- theta.zscore.df %>% mutate(K_Factor = paste0("K",K,"_",Factor), Gene = rownames(.)) %>% select(Gene, weight, K_Factor) %>% spread(key = "K_Factor", value = "weight")
## write.table(theta.zscore.df.wide, paste0(OUTDIR, "topic.zscore.Pearson.corr.txt"), row.names=F, quote=F, sep="\t")
## theta.zscore.df.wide.mtx <- theta.zscore.df.wide %>% `rownames<-`(.$Gene) %>% select(-Gene) %>% as.matrix()
## d <- cor(theta.zscore.df.wide.mtx, method="pearson")
## m <- as.matrix(d)

d <- cor(theta.zscore.df, method="pearson")
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
    cexCol=0.01 * ncol(mtx), cexRow=0.01 * ncol(mtx), #4.8/sqrt(nrow(mtx))
    ## cexCol=1/(ncol(mtx)^(1/3)), cexRow=1/(ncol(mtx)^(1/3)), #4.8/sqrt(nrow(mtx))
    main=title,
    ...
  )
}


pdf(file=paste0(FIGDIR,"/cluster.topic.zscore.by.Pearson.corr.pdf"), width=30, height=30)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic zscore clustering by Pearson Correlation"))
dev.off()
png(file=paste0(FIGDIR, "/cluster.topic.zscore.by.Pearson.corr.png"), width=3000, height=3000)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic zscore clustering by Pearson Correlation"))
dev.off()

## higher values of K
index <- which(theta.zscore.df %>% colnames() %>% strsplit(split="_") %>% sapply("[[",1) %>% gsub("K","",.) >= 30)
d <- cor(theta.zscore.df[index,index], method="pearson")
m <- as.matrix(d)
pdf(file=paste0(FIGDIR,"/cluster.topic.zscore.by.Pearson.corr.K30_higher.pdf"), width=30, height=30)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic zscore clustering by Pearson Correlation"))
dev.off()
png(file=paste0(FIGDIR, "/cluster.topic.zscore.by.Pearson.corr.K30_higher.png"), width=3000, height=3000)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic zscore clustering by Pearson Correlation"))
dev.off()

## set correlation < 0.5 to zero to expand the color range
m[m<0.5] <- 0
pdf(file=paste0(FIGDIR,"/cluster.topic.zscore.by.Pearson.corr.K30_higher.threshold_cor_0.5.pdf"), width=30, height=30)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic zscore clustering by Pearson Correlation"))
dev.off()
png(file=paste0(FIGDIR, "/cluster.topic.zscore.by.Pearson.corr.K30_higher.threshold_cor_0.5.png"), width=3000, height=3000)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic zscore clustering by Pearson Correlation"))
dev.off()



## topic defined by raw weights
d <- cor(theta.raw.df, method="pearson")
m <- as.matrix(d)

pdf(file=paste0(FIGDIR,"/cluster.topic.raw.by.Pearson.corr.pdf"), width=30, height=30)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic raw weight clustering by Pearson Correlation"))
dev.off()
png(file=paste0(FIGDIR, "/cluster.topic.raw.by.Pearson.corr.png"), width=3000, height=3000)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic raw weight clustering by Pearson Correlation"))
dev.off()


## higher values of K
index <- which(theta.raw.df %>% colnames() %>% strsplit(split="_") %>% sapply("[[",1) %>% gsub("K","",.) >= 30)
d <- cor(theta.raw.df[index,index], method="pearson")
m <- as.matrix(d)
pdf(file=paste0(FIGDIR,"/cluster.topic.raw.by.Pearson.corr.K30_higher.pdf"), width=30, height=30)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic raw weight clustering by Pearson Correlation"))
dev.off()
png(file=paste0(FIGDIR, "/cluster.topic.raw.by.Pearson.corr.K30_higher.png"), width=3000, height=3000)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic raw weight clustering by Pearson Correlation"))
dev.off()

## set correlation < 0.5 to zero to expand the color range
m[m<0.5] <- 0
pdf(file=paste0(FIGDIR,"/cluster.topic.raw.by.Pearson.corr.K30_higher.threshold_cor_0.5.pdf"), width=30, height=30)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic raw weight clustering by Pearson Correlation"))
dev.off()
png(file=paste0(FIGDIR, "/cluster.topic.raw.by.Pearson.corr.K30_higher.threshold_cor_0.5.png"), width=3000, height=3000)
plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, topic raw weight clustering by Pearson Correlation"))
dev.off()





## ```{r "lenient, K=14, scatter plot between number of no_IL1B edgeR DE genes and number of significant topics"}
## threshold=0.1
## ## adjusted p-value
## p <- edgeR.test.df %>% subset(K == 14 & cutoff.type == "lenient") %>% ggplot(aes(x=NoI.num.genes, y=topic.count, color=test.type)) + geom_point() + mytheme +
##   ylab(paste0("Number of Significant Topics (adjusted p-value < ", threshold, ")")) +
##   xlab(paste0("p.01_lfc1.15 number of DE genes"))
## print(p)
## ```

## ```{r "strict K=14 scatter plot between number of no_IL1B edgeR DE genes and number of significant topics"}

## p <- edgeR.test.df %>% subset(K == 14 & cutoff.type == "strict") %>% ggplot(aes(x=NoI.num.genes, y=topic.count, color=test.type)) + geom_point() + mytheme +
##   ylab(paste0("Number of Significant Topics (adjusted p-value < ", threshold, ")")) +
##   xlab(paste0("p.001_lfc1.2 number of DE genes"))
## print(p)
## ## fdr
## ```


## <br>
## Number of strict DE genes in each perturbation CDF
## Then subset to perturbation with <150 DE genes to see the distribution more closely
## <br>
## ```{r number of strice DE genes in each perturbation CDF, warning=F,  fig.show="hold", out.width="50%"}
## p <- edgeR.DE %>% subset(cutoff.type == "strict") %>% ggplot(aes(x=NoI.num.genes)) + stat_ecdf() + mytheme +
##   xlab("Number of DE genes") + ylab("Fraction of Perturbations")
## print(p)


## p <- edgeR.DE %>% subset(cutoff.type == "strict" & NoI.num.genes < 150) %>% ggplot(aes(x=NoI.num.genes)) + stat_ecdf() + mytheme +
##   xlab("Number of DE genes") + ylab("Fraction of Perturbations")
## print(p)

## ```
## <br>
## Number of lenient DE genes in each perturbation CDF
## Subset to perturbation with <500 DE genes to see the distribution more closely 
## <br>
## ```{r number of lenient DE genes in each perturbation CDF, warning=F, fig.show="hold", out.width="50%"}
## p <- edgeR.DE %>% subset(cutoff.type == "lenient") %>% ggplot(aes(x=NoI.num.genes)) + stat_ecdf() + mytheme +
##   xlab("Number of DE genes") + ylab("Fraction of Perturbations")
## print(p)

## p <- edgeR.DE %>% subset(cutoff.type == "lenient" & NoI.num.genes < 500) %>% ggplot(aes(x=NoI.num.genes)) + stat_ecdf() + mytheme +
##   xlab("Number of DE genes") + ylab("Fraction of Perturbations")
## print(p)

## ```


## <br>
## pick the cut-off at 30 DE genes for the strict p-value and logFC condition
## <br>
## ```{r pickCutOff}
## cut.off <- 30
## ```

## <br>
## Number of perturbation with large number of edgeR DE genes (p < 0.001, logFC > 1.2) that also has significant topics
## <br>
## ```{r percent strict DE + topic model matches versus K}
## total <- edgeR.DE$Gene %>% unique() %>% length() # total number of perturbations
## thresholded.DE <- edgeR.DE %>% subset(cutoff.type == "strict" & NoI.num.genes > cut.off) 

## # significant by adjusted p-value
## all.test.df.sig <- all.test.df %>% subset(adjusted.p.value < threshold)
## match.DE.test <- merge(thresholded.DE, all.test.df.sig, by="Gene")

## # how many perturbations are signficiant by EdgeR and past threshold, and at the same time has at least one DE topic?
## match.DE.test.summary <- match.DE.test %>% group_by(test.type, K) %>% summarize(num.matched.genes = n()) %>% mutate(fraction.matched.genes = num.matched.genes / total)

## # plot 
## K <- match.DE.test.summary %>% pull(K) %>% unique() # K for x tick labels
## p <- match.DE.test.summary %>% ggplot(aes(x=K, y=num.matched.genes, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Number of perturbations")) + mytheme +
##   ggtitle(paste0("Perturbations with > ", cut.off, " DE genes (p < 0.001, logFC > 1.2) \nand have at least one DE topic (adjusted p-value < ", threshold, ")")) +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## ```
## <br>
## Fraction of perturbation with large number of edgeR DE genes (p < 0.001, logFC > 1.2) that also has significant topics
## <br>
## ```{r Fraction of perturbations with > 30 DE gene and at least one DE topic}
## p <- match.DE.test.summary %>% ggplot(aes(x=K, y=fraction.matched.genes, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Fraction of perturbations")) + mytheme +
##   ggtitle(paste0("Perturbations with > ", cut.off, " DE genes (p < 0.001, logFC > 1.2) \nand have at least one DE topic (adjusted p-value < ", threshold, ")")) +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## ```


## <br> 
## Cut off at 20 DE genes from EdgeR
## <br>
## ```{r cut.off at > 20 strict category DE genes }
## cut.off <- 20
## total <- edgeR.DE$Gene %>% unique() %>% length() # total number of perturbations
## thresholded.DE <- edgeR.DE %>% subset(cutoff.type == "strict" & NoI.num.genes > cut.off) 

## # significant by adjusted p-value
## all.test.df.sig <- all.test.df %>% subset(adjusted.p.value < threshold)
## match.DE.test <- merge(thresholded.DE, all.test.df.sig, by="Gene")

## # how many perturbations are signficiant by EdgeR and past threshold, and at the same time has at least one DE topic?
## match.DE.test.summary <- match.DE.test %>% group_by(test.type, K) %>% summarize(num.matched.genes = n()) %>% mutate(fraction.matched.genes = num.matched.genes / total)
## K <- match.DE.test.summary %>% pull(K) %>% unique() # K for x tick labels
## p <- match.DE.test.summary %>% ggplot(aes(x=K, y=fraction.matched.genes, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Fraction of perturbations")) + mytheme +
##   ggtitle(paste0("Perturbations with > ", cut.off, " DE genes (p < 0.001, logFC > 1.2) \nand have at least one DE topic (adjusted p-value < ", threshold, ")")) +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## ```


## <br>
## Percent of perturbation with large number of edgeR DE genes (p < 0.01, logFC > 1.15) that also has significant topics
## <br>
## ```{r percent lenient DE + topic model matches versus K}
## cut.off <- 200
## total <- edgeR.DE$Gene %>% unique() %>% length() # total number of perturbations
## thresholded.DE <- edgeR.DE %>% subset(cutoff.type == "lenient" & NoI.num.genes > cut.off) 

## # significant by adjusted p-value
## all.test.df.sig <- all.test.df %>% subset(adjusted.p.value < threshold)
## match.DE.test <- merge(thresholded.DE, all.test.df.sig, by="Gene")

## # how many perturbations are signficiant by EdgeR and past threshold, and at the same time has at least one DE topic?
## match.DE.test.summary <- match.DE.test %>% group_by(test.type, K) %>% summarize(num.matched.genes = n()) %>% mutate(fraction.matched.genes = num.matched.genes / total)

## # plot 
## K <- match.DE.test.summary %>% pull(K) %>% unique() # K for x tick labels
## p <- match.DE.test.summary %>% ggplot(aes(x=K, y=fraction.matched.genes, color=test.type)) + geom_line() + geom_point() +
##   xlab("K") + ylab(paste0("Fraction of perturbations")) + mytheme +
##   ggtitle(paste0("Perturbations with > ", cut.off, " DE genes (p < 0.01, logFC > 1.15) \nand have at least one DE topic (adjusted p-value < ", threshold, ")")) +
##   scale_x_continuous("K", labels = as.character(K), breaks = K)
## print(p)
## ```


## <br>
## Percent of perturbation with at least one Seurat DE gene (p.adj < 0.1, FC > 1.5) that also has significant topics
## <br>
## ```{r percent lenient Seurat DE + topic model matches versus K}
## de.markers.sig.strict <- de.markers %>% subset(p_val_adj < 0.1 & abs(avg_logFC) > 0.69)
## de.markers.sig.lenient <- de.markers %>% subset(p_val_adj < 0.15 & abs(avg_logFC) > 0.4)



## ```



## ## ## FGSEA top pathways for each topic in each K
## ## ## top fgsea pathways for raw score
## ## threshold <- opt$p.adj.threshold
## ## fgsea.results.df.edited <- fgsea.results.df %>% group_by(K, database) %>% mutate(padj.over.K = p.adjust(pval))
## ## sig.fgsea.results.raw.df <- fgsea.results.df.edited %>% subset(padj.over.topics < threshold & type == "raw.score") 
## ## sig.fgsea.results.zscore.df <- fgsea.results.df.edited %>% subset(padj.over.topics < threshold & type == "z.score")
## ## toPlot.raw <- sig.fgsea.results.raw.df %>% group_by(K, topic, database) %>% arrange(padj.over.topics) %>% slice(1:10)
## ## toPlot.raw$padj.over.K <- as.numeric(toPlot.raw$padj.over.topics)
## ## toPlot.zscore <- sig.fgsea.results.zscore.df %>% group_by(K, topic, database) %>% arrange(padj.over.topics) %>% slice(1:10)

## ## # make plots for all K, topic, and database
## ## pdf(file=paste0(FIGDIR,"/Top.FGSEA.MSigDB.pathway.list_p.adj.", threshold, ".pdf"))
## ## for (k in toPlot.raw$K %>% unique()) {
## ## # for (k in c(14)) {
## ##   for (t in toPlot.raw$topic %>% unique() %>% sort()) {
## ##     for (db in toPlot.raw$database %>% unique()) {
## ##       if(nrow(toPlot.raw %>% subset(K == k & topic == t & database == db)) > 0){
## ##         p <- toPlot.raw %>% subset(K == k & topic == t & database == db) %>% arrange(padj.over.K) %>% ggplot(aes(x=pathway, y=padj.over.topics)) + geom_col() + mytheme + coord_flip() +
## ##           ggtitle(paste0("K = ", k, ", Topic ", t, ", database: ", db, "\nFGSEA from raw weights"))
## ##         print(p)
## ##       }
## ##       if(nrow(toPlot.zscore %>% subset(K == k & topic == t & database == db) > 0)){
## ##         p <- toPlot.zscore %>% subset(K == k & topic == t & database == db) %>% arrange(padj.over.K) %>% ggplot(aes(x=pathway, y=padj.over.topics)) + geom_col() + mytheme + coord_flip()  +
## ##           ggtitle(paste0("K = ", k, ", Topic ", t, ", database: ", db, "\nFGSEA from z-score (set of specific genes)"))
## ##         print(p)
## ##       }
## ##     }
## ##   }
## ## }
## ## dev.off()



## ## ## FGSEA top pathways for each topic in each K
## ## <br>
## ## Number of Enriched and Significant GSEA pathways per K
## ## <br>
## ## ```{r top fgsea pathways for raw score, fig.show="hold", out.width="50%"}

## ## ```






## ## ## 
## ## <br>
## ## <br>
## ## ```{r }
## ## ```





