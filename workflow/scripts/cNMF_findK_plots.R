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
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggpubr)) ## annotate arranged figure

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

## ## for all genes 210707 folder
## ## ## all genes directories (for sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/2kG.library/acrossK/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/acrossK/"
## opt$sampleName <- "2kG.library"

## ## for all genes (in sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/figures/all_genes/Perturb_2kG_dup4/acrossK/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/Perturb_2kG_dup4/acrossK/"
## opt$aggregated.data <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/Perturb_2kG_dup4/acrossK//aggregated.outputs.findK.RData"
## opt$sampleName <- "Perturb_2kG_dup4"

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

# ## for testing findK_plots for control only cells
# opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/figures/all_genes/2kG.library.ctrl.only/acrossK/"
# opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/analysis/all_genes/2kG.library.ctrl.only/acrossK/"
# opt$sampleName <- "2kG.library.ctrl.only"
# opt$aggregated.data <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/analysis/all_genes/2kG.library.ctrl.only/acrossK/aggregated.outputs.findK.RData"


## Directories and Constants
SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
threshold <- opt$p.adj.threshold
DATADIR=opt$datadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
FIGDIR=opt$figdir
check.dir <- c(OUTDIR, FIGDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))


## load("/Volumes/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2105_findK/analysis/no_IL1B/aggregated.outputs.findK.RData")
## load("/Volumes/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/differential_expression/210526_SeuratDE/outputs/no_IL1B/de.markers.RData")
mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))
mytheme <- theme_classic() + theme(axis.text = element_text(size = 5),
                                   axis.title = element_text(size = 6),
                                   plot.title = element_text(hjust = 0.5, face = "bold", size=7),
                                   axis.line = element_line(color = "black", size = 0.25),
                                   axis.ticks = element_line(color = "black", size = 0.25))
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)
palette.log2FC = colorRampPalette(c("blue", "white", "red"))(n = 100)
palette.pos = colorRampPalette(c("white", "red"))(n = 100)
adjusted.p.thr <- 0.05
fdr.thr <- 0.05

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

##########################################################################################
## GSEA plots
##################################################
## labels and parameters for the for loops
ranking.types <- c("zscore", "raw")
GSEA.types <- c("GOEnrichment", "ByWeightGSEA", "GSEA")
GSEA.type.labels <- c("GO Pathway Enrichment\n(Top 300 Program Gene by Hypergeometric Test)", "Gene Set Enrichment\n(All Genes)", "Gene Set Enrichment\n(Top 300 Program Gene by Hypergeomteric Test)")
##################################################
## function to output plot in the same format
plotGSEA <- function(toplot, nPathwayMetric, nPathwayMetricLabel) {
    p <- toplot %>% ggplot(aes(x=K, y=get(nPathwayMetric))) + geom_line(size=0.5) + geom_point(size=0.5) + mytheme +
        xlab("K") + ylab(nPathwayMetricLabel) #+ #ggtitle(paste0(GSEA.type.label)) +
        ## scale_x_continuous("K", labels = as.character(K), breaks = K)
    print(p)
    return(p)
}

##################################################
for (GSEA.type.i in 1:length(GSEA.types)) {
    GSEA.type <- GSEA.types[GSEA.type.i]
    GSEA.type.label <- GSEA.type.labels[GSEA.type.i]

##################################################
    ## process the GSEA data here
    tmp.df <- get(paste0("clusterProfiler.", GSEA.type, ".df")) %>%
        subset(p.adjust < fdr.thr) %>%
        group_by(K, type) %>%
        mutate(nPathways = n()) %>%
        select(K, type, nPathways, ID, Description) %>%
        unique %>%
        mutate(nUniquePathways = n(),
               normalizedNPathways = nPathways / K,
               normalizedNUniquePathways = nUniquePathways / K) %>%
        as.data.frame
    summary.df <- tmp.df %>%
        select(-ID, -Description) %>%
        unique %>%
        `rownames<-`(paste0(.$type, "_K", .$K))
    K <- summary.df %>% pull(K) %>% unique() # K for x tick labels

##################################################
    for (ranking.type in ranking.types) {
        toplot <- summary.df %>%
            subset(type == ranking.type) ## subset to selected program gene ranking type (zscore or raw)

        nPathwayMetrics <- toplot %>% select(-K, -type) %>% colnames
        nPathwayMetricLabels <- c("# Total Pathways", "# Unique Pathways", "Average Pathways\nper Program", "Aveage Unique Pathways\nper Program")

        plotFilename <- paste0(FIGDIR, GSEA.type, "_", ranking.type)
        pdf(paste0(plotFilename, ".pdf"), width=3, height=3)
        plot.list <- list() ## create a list to store all metrics
        for (i in 1:length(nPathwayMetrics)) { ## make a line plot 
            nPathwayMetric <- nPathwayMetrics[i]
            nPathwayMetricLabel <- nPathwayMetricLabels[i]
            plot.list[[nPathwayMetric]] <- plotGSEA(toplot, nPathwayMetric, nPathwayMetricLabel)
        }
        p <- plot_grid(plotlist = plot.list, nrow=length(plot.list), align="v", axis="lr") 
        p <- annotate_figure(p, top = text_grob(paste0(GSEA.type.label, "\n", ranking.type), size=8))
        print(p)
        eval(parse(text = paste0("p.", GSEA.type, ".", ranking.type, " <- p")))
        dev.off()
    }
}

## combine all plots in one panel
p.all.GSEA <- plot_grid(p.GOEnrichment.zscore, p.GOEnrichment.raw,
               p.GSEA.zscore, p.GSEA.raw,
               p.ByWeightGSEA.zscore, p.ByWeightGSEA.raw,
               nrow = 3, ncol = 2, align="hv", axis="tblr")
plotFilename <- paste0(FIGDIR, "/All_GSEA")
pdf(paste0(plotFilename, ".pdf"), width=6, height=10)
print(p.all.GSEA)
ggsave(paste0(plotFilename, ".eps"))
dev.off()

##################################################
## End of GSEA plots
##################################################


## Notes:
## why is there a sharp change from K=19 to K=21?


## Fraction of topics that has at least one significant (metric) versus K
num.program.genes <- 300
enrichment.thr <- 1
## plot number of TF enriched per K for {enhancers, promoters}
for (ep.type in c("promoter", "enhancer")) {
    ep.type.label <- ifelse(ep.type == "promoter", "Promoter", "Enhancer")
    pdf(file=paste0(FIGDIR, "/TF.motif.enrichment.", ep.type, ".fdr.thr", as.character(fdr.thr), ".pdf"), width=3, height=3)
    ## promoters
    toplot <- get(paste0("all.", ep.type, ".ttest.df")) %>%
        mutate(significant = two.sided.p.adjust < fdr.thr & enrichment > enrichment.thr) %>%
        group_by(K) %>%
        summarise(total=significant %>% as.numeric %>% sum) %>%
        mutate(average.per.topic = total / K) %>%
        as.data.frame
    ## K <- toplot %>% pull(K) %>% unique() # K for x tick labels
    title <- paste0("Transcription Factors Enriched in \n", ep.type.label, " of Program Genes ")
    ## total number of significant TF plots
    p1 <- toplot %>% ggplot(aes(x=K, y=total)) + geom_line(size=0.5) + geom_point(size=0.5) +
    ## p1 <- toplot %>% ggplot(aes(x=K, y=total)) + geom_col(color="gray35") +
        xlab("K") + ylab(paste0("# Transcription Factors\n (FDR < ", fdr.thr, ")")) + mytheme +
        ## scale_x_continuous("K", labels = as.character(K), breaks = K) +
        ggtitle(title)
    print(p1)
    p2 <- toplot %>% ggplot(aes(x=K, y=average.per.topic)) + geom_line(size=0.5) + geom_point(size=0.5) +
    ## p2 <- toplot %>% ggplot(aes(x=K, y=average.per.topic)) + geom_col(width=0.5, color="gray35") +
        xlab("K") + ylab(paste0("Average # Transcription Factors per Program\n(FDR < ", fdr.thr, ")")) + mytheme +
        ## scale_x_continuous("K", labels = as.character(K), breaks = K) +
        ggtitle(paste0("Transcription Factors Enriched in \nPromoter of the Top ", num.program.genes, " Genes of Each Program"))
    print(p2)

    ## plot number of unique TF enriched per K for promoters
    toplot.unique <- get(paste0("all.", ep.type, ".ttest.df")) %>%
        mutate(significant = two.sided.p.adjust < fdr.thr & enrichment > enrichment.thr) %>%
        select(K, motif, significant) %>%
        unique() %>%
        group_by(K) %>%
        summarise(total=significant %>% as.numeric %>% sum) %>%
        mutate(average.per.topic = total / K) %>%
        as.data.frame
    K <- toplot.unique %>% pull(K) %>% unique() # K for x tick labels
    ## total number of significant TF plots
    p3 <- toplot.unique %>% ggplot(aes(x=K, y=total)) + geom_line(size=0.5) + geom_point(size=0.5) +
        xlab("K") + ylab(paste0("Number of Transcription Factors with adjusted p-value < ", fdr.thr)) + mytheme #+
        ## scale_x_continuous("K", labels = as.character(K), breaks = K)
    print(p3)
    p4 <- toplot.unique %>% ggplot(aes(x=K, y=average.per.topic)) + geom_line(size=0.5) + geom_point(size=0.5) +
        xlab("K") + ylab(paste0("Average Number of Unique Transcription Factors per Program \n with adjusted p-value < ", fdr.thr)) + mytheme #+
        ## scale_x_continuous("K", labels = as.character(K), breaks = K) +
        ggtitle(paste0("Unique Transcription Factors Enriched in \nPromoters of the Top 100 Genes of Each Program"))
    print(p4)

    p <- plot_grid(p1 + ggtitle("") + ylab("# TFs"), p2 + ggtitle("") + ylab("# TF per\nProgram"),
                   p3 + ggtitle("") + ylab("# Unique TFs"), p4 + ggtitle("") + ylab("# Unique TF\nper Program"), nrow=4, align = "hv", axis="tblr")
    ## title = ep.type.label
    eval(parse(text = paste0("p.", ep.type, " <- annotate_figure(p, top = text_grob(label=title, face='bold', size=8))")))
    print(get(paste0("p.", ep.type)))    
    dev.off()
}

## put together all motif enrichment plot in one panel
p.all.TFMotifEnrichment <- plot_grid(p.promoter, p.enhancer, nrow=1)
plotFilename <- paste0(FIGDIR, "/All_TFMotifEnrichment")
pdf(paste0(plotFilename, ".pdf"), width=4, height=4)
print(p.all.TFMotifEnrichment)
ggsave(paste0(plotFilename, ".eps"))
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





