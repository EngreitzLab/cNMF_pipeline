## plot PoPS analysis results
## Helen Kang
## 211109

library(conflicted)
conflict_prefer("first", "dplyr")
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("collapse", "dplyr")
conflict_prefer("combine", "dplyr")
packages <- c("optparse","dplyr", "data.table", "reshape2", "ggplot2",
              "tidyr", "textshape","readxl", "gplots", "AnnotationDbi",
              "org.Hs.eg.db", "ggrepel", "gplots")
xfun::pkg_attach(packages)
conflict_prefer("select", "dplyr")


option.list <- list(
    make_option("--project", type="character", default = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211130_test_PoPS.plots/", help="project directory"),
    make_option("--sampleName", type="character", default="2kG.library", help="project name"),
    make_option("--output", type="character", default = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211130_test_PoPS.plots/outputs/", help="output directory"),
    make_option("--figure", type="character", default = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211130_test_PoPS.plots/figures/", help="figure directory"),
    make_option("--scratch.output", type="character", default="", help="output directory for large files"),
    make_option("--k.val", type="numeric", default=60, help="the value of K in this run"),
    make_option("--PoPS_outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/", help="PoPS output directory"),
    make_option("--coefs_with_cNMF", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60.coefs", help=""),
    make_option("--marginals_with_cNMF", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60.marginals", help=""),
    make_option("--preds_with_cNMF", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60.preds", help=""),
    make_option("--coefs_without_cNMF", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6.coefs", help=""),
    make_option("--marginals_without_cNMF", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6.marginals", help=""),
    make_option("--preds_without_cNMF", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6.preds", help=""),
    make_option("--prefix", type="character", default="CAD_aug6_cNMF60", help="magma file name (before genes.raw)"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--cNMF.features", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/data/features/pops_features_raw/topic.zscore.ensembl.scaled_k_60.dt_0_2.txt", help="normalized cNMF weights, unit variance and zero mean"),
    make_option("--all.features", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211101_normalized_features/outputs/full_features_with_cNMF.RDS", help=".RDS file with all features input into PoPS"),
    make_option("--external.features.metadata", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/metadata/metadata_jul17.txt", help="annotations for each external features"),
    make_option("--combined.preds", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60.combined.preds", help="preds file with results from no_cNMF run and with_cNMF run"),
    make_option("--coefs.defining.top.topic.RDS", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60_coefs.defining.top.topic.RDS", help=""),
    make_option("--preds.importance.score.key.columns", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60_preds.importance.score.key.columns.txt", help="")
)
opt <- parse_args(OptionParser(option_list=option.list))


## debug PoPS.plots.R using scRNAseq_11AMDox_1 sample
opt$sampleName <- "scRNAseq_2kG_11AMDox_1"
opt$output <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/"
opt$figure <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/figures/all_genes/scRNAseq_2kG_11AMDox_1/K5/"
opt$scratch.output <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/211101_20sample_snakemake/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/"
opt$prefix <- "CAD_aug6_cNMF5"
opt$k.val <- 5
opt$coefs_with_cNMF <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/CAD_aug6_cNMF5.coefs"
opt$preds_with_cNMF <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/CAD_aug6_cNMF5.preds"
opt$marginals_with_cNMF <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/CAD_aug6_cNMF5.marginals"
opt$coefs_without_cNMF <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/pops/CAD_aug6.coefs"
opt$preds_without_cNMF <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/pops/CAD_aug6.preds"
opt$marginals_without_cNMF <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/CAD_aug6_cNMF5.marginals"
opt$cNMF.features <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/topic.zscore.ensembl.scaled_k_5.dt_0_2.txt"
opt$all.features <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/full_features_CAD_aug6_cNMF5.RDS"
opt$external.features.metadata <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/metadata/metadata_jul17.txt"
opt$combined.preds <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/CAD_aug6_cNMF5.combined.preds"
opt$coefs.defining.top.topic.RDS <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/211101_20sample_snakemake/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/CAD_aug6_cNMF5_coefs.defining.top.topic.RDS"
opt$preds.importance.score.key.columns <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/K5/threshold_0_2/pops/CAD_aug6_cNMF5_PoPS_preds.importance.score.key.columns.txt"


SAMPLE=opt$sampleName
OUTDIR=opt$output
SCRATCH.OUTDIR=opt$scratch.output
FIGDIR=opt$figure
PREFIX=opt$prefix
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
DENSITY.THRESHOLD.FILENAME =paste0("dt_", DENSITY.THRESHOLD)
k <- opt$k.val

check.dir <- c(OUTDIR, FIGDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))



## graphing constants
mytheme <- theme_classic() + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(hjust = 0.5, face = "bold", size=14))
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)


## load metadata
meta.data.path <- opt$external.features.metadata
metadata <- read.delim(meta.data.path, stringsAsFactors=F)

## load data
preds <- read.table(file=opt$preds_with_cNMF,header=T, stringsAsFactors=F, sep="\t")
colnames(preds) <- paste0(colnames(preds), "_with.cNMF")
colnames(preds)[1] <- "ENSGID"
preds.before <- read.table(file=paste0(opt$preds_without_cNMF), header=T, stringsAsFactors=F, sep="\t")
colnames(preds.before) <- paste0(colnames(preds.before), "_without.cNMF")
colnames(preds.before)[1] <- "ENSGID"
preds.combined <- merge(preds, preds.before, by="ENSGID")
marginals <- read.table(file=opt$marginals_with_cNMF,header=T, stringsAsFactors=F, sep="\t") %>% merge(metadata, by="X")
coefs <- read.table(file=opt$coefs_with_cNMF,header=T, stringsAsFactors=F, sep="\t")
coefs.df <- coefs[4:nrow(coefs),] %>% merge(metadata, by.x="parameter", by.y="X") %>% arrange(desc(beta)) 
coefs.df$beta <- coefs.df$beta %>% as.numeric


## map ids
x <- org.Hs.egENSEMBL 
mapped_genes <- mappedkeys(x)
xx.entrez.to.ensembl <- as.list(x[mapped_genes]) # EntrezID to Ensembl
xx.ensembl.to.entrez <- as.list(org.Hs.egENSEMBL2EG) # Ensembl to EntrezID


y <- org.Hs.egGENENAME
y_mapped_genes <- mappedkeys(y)
entrez.to.genename <- as.list(y[y_mapped_genes])
genename.to.entrez <- as.list(org.Hs.egGENENAME)


z <- org.Hs.egSYMBOL
z_mapped_genes <- mappedkeys(z)
entrez.to.symbol <- as.list(z[z_mapped_genes])
symbol.to.entrez <- as.list(org.Hs.egSYMBOL)


## preds.df
preds.df <- preds %>% mutate(EntrezID = xx.ensembl.to.entrez[preds$ENSGID %>% as.character] %>% sapply("[[",1)) %>%
    mutate(Gene.name = entrez.to.genename[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character,
           Gene = entrez.to.symbol[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character)

## load all preds
preds.combined.df <- read.delim(file=opt$combined.preds, stringsAsFactors=F)

## load cNMF features
features <- read.delim(opt$cNMF.features, stringsAsFactors=F)

## load all features
all.features <- readRDS(opt$all.features)
features <- read.delim(opt$cNMF.features, stringsAsFactors=F)

# ## load gene x feature importance score
# load(file=paste0(OUTDIR, "/coefs.marginals.feature.outer.prod.RDS"))

## load coefs.defining.top.topic.df
coefs.defining.top.topic.df <- readRDS(file=opt$coefs.defining.top.topic.RDS) # paste0(SCRATCH.OUTDIR, "/", PREFIX, "_coefs.defining.top.topic.RDS"))

## load PoPS importance score with key columns table
PoPS_preds.importance.score.key <- read.delim(file=opt$preds.importance.score.key.columns, stringsAsFactors=F) #paste0(OUTDIR, "/PoPS_preds.importance.score.key.columns.txt"), stringsAsFactors = F)


##################################################
## Plots

## plot the list of topics and their PoPS component scores for gene of interest
gene.set <- c("GOSR2", "TLNRD1", "EDN1", "NOS3", "KLF2", "ERG", "CCM2", "KRIT")
pdf(paste0(FIGDIR, "/", PREFIX, "_", DENSITY.THRESHOLD.FILENAME, "_feature_x_gene.component.importance.score.coefs.pdf"), width=4, height=6)
for (gene.here in gene.set) {
    toPlot <- coefs.defining.top.topic.df %>% subset(grepl(paste0("^",gene.here,"$"), Gene)) %>% select(topic, gene.feature_x_beta)
    p <- toPlot %>% ggplot(aes(x=reorder(topic, gene.feature_x_beta), y=gene.feature_x_beta)) + geom_col(fill="#38b4f7") + theme_minimal() +
        coord_flip() + xlab("Feature (Topic)") + ylab("Feature x Gene\nImportance Score") + ggtitle(paste0(gene.here)) + mytheme
    print(p)

    toPlot <- PoPS_preds.importance.score.key %>% subset(grepl(paste0("^",gene.here,"$"), Gene)) %>% select(Long_Name, gene.feature_x_beta) %>% slice(1:15)
    p <- toPlot %>% ggplot(aes(x=reorder(Long_Name, gene.feature_x_beta), y=gene.feature_x_beta)) + geom_col(fill="#38b4f7") + theme_minimal() +
        coord_flip() + xlab("Features") + ylab("Feature x Gene\nImportance Score") + ggtitle(paste0(gene.here)) + mytheme
}
dev.off()

 
pdf(paste0(FIGDIR, "/", PREFIX, "_", DENSITY.THRESHOLD.FILENAME, "_Topic.coef.beta.pdf"), width=4, height=6) ## double check figures
toPlot <- coefs.df %>% subset(grepl("zscore",parameter)) %>% arrange(desc(beta))
p <- toPlot %>% ggplot(aes(x=reorder(parameter, beta), y=beta)) + geom_col(fill="#38b4f7") + theme_minimal() +
    coord_flip() + xlab("Features (topic)") + ylab("Beta Score") + ggtitle(paste0("K = ", k, " Topics")) + mytheme
print(p)
dev.off()
pdf(paste0(FIGDIR, "/all.coef.beta.pdf"))
toPlot <- coefs.df %>% arrange(desc(beta)) %>% slice(1:15)
p <- toPlot %>% ggplot(aes(x=reorder(Long_Name, beta), y=beta)) + geom_col(fill="#38b4f7") + theme_minimal() +
    coord_flip() + xlab("Features") + ylab("Beta Score") + ggtitle(paste0("K = ", k, " Topics")) + mytheme
print(p)
dev.off()
pdf(paste0(FIGDIR, "/all.marginals.beta.pdf"))
toPlot <- marginals %>% arrange(desc(beta)) %>% slice(1:15)
p <- toPlot %>% ggplot(aes(x=reorder(Long_Name, beta), y=beta)) + geom_col(fill="#38b4f7") + theme_minimal() +
    coord_flip() + xlab("Features") + ylab("Beta Score") + ggtitle(paste0("K = ", k, " Topics, Ranked by marginals beta score")) + mytheme
print(p)
dev.off()



## figure
pdf(paste0(FIGDIR, "/", PREFIX, "_", DENSITY.THRESHOLD.FILENAME, "_PoPS_score_list.pdf"), width=4, height=6)
top.PoPS.genes <- preds.df %>% slice(1:20)
toPlot <- data.frame(Gene= top.PoPS.genes %>% pull(Gene),
                     Score=top.PoPS.genes %>% pull(PoPS_Score_with.cNMF))
p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(fill="#38b4f7") + theme_minimal()
p <- p + coord_flip() + xlab("Top 50 Genes") + ylab("PoPS Score") + ggtitle(paste0(SAMPLE, ", ", PREFIX)) + mytheme
print(p)
dev.off()


## PoPS before vs after score
pdf(paste0(FIGDIR, "/", PREFIX, "_", DENSITY.THRESHOLD.FILENAME, "_before.vs.after.cNMF.pdf"))
labels <- preds.combined.df %>% subset((PoPS_Score_with.cNMF > (2 * PoPS_Score_without.cNMF)) & PoPS_Score_with.cNMF > 1)
p <- preds.combined.df %>% ggplot(aes(x=PoPS_Score_without.cNMF, y=PoPS_Score_with.cNMF)) + geom_point(size=0.5) + mytheme + 
xlab("PoPS Score (without cNMF features)") + ylab("PoPS Score(with cNMF features)") + geom_abline(slope=1, color="red") + geom_text_repel(data=labels, box.padding = 0.5,
                                                                                                     max.overlaps=30,
                                                                                                     aes(label=Gene), size=4,
                                                                                                     color="blue")
print(p)
dev.off()



# ## ## ## ## ## ##
# ## consider deleting:
# ## PoPS after - before for score, zscore, and ranking
# zscore <- function(array) (array - mean(array)) / sd(array)
# preds.all.df <- preds.all.df %>%
#     mutate(PoPS_Score_after_zscore = zscore(PoPS_Score_after),
#            PoPS_Score_before_zscore = zscore(PoPS_Score_before),
#            PoPS_Score_after.diff.before.zscore = PoPS_Score_after_zscore - PoPS_Score_before_zscore, ## check difference in zscore between with cNMF and no cNMF
#            PoPS_Score_after.diff.before = PoPS_Score_after - PoPS_Score_before) %>%
#     arrange(desc(PoPS_Score_before)) %>%
#     mutate(PoPS_Score_before_ranking = 1:n()) %>%
#     arrange(desc(PoPS_Score_after)) %>%
#     mutate(PoPS_Score_after_ranking = 1:n(),
#            PoPS_Score_before.diff.after.ranking = PoPS_Score_before_ranking - PoPS_Score_after_ranking) %>% ## check difference in ranking (before - after to check for "growth" in ranking)
#     arrange(desc(PoPS_Score_before.diff.after.ranking))


# ## plot after.diff.before results
# pdf(paste0(FIGDIR, "after.diff.before.cNMF.pdf"))
# top.PoPS.genes <- preds.all.df %>% slice(1:20) ## sorted by the change in ranking
# toPlot <- data.frame(Gene= top.PoPS.genes %>% pull(Gene),
#                      Score=top.PoPS.genes %>% pull(PoPS_Score_before.diff.after.ranking))
# p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(fill="#38b4f7") + theme_minimal()
# p <- p + coord_flip() + xlab("Top 20 Genes") + ylab("PoPS Ranking Increase") + ggtitle(paste0(SAMPLE)) + mytheme
# print(p)

# top.PoPS.genes <- preds.all.df %>% arrange(desc(PoPS_Score_after.diff.before.zscore)) %>% slice(1:20) ## sorted by the change in ranking
# toPlot <- data.frame(Gene= top.PoPS.genes %>% pull(Gene),
#                      Score=top.PoPS.genes %>% pull(PoPS_Score_after.diff.before.zscore))
# p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(fill="#38b4f7") + theme_minimal()
# p <- p + coord_flip() + xlab("Top 20 Genes") + ylab("PoPS z-score difference (with cNMF features - without)") + ggtitle(paste0(SAMPLE)) + mytheme
# print(p)

# top.PoPS.genes <- preds.all.df %>% arrange(desc(PoPS_Score_after.diff.before)) %>% slice(1:20) ## sorted by the change in ranking
# toPlot <- data.frame(Gene= top.PoPS.genes %>% pull(Gene),
#                      Score=top.PoPS.genes %>% pull(PoPS_Score_after.diff.before))
# p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(fill="#38b4f7") + theme_minimal()
# p <- p + coord_flip() + xlab("Top 20 Genes") + ylab("PoPS Score difference (with cNMF features - without)") + ggtitle(paste0(SAMPLE)) + mytheme
# print(p)

# ## waterfall plot for the difference ranking
# mytheme <- theme_classic() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, face = "bold"))
# labels <- preds.all.df %>% arrange(desc(PoPS_Score_after.diff.before)) %>% slice(1:20)

# p <- preds.all.df %>% ggplot(aes(x=reorder(Gene, PoPS_Score_after.diff.before), y=PoPS_Score_after.diff.before)) + geom_point(size=0.75) + mytheme +
#     theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + scale_color_manual(values = c("#E0E0E0", "#38b4f7")) +
#     xlab("Genes") + ylab(paste0("PoPS Score Difference (with cNMF features - without)")) +
#     geom_text_repel(data=labels, box.padding = 0.35,
#                     aes(label=Gene), size=5,
#                     color="black") +
#     theme(legend.position = "none") # legend at the bottom?
# print(p)
# dev.off()

