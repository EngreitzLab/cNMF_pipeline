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
    make_option("--output", type="character", default = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211130_test_PoPS.plots/outputs/", help="output directory"),
    make_option("--figure", type="character", default = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211130_test_PoPS.plots/figures/", help="figure directory"),
    make_option("--scratch.output", type="character", default="", help="output directory for large files"),
    make_option("--PoPS_outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/", help="PoPS output directory"),
    make_option("--coefs", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60.coefs", help=""),
    make_option("--marginals", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60.marginals", help=""),
    make_option("--preds", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60.preds", help=""),
    make_option("--prefix", type="character", default="CAD_aug6_cNMF60", help="magma file name (before genes.raw)"),
    make_option("--sampleName", type="character", default="2kG.library", help="Sample name for this run to output"),
    make_option("--all.features", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211101_normalized_features/outputs/full_features_with_cNMF.RDS", help=".RDS file with all features input into PoPS"),
    make_option("--external.features.metadata", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/metadata/metadata_jul17.txt", help="annotations for each external features"),
    make_option("--coefs.defining.top.topic.RDS", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60_coefs.defining.top.topic.RDS", help=""),
    make_option("--preds.importance.score.key.columns", type="character", defualt="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60_preds.importance.score.key.columns.txt", help="")
)
opt <- parse_args(OptionParser(option_list=option.list))


OUTDIR=opt$output
SCRATCH.OUTDIR=opt$scratch.output
FIGDIR=opt$figure
PREFIX=opt$prefix
SAMPLE=opt$sampleName

check.dir <- c(OUTDIR, FIGDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))



## graphing constants
mytheme <- theme_classic() + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(hjust = 0.5, face = "bold", size=14))
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)


## load metadata
meta.data.path <- opt$external.features.metadata
metadata <- read.delim(meta.data.path, stringsAsFactors=F)

## load data
print("loading PoPS results")
preds <- read.table(file=opt$preds,header=T, stringsAsFactors=F, sep="\t")
colnames(preds)[1] <- "ENSGID"
marginals <- read.table(file=opt$marginals,header=T, stringsAsFactors=F, sep="\t")
coefs <- read.table(file=opt$coefs,header=T, stringsAsFactors=F, sep="\t")
coefs.df <- coefs[4:nrow(coefs),] %>% arrange(desc(beta))
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


## function for adding name to df
add_gene_name_to_df <- function(df) {
    out <- df %>% mutate(EntrezID = xx.ensembl.to.entrez[df$ENSGID %>% as.character] %>% sapply("[[",1)) %>%
    mutate(Gene.name = entrez.to.genename[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character,
           Gene = entrez.to.symbol[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character)
    return(out)
}


## preds.df
preds.df <- preds %>% add_gene_name_to_df


## load all features
all.features <- loadRDS(opt$all.features)

# ## load gene x feature importance score
# load(file=paste0(OUTDIR, "/coefs.marginals.feature.outer.prod.RDS"))

## load coefs.defining.top.topic.df
coefs.defining.top.topic.df <- readRDS(file=opt$coefs.defining.top.topic.RDS) # paste0(SCRATCH.OUTDIR, "/", PREFIX, "_coefs.defining.top.topic.RDS"))

## load PoPS importance score with key columns table
PoPS_preds.importance.score.key <- read.delim(file=opt$preds.importance.score.key.columns, stringsAsFactors=F) #paste0(OUTDIR, "/PoPS_preds.importance.score.key.columns.txt"), stringsAsFactors = F)


##################################################
## Plots

# ## plot the list of topics and their PoPS component scores for gene of interest
# gene.set <- c("GOSR2", "TLNRD1", "EDN1", "NOS3", "KLF2", "ERG", "CCM2", "KRIT")
# pdf(paste0(FIGDIR, "/", PREFIX, "_", DENSITY.THRESHOLD.FILENAME, "_feature_x_gene.component.importance.score.coefs.pdf"), width=4, height=6)
# for (gene.here in gene.set) {
#     toPlot <- coefs.defining.top.topic.df %>% subset(grepl(paste0("^",gene.here,"$"), Gene)) %>% select(topic, gene.feature_x_beta)
#     p <- toPlot %>% ggplot(aes(x=reorder(topic, gene.feature_x_beta), y=gene.feature_x_beta)) + geom_col(fill="#38b4f7") + theme_minimal() +
#         coord_flip() + xlab("Feature (Topic)") + ylab("Feature x Gene\nImportance Score") + ggtitle(paste0(gene.here)) + mytheme
#     print(p)

#     toPlot <- PoPS_preds.importance.score.key %>% subset(grepl(paste0("^",gene.here,"$"), Gene)) %>% select(Long_Name, gene.feature_x_beta) %>% slice(1:15)
#     p <- toPlot %>% ggplot(aes(x=reorder(Long_Name, gene.feature_x_beta), y=gene.feature_x_beta)) + geom_col(fill="#38b4f7") + theme_minimal() +
#         coord_flip() + xlab("Features") + ylab("Feature x Gene\nImportance Score") + ggtitle(paste0(gene.here)) + mytheme
# }
# dev.off()

 
## top selected features
pdf(paste0(FIGDIR, "/", PREFIX, "_", SAMPLE, "_all.coef.beta.pdf"))
toPlot <- coefs.df %>% arrange(desc(beta)) %>% slice(1:15)
p <- toPlot %>% ggplot(aes(x=reorder(Long_Name, beta), y=beta)) + geom_col(fill="#38b4f7") + theme_minimal() +
    coord_flip() + xlab("Features") + ylab("Beta Score") + ggtitle(paste0("K = ", k, " Topics")) + mytheme
print(p)
dev.off()
 

## figure
pdf(paste0(FIGDIR, "/", PREFIX, "_", SAMPLE, "_PoPS_score_list.pdf"), width=4, height=6)
top.PoPS.genes <- preds.df %>% slice(1:20)
toPlot <- data.frame(Gene= top.PoPS.genes %>% pull(Gene),
                     Score=top.PoPS.genes %>% pull(PoPS_Score))
p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score) ) + geom_col(fill="#38b4f7") + theme_minimal()
p <- p + coord_flip() + xlab("Top 50 Genes") + ylab("PoPS Score") + ggtitle(paste0(SAMPLE)) + mytheme
print(p)
dev.off()


# ## PoPS before vs after score
# pdf(paste0(FIGDIR, "/", PREFIX, "_", DENSITY.THRESHOLD.FILENAME, "_before.vs.after.cNMF.pdf"))
# labels <- preds.combined.df %>% subset(PoPS_Score_after > 2 * PoPS_Score_before & PoPS_Score_after > 1)
# p <- preds.combined.df %>% ggplot(aes(x=PoPS_Score_before, y=PoPS_Score_after)) + geom_point(size=0.5) + mytheme + 
# xlab("PoPS Score (without cNMF features)") + ylab("PoPS Score(with cNMF features)") + geom_abline(slope=1, color="red") + geom_text_repel(data=labels, box.padding = 0.5,
#                                                                                                      max.overlaps=30,
#                                                                                                      aes(label=Gene), size=4,
#                                                                                                      color="blue")
# print(p)
# dev.off()



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

