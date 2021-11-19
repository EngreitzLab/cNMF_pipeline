## Helen Kang
## Script to prepare PoPS data for generating figures
## 211118 adapted from PoPS.data.processing.211115.R

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
    make_option("--project", type="character", default = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/", help="project directory"),
    make_option("--output", type="character", default = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/", help="output directory"),
    make_option("--figure", type="character", default = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/figures/", help="figure directory"),
    make_option("--coefs_with_cNMF", type="character", default="", help=""),
    make_option("--marginals_with_cNMF", type="character", default="", help=""),
    make_option("--preds_with_cNMF", type="character", default="", help=""),
    make_option("--coefs_without_cNMF", type="character", default="", help=""),
    make_option("--marginals_without_cNMF", type="character", default="", help=""),
    make_option("--preds_without_cNMF", type="character", default="", help=""),
    make_option("--prefix", type="character", default="CAD_aug6_cNMF60", help="magma file name (before genes.raw)"),
    make_option("--cNMF.features", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/data/features/pops_features_raw/topic.zscore.ensembl.scaled_k_60.dt_0_2.txt", help="normalized cNMF weights, unit variance and zero mean"),
    make_option("--all.features", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211101_normalized_features/outputs/full_features_with_cNMF.RDS", help=".RDS file with all features input into PoPS"),
    make_option("--recompute", type="logical", default=F, help="T for rerunning the entire script, F for only outputting the missing data")
)
opt <- parse_args(OptionParser(option_list=option.list))

SAMPLE=opt$prefix
OUTDIR=opt$output
FIGDIR=opt$figure
PREFIX=opt$prefix

check.dir <- c(OUTDIR, FIGDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))




## graphing constants
mytheme <- theme_classic() + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(hjust = 0.5, face = "bold", size=14))
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)


## load metadata
meta.data.path <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/metadata/metadata_jul17.txt"
metadata <- read.delim(meta.data.path, stringsAsFactors=F)


## load data
preds <- read.table(file=opt$preds_with_cNMF,header=T, stringsAsFactors=F, sep="\t")
colnames(preds) <- paste0(colnames(preds), "_with.cNMF")
colnames(preds)[1] <- "ENSGID"
preds.before <- read.table(file=paste0(opt$preds_without_cNMF), header=T, stringsAsFactors=F, sep="\t")
colnames(preds.before) <- paste0(colnames(preds.before), "_without.cNMF")
colnames(preds.before)[1] <- "ENSGID"
preds.combined <- merge(preds, preds.before, by="ENSGID")
marginals <- read.table(file=opt$marginals_with_cNMF,header=T, stringsAsFactors=F, sep="\t")
coefs <- read.table(file=opt$coefs_with_cNMF,header=T, stringsAsFactors=F, sep="\t")
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

## add Gene name to preds.df
preds.df <- preds %>% mutate(EntrezID = xx.ensembl.to.entrez[preds$ENSGID %>% as.character] %>% sapply("[[",1)) %>%
    mutate(Gene.name = entrez.to.genename[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character,
           Gene = entrez.to.symbol[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character)
file.name <- paste0(OUTDIR, "/", SAMPLE, ".combined.preds")
if( !file.exists(file.name) | opt$recompute ) {
preds.combined.df <- preds.combined %>% mutate(EntrezID = xx.ensembl.to.entrez[preds.combined$ENSGID %>% as.character] %>% sapply("[[",1)) %>%
    mutate(Gene.name = entrez.to.genename[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character,
           Gene = entrez.to.symbol[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character)
write.table(preds.combined.df %>% apply(2, as.character), file=file.name, row.names=F, quote=F, sep="\t")
} else {
    preds.combined.df <- read.delim(file.name, stringsAsFactors = F)
}



file.name <- paste0(OUTDIR, "/", PREFIX, "_coefs.marginals.feature.outer.prod.RDS")
if( !file.exists(file.name) | opt$recompute ) {

## load cNMF features
## features <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/data/features/pops_features_raw/topic.zscore.ensembl.scaled_k_60.dt_0_2.txt", stringsAsFactors=F)
features <- read.delim(opt$cNMF.features, stringsAsFactors=F)

## load all features
all.features.cNMF <- readRDS(opt$all.features)


## sort features and subset coefs, take a product ((gene x features) x (features x beta scalar))
coefs.cnmf <- coefs %>% subset(grepl("zscore", parameter))
coefs.cnmf.mtx <- coefs.cnmf %>% `rownames<-`(.$parameter) %>% select(-parameter)
coefs.cnmf.mtx$beta <- coefs.cnmf.mtx$beta %>% as.numeric
coefs.cnmf.mtx <- coefs.cnmf.mtx %>% as.matrix

## subset features
features.names.tokeep <- coefs.cnmf.mtx %>% rownames
features.tokeep <- features %>% `rownames<-`(.$ENSGID) %>% select(all_of(features.names.tokeep)) %>% as.matrix

## also test marginals
features.names.tokeep.df <- data.frame(X=features.names.tokeep)
marginals.tokeep <- inner_join(features.names.tokeep.df, marginals, by="X") %>% `rownames<-`(.$X) %>% select(beta) %>% as.matrix


coefs.mtx <- coefs.df %>% `rownames<-`(.$parameter) %>% select(-parameter)
coefs.mtx$beta <- coefs.mtx$beta %>% as.numeric
coefs.mtx <- coefs.mtx %>% as.matrix
all.features.cNMF.names.tokeep <- coefs.mtx %>% rownames
all.features.cNMF.tokeep <- all.features.cNMF %>% `rownames<-`(.$ENSGID) %>% select(all_of(all.features.cNMF.names.tokeep)) %>% as.matrix
all.features.cNMF.tokeep[all.features.cNMF.tokeep=="True"] <- 1
all.features.cNMF.tokeep[all.features.cNMF.tokeep=="False"] <- 0
storage.mode(all.features.cNMF.tokeep) <- "numeric"
all.features.cNMF.names.tokeep.df <- data.frame(X=all.features.cNMF.names.tokeep) ## to subset marginals
all.marginals.cNMF.tokeep <- inner_join(all.features.cNMF.names.tokeep.df, marginals, by="X") %>% `rownames<-`(.$X) %>% select(beta) %>% as.matrix

    saveRDS(all.features.cNMF.names.tokeep, file=paste0(OUTDIR, "/all.features.cNMF.keep.prioritized.mtx.RDS"))

## ## check if PoPS_Score = coefs * features
## ## multiply features with beta
## PoPS_Score.coefs.manual <- features.tokeep %*% coefs.cnmf.mtx %>% `colnames<-`("PoPS_Score.coefs.manual")
## PoPS_Score.marginals.manual <- features.tokeep %*% marginals.tokeep %>% `colnames<-`("PoPS_Score.marginals.manual")
## PoPS_Score.coefs.all <- all.features.cNMF.tokeep %*% coefs.mtx %>% `colnames<-`("PoPS_Score_all.coefs")
## PoPS_Score.marginals.all <- all.features.cNMF.tokeep %*% marginals.tokeep %>% `colnames<-`("PoPS_Score.all.marginals")
## ## PoPS_Score.marginals.manual.all <- features %>% `rownames<-`(.$ENSGID) %>% select(-ENSGID) %*% marginals


## get outer products of features and marginals
PoPS_Score.coefs.manual.outer.ENSG <- sweep(features.tokeep, 2, (coefs.cnmf.mtx %>% t), `*`) ### store this matrix
PoPS_Score.marginals.manual.outer.ENSG <- sweep(features.tokeep, 2, (marginals.tokeep %>% t), `*`) ### save this matrix
PoPS_Score.coefs.all.outer.ENSG <- sweep(all.features.cNMF.tokeep, 2, (coefs.mtx %>% t), `*`)
PoPS_Score.marginals.all.outer.ENSG <- sweep(all.features.cNMF.tokeep, 2, (all.marginals.cNMF.tokeep %>% t), `*`)


## function to convert the outer product ENSG name to Gene name
add_gene_name <- function(df) {
    out <- df %>% as.data.frame %>% mutate(EntrezID = xx.ensembl.to.entrez[df %>% rownames %>% as.character] %>% sapply("[[",1)) %>%
    mutate(Gene.name = entrez.to.genename[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character,
           Gene = entrez.to.symbol[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character)
    return(out)
}

PoPS_Score.coefs.manual.outer <- add_gene_name(PoPS_Score.coefs.manual.outer.ENSG)
PoPS_Score.marginals.manual.outer <- add_gene_name(PoPS_Score.marginals.manual.outer.ENSG)
PoPS_Score.coefs.all.outer <- add_gene_name(PoPS_Score.coefs.all.outer.ENSG)
PoPS_Score.marginals.all.outer <- add_gene_name(PoPS_Score.marginals.all.outer.ENSG)

## save the results
write.table(PoPS_Score.coefs.manual.outer %>% apply(2, as.character), file=paste0(OUTDIR, "/coefs.feature.outer.prod.txt"), quote=F, sep="\t")
write.table(PoPS_Score.marginals.manual.outer %>% apply(2, as.character), file=paste0(OUTDIR, "/marginals.feature.outer.prod.txt"), quote=F, sep="\t")
write.table(PoPS_Score.coefs.all.outer %>% apply(2, as.character), file=paste0(OUTDIR, "/coefs.all.feature.outer.prod.txt"), quote=F, sep="\t")
write.table(PoPS_Score.marginals.all.outer %>% apply(2, as.character), file=paste0(OUTDIR, "/marginals.all.feature.outer.prod.txt"), quote=F, sep="\t")



## find top topic that define PoPS for each gene
## function for sorting the feature x gene importance value
sort_feature_x_gene_importance <- function(df) {
    out <- df %>% mutate(ENSGID=rownames(.)) %>% melt(id.vars = c("Gene.name", "EntrezID", "Gene", "ENSGID"), variable.name="topic", value.name="gene.feature_x_beta") %>% group_by(Gene) %>% arrange(desc(gene.feature_x_beta)) %>% as.data.frame
    return(out)
}

coefs.defining.top.topic.df <- PoPS_Score.coefs.manual.outer %>% sort_feature_x_gene_importance
marginals.defining.top.topic.df <- PoPS_Score.marginals.manual.outer %>% sort_feature_x_gene_importance
all.coefs.defining.top.topic.df <- PoPS_Score.coefs.all.outer %>% sort_feature_x_gene_importance
all.marginals.defining.top.topic.df <- PoPS_Score.marginals.all.outer %>% sort_feature_x_gene_importance
save(PoPS_Score.coefs.manual.outer, PoPS_Score.marginals.manual.outer, PoPS_Score.coefs.all.outer, PoPS_Score.marginals.all.outer,
     coefs.defining.top.topic.df, marginals.defining.top.topic.df, all.coefs.defining.top.topic.df, all.marginals.defining.top.topic.df,
     file=file.name)
} else {
    load(file.name)
}


## output a table, one row per gene, with columns for gene symbol, PoPS score without cNMF, PoPS score with cNMF, top features from any source important for that gene, top topic features important for that gene
coefs.defining.top.topic.df.subset <- coefs.defining.top.topic.df %>% subset(grepl("^ENSG",ENSGID)) %>% group_by(ENSGID) %>% arrange(desc(gene.feature_x_beta)) %>% slice(1:10) %>% as.data.frame
all.coefs.defining.top.topic.df.subset <- all.coefs.defining.top.topic.df %>% subset(grepl("^ENSG",ENSGID)) %>% group_by(ENSGID) %>% arrange(desc(gene.feature_x_beta)) %>% slice(1:10) %>% as.data.frame


PoPS_preds.importance.score <- merge(preds.combined.df, all.coefs.defining.top.topic.df.subset %>% select(-Gene.name, -Gene, -EntrezID), by="ENSGID") %>% merge(., metadata, by.x="topic", by.y="X", all.x=T)
colnames(PoPS_preds.importance.score)[which(colnames(PoPS_preds.importance.score)=="topic")] <- "pathway"
write.table(PoPS_preds.importance.score %>% apply(2, as.character), file=paste0(OUTDIR, "/PoPS_preds.importance.score.all.columns.txt"), sep="\t", quote=F, row.names=F)
PoPS_preds.importance.score.key <- PoPS_preds.importance.score %>% select(Gene, Gene.name, PoPS_Score_with.cNMF, PoPS_Score_without.cNMF, pathway, Long_Name, gene.feature_x_beta)
write.table(PoPS_preds.importance.score.key %>% apply(2, as.character), file=paste0(OUTDIR, "/PoPS_preds.importance.score.key.columns.txt"), sep="\t", quote=F, row.names=F)

 


## get top features' top genes
top.features.names <- coefs.df %>% arrange(desc(beta)) %>% slice(1:10) %>% pull(parameter)
top.features.definition <- all.features.cNMF %>% select(all_of(c("ENSGID",top.features.names))) %>% add_gene_name_to_df
slice_top_genes_from_features <- function(df) {
    out <- df %>% melt(value.name="relative.importance", variable.name="features", id.vars=c("ENSGID", "Gene", "Gene.name", "EntrezID")) %>% group_by(features) %>% arrange(desc(relative.importance)) %>% slice(1:10)
}
top.genes.in.top.features <- top.features.definition %>% slice_top_genes_from_features %>% as.data.frame %>% merge(coefs.df, by.x="features", by.y="parameter") %>% arrange(desc(beta))
write.table(top.genes.in.top.features, file=paste0(OUTDIR, "/top.genes.in.top.features.coefs.txt"), row.names=F, quote=F, sep="\t")
