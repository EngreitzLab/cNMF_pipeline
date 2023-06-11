## Helen Kang
## Script to output all supplementary tables
## 220630

## setwd("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/SuppTables/")

suppressPackageStartupMessages(library(conflicted))
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("first", "dplyr")

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(reshape2)
    ## library(ggplot2)
    ## library(cowplot)
    ## library(ggpubr) ## ggarrange
    ## library(gplots) ## heatmap.2
    ## library(scales) ## geom_tile gradient rescale
    ## library(ggrepel)
    library(stringr)
    library(stringi)
    library(svglite)
    library(Seurat)
    library(SeuratObject)
    library(xlsx)
})


##########################################################################################
## Constants and Directories
option.list <- list(
    make_option("--sampleName", type="character", default="2kG.library", help="Name of the sample"),
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/analysis/all_genes/2kG.library.ctrl.only/K25/threshold_0_2/", help="Output directory"),
    make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--raw.mtx.RDS.dir",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.RDS", help="input matrix to cNMF pipeline"),
    make_option("--perturbSeq", type="logical", default=TRUE, help="Whether this is a Perturb-seq experiment")
)
opt <- parse_args(OptionParser(option_list=option.list))


## ## K562 gwps 2k overdispersed genes
## ## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/WeissmanK562gwps/K90/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/"
## opt$K.val <- 90
## opt$sampleName <- "WeissmanK562gwps"
## opt$perturbSeq <- TRUE

## ## ENCODE mouse heart
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/analysis/top2000VariableGenes/mouse_ENCODE_heart/K10/threshold_0_2/"
## opt$sampleName <- "mouse_ENCODE_heart"
## opt$K.val <- 10
## opt$perturbSeq <- FALSE


OUTDIR <- opt$outdir
SAMPLE <- opt$sampleName
k <- opt$K.val
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
check.dir <- c(OUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))
fdr.thr <- 0.05

db <- ifelse(grepl("mouse", SAMPLE), "org.Mm.eg.db", "org.Hs.eg.db")
suppressPackageStartupMessages(library(!!db)) ## load the appropriate database

## Load Data
## Batch Correlation table
load(paste0(OUTDIR, "/batch.correlation.RDS"))
if("topic" %in% colnames(max.batch.correlation.df)){
    max.batch.correlation.df <- max.batch.correlation.df %>%
        mutate(ProgramID = paste0("K", k, "_", gsub("topic_", "", topic))) %>%
        as.data.frame
}
if ("ProgramID" %in% colnames(max.batch.correlation.df)) {
    max.batch.correlation.df <- max.batch.correlation.df %>%
        mutate(ProgramID = paste0("K", k, "_", gsub("topic_", "", ProgramID))) %>%
        as.data.frame
}


##################################################
## load cNMF results
cNMF.result.file <- paste0(OUTDIR,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
if(file.exists(cNMF.result.file)) {
    message(paste0("loading cNMF result file: \n", cNMF.result.file))
    load(cNMF.result.file)
} else {
	print(paste0(cNMF.result.file, " does not exist"))
}


## helper function to map between ENSGID and SYMBOL
map.ENSGID.SYMBOL <- function(df) {
    ## need column `Gene` to be present in df
    ## detect gene data type (e.g. ENSGID, Entrez Symbol)
    if(!("Gene" %in% colnames(df))) df$Gene = df$ENSGID
    gene.type <- ifelse(nrow(df) == sum(as.numeric(grepl("^ENS", df$Gene))),
                        "ENSGID",
                        "Gene")
    if(gene.type == "ENSGID") {
        mapped.genes <- mapIds(get(db), keys=df$Gene, keytype = "ENSEMBL", column = "SYMBOL")
        df <- df %>% mutate(ENSGID = Gene, Gene = mapped.genes)
    } else {
        mapped.genes <- mapIds(get(db), keys=df$Gene, keytype = "SYMBOL", column = "ENSEMBL")
        df <- df %>% mutate(ENSGID = mapped.genes)
    }
    df <- df %>% mutate(Gene = ifelse(is.na(Gene), ENSGID, Gene))
    return(df)
}

if(sum(c("ENSGID", "Gene") %in% colnames(median.spectra.zscore.df)) < 2) {
    median.spectra.zscore.df <- median.spectra.zscore.df %>% map.ENSGID.SYMBOL
}
    
## get list of topic defining genes by z-score coefficients
if(!("theta.zscore.rank.df" %in% ls())) {
    theta.zscore.rank.list <- vector("list", ncol(theta.zscore))## initialize storage list
    for(i in 1:ncol(theta.zscore)) {
        topic <- paste0("topic_", colnames(theta.zscore)[i])
        theta.zscore.rank.list[[i]] <- theta.zscore %>%
            as.data.frame %>%
            select(all_of(i)) %>%##here
            `colnames<-`("topic.zscore") %>%
            mutate(Gene = rownames(.)) %>%
            arrange(desc(topic.zscore), .before="topic.zscore") %>%
            mutate(zscore.specificity.rank = 1:n()) %>% ## add rank column
            mutate(Topic = topic) ## add topic column
    }
    theta.zscore.rank.df <- do.call(rbind, theta.zscore.rank.list) %>%  ## combine list to df
        `colnames<-`(c("topic.zscore", "Gene", "zscore.specificity.rank", "ProgramID")) %>%
        mutate(ProgramID = gsub("topic_", paste0("K", k, "_"), ProgramID)) %>%
        as.data.frame %>% map.ENSGID.SYMBOL
}


if(!("theta.tpm.rank.df" %in% ls())) {
    theta.tpm.rank.list <- vector("list", ncol(theta.raw))## initialize storage list
    for(i in 1:ncol(theta.zscore)) {
        topic <- paste0("topic_", colnames(theta.raw)[i])
        theta.tpm.rank.list[[i]] <- theta %>%
            as.data.frame %>%
            select(all_of(i)) %>%
            `colnames<-`("topic.zscore") %>%
            mutate(Gene = rownames(.)) %>%
            arrange(desc(topic.zscore), .before="topic.zscore") %>%
            mutate(zscore.specificity.rank = 1:n()) %>% ## add rank column
            mutate(Topic = topic) ## add topic column
    }
    theta.tpm.rank.df <- do.call(rbind, theta.tpm.rank.list) %>%  ## combine list to df
        `colnames<-`(c("program.tpm.coef", "Gene", "tpm.coef.rank", "ProgramID")) %>%
        mutate(ProgramID = gsub("topic_", paste0("K", k, "_"), ProgramID)) %>%
        as.data.frame %>% map.ENSGID.SYMBOL %>% mutate(Gene = ifelse(is.na(Gene), ENSGID, Gene))
}

## get list of topic genes by median spectra weight
if(!("median.spectra.rank.df" %in% ls())) {
    median.spectra.rank.list <- vector("list", ncol(median.spectra))## initialize storage list
    for(i in 1:ncol(median.spectra)) {
        topic <- paste0("topic_", colnames(median.spectra)[i])
        median.spectra.rank.list[[i]] <- median.spectra %>%
            as.data.frame %>%
            select(all_of(i)) %>%
            `colnames<-`("median.spectra") %>%
            mutate(Gene = rownames(.)) %>%
            arrange(desc(median.spectra), .before="median.spectra") %>%
            mutate(median.spectra.rank = 1:n()) %>% ## add rank column
            mutate(Topic = topic) ## add topic column
    }
    median.spectra.rank.df <- do.call(rbind, median.spectra.rank.list) %>%  ## combine list to df
        `colnames<-`(c("median.spectra", "Gene", "median.spectra.rank", "ProgramID")) %>%
        mutate(ProgramID = gsub("topic_", paste0("K", k, "_"), ProgramID)) %>%
        as.data.frame %>% map.ENSGID.SYMBOL %>% mutate(Gene = ifelse(is.na(Gene), ENSGID, Gene))
    ## median.spectra.zscore.df <- median.spectra.zscore.df %>% mutate(Gene = ENSGID) ## quick fix, need to add "Gene" column to this dataframe in analysis script
}

## Load Regulator Data
if(opt$perturbSeq) {
    MAST.file.name <- paste0(OUTDIR, "/", SAMPLE, "_MAST_DEtopics.txt")
    message(paste0("loading ", MAST.file.name))
    MAST.df <- read.delim(MAST.file.name, stringsAsFactors=F, check.names=F)
    if(grepl("topic", MAST.df$primerid) %>% sum > 0) MAST.df <- MAST.df %>% mutate(ProgramID = paste0("K", k, "_", gsub("topic_", "", primerid))) %>% as.data.frame    
}

## Load Promoter and Enhancer TF Motif Enrichment Data
add.ProgramID <- function(df) {
    if("topic" %in% c(df %>% colnames)) {
        return(df %>% mutate(ProgramID = paste0('K', k, '_', gsub('topic_', '', topic))) %>% as.data.frame)
    }
}
num.top.genes <- 300
## all.ttest.df.path <- paste0(OUTDIRSAMPLE,"/", ep.type, ".topic.top.", num.top.genes, ".zscore.gene_motif.count.ttest.enrichment_motif.thr.", motif.match.thr.str, "_", SUBSCRIPT.SHORT,".txt")
promoter.ttest.df.path <- paste0(OUTDIR,"/", "promoter", ".topic.top.", num.top.genes, ".zscore.gene_motif.count.ttest.enrichment_motif.thr.", "pval1e-4", "_", SUBSCRIPT.SHORT,".txt")
enhancer.ttest.df.path <- paste0(OUTDIR,"/", "promoter", ".topic.top.", num.top.genes, ".zscore.gene_motif.count.ttest.enrichment_motif.thr.", "pval1e-6", "_", SUBSCRIPT.SHORT,".txt")
promoter.ttest.df <- read.delim(promoter.ttest.df.path, stringsAsFactors=F) %>% add.ProgramID
enhancer.ttest.df <- read.delim(enhancer.ttest.df.path, stringsAsFactors=F) %>% add.ProgramID


## GO terms
## file.name <- paste0(OUTDIRSAMPLE, "/clusterProfiler_GeneRankingType", ranking.type.here, "_EnrichmentType", GSEA.type,".txt")
file.name <- paste0(OUTDIR, "/clusterProfiler_GeneRankingType", "zscore", "_EnrichmentType", "GOEnrichment",".txt")
theta.zscore.GO.df <- read.delim(file.name, stringsAsFactors=F)
file.name <- paste0(OUTDIR, "/clusterProfiler_GeneRankingType", "median_spectra_zscore", "_EnrichmentType", "GOEnrichment",".txt")
median_spectra_zscore.GO.df <- read.delim(file.name, stringsAsFactors=F)


####################################################################################################
## Gather all columns
## MaxBatchCorrelation
MaxBatchCorrelation <- max.batch.correlation.df %>%
    select(ProgramID, maxPearsonCorrelation) %>%
    `colnames<-`(c("ProgramID", "MaxBatchCorrelation")) %>%
    as.data.frame

## nSigPerturbationsProgramUp
SigPerturbations.df <- MAST.df %>%
    subset(fdr.across.ptb < fdr.thr) %>%
    as.data.frame
nSigPerturbationsProgramUp <- SigPerturbations.df %>%
    subset(coef > log(1.1)) %>%
    group_by(ProgramID) %>%
    summarize(nSigPerturbationsProgramUp = n()) %>%
    as.data.frame
nSigPerturbationsProgramDown <- SigPerturbations.df %>%
    subset(coef < log(0.9)) %>%
    group_by(ProgramID) %>%
    summarize(nSigPerturbationsProgramDown = n()) %>%
    as.data.frame

## nSigMotifsPromoter
nSigMotifsPromoter <- promoter.ttest.df %>%
    subset(one.sided.p.adjust < fdr.thr) %>%
    group_by(ProgramID) %>%
    summarize(nSigMotifsPromoter = n()) %>%
    as.data.frame

## nSigMotifsEnhancer
nSigMotifsEnhancer <- enhancer.ttest.df %>%
    subset(one.sided.p.adjust < fdr.thr) %>%
    group_by(ProgramID) %>%
    summarize(nSigMotifsEnhancer = n()) %>%
    as.data.frame

## ProgramGenesZScoreCoefficientTop10
ProgramGenesZScoreCoefficientTop10 <- theta.zscore.rank.df %>%
    mutate(Gene = ifelse(is.na(Gene), ENSGID, Gene)) %>%
    subset(zscore.specificity.rank < 10) %>%
    group_by(ProgramID) %>%
    arrange(desc(topic.zscore)) %>%
    summarize(ProgramGenesZScoreCoefficientTop10 = paste0(Gene, collapse=",")) %>%
    as.data.frame

## ProgramGenesTPMCoefficientTop10
ProgramGenesTPMCoefficientTop10 <- theta.tpm.rank.df %>%
    subset(tpm.coef.rank < 10) %>%
    group_by(ProgramID) %>%
    arrange(desc(program.tpm.coef)) %>%
    summarize(ProgramGenesTPMCoefficientTop10 = paste0(Gene, collapse=",")) %>%
    as.data.frame

## ProgramGenesMedianSpectraTop10
ProgramGenesMedianSpectraTop10 <- median.spectra.rank.df %>%
    subset(median.spectra.rank < 10) %>%
    mutate(Gene = ifelse(is.na(Gene), ENSGID, Gene)) %>%
    group_by(ProgramID) %>%
    arrange(desc(median.spectra)) %>%
    summarize(ProgramGenesMedianSpectraTop10 = paste0(Gene, collapse = ",")) %>%
    as.data.frame
    
## ProgramGenesMedianSpectraZScoreTop10
ProgramGenesMedianSpectraZScoreTop10 <- median.spectra.zscore.df %>%
    subset(median.spectra.zscore.rank < 10) %>%
    group_by(ProgramID) %>%
    arrange(desc(median.spectra.zscore)) %>%
    summarize(ProgramGenesMedianSpectraZScoreTop10 = paste0(Gene, collapse = ",")) %>%
    as.data.frame

## ProgramGenesMotifsPromoter
ProgramGenesMotifsPromoter <- promoter.ttest.df %>%
    subset(one.sided.p.adjust < fdr.thr) %>%
    group_by(ProgramID) %>%
    arrange(one.sided.p.adjust) %>%
    summarize(ProgramGenesMotifsPromoter = paste0(motif, collapse = ",")) %>%
    as.data.frame

## ProgramGenesMotifsEnhancer
ProgramGenesMotifsEnhancer <- enhancer.ttest.df %>%
    subset(one.sided.p.adjust < fdr.thr) %>%
    group_by(ProgramID) %>%
    arrange(one.sided.p.adjust) %>%
    summarize(ProgramGenesMotifsEnhancer = paste0(motif, collapse = ",")) %>%
    as.data.frame

## ProgramGenesZScoreCoefficientGOTermsTop10
ProgramGenesZScoreCoefficientGOTermsTop10 <- theta.zscore.GO.df %>%
    group_by(ProgramID) %>%
    arrange(fdr.across.ont) %>%
    slice(1:10) %>%
    mutate(GOTerm = paste0(ONTOLOGY, ":", ID, ":", Description)) %>%
    summarize(ProgramGenesZScoreCoefficientGOTermsTop10 = paste0(GOTerm, collapse = ",")) %>%
    as.data.frame

## ProgramGenesMedianSpectraZScoreGOTermsTop10
ProgramGenesMedianSpectraZScoreGOTermsTop10 <- median_spectra_zscore.GO.df %>%
    group_by(ProgramID) %>%
    arrange(fdr.across.ont) %>%
    slice(1:10) %>%
    mutate(GOTerm = paste0(ONTOLOGY, ":", ID, ":", Description)) %>%
    summarize(ProgramGenesMedianSpectraZScoreGOTermsTop10 = paste0(GOTerm, collapse = ",")) %>%
    as.data.frame
    


## ProgramGenesTF*


## Combine all columns
ProgramSummary.list <- list(
    MaxBatchCorrelation = MaxBatchCorrelation,
    nSigPerturbationsProgramUp = nSigPerturbationsProgramUp,
    nSigPerturbationsProgramDown = nSigPerturbationsProgramDown,
    nSigMotifsPromoter = nSigMotifsPromoter,
    nSigMotifsEnhancer = nSigMotifsEnhancer,
    ProgramGenesZScoreCoefficientTop10 = ProgramGenesZScoreCoefficientTop10,
    ProgramGenesTPMCoefficientTop10 = ProgramGenesTPMCoefficientTop10,
    ProgramGenesMedianSpectraTop10 = ProgramGenesMedianSpectraTop10,
    ProgramGenesMedianSpectraZScoreTop10 = ProgramGenesMedianSpectraZScoreTop10,
    ProgramGenesMotifsPromoter = ProgramGenesMotifsPromoter,
    ProgramGenesMotifsEnhancer = ProgramGenesMotifsEnhancer,
    ProgramGenesZScoreCoefficientGOTermsTop10 = ProgramGenesZScoreCoefficientGOTermsTop10,
    ProgramGenesMedianSpectraZScoreGOTermsTop10 = ProgramGenesMedianSpectraZScoreGOTermsTop10
)
    
ProgramSummary.df <- Reduce(function(x, y, ...) full_join(x, y, by = "ProgramID", ...), ProgramSummary.list)

## ## Write Table to Text File
fileName <- paste0(SAMPLE, "_ProgramSummary_", SUBSCRIPT.SHORT)
write.table(ProgramSummary.df, file=paste0(OUTDIR, "/", fileName, ".xlsx"), row.names=F, quote=F, sep="\t")
## Write Excel File
write.xlsx(ProgramSummary.df, file=paste0(OUTDIR, "/", fileName, ".xlsx"), sheetName="Gene Selection Table (For experimental design)", row.names=F, showNA=F)

