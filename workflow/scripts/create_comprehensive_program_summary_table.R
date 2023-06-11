## Helen Kang
## 220316
## Make Summary Table

library(conflicted)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("first", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("rename", "dplyr")

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(reshape2)
    library(ggplot2)
    library(ggpubr) ## ggarrange
    library(gplots) ## heatmap.2
    library(ggrepel)
    library(readxl)
    library(xlsx) ## might not need this package
    library(writexl)
    library(org.Hs.eg.db)
})


option.list <- list(
    make_option("--sampleName", type="character", default="2kG.library", help="Name of the sample"),
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/analysis/all_genes/2kG.library.ctrl.only/K25/threshold_0_2/", help="Output directory"),
    make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--perturbSeq", type="logical", default=TRUE, help="Whether this is a Perturb-seq experiment")
)
opt <- parse_args(OptionParser(option_list=option.list))



## K562 gwps 2k overdispersed genes
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/WeissmanK562gwps/K90/"
opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/"
opt$K.val <- 90
opt$sampleName <- "WeissmanK562gwps"
opt$perturbSeq <- TRUE
opt$scratch.outdir <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/K90/analysis/comprehensive_program_summary/"
opt$barcodeDir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/data/K562_gwps_raw_singlecell_01_metadata.txt"

## OUTDIR <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220316_regulator_topic_definition_table/outputs/"
## FIGDIR <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220316_regulator_topic_definition_table/figures/"
## SCRATCHOUTDIR <- "/scratch/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220316_regulator_topic_definition_table/outputs/"
OUTDIR <- opt$outdir
SCRATCHOUTDIR <- opt$scratch.outidr
check.dir <- c(OUTDIR, SCRATCHOUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

mytheme <- theme_classic() + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, face = "bold"))
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)



## parameters
## OUTDIRSAMPLE <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/"
OUTDIRSAMPLE <- opt$outdir
k <- opt$K.val 
SAMPLE <- opt$sampleName
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)

## load ref table (for Perturbation distance to GWAS loci annotation in Perturb_plus column)
ref.table <- read.delim(file="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/ref.table.txt", header=T, check.names=F, stringsAsFactors=F)

## ## load test results
## all.test.combined.df <- read.delim(paste0("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220312_compare_statistical_test/outputs/all.test.combined.df.txt"), stringsAsFactors=F)
## all.test.MAST.df <- all.test.combined.df %>% subset(test.type == "batch.correction")

## MAST.df.4n.input <- read.delim(paste0("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/2kG.library4n3.99x_MAST.txt"), stringsAsFactors=F, check.names=F)
## MAST.df.4n <- MAST.df.4n.input %>% ## remove multiTarget entries
##     subset(!grepl("multiTarget", perturbation)) %>%
##     group_by(zlm.model.name) %>%
##     mutate(fdr.across.ptb = p.adjust(`Pr(>Chisq)`, method="fdr")) %>%
##     ungroup() %>%
##     subset(zlm.model.name == "batch.correction") %>%
##     select(-zlm.model.name) %>%
##     as.data.frame
## colnames(MAST.df.4n) <- c("Topic", "p.value", "log2FC", "log2FC.ci.hi", "log2fc.ci.lo", "fdr", "Perturbation", "fdr.across.ptb")

## Load Regulator Data
if(opt$perturbSeq) {
    MAST.file.name <- paste0(OUTDIR, "/", SAMPLE, "_MAST_DEtopics.txt")
    message(paste0("loading ", MAST.file.name))
    MAST.df <- read.delim(MAST.file.name, stringsAsFactors=F, check.names=F) %>% rename("Perturbation" = "perturbation") %>%
        rename("log2FC" = "coef", "log2FC.ci.hi" = ci.hi, "log2FC.ci.lo" = ci.lo, "p.value" = "Pr(>Chisq)")
    if(grepl("topic", MAST.df$primerid) %>% sum > 0) MAST.df <- MAST.df %>% mutate(ProgramID = paste0("K", k, "_", gsub("topic_", "", primerid))) %>% as.data.frame    
}



## ## 2n MAST
## file.name <- paste0("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes//Perturb_2kG_dup4/acrossK/aggregated.outputs.findK.perturb-seq.RData")
## load(file.name)
## MAST.df.2n <- MAST.df %>%
##     filter(K == 60) %>%
##     select(-K) %>%
##     select(-zlm.model.name)
## colnames(MAST.df.2n) <- c("Topic", "p.value", "log2FC", "log2FC.ci.hi", "log2fc.ci.lo", "fdr", "Perturbation", "fdr.across.ptb")



## load gene annotations (Refseq + Uniprot)
gene.summaries.path <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/combined.gene.summaries.txt"
gene.summary <- read.delim(gene.summaries.path, stringsAsFactors=F)


## load topic model results
cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
}



## modify theta.zscore if Gene is in ENSGID
db <- ifelse(grepl("mouse", SAMPLE), "org.Mm.eg.db", "org.Hs.eg.db")
gene.ary <- theta.zscore %>% rownames
if(grepl("^ENSG", gene.ary) %>% as.numeric %>% sum == nrow(theta.zscore)) {
    GeneSymbol.ary <- mapIds(get(db), keys=gene.ary, keytype = "ENSEMBL", column = "SYMBOL")
    GeneSymbol.ary[is.na(GeneSymbol.ary)] <- row.names(theta.zscore)[is.na(GeneSymbol.ary)]
    rownames(theta.zscore) <- GeneSymbol.ary
}

## omega.4n <- omega
## theta.zscore.4n <- theta.zscore
## theta.raw.4n <- theta.raw
meta_data <- read.delim(opt$barcodeDir, stringsAsFactors=F) 

## ## batch topics
## batch.topics.4n <- read.delim(file=paste0(OUTDIRSAMPLE, "/batch.topics.txt"), stringsAsFactors=F) %>% as.matrix %>% as.character

ann.omega <- merge(meta_data, omega, by.x="CBC", by.y=0, all.T)


##########################################################################################
## create table
create_topic_definition_table <- function(theta.zscore, t) {
    out <- theta.zscore[,t] %>%
        as.data.frame %>%
        `colnames<-`(c("zscore")) %>%
        mutate(Perturbation = rownames(theta.zscore), .before="zscore") %>%
        merge(gene.summary, by.x="Perturbation", by.y="Gene", all.x=T) %>%
        arrange(desc(zscore)) %>%
        mutate(Rank = 1:n(), .before="Perturbation") %>%
        mutate(ProgramID = paste0("K", k, "_", t), .before="zscore") %>%
        arrange(Rank) %>%
        mutate(My_summary = "", .after = "zscore") %>%
        select(Rank, ProgramID, Perturbation, zscore, My_summary, FullName, Summary)
}

create_topic_regulator_table <- function(all.test, program.here, fdr.thr = 0.1) {
    out <- MAST.df %>%
        subset(ProgramID == program.here &
               fdr.across.ptb < fdr.thr) %>%
        select(Perturbation, fdr.across.ptb, log2FC, log2FC.ci.hi, log2FC.ci.lo, fdr, p.value) %>%
        merge(gene.summary, by.x="Perturbation", by.y="Gene", all.x=T) %>%
        arrange(fdr.across.ptb, desc(log2FC)) %>%
        mutate(Rank = 1:n(), .before="Perturbation") %>%
        mutate(My_summary = "", .after="Perturbation") %>%
        mutate(ProgramID = program.here, .after="Rank") %>%
        ## merge(., ref.table %>% select("Symbol", "TSS.dist.to.SNP", "GWAS.classification"), by.x="Perturbation", by.y="Symbol", all.x=T) %>%
        ## mutate(EC_ctrl_text = ifelse(.$GWAS.classification == "EC_ctrls", "(+)", "")) %>%
        ## mutate(GWAS.class.text = ifelse(grepl("CAD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb"),
        ##                          ifelse(grepl("IBD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb_IBD"), ""))) %>%
        ## mutate(Perturb_plus = paste0(Perturbation, GWAS.class.text, EC_ctrl_text)) %>%
    select(Rank, ProgramID, Perturbation, fdr.across.ptb, log2FC, My_summary, FullName, Summary, log2FC.ci.hi, log2FC.ci.lo, fdr, p.value) %>% ## removed Perturb_plus
        arrange(Rank)
}



create_summary_table <- function(ann.omega, theta.zscore, all.test, meta_data) {
    df.list <- vector("list", k)
    for (t in 1:k) {
        program.here <- paste0("K", k, "_", t)

        ## topic defining genes
        ann.theta.zscore <- theta.zscore %>% create_topic_definition_table(t)
        ann.top.theta.zscore <- ann.theta.zscore %>% subset(Rank <= 100) ## select the top 100 topic defining genes to output
 
        ## regulators
        if(opt$perturbSeq) regulator.MAST.df <- all.test %>% create_topic_regulator_table(program.here, 0.3)
 
        ## write table to scratch dir
        file.name <- paste0(SCRATCHOUTDIR, program.here, "_table.csv")
        sink(file=file.name) ## open the document
        ## cat("Author,PERTURBATIONS SUMMARIES\n\"\n\"\n\"\n\"\n\"\n\"\n\"\n\n\n\nAuthor,TOPIC SUMMARIES\n\"\n\"\n\"\n\"\n\"\n\"\n\"\n\nAuthor,TESTABLE HYPOTHESIS IDEAS:\n\"\n\"\n\"\n\"\n\"\n\"\n\"\n\nAuthor,OTHER THOUGHTS:\n\"\n\"\n\"\n\"\n\"\n\"\n\"\n\nTOPIC DEFINING GENES (TOP 100),\n") ## headers
        cat("Author,PERTURBATIONS SUMMARIES,,,,,,,,,,,\n\"\n\"\n\"\n\"\n\"\n\"\n\"\n\"\nAuthor,TOPIC SUMMARIES\n\"\n\"\n\"\n\"\n\"\n\"\n\"\n\"\nAuthor,TESTABLE HYPOTHESIS IDEAS:\n\"\n\"\n\"\n\"\n\"\n\"\n\"\n\"\nAuthor,OTHER THOUGHTS:\n\"\n\"\n\"\n\"\n\"\n\"\n\"\n\"\n\nTOPIC DEFINING GENES (TOP 100),\n") ## headers
        write.csv(ann.top.theta.zscore, row.names=F) ## topic defining genes
        cat("\n\"\n\"\nPERTURBATIONS REGULATING TOPIC AT FDR < 0.3 (most significant on top),,,,,,,,,,,\n") ## headers

        if(opt$perturbSeq) write.csv(regulator.MAST.df, row.names=F) ## regulators
        sink() ## close the document
         
        ## read the assembled table to save to list
        df.list[[t]] <- read.delim(file.name, stringsAsFactors=F, check.names=F, sep=",")
    }


    ## output to xlsx
    names(df.list) <- paste0("Program ", 1:k)
    write_xlsx(df.list, paste0(OUTDIR, "/", SAMPLE, "_k_", k, ".dt_", DENSITY.THRESHOLD, "_ComprehensiveProgramSummary.xlsx"))

    return(df.list)
}



df <- create_summary_table(ann.omega, theta.zscore, MAST.df, meta_data)
