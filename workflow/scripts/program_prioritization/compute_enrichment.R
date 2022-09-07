# Q:  Among the GWAS signals likely to act in endothelial cells, 
#     are the likely genes (1-2 closest to EC peak) enriched for regulating or being expressed in certain topics? 
#
# To formalize this last analysis, here’s what I think we would need to do: 
# 1. Define GWAS signals that may act in endothelial cells:
#    - Filter to credible sets without association to lipids 
#      (e.g. based on PheWAS from Aragam table for now; ideally based on coloc analysis)
#    - Filter to credible sets with noncoding variant in EC peak or ABC enhancer link to ECs 
#      (i.e., where RankOfDistanceToECPeakWithVariant is not NA or MaxABC.Rank.ECOnly is not NA)
#    - Add signals with coding or splice site variants? (e.g. CCM2)
# 2. Test for enrichment of each topic among the top genes in these GWAS signals: 
#    - For each topic, get the set of genes that are either IN the topic (top 100 by z-score for now; may also want to try top 300) or regulate the topic (FDR < 0.3 for now; may also want to try <0.1).  
#    - Then, test for enrichment among the genes closest to EC peak (e.g. top 2) — i.e., where there is convergent evidence between genes in a pathway and regulatory links from endothelial cells enhancers — either using Fisher’s exact test, or by permuting the topic-gene assignments (e.g., permute the “TopicsRegulatedByThisGene” and “TopicInWhichGeneInTop100ZScoreSpecificGenes” columns together to new genes.  
#    - As the background set, consider all (perturbed?) genes within 1 Mb of these EC loci, or all perturbed genes in GWAS loci (CAD or IBD) 
#    - Consider how to deal with loci with multiple independent signals that point to the same gene (e.g. LDLR, COL4A1), which might lead to ‘double-counting’ some genes depending on how we do the test


## Possible background models:
## 1.  All expressed genes in TeloHAECs (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to a random expressed gene in TeloHAEC?)
## 2.  Other genes in EC CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 3.  Other genes in CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 4.  Other expressed genes in EC CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 5.  Other expressed genes in CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 6.  Other expressed, 2000Kperturbed genes in EC CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 7.  Other expressed, 2000Kperturbed genes in CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)

## Possible tests:
## 	- Fisher's exact test
## 	- Permute the gene labels (e.g. among expressed perturbed genes in EC CAD GWAS loci) with respect to all of the other data
##        Make a [Gene,GeneInTopic,GeneRegulatesTopic] %>% unique() data frame, and shuffle it
##        Merge the shuffled data frame back into [CredibleSet, Gene, RankInGWASLocus] %>% unique()
##        Count the number of overlaps sum((RankInGWASLocus <= 2 & (topic %in% GeneRegulatesTopic | topic %in% GeneInTopic)))
##  	  Do that 1000 times, create a distribution.  P-value = fraction of permutations that equal or exceed the observed value

## Resource == "Aragam2021"
##  keep Harst where there is not a Aragam2021 lead SNP within 20Kb
##    /oak/stanford/groups/engreitz/Papers/ECPerturbSeq/Analysis/211219_GenePrioritization/ECPerturbSeq2021-Analysis/GWAS_tables/




## Libraries
library(conflicted)
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("first", "dplyr")
packages <- c("optparse","dplyr","data.table","reshape2","ggplot2","ggpubr","conflicted",
              "cluster","textshape","readxl","writexl","tidyr","org.Hs.eg.db","stats",
              "gplots", "stringi" # heatmap.2
              ) 
xfun::pkg_attach(packages)
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("replace_na", "dplyr")
conflict_prefer("filter", "dplyr")

## optparse
option.list <- list(
    make_option("--input.GWAS.table", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/GWAS_tables/CAD_GWAS_gene_incl_ubq_genes.txt", help="Input CAD GWAS Table"),
    make_option("--cNMF.PoPS.input.table", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/analysis/all_genes/2kG.library.ctrl.only/K60/threshold_0_2/cNMF.PoPS.input.to.large.table.pipeline_k_60.dt_0_2.CAD_aug6.txt", help="Input table with (Regulator) and program gene information"),
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/analysis/all_genes/2kG.library.ctrl.only/K60/threshold_0_2/", help="Output directory"),
    make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
    # make_option("--prefix", type="character", default="", help="PoPS prefix"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--num.tests", type="numeric", default=2, help="number of statistical test to do from the top of statistical.test.df.txt"),
    make_option("--statistical.test.list", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/GWAS_tables/statistical.test.list.txt", help="table to specify program prioritization tests to perform"),
    make_option("--perturbseq", type="logical", default=F)
)
opt <- parse_args(OptionParser(option_list=option.list))


## Directories
## FIGDIR <- "./figures/"
OUTDIR <- opt$outdir
check.dir <- c(OUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))


DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
# PREFIX <- opt$prefix

## load batch program list
batch.topics <- read.delim(opt$batch.topics, stringsAsFactors=F)


## Perturb-seq vs 10X names (from Gavin)
## ptb10xNames.df <- read_xlsx("../data/Perturbation 10X names.xlsx") %>% as.data.frame
ptb10xNames.df <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/220627_add_Perturbation 10X names.txt", stringsAsFactors=F, check.names=F)
## write.table(ptb10xNames.df,"../data/220627_add_Perturbation 10X names.220628.txt", quote=F, row.names=F, sep="\t")
perturbseq.gene.names.to10X <- function(Gene) stri_replace_all_regex(Gene, pattern = ptb10xNames.df$Symbol, replace = ptb10xNames.df$`Name used by CellRanger`, vectorize=F)
tenX.gene.names.toperturbseq <- function(Gene) stri_replace_all_regex(Gene, pattern = ptb10xNames.df$`Name used by CellRanger`, replace = ptb10xNames.df$Symbol, vectorize=F) 


## Load Data
gtf.10X.df <- readRDS("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/refdata-cellranger-arc-GRCh38-2020-A_genes.gtf_df.RDS") ## load 10X gtf file

## PP.GWAS.df <- read.delim("PP_GWAS_gene_incl_ubq_genes.txt", header=T, stringsAsFactors=F)
CAD.GWAS.df <- CAD.GWAS.df.original <- read.delim(opt$input.GWAS.table, header=T, stringsAsFactors=F, sep="\t") %>% ## Rosa's table
    mutate(original.gene = gene) %>%
    mutate(gene = gene %>% perturbseq.gene.names.to10X)
CAD.GWAS.df <- CAD.GWAS.df %>% merge(gtf.10X.df %>% select(gene_name, gene_type) %>% unique, by.x="gene", by.y="gene_name", all.x=T)
## write.table(CAD.GWAS.df %>% subset(gene_type != "protein_coding"), paste0(OUTDIR, "/nonProteinCodingGene.CAD.GWAS.df.txt"), sep="\t", quote=F, row.names=F)
CAD.GWAS.df.duplicatedRows <- CAD.GWAS.df.original %>% subset(gene %in% CAD.GWAS.df.original$gene[CAD.GWAS.df.original$gene %>% duplicated])
## write.table(CAD.GWAS.df.duplicatedRows, paste0(OUTDIR, "/duplicated_CAD_GWAS_gene_incl_ubq_genes.txt"), sep="\t", quote=F, row.names=F) ## check for duplicated genes
IBD.GWAS.df <- read_excel("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/200511_2000_Gene_Library_partial_IBD_loci_selection.xlsx", skip=22) %>% select(c(29:72)) %>% as.data.frame %>%
    rowwise %>%
    mutate(Gene = strsplit(Symbol_txptIDs, split=";") %>% sapply(`[[`, 1) %>% unlist,
           DistanceToTSS = TSS_start - BestSNPPos) %>%
    select(Gene, CredibleSet, BestSNP, TSS_chr, TSS_start, DistanceToTSS) %>%
    setnames(old = c("Gene", "TSS_chr", "TSS_start", "DistanceToTSS"), new = c("gene", "chr", "TSS_position", "SNP_to_TSS_distance")) %>%
    group_by(CredibleSet) %>%
    arrange(abs(SNP_to_TSS_distance)) %>%
    mutate(rank_SNP_to_TSS = 1:n()) %>%
    subset(abs(SNP_to_TSS_distance) < 500000) %>%
    as.data.frame


## Load data for all genes relation to the topics (not limited to CAD GWAS genes)
narrow.df <- read.delim(opt$cNMF.PoPS.input.table, header=T, stringsAsFactors=F)## %>% select(-Top5FeaturesThatContributeToPoPSScore)
narrow.df <- narrow.df[narrow.df$Gene!="NULL", ]

## Gavin's full bulk RNA-seq TPM table
## RNAseq_CITV <- read_xlsx(paste0("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/200815_subset_RNAseq_CITV1.xlsx"), sheet="RNAseq") ## old
## TPM.df <- RNAseq_CITV %>% select(`...1`, Telo_C_1, Telo_C_2, Telo_C_3) %>% `colnames<-`(c("Gene", "TeloC1", "TeloC2", "TeloC3")) %>% mutate(TeloHAEC_ctrl_TPM = (TeloC1 + TeloC2 + TeloC3) / 3) %>% as.data.frame ## old
RNAseq_CITV <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/CAD_SNP_INFO/RNASeq_Telo_Eahy_pm_IL1b_TNF_VEGF.txt", stringsAsFactors=F)
TPM.df <- RNAseq_CITV %>% select(Gene_symbol, TeloHAEC.Ctrl_avg) %>% `colnames<-`(c("gene", "TeloHAEC_ctrl_TPM"))
narrow.df.TPM <- merge(TPM.df, narrow.df %>% select(-ENSGID, -EntrezID) %>% unique, by.x=c("gene"), by.y=c("Gene"), all.y=T)
narrow.df.TPM$gene[narrow.df.TPM$gene == "MESDC1"] <- "TLNRD1" ## one time fix MESDC1 -> TLNRD1



## ## debug genes without PoPS_Score (add the correct information from narrow.df)
## debug.CAD.GWAS.df <- CAD.GWAS.df %>% subset(!(CAD.GWAS.df$gene %in% narrow.df$Gene))
## debug.genes <- debug.CAD.GWAS.df %>% pull(gene) %>% unique %>% sort
## debug.genes.10X.format <- debug.genes %>% perturbseq.gene.names.to10X ## in 10X format
## debug.genes.perturbseq.format <- debug.genes %>% tenX.gene.names.toperturbseq ## in perturbseq.format
## debug.genes.all.formats <- c(debug.genes.10X.format, debug.genes.perturbseq.format) %>% unique %>% sort  ## 218 genes (with duplicates of names in different format)##### ## merge narrow.df into CAD.GWAS.df for debug genes
## debug.narrow.df <- narrow.df %>% subset(Gene %in% debug.genes.all.formats)

## debug.narrow.df.genes <- debug.narrow.df %>% pull(Gene) %>% unique
## debug.narrow.df.genes.10X.format <- debug.narrow.df.genes %>% perturbseq.gene.names.to10X
## debug.narrow.df.genes.perturbseq.format <- debug.narrow.df.genes %>% tenX.gene.names.toperturbseq
## debug.narrow.df.genes.all.formats <- c(debug.narrow.df.genes.10X.format, debug.narrow.df.genes.perturbseq.format) %>% unique ## all formats for genes that can be found in narrow.df (has cNMF analysis results)


## ##########################################################################################
## ## load cNMF results to get the list of input genes to cNMF
## k <- 60
## DENSITY.THRESHOLD <- gsub("\\.","_", "0.2")
## SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
## OUTDIRSAMPLE <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/"
## cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
## print(cNMF.result.file)
## if(file.exists(cNMF.result.file)) {
##     print("loading cNMF result file")
##     load(cNMF.result.file)
## }


########################################################################################
## Load in information about coding variants, and merge into main table
## To do: Move this into Rosa's script combine_tables_ubq_gene.R
all.cs.aragam <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/CAD_SNP_INFO/Aragam2021/all.cs.txt", stringsAsFactors=F)
all.cs.harst <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/CAD_SNP_INFO/Harst2017/all.cs.txt", stringsAsFactors=F)

addCodingGenes <- function(df, all.cs) {
    codingGenes <- all.cs %>% filter(AnyCoding) %>% select(CredibleSet, CodingVariantGene)
    df <- df %>% merge(codingGenes, all.x=TRUE) %>%
            rowwise() %>%
            mutate(GeneContainsCodingVariant=ifelse(is.na(CodingVariantGene), FALSE, (gene %in% strsplit(CodingVariantGene,";")[[1]]))) %>%
            ungroup() %>%
            as.data.frame()
    return(df)
}

CAD.GWAS.df <- CAD.GWAS.df %>% 
    ## mutate(Resource=ifelse(Resource=="Aragam_2021","Aragam2021","Harst2017"),
    ##        CredibleSet=paste0(Resource,"-",Lead_SNP_rsID)) %>%
    addCodingGenes(rbind(all.cs.aragam,all.cs.harst)) %>%
    mutate(CredibleSet = ifelse(CredibleSet %in% c("Aragam2021-rs768453105", "Harst2017-rs138120077"), "Aragam2021-rs768453105_Harst2017-rs138120077", CredibleSet))



########################################################################################
## Define key variables for analysis
CAD.GWAS.df <- CAD.GWAS.df %>% mutate(
    ExpressedInTeloHAEC = (TeloHAEC_ctrl_TPM >= 1),
    InPlausibleECLocus = ( !LipidLevelsAssociated & (!is.na(RankOfDistanceToECPeakWithVariant) | !is.na(MaxABC.Rank.ECOnly) | !is.na(CodingVariantGene)) ), ## Also require that locus is not associated with lipid levels
    InPlausibleECLocus_includeLipid = (!is.na(RankOfDistanceToECPeakWithVariant) | !is.na(MaxABC.Rank.ECOnly) | !is.na(CodingVariantGene)),
    InPlausibleLocus = ( !LipidLevelsAssociated & (!is.na(rank_SNP_to_TSS) | !is.na(MaxABC.Rank) | !is.na(CodingVariantGene)) ),
    InPlausibleLocus_includeLipid = !is.na(rank_SNP_to_TSS) | !is.na(MaxABC.Rank) | !is.na(CodingVariantGene),
    TopCandidate = (RankOfDistanceToECPeakWithVariant <= 2 | MaxABC.Rank.ECOnly <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE),
    TopCandidateInECLocus = TopCandidate & InPlausibleECLocus,
    ExpressedTopCandidateInECLocus = ((RankOfDistanceToECPeakWithVariant <= 2 | MaxABC.Rank.ECOnly <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE)) & InPlausibleECLocus & ExpressedInTeloHAEC,
    TopCandidateInECLocus_includeLipid = ((RankOfDistanceToECPeakWithVariant <= 2 | MaxABC.Rank.ECOnly <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE)) & InPlausibleECLocus_includeLipid
    ## TopCandidate = ((rank_SNP_to_TSS <= 2 | MaxABC.Rank <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE)) & InPlausibleLocus,
    ## TopCandidate_includeLipid = ((rank_SNP_to_TSS <= 2 | MaxABC.Rank <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE)) & InPlausibleLocus_includeLipid
    )

CAD.GWAS.df <- CAD.GWAS.df %>% 
    rowwise() %>%
    mutate(ProgramsLinkedToGene= c(strsplit(ProgramsRegulatedByThisGene,",")[[1]], strsplit(ProgramsInWhichGeneIsInTop300ZScoreSpecificGenes, ",")[[1]]) %>% unique() %>% paste(collapse=',')) %>%
    as.data.frame()

## merge CAD.GWAS.df into narrow.df.TPM to get all expressed gene's information
remove_na <- function(x) x[is.na(x)] <- NULL
narrow.df.TPM <- merge(narrow.df.TPM, # %>%
                       CAD.GWAS.df %>%
                       select(-one_of(colnames(narrow.df.TPM)[!(colnames(narrow.df.TPM) %in% c("gene"))])),
                       by=c("gene"), all.x=T
                       ) %>%
    rowwise() %>%
    ## mutate(ProgramsLinkedToGene= c(strsplit(ProgramsRegulatedByThisGene,",")[[1]], strsplit(ProgramsInWhichGeneIsInTop300ZScoreSpecificGenes, ",")[[1]]) %>% unique() %>% paste(collapse=',')) %>%
    as.data.frame ## because CAD.GWAS.df has dupicated genes, TopicsInWhichGeneIsInTop100ZScoreSpecificGenes sum up to > 100
narrow.df.TPM[is.na(narrow.df.TPM)] <- FALSE
narrow.df.TPM <- narrow.df.TPM %>%
    rowwise %>%
    mutate(AllGenesIncNMFInput = gene %in% rownames(theta)) %>%
    as.data.frame

## merge IBD.GWAS.df into narrow.df.TPM to get all TeloHAEC expressed gene's information
CAD.IBD.colnames.intersect <- intersect(CAD.GWAS.df %>% colnames, IBD.GWAS.df %>% colnames)
CAD.IBD.colnames.diff <- setdiff(CAD.GWAS.df %>% colnames, IBD.GWAS.df %>% colnames) %>% append("gene")
IBD.GWAS.df <- merge(narrow.df.TPM %>% select(all_of(CAD.IBD.colnames.diff)), IBD.GWAS.df %>% select(all_of(CAD.IBD.colnames.intersect)), by="gene", all.y=T) %>%
    mutate(ExpressedInTeloHAEC = (TeloHAEC_ctrl_TPM >= 1))
IBD.GWAS.df[is.na(IBD.GWAS.df)] <- FALSE

#######################################################################################
## Output some statistics


## function to combine signals when they are within 10kb of each other
group_loci <- function(df) {
    tmp <- df %>% select(SNP_position) %>% as.matrix %>% as.numeric
    threshold <- 1e5
    group <- rep(1, length(tmp))
    for (i in 1:(length(tmp) - 1)) {
        for (j in i:length(tmp)) {
            difference <- abs(tmp[i] - tmp[j])
            if ( difference < threshold ) {
                group[j] <- group[i]
            } else {
                group[j] <- group[j-1] + 1
            }
        }
    }
    group <- paste0(df$chr %>% unique, "_", group)
    return(group)
}

CAD.GWAS.df <- CAD.GWAS.df %>%
    group_by(chr) %>%
    mutate(group_loci = tibble(chr, SNP_position) %>% group_loci) %>%
    as.data.frame  %>%
    arrange(chr, SNP_position)

IBD.GWAS.df <- IBD.GWAS.df %>%
    group_by(chr) %>%
    mutate(group_loci = tibble(chr, SNP_position) %>% group_loci) %>%
    as.data.frame  %>%
    arrange(chr, SNP_position)


## Count of credible sets that count as "plausibleECLocus":
cat("Number of credible sets that count as 'in a plausible EC locus': ", 
    CAD.GWAS.df %>% filter(InPlausibleECLocus) %>% pull(CredibleSet) %>% unique() %>% length(), 
    " out of ",
    CAD.GWAS.df %>% pull(CredibleSet) %>% unique() %>% length(),
    "\n")

## Count of loci that count as "plausibleECLocus":
cat("Number of loci that count as 'in a plausible EC locus': ", 
    CAD.GWAS.df %>% filter(InPlausibleECLocus) %>% pull(group_loci) %>% unique() %>% length(), 
    " out of ",
    CAD.GWAS.df %>% pull(group_loci) %>% unique() %>% length(),
    "\n")


## Count of credible sets that count as "plausibleECLocus_includeLipid":
cat("Number of credible sets that count as 'in a plausible EC locus (including Lipid Loci)': ", 
    CAD.GWAS.df %>% filter(InPlausibleECLocus_includeLipid) %>% pull(CredibleSet) %>% unique() %>% length(), 
    " out of ",
    CAD.GWAS.df %>% pull(CredibleSet) %>% unique() %>% length(),
    "\n")

## Count of loci that count as "plausibleECLocus_includeLipid":
cat("Number of loci that count as 'in a plausible EC locus (including Lipid Loci)': ", 
    CAD.GWAS.df %>% filter(InPlausibleECLocus_includeLipid) %>% pull(group_loci) %>% unique() %>% length(), 
    " out of ",
    CAD.GWAS.df %>% pull(group_loci) %>% unique() %>% length(),
    "\n")

## Count of credible sets that count as "plausibleLocus":
cat("Number of credible sets that count as 'in a plausible locus': ", 
    CAD.GWAS.df %>% filter(InPlausibleLocus) %>% pull(CredibleSet) %>% unique() %>% length(), 
    " out of ",
    CAD.GWAS.df %>% pull(CredibleSet) %>% unique() %>% length(),
    "\n")

## Count of loci that count as "plausibleLocus":
cat("Number of loci that count as 'in a plausible locus': ", 
    CAD.GWAS.df %>% filter(InPlausibleLocus) %>% pull(group_loci) %>% unique() %>% length(), 
    " out of ",
    CAD.GWAS.df %>% pull(group_loci) %>% unique() %>% length(),
    "\n")

## Count of credible sets that count as "plausibleLocus_includeLipid":
cat("Number of credible sets that count as 'in a plausible locus (including Lipid Loci)': ", 
    CAD.GWAS.df %>% filter(InPlausibleLocus_includeLipid) %>% pull(CredibleSet) %>% unique() %>% length(), 
    " out of ",
    CAD.GWAS.df %>% pull(CredibleSet) %>% unique() %>% length(),
    "\n")

## Count of loci that count as "plausibleLocus_includeLipid":
cat("Number of loci that count as 'in a plausible locus (including Lipid Loci)': ", 
    CAD.GWAS.df %>% filter(InPlausibleLocus_includeLipid) %>% pull(group_loci) %>% unique() %>% length(), 
    " out of ",
    CAD.GWAS.df %>% pull(group_loci) %>% unique() %>% length(),
    "\n")


## Table of genes with coding variants and their associations with lipids or BP
CAD.GWAS.df %>% filter(GeneContainsCodingVariant) %>% 
    select(CredibleSet,CodingVariantGene,LipidLevelsAssociated,BloodPressureAssociated) %>% 
    unique() %>% 
    write.table(file=paste0(OUTDIR, "CodingVariantSummary.tsv"), sep='\t', quote=F, col.names=T, row.names=F)


#######################################################################################
## Combine CAD and IBD GWAS data frame
## ## find the genes in IBD GWAS loci but not in CAD GWAS loci
## IBD.CAD.diff.genes <- setdiff(IBD.GWAS.df$TargetGene %>% unique, CAD.GWAS.df$gene %>% unique)
## IBD.unique.genes.in.CAD.GWAS.df.format <- data.frame(matrix(ncol = ncol(CAD.GWAS.df), nrow=length(IBD.CAD.diff.genes)))
## colnames(IBD.unique.genes.in.CAD.GWAS.df.format) <- colnames(CAD.GWAS.df)
## IBD.unique.genes.in.CAD.GWAS.df.format$gene <- IBD.CAD.diff.genes
## ## need to also sync in columns such as SNP_position, Lead_SNP_rsID, CredibleSet, TSS_position, SNP_to_TSS_distance, etc
CAD.IBD.GWAS.df <- rbind(IBD.GWAS.df, CAD.GWAS.df) %>% unique


## save CAD.GWAS.df
write.table(CAD.GWAS.df, paste0(OUTDIR, "CAD.GWAS.df.txt"), sep="\t", quote=F, row.names=F)
write.table(CAD.GWAS.df %>% filter(TopCandidateInECLocus), paste0(OUTDIR, "CAD.GWAS.df.TopCandidateInECLocus.txt"), sep="\t", quote=F, row.names=F)
write.table(narrow.df.TPM, paste0(OUTDIR, "narrow.df.TPM.txt"), sep="\t", quote=F, row.names=F)
write.table(CAD.IBD.GWAS.df, paste0(OUTDIR, "CAD.IBD.GWAS.df.txt"), sep="\t", quote=F, row.names=F)

candidateCADGenes.df <- CAD.GWAS.df %>% filter(TopCandidateInECLocus) %>% unique()
write.table(candidateCADGenes.df, paste0(OUTDIR, "candidateCADGenes.df.txt"), sep="\t", quote=F, row.names=F)
write.table(candidateCADGenes.df %>% pull(gene) %>% unique %>% sort, paste0(OUTDIR, "candidateCADGenes.list.txt"), sep="\t", quote=F, row.names=F, col.names=F)

## ## check for missing genes ## not missing any genes after Rosa's update 220422
## missing.genes <- setdiff(CAD.GWAS.df.original$gene %>% unique, CAD.GWAS.df$gene %>% unique)
## CAD.GWAS.df.missing.genes.tocheck <- CAD.GWAS.df.original %>% subset(gene %in% missing.genes)
## write.table(CAD.GWAS.df.missing.genes.tocheck, paste0(OUTDIR, "CAD.GWAS.df.missing.genes.tocheck.txt"), sep="\t", quote=F, row.names=F)


## subset CAD.GWAS.df to only the necessary columns
CAD.GWAS.df.all.columns <- CAD.GWAS.df
CAD.GWAS.df <- CAD.GWAS.df.all.columns %>% select(gene, original.gene, TopCandidateInECLocus, on_2kG_lib, ProgramsInWhichGeneIsInTop300ZScoreSpecificGenes, ProgramsRegulatedByThisGene, TopCandidate, InPlausibleECLocus, LipidLevelsAssociated, rank_SNP_to_TSS, MaxABC.Rank, CodingVariantGene, gene_type) %>% unique


# ## output debug.genes that are not in narrow.df
# theta.zscore.rank.df <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/theta.zscore.rank.df.txt", stringsAsFactors=F) ## need theta.zscore.rank.df
# debug.genes.not.in.narrow.df <- setdiff(debug.genes.10X.format, debug.narrow.df.genes.10X.format) %>% perturbseq.gene.names.to10X %>% unique %>% tenX.gene.names.toperturbseq %>% unique %>%
#     as.data.frame %>%
#     `colnames<-`("Gene") %>%
#     mutate(GenesIncNMF = Gene %in% (theta.zscore.rank.df %>% pull(Gene) %>% unique),
#            IsProgramGene = Gene %in% (theta.zscore.rank.df %>% subset(zscore.specificity.rank <= 300) %>% pull(Gene))) %>%
#     merge(gtf.10X.df %>% select(gene_name, gene_type) %>% unique, by.x="Gene", by.y="gene_name", all.x=T) %>%
#     merge(CAD.GWAS.df, by.x="Gene", by.y="gene") %>%
#     unique
# write.table(debug.genes.not.in.narrow.df, paste0(OUTDIR, "genes.not.found.in.cNMF.analysis.results.txt"), sep="\t", row.names=F, quote=F)


## CAD.GWAS.df <- CAD.GWAS.df %>% subset(gene_type != "lncRNA")

## ## load gene annotations (Refseq + Uniprot)
## gene.summaries.path <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/combined.gene.summaries.txt"
## gene.summary <- read.delim(gene.summaries.path, stringsAsFactors=F)
## ## save TopCandidateGene list
## topic.list <- c(8, 15, 39, 48)
## topic.name.list <- paste0("K60_", topic.list)
## TopCandidateGene.df <- CAD.GWAS.df %>%
##     subset(ExpressedInTeloHAEC & TopCandidateInECLocus) %>%
##     select(gene, ProgramsInWhichGeneIsInTop300ZScoreSpecificGenes) %>%
##     rowwise %>%
##     subset(topic.name.list %in% (strsplit(ProgramsInWhichGeneIsInTop300ZScoreSpecificGenes, split="\\|") %>% unlist %>% as.character)) ## debug strsplit and %in% 
    
## TopCandidateGene <- CAD.GWAS.df %>%
##     subset(ExpressedInTeloHAEC & TopCandidateInECLocus) %>%
##     select(gene) %>%
##     unique %>%
##     merge(gene.summary, by.x="gene", by.y="Gene", all.x=T)
## write.table(TopCandidateGene, "CAD.GWAS.TopCandidateGene.txt", sep="\t", quote=F, row.names=F)


########################################################################################
## Run statistical Tests
## Create subset dfs based on background conditions
# k <- 60
k <- opt$K.val 

## read table that specifies which gene set is for alternative hypothesis and which set of topics to use
statistical.test.list.df <- read.delim(opt$statistical.test.list, stringsAsFactors=F)

## num.tests <- opt$num.tests ## actual number of tests requested by the user
## statistical.test.list.df <- statistical.test.list.df[1:num.tests,]

num.tests <- nrow(statistical.test.list.df) ## specify number of statistical tests to conduct
df.list.to.test <- vector("list", num.tests) ## initiate background data storage list
names(df.list.to.test) <- statistical.test.list.df$test.name

## remove batch topics from test
## load topic summary
# TopicSummary.df <- read.delim("../data/210730_cNMF_topic_model_analysis.xlsx - TopicCatalogNEW.tsv", stringsAsFactors=F)

# ## batch effect topic
# batch.topics <- TopicSummary.df %>%
#     filter(ProgramCategoryLabel == "Batch") %>%
#     pull(ProgramID) %>%
#     unique

# topics_to_test <- setdiff(c(1:k), batch.topics %>% gsub("K60_", "", .))

## do not remove batch topic from the test
topics_to_test = c(1:k)

all.test.results.background.list <- vector("list", num.tests)
for(i in 1:num.tests) {
    background.name <- names(df.list.to.test)[i] ## get the name of the background filter
    subset.df <- eval(parse(text = paste0("subset.df <- ", statistical.test.list.df$background.subset.command[i])))

    all.test.results.list <- vector("list", k) ## initialize variable
    ## loop over every topic t
    for(t in topics_to_test) {
        topic <- paste0("K60_", t)

        ## new scratch for test that only looks at 'expressed genes in a topic'
        genes.to.test <- statistical.test.list.df$genes.to.test[i]
        topics.to.test <- statistical.test.list.df$topics.to.test[i]

        list.to.test <- subset.df %>%
            rowwise %>%
            mutate(hypothesis.gene = ifelse(eval(parse(text = genes.to.test)), 1, 0),
                   LinkedToTopic.gene = ifelse(topic %in% (get(topics.to.test) %>% strsplit("\\|") %>% unlist %>% as.character), 1, 0)) %>% ## debug
            select(gene, hypothesis.gene, LinkedToTopic.gene) %>%
            unique %>%
            group_by(gene) %>%
            summarize(hypothesis.gene = ifelse(hypothesis.gene %>% sum > 0, "CandidateGene", "NotCandidateGene"),
                      LinkedToTopic.gene = ifelse(LinkedToTopic.gene %>% sum > 0, "LinkedToTopic", "NotLinkedToTopic")) %>%
            as.data.frame


        table.to.test <- table(list.to.test$hypothesis.gene, list.to.test$LinkedToTopic.gene) ## count frequencies
        if(ncol(table.to.test) == 1) {
            columns.needed <- c("LinkedToTopic", "NotLinkedToTopic")
            column.to.add <- columns.needed[!(columns.needed %in% colnames(table.to.test))]
            table.to.add <- matrix(c(0,0), nrow=2, dimnames=list(rownames(table.to.test),c(column.to.add)))
            if(column.to.add == "NotLinkedToTopic") table.to.test <- cbind(table.to.test, table.to.add) else table.to.test <- cbind(table.to.add, table.to.test)
        } else if (nrow(table.to.test) == 1) {
            rows.needed <- c("CandidateGene", "NotCandidateGene")
            row.to.add <- row.needed[!(row.needed %in% rownames(table.to.test))]
            table.to.add <- matrix(c(0,0), ncol=2, dimnames=list(c(row.to.add), colnames(table.to.test)))
            if(row.to.add == "NotCandidateGene") table.to.test <- rbind(table.to.test, table.to.add) else table.to.test <- rbind(table.to.add, table.to.test)
        }


        flattened.table.to.store <- table.to.test %>%
            matrix(nrow=1, ncol=4) %>% ## flatten the table
            as.data.frame %>%
            `colnames<-`(c("LinkedToTopic_CandidateGene",
                           "LinkedToTopic_NotCandidateGene",
                           "NotLinkedToTopic_CandidateGene",
                           "NotLinkedToTopic_NotCandidateGene")) %>% ## assign the table values to the corresponding column names
            mutate(enrichment = ( LinkedToTopic_CandidateGene / (LinkedToTopic_CandidateGene + NotLinkedToTopic_CandidateGene) ) / ( LinkedToTopic_NotCandidateGene / ( LinkedToTopic_NotCandidateGene + NotLinkedToTopic_NotCandidateGene) ), ## calculate enrichment
                   enrichment.log2fc = log2(enrichment), ## log2fc of the enrichment
                   Gene_LinkedToTopic_CandidateGene = list.to.test %>% ## list the genes that fall into each of the four categories in the contingency table
                       subset(hypothesis.gene == "CandidateGene" & LinkedToTopic.gene == "LinkedToTopic") %>%
                       pull(gene) %>%
                       paste0(collapse=","),
                   Gene_LinkedToTopic_NotCandidateGene = list.to.test %>%
                       subset(hypothesis.gene == "NotCandidateGene" & LinkedToTopic.gene == "LinkedToTopic") %>%
                       pull(gene) %>%
                       paste0(collapse=","),
                   Gene_NotLinkedToTopic_CandidateGene = list.to.test %>%
                       subset(hypothesis.gene == "CandidateGene" & LinkedToTopic.gene == "NotLinkedToTopic") %>%
                       pull(gene) %>%
                       paste0(collapse=",")#,
                   ## Gene_NotLinkedToTopic_NotCandidateGene = list.to.test %>%
                   ##     subset(hypothesis.gene == "NotCandidateGene" & LinkedToTopic.gene == "NotLinkedToTopic") %>%
                   ##     pull(gene) %>%
                   ##     paste0(collapse=",")
                   )

        ## binomial test
        ## example: given a group of genes (for null hypothesis), the percentage of these genes linked to a topic is p, plug in this probability into the following test.

        ## binomial test
        ## binom.test(x, n, p = 0.5)
        ## x: number of successes (number of orange genes linked to the topic)
        ## n: number of trials (total number of orange genes)
        ## p: hypothesized probability of sucess (from the null hypothesis)

        ## null distribution probability
        p.binomial.null <- (list.to.test %>% subset(LinkedToTopic.gene == "LinkedToTopic") %>% nrow) / (list.to.test %>% nrow)
        p.binomial.test <- binom.test(list.to.test %>% subset(hypothesis.gene == "CandidateGene" & LinkedToTopic.gene == "LinkedToTopic") %>% nrow, ## number of candidate genes that link to the topic
                                      list.to.test %>% subset(hypothesis.gene == "CandidateGene") %>% nrow, ## number of candidate genes
                                      p = p.binomial.null)$p.value ## binomial test

        if (dim(table.to.test)[1] == 2 & dim(table.to.test) [2] == 2) {
            fisher.p.val <- fisher.test(table.to.test, alternative="greater")$p.value
        } else {
            fisher.p.val = NA ## add a row of 0?
        }
        all.test.results.list[[t]] <- data.frame(Topic = topic,
                                                 fisher.p.value = fisher.p.val,
                                                 binomial.p.value = p.binomial.test) %>%
            cbind(flattened.table.to.store) ## put all results for one topic in one line to store
    }
    all.test.results.background.list[[i]] <- do.call(rbind, all.test.results.list) %>% mutate(fisher.p.adjust = p.adjust(fisher.p.value), binomial.p.adjust = p.adjust(binomial.p.value), background = background.name, .after="binomial.p.value")
}
all.test.results <- do.call(rbind, all.test.results.background.list) %>%
    as.data.frame %>%
    ## mutate(Gene_LinkedToTopic_NotCandidateGene = as.character(Gene_LinkedToTopic_NotCandidateGene)) %>%
    mutate(Gene_NotLinkedToTopic_CandidateGene = ifelse(NotLinkedToTopic_CandidateGene > 300, "", Gene_NotLinkedToTopic_CandidateGene),
           Gene_LinkedToTopic_CandidateGene = ifelse(LinkedToTopic_CandidateGene > 150, "", Gene_LinkedToTopic_CandidateGene)) #Gene_LinkedToTopic_NotCandidateGene))

all.test.results.nonbatch <- all.test.results %>%
    subset(Topic %in% paste0("K60_", setdiff(c(1:k), batch.topics %>% gsub("K60_", "", .)))) %>%
    group_by(background) %>%
    mutate(fisher.p.adjust = p.adjust(fisher.p.value, method="fdr"),
           binomial.p.adjust = p.adjust(binomial.p.value, method="fdr"),
           .after="binomial.p.value") %>%
    as.data.frame

## output statistical test results
write.table(all.test.results, file=paste0(OUTDIR, "program_prioritization.txt"), sep="\t", quote=F, row.names=F)
write.table(all.test.results.nonbatch, file=paste0(OUTDIR, "program_prioritization.nonbatch.txt"), sep="\t", quote=F, row.names=F)

if(opt$perturbseq) {
## combine regulator and program gene's p-value, then do multiple hypothesis correction ## regulators and program genes should be statistically independent
regulator.programGene.test.pairs <- data.frame(regulator = statistical.test.list.df$test.name[seq(2, num.tests, by=2)],
                                               programGene = statistical.test.list.df$test.name[seq(1, num.tests-1, by=2)],
                                               name = statistical.test.list.df %>% separate(col="test.name", sep="_geneSet.", into = c("regulator.or.program.gene", "toProcess")) %>% separate(col="toProcess", sep="_", into=c("test.handle", "background")) %>% pull(test.handle) %>% unique)
regulator.programGene.combined.pval.nonbatch.list <- vector("list", nrow(regulator.programGene.test.pairs))
for(i in 1:nrow(regulator.programGene.test.pairs)) {
    regulator.results.nonbatch <- all.test.results.nonbatch %>%
        subset(background == regulator.programGene.test.pairs$regulator[i]) %>%
        `colnames<-`(paste0("Regulator_", colnames(.))) %>%
        dplyr::rename("ProgramID" = "Regulator_Topic")
    programGene.results.nonbatch <- all.test.results.nonbatch %>% subset(background == regulator.programGene.test.pairs$programGene[i]) %>%
        `colnames<-`(paste0("ProgramGene_", colnames(.))) %>%
        dplyr::rename("ProgramID" = "ProgramGene_Topic")

    regulator.programGene.combined.pval.nonbatch.list[[i]] <- merge(regulator.results.nonbatch,
                                                                    programGene.results.nonbatch,
                                                                    by="ProgramID") %>%
        mutate(LinkedGenes_fisher.p.value = Regulator_fisher.p.value * ProgramGene_fisher.p.value,
               LinkedGenes_binomial.p.value = Regulator_binomial.p.value * ProgramGene_binomial.p.value,
               LinkedGenes_fisher.p.adjust = p.adjust(LinkedGenes_fisher.p.value, method="fdr"),
               LinkedGenes_binomial.p.adjust = p.adjust(LinkedGenes_binomial.p.value, method="fdr"),
               LinkedGenes_fisherNegLog10FDR = -log10(LinkedGenes_fisher.p.adjust),
               LinkedGenes_binomialNegLog10FDR = -log10(LinkedGenes_binomial.p.adjust),
               LinkedGenes_LinkedToTopic_CandidateGene = Regulator_LinkedToTopic_CandidateGene + ProgramGene_LinkedToTopic_CandidateGene,
               .after = "ProgramID") %>%
        arrange(LinkedGenes_binomial.p.adjust) %>% mutate(test.name = regulator.programGene.test.pairs$name[i], .after = "ProgramID")
}
regulator.programGene.combined.pval.nonbatch <- do.call(rbind, regulator.programGene.combined.pval.nonbatch.list) %>%
    arrange(LinkedGenes_fisher.p.adjust, desc(LinkedGenes_LinkedToTopic_CandidateGene))

write.table(regulator.programGene.combined.pval.nonbatch, file=paste0(OUTDIR, "combinedRegulatorProgramGeneEnrichmentTest.to.prioritize.topics.nonbatch.txt"), sep="\t", quote=F, row.names=F)

regulator.programGene.combined.pval.nonbatch.toSuppTable <- regulator.programGene.combined.pval.nonbatch %>%
    select(ProgramID,
           LinkedGenes_fisher.p.value,
           LinkedGenes_fisher.p.adjust,
           LinkedGenes_fisherNegLog10FDR,
           LinkedGenes_LinkedToTopic_CandidateGene,
           Regulator_fisher.p.value,
           Regulator_fisher.p.adjust,
           Regulator_LinkedToTopic_CandidateGene,
           Regulator_LinkedToTopic_NotCandidateGene,
           Regulator_NotLinkedToTopic_CandidateGene,
           Regulator_NotLinkedToTopic_NotCandidateGene,
           Regulator_enrichment,
           Regulator_enrichment.log2fc,
           Regulator_Gene_LinkedToTopic_CandidateGene,
           ProgramGene_fisher.p.value,
           ProgramGene_fisher.p.adjust,
           ProgramGene_LinkedToTopic_CandidateGene,
           ProgramGene_LinkedToTopic_NotCandidateGene,
           ProgramGene_NotLinkedToTopic_CandidateGene,
           ProgramGene_NotLinkedToTopic_NotCandidateGene,
           ProgramGene_enrichment,
           ProgramGene_Gene_LinkedToTopic_CandidateGene
           )
           
write.table(regulator.programGene.combined.pval.nonbatch.toSuppTable, file=paste0(OUTDIR, "SuppTable.combinedRegulatorProgramGeneEnrichmentTestToPrioritizeCADPrograms.nonbatch.txt"), sep="\t", quote=F, row.names=F)


}

## ## output CAD.GWAS.df that creates the program priorization table
## write.table(CAD.GWAS.df, "../SuppTables/CAD.GWAS.GenePrioritzationTable.txt", sep="\t", quote=F, row.names=F)


## ########################################################################################
## ## Plot statistical test results
## ## volcano plot ## x-axis: enrichment, y-axis: -log10(p-value)

## ## known CAD gene
## knownCADGene <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/known_CAD_gene_set.txt", stringsAsFactors=F) %>% as.matrix %>% as.character %>% sort

## toplot <- all.test.results %>%
##     select(-Gene_NotLinkedToTopic_CandidateGene, -Gene_LinkedToTopic_NotCandidateGene) %>%
##     rowwise %>% 
##     mutate(KnownCADGene = knownCADGene %in% strsplit(Gene_LinkedToTopic_CandidateGene, split=",")[[1]] %>% as.numeric %>% sum,
##            NovelCADGene = LinkedToTopic_CandidateGene - KnownCADGene) %>%
##     melt(id.vars=colnames(.)[!colnames(.) %in% c("KnownCADGene", "NovelCADGene")], value.name = "LinkedToTopic_CandidateGene_count", variable.name = "CandidateGeneType") %>% ## keep original columns the same and melt the new columns
##     as.data.frame
## write.table(toplot, file="all.results.table.to.plot.txt", quote=F, row.names=F, sep="\t")


## ## bar plot ## rank the topics by p-value and by the number of LinkedToTopic_CosestToECPeak gene, color by if the gene is a known CAD gene
## pdf("figures/statistical.test.result.barplot.pdf")
## for ( test.type in names(df.list.to.test) ) {
##     toplot.here <- toplot %>% subset(background == test.type)
##     topics.to.keep <- toplot.here %>% select(Topic, LinkedToTopic_CandidateGene) %>% arrange(desc(LinkedToTopic_CandidateGene)) %>% unique %>% slice(1:10) %>% pull(Topic) %>% as.character
##     toplot.here <- toplot.here %>% subset(Topic %in% topics.to.keep)
##     toplot.here.pval <- toplot.here %>% select(Topic, fisher.p.adjust, binomial.p.adjust) %>% unique %>% mutate(CandidateGeneType="-log10 adjusted p-value")

##     p <- toplot.here %>% ggplot(aes(x = reorder(Topic, LinkedToTopic_CandidateGene), y = LinkedToTopic_CandidateGene_count, fill = CandidateGeneType), color="black") + geom_col(width=0.7) + geom_col(data=toplot.here.pval, aes(x=Topic, y=-log10(binomial.p.adjust)), width=0.3) + mytheme +
##         geom_hline(yintercept=1, linetype="dotted", color = "red") +
##         coord_flip() +
##         xlab("Topic") + ylab("Gene Count") + ggtitle(test.type)
##     print(p)
## }
## dev.off()


## ## plot to compare different statistitcal hypothesis
## pdf(paste0(FIGDIR, "compare.statistical.test.pdf"), width=10, height=10)
## for (t in topics_to_test) {
##     topic <- paste0("topic_", t)
##     ## x-axis: test.name
##     ## y-axis: -log10 p.adjust
##     test.result.topic <- all.test.results %>% subset(Topic == topic)
##     toplot <- test.result.topic %>%
##         select(background, binomial.p.adjust, enrichment) %>%
##         mutate(negative.log10.binomial.p.adjust = -log10(binomial.p.adjust)) %>%
##         select(-binomial.p.adjust) %>%
##         melt(id.vars="background", value.name = "value", variable.name = "data.type") %>%
##         as.data.frame
##     p <- toplot %>% ggplot(aes(x=background, y=value, fill = data.type)) + geom_bar(stat="identity", position = "dodge", width=0.5) + mytheme +
##         geom_hline(yintercept = 1, linetype = "dotted") + 
##         xlab("Statistical Test") + ylab("") + ggtitle(paste0("Topic ", t)) +
##         theme(axis.text.x = element_text(angle = 90, hjust=1))
##     print(p)
## }
## dev.off()






## develop null hypothesis
## example: given a bag of 200 orange marbles and 800 blue margles, if you pick one at a time and repeat for 200 times, what is the chance that the marbles you picked are orange? (doesn't sound very right to Helen)
## another example: given a group of genes (for null hypothesis), the percentage of these genes linked to a topic is p, plug in this probability into the following test.

## binomial test
## binom.test(x, n, p = 0.5)
## x: number of successes (number of orange genes linked to the topic)
## n: number of trials (total number of orange genes)
## p: hypothesized probability of sucess (from the null hypothesis)


## save the results


## Possible background models:
## 1.  All expressed genes in TeloHAECs (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to a random expressed gene in TeloHAEC?)
## 2.  Other genes in EC CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 3.  Other genes in CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 4.  Other expressed genes in EC CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 5.  Other expressed genes in CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 6.  Other expressed, 2000Kperturbed genes in EC CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)
## 7.  Other expressed, 2000Kperturbed genes in CAD GWAS loci.  (Are the closest genes to EC CAD GWAS loci enriched for being in Topic 15, compared to other genes in those loci?)





## #########################################################
## ## A few one-off tests:

## ## Are the nearby EC genes enriched for regulating ANY topic?
## CAD.GWAS.df <- CAD.GWAS.df %>% 
##     rowwise() %>%
##     mutate(
##         RegulatesAnyTopic=grepl("opic",TopicsRegulatedByThisGene),
##         ExpressedInAnyTopic=grepl("opic",TopicsInWhichGeneIsInTop300ZScoreSpecificGenes),
##         LinkedToAnyTopic=grepl("opic",TopicsLinkedToGene),
##         LinkedToAnyEnrichedTopic=any(paste0("topic_",c(8,39,36,15,51,47,48,29,31,2,25,35)) %in% strsplit(TopicsLinkedToGene,",")[[1]]),
##         LinkedToAnyCCMTopic=any(paste0("topic_",c(8,39,15,48,35)) %in% strsplit(TopicsLinkedToGene,",")[[1]]),
##         LinkedToPrioritizedTopic=any(paste0("topic_",c(8,39,15,48)) %in% strsplit(TopicsLinkedToGene,",")[[1]])
##     ) %>%
##     ungroup() %>% 
##     as.data.frame()

## for (col in c("RegulatesAnyTopic","ExpressedInAnyTopic","LinkedToAnyTopic","LinkedToAnyEnrichedTopic","LinkedToAnyCCMTopic","LinkedToPrioritizedTopic")) {
##     print(col)
##     testTable <- with(CAD.GWAS.df %>% filter(ExpressedInTeloHAEC), 
##             c(
##                 gene[TopCandidateInECLocus & get(col)] %>% unique() %>% length(),
##                 gene[!TopCandidateInECLocus & get(col)] %>% unique() %>% length(),
##                 gene[TopCandidateInECLocus & !get(col)] %>% unique() %>% length(),
##                 gene[!TopCandidateInECLocus & !get(col)] %>% unique() %>% length()
##             )) %>%
##         matrix(nrow=2, ncol=2)
##     print(testTable)
##     fisher.test(testTable) %>% print()
## }


## ## scratch 220304 Helen
## CAD.GWAS.df %>%
##     select(gene, LinkedToPrioritizedTopic, ExpressedInTeloHAEC, TopCandidateInECLocus) %>%
##     subset(TopCandidateInECLocus & ExpressedInTeloHAEC) %>%
##     unique %>%
##     group_by(LinkedToPrioritizedTopic) %>%
##     summarize(numGenes = n())



## ################################################################################
## ## scratch 220107 Helen

## ########################################################################################
## ## Run Permutation Test
## ## permute rank_SNP_to_TSS within each GWAS loci
## permuted.rank_SNP_to_TSS.CAD.GWAS.df <- CAD.GWAS.df %>% group_by(CredibleSet) %>% mutate(rank_SNP_to_TSS = sample(rank_SNP_to_TSS)) %>% as.data.frame
## t <- 15
## topic <- paste0("topic_", t)
## ## CAD.GWAS.df %>% mutate(subset.topic = topic %in% (TopicsLinkedToGene %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(subset.topic)
## original.totest <- CAD.GWAS.df %>% rowwise %>% mutate(find.topic = topic %in% (TopicsLinkedToGene %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(find.topic) %>% as.data.frame
## permuted.null.totest <- permuted.rank_SNP_to_TSS.CAD.GWAS.df %>% rowwise %>% mutate(find.topic = topic %in% (TopicsLinkedToGene %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(find.topic) %>% as.data.frame

## wilcox.test(permuted.null.totest$rank_SNP_to_TSS, original.totest$rank_SNP_to_TSS)


## ## find how the genes linked to a topic are related to the GWAS loci
## ## Replot number of genes linked to topic T versus distance ranking histogram
## pdf(file=paste0(FIGDIR, "number.of.TopicsLinkedToGene.vs.distance.rank.pdf"))
## for (t in topics_to_test) {
##     topic <- paste0("topic_", t)
##     ## CAD.GWAS.df %>% mutate(subset.topic = topic %in% (TopicsLinkedToGene %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(subset.topic)
##     toplot <- CAD.GWAS.df %>% rowwise %>% mutate(find.topic = topic %in% (TopicsLinkedToGene %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(find.topic) %>% as.data.frame
##     p <- toplot %>% ggplot(aes(x=rank_SNP_to_TSS)) + geom_histogram(bins=max(toplot$rank_SNP_to_TSS)) + mytheme +
##         ggtitle(paste0("Topic ", t))
##     print(p)
## }
## dev.off()

## ##########################################################################################
## ## plot enrichment
## pdf("rank_SNP_to_TSS.enrichment.pdf")
## ## Plot observed genes at rank X over expected genes at rank X (enrichment) for each topic
## enrichment.short.df.list <- vector("list", length(topics_to_test))
## gene.choice.criteria.count.list <- c(100, 300, 500)
## gene.choice.criteria.count <- 100 # iterate on this
## gene.choice.critera <- paste0("TopicsInWhichGeneIsInTop", gene.choice.criteria.count, "ZScoreSpecificGenes")
## for (i in 1:length(topics_to_test)) {
##     t <- topics_to_test[i]
##     topic <- paste0("topic_", t)
##     observed.df <- CAD.GWAS.df %>%
##         select(gene, CredibleSet, rank_SNP_to_TSS, all_of(gene.choice.criteria)) %>%
##         unique %>%
##         rowwise %>%
##         mutate(find.topic = topic %in% (get(gene.choice.criteria) %>% tolower %>% strsplit(",") %>% unlist)) %>%
##         filter(find.topic) %>%
##         pull(rank_SNP_to_TSS) %>%
##         table %>%
##         as.data.frame %>%
##         `colnames<-`(c("rank_SNP_to_TSS", "count"))
##     expected.df <- CAD.GWAS.df %>%
##         select(gene, CredibleSet, rank_SNP_to_TSS) %>%
##         unique %>%
##         pull(rank_SNP_to_TSS) %>%
##         table %>%
##         as.data.frame %>%
##         `colnames<-`(c("rank_SNP_to_TSS", "total.genes.at.rank")) %>%
##         mutate(expected.count = total.genes.at.rank * gene.choice.criteria.count / nrow(theta)) ## adjust expected count based on gene selection criteria
##     enrichment.df <- merge(observed.df, expected.df, by="rank_SNP_to_TSS", all=T) %>%
##         mutate(count = replace_na(count, replace=0)) %>%
##         mutate(enrichment = count / expected.count) %>%
##         arrange(rank_SNP_to_TSS)
##     enrichment.df$rank_SNP_to_TSS <- enrichment.df$rank_SNP_to_TSS %>% as.character %>% as.numeric
##     enrichment.df <- enrichment.df %>%
##         mutate(gene.selection.criteria = gene.choice.criteria)
##     enrichment.short.df <- rbind(enrichment.df[which(enrichment.df$rank_SNP_to_TSS < 5),c("rank_SNP_to_TSS", "count", "expected.count", "total.genes.at.rank")], enrichment.df[which(enrichment.df$rank_SNP_to_TSS >= 5),c("rank_SNP_to_TSS", "count", "expected.count", "total.genes.at.rank")] %>% apply(2, sum)) %>%
##         mutate(enrichment = count / expected.count,
##                gene.selection.criteria = gene.choice.criteria)
##     enrichment.short.df$rank_SNP_to_TSS[enrichment.short.df$rank_SNP_to_TSS > 5] <- "5+"

##     ## ## binomial test on enrichment
##     ## p.binomial.null <- enrichment.short.df$expected.count[rank.index] / nrow(theta) ## expected number of genes linked to topic (or should we use number of genes at rank?) / number of genes in the input of cNMF pipeline
##     ## p.binomial.test <- binom.test(enrichment.short.df$count[rank.index], ## number of genes linked to topic,
##     ##                               gene.choice.criteria.count, ## gene selection criteria (e.g. 100, 300, 500)
##     ##                               p = p.binomial.null)$p.value
##     ## ### Example
##     ## ## p.binomial.null <- (list.to.test %>% subset(LinkedToTopic.gene == "LinkedToTopic") %>% nrow) / (list.to.test %>% nrow)
##     ## ## p.binomial.test <- binom.test(list.to.test %>% subset(hypothesis.gene == "CandidateGene" & LinkedToTopic.gene == "LinkedToTopic") %>% nrow, ## number of candidate genes that link to the topic
##     ## ##                               list.to.test %>% subset(hypothesis.gene == "CandidateGene") %>% nrow, ## number of candidate genes
##     ## ##                               p = p.binomial.null)$p.value ## binomial test

##     ## plot
##     p.enrichment <- enrichment.short.df %>% ggplot(aes(x=rank_SNP_to_TSS, y=enrichment)) + geom_col(fill="#38b4f7") + LARGEFONT +
##         xlab("Rank of Distance from SNP to TSS") + ylab("Enrichment")
##     p.count <- enrichment.short.df %>%
##         ggplot(aes(x=rank_SNP_to_TSS, y=count)) + geom_col(position="dodge", width=0.5) + LARGEFONT +
##         xlab("Rank of Distance form SNP to TSS") + ylab("Count")
##     p <- ggarrange(p.count, p.enrichment, nrow=1)
##     annotate_figure(p, top = paste0("Topic ", t, ", Enrichment for Top ", gene.choice.criteria.count, " Topic Defining Genes"), face = "bold")
##     enrichment.short.df.list[[i]] <- enrichment.short.df %>% mutate(Topic = topic)
## }
## dev.off()
## enrichment.all.df <- do.call(rbind, enrichment.short.df.list)
## enrichment.df.toheatmap <- enrichment.all.df %>%
##     select(rank_SNP_to_TSS, enrichment, Topic) %>%
##     spread(key = Topic, value = enrichment, fill=0) %>%
##     `rownames<-`(.$rank_SNP_to_TSS %>% as.character) %>%
##     select(-rank_SNP_to_TSS) %>%
##     as.matrix

## ##########################################################################################
## ## heatmap of rank_SNP_to_TSS enrichment
## ## heatmap function
## palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)
## plotHeatmap <- function(mtx, clusterRow = T, clusterCol = T, cexRow=2.5/(nrow(mtx)^(1/3)), cexCol=1.7/(ncol(mtx)^(1/3)), title="", xlab="", ylab=""){
##     heatmap.2(
##     mtx, 
##     Rowv=clusterRow, 
##     Colv=clusterCol,
##     trace='none',
##     key=T,
##     col=palette,
##     labCol=colnames(mtx),
##     margins=c(5,5), 
##     cex.main=0.1, 
##     cexRow=cexRow,
##     cexCol=cexCol,
##     main=title,
##     xlab=xlab,
##     ylab=ylab
## )
## }
## ## output pdf
## pdf("rank_SNP_to_TSS.enrichment.in.topic.heatmap.pdf")
## plotHeatmap(enrichment.df.toheatmap %>% t(), clusterCol=F, xlab="Distance Rank from SNP to TSS", ylab="Topic")
## dev.off()


## ## with permutation test
## permuted.rank_SNP_to_TSS.CAD.GWAS.df <- CAD.GWAS.df %>% group_by(CredibleSet) %>% mutate(rank_SNP_to_TSS = sample(rank_SNP_to_TSS)) %>% as.data.frame
## pdf(file=paste0(FIGDIR, "number.of.TopicsLinkedToGene.vs.distance.rank.with.permutation.test.pdf"))
## for (t in topics_to_test) {
##     topic <- paste0("topic_", t)
##     gene.choice.criteria <- "TopicsRegulatedByThisGene"
##     ## CAD.GWAS.df %>% mutate(subset.topic = topic %in% (TopicsLinkedToGene %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(subset.topic)
##     original.totest <- CAD.GWAS.df %>% rowwise %>% mutate(find.topic = topic %in% (get(gene.choice.criteria) %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(find.topic) %>% as.data.frame
##     permuted.null.totest <- permuted.rank_SNP_to_TSS.CAD.GWAS.df %>% rowwise %>% mutate(find.topic = topic %in% (get(gene.choice.criteria) %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(find.topic) %>% as.data.frame
##     p.original <- original.totest %>% ggplot(aes(x=rank_SNP_to_TSS)) + geom_histogram(bins=max(original.totest$rank_SNP_to_TSS)) + mytheme +
##         ggtitle(paste0("Original Data"))
##     p.permutation <- permuted.null.totest %>% ggplot(aes(x=rank_SNP_to_TSS)) + geom_histogram(bins=max(permuted.null.totest$rank_SNP_to_TSS)) + mytheme +
##         ggtitle(paste0("Permuted null distribution"))
##     p.value <- wilcox.test(permuted.null.totest$rank_SNP_to_TSS, original.totest$rank_SNP_to_TSS)$p.val
##     p <- ggarrange(p.original, p.permutation, nrow=1)
##     p <- annotate_figure(p, top = text_grob(paste0("Topic ", t, ", p-value: ", p.value)))
##     print(p)
## }
## dev.off()


## ## notes for permutation test: ## referred to http://www.stat.uchicago.edu/~yibi/teaching/stat222/2017/Lectures/nonpar.pdf
## ## another useful notes on one sample test: http://users.stat.umn.edu/~helwig/notes/perm-Notes.pdf
## ## 1. sample permuted.null.totest, find the mean difference between permuted.null.totest$rank_SNP_to_TSS and the original
## ## 2. repeat for M times, with M a large number (e.g. 1000)
## ## 3. count the number k of repetitions that produce mean at least as extreme (here, close to 0 and smaller than our original mean) as the mean difference of the original
## ## 4. when M is large enough, k/M should approximate p-value
## ## bonus: extrapolate k/M by varying M.

## nperm <- 10000
## mean.perm <- matrix(1, nrow=length(topics_to_test), ncol=nperm) ## initialize matrix to store the mean of permutated distance_to_top_SNP
## p.val.toExtrapolate <- matrix(NA, nrow=length(topics_to_test), ncol = log10(nperm)) ## want to extrapolate at nperm = {10, 100, 1000, 10000}
## for (topic.index in 2:length(topics_to_test)) {
##     t <- topics_to_test[topic.index]
##     topic <- paste0("topic_", t)
##     print(topic) ## report topic status
##     original.totest <- CAD.GWAS.df %>% rowwise %>% mutate(find.topic = topic %in% (TopicsLinkedToGene %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(find.topic) %>% as.data.frame
##     for (n in 1:nperm) {
##         if(n %% 100 == 0) print(n) ## report the status of the iterations
##         permuted.rank_SNP_to_TSS.CAD.GWAS.df <- CAD.GWAS.df %>% group_by(CredibleSet) %>% mutate(rank_SNP_to_TSS = sample(rank_SNP_to_TSS)) %>% as.data.frame ## permutation
##         permuted.null.totest <- permuted.rank_SNP_to_TSS.CAD.GWAS.df %>% rowwise %>% mutate(find.topic = topic %in% (TopicsLinkedToGene %>% tolower %>% strsplit(",") %>% unlist)) %>% filter(find.topic) %>% as.data.frame ## subset to genes linked to the topic
##         ## p.value.perm[n] <- wilcox.test(permuted.null.totest$rank_SNP_to_TSS, original.totest$rank_SNP_to_TSS)$p.val
##         mean.perm[topic.index,n] <- mean(permuted.null.totest$rank_SNP_to_TSS) ## mean of the permutation result
##     }
##     ## count number of permutations close to zero
##     max.exponent <- log10(nperm)
##     for (exponent in 1:max.exponent) {
##         p.val.toExtrapolate[topic.index, exponent] <- sum(as.numeric(mean.perm[1:(10^exponent)] < mean(original.totest$rank_SNP_to_TSS))) / (10^exponent) ## at each exponent value, calculate a p-value (to be extrapolated)
##     }
## } ## this took 8 hours?
## ## format the matrix
## rownames(p.val.toExtrapolate) <- paste0("topic_", topics_to_test)
## colnames(p.val.toExtrapolate) <- 10^(1:log10(nperm))

## ## save the data
## saveRDS(p.val.toExtrapolate, file="p.val.scratch.RDS")

## ## plot permutation p-value results
## ## plot extrapolation of p-value for each topic
## for (t in topics_to_test) {
##     topic <- paste0("topic_", t)
##     toplot <- p.val.toExtrapolate %>% as.data.frame %>% subset(rownames(.) %in% topic) %>% melt(value.name = "p.value", variable.name = "num.iterations")
## }

## ## plot p-value for all topics
## pdf(file=paste0(FIGDIR, "permutation.test.rank_SNP_to_TSS.nperm10000.pdf"))
## toplot <- p.val.toExtrapolate %>% as.data.frame %>% mutate(Topic = rownames(.))
## p <- toplot %>% ggplot(aes(x=reorder(Topic, desc(`10000`)), y=`10000`)) + geom_col() + mytheme + coord_flip() +
##     xlab("Topic") + ylab("p-value")
## print(p)
## dev.off()


## ## debug difference between CAD.GWAS.df and narrow.df.TPM
## CAD.GWAS.genes.not.in.narrow.df.TPM <- setdiff(CAD.GWAS.df %>% pull(gene) %>% unique, narrow.df.TPM.cNMFInput %>% pull(gene) %>% unique)
## CAD.GWAS.df %>% subset(gene %in% CAD.GWAS.genes.not.in.narrow.df.TPM) %>% pull(TopCandidateInECLocus) %>% table
## duplicated.narrow.df.TPM.EntrezID <- narrow.df.TPM %>% select(gene, EntrezID, ENSGID) %>% filter(duplicated(EntrezID) & EntrezID!="NULL") %>% pull(EntrezID) %>% unique
## duplicated.narrow.df.TPM.ENSGID <- narrow.df.TPM %>% select(gene, EntrezID, ENSGID) %>% filter(duplicated(ENSGID) & ENSGID!="NULL") %>% pull(ENSGID) %>% unique
## narrow.df.TPM %>% subset(EntrezID %in% duplicated.EntrezID) %>% select(gene, EntrezID, ENSGID)
## narrow.df.TPM %>% subset(ENSGID %in% duplicated.ENSGID) %>% select(gene, EntrezID, ENSGID)

## duplicated.narrow.df.EntrezID <- narrow.df %>% select(Gene, EntrezID, ENSGID) %>% filter(duplicated(EntrezID) & EntrezID!="NULL") %>% pull(EntrezID) %>% unique
## duplicated.narrow.df.ENSGID <- narrow.df %>% select(Gene, EntrezID, ENSGID) %>% filter(duplicated(ENSGID) & ENSGID!="NULL") %>% pull(ENSGID) %>% unique
## narrow.df %>% subset(EntrezID %in% duplicated.narrow.df.EntrezID) %>% select(Gene, EntrezID, ENSGID)
## narrow.df %>% subset(ENSGID %in% duplicated.narrow.df.ENSGID) %>% select(Gene, EntrezID, ENSGID)


## narrow.df %>% group_by(Gene) %>% mutate(num.rows = n()) %>% subset(num.rows > 1) %>% as.data.frame %>% write.table("scratch.narrow.df.duplicated.txt", row.names=F, quote=F, sep="\t")


