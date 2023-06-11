## Helen Kang
## Script to do program enrichment for any GWAS traits
## 230119

## setwd("/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/")

## getwd()
## statistical.test.dir <- "workflow/scripts/program_prioritization/All_GWAS_traits.statistical.test.list.txt"
## print(statistical.test.dir)
## file.exists(statistical.test.dir)
## statistical.test.list.df <- read.delim(paste0("workflow/scripts/program_prioritization/All_GWAS_traits.statistical.test.list.txt"), stringsAsFactors=F)

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
    make_option("--input.GWAS.table", type="character", default="/oak/stanford/groups/engreitz/Users/rosaxma/2111_pipeline_output/UKB/overlap/MAP/MAP_GWAS_gene_incl_ubq_genes.txt", help="Input CAD GWAS Table"),
    make_option("--cNMF.table", type="character", default="", help="Table with cNMF program result"),
    make_option("--outdir", type="character", default="./outputs/MAP/", help="Output directory"),
    make_option("--figdir", type="character", default="./figures/", help="Output directory"),
    make_option("--sampleName", type="character", default="2kG.library", help="Name of the sample"),
    make_option("--celltype", type="character", default="EC", help="Cell type in GWAS table"),
    make_option("--K.val", type="numeric", default=60, help="The value of K"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--outdirsample", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/", help="path to cNMF analysis results"), ## or for 2n1.99x: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/Perturb_2kG_dup4/K60/threshold_0_2/"
    make_option("--num.tests", type="numeric", default=2, help="number of statistical test to do from the top of statistical.test.df.txt"),
    make_option("--trait.name", type="character", default="MAP", help="name of the trait"),
    make_option("--coding.variant.df", type="character", default="/oak/stanford/groups/engreitz/Users/rosaxma/2111_pipeline_output/UKB/SBP/SBPvariant.list.1.coordinate.txt", help="Data frame with coding variant information for GWAS trait"),
    make_option("--regulator.analysis.type", type="character", default="GWASWide", help="path to statistical test recipe"),
    make_option("--perturbSeq", type="logical", default=FALSE, help="Whether this is a Perturb-seq experiment"),
    make_option("--TPM.table", type="character", default="", help="Path to TPM table for the correct cell type")
)
opt <- parse_args(OptionParser(option_list=option.list))

## ## sdev for K562 gwps 2k overdispersed genes K=80
## opt$input.GWAS.table <- "/oak/stanford/groups/engreitz/Users/rosaxma/2111_pipeline_output/UKB/overlap/RBC/RBC_GWAS_gene_incl_ubq_genes.txt"
## opt$cNMF.table <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K80/threshold_0_2/prepare_compute_enrichment.txt"
## opt$sampleName <- "WeissmanK562gwps"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/WeissmanK562gwps/K80/program_prioritization_GenomeWide/RBC/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K80/threshold_0_2/program_prioritization_GenomeWide/RBC/"
## opt$trait.name <- "RBC"
## opt$celltype <- "K562"
## opt$K.val <- 80
## opt$coding.variant.df <- "/oak/stanford/groups/engreitz/Users/rosaxma/2111_pipeline_output/UKB/RBC/RBCvariant.list.1.coordinate.txt"
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K80/threshold_0_2/"
## opt$perturbSeq <- TRUE
## opt$regulator.analysis.type <- "GenomeWide"



## ## sdev for K562 gwps 2k overdispersed genes K=90 RBC
## opt$input.GWAS.table <- "/oak/stanford/groups/engreitz/Users/rosaxma/2111_pipeline_output/UKB/overlap/RBC/RBC_GWAS_gene_incl_ubq_genes.txt"
## opt$cNMF.table <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/prepare_compute_enrichment.txt"
## opt$sampleName <- "WeissmanK562gwps"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/WeissmanK562gwps/K90/program_prioritization_GenomeWide/RBC/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/program_prioritization_GenomeWide/RBC/"
## opt$trait.name <- "RBC"
## opt$celltype <- "K562"
## opt$K.val <- 90
## opt$coding.variant.df <- "/oak/stanford/groups/engreitz/Users/rosaxma/2111_pipeline_output/UKB/RBC/RBCvariant.list.1.coordinate.txt"
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/"
## opt$perturbSeq <- TRUE
## opt$regulator.analysis.type <- "GenomeWide"
## opt$TPM.table <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/230512_TPM/outputs/K562.TPM.df.txt"


## ## sdev for K562 gwps 2k overdispersed genes K=90 Plt
## opt$input.GWAS.table <- "/oak/stanford/groups/engreitz/Users/rosaxma/2111_pipeline_output/UKB/overlap/Plt/Plt_GWAS_gene_incl_ubq_genes.txt"
## opt$cNMF.table <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/prepare_compute_enrichment.txt"
## opt$sampleName <- "WeissmanK562gwps"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/WeissmanK562gwps/K90/program_prioritization_GenomeWide/Plt/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/program_prioritization_GenomeWide/Plt/"
## opt$trait.name <- "Plt"
## opt$celltype <- "K562"
## opt$K.val <- 90
## opt$coding.variant.df <- "/oak/stanford/groups/engreitz/Users/rosaxma/2111_pipeline_output/UKB/Plt/Pltvariant.list.1.coordinate.txt"
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/"
## opt$perturbSeq <- TRUE
## opt$regulator.analysis.type <- "GenomeWide"
## opt$TPM.table <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/230512_TPM/outputs/K562.TPM.df.txt"


## ## sdev for EC Perturb-seq K=60
## opt$input.GWAS.table <- "/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/GWAS_tables/CAD_GWAS_gene_incl_ubq_genes_aragam_all_harst.txt"
## opt$cNMF.table <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/prepare_compute_enrichment.txt"
## opt$sampleName <- "2kG.library"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/2kG.library/K60/program_prioritization_GWASWide/CAD/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/program_prioritization_GWASWide/CAD/"
## opt$trait.name <- "CAD"
## opt$celltype <- "EC"
## opt$K.val <- 60
## opt$coding.variant.df <- ""
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/"
## opt$perturbSeq <- TRUE
## opt$regulator.analysis.type <- "GWASWide"



## Directories
FIGDIR <- opt$figdir
OUTDIR <- opt$outdir
DATADIR <- "/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/" ## update this path
check.dir <- c(FIGDIR, OUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))


##########################################################################################
## Load Data
## load cNMF results to get the list of input genes to cNMF (need theta)
k <- opt$K.val
SAMPLE <- opt$sampleName
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
OUTDIRSAMPLE <- opt$outdirsample
## OUTDIRSAMPLE <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/"
celltype <- opt$celltype

cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
}


## load 10X reference
gtf.10X.df <- readRDS(paste0("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/refdata-cellranger-arc-GRCh38-2020-A_genes.gtf_df.RDS")) ## load 10X gtf file
Gene.ENSEMBL.10X.df <- gtf.10X.df %>% 
    mutate(ENSGID = gene_id,
           Gene10X = gene_name) %>%
    select(Gene10X, ENSGID) %>%
    unique
## gtf <- importGTF("/home/groups/engreitz/Software/cellranger-arc-1.0.1/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf")

db <- ifelse(grepl("mouse", SAMPLE), "org.Mm.eg.db", "org.Hs.eg.db")
library(!!db) ## load the appropriate database
## helper function to map between ENSGID and SYMBOL
map.ENSGID.SYMBOL <- function(df) {
    ## need column `Gene` to be present in df
    ## detect gene data type (e.g. ENSGID, Entrez Symbol)
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
    df <- df %>%
        mutate(Gene = ifelse(is.na(Gene), "NA", Gene),
               ENSGID = ifelse(is.na(ENSGID), "NA", ENSGID)) 
    df <- df %>% merge(Gene.ENSEMBL.10X.df, by="ENSGID", all.x=T)
    df <- df %>% mutate(Gene = ifelse(Gene == "NA", Gene10X, Gene))
    notMatched.df <- df %>% subset(is.na(Gene10X))
    if(nrow(notMatched.df) > 0) {
        notMatched.index <- notMatched.df %>% rownames
        match.df <- merge(notMatched.df %>% select(-ENSGID, -Gene10X), Gene.ENSEMBL.10X.df, by.x="Gene", by.y="Gene10X") %>% mutate(Gene10X = Gene) %>% select(all_of(df %>% colnames))
        matched.index <- notMatched.df %>% subset(Gene %in% c(match.df$Gene %>% unique)) %>% rownames %>% as.numeric
        df <- rbind(df[-matched.index,], match.df)
    }
    df <- df %>% mutate(OriginalGene = Gene,
                        Gene = Gene10X)
    ## toconvert <- df %>% subset(Gene != Gene10X) %>% select(ENSGID, Gene, Gene10X)
    ## toconvert.index <- toconvert %>% 
    ## df %>% subset(Gene %in% Gene.ENSEMBL.10X.df$Gene10X & is.na(ENSGID))    
    return(df)
}


if(grepl("2kG.library", SAMPLE)) {
    ## Perturb-seq vs 10X names (from Gavin)
    ## ptb10xNames.df <- read_xlsx("../data/Perturbation 10X names.xlsx") %>% as.data.frame
    ptb10xNames.df <- read.delim(paste0(DATADIR, "/data/220627_add_Perturbation 10X names.txt"), stringsAsFactors=F, check.names=F)
    ## write.table(ptb10xNames.df,"../data/220627_add_Perturbation 10X names.220628.txt", quote=F, row.names=F, sep="\t")
    perturbseq.gene.names.to10X <- function(Gene) stri_replace_all_regex(Gene, pattern = ptb10xNames.df$Symbol, replace = ptb10xNames.df$`Name used by CellRanger`, vectorize=F)
    tenX.gene.names.toperturbseq <- function(Gene) stri_replace_all_regex(Gene, pattern = ptb10xNames.df$`Name used by CellRanger`, replace = ptb10xNames.df$Symbol, vectorize=F) 
}


## gtf.10X.df <- readRDS(paste0(DATADIR, "/data/refdata-cellranger-arc-GRCh38-2020-A_genes.gtf_df.RDS")) ## load 10X gtf file
message(paste0("Loading input GWAS table file from ", opt$input.GWAS.table))
GWAS.df <- GWAS.df.original <- read.delim(opt$input.GWAS.table, stringsAsFactors=F) %>% mutate(original.gene = gene) %>%
    select(-ProgramsInWhichGeneIsInTop100ZScoreSpecificGenes, -ProgramsInWhichGeneIsInTop300ZScoreSpecificGenes, -ProgramsInWhichGeneIsInTop500ZScoreSpecificGenes, -PoPS_Score, -PoPS.Rank, -Top5ProgramsThatContributeToPoPSScore) %>%
    mutate(GeneInGWASTable = TRUE,
           Gene = gene) %>%
    map.ENSGID.SYMBOL

if(grepl("2kG.library", SAMPLE)){
    GWAS.df <- GWAS.df %>%
        ## mutate(original.gene = gene) %>%
        mutate(gene = gene %>% perturbseq.gene.names.to10X)
    
    ## add columns on_2kG_lib and TeloHAEC_ctrl_TPM
    genes.on.2kG.lib <- read_xlsx(paste0(DATADIR, "/data/Table.S1.1 revised 220722.xlsx")) %>%
        mutate(original.gene = `Symbol (for library design)`,
               gene = `Symbol (in CellRanger)`)
    RNAseq_CITV <- read.delim(paste0(DATADIR, "/CAD_SNP_INFO/RNASeq_Telo_Eahy_pm_IL1b_TNF_VEGF.txt"), stringsAsFactors=F)
    TPM.df <- RNAseq_CITV %>% select(Gene_symbol, TeloHAEC.Ctrl_avg) %>% `colnames<-`(c("gene", "TeloHAEC_ctrl_TPM"))


    ## GWAS.df <- GWAS.df %>%
    ##     mutate(on_2kG_lib = gene %in% genes.on.2kG.lib$gene) %>%
    ##     merge(TPM.df, by="gene", all.x=T)
}

perturbed_gene_ary <- barcode.names$Gene %>% unique


## Load data for all genes relation to the topics (not limited to CAD GWAS genes)
message(paste0("Loading input G2P table file from ", opt$cNMF.table))
narrow.df <- read.delim(opt$cNMF.table, header=T, stringsAsFactors=F)## %>% select(-Top5FeaturesThatContributeToPoPSScore)
## narrow.df <- narrow.df[narrow.df$Gene!="NULL", ]

if(grepl("2kG.library", SAMPLE)) {
    narrow.df.TPM <- merge(TPM.df, narrow.df %>% select(-ENSGID) %>% unique, by.x=c("gene"), by.y=c("Gene"), all.y=T)
    narrow.df.TPM$gene[narrow.df.TPM$gene == "MESDC1"] <- "TLNRD1" ## one time fix MESDC1 -> TLNRD1
} else {
    narrow.df.TPM <- narrow.df %>% mutate(gene = Gene)
}

if(!grepl("2kG.library", SAMPLE)) {
    coding_variant.df <- read.delim(opt$coding.variant.df, stringsAsFactors=F)

    GeneWithCodingVariant <- coding_variant.df$CodingVariantGene %>% unique
    GWAS.df <- GWAS.df %>%
        mutate(GeneContainsCodingVariant = (gene %in% GeneWithCodingVariant) | (original.gene %in% GeneWithCodingVariant),
               CodingVariantGene = ifelse(GeneContainsCodingVariant, gene, NA))
} else {
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

    GWAS.df <- GWAS.df %>% 
        ## mutate(Resource=ifelse(Resource=="Aragam_2021","Aragam2021","Harst2017"),
        ##        CredibleSet=paste0(Resource,"-",Lead_SNP_rsID)) %>%
        addCodingGenes(rbind(all.cs.aragam,all.cs.harst)) %>%
        mutate(CredibleSet = ifelse(CredibleSet %in% c("Aragam2021-rs768453105", "Harst2017-rs138120077"), "Aragam2021-rs768453105_Harst2017-rs138120077", CredibleSet))
}

## include TPM table
if(!("TPM" %in% colnames(GWAS.df))) {
    TPM.df <- read.delim(opt$TPM.table, stringsAsFactors=F)
    GWAS.tmp.df <- merge(GWAS.df, TPM.df, by.x="gene", by.y="Gene", all.x=T) 
    GWAS.df <- GWAS.tmp.df
}



########################################################################################
## Define key variables for analysis
GWAS.df <- GWAS.df %>%
    mutate(
        Expressed = (TPM >= 1),
        InPlausibleCellTypeSpecificLocus = ( !LipidLevelsAssociated & (!is.na(get(paste0("RankOfDistanceTo", celltype, "PeakWithVariant"))) | !is.na(get(paste0("MaxABC.Rank.", celltype, "Only"))) | !is.na(CodingVariantGene)) ), ## Also require that locus is not associated with lipid levels
        InPlausibleCellTypeSpecificLocus_includeLipid = (!is.na(get(paste0("RankOfDistanceTo", celltype, "PeakWithVariant"))) | !is.na(get(paste0("MaxABC.Rank.", celltype, "Only"))) | !is.na(CodingVariantGene)),
        InPlausibleLocus = ( !LipidLevelsAssociated & (!is.na(rank_SNP_to_TSS) | !is.na(MaxABC.Rank) | !is.na(CodingVariantGene)) ),
        InPlausibleLocus_includeLipid = !is.na(rank_SNP_to_TSS) | !is.na(MaxABC.Rank) | !is.na(CodingVariantGene),
        TopCandidate = (get(paste0("RankOfDistanceTo", celltype, "PeakWithVariant")) <= 2 | MaxABC.Rank.ECOnly <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE),
        TopCandidateInCellTypeSpecificLocus = TopCandidate & InPlausibleCellTypeSpecificLocus
        ## !!(paste0("ExpressedTopCandidateIn", celltype, "Locus")) := ((get(paste0("RankOfDistanceTo", celltype, "PeakWithVariant")) <= 2 | get(paste0("MaxABC.Rank.", celltype, "Only")) <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE)) & get(paste0("InPlausible", celltype, "Locus")) & Expressed,
        ## TopCandidateInCellTypeSpecificLocus_includeLipid = ((get(paste0("RankOfDistanceTo", celltype, "PeakWithVariant")) <= 2 | get(paste0("MaxABC.Rank.", celltype, "Only")) <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE)) & get(paste0("InPlausible", celltype, "Locus_includeLipid"))
        ## TopCandidate = ((rank_SNP_to_TSS <= 2 | MaxABC.Rank <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE)) & InPlausibleLocus,
        ## TopCandidate_includeLipid = ((rank_SNP_to_TSS <= 2 | MaxABC.Rank <= 2 | GeneContainsCodingVariant) %>% replace_na(FALSE)) & InPlausibleLocus_includeLipid
    ) %>%
    mutate(perturbed_gene = (gene %in% perturbed_gene_ary)) %>%
    as.data.frame()



##########################################################################################
## quick statistics on GWAS.df
## Count of V2G linked genes
cat("Number of genes with V2G links: ",
    GWAS.df %>% filter(TopCandidateInCellTypeSpecificLocus) %>% pull(gene) %>% unique %>% length,
    "\n")

## Count of credible sets that count as "plausibleCellTypeSpecificLocus":
cat("Number of credible sets that count as 'in a plausible CellTypeSpecific locus': ", 
    GWAS.df %>% filter(InPlausibleCellTypeSpecificLocus) %>% pull(CredibleSet) %>% unique() %>% length(), 
    " out of ",
    GWAS.df %>% pull(CredibleSet) %>% unique() %>% length(),
    "\n")


## Count of credible sets that count as "plausibleLocus":
cat("Number of credible sets that count as 'in a plausible locus': ", 
    GWAS.df %>% filter(InPlausibleLocus) %>% pull(CredibleSet) %>% unique() %>% length(), 
    " out of ",
    GWAS.df %>% pull(CredibleSet) %>% unique() %>% length(),
    "\n")



## merge CAD.GWAS.df into narrow.df.TPM to get all expressed gene's information
remove_na <- function(x) x[is.na(x)] <- NULL
overlapping.colnames <- intersect(narrow.df.TPM %>% colnames, GWAS.df %>% colnames)
narrow.df.TPM <- merge(narrow.df.TPM, # %>%
                       GWAS.df %>%
                       select(-one_of((overlapping.colnames[!(overlapping.colnames %in% c("gene"))]))),
                       by=c("gene"), all=T
                       ) %>%
    rowwise() %>%
    mutate(ProgramsLinkedToGene= c(strsplit(ProgramsRegulatedByThisGene,",")[[1]], strsplit(ProgramsInWhichGeneIsInTop300ZScoreSpecificGenes, ",")[[1]]) %>% unique() %>% paste(collapse=',')) %>%
    ## mutate(perturbed_gene = (gene %in% perturbed_gene_ary) | (OriginalGene %in% perturbed_gene_ary)) %>%
    ## ## mutate(ProgramsLinkedToGene= c(strsplit(ProgramsRegulatedByThisGene,",")[[1]], strsplit(ProgramsInWhichGeneIsInTop300ZScoreSpecificGenes, ",")[[1]]) %>% unique() %>% paste(collapse=',')) %>%
    as.data.frame ## because CAD.GWAS.df has dupicated genes, TopicsInWhichGeneIsInTop100ZScoreSpecificGenes sum up to > 100
narrow.df.TPM[is.na(narrow.df.TPM)] <- FALSE

## head(narrow.df.TPM)
## print(colnames(narrow.df.TPM))
## GenesIncNMF <- c(theta.zscore.rank.df$Gene %>% unique)
narrow.df.TPM <- narrow.df.TPM %>%
    ## rowwise %>%
    mutate(AllGenesIncNMFInput = IncNMFAnalysis) %>%
    as.data.frame

## write GWAS table to file
message("write GWAS.df to file")
write.table(GWAS.df, paste0(OUTDIR, opt$trait.name, ".GWAS.df.txt"), sep="\t", quote=F, row.names=F)
write.table(narrow.df.TPM, paste0(OUTDIR, opt$trait.name, ".narrow.df.TPM.txt"), sep="\t", quote=F, row.names=F)

# narrow.df.TPM <- read.delim(paste0(OUTDIR, opt$trait.name, ".narrow.df.TPM.txt"), stringsAsFactors=F)



########################################################################################
## Run statistical Tests
## Create subset dfs based on background conditions
 
## read table that specifies which gene set is for alternative hypothesis and which set of topics to use
message("reading statistical test list")
statistical.test.list.df <- read.delim(paste0("workflow/scripts/program_prioritization/All_GWAS_traits.statistical.test.list_", opt$regulator.analysis.type, ".txt"), stringsAsFactors=F)
## statistical.test.list.df <- read.delim(paste0("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230119_apply_V2G2P_to_more_GWAS_traits/All_GWAS_traits.statistical.test.list.txt"), stringsAsFactors=F)

message("double check statistical test criteria")
statistical.test.list.df[1:2,] ## double check test criteria

num.tests <- opt$num.tests ## actual number of tests requested by the user
statistical.test.list.df <- statistical.test.list.df[1:num.tests,]

## num.tests <- nrow(statistical.test.list.df) ## specify number of statistical tests to conduct
df.list.to.test <- vector("list", num.tests) ## initiate background data storage list
names(df.list.to.test) <- statistical.test.list.df$test.name

## remove batch topics from test
if(grepl("2kG.library", SAMPLE)) {
    ## load topic summary
    TopicSummary.df <- read.delim(paste0(DATADIR, "/data/210730_cNMF_topic_model_analysis.xlsx - TopicCatalogNEW.tsv"), stringsAsFactors=F)

    ## batch effect topic
    batch.topics <- TopicSummary.df %>%
        filter(ProgramCategoryLabel == "Batch") %>%
        pull(ProgramID) %>%
        unique
} else {
    message("reading batch topics")
    filename = paste0(opt$outdirsample, "/batch.topics.txt")
    empty = file.size(filename) == 0L
    if(empty) {
        message(paste0("There is no batch topics in K = ", k))
        batch.topics = character(0)
    } else {
        batch.topics <- read.delim(filename, stringsAsFactors=F, header=F) %>%
            mutate(ProgramID = paste0("K", k, "_", gsub("topic_", "", V1))) %>%
            pull(ProgramID) %>%
            unique
    }
}

## topics_to_test <- setdiff(c(1:k), batch.topics %>% gsub(paste0("K", k, "_"), "", .))

## do not remove batch topic from the test
topics_to_test = c(1:k)




all.test.results.background.list <- vector("list", num.tests)
for(i in 1:num.tests) {
    background.name <- names(df.list.to.test)[i] ## get the name of the background filter
    subset.df <- eval(parse(text = paste0("subset.df <- ", statistical.test.list.df$background.subset.command[i]))) 

    all.test.results.list <- vector("list", k) ## initialize variable
    ## loop over every topic t
    for(t in topics_to_test) {
        topic <- paste0("K", k, "_", t)

        ## new scratch for test that only looks at 'expressed genes in a topic'
        genes.to.test <- statistical.test.list.df$genes.to.test[i]
        topics.to.test <- statistical.test.list.df$topics.to.test[i]

        list.to.test <- subset.df %>%
            rowwise %>%
            mutate(hypothesis.gene = ifelse(eval(parse(text = genes.to.test)), 1, 0),
                   LinkedToTopic.gene = ifelse(topic %in% (get(topics.to.test) %>% strsplit("\\|") %>% unlist %>% as.character), 1, 0)) %>% ## debug
            ## mutate(gene = paste0(gene, "_", Gene, "_", OriginalGene)) %>%
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
    mutate(Gene_NotLinkedToTopic_CandidateGene = ifelse(NotLinkedToTopic_CandidateGene > 300, "", Gene_NotLinkedToTopic_CandidateGene),
           Gene_LinkedToTopic_CandidateGene = ifelse(LinkedToTopic_CandidateGene > 150, "", Gene_LinkedToTopic_CandidateGene)) #Gene_LinkedToTopic_NotCandidateGene))

all.test.results.nonbatch <- all.test.results %>%
    subset(Topic %in% paste0("K", k, "_", setdiff(c(1:k), batch.topics %>% gsub(paste0("K", k, "_"), "", .)))) %>%
    group_by(background) %>%
    mutate(fisher.p.adjust = p.adjust(fisher.p.value, method="fdr"),
           binomial.p.adjust = p.adjust(binomial.p.value, method="fdr"),
           .after="binomial.p.value") %>%
    as.data.frame

## output statistical test results
write.table(all.test.results, file=paste0(OUTDIR, opt$trait.name, ".fisher.exact.test.to.prioritize.topics.txt"), sep="\t", quote=F, row.names=F)
write.table(all.test.results.nonbatch, file=paste0(OUTDIR, opt$trait.name, ".fisher.exact.test.to.prioritize.topics.nonbatch.txt"), sep="\t", quote=F, row.names=F)

head(all.test.results)

## ## reload data
## all.test.results.nonbatch <- read.delim(paste0(OUTDIR, opt$trait.name, ".fisher.exact.test.to.prioritize.topics.nonbatch.txt"), stringsAsFactors=F)



regulator.programGene.test.pairs <- data.frame(regulator = statistical.test.list.df$test.name[seq(2, num.tests, by=2)],
                                               programGene = statistical.test.list.df$test.name[seq(1, num.tests-1, by=2)],
                                               name = statistical.test.list.df[1:num.tests,] %>%
                                                   separate(col="test.name", sep="_geneSet.", into = c("regulator.or.program.gene", "toProcess")) %>%
                                                   separate(col="toProcess", sep="_", into=c("test.handle", "background")) %>% pull(test.handle) %>% unique,
                                               command = statistical.test.list.df$genes.to.test[seq(1, num.tests-1, by=2)],
                                               regulator.selection.command = statistical.test.list.df$background.subset.command[seq(2, num.tests, by=2)])
regulator.programGene.combined.pval.nonbatch.list <- vector("list", nrow(regulator.programGene.test.pairs))

for(i in 1:nrow(regulator.programGene.test.pairs)) {
    regulator.results.nonbatch <- all.test.results.nonbatch %>%
        subset(background == regulator.programGene.test.pairs$regulator[i]) %>%
        `colnames<-`(paste0("Regulator_", colnames(.))) %>%
        dplyr::rename("ProgramID" = "Regulator_Topic")
    programGene.results.nonbatch <- all.test.results.nonbatch %>% subset(background == regulator.programGene.test.pairs$programGene[i]) %>%
        `colnames<-`(paste0("ProgramGene_", colnames(.))) %>%
        dplyr::rename("ProgramID" = "ProgramGene_Topic")

    ## calculate expected number of genes ## adapted from 221115_compute_enrichment.R
    mean.Regulators = regulator.results.nonbatch %>%
        mutate(nRegulators=Regulator_LinkedToTopic_CandidateGene + Regulator_LinkedToTopic_NotCandidateGene) %>%
        pull(nRegulators) %>%
        mean
    ## nPerturbedGeneInCADGWASLoci = GWAS.df %>% filter(perturbed_gene == 1) %>% pull(gene) %>% unique %>% length
    ## nPerturbedGeneCandidate = GWAS.df %>% filter(perturbed_gene == 1 & eval(parse(text = regulator.programGene.test.pairs$command[i] %>% as.character))) %>% pull(gene) %>% unique %>% length
    nPerturbedGene = eval(parse(text = regulator.programGene.test.pairs$regulator.selection.command[i] %>% as.character)) %>% pull(gene) %>% unique %>% length ## old: narrow.df.TPM %>% filter(perturbed_gene == 1) %>% pull(gene) %>% unique %>% length
    nPerturbedGeneCandidate = eval(parse(text = regulator.programGene.test.pairs$regulator.selection.command[i] %>% as.character)) %>% filter(eval(parse(text = regulator.programGene.test.pairs$command[i] %>% as.character))) %>% pull(gene) %>% unique %>% length ## old: narrow.df.TPM %>% filter(perturbed_gene == 1 & eval(parse(text = regulator.programGene.test.pairs$command[i] %>% as.character))) %>% pull(gene) %>% unique %>% length
    expected.Regulators = mean.Regulators * nPerturbedGeneCandidate / nPerturbedGene
    nProgramGeneBackground = narrow.df.TPM %>% filter(IncNMFAnalysis == 1) %>% pull(gene) %>% unique %>% length
    nProgramGeneCandidate = narrow.df.TPM %>% subset(IncNMFAnalysis == 1 & eval(parse(text = regulator.programGene.test.pairs$command[i] %>% as.character))) %>% pull(gene) %>% unique %>% length
    expected.ProgramGenes = 300 * nProgramGeneCandidate / nProgramGeneBackground
    expected.nLinkedGenes = expected.Regulators + expected.ProgramGenes

    regulator.programGene.combined.pval.nonbatch.list[[i]] <- merge(regulator.results.nonbatch,
                                                                    programGene.results.nonbatch,
                                                                    by="ProgramID") %>%
        mutate(LinkedGenes_ChiSquareTestStatistic = -2 * (log(Regulator_fisher.p.value) + log(ProgramGene_fisher.p.value)),
               LinkedGenes_ChiSquare.p.value = pchisq(LinkedGenes_ChiSquareTestStatistic, df = 4, lower.tail=F),
               LinkedGenes_fisher.p.value = Regulator_fisher.p.value * ProgramGene_fisher.p.value,
               LinkedGenes_binomial.p.value = Regulator_binomial.p.value * ProgramGene_binomial.p.value,
               LinkedGenes_fisher.p.adjust = p.adjust(LinkedGenes_fisher.p.value, method="fdr"),
               LinkedGenes_binomial.p.adjust = p.adjust(LinkedGenes_binomial.p.value, method="fdr"),
               LinkedGenes_ChiSquare.p.adjust = p.adjust(LinkedGenes_ChiSquare.p.value, method="fdr"),
               LinkedGenes_fisherNegLog10FDR = -log10(LinkedGenes_fisher.p.adjust),
               LinkedGenes_binomialNegLog10FDR = -log10(LinkedGenes_binomial.p.adjust),
               LinkedGenes_ChiSquareNegLog10FDR = -log10(LinkedGenes_ChiSquare.p.adjust),
               ## LinkedGenes_LinkedToTopic_CandidateGene = Regulator_LinkedToTopic_CandidateGene + ProgramGene_LinkedToTopic_CandidateGene, 
               .after = "ProgramID") %>%
        ## adapted from 211115_compute_enrichment.R
                rowwise %>% 
    mutate(LinkedGenes_LinkedToTopic_CandidateGene = c(Regulator_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character, ProgramGene_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character) %>% unique %>% length,
           LinkedGenes_Gene_LinkedToTopic_CandidateGene = c(Regulator_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character, ProgramGene_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character) %>% unique %>% paste0(collapse=","),
        InGeneSet_Regulator_and_ProgramGene = intersect(Regulator_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character, ProgramGene_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character) %>% unique %>% length,
           Gene_InGeneSet_Regulator_and_ProgramGene = intersect(Regulator_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character, ProgramGene_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character) %>% unique %>% paste0(collapse=","),
           Unique_Regulator_InGeneSet_CandidateGene = setdiff(Regulator_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character, ProgramGene_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character) %>% unique %>% length,
                   Unique_ProgramGene_InGeneSet_CandidateGene = setdiff(ProgramGene_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character, Regulator_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character) %>% unique %>% length,
                   Unique_Regulator_Gene_InGeneSet_CandidateGene = setdiff(Regulator_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character, ProgramGene_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character) %>% unique %>% paste0(collapse=","),
                   Unique_ProgramGene_Gene_InGeneSet_CandidateGene = setdiff(ProgramGene_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character, Regulator_Gene_LinkedToTopic_CandidateGene %>% strsplit(split=",") %>% unlist %>% as.matrix %>% as.character) %>% unique %>% paste0(collapse=","),

           .after = "LinkedGenes_binomialNegLog10FDR") %>%
        mutate(.after = "LinkedGenes_LinkedToTopic_CandidateGene",
               LinkedGenes_Expected = expected.nLinkedGenes,
               LinkedGenes_Enrichment = LinkedGenes_LinkedToTopic_CandidateGene / LinkedGenes_Expected) %>%
    mutate(.before = "Regulator_enrichment",
           Regulator_Expected = expected.Regulators) %>%
    mutate(.before = "ProgramGene_enrichment",
           ProgramGene_Expected = expected.ProgramGenes) %>%
        as.data.frame %>%
## end of adapted code ##
        arrange(LinkedGenes_binomial.p.adjust) %>% mutate(test.name = regulator.programGene.test.pairs$name[i], .after = "ProgramID")
}

regulator.programGene.combined.pval.nonbatch <- do.call(rbind, regulator.programGene.combined.pval.nonbatch.list) %>%
    arrange(LinkedGenes_ChiSquare.p.adjust, desc(LinkedGenes_LinkedToTopic_CandidateGene))

write.table(regulator.programGene.combined.pval.nonbatch, file=paste0(OUTDIR, "/", opt$trait.name, ".program_prioritization.txt"), sep="\t", quote=F, row.names=F)

## create a table with V2G2P prioritized genes for Program Genes
fdr.thr <- 0.05
prioritized.program.genes.df <- regulator.programGene.combined.pval.nonbatch %>%
    subset(ProgramGene_fisher.p.adjust < fdr.thr) %>%
    select(ProgramID, ProgramGene_Gene_LinkedToTopic_CandidateGene) %>%
    separate_rows(ProgramGene_Gene_LinkedToTopic_CandidateGene, sep=",", convert=F) %>%
    `colnames<-`(c("ProgramID", "V2G2P Program Gene")) %>%
    as.data.frame



    
####################################################################################################
## plots ## adapted from TeloHAEC_Perturb-seq_2kG/221115_stimulation_condition_V2G/221115_compute_enrichment.R
source('/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/figures/helper_scripts/plot_helper_functions.R')
## reload data
regulator.programGene.combined.pval.nonbatch <- read.delim(paste0(OUTDIR, "/", opt$trait.name, ".program_prioritization.txt"), stringsAsFactors=F)
## regulator.programGene.combined.pval.nonbatch <- read.delim(paste0(OUTDIR, opt$trait.name, ".combinedRegulatorProgramGeneEnrichmentTest.to.prioritize.topics.nonbatch.txt"), stringsAsFactors=F)
theme.here <- theme(legend.key.size=unit(0.1, unit="in"),
                    legend.text = element_text(size=5),
                    legend.position = "bottom",
                    legend.direction="vertical",
                    legend.margin = margin(1,1,1,1,unit="pt"),
                    panel.spacing = unit(1, units="pt"))
mytheme <- theme_classic() + theme(axis.text = element_text(size = 5),
                                   axis.title = element_text(size = 6),
                                   plot.title = element_text(hjust = 0.5, face = "bold", size=6),
                                   axis.line = element_line(color = "black", size = 0.25),
                                   axis.ticks = element_line(color = "black", size = 0.25))


if(opt$perturbSeq) {
    category_ary <- c("LinkedGenes", "ProgramGene", "Regulator")
} else {
    category_ary <- c("ProgramGene")
}

fdr.thr <- 0.05
for(category in category_ary) {
    expectedNumEnrichedGenes <- regulator.programGene.combined.pval.nonbatch %>% pull(get(paste0(category, "_Expected"))) %>% unique
    toplot <- regulator.programGene.combined.pval.nonbatch %>%
        mutate(significant = ifelse(get(paste0(category, "_", ifelse(category == "LinkedGenes", "fisher", "fisher"), ".p.adjust")) < fdr.thr,
                             ifelse(get(paste0(category, "_", ifelse(category == "LinkedGenes", "fisher", "fisher"), ".p.adjust")) < (fdr.thr / 10),
                             ifelse(get(paste0(category, "_", ifelse(category == "LinkedGenes", "fisher", "fisher"), ".p.adjust")) < (fdr.thr / 100), "***", "**"), ## change to "ChiSquare" in the true slot if you want to combine regulator and program genes by Fisher's method
                             "*"),
                             "")) %>%
        ## add.df.Program.name %>% ## todo (generalize)
        arrange(desc(get(paste0(category, "_LinkedToTopic_CandidateGene"))))
    ## numSignificantPrograms <- nrow(toplot)

    if(grepl("2kG.library", SAMPLE)) {
        toplot <- toplot %>% add.df.Program.name
    } else {
        toplot <- toplot %>% mutate(truncatedLabel = ProgramID)
    }

    if(category == "LinkedGenes") {
        toplot <- toplot %>%
            arrange(desc(LinkedGenes_LinkedToTopic_CandidateGene),
                    desc(InGeneSet_Regulator_and_ProgramGene),
                    desc(Unique_ProgramGene_InGeneSet_CandidateGene),
                    desc(Unique_Regulator_Gene_InGeneSet_CandidateGene))
    }

    label.order <- toplot$truncatedLabel %>% rev
    toplot <- toplot %>%
        mutate(truncatedLabel = factor(truncatedLabel, levels=label.order))
    
    if(category == "LinkedGenes") {
        maxNumGenes <- toplot$LinkedGenes_LinkedToTopic_CandidateGene %>% max
        toplot.here <- toplot %>%
            select(truncatedLabel, InGeneSet_Regulator_and_ProgramGene, Unique_Regulator_InGeneSet_CandidateGene, Unique_ProgramGene_InGeneSet_CandidateGene, significant) %>%
            melt(id.vars=c("truncatedLabel", "significant"), value.name="numGenes", variable.name="LinkType") %>%
            mutate(LinkTypeText = ifelse(grepl("_and_", LinkType), "Regulator and Co-regulated Genes",
                                  ifelse(grepl("Regulator", LinkType), "Regulators", "Co-regulated Genes")),
                   LinkTypeText = factor(LinkTypeText, levels=c("Co-regulated Genes", "Regulators", "Regulator and Co-regulated Genes") %>% rev))

        p <- toplot.here %>% ggplot(aes(x=truncatedLabel, y=numGenes, fill=LinkTypeText)) + geom_col() + coord_flip() + mytheme +
            geom_text(aes(label=significant, y=maxNumGenes*1.1), nudge_x=-0.5, size=3, color='gray') +
            scale_fill_manual(name = "", values = c("gray30", "#38b4f7", "#0141a8") %>% rev) +
            xlab("Programs") + ylab("# Genes Linked to Program") +
            ggtitle(paste0(opt$trait.name, " GWAS trait")) +
            theme.here +
            geom_hline(yintercept=expectedNumEnrichedGenes, linetype="dashed", color="gray")

    } else {
        toplot.here <- toplot %>%
            select(truncatedLabel, paste0(category, "_LinkedToTopic_CandidateGene"), significant)

        p <- toplot.here %>% ggplot(aes(x=truncatedLabel, y=get(paste0(category, "_LinkedToTopic_CandidateGene")))) + geom_col(fill="gray30") + coord_flip() + mytheme +
            geom_text(aes(label=significant, y=max(get(paste0(category, "_LinkedToTopic_CandidateGene")))*1.1), nudge_x=-0.5, size=3, color='gray') +
            xlab("Programs") + ylab("# Genes Linked to Program") +
            ggtitle(paste0(opt$trait.name, " GWAS trait")) +
            theme.here +
            geom_hline(yintercept=expectedNumEnrichedGenes, linetype="dashed", color="gray")
    }

    ## if(opt$perturbSeq) p <- p + scale_fill_manual(name = "", values = c("gray30", "#38b4f7", "#0141a8") %>% rev) 
    filename <- paste0(FIGDIR, "/", SAMPLE, "_K", k, "_dt_", DENSITY.THRESHOLD, "_", opt$trait.name, "_", category, "_GeneCountBarPlot")
    pdf(paste0(filename, ".pdf"), width=2.5, height=3/60*k+1)
    print(p)
    dev.off()

    ## Enrichment Plot
    column.here <- ifelse(category == "LinkedGenes", "LinkedGenes_Enrichment", paste0(category, "_enrichment"))
    toplot <- toplot %>% arrange(desc(get(column.here)))
    label.order <- toplot$truncatedLabel %>% rev
    toplot <- toplot %>%
        mutate(truncatedLabel = factor(truncatedLabel, levels=label.order))
        
    p <- toplot %>% ggplot(aes(x=truncatedLabel, y=get(column.here))) + geom_col(fill="gray30") + coord_flip() + mytheme +
        geom_text(aes(label=significant, y=max(get(column.here))*1.1), nudge_x=-0.5, size=3, color='gray') +
        xlab("Program") + ylab(paste0(category, "Enrichment")) +
        ggtitle(paste0(opt$trait.name, " GWAS trait")) +
        theme.here

    filename <- paste0(FIGDIR, "/", SAMPLE, "_K", k, "_dt_", DENSITY.THRESHOLD, "_", opt$trait.name, "_", category, "_EnrichmentBarPlot")
    pdf(paste0(filename, ".pdf"), width=2.5, height=3/60*k+1)
    print(p)
    dev.off()

}


## output V2G2P genes
fdr.thr <- 0.05
for(key in category_ary) {
    significantGenes.df <- regulator.programGene.combined.pval.nonbatch %>%
        subset(get(paste0(key, "_fisher.p.adjust")) < fdr.thr) %>%
        select(ProgramID, paste0(key, "_Gene_LinkedToTopic_CandidateGene")) %>%
        as.data.frame
    write.table(significantGenes.df, paste0(OUTDIR, "/significant", key, ".df.txt"), sep="\t", quote=F, row.names=F)

    significantGenes <- significantGenes.df %>%
        pull(get(paste0(key, "_Gene_LinkedToTopic_CandidateGene"))) %>%
        paste0(collapse=",") %>%
        strsplit(split=",") %>%
        unlist %>%
        unique
    write.table(significantGenes, paste0(OUTDIR, "/significant", key, ".txt"), sep="\n", quote=F, row.names=F, col.names=F)

    significantGenes.formatted.df <- significantGenes.df %>%
        separate_rows(paste0(key, "_Gene_LinkedToTopic_CandidateGene")) %>%
        mutate(t = gsub(paste0("K", k, "_"), "", ProgramID)) %>%
        arrange(t) %>%
        group_by(get(paste0(key, "_Gene_LinkedToTopic_CandidateGene"))) %>%
        summarize(ProgramID = paste0(ProgramID, collapse=",")) %>%
        `colnames<-`(c(paste0(key,"_Gene_LinkedToTopic_CandidateGene"), "ProgramID")) %>%
        as.data.frame
    write.table(significantGenes.formatted.df, paste0(OUTDIR, "/significant", key, ".formatted.df.txt"), sep="\t", quote=F, row.names=F)
}

## add program gene and regulator membership to linked genes table
if(opt$perturbSeq) {
    linkedSignificantGenes.formatted.df <- read.delim(paste0(OUTDIR, "/significantLinkedGenes.formatted.df.txt"), stringsAsFactors=F) %>% `colnames<-`(c("Gene", "LinkedPrograms"))
    for(key in c("ProgramGene", "Regulator")) {
        assign(paste0(key, "SignificantGenes.formatted.df"), regulator.programGene.combined.pval.nonbatch %>%
                                                             subset(get(paste0("LinkedGenes_fisher.p.adjust")) < fdr.thr) %>%
                                                             select(ProgramID, paste0(key, "_Gene_LinkedToTopic_CandidateGene")) %>%
                                                             separate_rows(paste0(key, "_Gene_LinkedToTopic_CandidateGene")) %>%
                                                             mutate(t = gsub(paste0("K", k, "_"), "", ProgramID)) %>%
                                                             arrange(t) %>%
                                                             group_by(get(paste0(key, "_Gene_LinkedToTopic_CandidateGene"))) %>%
                                                             summarize(ProgramID = paste0(ProgramID, collapse=",")) %>%
                                                             `colnames<-`(c("Gene", key)) %>%
                                                             as.data.frame)
    }
    significantGenes.formatted.df <- merge(linkedSignificantGenes.formatted.df, ProgramGeneSignificantGenes.formatted.df, by="Gene", all.x=T) %>% merge(RegulatorSignificantGenes.formatted.df, by="Gene", all.x=T) %>% `colnames<-`(c("Gene", "LinkedProgram", "PartOfProgram", "RegulatorOfProgram"))
    write.table(significantGenes.formatted.df, paste0(OUTDIR, "/significantLinkedGenes.formatted.df.txt"), sep="\t", quote=F, row.names=F)
}



## scratch 230316
## to find genes that are not regulator or program genes in the final narrow.df.TPM table, then check if they are missed
## not.in.df.genes <- c(narrow.df.TPM %>% filter(!perturbed_gene & ! IncNMFAnalysis) %>% pull(gene))
## narrow.df %>% subset(OriginalGene %in% not.in.df.genes[not.in.df.genes%in% c(narrow.df$Gene, narrow.df$Gene10X, narrow.df$OriginalGene)]) %>% dim
