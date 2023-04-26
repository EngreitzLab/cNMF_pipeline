## Curate a table
## Helen Kang
## 211119

## setwd("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/sample_cNMF_PoPS_input_table/")

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


packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
              "tidyr", "grid", "gtable", "gridExtra","ggrepel","ramify",
              "ggpubr","gridExtra",
              "org.Hs.eg.db","limma","fgsea", "conflicted",
              "cluster","textshape","readxl", 
              "ggdist", "gghalves", "Seurat", "writexl",
              "stringi") 
xfun::pkg_attach(packages)
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


# setwd("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/figures/helper_scripts/")
# source("../../figures/helper_scripts/plot_helper_functions.R")
# source("../../figures/helper_scripts/load_data.R")
# source("../../figures/helper_scripts/load_edgeR.R")
# setwd("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/sample_cNMF_PoPS_input_table/")

source("./workflow/scripts/helper_functions.R")

option.list <- list(
    make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/figures/all_genes/2kG.library.ctrl.only/K25/threshold_0_2/", help="Figure directory"),
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211206_ctrl_only_snakemake/analysis/all_genes/2kG.library.ctrl.only/K25/threshold_0_2/", help="Output directory"),
    ## make_option("--olddatadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/", help="Input 10x data directory"),
    make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/", help="Input 10x data directory"),
    make_option("--topic.model.result.dir", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/210625_snakemake_output/top3000VariableGenes_acrossK/2kG.library/", help="Topic model results directory"),
    make_option("--sampleName", type="character", default="2kG.library.ctrl.only", help="Name of Samples to be processed, separated by commas"),
    make_option("--K.val", type="numeric", default=25, help="K value to analyze"),
    ## make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    ## make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--perturbSeq", type="logical", default=F, help="T for Perturb-seq experiment, F for no perturbation"),

    ## PoPS results
    make_option("--preds_with_cNMF", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6_cNMF60.preds", help="PoPS Score with cNMF input"),
    make_option("--preds_without_cNMF", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/CAD_aug6.preds", help="PoPS Score with cNMF input"),

    
    ## summary plot parameters
    make_option("--test.type", type="character", default="per.guide.wilcoxon", help="Significance test to threshold perturbation results"),
    make_option("--adj.p.value.thr", type="numeric", default=0.1, help="adjusted p-value threshold"),
    make_option("--recompute", type="logical", default=F, help="T for recomputing statistical tests and F for not recompute")

)
opt <- parse_args(OptionParser(option_list=option.list))


## ## all genes directories (for sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/2kG.library/K60/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/"
## opt$K.val <- 60
## opt$sampleName <- "2kG.library"
## opt$perturbSeq <- TRUE

## ## K562 gwps 2k overdispersed genes
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/WeissmanK562gwps/K80/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K80/threshold_0_2/"
## opt$K.val <- 80
## opt$sampleName <- "WeissmanK562gwps"
## opt$perturbSeq <- TRUE


## ## K562 gwps 2k overdispersed genes
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/WeissmanK562gwps/K90/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/"
## opt$K.val <- 90
## opt$sampleName <- "WeissmanK562gwps"
## opt$perturbSeq <- TRUE


OUTDIRSAMPLE <- OUTDIR <- opt$outdir
DATADIR <- opt$datadir
SAMPLE=opt$sampleName
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
k <- opt$K.val
## OUTDIRSAMPLE=paste0(OUTDIR, SAMPLE, "/K",k,"/threshold_", DENSITY.THRESHOLD, "/")
SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)
SUBSCRIPT.SHORT=paste0("k_", k,".dt_",DENSITY.THRESHOLD)
if(!dir.exists(OUTDIR)) dir.create(OUTDIR)
fdr.thr <- 0.05
db <- ifelse(grepl("mouse", SAMPLE), "org.Mm.eg.db", "org.Hs.eg.db")
library(!!db) ## load the appropriate database

if(grepl("2kG.library", SAMPLE)) {
    
    ## map ids ## slow version
    x <- org.Hs.egENSEMBL 
    mapped_genes <- mappedkeys(x)
    entrez.to.ensembl <- as.list(x[mapped_genes]) # EntrezID to Ensembl
    ensembl.to.entrez <- as.list(org.Hs.egENSEMBL2EG) # Ensembl to EntrezID

    y <- org.Hs.egGENENAME
    y_mapped_genes <- mappedkeys(y)
    entrez.to.genename <- as.list(y[y_mapped_genes])
    genename.to.entrez <- as.list(org.Hs.egGENENAME)

    ## map between EntrezID and Gene Symbol
    z <- org.Hs.egSYMBOL
    z_mapped_genes <- mappedkeys(z)
    entrez.to.symbol <- as.list(z[z_mapped_genes])

    entrez.to.symbol <- as.list(org.Hs.egSYMBOL)
    symbol.to.entrez <- as.list(org.Hs.egSYMBOL2EG)
}

## load 10X reference
gtf.10X.df <- readRDS(paste0("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/refdata-cellranger-arc-GRCh38-2020-A_genes.gtf_df.RDS")) ## load 10X gtf file
Gene.ENSEMBL.10X.df <- gtf.10X.df %>% 
    mutate(ENSGID = gene_id,
           Gene10X = gene_name) %>%
    select(Gene10X, ENSGID) %>%
    unique
## gtf <- importGTF("/home/groups/engreitz/Software/cellranger-arc-1.0.1/refdata-cellranger-arc-GRCh38-2020-A/genes/genes.gtf")

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

## ## load Perturb-seq analysis results
## OUTDIRSAMPLE <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/"
## SUBSCRIPT <- "k_60.dt_0_2.minGuidePerPtb_1.minCellPerGuide_2"
file.name <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT,".RData")
print(file.name) 
if(file.exists((file.name))) { 
    print(paste0("loading ",file.name))
    load(file.name)
}

if(opt$perturbSeq) {
    MAST.file.name <- paste0(OUTDIRSAMPLE, "/", SAMPLE, "_MAST_DEtopics.txt")
    print(paste0("loading ", MAST.file.name))
    MAST.df <- read.delim(MAST.file.name, stringsAsFactors=F, check.names=F)
    if(grepl("2kG.library", SAMPLE)){
        MAST.df <- MAST.df %>%
            subset(zlm.model.name == "batch.correction")
        multiTarget.genes <- MAST.df %>%
            subset(grepl("multiTarget", perturbation)) %>% 
            pull(perturbation) %>%
            unique %>%
            gsub("_multiTarget", "", .)
        MAST.df <- MAST.df %>%
            mutate(ProgramID = gsub("topic_", "K60_", primerid)) %>%
            subset(!(perturbation %in% multiTarget.genes) & ## for multiTarget genes, swap in the set where cells have GeneA and GeneA-and-GeneB guides.
                   !grepl("-and-", perturbation)) %>% ## remove the set where cells only have GeneA-and-GeneB guides.
            mutate(perturbation = gsub("_multiTarget", "", perturbation)) %>%
            group_by(zlm.model.name) %>%
            mutate(fdr.across.ptb = p.adjust(`Pr(>Chisq)`, method="fdr")) %>%
            subset(zlm.model.name %in% c("batch.correction")) %>%
            select(-zlm.model.name) %>%
            group_by(ProgramID) %>%
            arrange(desc(coef)) %>%
            mutate(coef_rank = 1:n()) %>%
            as.data.frame
    }
}


## make theta.zscore.rank.df
## if(!("theta.zscore.rank.df" %in% ls())) {
theta.zscore.rank.df <- theta.zscore %>%
    as.data.frame %>%
    mutate(Gene = rownames(.)) %>%
    melt(id.vars="Gene", variable.name="ProgramID", value.name="zscore.specificity") %>%
    mutate(ProgramID = paste0("K", k, "_", ProgramID)) %>%
    group_by(ProgramID) %>%
    arrange(desc(zscore.specificity)) %>%
    mutate(zscore.specificity.rank = 1:n()) %>%
    ungroup %>%
    arrange(ProgramID, zscore.specificity.rank) %>%
    as.data.frame
## }
theta.zscore.rank.df <- theta.zscore.rank.df %>% map.ENSGID.SYMBOL

if(grepl("2kG.library", SAMPLE)) {

    ## load PoPS outputs
    preds <- read.table(file=opt$preds_without_cNMF,header=T, stringsAsFactors=F, sep="\t")

    ## load PoPS processed gene x feature score
    ## load("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/coefs.marginals.feature.outer.prod.RDS") ## takes a while
    OUTDIRPOPS <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/"
    ## PoPS_Score.coefs.all.outer <- read.delim(paste0(OUTDIRPOPS, "/coefs.all.feature.outer.prod.txt"), stringsAsFactors=F) ## also takes a long time to load

    ## top.genes.in.top.features <- read.delim(file=paste0(OUTDIRPOPS, "/top.genes.in.top.features.coefs.txt"), stringsAsFactors=F)
    ## load top features defining each gene's PoPS score
    PoPS_preds.importance.score.key <- read.delim(paste0("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211108_withoutBBJ/outputs/PoPS_preds.importance.score.key.columns.txt"), stringsAsFactors=F)

    ## load top topic defining each gene's PoPS score
    PREFIX <- paste0("CAD_aug6_cNMF", k)
    load(paste0(OUTDIRPOPS, "/", PREFIX, "_coefs.defining.top.topic.RDS"))

    coefs <- read.delim(paste0(OUTDIRPOPS, "/", PREFIX, ".coefs"), header=T, stringsAsFactors=F)
    coefs.df <- coefs[4:nrow(coefs),]
    coefs.cNMF.df <- coefs.df %>% subset(grepl("zscore",parameter))
    coefs.cNMF.names <- coefs.cNMF.df %>% pull(parameter)

    ## load EdgeR log2fcs and p-values
    log2fc.edgeR <- read.table(paste0("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/EdgeR/ALL_log2fcs_dup4_s4n3.99x.txt"), header=T, stringsAsFactors=F)
    p.value.edgeR <- read.table(paste0("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/EdgeR/ALL_Pvalues_dup4_s4n3.99x.txt"), header=T, stringsAsFactors=F)

    ## load known CAD gene set
    params.known.CAD.Gene.set <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/known_CAD_gene_set.txt", stringsAsFactors=F, header=F) %>% as.matrix %>% as.character
}


## Helper function to add EntrezID and ENSGID
add.EntrezID.ENSGID <- function(df) {
    return ( df %>%
             mutate(EntrezID = symbol.to.entrez[.$Gene %>% as.character] %>% sapply("[[",1) %>% as.character) %>%
            mutate(ENSGID = entrez.to.ensembl[.$EntrezID %>% as.character] %>% sapply("[[", 1) %>% as.character) 
            )
}
    
## ## MAST results
## MAST.df.4n.original <- read.delim(paste0("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/2kG.library4n3.99x_MAST.txt"), stringsAsFactors=F, check.names=F)
## MAST.df <- MAST.df.4n.original %>%
##     subset(!grepl("multiTarget", perturbation)) %>%
##     group_by(zlm.model.name) %>%
##     mutate(fdr.across.ptb = p.adjust(`Pr(>Chisq)`, method="fdr")) %>%
##     subset(zlm.model.name %in% c("batch.correction")) %>%
##     select(-zlm.model.name) %>%
##     mutate(ProgramID = gsub("topic_", "K60_", primerid)) %>%
##     group_by(ProgramID) %>%
##     arrange(desc(coef)) %>%
##     mutate(coef_rank = 1:n()) %>%
##     as.data.frame
## sigPerturbationsProgram <- MAST.df %>%
##     subset(fdr.across.ptb < fdr.thr) %>%
##     as.data.frame

if(grepl("2kG.library", SAMPLE)) {
    ## curated gene name conversion list
    geneNameConversion <- read.delim(paste0(DATADIR, "/heterogenous_geneName_conversion_toMatchENSGID.txt"), stringsAsFactors=F) %>%
        separate(col=Conversion, into=c("fromGene", "toGene"), sep=" -> ", remove=F)
    geneNameConversionBackToPerturbseq <- geneNameConversion %>% filter(!KeepForGWASTable)

    ## helper function to conver gene names
    convert.gene.names.toMatchENSGID <- function(Gene) stri_replace_all_regex(Gene, pattern = geneNameConversion$fromGene, replace = geneNameConversion$toGene, vectorize=F)

    convert.gene.names.toMatchPerturbseq <- function(Gene) stri_replace_all_regex(Gene, pattern = geneNameConversionBackToPerturbseq$toGene, replace = geneNameConversionBackToPerturbseq$fromGene, vectorize=F)
}

## new gene name conversion method ## 220628
ptb10xNames.df <- read.delim(paste0(DATADIR, "/220627_add_Perturbation 10X names.txt"), stringsAsFactors=F, check.names=F)
perturbseq.gene.names.to10X <- function(Gene) {
    if(Gene %in% ptb10xNames.df$Symbol) {
        out <- stri_replace_all_regex(Gene, pattern = ptb10xNames.df$Symbol, replace = ptb10xNames.df$`Name used by CellRanger`, vectorize=F)
    } else {
        out <- Gene
    }
    return(out)
}
tenX.gene.names.toperturbseq <- function(Gene) stri_replace_all_regex(Gene, pattern = ptb10xNames.df$`Name used by CellRanger`, replace = ptb10xNames.df$Symbol, vectorize=F) 

print("loaded all prerequisite data")


##########################################################################################
ptb.topic.thr <- 0.05 ## FDR threshold
## Column: ProgramsRegulatedByThisGene
if(opt$perturbSeq) {
    ProgramsRegulatedByThisGene.df <- MAST.df %>%
        mutate(OriginalGene = perturbation,
               Gene = perturbseq.gene.names.to10X(perturbation),
               t = gsub("topic_", "", primerid) %>% as.numeric,
               ProgramID = paste0("K", k, "_",  t)) %>%
        mutate(significant = fdr.across.ptb < ptb.topic.thr &
                   (coef > log(1.1) | coef < log(0.9)),
               ProgramID = ifelse(significant, ProgramID, "NA")) %>%
        ## subset(fdr.across.ptb < ptb.topic.thr) %>%
        ## mutate(Gene = convert.gene.names.toMatchENSGID(perturbation),         
        map.ENSGID.SYMBOL
    ## mutate(EntrezID = symbol.to.entrez[.$Gene %>% as.character] %>% sapply("[[",1) %>% as.character) %>%
    ## mutate(ENSGID = entrez.to.ensembl[.$EntrezID %>% as.character] %>% sapply("[[", 1) %>% as.character)
    ProgramsRegulatedByThisGene.tokeep <- ProgramsRegulatedByThisGene.df %>%
        arrange(t) %>%
        ## select(Gene, OriginalGene, ENSGID, EntrezID, ProgramID) %>%
        select(Gene, OriginalGene, Gene10X, ENSGID, ProgramID) %>%
        unique %>% 
        ## group_by(ENSGID, Gene, OriginalGene, EntrezID) %>%
        group_by(ENSGID, Gene, Gene10X, OriginalGene) %>%        
        summarize(ProgramsRegulatedByThisGene=paste0(ProgramID, collapse="|")) %>%
        mutate(ProgramsRegulatedByThisGene = gsub("NA[|]|[|]NA", "", ProgramsRegulatedByThisGene),
               ProgramsRegulatedByThisGene = gsub("[|]NA[|]", "|", ProgramsRegulatedByThisGene),
               perturbed_gene = TRUE) %>% ## clean up "NA"s
        ## mutate(Gene = convert.gene.names.toMatchENSGID(rownames(.))) %>%
        ## mutate(Gene = perturbseq.gene.names.to10X(Gene)) %>%
        as.data.frame
    ## if(grepl("2kG.library", SAMPLE)) {
    ##     ProgramsRegulatedByThisGene.df <- ProgramsRegulatedByThisGene.df %>% mutate(PreviousConversionGene = convert.gene.names.toMatchENSGID(perturbation))
    ##     ProgramsRegulatedByThisGene.tokeep <- ProgramsRegulatedByThisGene.df %>%
    ##         mutate(PreviousConversionGene = convert.gene.names.toMatchENSGID(Gene)) %>%
    ##         as.data.frame
    ## } 
}

## ## ## ## Column: ProgramsInWhichGeneIsInTop100ZScoreSpecificGenes
## ## ## params.top.n.genes.in.topic <- 100
## if(grepl("2kG.library", SAMPLE)) {
##     ProgramsInWhichGeneIsExpressed.df <- theta.zscore.rank.df %>%
##         mutate(OriginalGene = Gene,
##                PreviousConversionGene = convert.gene.names.toMatchENSGID(Gene),
##                Gene = perturbseq.gene.names.to10X(Gene)) %>%
##         ## add.EntrezID.ENSGID %>%
##         as.data.frame
## } else {
##     ProgramsInWhichGeneIsExpressed.df <- theta.zscore.rank.df %>%
##         ## mutate(OriginalGene = Gene,
##         ##        PreviousConversionGene = convert.gene.names.toMatchENSGID(Gene),
##         ##        Gene = perturbseq.gene.names.to10X(Gene)) %>%
##         ## add.EntrezID.ENSGID %>%
##         as.data.frame
## }

IncNMFAnalysis.tokeep <- theta.zscore.rank.df %>%
    select(Gene, Gene10X, ENSGID) %>%
    unique %>%
    mutate(IncNMFAnalysis = TRUE) %>%
    as.data.frame

params.top.n.genes.in.topic.list <- c(100, 300, 500)
ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes.tokeep <- Reduce(
    function(x, y, ...) full_join(x, y, by = c("Gene", "Gene10X", "ENSGID"), ...), ## want to join all the data frames that has columns c("Gene", "ProgramsInWhichGeneIsInTopNZScoreSpecificGenes")
    out <- lapply(params.top.n.genes.in.topic.list, function (params.top.n.genes.in.topic) { # for each selected number of top genes we want to include
        column.name <- paste0("ProgramsInWhichGeneIsInTop", params.top.n.genes.in.topic, "ZScoreSpecificGenes") # store the column name in a variable
        ## out <- ProgramsInWhichGeneIsExpressed.df %>%
        ##     group_by(ProgramID) %>% 
        ## arrange(desc(Topic.zscore)) %>% # sort the z-score within each topic
        ## slice(1:params.top.n.genes.in.topic) %>% # select the top N genes
        out <- theta.zscore.rank.df %>%
            mutate(IsProgramGene = zscore.specificity.rank <= params.top.n.genes.in.topic,
                   ProgramID = ifelse(IsProgramGene, ProgramID, "NA")) %>%
            ## subset(zscore.specificity.rank <= params.top.n.genes.in.topic) %>%
            ## ungroup %>%
            group_by(Gene, Gene10X, ENSGID) %>% # for each gene
            summarize(!!column.name := paste0(ProgramID, collapse="|")) %>% # paste the topics each gene links to by ","
            mutate(!!column.name := gsub("NA[|]|[|]NA", "", get(column.name)),
                   !!column.name := gsub("[|]NA[|]", "|", get(column.name))) %>% ## clean up "NA"s
            as.data.frame
    }))
## if(grepl("2kG.library", SAMPLE)) {
##     ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes.tokeep <- ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes.tokeep %>%     
##         ## mutate(Gene = convert.gene.names.toMatchENSGID(Gene)) %>%
##         mutate(OriginalGene = Gene) %>%
##         mutate(PreviousConversionGene = convert.gene.names.toMatchENSGID(Gene)) %>%
##         mutate(Gene = perturbseq.gene.names.to10X(Gene)) %>%
##         add.EntrezID.ENSGID %>%
##         as.data.frame
## }

## ## Column: PoPSEnrichedProgramsRegulatedByThisGene ## FDR < 0.05 and in PoPS prioritized features
## if(opt$perturbSeq & grepl("2kG.library", SAMPLE)){
## cNMF.PoPS.enriched.topics <- coefs.cNMF.names %>% gsub("zscore_K60_topic", "K60_", .) ## get the list of PoPS prioritized cNMF topics and reformat
## PoPSEnrichedProgramsRegulatedByThisGene.tokeep <- ProgramsRegulatedByThisGene.df %>%
##     select(Gene, OriginalGene, PreviousConversionGene, ENSGID, ProgramID) %>%
##     subset(ProgramID %in% cNMF.PoPS.enriched.topics) %>%
##     group_by(ENSGID, Gene, OriginalGene, PreviousConversionGene) %>%
##     summarize(PoPSEnrichedProgramsRegulatedByThisGene=paste0(ProgramID, collapse="|")) %>%
##     ## mutate(Gene = convert.gene.names.toMatchENSGID(Gene)) %>%
##     ## mutate(PreviousConversionGene = convert.gene.names.toMatchENSGID(Gene)) %>%
##     ## mutate(Gene = perturbseq.gene.names.to10X(Gene)) %>%
##     as.data.frame
## }

## ## Column: PoPSEnrichedProgramsInWhichGeneIsInTop100ZScoreSpecificGenes
## cNMF.PoPS.enriched.topics <- coefs.cNMF.names %>% gsub("zscore_K60_topic", "K60_", .) ## get the list of PoPS prioritized cNMF topics and reformat
## params.top.n.genes.in.topic <- 300
## PoPSEnrichedProgramsInWhichGeneIsInTop300ZScoreSpecificGenes.tokeep <- ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes.tokeep %>% ## ProgramsInWhichGeneIsExpressed.df %>%
##     subset(ProgramID %in% cNMF.PoPS.enriched.topics) %>%
##     ## group_by(ProgramID) %>%
##     ## arrange(desc(Topic.zscore)) %>%
##     ## slice(1:params.top.n.genes.in.topic) %>%
##     ## ungroup %>%
##     subset(zscore.specificity.rank <= params.top.n.genes.in.topic) %>%
##     group_by(Gene, PreviousConversionGene, OriginalGene, ENSGID) %>%
##     summarize(PoPSEnrichedProgramsInWhichGeneIsInTop300ZScoreSpecificGenes = paste0(ProgramID, collapse="|")) %>%
##     ## mutate(Gene = convert.gene.names.toMatchENSGID(Gene)) %>%
##     ## mutate(PreviousConversionGene = convert.gene.names.toMatchENSGID(Gene)) %>%
##     ## mutate(Gene = perturbseq.gene.names.to10X(Gene)) %>%
##     ## add.EntrezID.ENSGID %>%
##     as.data.frame


## ## PoPS.Score
## PoPS.Score.tokeep <- preds %>% select(ENSGID, PoPS_Score) %>%
##     mutate(EntrezID = ensembl.to.entrez[.$ENSGID %>% as.character] %>% sapply("[[",1) %>% as.character) %>%
##     mutate(Gene = entrez.to.symbol[.$EntrezID %>% as.character] %>% sapply("[[",1) %>% as.character,
##            ## Gene = convert.gene.names.toMatchENSGID(Gene)) %>%
##            PreviousConversionGene = convert.gene.names.toMatchENSGID(Gene)) %>%
##     mutate(OriginalGene = Gene) %>%
##     mutate(Gene = perturbseq.gene.names.to10X(Gene)) %>%
##     subset(!(EntrezID == "NULL" & Gene == "NULL" & PreviousConversionGene == "NULL" & OriginalGene == "NULL")) %>%
##     as.data.frame


## ## ## PoPS.Rank [rank of the score among genes near this CredibleSet/GWAS signal]
## ## create later when we merge in CredibleSet/GWAS signal

## ## ## Top5ProgramsThatContributeToPoPSScore
## ## Top5ProgramsThatContributeToPoPSScore.tokeep <- coefs.defining.top.topic.df %>%
## ##     mutate(ProgramID = gsub("zscore_K60_topic", "K60_", topic),
## ##            Gene = convert.gene.names.toMatchENSGID(Gene)) %>%
## ##     group_by(Gene) %>%
## ##     arrange(desc(gene.feature_x_beta)) %>%
## ##     slice(1:5) %>%
## ##     summarize(Top5ProgramsThatContributeToPoPSScore = paste0(ProgramID, collapse="|")) %>%
## ##     add.EntrezID.ENSGID %>%
## ##     as.data.frame

## ## Top5FeaturesThatContributeToPoPSScore
## Top5FeaturesThatContributeToPoPSScore.tokeep <- PoPS_preds.importance.score.key %>%
##     group_by(Gene) %>%
##     arrange(desc(gene.feature_x_beta)) %>%
##     slice(1:5) %>%
##     mutate(to.display=paste0(pathway, ":", Long_Name)) %>%
##     summarize(Top5FeaturesThatContributeToPoPSScore = paste0(to.display, collapse="|")) %>%
##     ## mutate(Gene = convert.gene.names.toMatchENSGID(Gene)) %>%
##     mutate(OriginalGene = Gene) %>% 
##     mutate(PreviousConversionGene = convert.gene.names.toMatchENSGID(Gene)) %>%
##     mutate(Gene = perturbseq.gene.names.to10X(Gene)) %>%
##     add.EntrezID.ENSGID %>%
##     as.data.frame
## }

## ## DoesThisGeneWhenPerturbedRegulateAKnownCADGene 
## if(opt$perturbSeq & grepl("2kG.library", SAMPLE)){
##     CAD.gene.regulator.log2fc <- log2fc.edgeR %>%
##         subset(grepl(paste0(params.known.CAD.Gene.set, ":ENSG") %>%
##                      paste0(collapse="|"), .$gene))
##     CAD.gene.regulator.p.value <- p.value.edgeR %>%
##         subset(grepl(paste0(params.known.CAD.Gene.set, ":ENSG") %>%
##                      paste0(collapse="|"), .$gene))
##     params.EdgeR.p.value.thr <- "0.05"
##     CAD.gene.regulator.df <- CAD.gene.regulator.p.value %>%
##         melt(id.vars="genes", value.name="p.value", variable.name="perturbation") %>%
##         subset(p.value < params.EdgeR.p.value.thr)
##     DoesThisGeneWhenPerturbedRegulateAKnownCADGene.tokeep <- CAD.gene.regulator.df %>%
##         separate(genes, into=c("Gene", "ENSGID"), sep=":") %>%
##         group_by(perturbation) %>%
##         summarize(DoesThisGeneWhenPerturbedRegulateAKnownCADGene = paste0(Gene, collapse="|"))
##     colnames(DoesThisGeneWhenPerturbedRegulateAKnownCADGene.tokeep)[colnames(DoesThisGeneWhenPerturbedRegulateAKnownCADGene.tokeep) == "perturbation"] <- "Gene"
##     DoesThisGeneWhenPerturbedRegulateAKnownCADGene.tokeep <- DoesThisGeneWhenPerturbedRegulateAKnownCADGene.tokeep %>%
##         ## mutate(Gene = convert.gene.names.toMatchENSGID(Gene)) %>%
##         mutate(OriginalGene = Gene) %>%
##         mutate(PreviousConversionGene = convert.gene.names.toMatchENSGID(Gene)) %>%
##         mutate(Gene = perturbseq.gene.names.to10X(Gene)) %>%
##         add.EntrezID.ENSGID %>%
##         as.data.frame
## }

## put together all
if(opt$perturbSeq){
    ## if(grepl("2kG.library", SAMPLE)) {
    ##     list.all <- list(ProgramsRegulatedByThisGene = ProgramsRegulatedByThisGene.tokeep,
    ##                      ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes = ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes.tokeep,
    ##                      PoPSEnrichedProgramsRegulatedByThisGene = PoPSEnrichedProgramsRegulatedByThisGene.tokeep,
    ##                      PoPSEnrichedProgramsInWhichGeneIsInTop300ZScoreSpecificGenes = PoPSEnrichedProgramsInWhichGeneIsInTop300ZScoreSpecificGenes.tokeep,
    ##                      PoPS.Score = PoPS.Score.tokeep,
    ##                      DoesThisGeneWhenPerturbedRegulateAKnownCADGene = DoesThisGeneWhenPerturbedRegulateAKnownCADGene.tokeep,
    ##                      IncNMFAnalysis = IncNMFAnalysis.tokeep) ## caveat: Not all perturbations are in common Gene Symbol format (e.g. FAM212A, Icam2, Slc9a3r2)
    ## } else {
        list.all <- list(ProgramsRegulatedByThisGene = ProgramsRegulatedByThisGene.tokeep,
                         ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes = ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes.tokeep,
                         IncNMFAnalysis = IncNMFAnalysis.tokeep) ## caveat: Not all perturbations are in common Gene Symbol format (e.g. FAM212A, Icam2, Slc9a3r2)
    ## }
} else {

    list.all <- list(ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes = ProgramsInWhichGeneIsInTopNListZScoreSpecificGenes.tokeep,
                     ## PoPSEnrichedProgramsInWhichGeneIsInTop300ZScoreSpecificGenes = PoPSEnrichedProgramsInWhichGeneIsInTop300ZScoreSpecificGenes.tokeep,
                     ## PoPS.Score = PoPS.Score.tokeep,
                     IncNMFAnalysis = IncNMFAnalysis.tokeep) 
}
## if(grepl("2kG.library", SAMPLE)) {
##     df.all <- Reduce(function(x, y, ...) full_join(x, y %>% select(-OriginalGene, -PreviousConversionGene), by = c("Gene", "ENSGID", "EntrezID"), ...), list.all)
## } else {
    df.all <- Reduce(function(x, y, ...) full_join(x, y, by = c("Gene", "Gene10X", "ENSGID"), ...), list.all)
## }
df.all[df.all == "NULL"] <- NA
na.row.index <- which(df.all$Gene %>% is.na & df.all$ENSGID %>% is.na & df.all$PoPS_Score %>% is.na)

## Perturb-seq specific trimming of duplicated Gene names
if(opt$perturbSeq) {
    ## if(grepl("2kG.library", SAMPLE)) {
    ##     df.all <- df.all %>% subset(!(is.na(Gene) & is.na(ProgramsRegulatedByThisGene) & !(Gene %in% c(theta.zscore %>% rownames, MAST.df$perturbation %>% unique)) & is.na(ENSGID))) %>% ## combine duplicated Gene names' rows
    ##         group_by(Gene, ENSGID) %>% ## consider implement match with 10X
    ##         summarize_all(., function(x) {
    ##             out <- paste0(x %>% unique, collapse="|")
    ##             if(grepl("TRUE", out %>% as.character)) return("TRUE") else return(out)
    ##         }) %>%
    ##         as.data.frame
    ## } else {
        df.all <- df.all %>% subset(!(is.na(Gene) & is.na(ProgramsRegulatedByThisGene) & !(Gene %in% c(theta.zscore %>% rownames, MAST.df$perturbation %>% unique)) & is.na(ENSGID))) %>% ## combine duplicated Gene names' rows
            group_by(Gene, Gene10X, ENSGID) %>%
            summarize_all(., function(x) {
                out <- paste0(x %>% unique, collapse="|")
                if(grepl("TRUE", out %>% as.character)) return("TRUE") else return(out)
            }) %>%
            mutate(IncNMFAnalysis = ifelse(IncNMFAnalysis == "TRUE", TRUE, FALSE),
                   perturbed_gene = ifelse(perturbed_gene == "TRUE", TRUE, FALSE)) %>%
            as.data.frame
    ## }        
}


## write table
print("writing prepared table to file")
write.table(df.all, file=paste0(OUTDIRSAMPLE, "/prepare_compute_enrichment.txt"), row.names=F, quote=F, sep="\t")

