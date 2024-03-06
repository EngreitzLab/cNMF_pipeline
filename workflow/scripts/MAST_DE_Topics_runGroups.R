## Helen Kang
## run MAST with batch correction in groups of 200 perturbations
## 230128 (adapted from 220217)




library(conflicted)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("first", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("melt", "reshape2")
conflict_prefer("filter", "dplyr")

packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
              "tidyr", "grid", "gtable", "gridExtra","ggrepel",#"ramify",
              "ggpubr","gridExtra", "parallel", "future",
              "org.Hs.eg.db","limma","conflicted", #"fgsea", 
              "cluster","textshape","readxl", 
              "ggdist", "gghalves", "Seurat", "writexl", "SingleCellExperiment", "MAST") #              "GGally","RNOmni","usedist","GSEA","clusterProfiler","IsoplotR","wesanderson",
xfun::pkg_attach(packages)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")


##########################################################################################
## Constants and Directories
option.list <- list(
    make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
    make_option("--sampleName", type="character", default="2kG.library", help="sample name"),
    make_option("--barcode.names", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/2kG.library.barcodes.tsv", help="barcodes.tsv for all cells"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--outdirsample", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/", help="path to cNMF analysis results"), ## or for 2n1.99x: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/Perturb_2kG_dup4/K60/threshold_0_2/"
    make_option("--scatteroutput", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/MAST/K80/threshold_0_2/", help="path to gene breakdown table output"),
    make_option("--gene.group.list", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/MAST/GeneNames_Group43.txt"),
    make_option("--scatter.gene.group", type="numeric", default=50, help="Gene group index"),
    
    ## script dir
    make_option("--scriptdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts/", help="location for this script and functions script")

)
opt <- parse_args(OptionParser(option_list=option.list))


## ## sdev debug K562 gwps
## opt$sampleName <- "WeissmanK562gwps"
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/data/K562_gwps_raw_singlecell_01_metadata.txt"
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K35/threshold_0_2/"
## opt$scatteroutput <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/WeissmanK562gwps/MAST/K80/threshold_0_2/"
## opt$gene.group.list <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/MAST/GeneNames_Group43.txt"
## opt$scatter.gene.group <- 43
## opt$scriptdir <- "/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts"
## opt$K.val <- 80

## ## no IL1B
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/data/no_IL1B.barcodes.txt"
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/analysis/top2000VariableGenes/no_IL1B/K15/threshold_0_2/"
## opt$scatteroutput <- "/scratch/groups/engreitz/Users/kangh/cNMF_pipeline/tutorials/2306_V2G2P_prep/top2000VariableGenes/no_IL1B/MAST/K15/threshold_0_2/"
## opt$gene.group.list <- "/scratch/groups/engreitz/Users/kangh/cNMF_pipeline/tutorials/2306_V2G2P_prep/top2000VariableGenes/no_IL1B/MAST/GeneNames_Group12.txt"
## opt$scatter.gene.group <- 12
## opt$sampleName <- "no_IL1B"
## opt$K.val <- 15
## opt$density.thr <- 0.2
## opt$scriptdir <- "/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts"

## ## ish SCHEMA RNA-seq data
## opt$barcode.namse <- "/stanley/sheng-lab/helenkang/SCHEMA_bulkRNAseq_analysis/240223_initial_analysis/outputs/SCHEMA_integration_metadata.txt"
## opt$outdirsample <- "/stanley/sheng-lab/helenkang/SCHEMA_bulkRNAseq_analysis/240301_snakemake/analysis/top2000VariableGenes/SCHEMA_integration/K3/threshold_0_2/"
## opt$scatteroutput <- "/broad/hptmp/Helen/240301_snakemake_240303Updated/top2000VariableGenes/SCHEMA_integration/MAST/K3/threshold_0_2/"
## opt$gene.group.list <- "/broad/hptmp/Helen/240301_snakemake_240303Updated/top2000VariableGenes/SCHEMA_integration/MAST/GeneNames_Group8.txt"
## opt$scatter.gene.group <- 8
## opt$sampleName <- "SCHEMA_integration"
## opt$K.val <- 3
## opt$density.thr <- 0.2
## opt$scriptdir <- "/seq/lincRNA/Helen/cNMF_pipeline/Perturb-seq/workflow/scripts"


k <- opt$K.val 
SAMPLE <- opt$sampleName
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
## OUTDIRSAMPLE <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K31/threshold_0_2/"
OUTDIRSAMPLE <- opt$outdirsample
SCATTEROUTDIR <- opt$scatteroutput
SCATTERINDEX <- opt$scatter.gene.group

## INPUTDIR <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/inputs/"
OUTDIR <- OUTDIRSAMPLE
check.dir <- c(OUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

## mytheme <- theme_classic() + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, face = "bold"))
## palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)
## p.adjust.thr <- 0.1



##########################################################################################
## load data
## ## load known CAD genes
## CAD.genes <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/known_CAD_gene_set.txt", header=F, stringsAsFactors=F) %>% as.matrix %>% as.character

## ## load topic model results
## cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
## print(cNMF.result.file)
## if(file.exists(cNMF.result.file)) {
##     print("loading cNMF result file")
##     load(cNMF.result.file)
## }

## file.name <- paste0(OUTDIRSAMPLE,"/cNMFAnalysis.",SUBSCRIPT,".RData")
## print(file.name) 
## if(file.exists((file.name))) { 
##     print(paste0("loading ",file.name))
##     load(file.name) 
## }

## ## add UMI / cell information
## umiPerCell.df <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/outputs/UMIPerCell.txt", stringsAsFactors=F, row.names=1) %>% mutate(long.CBC = rownames(.)) %>%
##     separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
##     separate(col="CBC", into=c("CBC", "sample"), sep="-scRNAseq_2kG_", remove=F)

## sample.to.10X.lane <- data.frame(sample_num = 1:20,
##                                  sample = umiPerCell.df$sample %>% unique %>% sort)
## umiPerCell.df <- merge(umiPerCell.df, sample.to.10X.lane, by="sample") %>%
##     mutate(CBC_10x = paste0(CBC, "-", sample_num))  

## ## add guide / cell information
## guidePerCell.df <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/inputs/_UPDATED_ALL_SAMPLES_dup4_NON_NA_CALLS_FOR_EA_CBC_FULL_INFO.txt", stringsAsFactors=F) 

## ## add genes detected / cell information
## geneDetectedPerCell.df <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/outputs/geneDetectedPerCell.txt", stringsAsFactors=F) %>%
##     mutate(long.CBC = rownames(.)) %>%
##     separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
##     separate(col="CBC", into=c("CBC", "sample"), sep="-scRNAseq_2kG_", remove=F) %>%
##     merge(sample.to.10X.lane, by="sample") %>%
##     mutate(CBC_10x = paste0(CBC, "-", sample_num))  


## load gene list
gene.ary <- read.delim(opt$gene.group.list, header=F, stringsAsFactors=F) %>% unlist %>% as.character

## load log2(TPM + 1)
load(paste0(opt$scatteroutput, "/", SAMPLE, "_MAST_log2TPM_barcodes.RDS"))

print("organize meta data")
## get metadata from log2.X.full rownames
# if(SAMPLE %in% c("2kG.library", "Perturb_2kG_dup4")) {
if( grepl("2kG.library|Perturb_2kG_dup4", SAMPLE) ) {
    barcode.names <- read.table(opt$barcode.names, header=F, stringsAsFactors=F) %>% `colnames<-`("long.CBC")
    ## rownames(omega) <- barcode.names %>% pull(long.CBC) %>% gsub("CSNK2B-and-CSNK2B", "CSNK2B",.)
    meta_data <- barcode.names %>%
        rownames %>%
        as.data.frame %>%
        `colnames<-`("long.CBC") %>%
        separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
        ## separate(col="CBC", into=c("CBC", "sample"), sep="-", remove=F) %>%
        separate(col="CBC", into=c("CBC", "sample"), sep="-scRNAseq_2kG_", remove=F) %>%
        mutate(Gene = gsub("-TSS2", "", Gene.full.name),
               CBC = gsub("RHOA-and-", "", CBC)) %>%
        filter(Gene %in% c(gene.ary, "negative-control", "safe-targeting")) %>%
        as.data.frame
    sample.to.10X.lane <- data.frame(sample_num = 1:20,
                                     sample = meta_data$sample %>% unique %>% sort)

    ## meta_data <- merge(meta_data, sample.to.10X.lane, by="sample") %>%
    ##     mutate(CBC_10x = paste0(CBC, "-", sample_num)) %>%
    ##     merge(guidePerCell.df %>% select(CBC_10x, guides_per_cbc, max_umi_ct), by="CBC_10x") %>%
    ##     merge(umiPerCell.df, by="long.CBC") %>%
    ##     merge(geneDetectedPerCell.df, by.x="long.CBC", by.y=0)
} else {
    ## barcode.names <- read.table(opt$barcode.names, header=T, stringsAsFactors=F) ## %>% `colnames<-`("long.CBC")
    ## print("finished loading barcode names")
    ## print(paste0("omega dimensions: ", dim(omega)))
    ## print(paste0("barcode names dimensions: ", dim(barcode.names)))
    meta_data <- barcode.names.subset %>% filter(Gene %in% c(gene.ary, "negative-control"))
    meta_data <- meta_data %>% mutate(sample = factor(sample))
}    



## load model fitting fomulas
## test.cmd.df <- read.delim(paste0(INPUTDIR, "MAST_model_formulas.txt"), stringsAsFactors=F)
test.cmd.df <- data.frame(zlm.model.command = "~condition + lane",
                          zlm.model.name = "batch.correction")

## ## 220222 scratch
## if ( !( "ann.omega" %in% ls()) ) {
##     if( grepl("2kG.library|Perturb_2kG_dup4", SAMPLE) ) {
##         ann.omega <- omega %>%
##             as.data.frame %>%
##             mutate(long.CBC = rownames(.)) %>%
##             separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
##             mutate(Gene = gsub("-TSS2$", "", Gene.full.name))
##     } else {
##         ann.omega <- omega %>% as.data.frame %>%
##             mutate(CBC = rownames(.)) %>%
##             merge(barcode.names, by="CBC", all.x=T)
##         ## meta_data <- barcode.names
##     }
## } else {                        
##     ann.omega <- ann.omega %>%
##         mutate(Gene = gsub("-TSS2$", "", Gene.full.name)) ## remove TSS2 annotation
## }
## test

## omega.tpm <- omega %>% as.matrix %>% apply(2, function(x) x / sum(x) * 1000000) ## convert to TPM
## ## end of test
## ## omega.tpm <- ann.omega[,1:k] %>% as.matrix %>% apply(2, function(x) x / sum(x) * 1000000) ## convert to TPM
## log2.omega <- (omega.tpm + 1) %>% log2 ## log2(TPM + 1)
## ## ann.omega.original <- ann.omega

## ## load log2(TPM + 1)
## log2.omega <- readRDS(paste0(opt$scatteroutput, "/", SAMPLE, "_MAST_log2TPM.RDS"))

if( grepl("2kG.library|Perturb_2kG_dup4", SAMPLE) ) {
    df <- merge(log2.omega %>% as.data.frame %>% mutate(long.CBC = rownames(.)), meta_data, by="long.CBC")
} else {
    df <- merge(log2.omega %>% as.data.frame %>% mutate(CBC = rownames(.)), meta_data, by="CBC", all.y=T)
}


## gene.here <- "MESDC1"
gene.list <- gene.ary
num.ptb <- length(gene.list)
cat(paste0("number of genes: ", num.ptb, "\n"))

## MAST.list <- vector("list", num.ptb)
MAST.list <- mclapply(1:num.ptb, function(i) {
    gene.here <- gene.list[i]
    out <- tryCatch(
    { suppressWarnings({ 
        totest.df <- df %>%
            mutate(Gene = gsub("safe-targeting","negative-control", Gene)) %>% ## combine safe targeting and negative control guide to call them "negative-control"
            subset(Gene %in% c("negative-control", gene.here)) ## take a small slice of data (one perturbation one control)

        scaRaw <- FromMatrix(totest.df %>% select(-any_of(colnames(meta_data))) %>% t)
        colData(scaRaw)$perturb_status <- totest.df %>% pull(Gene)
        colData(scaRaw)$lane <- totest.df %>% 
                          pull(sample)
        cond <- factor(colData(scaRaw)$perturb_status)
        cond <- relevel(cond, "negative-control")
        colData(scaRaw)$condition <- cond
        
        format.zlm.result <- function(zlmCond) {
            condition.str <- paste0('condition', gene.here)
            summaryCond <- summary(zlmCond, doLRT=condition.str, parallel = T)

            summaryDt <- summaryCond$datatable
            fcHurdle <- merge(summaryDt[contrast==condition.str & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                              summaryDt[contrast==condition.str & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') %>% #logFC coefficients
                mutate(fdr:=p.adjust(`Pr(>Chisq)`, 'fdr'),
                       perturbation = gene.here)
            return(fcHurdle)
        }

        fcHurdle <- do.call(rbind, lapply(1:nrow(test.cmd.df), function(i) {
            zlmCond <- eval(parse(text = paste0("zlmCond <- zlm(", test.cmd.df$zlm.model.command[i], ", scaRaw)")))
            fcHurdle.here  <- format.zlm.result(zlmCond) %>% mutate(zlm.model.name = test.cmd.df$zlm.model.name[i])
            return(fcHurdle.here)
        }))
        cat(gene.here)
        
        return(fcHurdle)
    })},
    warning = function(cond) {
        return(fcHurdle)
    },
    error = function(cond) {
        return(data.frame(primerid = NA,
                          `Pr(>Chisq)`= NA,
                          coef = NA,
                          ci.hi = NA,
                          ci.lo = NA,
                          fdr = NA,
                          perturbation = gene.here, check.names=F,
                          zlm.model.name=NA)
               )
    }
    )
    ## MAST.list[[i]] <- fcHurdle %>% mutate(perturbation = gene.here)
}, mc.cores = max(1, floor(availableCores() - 1))) ## 64G, K=35, 3 perturbations took 30 minutes


MAST.df <- do.call(rbind, MAST.list) %>%
    ## group_by(zlm.model.name) %>%
    ## mutate(fdr.across.ptb = p.adjust(`Pr(>Chisq)`, method='fdr')) %>%
    as.data.frame

write.table(MAST.df, file=paste0(SCATTEROUTDIR, "/", SAMPLE, "_MAST_DEtopics_Group", opt$scatter.gene.group, ".txt"), quote=F, row.names=F, sep="\t")

## MAST.df <- read.delim(paste0(OUTDIRSAMPLE, "/", SAMPLE, "_MAST_DEtopics.txt"), stringsAsFactors=F)

## ## end of 220222 scratch


## perturb_status <- meta_data %>%
##     as.data.frame %>%
##     mutate(Gene = gsub("-TSS2$", "", Gene.full.name)) %>%
##     select(Gene)
## colData(scaRaw)$perturb_status <- perturb_status %>% pull(Gene)
## colData(scaRaw)$sample <- meta_data %>%
##                   separate(col="CBC", into=c("CBC", "sample"), sep="-scRNAseq_2kG_") %>%
##                   pull(sample)

## cond <- factor(gsub("safe-targeting", "negative-control", perturb_status$Gene))
## cond <- relevel(cond, "negative-control")
## colData(scaRaw)$condition <- cond
## zmCond <- zlm(~condition + sample, scaRaw)

## save.image("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/zmCond.RDS")
