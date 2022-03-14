## Helen Kang
## run MAST with batch correction
## 220217




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
conflict_prefer("select", "dplyr")
conflict_prefer("melt", "reshape2")
conflict_prefer("filter", "dplyr")

packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
              "tidyr", "grid", "gtable", "gridExtra","ggrepel",#"ramify",
              "ggpubr","gridExtra",
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
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--outdirsample", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/", help="path to cNMF analysis results"), ## or for 2n1.99x: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/Perturb_2kG_dup4/K60/threshold_0_2/"

    ## script dir
    make_option("--scriptdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts/", help="location for this script and functions script")

)
opt <- parse_args(OptionParser(option_list=option.list))

k <- opt$K.val 
SAMPLE <- opt$sampleName
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
## OUTDIRSAMPLE <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/"
OUTDIRSAMPLE <- opt$outdirsample

INPUTDIR <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/inputs/"
OUTDIR <- OUTDIRSAMPLE
check.dir <- c(INPUTDIR, OUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

mytheme <- theme_classic() + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, face = "bold"))
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)
p.adjust.thr <- 0.1



##########################################################################################
## load data
## load known CAD genes
CAD.genes <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/data/known_CAD_gene_set.txt", header=F, stringsAsFactors=F) %>% as.matrix %>% as.character

## load topic model results
cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
}

file.name <- paste0(OUTDIRSAMPLE,"/cNMFAnalysis.",SUBSCRIPT,".RData")
print(file.name) 
if(file.exists((file.name))) { 
    print(paste0("loading ",file.name))
    load(file.name) 
}

## add UMI / cell information
umiPerCell.df <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/outputs/UMIPerCell.txt", stringsAsFactors=F, row.names=1) %>% mutate(long.CBC = rownames(.)) %>%
    separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
    separate(col="CBC", into=c("CBC", "sample"), sep="-scRNAseq_2kG_", remove=F)

sample.to.10X.lane <- data.frame(sample_num = 1:20,
                                 sample = umiPerCell.df$sample %>% unique %>% sort)
umiPerCell.df <- merge(umiPerCell.df, sample.to.10X.lane, by="sample") %>%
    mutate(CBC_10x = paste0(CBC, "-", sample_num))  

## add guide / cell information
guidePerCell.df <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/inputs/_UPDATED_ALL_SAMPLES_dup4_NON_NA_CALLS_FOR_EA_CBC_FULL_INFO.txt", stringsAsFactors=F) 

## add genes detected / cell information
geneDetectedPerCell.df <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/outputs/geneDetectedPerCell.txt", stringsAsFactors=F) %>%
    mutate(long.CBC = rownames(.)) %>%
    separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
    separate(col="CBC", into=c("CBC", "sample"), sep="-scRNAseq_2kG_", remove=F) %>%
    merge(sample.to.10X.lane, by="sample") %>%
    mutate(CBC_10x = paste0(CBC, "-", sample_num))  

## get metadata from log2.X.full rownames
if(SAMPLE != "Perturb_2kG_dup4") {
    meta_data <- omega %>%
        rownames %>%
        as.data.frame %>%
        `colnames<-`("long.CBC") %>%
        separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
        ## separate(col="CBC", into=c("CBC", "sample"), sep="-", remove=F) %>%
        separate(col="CBC", into=c("CBC", "sample"), sep="-scRNAseq_2kG_", remove=F) %>%
        mutate(Gene = gsub("-TSS2$", "", Gene.full.name),
               CBC = gsub("RHOA-and-", "", CBC)) %>%
        as.data.frame
    sample.to.10X.lane <- data.frame(sample_num = 1:20,
                                     sample = meta_data$sample %>% unique %>% sort)

    meta_data <- merge(meta_data, sample.to.10X.lane, by="sample") %>%
        mutate(CBC_10x = paste0(CBC, "-", sample_num)) %>%
        merge(guidePerCell.df %>% select(CBC_10x, guides_per_cbc, max_umi_ct), by="CBC_10x") %>%
        merge(umiPerCell.df, by="long.CBC") %>%
        merge(geneDetectedPerCell.df, by.x="long.CBC", by.y=0)
} else {
    meta_data <- omega %>%
        rownames %>%
        as.data.frame %>%
        `colnames<-`("long.CBC") %>%
        separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC_10x"), sep=":", remove=F) %>%
        separate(col="CBC_10x", into=c("CBC", "sample"), sep="-", remove=F) %>%
        mutate(Gene = gsub("-TSS2$", "", Gene.full.name),
               CBC = gsub("RHOA-and-", "", CBC)) %>%
        as.data.frame

    ## meta_data <- meta_data %>%
    ##     merge(guidePerCell.df %>% select(CBC_10x, guides_per_cbc, max_umi_ct), by="CBC_10x") %>%
    ##     merge(umiPerCell.df %>% select(CBC_10x, num.UMI), by="CBC_10x") %>%
    ##     merge(geneDetectedPerCell.df %>% select(CBC_10x, num.geneDetected), by="CBC_10x")
    meta_data$sample <- as.factor(meta_data$sample)
}    
 

## load model fitting fomulas
test.cmd.df <- read.delim(paste0(INPUTDIR, "MAST_model_formulas.txt"), stringsAsFactors=F)


## ## follow MAST tutorial to create the approapriate objects
## ## create ann.omega
## if ( !( "ann.omega" %in% ls()) ) {
##     ann.omega <- omega %>%
##         as.data.frame %>%
##         mutate(long.CBC = rownames(.)) %>%
##         separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
##         mutate(Gene = gsub("-TSS2$", "", Gene.full.name))
## }                        
## meta_data <- ann.omega[,61:65]

## 220222 scratch
if ( !( "ann.omega" %in% ls()) ) {
    ann.omega <- omega %>%
        as.data.frame %>%
        mutate(long.CBC = rownames(.)) %>%
        separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>%
        mutate(Gene = gsub("-TSS2$", "", Gene.full.name))
} else {                        
    ann.omega <- ann.omega %>%
        mutate(Gene = gsub("-TSS2$", "", Gene.full.name)) ## remove TSS2 annotation
}
omega.tpm <- ann.omega[,1:k] %>% as.matrix %>% apply(2, function(x) x / sum(x) * 1000000) ## convert to TPM
log2.omega <- (omega.tpm + 1) %>% log2 ## log2(TPM + 1)
## ann.omega.original <- ann.omega
df <- merge(log2.omega %>% as.data.frame %>% mutate(long.CBC = rownames(.)), meta_data, by="long.CBC")


gene.here <- "MESDC1"
gene.list <- df$Gene %>% unique
num.ptb <- length(gene.list)

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
        if(SAMPLE != "Perturb_2kG_dup4") {
        colData(scaRaw)$umiPerCell <- totest.df %>% pull(num.UMI)
        colData(scaRaw)$topGuideUMI <- totest.df %>% pull(max_umi_ct)
        colData(scaRaw)$geneDetectedPerCell <- totest.df %>% pull(num.geneDetected)
        }
        
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

        ## scratch 220303
 
        ## zlmCond <- zlm(~condition + lane, scaRaw)
        ## fcHurdle.list <- vector("list", nrow(test.cmd.df)) ## create storage list for the different models we want to fit the data to
        ## system.time(
        fcHurdle <- do.call(rbind, lapply(1:nrow(test.cmd.df), function(i) {
            zlmCond <- eval(parse(text = paste0("zlmCond <- zlm(", test.cmd.df$zlm.model.command[i], ", scaRaw)")))
            fcHurdle.here  <- format.zlm.result(zlmCond) %>% mutate(zlm.model.name = test.cmd.df$zlm.model.name[i])
            ## write.table(fcHurdle.here, paste0("./outputs/", gene.here, "_", test.cmd.df$zlm.model.name[i], ".txt"), quote=F, sep="\t", row.names=F)
            return(fcHurdle.here)
        }))
        print(gene.here)
        ## write.table(fcHurdle, paste0("./outputs/", gene.here, "_DEtopics.txt"), quote=F, sep="\t", row.names=F)
        ## )
        ## fcHurdle <- do.call(rbind, fcHurdle.list)

        ## system.time({
        ## fcHurdle.list <- vector("list", nrow(test.cmd.df)) ## create storage list for the different models we want to fit the data to
        ## for(i in 1:nrow(test.cmd.df)) {
        ## zlmCond <- eval(parse(text = paste0("zlmCond <- zlm(", test.cmd.df$zlm.model.command[i], ", scaRaw)")))
        ##     fcHurdle.list[[i]] <- format.zlm.result(zlmCond) %>% mutate(zlm.model.name = test.cmd.df$zlm.model.name[i])
        ## }
        ## fcHurdle <- do.call(rbind, fcHurdle.list)
        ## })
        
        return(fcHurdle)
    })},
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
}, mc.cores = detectCores() - 2)
MAST.df <- do.call(rbind, MAST.list) %>%
    group_by(zlm.model.name) %>%
    mutate(fdr.across.ptb = p.adjust(`Pr(>Chisq)`, method='fdr')) %>%
    as.data.frame
write.table(MAST.df, file=paste0(OUTDIRSAMPLE, "/", SAMPLE, "_MAST_DEtopics.txt"), quote=F, row.names=F, sep="\t")

MAST.df <- read.delim(paste0(OUTDIRSAMPLE, "/", SAMPLE, "_MAST_DEtopics.txt"), stringsAsFactors=F)

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
