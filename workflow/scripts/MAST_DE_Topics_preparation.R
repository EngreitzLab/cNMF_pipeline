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
    ## make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
    ## make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
    make_option("--outdirsample", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/", help="path to cNMF analysis results"), ## or for 2n1.99x: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/Perturb_2kG_dup4/K60/threshold_0_2/"
    make_option("--scatteroutput", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/MAST/", help="path to gene breakdown table output"),
    make_option("--numCtrl", type="numeric", default=5000, help="number of control cells to use for MAST")
    
    ## ## script dir
    ## make_option("--scriptdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts/", help="location for this script and functions script")

)
opt <- parse_args(OptionParser(option_list=option.list))


## ## sdev debug K562 gwps
## opt$sampleName <- "WeissmanK562gwps"
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/data/K562_gwps_raw_singlecell_01_metadata.txt"
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K35/threshold_0_2/"
## opt$scatteroutput <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/MAST/K35/threshold_0_2/"
## opt$scatter.gene.group <- 496
## opt$scriptdir <- "/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts"
## opt$K.val <- 35


k <- opt$K.val 
SAMPLE <- opt$sampleName
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
## SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
## OUTDIRSAMPLE <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K31/threshold_0_2/"
OUTDIRSAMPLE <- opt$outdirsample
SCATTEROUTDIR <- opt$scatteroutput
SCATTERINDEX <- opt$scatter.gene.group

INPUTDIR <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220217_MAST/inputs/"
OUTDIR <- OUTDIRSAMPLE
check.dir <- c(INPUTDIR, OUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))



##########################################################################################
## load data

## load topic model results
cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
}

barcode.names <- read.delim(opt$barcode.names, stringsAsFactors=F)
## separate out control cells
omega.ctrl.index <- barcode.names %>% mutate(rowindex = 1:n()) %>% filter(Gene == "negative-control") %>% pull(rowindex)
omega.ptb.index <- barcode.names %>% mutate(rowindex = 1:n()) %>% filter(Gene != "negative-control") %>% pull(rowindex)
omega.ctrl <- omega[omega.ctrl.index,]
omega.ptb <- omega[omega.ptb.index,]
barcode.names.ctrl <- barcode.names[omega.ctrl.index,]
barcode.names.ptb <- barcode.names[omega.ptb.index,]
## randomly subset to 5000 cells
ctrl.subset.index <- sample(1:length(omega.ctrl.index), min(length(omega.ctrl.index), opt$numCtrl), replace=FALSE)
omega.ctrl.subset <- omega.ctrl[ctrl.subset.index,]
barcode.names.ctrl.subset <- barcode.names.ctrl[ctrl.subset.index,]

omega.new <- rbind(omega.ptb, omega.ctrl.subset)
barcode.names.subset <- rbind(barcode.names.ptb, barcode.names.ctrl.subset)

omega.tpm <- omega.new %>% as.matrix %>% apply(2, function(x) x / sum(x) * 1000000) ## convert to TPM
log2.omega <- (omega.tpm + 1) %>% log2 ## log2(TPM + 1)

## output log2(TPM + 1)
save(log2.omega, barcode.names.subset, file=paste0(opt$scatteroutput, "/", SAMPLE, "_MAST_log2TPM_barcodes.RDS"))
