## Helen Kang
## Script to gather all MAST runs
## 230128


library(conflicted)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("first", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("melt", "reshape2")
conflict_prefer("filter", "dplyr")

packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
              "tidyr") #, "grid", "gtable", "gridExtra","ggrepel",#"ramify",
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
    make_option("--scatteroutput", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/MAST/K35/threshold_0_2/", help="path to gene breakdown table output"),
    make_option("--total.scatter.gene.group", type="numeric", default=50, help="Total number of groups for running MAST"),
    
    ## script dir
    make_option("--scriptdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts/", help="location for this script and functions script")

)
opt <- parse_args(OptionParser(option_list=option.list))


## ## sdev debug K562 gwps
## opt$sampleName <- "WeissmanK562gwps"
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/data/K562_gwps_raw_singlecell_01_metadata.txt"
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K35/threshold_0_2/"
## opt$scatteroutput <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/230104_snakemake_WeissmanLabData/top2000VariableGenes/WeissmanK562gwps/MAST/K35/threshold_0_2/"
## opt$total.scatter.gene.group <- 494
## opt$scriptdir <- "/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts"
## opt$K.val <- 35


k <- opt$K.val 
SAMPLE <- opt$sampleName
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
## SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
OUTDIRSAMPLE <- opt$outdirsample
SCATTEROUTDIR <- opt$scatteroutput
NUMGROUPS <- opt$total.scatter.gene.group



OUTDIR <- OUTDIRSAMPLE
check.dir <- c(OUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

## load data
MAST.df.list <- vector("list", NUMGROUPS)
for (i in 1:NUMGROUPS) {
    MAST.df.list[[i]] <- read.delim(paste0(SCATTEROUTDIR, "/", SAMPLE, "_MAST_DEtopics_Group", i, ".txt"), stringsAsFactors=F, check.names=F)
}
MAST.df <- do.call(rbind, MAST.df.list) %>%
    group_by(zlm.model.name) %>%
    mutate(fdr.across.ptb = p.adjust(`Pr(>Chisq)`, method='fdr')) %>%
    as.data.frame

write.table(MAST.df, file=paste0(OUTDIRSAMPLE, "/", SAMPLE, "_MAST_DEtopics.txt"), quote=F, row.names=F, sep="\t")

