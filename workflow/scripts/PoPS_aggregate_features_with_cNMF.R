## aggregate PoPS raw features txt files to a single RDS file
## Helen Kang
## 211118 referenced 211101 log.R

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
              "tidyr", "textshape","readxl", "AnnotationDbi")
xfun::pkg_attach(packages)
conflict_prefer("select", "dplyr")


option.list <- list(
    make_option("--feature.RDS", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/data/features/pops_features_raw/"),
    make_option("--cNMF.features", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211101_normalized_features/outputs/", help="cNMF features in ENSGID to add to all features"),
    make_option("--output", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211101_normalized_features/outputs/", help="output directory"),
    make_option("--prefix", type="character", default="CAD_aug6_cNMF60", help="use a format of MAGMA_{with, without}cNMF to specify which features are included")
)
opt <- parse_args(OptionParser(option_list=option.list))


## ## for sdev
## opt$feature.dir <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/211101_20sample_snakemake/pops/features/pops_features_raw"
## opt$output <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/pops/"
## opt$prefix <- "CAD_aug6_cNMF30"


OUTDIR=opt$output
PREFIX=opt$prefix



## load cNMF features
features <- read.delim(opt$cNMF.features, stringsAsFactors=F)

## load external features RDS
all.features <- readRDS(opt$feature.RDS) ## features without cNMF

## combine cNMF features and external features
all.features.cNMF <- merge(all.features, features, by="ENSGID") ## all features with cNMF
saveRDS(all.features.cNMF, file=paste0(OUTDIR, "/full_features_", PREFIX, ".RDS"))
write.table(all.features.cNMF, file=paste0(OUTDIR, "/full_features_", PREFIX, ".txt", row.names=F, quote=F, sep="\t"))
