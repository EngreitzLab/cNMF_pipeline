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
    make_option("--feature.dir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/data/features/pops_features_raw/"),
    make_option("--output", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/211101_normalized_features/outputs/", help="output directory"),
    make_option("--magma_prefix", type="character", default="CAD_aug6", help="MAGMA prefix"),
    make_option("--sampleName", type="character", default="2kG.library", help="Sample name for this run to output")
)
opt <- parse_args(OptionParser(option_list=option.list))


## ## for sdev
## opt$feature.dir <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/211101_20sample_snakemake/pops/features/pops_features_raw"
## opt$output <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/scRNAseq_2kG_11AMDox_1/pops/"
## opt$prefix <- "CAD_aug6_cNMF30"


OUTDIR=opt$output
PREFIX=opt$prefix
SAMPLE=opt$sampleName


## ## load cNMF features
## features <- read.delim(opt$, stringsAsFactors=F)

## get files in feature directory
feature.file.paths <- dir(opt$feature.dir)

## load all features
num_feature_files = length(feature.file.paths)
all.features.list <- vector("list", num_feature_files)
for (feature.raw.index in 1:num_feature_files) {
    all.features.list[[feature.raw.index]] <- read.delim(paste0(opt$feature.dir, "/", feature.file.paths[[feature.raw.index]]), stringsAsFactors=F)
}

all.features <- Reduce(function(x,y) merge(x,y,by="ENSGID"), all.features.list) ## features without cNMF
## all.features.cNMF <- merge(all.features, features, by="ENSGID") ## all features with cNMF
print(paste0("Saving all.features to ", OUTDIR, " as an RDS file"))
saveRDS(all.features, file=paste0(OUTDIR, "/full_external_features_", PREFIX, "_", SAMPLE, ".RDS"))
print(paste0("Saving all.features to ", OUTDIR, " as a text file"))
write.table(all.features, file=paste0(OUTDIR, "/full_external_features_", PREFIX, "_", SAMPLE, ".txt"), row.names=F, quote=F, sep="\t")
print("done!")