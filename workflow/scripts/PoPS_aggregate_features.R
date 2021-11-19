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
    make_option("--prefix", type="character", default="CAD_aug6_cNMF60", help="use a format of MAGMA_{with, without}cNMF to specify which features are included")
)
opt <- parse_args(OptionParser(option_list=option.list))
OUTDIR=opt$output
PREFIX=opt$prefix



# ## load cNMF features
# features <- read.delim("/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210831_PoPS/data/features/pops_features_raw/topic.zscore.ensembl.scaled_k_60.dt_0_2.txt", stringsAsFactors=F)

## get files in feature directory
feature.file.paths <- dir(opt$feature.dir)

## load all features
num_feature_files = length(feature.file.paths)
all.features.list <- vector("list", num_feature_files)
for (feature.raw.index in 1:num_feature_files) {
    all.features.list[[feature.raw.index]] <- read.delim(feature.file.paths[[feature.raw.index]], stringsAsFactors=F)
}

all.features <- Reduce(function(x,y) merge(x,y,by="ENSGID"), all.features.list) ## features without cNMF
all.features.cNMF <- merge(all.features, features, by="ENSGID") ## all features with cNMF
saveRDS(all.features, file=paste0(OUTDIR, "/full_features_", PREFIX, ".RDS"))
write.table(all.features, file=paste0(OUTDIR, "/full_features_", PREFIX, ".txt", row.names=F, quote=F, sep="\t")
