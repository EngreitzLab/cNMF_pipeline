## Make Seurat Object
## Helen Kang
## 210613

.libPaths("/home/groups/engreitz/Software/R_3.6.1")

suppressPackageStartupMessages(library(optparse))

option.list <- list(
    # make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/heart_atlas/2106_FT007_Analysis/outputs/", help="Output directory"),
    # make_option("--datadir",type="character", default="/oak/stanford/groups/engreitz/Users/kangh/process_sequencing_data/210611_FT007_CM_CMO/gex_FT007_50k/outs/filtered_feature_bc_matrix/", help="Data directory"),
    # make_option("--project",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/heart_atlas/2106_FT007_Analysis/",help="Project Directory"),
    # make_option("--sampleName",type="character",default="gex_FT007_50k", help="Sample name"),
    make_option("--inputSeuratObject", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/data/FT010_fresh_3min.SeuratObject.RDS", help="Path to the Seurat Object"),
    make_option("--output_h5ad", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/data/FT010_fresh_3min.h5ad"),
    make_option("--output_gene_name_txt", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/data/FT010_fresh_3min.h5ad.all.genes.txt")
    # make_option("--recompute", type="logical", default=F, help="T for recomputing UMAP from 10x count matrix")
)
opt <- parse_args(OptionParser(option_list=option.list))

suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
# suppressPackageStartupMessages(library(readxl))
# suppressPackageStartupMessages(library(ggrepel))

# mytheme <- theme_classic() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5))

# source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")


#######################################################################
## Constants
# PROJECT=opt$project
# DATADIR=opt$datadir
# OUTDIR=opt$outdir
# FIGDIR= paste0(PROJECT, "/figures/")
# SAMPLE=opt$sampleName
# OUTDIRSAMPLE=paste0(OUTDIR,"/",SAMPLE,"/")
# FIGDIRSAMPLE=paste0(FIGDIR,SAMPLE,"/")
# palette = colorRampPalette(c("#38b4f7", "white", "red"))(n=100)
# # create dir if not already
# check.dir <- c(OUTDIR, OUTDIRSAMPLE)
# invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x) }))


## convert Seurat Object to h5ad
anndata <- import("anndata", convert = FALSE) ## load AnnData module

## load Seruat Object
s <- readRDS(opt$inputSeuratObject)



adata <- anndata$AnnData(
    X = t(GetAssayData(object = s) %>% as.matrix),
    obs = data.frame(s@meta.data),
    var = s %>% rownames %>% as.data.frame %>% `rownames<-`(s %>% rownames) %>% `colnames<-`("Gene")
)

anndata$AnnData$write(adata, opt$output_h5ad)

## write gene names
gene.names <- s %>% rownames
write.table(gene.names, file=opt$output_gene_name_txt, col.names=F, row.names=F, quote=F, sep="\t")

