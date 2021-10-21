## calculate UMAP given already available Seurat Objects
## Helen Kang
## 210519

.libPaths("/home/groups/engreitz/Software/R_3.6.1")

suppressPackageStartupMessages(library(optparse))

 
option.list <- list(
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/heart_atlas/2105_FT005_Analysis/outputs/", help="Output directory"),
    make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/heart_atlas/2105_FT005_Analysis/outputs/", help="Output directory"),
    make_option("--sampleName",type="character",default="FT005_gex_new_pipeline", help="Sample name"),
    make_option("--project",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/heart_atlas/2105_FT005_Analysis/",help="Project Directory"),
    make_option("--inputSeuratObject", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/heart_atlas/2105_FT005_Analysis/outputs/FT005_gex/withUMAP.SeuratObject.RDS", help="Path to the Seurat Object"),
    make_option("--compareSeuratObject", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/heart_atlas/2105_FT005_Analysis/outputs/FT005_gex_new_pipeline/withUMAP.SeruatObject.RDS", help="Path to the Seurat Object"),
    make_option("--maxMt", type="numeric", default=50, help="filter out cells with percent mitochondrial gene higher than this threhsold"),
    make_option("--maxCount", type="numeric", default=25000, help="filter out cells with UMI count more than this threshold"),
    make_option("--minUniqueGenes", type="numeric", default=0, help="filter out cells with unique gene detected less than this threshold"),
    make_option("--UMAP.resolution", type="numeric", default=0.06, help="UMAP resolution. The default is 0.06")
)
opt <- parse_args(OptionParser(option_list=option.list))

library(SeuratObject)
library(Seurat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggrepel))

source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")

mytheme <- theme_classic() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5)) 


#######################################################################
## Constants
PROJECT=opt$project
OUTDIR=opt$outdir
FIGDIR=opt$figdir
SAMPLE=opt$sampleName
# OUTDIRSAMPLE=paste0(OUTDIR,"/",SAMPLE,"/")
FIGDIRSAMPLE=paste0(FIGDIR,"/",SAMPLE,"/")
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n=100)
# create dir if not already
# check.dir <- c(OUTDIR, OUTDIRSAMPLE, FIGDIR, FIGDIRSAMPLE)
check.dir <- c(OUTDIR, FIGDIR, FIGDIRSAMPLE)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x) }))



## function for calculating UMAP
calcUMAP <- function(s) {
    s <- SCTransform(s)
    s <- RunPCA(s, verbose = FALSE)
    s <- FindNeighbors(s, dims = 1:10)
    s <- FindClusters(s, resolution = opt$UMAP.resolution)
    s <- RunUMAP(s, dims = 1:10)
}



## load data
s <- readRDS(opt$inputSeuratObject)
s <- calcUMAP(s)

## Choose filters, use cells with Good_singlet, and filter MT/RP
s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
s[["percent.ribo"]] <- PercentageFeatureSet(s, pattern = "^RPS|^RPL")


## plot QC
s.meta <- SeuratObject::FetchData(s, colnames(s[[]]))  ## This weird Seurat syntax gets the list of all metadata vars and then fetches a matrix of the data
# plotSingleCellStats(s.meta, mtMax=NULL, nCountMax=NULL, paste0(FIGDIRSAMPLE,"/QC.single.cell.stats.UMAPres.", opt$UMAP.resolution, ".pdf"))

# saveRDS(s, paste0(OUTDIRSAMPLE,"/", SAMPLE, ".withUMAP.", opt$UMAP.resolution, ".SeuratObject.RDS")) 
# saveRDS(s, paste0(OUTDIRSAMPLE,"/", SAMPLE, ".withUMAP_SeuratObject.RDS")) 
saveRDS(s, paste0(OUTDIR,"/", SAMPLE, ".withUMAP_SeuratObject.RDS")) 

# plotUmap(s, title="", split=NULL, paste0(FIGDIRSAMPLE,"/QC.unfiltered.umap.pdf"))
# # plotUmap(s, title="", split=NULL, paste0(FIGDIR,"/",SAMPLE, "_QC.unfiltered.umap.pdf"))


# ## Filters
# maxMt <- opt$maxMt
# maxCount <- opt$maxCount
# minUniqueGenes <- opt$minUniqueGenes

# sf <- subset(s, subset = percent.mt < maxMt & nCount_RNA < maxCount & nFeature_RNA > minUniqueGenes)

# ## plot QC
# sf.meta <- SeuratObject::FetchData(sf, colnames(sf[[]]))  ## This weird Seurat syntax gets the list of all metadata vars and then fetches a matrix of the data
# plotSingleCellStats(sf.meta, mtMax=NULL, nCountMax=NULL, paste0(FIGDIRSAMPLE,"/QC.filtered.single.cell.stats.pdf"))

# plotUmap(sf, title="", split=NULL, paste0(FIGDIRSAMPLE,"/QC.filtered.umap.pdf"))

# saveRDS(sf, paste0(OUTDIRSAMPLE,"/", SAMPLE, ".withUMAP.", opt$UMAP.resolution, ".filtered.SeuratObject.RDS")) 


