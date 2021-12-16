## Make Seurat Object
## Helen Kang
## 210613

## .libPaths("/home/groups/engreitz/Software/R_3.6.1")

suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/heart_atlas/2106_FT007_Analysis/outputs/", help="Output directory"),
    make_option("--datadir",type="character", default="/oak/stanford/groups/engreitz/Users/kangh/process_sequencing_data/210611_FT007_CM_CMO/gex_FT007_50k/outs/filtered_feature_bc_matrix/", help="Data directory"),
    # make_option("--project",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/heart_atlas/2106_FT007_Analysis/",help="Project Directory"),
    make_option("--sampleName",type="character",default="gex_FT007_50k", help="Sample name"),
    make_option("--recompute", type="logical", default=F, help="T for recomputing UMAP from 10x count matrix")
)
opt <- parse_args(OptionParser(option_list=option.list))

## opt$outdir <- "/oak/stanford/groups/engreitz/Users/katherine/perturb-seq/lowmoi_cnmf/analysis/data/"
## opt$datadir <- "/oak/stanford/groups/engreitz/Users/katherine/perturb-seq/data/pilot_low_moi"
## opt$sampleName <- "lowmoi"

## library(SeuratObject)
library(Seurat)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggrepel))

mytheme <- theme_classic() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5))

## source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")


#######################################################################
## Constants
# PROJECT=opt$project
DATADIR=opt$datadir
OUTDIR=opt$outdir
# FIGDIR= paste0(PROJECT, "/figures/")
SAMPLE=opt$sampleName
# OUTDIRSAMPLE=paste0(OUTDIR,"/",SAMPLE,"/")
# FIGDIRSAMPLE=paste0(FIGDIR,SAMPLE,"/")
palette = colorRampPalette(c("#38b4f7", "white", "red"))(n=100)
# create dir if not already
check.dir <- c(OUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))


outSeuratObjectPath <- paste0(OUTDIR,"/", SAMPLE,".SeuratObject.RDS")
if(opt$recompute | !file.exists(outSeuratObjectPath)) {
#######################################################################
    ## Read 10X Matrix
    features <- read.delim(file=paste0(DATADIR, "/features.tsv.gz"), header=F, stringsAsFactors=F)
    if(features %>% ncol == 1) {
        my.data <- Read10X(data.dir = DATADIR, gene.column=1) ## for manually assembled features (genes) that only has one column
    } else {
        my.data <- Read10X(data.dir = paste0(DATADIR)) # two matrices: "Gene Expression", "Peaks"
    }
    if(class(my.data) == "list") {
			message("loading Gene Expression Matrix from Multiome data")
			gex.mtx <- my.data[["Gene Expression"]]
			} else {
			message("loading Gene Expression Matrix")
			gex.mtx <- my.data
		}


    ## Create Seurat Object
    message("creating Seurat object")
    s <- CreateSeuratObject(
        counts = gex.mtx,
        project = SAMPLE
    )

    ## Choose filters, use cells with Good_singlet, and filter MT/RP
    s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
    s[["percent.ribo"]] <- PercentageFeatureSet(s, pattern = "^RPS|^RPL")

    # ## plot QC
    # s.meta <- SeuratObject::FetchData(s, colnames(s[[]]))  ## This weird Seurat syntax gets the list of all metadata vars and then fetches a matrix of the data
    # plotSingleCellStats(s.meta, mtMax=NULL, nCountMax=NULL, paste0(FIGDIRSAMPLE,"/QC.single.cell.stats.pdf"))

    ## remove non-protein coding genes and genes detected in fewer than 10 cells
    tokeep <- which(!(grepl("^LINC|^[A-Za-z][A-Za-z][0-9][0-9][0-9][0-9][0-9][0-9]\\.", s %>% rownames))) ## remove antisense genes, divergent transcripts, non coding genes, and intronic transcripts
    s.subset <- s[tokeep,]
    s.subset <- subset(s.subset, subset= nCount_RNA > 200 & nFeature_RNA > 200) # remove cells with less than 200 UMIs and less than 200 genes
    ## tokeep <- which(s.subset@assays$RNA@counts %>% apply(1, sum) > 10) # keep genes detected in more than 10 UMIs
    tokeep <- which(s.subset@assays$RNA@counts %>% apply(1, function(x) ((x > 0) %>% as.numeric %>% sum > 10)))
    s.subset <- s.subset[tokeep,]

    s <- s.subset

    saveRDS(s, outSeuratObjectPath)

} else { message(paste0("File already exists at ", outSeuratObjectPath)) }
