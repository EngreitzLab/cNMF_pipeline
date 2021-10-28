
## Helen Kang
## Topic Model Analysis Only (no plot output)
## 210503

## .libPaths("/home/groups/engreitz/Software/R_3.6.1")


## packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
##               "CountClust", "Hmisc", "tidyr", "grid", "gtable", "gridExtra","ggrepel","ramify",
##               "GGally","RNOmni","usedist","ggpubr","gridExtra","GSEA",
##               "org.Hs.eg.db","limma","clusterProfiler","fgsea", "conflicted",
##               "cluster","textshape","readxl", "IsoplotR", "wesanderson", 
##               "ggdist", "gghalves", "Seurat", "writexl")
## library(textshape)
## library(readxl)
## # library(Seurat)
packages <- c("optparse", "data.table", "reshape2", "fgsea", "conflicted", "readxl", "writexl", "org.Hs.eg.db", "tidyr", "dplyr")
xfun::pkg_attach(packages)
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")


## source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")

option.list <- list(
  make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/figures/all_genes/", help="Figure directory"),
  make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/all_genes/", help="Output directory"),
  # make_option("--olddatadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/", help="Input 10x data directory"),
  make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/", help="Input 10x data directory"),
  make_option("--topic.model.result.dir", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/210707_snakemake_maxParallel/all_genes_acrossK/2kG.library/", help="Topic model results directory"),
  make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
  # make_option("--sep", type="logical", default=F, help="Whether to separate replicates or samples"),
  make_option("--K.list", type="character", default="2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25", help="K values available for analysis"),
  make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
  make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
  make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
  make_option("--ABCdir",type="character", default="/oak/stanford/groups/engreitz/Projects/ABC/200220_CAD/ABC_out/TeloHAEC_Ctrl/Neighborhoods/", help="Path to ABC enhancer directory"),
  make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
  make_option("--raw.mtx.dir",type="character",default="stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/data/no_IL1B_filtered.normalized.ptb.by.gene.mtx.filtered.txt", help="input matrix to cNMF pipeline"),
  make_option("--raw.mtx.RDS.dir",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.RDS", help="input matrix to cNMF pipeline"), # the first lane: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.expandedMultiTargetGuide.RDS"
  make_option("--subsample.type", type="character", default="", help="Type of cells to keep. Currently only support ctrl"),
  # make_option("--barcode.names", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/barcodes.tsv", help="barcodes.tsv for all cells"),
  make_option("--reference.table", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/210702_2kglib_adding_more_brief_ca0713.xlsx"),
  
  ## fisher motif enrichment
  ## make_option("--outputTable", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.k_14.df_0_2.txt", help="Output directory"),
  ## make_option("--outputTableBinary", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210607_snakemake_output/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.binary.k_14.df_0_2.txt", help="Output directory"),
  ## make_option("--outputEnrichment", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210607_snakemake_output/outputs/no_IL1B/topic.top.100.zscore.gene.motif.fisher.enrichment.k_14.df_0_2.txt", help="Output directory"),
  make_option("--motif.promoter.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv", help="All promoter's motif matches"),
  make_option("--motif.enhancer.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv", help="All enhancer's motif matches specific to {no,plus}_IL1B"),
  make_option("--enhancer.fimo.threshold", type="character", default="1.0E-4", help="Enhancer fimo motif match threshold"),

  #summary plot parameters
  make_option("--test.type", type="character", default="per.guide.wilcoxon", help="Significance test to threshold perturbation results"),
  make_option("--adj.p.value.thr", type="numeric", default=0.1, help="adjusted p-value threshold"),
  make_option("--recompute", type="logical", default=F, help="T for recomputing statistical tests and F for not recompute")
  
)
opt <- parse_args(OptionParser(option_list=option.list))

## ## all genes directories (for sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/"
## opt$K.val <- 60

## ## debug ctrl
## opt$topic.model.result.dir <- "/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/210810_snakemake_ctrls/all_genes_acrossK/2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation"
## opt$sampleName <- "2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210810_snakemake_ctrls/figures/2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210810_snakemake_ctrls/analysis/2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation/all_genes/"
## opt$barcode.names <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210806_curate_ctrl_mtx/outputs/2kG.library.no.DE.gene.with.FDR.less.than.0.1.perturbation.barcodes.tsv"
## opt$K.val <- 60


## mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))

SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
DATADIR=opt$olddatadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
TMDIR=opt$topic.model.result.dir
# SEP=opt$sep
# K.list <- strsplit(opt$K.list,",") %>% unlist() %>% as.numeric()
k <- opt$K.val
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
FIGDIR=opt$figdir
FIGDIRSAMPLE=paste0(FIGDIR, "/", SAMPLE, "/K",k,"/")
FIGDIRTOP=paste0(FIGDIRSAMPLE,"/",SAMPLE,"_K",k,"_dt_", DENSITY.THRESHOLD,"_")
OUTDIRSAMPLE=paste0(OUTDIR, "/", SAMPLE, "/K",k,"/threshold_", DENSITY.THRESHOLD, "/")
FGSEADIR=paste0(OUTDIRSAMPLE,"/fgsea/")
FGSEAFIG=paste0(FIGDIRSAMPLE,"/fgsea/")

## subscript for files
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
# SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr

## ## directories for factor motif enrichment
## FILENAME=opt$filename


## ## modify motif.enhancer.background input directory ##HERE: perhaps do a for loop for all the desired thresholds (use strsplit on enhancer.fimo.threshold)
## opt$motif.enhancer.background <- paste0(opt$motif.enhancer.background, opt$enhancer.fimo.threshold, "/fimo.formatted.tsv")

 
# create dir if not already
check.dir <- c(OUTDIR, FIGDIR, paste0(FIGDIR,SAMPLE,"/"), paste0(FIGDIR,SAMPLE,"/K",k,"/"), paste0(OUTDIR,SAMPLE,"/"), OUTDIRSAMPLE, FIGDIRSAMPLE, FGSEADIR, FGSEAFIG)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)


######################################################################
## Load topic model results
    
cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
}


##########################################################################
## GSEA

## check if files already exist
check.file <- paste0(FGSEADIR,"/fgsea_all_pathways_df_", c("raw.score", "z.score"), "_", SUBSCRIPT.SHORT, ".RData") 
if(file.exists(check.file) %>% as.numeric() %>% sum() == length(check.file) & !opt$recompute ) {
    for (i in 1:length(check.file)) {
        print(paste0(check.file[i]), " already exists")
        # load(check.file[i])
    }
} else {

## import gene set annotation data
msigdb.all <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data//msigdb.v7.2.entrez.gmt")) # downloaded from https://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb
msigdb.H <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data//h.all.v7.2.entrez.gmt"))
msigdb.c2.cp <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c2.cp.v7.2.entrez.gmt"))
msigdb.c2 <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c2.all.v7.2.entrez.gmt"))
msigdb.c3 <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c3.all.v7.2.entrez.gmt"))
msigdb.c5 <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c5.all.v7.2.entrez.gmt"))

msigdb.names <- paste0("msigdb.", c("all", "H", "c2.cp", "c2", "c3", "c5"))

## get conversion list (gene name to entrez id)
convertion.list <- as.list(org.Hs.egALIAS2EG) %>% lapply(., `[[`, 1)
gene.to.genbank.id <- do.call(rbind, convertion.list) %>% as.data.frame() %>% `colnames<-`("genbank_id")
gene.to.genbank.id$Gene <- rownames(gene.to.genbank.id)
gene.to.genbank.id$genbank_id <- gene.to.genbank.id$genbank_id %>% as.character() %>% as.numeric()

## loop over all pathways
for (type in c("raw.score", "z.score")){ # "raw.score", "z.score"
    if (type == "z.score") input.data <- theta.zscore else input.data <- theta
    input.data <- input.data %>% as.data.frame() %>% mutate(Gene = rownames(.)) %>% melt(id.vars="Gene", value.name="weight", variable.name="topic") ## convert to long table
    fgsea.df.list <- vector("list", length(msigdb.names))
    for (msigdb.index in 1:length(msigdb.names)) {
        msigdb.pathway <- msigdb.names[msigdb.index]
        ## specify data set to use
        gsea.pathways=get(msigdb.pathway)
        ## initialize output storage variables
        ## out.gsea <- data.frame(ID=character(), Description=character(), setSize=double(), enrichmentScore=double(), NES=character(), pvalue=double(), p.adjust=double(), qvalues=double())
        fgseaRes <- vector("list", k)
        fgseaRes.df <- vector("list", k) # store fgseaRes that's converted to df format (remove leading Edge)
        for(t in c(1:k)){
            input.data.id <- input.data %>% subset(topic==t) %>% merge(., gene.to.genbank.id, by="Gene") %>% arrange(desc(weight))

            ## fgsea
            gene.list <- input.data.id$weight
            names(gene.list) <- input.data.id$genbank_id
            ## if (short.list) {
            ##   names(gene.list)[(number+1):length(gene.list)] <- 0
            ## }
            fgseaRes[[t]] <- fgsea(pathways = gsea.pathways, 
                                   stats    = gene.list,
                                   nperm    = 10000,
                                   ## eps      = 0.0,
                                   minSize  = 3,
                                   maxSize  = 800) %>% mutate(topic=t) %>% arrange(pval) 
            fgseaRes[[t]]$leadingEdge <- sapply(fgseaRes[[t]]$leadingEdge, paste, collapse=",") # convert list in leadingEdge to string
            fgseaRes.df[[t]] <- data.frame(pathway = fgseaRes[[t]]$pathway,
                                           pval = fgseaRes[[t]]$pval,
                                           padj = fgseaRes[[t]]$padj,
                                           ES = fgseaRes[[t]]$ES,
                                           NES = fgseaRes[[t]]$NES,
                                           nMoreExtreme = fgseaRes[[t]]$nMoreExtreme,
                                           numGenesInList = fgseaRes[[t]]$size,
                                           topic = fgseaRes[[t]]$topic)
        }
        fgsea.df.list[[msigdb.index]] <- do.call(rbind, fgseaRes.df) %>% mutate(database = msigdb.pathway, padj.over.topics = p.adjust(pval))
        ## save fgseaRes list
        save(fgseaRes, file=paste0(FGSEADIR,"/fgsea_",msigdb.pathway,"_",type,"_", SUBSCRIPT.SHORT, ".RData"))
    }
    fgsea.df <- do.call(rbind, fgsea.df.list)
    save(fgsea.df, file=paste0(FGSEADIR,"/fgsea_all_pathways_df_",type, "_", SUBSCRIPT.SHORT, ".RData"))
    write.table(fgsea.df, file=paste0(FGSEADIR, "/fgsea_", type, "_", SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep="\t")
    write.table(fgsea.df %>% subset(padj < 0.1), file=paste0(FGSEADIR, "/fgsea_", type, "_p.adj.0.1_", SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep="\t")
}



}




