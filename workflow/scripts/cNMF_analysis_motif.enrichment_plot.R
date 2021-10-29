## Helen Kang
## Topic Model Analysis 
## 210503

## .libPaths("/home/groups/engreitz/Software/R_3.6.1")


## packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
##               "CountClust", "Hmisc", "tidyr", "grid", "gtable", "gridExtra","ggrepel","ramify",
##               "GGally","RNOmni","usedist","ggpubr","gridExtra","GSEA",
##               "org.Hs.eg.db","limma","clusterProfiler","fgsea", "conflicted",
##               "cluster","textshape","readxl", "IsoplotR", "wesanderson", 
##               "ggdist", "patchwork", "gghalves", "Seurat", "writexl")
## library(textshape)
## library(readxl)
packages <- c("optparse","dplyr", "ggplot2", "reshape2", "ggrepel", "conflicted")
## library(Seurat)
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
  make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/", help="Output directory"),
  # make_option("--olddatadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/", help="Input 10x data directory"),
  make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/", help="Input 10x data directory"),
  # make_option("--topic.model.result.dir", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/210625_snakemake_output/top3000VariableGenes_acrossK/2kG.library/", help="Topic model results directory"),
  make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
  # make_option("--sep", type="logical", default=F, help="Whether to separate replicates or samples"),
  # make_option("--K.list", type="character", default="2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25", help="K values available for analysis"),
  make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
  # make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
  # make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
  # make_option("--ABCdir",type="character", default="/oak/stanford/groups/engreitz/Projects/ABC/200220_CAD/ABC_out/TeloHAEC_Ctrl/Neighborhoods/", help="Path to ABC enhancer directory"),
  make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
  # make_option("--raw.mtx.dir",type="character",default="stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/data/no_IL1B_filtered.normalized.ptb.by.gene.mtx.filtered.txt", help="input matrix to cNMF pipeline"),
  # make_option("--raw.mtx.RDS.dir",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.RDS", help="input matrix to cNMF pipeline"), # the first lane: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.expandedMultiTargetGuide.RDS"
  # make_option("--subsample.type", type="character", default="", help="Type of cells to keep. Currently only support ctrl"),
  # make_option("--barcode.names", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/barcodes.tsv", help="barcodes.tsv for all cells"),
  # make_option("--reference.table", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/210702_2kglib_adding_more_brief_ca0713.xlsx"),
  
  ## fisher motif enrichment
  ## make_option("--outputTable", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.k_14.df_0_2.txt", help="Output directory"),
  ## make_option("--outputTableBinary", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210607_snakemake_output/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.binary.k_14.df_0_2.txt", help="Output directory"),
  ## make_option("--outputEnrichment", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210607_snakemake_output/outputs/no_IL1B/topic.top.100.zscore.gene.motif.fisher.enrichment.k_14.df_0_2.txt", help="Output directory"),
  # make_option("--motif.promoter.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv", help="All promoter's motif matches"),
  # make_option("--motif.enhancer.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv", help="All enhancer's motif matches specific to {no,plus}_IL1B"),
  # make_option("--enhancer.fimo.threshold", type="character", default="1.0E-4", help="Enhancer fimo motif match threshold"),

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

mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))

SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
# STATIC.SAMPLE=c("Telo_no_IL1B_T200_1", "Telo_no_IL1B_T200_2", "Telo_plus_IL1B_T200_1", "Telo_plus_IL1B_T200_2", "no_IL1B", "plus_IL1B",  "pooled")
# DATADIR=opt$olddatadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
## TMDIR=opt$topic.model.result.dir
## SEP=opt$sep
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
## selected.gene <- c("EDN1", "NOS3", "TP53", "GOSR2", "CDKN1A")
# ABC genes
# gene.set <- c("INPP5B", "SF3A3", "SERPINH1", "NR2C1", "FGD6", "VEZT", "SMAD3", "AAGAB", "GOSR2", "ATP5G1", "ANGPTL4", "SRBD1", "PRKCE", "DAGLB") # ABC_0.015_CAD_pp.1_genes #200 gene library

# # cell cycle genes
# ## need to update these for 2kG library
# gene.list.three.groups <- read.delim(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/ptbd.genes_three.groups.txt"), header=T, stringsAsFactors=F)
# enhancer.set <- gene.list.three.groups$Gene[grep("E_at_", gene.list.three.groups$Gene)]
# CAD.focus.gene.set <- gene.list.three.groups %>% subset(Group=="CAD_focus") %>% pull(Gene) %>% append(enhancer.set)
# EC.pos.ctrl.gene.set <- gene.list.three.groups %>% subset(Group=="EC_pos._ctrls") %>% pull(Gene)

# cell.count.thr <- opt$cell.count.thr # greater than this number, filter to keep the guides with greater than this number of cells
# guide.count.thr <- opt$guide.count.thr # greater than this number, filter to keep the perturbations with greater than this number of guides

# guide.design = read.delim(file=paste0(DATADIR, "/200607_ECPerturbSeqMiniPool.design.txt"), header=T, stringsAsFactors = F)


# ## add GO pathway log2FC
# GO <- read.delim(file=paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/GO.Pathway.table.brief.txt"), header=T, check.names=FALSE)
# GO.list <- read.delim(file="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/GO.Pathway.list.brief.txt", header=T, check.names=F)
# colnames(GO)[1] <- "Gene"
# colnames(GO.list)[1] <- "Gene"
# ## load all sample, K, topic's top 100 genes (by TopFeatures() KL-score measure)
# ## allGeneKtopic100 <- read.delim(paste0(TMDIR, "no.plus.pooled.top100.topicStats.txt"), header=T)
# # load non-expressed control gene list
# non.expressed.genes <- read.delim(file="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/non.expressed.ctrl.genes.txt", header=F, stringsAsFactors=F) %>% unlist %>% as.character() %>% sort()

# # perturbation type list
# gene.set.type.df <- data.frame(Gene=guide.design %>% pull(guideSet) %>% unique(),
#                                type=rep("other", guide.design %>% pull(guideSet) %>% unique() %>% length())) 
# gene.set.type.df$Gene <- gene.set.type.df$Gene %>% as.character()
# gene.set.type.df$type <- gene.set.type.df$type %>% as.character()
# gene.set.type.df$type[which(gene.set.type.df$Gene %in% non.expressed.genes)] <- "non-expressed"
# gene.set.type.df$type[which(gene.set.type.df$Gene %in% CAD.focus.gene.set)] <- "CAD focus"
# gene.set.type.df$type[grepl("^safe|^negative", gene.set.type.df$Gene)] <- "negative-control"
# gene.set.type.df$Gene[which(gene.set.type.df$Gene == "negative_control")] <- "negative-control"
# gene.set.type.df$Gene[which(gene.set.type.df$Gene == "safe_targeting")] <- "safe-targeting"
# # gene.set.type.df$type[which(gene.set.type.df$Gene %in% gene.set)] <- "ABC"

# gene.set.type.df.200 <- gene.set.type.df

# # reference table
# ref.table <- read_xlsx(opt$reference.table, sheet="2000_gene_library_annotated") 
# gene.set.type.df <- ref.table %>% select(Symbol, `Class(es)`) %>% `colnames<-`(c("Gene", "type"))
# gene.set.type.df$type[grepl("EC_ctrls", gene.set.type.df$type)] <- "EC_ctrls"
# gene.set.type.df$type[grepl("NonExpressed", gene.set.type.df$type)] <- "non-expressed"
# gene.set.type.df$type[grepl("abc.015", gene.set.type.df$type)] <- "ABC"
# gene.set.type.df <- rbind(gene.set.type.df, c("negative-control", "negative-control"), c("safe-targeting", "safe-targeting"))
# non.expressed.genes <- gene.set.type.df %>% subset(type == "non-expressed") %>% pull(Gene)
# # ABC genes
# gene.set <- gene.set.type.df %>% subset(grepl("ABC", type)) %>% pull(Gene)

# ## add GWAS classification
# modified.ref.table <- ref.table %>% mutate(GWAS.classification="")
# CAD.index <- which(grepl("CAD_Loci",ref.table$`Class(es)`))
# EC_ctrls.index <- which(grepl("^EC_ctrls",ref.table$`Class(es)`))
# ABC_linked.index <- which(grepl("MIG_etc",ref.table$`Class(es)`))
# IBD.index <- which(grepl("Non-CAD_loci_IBD",ref.table$`Class(es)`))
# non.expressed.index <- which(grepl("NonExpressed",ref.table$`Class(es)`))
# poorly.annotated.9p21.index <- which(grepl("9p21",ref.table$`Class(es)`))
#                                         # length(CAD.index) + length(EC_ctrls.index) + length(ABC_linked.index) + length(IBD.index) + length(non.expressed.index) + length(poorly.annotated.9p21.index)
# modified.ref.table$GWAS.classification[ABC_linked.index] <- "ABC"
# modified.ref.table$GWAS.classification[IBD.index] <- "IBD"
# modified.ref.table$GWAS.classification[non.expressed.index] <- "NonExpressed"
# modified.ref.table$GWAS.classification[poorly.annotated.9p21.index] <- "9p21.poorly.annotated"
# modified.ref.table$GWAS.classification[EC_ctrls.index] <- "EC_ctrls"
# modified.ref.table$GWAS.classification[CAD.index] <- "CAD"

# modified.ref.table <- modified.ref.table %>% group_by(GWAS.classification) %>% mutate(gene.count.per.GWAS.category = n())
# ref.table <- modified.ref.table

# ## add TSS distance to SNP
# modified.ref.table <- ref.table %>% mutate(TSS.dist.to.SNP = abs(`TSS v. SNP loc`))
# not.in.SNP.index <- which(is.na(modified.ref.table$`TSS v. SNP loc`))
# modified.ref.table$TSS.dist.to.SNP[not.in.SNP.index] <- NA
# ref.table <- modified.ref.table %>% ungroup()

# ## add closest gene to top GWAS loci ranking
# modified.ref.table <- ref.table %>%
#     group_by(`Top SNP ID`) %>% # per SNP metrics
#     arrange(abs(`TSS v. SNP loc`)) %>%
#     mutate(TSS.v.SNP.ranking = 1:n(),
#            total.gene.in.this.loci = n()) %>% ungroup() %>%
#     group_by(`Top SNP ID`, GWAS.classification) %>% # per SNP per GWAS class (CAD, IBD, NonExpressed, ABC, 9p21.poorly.annotated) 
#     arrange(abs(`TSS v. SNP loc`)) %>%
#     mutate(TSS.v.SNP.ranking.in.GWAS.category = 1:n(),
#            total.gene.in.this.loci.in.GWAS.category = n()) %>% ungroup()
# not.in.SNP.index <- which(is.na(modified.ref.table$`TSS v. SNP loc`))
# modified.ref.table$TSS.v.SNP.ranking.in.GWAS.category[not.in.SNP.index] <- NA
# modified.ref.table$TSS.v.SNP.ranking[not.in.SNP.index] <- NA
# ref.table <- modified.ref.table

# ## add gene count per distance ranking per GWAS loci
# modified.ref.table <- ref.table
# modified.ref.table <- modified.ref.table %>%
#     group_by(TSS.v.SNP.ranking) %>% # per ranking, not considering which GWAS category the gene is from
#     mutate(total.TSS.v.SNP.ranking.count = n()) %>% ungroup() %>%
#     group_by(GWAS.classification, TSS.v.SNP.ranking.in.GWAS.category) %>% # per GWAS category and per ranking
#     mutate(total.TSS.v.SNP.ranking.count.per.GWAS.classification = n()) %>% ungroup()
# not.in.SNP.index <- which(is.na(modified.ref.table$TSS.v.SNP.ranking))
# modified.ref.table$total.TSS.v.SNP.ranking.count[not.in.SNP.index] <- NA
# modified.ref.table$total.TSS.v.SNP.ranking.count.per.GWAS.classification[not.in.SNP.index] <- NA
# ref.table <- modified.ref.table

# write.table(ref.table, file=paste0(opt$datadir, "/ref.table.txt"), row.names=F, quote=F, sep="\t")

# ## ref.table ranking count summary table
# ref.table.gene.to.SNP.dist.ranking.count.summary.allGWAS <- ref.table %>% select(TSS.v.SNP.ranking, total.TSS.v.SNP.ranking.count) %>% mutate(GWAS.classification="all") %>% unique()
# ref.table.gene.to.SNP.dist.ranking.count.summary.indGWAS <- ref.table %>% select(TSS.v.SNP.ranking.in.GWAS.category, total.TSS.v.SNP.ranking.count.per.GWAS.classification, GWAS.classification) %>% `colnames<-`(c("TSS.v.SNP.ranking", "total.TSS.v.SNP.ranking.count", "GWAS.classification")) %>% unique()
# ref.table.gene.to.SNP.dist.ranking.count.summary <- rbind(ref.table.gene.to.SNP.dist.ranking.count.summary.allGWAS, ref.table.gene.to.SNP.dist.ranking.count.summary.indGWAS)
# ref.table.summary.na.index <- which(is.na(ref.table.gene.to.SNP.dist.ranking.count.summary$TSS.v.SNP.ranking))
# ref.table.gene.to.SNP.dist.ranking.count.summary <- ref.table.gene.to.SNP.dist.ranking.count.summary[-ref.table.summary.na.index,]
# rm(ref.table.summary.na.index)


# # convert enhancer SNP rs number to enhancer target gene name # need 2kG library version
# enh.snp.to.gene <- read.delim(paste0(DATADIR, "/enhancer.SNP.to.gene.name.txt"), header=T, stringsAsFactors = F) %>% mutate(Enhancer_name=gsub("_","-", Enhancer_name))

# # gene corresponding pathway
# gene.def.pathways <- read_excel(paste0(DATADIR,"topic.gene.definition.pathways.xlsx"), sheet="Gene_Pathway")

# ## Gavin's new list
# gene.classes.ranked <- read.table(paste0(opt$datadir, "Gene_Classes_Ranked_for_CAD_n_EC.txt"), header=T, stringsAsFactors = F)
# summaries <- read.delim(paste0(opt$datadir, "Gene_Summaries_n_Classes.txt"), sep="\t", header=T, stringsAsFactors = F)
# gene.summaries <- read_xlsx(paste0(opt$datadir, "Gene_Summaries.xlsx"), sheet="uniprot_summaries")


# ## Perturbation name and 10X gene name conversion table
# ptb.10X.name.conversion <- read_xlsx(paste0(opt$datadir, "Perturbation 10X names.xlsx"))

# ## EdgeR log2fcs and p-values
# log2fc.edgeR <- read.table(paste0(opt$datadir, "/EdgeR/ALL_log2fcs_dup4_s4n3.99x.txt"), header=T, stringsAsFactors=F)
# p.value.edgeR <- read.table(paste0(opt$datadir, "/EdgeR/ALL_Pvalues_dup4_s4n3.99x.txt"), header=T, stringsAsFactors=F)

# print("loaded all prerequisite data")




######################################################################
## Process topic model results
cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
} else {
    warning(paste0(cNMF.result.file, " does not exist"))
}

# ## load ann.omega
# file.name <- ifelse(SEP,
#                     paste0(OUTDIRSAMPLE,"/cNMFAnalysis.",SUBSCRIPT,".sep.RData"),
#                     paste0(OUTDIRSAMPLE,"/cNMFAnalysis.",SUBSCRIPT,".RData"))
# print(file.name) 
# if(file.exists((file.name))) { 
#     print(paste0("loading ",file.name))
#     load(file.name) 
# }



# ## load statistical test data
# toSave.features <- read.delim(paste0(OUTDIRSAMPLE, "/topic.KL.score_K", k, ".dt_", DENSITY.THRESHOLD, ".txt"), header=T, stringsAsFactors=F)
# all.test <- read.delim(file=paste0(OUTDIRSAMPLE, "/all.test.", SUBSCRIPT, ".txt"), header=T, stringsAsFactors=F) ##here
# realPvals.df <- read.delim(file=paste0(OUTDIRSAMPLE, "/all.expressed.genes.pval.fdr.",SUBSCRIPT,".txt"), header=T, stringsAsFactors=F)
# all.test.guide.w <- all.test %>% subset(test.type==opt$test.type)
# realPvals.df.guide.w <- realPvals.df %>% subset(test.type==opt$test.type)
# fdr.thr <- 0.1
# topFeatures.raw.weight <- theta.zscore %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene", variable.name="topic", value.name="scores") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:50)


## load Seurat Object with UMAP


# ## load full matrix
# file.name=paste0(opt$raw.mtx.RDS.dir, "_modified.multiTargetGuide.cell.names.RDS")
# fc.file.name=paste0(opt$raw.mtx.RDS.dir, "_FC_modified.multiTargetGuide.cell.names.RDS")
# log2fc.file.name=paste0(opt$raw.mtx.RDS.dir, "_log2FC_modified.multiTargetGuide.cell.names.RDS")
# if(file.exists(file.name) & file.exists(fc.file.name) & file.exists(log2fc.file.name)) {
#     print(paste0("Loading ", file.name))
#     X.full = readRDS(file.name)
#     print(paste0("Loading ", fc.file.name))
#     fc.X.full = readRDS(fc.file.name)
#     print(paste0("Loading ", log2fc.file.name))
#     log2fc.X.full = readRDS(log2fc.file.name)##here210813
# } else {
#     warning(paste0("full scRNA-seq matrix does not exist"))
# }


# ## modify X.full
# ## old 210819
# tokeep.index <- which(rownames(X.full) %in% ann.omega.filtered$long.CBC) ## need to expand X.full to match ann.omega.filtered due to guides that target two promoters
# fc.X.full <- fc.X.full[tokeep.index,] 
# X.full <- X.full[tokeep.index,] 
# colnames(X.full) <- colnames(fc.X.full) <- colnames(fc.X.full)  %>% strsplit(., split=":") %>% sapply("[[",1)

# ## ## adjust colnames, remove ENSG number
# ann.X.full.filtered <- X.full
# ## add back ENSG names?
# ## tmp <- colnames(ann.X.full.filtered) %>% strsplit(., split=":") %>% sapply("[[",1)
# ## tmpp <- data.frame(table(tmp)) %>% subset(Freq > 1)  # keep row names that have duplicated gene names but different ENSG names
# ## tmp.copy <- tmp
# ## tmp.copy[grepl(paste0(tmpp$tmp,collapse="|"),tmp)] <- colnames(ann.X.full.filtered)[grepl(paste0(tmpp$tmp,collapse="|"), colnames(ann.X.full.filtered))]
# ## colnames(ann.X.full.filtered) <- tmp.copy # the above section takes a while
# ## end of old 210819

# ## get ctrl log2 transformed expression
# tokeep.index <- which(grepl("control|targeting",rownames(X.full)))
# ## ctrl.X <- ann.X.full.filtered %>% subset(grepl("control|targeting",Gene))
# ctrl.X <- X.full[tokeep.index,]
# fc.ctrl.X <- fc.X.full[tokeep.index,]
# ## get ctrl topic weights
# ctrl.ann.omega <- ann.omega.filtered %>% subset(grepl("control|targeting",Gene)) %>% `rownames<-`(.$long.CBC)
# X.gene.names <- rownames(X.full) %>% strsplit(., split=":") %>% sapply("[[",1) %>% gsub("_multiTarget|-TSS2","",.)

## load motif enrichment results
file.name <- paste0(OUTDIRSAMPLE,"/cNMFAnalysis.factorMotifEnrichment.",SUBSCRIPT.SHORT,".RData")
print(file.name)
if(file.exists((file.name))) { 
    load(file.name)
    print(paste0("loading ", file.name))
}
motif.enrichment.variables <- c("all.enhancer.fisher.df", "all.promoter.fisher.df", 
                                "promoter.wide", "enhancer.wide", "promoter.wide.binary", "enhancer.wide.binary",
                                "enhancer.wide.10en6", "enhancer.wide.binary.10en6", "all.enhancer.fisher.df.10en6",
                                "promoter.wide.10en6", "promoter.wide.binary.10en6", "all.promoter.fisher.df.10en6",
                                "all.promoter.ttest.df", "all.promoter.ttest.df.10en6", "all.enhancer.ttest.df", "all.enhancer.ttest.df.10en6")
motif.enrichment.variables.missing <- (!(motif.enrichment.variables %in% ls())) %>% as.numeric %>% sum 
if ( motif.enrichment.variables.missing > 0 ) {
    warning(paste0(motif.enrichment.variables[!(motif.enrichment.variables %in% ls())], " not available"))
}

# ## load count.by.GWAS
# count.by.GWAS <- read.delim(file=paste0(OUTDIRSAMPLE,"/count.by.GWAS.classes_p.adj.",p.value.thr %>% as.character,"_",SUBSCRIPT,".txt"), header=T, stringsAsFactors = F)
# count.by.GWAS.withTopic <- read.delim(file=paste0(OUTDIRSAMPLE,"/count.by.GWAS.classes.withTopic_p.adj.",p.value.thr %>% as.character,"_",SUBSCRIPT,".txt"), header=T, stringsAsFactors=F)


# ## load UMAP data
# file.name <- paste0(OUTDIRSAMPLE, "gene.score.SeuratObject.RDS")
# if(file.exists(file.name)) {
#     s.gene.score <- readRDS(file.name)
# } else {
#     warning(paste0(file.name, " does not exist"))
# } 

# ## full UMAP
# file.name <- paste0(opt$datadir,"/", SAMPLE, "_SeuratObjectUMAP.rds") ## todo: change to calcUMAP output
# ## subset.file.name <- paste0(opt$datadir,"/", SAMPLE, "_subset_SeuratObjectUMAP.rds")
# options(future.globals.maxSize=1000*1024^2)
# s <- readRDS(file.name)
# ## s@meta.data <- s@meta.data[,-which(grepl("topic",colnames(s@meta.data)))]
# s <- AddMetaData(s, metadata = omega, col.name = paste0("K",k,"_",colnames(omega)))


## End of data loading



##########################################################################
## Plots


## volcano plots
volcano.plot <- function(toplot, ep.type, ranking.type, label.type="") {
    if( label.type == "pos") {
        label <- toplot %>% subset(-log10(p.adjust) > 1 & enrichment.log2fc > 0) %>% mutate(motif.toshow = gsub("HUMAN.H11MO.", "", motif))
    } else {
        label <- toplot %>% subset(-log10(p.adjust) > 1) %>% mutate(motif.toshow = gsub("HUMAN.H11MO.", "", motif))
    }
    t <- gsub("topic_", "", toplot$topic[1])
    p <- toplot %>% ggplot(aes(x=enrichment.log2fc, y=-log10(p.adjust))) + geom_point(size=0.5) + mytheme +
        ggtitle(paste0(SAMPLE[1], " Topic ", t, " Top 100 ", ranking.type," ", ifelse(ep.type=="promoter", "Promoter", "Enhancer"), " Motif Enrichment")) + xlab("Motif Enrichment (log2FC)") + ylab("-log10(adjusted p-value)") +
        geom_text_repel(data=label, box.padding = 0.5,
                        aes(label=motif.toshow), size=5,
                        color="black") + theme(text=element_text(size=16), axis.title=element_text(size=16), axis.text=element_text(size=16), plot.title=element_text(size=14))
    print(p)
    p <- toplot %>% ggplot(aes(x=enrichment.log2fc, y=-log10(p.value))) + geom_point(size=0.5) + mytheme +
        ggtitle(paste0(SAMPLE[1], " Topic ", t, " Top 100 ", ranking.type," ", ifelse(ep.type=="promoter", "Promoter", "Enhancer"), " Motif Enrichment")) + xlab("Motif Enrichment (log2FC)") + ylab("-log10(p-value)") +
        geom_text_repel(data=label, box.padding = 0.5,
                        aes(label=motif.toshow), size=5,
                        color="black") + theme(text=element_text(size=16), axis.title=element_text(size=16), axis.text=element_text(size=16), plot.title=element_text(size=14))
    return(p)
}

## function for all volcano plots
all.volcano.plots <- function(all.fisher.df, ep.type, ranking.type, label.type="") {
    for ( t in 1:k ){
        toplot <- all.fisher.df %>% subset(topic==paste0("topic_",t))
        volcano.plot(toplot, ep.type, ranking.type, label.type) %>% print()
    }
}


##########################################################################
## motif enrichment plot
ep.names <- c("enhancer", "promoter")
for (ep.type in ep.names) {
    pdf(file=paste0(FIGDIRTOP,"zscore.",ep.type,".motif.enrichment.pdf"))
    all.volcano.plots(get(paste0("all.",ep.type,".fisher.df")), ep.type, ranking.type="z-score")
    dev.off()
    pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment_motif.thr.10e-6.pdf"))
    all.volcano.plots(get(paste0("all.",ep.type,".fisher.df.10en6")), ep.type, ranking.type="z-score")
    dev.off()
    pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest.pdf"))
    all.volcano.plots(get(paste0("all.",ep.type,".ttest.df")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score")
    dev.off() 
    pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest_motif.thr.10e-6.pdf"))
    all.volcano.plots(get(paste0("all.",ep.type,".ttest.df.10en6")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score")
    dev.off() 
    pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest.labelPos.pdf"), width=6, height=6)
    all.volcano.plots(get(paste0("all.",ep.type,".ttest.df")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score", label.type="pos")
    dev.off() 
    pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest_motif.thr.10e-6.labelPos.pdf"))
    all.volcano.plots(get(paste0("all.",ep.type,".ttest.df.10en6")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score", label.type="pos")
    dev.off() 
}

