
## Helen Kang
## Topic Model Analysis Only (no plot output)
## 210503

## .libPaths("/home/groups/engreitz/Software/R_3.6.1")
library(conflicted)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")


packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
              "tidyr", "grid", "gtable", "gridExtra","ggrepel","ramify",
              "ggpubr","gridExtra","RNOmni",
              "org.Hs.eg.db","conflicted",
              "textshape","readxl", 
              "ggdist", "gghalves", "writexl") #              "GGally","RNOmni","usedist","GSEA","clusterProfiler","IsoplotR","wesanderson",
xfun::pkg_attach(packages)
conflict_prefer("combine", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")



## source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")
## source(file.path(snakemake@scriptdir, "topicModelAnalysis.functions.R"))

option.list <- list(
    make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/figures/all_genes/", help="Figure directory"),
    make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/", help="Output directory"),
                                        # make_option("--olddatadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/", help="Input 10x data directory"),
    make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/", help="Input 10x data directory"),
                                        # make_option("--topic.model.result.dir", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/211116_snakemake_dup4_cells/all_genes_acrossK/2kG.library/", help="Topic model results directory"),
    make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
    make_option("--sep", type="logical", default=F, help="Whether to separate replicates or samples"),
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
    make_option("--recompute", type="logical", default=F, help="T for recomputing statistical tests and F for not recompute"),
    
    ## script dir
    make_option("--scriptdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/cNMF_pipeline/Perturb-seq/workflow/scripts/", help="location for this script and functions script")
)
opt <- parse_args(OptionParser(option_list=option.list))


source(paste0(opt$scriptdir, "/topicModelAnalysis.functions.R"))

## ## dup4_2n1.99x
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/figures/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/"
## opt$K.val <- 5
## opt$sampleName <- "Perturb_2kG_dup4"

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

## ## debug scRNAseq_2kG_11AMDox_1
## opt$sampleName <- "scRNAseq_2kG_11AMDox_1"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/figures/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes/"
## opt$K.val <- 14
## opt$topic.model.result.dir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/analysis/all_genes_acrossK/scRNAseq_2kG_11AMDox_1"


mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))

SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
DATADIR=opt$olddatadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
                                        # TMDIR=opt$topic.model.result.dir
SEP=opt$sep
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
SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr

## ## directories for factor motif enrichment
## FILENAME=opt$filename


## ## modify motif.enhancer.background input directory ##HERE: perhaps do a for loop for all the desired thresholds (use strsplit on enhancer.fimo.threshold)
## opt$motif.enhancer.background <- paste0(opt$motif.enhancer.background, opt$enhancer.fimo.threshold, "/fimo.formatted.tsv")


                                        # create dir if not already
check.dir <- c(OUTDIR, FIGDIR, OUTDIRSAMPLE, FIGDIRSAMPLE, FGSEADIR, FGSEAFIG)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

## palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)
                                        # selected.gene <- c("EDN1", "NOS3", "TP53", "GOSR2", "CDKN1A")
                                        # # ABC genes
                                        # gene.set <- c("INPP5B", "SF3A3", "SERPINH1", "NR2C1", "FGD6", "VEZT", "SMAD3", "AAGAB", "GOSR2", "ATP5G1", "ANGPTL4", "SRBD1", "PRKCE", "DAGLB") # ABC_0.015_CAD_pp.1_genes #200 gene library

                                        # # cell cycle genes
                                        # ## need to update these for 2kG library
                                        # gene.list.three.groups <- read.delim(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/ptbd.genes_three.groups.txt"), header=T, stringsAsFactors=F)
                                        # enhancer.set <- gene.list.three.groups$Gene[grep("E_at_", gene.list.three.groups$Gene)]
                                        # CAD.focus.gene.set <- gene.list.three.groups %>% subset(Group=="CAD_focus") %>% pull(Gene) %>% append(enhancer.set)
                                        # EC.pos.ctrl.gene.set <- gene.list.three.groups %>% subset(Group=="EC_pos._ctrls") %>% pull(Gene)

cell.count.thr <- opt$cell.count.thr # greater than this number, filter to keep the guides with greater than this number of cells
guide.count.thr <- opt$guide.count.thr # greater than this number, filter to keep the perturbations with greater than this number of guides

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

                                        # reference table
                                        ref.table <- read_xlsx(opt$reference.table, sheet="2000_gene_library_annotated") 
                                        gene.set.type.df <- ref.table %>% select(Symbol, `Class(es)`) %>% `colnames<-`(c("Gene", "type"))
                                        gene.set.type.df$type[grepl("EC_ctrls", gene.set.type.df$type)] <- "EC_ctrls"
                                        gene.set.type.df$type[grepl("NonExpressed", gene.set.type.df$type)] <- "non-expressed"
                                        gene.set.type.df$type[grepl("abc.015", gene.set.type.df$type)] <- "ABC"
                                        gene.set.type.df <- rbind(gene.set.type.df, c("negative-control", "negative-control"), c("safe-targeting", "safe-targeting"))
                                        non.expressed.genes <- gene.set.type.df %>% subset(type == "non-expressed") %>% pull(Gene)
                                        # ABC genes
                                        gene.set <- gene.set.type.df %>% subset(grepl("ABC", type)) %>% pull(Gene)

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


                                        # print("loaded all prerequisite data")



                                        # ## for the guides that target multiple genes, we will split the gene annotation and duplicate the cell entry, so that each gene will get a cell read out.
                                        # adjust.multiTargetGuide.rownames <- function(omega) {
                                        #     ## duplicate cells with guide that targets multiple genes
                                        #     cells.with.multiTargetGuide <- rownames(omega)[grepl("and",rownames(omega))]
                                        #     if(length(cells.with.multiTargetGuide) > 0) {
                                        #         ## split index by number of guide targets
                                        #         cells.with.multiTargetGuide.index <- which(grepl("and",rownames(omega)))
                                        #         cells.with.singleTargetGuide.index <- which(!grepl("and",rownames(omega)))
                                        #         ## get multi target gene names
                                        #         multiTarget.names <- cells.with.multiTargetGuide %>% strsplit(., split=":") %>% sapply("[[",1) ## full names: GeneA-and-GeneB
                                        #         multiTarget.Guide.CBC <- cells.with.multiTargetGuide %>% strsplit(., split=":") %>% sapply( function(x) paste0(x[[2]], ":", x[[3]]) )
                                        #         multiTarget.names.1 <- multiTarget.names %>% strsplit(., split="-and-") %>% sapply ("[[",1)
                                        #         multiTarget.names.2 <- multiTarget.names %>% strsplit(., split="-and-") %>% sapply ("[[",2)
                                        #         multiTarget.names.all <- multiTarget.names.1 %>% append(multiTarget.names.2) %>% unique() ## get all the genes/enhancers that have guides targeting other gene/enhancer at the same time
                                        #         cells.with.multiTarget.gene.names.index <- which(grepl(paste0(multiTarget.names.all,collapse="|"), rownames(omega)))

                                        #         multiTarget.long.CBC.1 <- sapply(1:length(multiTarget.names), function(i) {
                                        #             paste0(multiTarget.names.1[i], "_multiTarget:", multiTarget.Guide.CBC[i])
                                        #         })
                                        #         multiTarget.long.CBC.2 <- sapply(1:length(multiTarget.names), function(i) {
                                        #             paste0(multiTarget.names.2[i], "_multiTarget:", multiTarget.Guide.CBC[i])
                                        #         })
                                        #         ## change original df's rownames
                                        #         multiTargetGuide.mtx <- omega[cells.with.multiTargetGuide.index,]
                                        #         multiTargetGuide.mtx.1 <- multiTargetGuide.mtx %>% `rownames<-`(multiTarget.long.CBC.1)
                                        #         multiTargetGuide.mtx.2 <- multiTargetGuide.mtx %>% `rownames<-`(multiTarget.long.CBC.2)
                                        #         ## pull cells with guides that has a single target, but the perturbed gene has multiTarget guide
                                        #         expanded.gene.name.df <- do.call(rbind, lapply(1:length(multiTarget.names.all), function(i) {
                                        #             gene.name.here <- multiTarget.names.all[i]
                                        #             toPaste.gene.name.here <- multiTarget.names[which(grepl(gene.name.here, multiTarget.names))] %>% gsub("-TSS2","",.) %>% unique()
                                        #             out <- do.call(rbind, lapply(1:length(toPaste.gene.name.here), function(j) {
                                        #                 singleTarget.cell.index.here <- which(grepl(gene.name.here,rownames(omega)) & !grepl("and",rownames(omega)))
                                        #                 singleTarget.cell.df <- omega[singleTarget.cell.index.here,]
                                        #                 rownames(singleTarget.cell.df) <- gsub(gene.name.here,toPaste.gene.name.here[j],rownames(singleTarget.cell.df))
                                        #                 return(singleTarget.cell.df)
                                        #             }))
                                        #             return(out)
                                        #         })) # takes two minutes
                                        #         omega <- rbind(omega[cells.with.singleTargetGuide.index,], multiTargetGuide.mtx.1, multiTargetGuide.mtx.2, expanded.gene.name.df)
                                        #     }
                                        #     return(omega)
                                        # }    


######################################################################
## Process topic model results
## for ( n in 1:length(SAMPLE) ) {

## if (SEP) {
##     guideCounts <- loadGuides(n, sep=T) %>% mutate(Gene=Gene.marked)
##     tmp.labels <- guideCounts$Gene %>% unique() %>% strsplit("-") %>% sapply("[[",2) %>% unique()
##     tmp.labels <- tmp.labels[!(tmp.labels %in% c("control","targeting"))]
##     rep1.label <- paste0("-",tmp.labels[1])
##     rep2.label <- paste0("-",tmp.labels[2])
## } else guideCounts <- loadGuides(n) %>% mutate(Gene=Gene.marked)


cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
print(cNMF.result.file)
if(file.exists(cNMF.result.file)) {
    print("loading cNMF result file")
    load(cNMF.result.file)
} else {
    theta.path <- paste0(TMDIR, "/", SAMPLE, ".gene_spectra_tpm.k_", k, ".dt_", DENSITY.THRESHOLD,".txt")
    theta.zscore.path <- paste0(TMDIR, "/", SAMPLE, ".gene_spectra_score.k_", k, ".dt_", DENSITY.THRESHOLD,".txt")
    print(theta.path)
    theta.raw <- read.delim(theta.path, header=T, stringsAsFactors=F, check.names=F, row.names=1)
    ## theta.raw <- read.delim(theta.path, header=T, stringsAsFactors=F, check.names=F) %>% select(-``)
    print("finished reading raw weights for topics")
    tmp.theta <- theta.raw
    tmp.theta[tmp.theta==0] <- min(tmp.theta[tmp.theta > 0])/100
    theta <- tmp.theta %>% apply(1, function(x) x/sum(x)) %>% `colnames<-`(c(1:k))
    theta.raw <- theta.raw %>% t() %>% as.data.frame() %>% `colnames<-`(c(1:k))
    print("loading topic z-score (specificity score)")
    theta.zscore <- read.delim(theta.zscore.path, header=T, stringsAsFactors=F, check.names=F, row.names=1) %>% t() %>% `colnames<-`(c(1:k))
    tmp <- rownames(theta) %>% strsplit(., split=":") %>% sapply("[[",1)
    tmpp <- data.frame(table(tmp)) %>% subset(Freq > 1)  # keep row names that have duplicated gene names but different ENSG names
    tmp.copy <- tmp
    tmp.copy[grepl(paste0(tmpp$tmp,collapse="|"),tmp)] <- rownames(theta)[grepl(paste0(tmpp$tmp,collapse="|"),rownames(theta))]
    rownames(theta) <- rownames(theta.raw) <- rownames(theta.zscore) <- tmp.copy
    ## truncate.theta.names <- function(theta) {
    ##     theta.gene.names <- rownames(theta) %>% strsplit(., split=":") %>% sapply("[[",1) # remove ENSG names
    ##     rownames(theta) <- theta.gene.names
    ##     return(theta)
    ## }
    ## theta.raw <- truncate.theta.names(theta.raw)
    ## theta.zscore <- truncate.theta.names(theta.zscore)
    omega.path <- paste0(TMDIR, "/", SAMPLE, ".usages.k_", k, ".dt_", DENSITY.THRESHOLD, ".consensus.txt")
    print(omega.path)
    omega.original <- omega <- read.delim(omega.path, header=T, stringsAsFactors=F, check.names=F, row.names = 1)  %>% apply(1, function(x) x/sum(x)) %>% t()
    colnames(omega) <- paste0("topic_",colnames(omega))
                                        # barcode.names <- read.table(opt$barcode.names, header=F, stringsAsFactors=F) %>% `colnames<-`("long.CBC")
                                        # rownames(omega) <- rownames(omega.original) <- barcode.names %>% pull(long.CBC) %>% gsub("CSNK2B-and-CSNK2B", "CSNK2B",.)
                                        # omega <- adjust.multiTargetGuide.rownames(omega)
                                        # barcode.names <- data.frame(long.CBC=rownames(omega)) %>% separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC"), sep=":", remove=F) %>% mutate(Gene = gsub("_multiTarget|-TSS", "", Gene.full.name))

    print("save the data")
    ensembl.theta.zscore.names <- mapIds(org.Hs.eg.db, keys = rownames(theta.zscore), keytype = "SYMBOL", column="ENSEMBL")
    ensembl.theta.zscore.names[ensembl.theta.zscore.names %>% is.na] <- rownames(theta.zscore)[ensembl.theta.zscore.names %>% is.na]
    theta.zscore.ensembl <- theta.zscore
    colnames(theta.zscore.ensembl) <- paste0("zscore_K", k, "_topic", colnames(theta.zscore.ensembl))
    theta.zscore.ensembl <- theta.zscore.ensembl %>% as.data.frame %>% mutate(ENSGID=ensembl.theta.zscore.names,.before=paste0("zscore_K",k,"_topic1"))

    ensembl.theta.raw.names <- mapIds(org.Hs.eg.db, keys = rownames(theta.raw), keytype = "SYMBOL", column="ENSEMBL")
    ensembl.theta.raw.names[ensembl.theta.raw.names %>% is.na] <- rownames(theta.raw)[ensembl.theta.raw.names %>% is.na]
    theta.raw.ensembl <- theta.raw
    colnames(theta.raw.ensembl) <- paste0("raw_K", k, "_topic", colnames(theta.raw.ensembl))
    theta.raw.ensembl <- theta.raw.ensembl %>% as.data.frame %>% mutate(ENSGID=ensembl.theta.raw.names,.before=paste0("raw_K",k,"_topic1"))

    ## normalize to zero mean + unit variance
    theta.raw.ensembl.scaled <- theta.raw.ensembl %>% select(-ENSGID) %>% apply(2, scale)  %>% as.data.frame %>% mutate(ENSGID=ensembl.theta.zscore.names,.before=paste0("raw_K",k,"_topic1"))
    theta.zscore.ensembl.scaled <- theta.zscore.ensembl %>% select(-ENSGID) %>% apply(2, scale) %>% as.data.frame %>% mutate(ENSGID=ensembl.theta.zscore.names,.before=paste0("zscore_K",k,"_topic1"))
    
    
    save(theta, theta.raw, theta.zscore, omega, theta.path, omega.path, # barcode.names,
         file=cNMF.result.file)

    write.table(theta.zscore, file=paste0(OUTDIRSAMPLE, "/topic.zscore_",SUBSCRIPT.SHORT, ".txt"), row.names=T, quote=F, sep="\t")
    write.table(theta.raw, file=paste0(OUTDIRSAMPLE, "/topic.raw.score_",SUBSCRIPT.SHORT, ".txt"), row.names=T, quote=F, sep="\t")
    write.table(theta.zscore.ensembl, file=paste0(OUTDIRSAMPLE, "/topic.zscore.ensembl_",SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep="\t")
    write.table(theta.raw.ensembl, file=paste0(OUTDIRSAMPLE, "/topic.raw.ensembl_",SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep="\t")
    write.table(theta.zscore.ensembl.scaled, file=paste0(OUTDIRSAMPLE, "/topic.zscore.ensembl.scaled_", SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep = "\t")
    write.table(theta.raw.ensembl.scaled, file=paste0(OUTDIRSAMPLE, "/topic.raw.ensembl.scaled_", SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep = "\t")
    
}  


                                        # # modify GO.list if "pooled"
                                        # if (SEP) {
                                        #   tmp.plus <- GO.list %>% mutate(Gene = paste0(GO.list$Gene,rep1.label), Pathway = paste0(GO.list$Pathway,rep1.label))
                                        #   tmp.no <- GO.list %>% mutate(Gene = paste0(GO.list$Gene,rep2.label), Pathway = paste0(GO.list$Pathway,rep2.label))
                                        #   GO.list <- rbind(tmp.no, tmp.plus)
                                        # }

                                        # for ( k in K.list ) {
file.name <- ifelse(SEP,
                    paste0(OUTDIRSAMPLE,"/cNMFAnalysis.",SUBSCRIPT,".sep.RData"),
                    paste0(OUTDIRSAMPLE,"/cNMFAnalysis.",SUBSCRIPT,".RData"))
print(file.name) 
if(file.exists((file.name))) { 
    print(paste0("loading ",file.name))
    load(file.name) 
}
if( (opt$subsample.type != "ctrl") & ( !("fc.ann.omega" %in% ls()) | opt$recompute) ) { # if the variable fc.ann.omega is present, then it means the below calculations are done previously

    barcode.names <- data.frame(long.CBC = rownames(omega)) %>%
        separate(col="long.CBC", into=c("Gene.full.name", "Guide", "CBC-sample"), sep=":", remove=F) %>%
        separate(col="CBC-sample", into=c("CBC", "sample"), sep="-", remove=F) %>%
        mutate(Gene = gsub("-TSS$", "", Gene.full.name))
    
                                        # get topic colnames
    topic.names <- colnames(omega)[which(grepl("topic",colnames(omega)))]
    
    ann.omega <- cbind(omega, barcode.names)  # annotated.omega
    
                                        # remove guides with less than 5 cells and remove genes with only one resulting guide
    ann.omega.filtered <- ann.omega %>%
        group_by(Gene, Guide) %>%
        mutate(nCount=n()) %>%
        subset(nCount > cell.count.thr) %>%
        select(-nCount) %>%
        ungroup() %>%
        group_by(Gene) %>%
        mutate(nCount = length(unique(Guide))) %>%
        subset(nCount > guide.count.thr) %>%
        select(-nCount) %>%
        as.data.frame()
    rownames(ann.omega.filtered) <- ann.omega.filtered$long.CBC
    
    ## neg.omega <- ann.omega.filtered %>% subset(grepl("^negative|^safe", Gene)) 
    
    avg.gene <- ann.omega.filtered %>%
        group_by(Gene) %>%
        summarise_at(colnames(omega), mean) %>%
        as.data.frame()
    rownames(avg.gene) <- avg.gene$Gene
    avg.neg <- avg.gene %>%
        subset(grepl("^negative|^safe", Gene)) %>%
        mutate(Gene = "combined-control") %>%
        group_by(Gene) %>%
        summarise_at(colnames(omega), mean) %>%
        as.data.frame()
    
    ## wilcoxon section ##hiGavin!
    
##################################################
    ## local function for computing p.values
    get.p.value <- function(df.ann.filtered, INT=F) {
        
        if(SEP) {
            num.reps=2
            labels <- c(rep1.label, rep2.label) # expand to accomodate multiple samples
            perguide.gene.test.list <- vector("list", ncol(omega)*num.reps) # iterate through topics + reps
            for (label.index in 1:num.reps){
                label.here <- labels[label.index]
                ctrl.this.rep <- df.ann.filtered %>% subset(grepl("^negative|^safe", Gene) & grepl(label.here, Gene, fixed=TRUE)) 
                df.ann.filtered.this.rep <- df.ann.filtered %>% subset(grepl(label.here, Gene, fixed=TRUE)) 
                for ( t in 1:length(colnames(omega)) ) {
                    col <- colnames(omega)[t] # current topic
                    ctrl <- ctrl.this.rep %>% select(all_of(col)) %>% as.matrix() %>% as.numeric()
                    if(INT) {
                        result <- df.ann.filtered.this.rep %>% group_by(Gene) %>% summarise(ttest.p = t.test(x = get(col) %>% as.numeric(), y = ctrl)$p.value)
                    } else {
                        result <- df.ann.filtered.this.rep %>% group_by(Gene) %>% 
                            summarise(wilcox.p = wilcox.test(x = get(col) %>% as.numeric(), y = ctrl)$p.value,
                                      ttest.p = t.test(x = get(col) %>% as.numeric(), y = ctrl)$p.value)
                    }
                    result$Topic <- rep(col, nrow(result))
                    perguide.gene.test.list[[label.index * t]] <- result
                }
            }
        } else {
            perguide.gene.test.list <- vector("list", ncol(omega)) # iterate through topics
            for ( t in 1:length(colnames(omega)) ) {
                col <- colnames(omega)[t] # current topic
                
                ctrl <- df.ann.filtered %>% select(all_of(col), Gene) %>% subset(grepl("control|targeting",Gene)) %>% select(-Gene) %>% as.matrix() %>% as.numeric()
                if(INT) {
                    result <- df.ann.filtered.this.rep %>% group_by(Gene) %>% summarise(ttest.p = t.test(x = get(col) %>% as.numeric(), y = ctrl)$p.value)
                } else {
                    result <- df.ann.filtered %>% group_by(Gene) %>% summarise(wilcox.p = wilcox.test(x = get(col) %>% as.numeric(), y = ctrl)$p.value,
                                                                               ttest.p = t.test(x = get(col) %>% as.numeric(), y = ctrl)$p.value)
                }
                result$Topic <- rep(col, nrow(result))
                perguide.gene.test.list[[t]] <- result
            }
        }
        perguide.gene.test <- do.call(rbind, perguide.gene.test.list)
        
        if(INT) {
            perguide.gene.test <- perguide.gene.test %>% mutate(adjusted.ttest.p = p.adjust(ttest.p, method="fdr"))
        } else { 
            perguide.gene.test <- perguide.gene.test %>% mutate(adjusted.wilcox.p = p.adjust(wilcox.p, method="fdr"), adjusted.ttest.p = p.adjust(ttest.p, method="fdr"))
        }
        return(perguide.gene.test)
    }
##################################################
    
    
    if(!SEP){
        ## avg cell and avg ctrl before comparing (then log2() + heatmap)
        avg.gene <- ann.omega.filtered %>% group_by(Gene) %>% summarise_at(colnames(omega), mean)
        avg.neg <- ann.omega.filtered %>% subset(grepl("^negative|^safe", Gene)) %>% mutate(Gene = "combined-control") %>% group_by(Gene) %>% summarise_at(colnames(omega), mean)
        mat <- data.frame(avg.gene[,-1], row.names = avg.gene$Gene) %>% as.matrix()
        vec <- avg.neg %>% select(colnames(omega)) %>% as.matrix()
        sweep.wrapper <- function(mat, vec, diff=3) {
            output <- sweep(mat, 2, vec, `/`) %>% log2()
            output[!is.finite(output)] <- floor(min(output[is.finite(output)]) - diff)
            output <- output %>% as.data.frame()
            return(output)
        }
        gene.score <- sweep.wrapper(mat,vec)
        write.table(gene.score, file=paste0(OUTDIRSAMPLE, "/avg.gene_x_topic.mtx_", SUBSCRIPT, ".txt"), row.names=T, col.names=T, quote=F, sep="\t")
        ##   gene.score <- sweep(mat, 2, vec, `/`) %>% log2()
        ##   gene.score[!is.finite(gene.score)] <- floor(min(gene.score[is.finite(gene.score)]) - 3)
        ## gene.score <- gene.score %>% as.data.frame()
        
        ## avg guide and avg ctrl before comparing (then log2() + heatmap)
        avg.guide <- ann.omega.filtered %>% mutate(Guide = paste0(Guide,"::",Gene)) %>% group_by(Guide) %>% mutate(cell.count = n()) %>% summarise_at(c(colnames(omega),"cell.count"), mean)
        mat <- data.frame(avg.guide %>% select(-Guide, -cell.count), row.names = avg.guide$Guide) %>% as.matrix()
        ## guide.score <- sweep(mat, 2, vec, `/`) %>% log2() %>% as.data.frame()
        guide.score <- sweep.wrapper(mat,vec)
        
                                        # ## asinh on control topics
                                        # asinh.neg <- neg.omega %>% select(colnames(omega)) %>% asinh() %>% `rownames<-`(neg.omega$Row.names)
                                        # asinh.all <- ann.omega %>% select(colnames(omega)) %>% asinh() %>% `rownames<-`(ann.omega$Row.names)
        
        ## Wilcoxon Test
        gene.test <- get.p.value(ann.omega.filtered)
        wilcox.gene.score <- merge.score.with.test(gene.score, gene.test, "wilcox.p", p.value.thr=0.05, adj.p.value.thr=0.1)
        
        
                                        # test gene on per-guide averages [more conservative, 201209.7]
                                        # perguide.gene.score is basically the same as guide.score except it's annotated
        ## might not need this anymore! [210628]
        guide.score.ann <- guide.score %>% as.data.frame() %>% mutate(Guide=rownames(guide.score)) %>%
            merge(., avg.guide %>% select(Guide,cell.count), by="Guide")
        
                                        # # remove genes with only one guide (and guide with less than 5 cells [201216.2])
                                        # guide.score.ann.filtered <- guide.score.ann %>% subset(cell.count >=5) %>% group_by(Gene) %>% mutate(nCount = n()) %>% subset(nCount > 1) %>% select(-nCount) %>% as.data.frame()
                                        # avg.guide.ann <- merge(avg.guide %>% separate(col="Guide", into=c("Guide","Gene"), sep="::") %>% select(-Gene), barcode.names %>% select(Gene, Guide) %>% unique() %>% mutate(Guide.Gene = paste0(Guide,"::",Gene)), by.x="Guide", by.y="Guide.Gene") [210628, old]
        avg.guide.ann <- merge(avg.guide %>% mutate(Guide.Gene=Guide) %>% select(-Guide), barcode.names %>% select(Gene, Guide) %>% unique() %>% mutate(Guide.Gene = paste0(Guide,"::",Gene)), by="Guide.Gene")
                                        # avg.guide.ann.filtered <- avg.guide.ann %>% subset(cell.count >=5) %>% group_by(Gene) %>% mutate(nCount = n()) %>% subset(nCount > 1) %>% select(-nCount) %>% as.data.frame()
        
                                        # wilcoxon and t.test on normalized and unnormalized per guide average data
        perguide.gene.test <- get.p.value(avg.guide.ann)
        wilcox.perguide.gene.score <- merge.score.with.test(gene.score, perguide.gene.test, "wilcox.p", p.value.thr=0.05, adj.p.value.thr=0.1)
        
                                        # end of Wilcoxon Test section
        
        
        ## # for testing purposes 210105
        ## here.perguide.gene.test <- perguide.gene.test 
        ## here.gene.test <- gene.test
        ## here.wilcox.gene.score <- wilcox.gene.score
        ## here.wilcox.perguide.gene.score <- wilcox.perguide.gene.score 
        
        ## Perturbation z-score [210809]
        ptb.zscore <- apply(gene.score, 1, function(x) (x - mean(x)) / sd(x)) %>% t
        write.table(ptb.zscore, file=paste0(OUTDIRSAMPLE, "Perturbation.zscore.txt"), row.names=F, quote=F, sep="\t")

        save(omega, theta, ann.omega, ann.omega.filtered, avg.gene, avg.neg, gene.score,
             avg.guide, guide.score, gene.test, 
             wilcox.gene.score, perguide.gene.test, wilcox.perguide.gene.score, ptb.zscore,
             file=file.name) # 916 MB for n=1, K=2 # maybe put these in scratch space?
        
        
    } else {
        ## avg cell and avg ctrl before comparing (then log2() + heatmap)
        tmp <- avg.neg %>% group_by(grepl(all_of(rep1.label),Gene, fixed=TRUE) ) %>% summarise_at(colnames(omega), mean)
        colnames(tmp) <- colnames(avg.gene)
        avg.neg <- rbind(avg.neg,tmp)
        mat <- avg.gene %>% subset(grepl(rep1.label,Gene, fixed=TRUE)) %>% select(colnames(omega)) %>% as.matrix()
        vec <- avg.neg %>% subset(Gene==TRUE) %>% select(colnames(omega)) %>% as.matrix()
        ## gene.score <- sweep(mat, 2, vec, `/`) %>% log2() %>% as.data.frame()
        gene.score <- sweep.wrapper(mat,vec)
        mat <- avg.gene %>% subset(grepl(rep2.label,Gene, fixed=TRUE)) %>% select(colnames(omega)) %>% as.matrix()
        vec <- avg.neg %>% subset(Gene==FALSE) %>% select(colnames(omega)) %>% as.matrix()
        ## gene.score <- rbind(gene.score, sweep(mat, 2, vec, `/`) %>% log2() %>% as.data.frame())
        gene.score <- rbind(gene.score, sweep.wrapper(mat,vec))

        
        ## avg guide and avg ctrl before comparing (then log2() + heatmap)
        avg.guide <- ann.omega.filtered %>% mutate(Guide = paste0(Guide,"::",Gene)) %>% group_by(Guide) %>% mutate(cell.count = n()) %>% summarise_at(c(colnames(omega),"cell.count"), mean)
        mat <- avg.guide %>% column_to_rownames("Guide") %>% subset(grepl(rep2.label,rownames(.), fixed=TRUE)) %>% select(-cell.count) %>% as.matrix()
        guide.score <- sweep.wrapper(mat,vec)  # sweep(mat, 2, vec, `/`) %>% log2() %>% as.data.frame()
        mat <- avg.guide %>% column_to_rownames("Guide") %>% subset(grepl(rep1.label,rownames(.), fixed=TRUE)) %>% select(-cell.count) %>% as.matrix()
        vec <- avg.neg %>% subset(Gene==TRUE) %>% select(colnames(omega)) %>% as.matrix()
        guide.score <- rbind(guide.score, sweep.wrapper(mat,vec)) ## sweep(mat, 2, vec, `/`) %>% log2() %>% as.data.frame())
        
                                        # ## asinh on control topics
                                        # asinh.neg <- neg.omega %>% select(colnames(omega)) %>% asinh() %>% `rownames<-`(neg.omega$Row.names)
                                        # asinh.all <- ann.omega %>% select(colnames(omega)) %>% asinh() %>% `rownames<-`(ann.omega$Row.names)
        
        ## Wilcoxon Test      
        gene.test <- get.p.value(ann.omega.filtered)
        wilcox.gene.score <- merge.score.with.test(gene.score, gene.test, "wilcox.p")
                                        # test gene on per-guide averages [more conservative, 201209.7]
                                        # perguide.gene.score is basically the same as guide.score except it's annotated
        guide.score.ann <- guide.score %>% as.data.frame() %>% #mutate(Guide=rownames(guide.score)) %>%
            merge(avg.guide %>% select(Guide,cell.count), ., by.x="Guide", by.y=0) %>%
            separate(., col="Guide", into=c("Guide", "Gene"), sep="::")
        
                                        # # remove genes with only one guide (and guide with less than 5 cells [201216.2])
                                        # guide.score.ann.filtered <- guide.score.ann %>% subset(cell.count >=5) %>% group_by(Gene) %>% mutate(nCount = n()) %>% subset(nCount > 1) %>% select(-nCount) %>% as.data.frame()
        avg.guide.ann <- avg.guide %>% separate(col="Guide", into=c("Guide","Gene"), sep="::")
                                        # avg.guide.ann.filtered <- avg.guide.ann %>% subset(cell.count >=5) %>% group_by(Gene) %>% mutate(nCount = n()) %>% subset(nCount > 1) %>% select(-nCount) %>% as.data.frame()
        
                                        # wilcoxon and t.test on normalized and unnormalized per guide average data
        perguide.gene.test <- get.p.value(avg.guide.ann)
        wilcox.perguide.gene.score <- merge.score.with.test(gene.score, perguide.gene.test, "wilcox.p", p.value.thr=0.05, adj.p.value.thr=0.1)
        
                                        # end of Wilcoxon Test section
        ## Perturbation z-score [210809]
        ptb.zscore <- apply(gene.score, 1, function(x) (x - mean(x)) / sd(x)) %>% t
        write.table(ptb.zscore, file=paste0(OUTDIRSAMPLE, "Perturbation.zscore.txt"), row.names=F, quote=F, sep="\t")
        
        save(omega, theta, ann.omega, ann.omega.filtered, avg.gene, avg.neg, gene.score,
             avg.guide, guide.score, gene.test,
             wilcox.gene.score, perguide.gene.test, wilcox.perguide.gene.score, ptb.zscore,
             file=file.name) 
        
        
    }
    
    print("starting INT")
##### INT section
    ## rank norm on all cells # 210105 using filtered omega
    rankNorm.omega <- ann.omega.filtered %>% `rownames<-`(ann.omega.filtered$long.CBC) %>% select(all_of(topic.names)) %>% apply(2,rankNorm) %>% as.data.frame()
                                        # old <210105
                                        # rankNorm.omega <- omega %>% apply(2,RankNorm) %>% as.data.frame()
    
    ann.rankNorm.omega <- merge(barcode.names,rankNorm.omega,by.x="long.CBC", by.y=0)
                                        # old <210105
                                        # ann.rankNorm.omega.filtered <- ann.rankNorm.omega %>% group_by(Gene) %>% mutate(nCount = n()) %>% subset(nCount > 1) %>% select(-nCount) %>% as.data.frame()
    
    ## rank norm on all guide averages
    rankNorm.avg.guide <- avg.guide %>% column_to_rownames("Guide") %>% select(-cell.count) %>% apply(2,rankNorm) %>% as.data.frame()
    ann.rankNorm.avg.guide <- merge(avg.guide %>% select(Guide,cell.count) %>% unique(), rankNorm.avg.guide, by.x="Guide", by.y=0) %>%
        separate(., col="Guide", into=c("Guide", "Gene"), sep="::")
                                        # old <210105
                                        # ann.rankNorm.avg.guide.filtered <- ann.rankNorm.avg.guide %>% subset(cell.count >= 5) %>% group_by(Gene) %>% mutate(nCount = n()) %>% subset(nCount > 1) %>% select(-nCount) %>% as.data.frame()
    
    ## branching off for pooled #TODO: how to merge the two?
    if (SEP) {  
        INT.gene.test <- get.p.value(ann.rankNorm.omega, INT=T)
        INT.avg.guide.test <- get.p.value(ann.rankNorm.avg.guide, INT=T)
        
    } else {
        ## # get negative controls
        ## rankNorm.neg <- ann.rankNorm.omega %>% subset(Gene %in% c("negative-control", "safe-targeting")) %>% as.data.frame() # there must be a away to subset rows that contains the keywords and doens't have to be exact.
                                        # rankNorm.neg <- ann.rankNorm.omega.filtered %>% subset(Gene %in% c("negative-control", "safe-targeting")) %>% as.data.frame() # there must be a away to subset rows that contains the keywords and doens't have to be exact.
        
                                        # t-test on INT
##################################################
        ## local function for INT.ttest
        get.p.value <- function(ann.rankNorm.omega.filtered){
            INT.gene.test.list <- vector("list", k) # iterate through topics
            for ( t in 1:length(topic.names) ) {
                col <- topic.names[t] # current topic
                ctrl <- ann.rankNorm.omega.filtered %>% subset(grepl("control$|targeting$", Gene)) %>% select(all_of(col)) %>% as.matrix() %>% as.numeric()
                result <- ann.rankNorm.omega.filtered %>% group_by(Gene) %>% summarise(ttest.p = t.test(x = get(col) %>% as.numeric(), y = ctrl)$p.value)
                result$Topic <- rep(col, nrow(result))
                INT.gene.test.list[[t]] <- result
            }
            INT.gene.test <- do.call(rbind, INT.gene.test.list)
        }
##################################################
        INT.gene.test <- get.p.value(ann.rankNorm.omega)
        INT.avg.guide.test <- get.p.value(ann.rankNorm.avg.guide)
                                        # INT.gene.test <- get.p.value(ann.rankNorm.omega.filtered)
                                        # INT.avg.guide.test <- get.p.value(ann.rankNorm.avg.guide.filtered)
        
    }
    
    print("done with INT")
    
    ## merge "pooled" back in
                                        # INT.gene.test.mtx <- INT.gene.test %>% spread(key=Topic, value=ttest.p) # stretch this out into a matrix of ( gene x topic )
    INT.gene.test <- INT.gene.test %>% mutate(adjusted.ttest.p = p.adjust(ttest.p, method = "fdr"))
    INT.avg.guide.test <- INT.avg.guide.test %>% mutate(adjusted.ttest.p = p.adjust(ttest.p, method = "fdr"))
    
    INT.ttest.gene.score <- merge.score.with.test(gene.score, INT.gene.test, "ttest.p")
    INT.ttest.avg.guide.score <- merge.score.with.test(gene.score, INT.avg.guide.test, "ttest.p", p.value.thr=0.05, adj.p.value.thr=0.1)
                                        # end of INT t-test section
    
    
                                        # # for testing purposes 210105
                                        # here.INT.gene.test <- INT.gene.test
                                        # here.INT.avg.guide.test <- INT.avg.guide.test
                                        # here.INT.ttest.gene.score <- INT.ttest.gene.score 
                                        # here.INT.ttest.avg.guide.score <- INT.ttest.avg.guide.score
                                        # here.ann.rankNorm.omega <- ann.rankNorm.omega
                                        # here.ann.rankNorm.avg.guide <- ann.rankNorm.avg.guide 
    
    
    ## add individual log2fc against average of control ## add sep
    
    if(SEP) {
        num.reps=2
        labels <- c(rep1.label, rep2.label) # expand to accomodate multiple samples
        fc.omega.list <- vector("list", num.reps) # iterate through topics + reps
        for (label.index in 1:num.reps){
            label.here <- labels[label.index]
            ctrl.vec <- ann.omega.filtered %>% subset(grepl(label.here, Gene, fixed=T) & grepl(paste0("^safe-targeting|^negative-control"), Gene)) %>% select(all_of(topic.names)) %>% apply(2,mean) %>% as.array()
            omega.filtered <- ann.omega.filtered %>% `rownames<-`(ann.omega.filtered$long.CBC) %>% subset(grepl(label.here, Gene, fixed=T)) %>% select(all_of(topic.names)) #210105
            fc.omega.list[[label.index]] <- sweep(omega.filtered,2,ctrl.vec,"/") # use this to plot # and save this!
        }
        fc.omega <- do.call(rbind, fc.omega.list)
    } else {
        ctrl.vec <- ann.omega.filtered %>% subset(grepl(paste0("^safe-targeting|^negative-control"), Gene)) %>% select(all_of(topic.names)) %>% apply(2,mean) %>% as.array()
        omega.filtered <- ann.omega.filtered %>% `rownames<-`(ann.omega.filtered$long.CBC) %>% select(all_of(topic.names)) #210105
        fc.omega <- sweep(omega.filtered,2,ctrl.vec,"/") # use this to plot # and save this!
                                        # old <210105
                                        # fc.omega <- sweep(omega,2,ctrl.vec,"/") # use this to plot # and save this!
    }
    
    log2fc.omega <- log2(fc.omega)
    tmp <- log2fc.omega %>% `colnames<-`(topic.names) # keep log2fc.omega colnames same as omega, alter tmp for merging with guide information
    log2fc.ann.omega <- merge(ann.omega.filtered %>% select(long.CBC, CBC, Gene, Guide), tmp, by.x="long.CBC", by.y=0)
    fc.ann.omega <- merge(ann.omega.filtered %>% select(long.CBC, CBC, Gene, Guide), fc.omega %>% `colnames<-`(topic.names), by.x="long.CBC", by.y=0) %>% `rownames<-`(.$long.CBC)
                                        # old <210105
                                        # log2fc.ann.omega <- merge(ann.omega %>% select(Row.names,Gene, Guide), tmp, by.x="Row.names", by.y=0)
                                        # fc.ann.omega <- merge(ann.omega %>% select(Row.names, Gene, Guide), fc.omega %>% `colnames<-`(topic.names), by.x="Row.names", by.y=0)
    
    
                                        # # for testing purposes 210105
                                        # here.log2fc.ann.omega <-log2fc.ann.omega
                                        # here.fc.ann.omega <- fc.ann.omega
    

    ## save the results
    if (SEP) {  # removed rankNorm.neg
        save(omega, theta, theta.zscore, ann.omega, ann.omega.filtered, avg.gene, avg.neg, gene.score, avg.guide, 
             guide.score, gene.test, wilcox.gene.score,
             perguide.gene.test, wilcox.perguide.gene.score, ptb.zscore,
             ## ann.GO.omega, topic.names, avg.GO, GO.score, GO.test, GO.test.mtx, wilcox.GO.score,
             ann.rankNorm.omega, INT.gene.test, INT.ttest.gene.score, 
             INT.avg.guide.test, INT.ttest.avg.guide.score,
             log2fc.omega, log2fc.ann.omega, fc.omega, fc.ann.omega,
             file=file.name) 
    } else {
        save(omega, theta, theta.zscore,
             ann.omega, ann.omega.filtered, avg.gene, avg.neg, gene.score, avg.guide, 
             guide.score, gene.test, wilcox.gene.score,
             perguide.gene.test, wilcox.perguide.gene.score, ptb.zscore,
             ## ann.GO.omega, topic.names, avg.GO, GO.score, GO.test, GO.test.mtx, wilcox.GO.score,
             ann.rankNorm.omega, INT.gene.test, INT.ttest.gene.score, 
             INT.avg.guide.test, INT.ttest.avg.guide.score,
             log2fc.omega, log2fc.ann.omega, fc.omega, fc.ann.omega,
             file=file.name)
    }
}




# if(opt$subsample.type!="ctrl") {

##########################################################################
## aggregate all statistical test results
all.test <- rbind(gene.test %>% select(Gene, Topic, adjusted.wilcox.p, wilcox.p) %>% `colnames<-`(c("Gene","Topic","adjusted.p.value","p.value")) %>% mutate(test.type="per.cell.wilcoxon"),
                  INT.gene.test %>% select(Gene, Topic, adjusted.ttest.p, ttest.p) %>% `colnames<-`(c("Gene","Topic","adjusted.p.value","p.value")) %>% mutate(test.type="per.cell.INT.ttest"),
                  perguide.gene.test %>% select(Gene, Topic, adjusted.wilcox.p, wilcox.p) %>% `colnames<-`(c("Gene","Topic","adjusted.p.value","p.value")) %>% mutate(test.type="per.guide.wilcoxon"),
                  INT.avg.guide.test %>% select(Gene, Topic, adjusted.ttest.p, ttest.p) %>% `colnames<-`(c("Gene","Topic","adjusted.p.value","p.value")) %>% mutate(test.type="per.guide.INT.ttest")) %>%
  mutate(gene.type=ifelse(grepl(paste0(c(non.expressed.genes, "targeting$|control$"), collapse="|"), Gene), "non-expressed", "expressed")) %>%
  mutate(gene.type=ifelse(grepl("targeting$|control$", Gene), "control", gene.type))
write.table(all.test, file=paste0(OUTDIRSAMPLE, "/all.test.", SUBSCRIPT, ".txt"),
            sep="\t", row.names=F, col.names=T, quote=F)




# ##########################################################################
# ## empirical fdr using the non-expressed perturbation results
# # [false discovery rate = # events above your p-value threshold that are not real signal / # events above your p-value threshold]

#     fdr = function(pval, realPvals, nullPvals) min(1, (sum(nullPvals<=pval)/length(nullPvals)) / (sum(realPvals<=pval)/length(realPvals)))
#     for ( n in 1:length(SAMPLE) ) {
# # all.test <- read.delim(file=paste0(OUTDIR,SAMPLE[n], "/all.test.K", k, ".txt"), header=T, stringsAsFactors=F)
#         realPvals.df.list <- vector("list", length(unique(all.test$test.type)))
#         for (j in 1:length(unique(all.test$test.type))) {
#             test.name<-unique(all.test$test.type)[j]
#             df.here <- all.test %>% subset(test.type==test.name)
#             realPvals <- df.here %>% subset(gene.type=="expressed") %>% select(p.value) %>% unlist() %>% as.numeric()
#             nullPvals <- df.here %>% subset(gene.type=="non-expressed") %>% select(p.value) %>% unlist() %>% as.numeric()
#             fdr.output <- do.call(rbind, lapply(1:length(realPvals), function(i) fdr(realPvals[i], realPvals, nullPvals)))
#             realPvals.df.list[[j]] <- df.here %>% subset(gene.type=="expressed") %>% mutate(fdr=fdr.output)
#         }
#         realPvals.df <- do.call(rbind, realPvals.df.list)
#         write.table(realPvals.df, file=paste0(OUTDIRSAMPLE, "/all.expressed.genes.pval.fdr.",SUBSCRIPT,".txt"), sep="\t", row.names=F, col.names=T, quote=F)
#     }




# ##########################################################################
# ## CAD GWAS gene annotation
# ## GWAS.annotation.df <- ref.table %>% select(`Symbol`, chromosome, `TSS position (strongest)`, `Top SNP ID`, `TSS v. SNP loc`, GWAS.classification, gene.count.per.GWAS.category, TSS.v.SNP.ranking, TSS.v.SNP.ranking.in.GWAS.category, total.gene.in.this.loci.in.GWAS.category, total.TSS.v.SNP.ranking.count, total.TSS.v.SNP.ranking.count.per.GWAS.classification)
# GWAS.annotation.df <- ref.table
# all.test.GWAS <- merge(all.test, GWAS.annotation.df, by.x="Gene", by.y="Symbol", all=T) ## HERE: add total rank count even if that rank doesn't have significant gene
# test.condition.names <- all.test$test.type %>% unique()
# ## count.by.GWAS.list <- count.by.GWAS.withTopic.list <- vector("list",length(test.condition.names))
# ## for (i in 1:length(test.condition.names)) {
#     ## test.condition <- test.condition.names[i]
#     p.value.thr <- 0.1
#     passed.filter.GWAS <- all.test.GWAS %>% subset(adjusted.p.value < p.value.thr)
#     count.by.GWAS <- passed.filter.GWAS %>% select(-Topic, -adjusted.p.value, -p.value) %>% unique() %>%
#         group_by(test.type, GWAS.classification) %>%
#         mutate(pass.filter.GWAS.classification.count = n(),
#                GWAS.class.enrichment = pass.filter.GWAS.classification.count / gene.count.per.GWAS.category) %>% ungroup() %>%
#         group_by(test.type, GWAS.classification, TSS.v.SNP.ranking.in.GWAS.category) %>%
#         mutate(passed.filter.ranking.count = n()) %>% ungroup()
#     count.by.GWAS.withTopic <- passed.filter.GWAS %>% select(-adjusted.p.value, -p.value) %>% unique() %>%
#         group_by(test.type, GWAS.classification) %>%
#         mutate(pass.filter.GWAS.classification.count = n(),
#                GWAS.class.enrichment = pass.filter.GWAS.classification.count / gene.count.per.GWAS.category) %>% ungroup() %>%
#         group_by(test.type, GWAS.classification, TSS.v.SNP.ranking.in.GWAS.category) %>%
#         mutate(passed.filter.ranking.count = n()) %>% ungroup()
#         ## count.by.GWAS.list[[i]] <- passed.filter.GWAS %>% select(-Topic, -adjusted.p.value, -p.value) %>% unique() %>% group_by(GWAS.classification,TSS.v.SNP.ranking.in.GWAS.category) %>% mutate(ranking.count = n()) %>% ungroup() %>% group_by(GWAS.classification) %>% mutate(pass.filter.GWAS.classification.count = n(), GWAS.class.enrichment = pass.filter.GWAS.classification.count / total.gene.in.this.loci.in.GWAS.category)
#     ## no.ranking.index <- which(count.by.GWAS.list[[i]]$`TSS.v.SNP.ranking` %>% is.na())
#     ## count.by.GWAS.list[[i]]$ranking.count[no.ranking.index] <- NA
#     ## tmp <- data.frame(table(passed.filter.GWAS$GWAS.classification)) %>% `colnames<-`(c("GWAS.classification","count")) %>% mutate(test.type=test.condition)
# # }
# ## count.by.GWAS <- do.call(rbind, count.by.GWAS.list)
# ## count.by.GWAS.withTopic <- do.call(rbind, count.by.GWAS.withTopic.list)
# write.table(count.by.GWAS,file=paste0(OUTDIRSAMPLE,"/count.by.GWAS.classes_p.adj.",p.value.thr %>% as.character,"_",SUBSCRIPT,".txt"), sep="\t", col.names=T, row.names=F, quote=F)
# write.table(count.by.GWAS.withTopic,file=paste0(OUTDIRSAMPLE,"/count.by.GWAS.classes.withTopic_p.adj.",p.value.thr %>% as.character,"_",SUBSCRIPT,".txt"), sep="\t", col.names=T, row.names=F, quote=F)

# } ## if this sample does not contain only the controls


# ##########################################################################
# ## KL score for genes specific to topic
# geneName <- theta %>% row.names()
# CBC <- omega %>% row.names()

# # select genes that are most critical for separating topics
# features <- ExtractTopFeatures(theta, top_features = length(geneName), method = "poisson", options = "min", shared=TRUE)
# features$genes <- apply(features$indices, c(1,2), function(x) geneName[x])

#                                         # topFeatures <- features
#                                         # topFeatures$scores <- features$scores[,1:50]
#                                         # topFeatures$indices <- features$indices[,1:50]
#                                         # topFeatures$genes <- features$genes[,1:50]
#                                         #
#                                         # # 100 * x for topFeatures$scores
#                                         # topFeatures$scores <- topFeatures$scores * 100

#                                         # store to file
# toSave.features.list <- vector("list", nrow(features$genes))
# for (i in 1:nrow(features$genes)) {
#     tmp <- data.frame(genes = features$genes[i,],
#                       scores = features$scores[i,],
#                       topic = rep(i, ncol(features$genes)),
#                       indices = features$indices[i,],
#                       sample = rep(SAMPLE[n],ncol(features$genes)),
#                       K = rep(k, ncol(features$genes)))
#     toSave.features.list[[i]] <- tmp
# }
# toSave.features <- do.call(rbind, toSave.features.list)
# write.table(toSave.features, file = paste0(OUTDIRSAMPLE, "/topic.KL.score_K", k, ".dt_", DENSITY.THRESHOLD, ifelse(SEP, ".sep", ""), ".txt"),
#             sep='\t', quote=F, col.names=T, row.names=F)




# ##########################################################################
# ## Plots

# # replicates on labels 
# if(SEP){
#   # labels
#   tmp.labels <- rownames(gene.score) %>% unique() %>% strsplit("-") %>% sapply("[[",2) %>% unique()
#   tmp.labels <- tmp.labels[!(tmp.labels %in% c("control","targeting"))]
#   rep.labels <- paste0("-", tmp.labels)
#   n.rep <- length(rep.labels)
#   rep1.label <- paste0("-",tmp.labels[1])
#   rep2.label <- paste0("-",tmp.labels[2])
  
#   add.tag <- function(df, rep.labels, col="Gene"){
#     n.rep <- length(rep.labels)
#     if(class(df) == "character") {
#       df.list <- vector("list", n.rep)
#       for(rep.index in 1:n.rep) {
#         df.list[[rep.index]] <- paste0(df, rep.labels[rep.index])
#       }
#       return(do.call(append, df.list))
#     } else {
#       df.list <- vector("list", n.rep)
#       for(rep.index in 1:n.rep) {
#         tmp.df <- df
#         tmp.df[[col]] <- paste0(tmp.df[[col]], rep.labels[rep.index])
#         df.list[[rep.index]] <- tmp.df #df %>% mutate(!!column:=col = paste0(!!as.symbol(col), rep.labels[rep.index]))
#       }
#       return(do.call(rbind, df.list))
#     }
#   }
  
#   # label reference files
#   gene.set.type.df <- add.tag(gene.set.type.df, rep.labels)
#   gene.set <- add.tag(gene.set, rep.labels)
#   EC.pos.ctrl.gene.set <- add.tag(EC.pos.ctrl.gene.set, rep.labels)
#   CAD.focus.gene.set <- add.tag(CAD.focus.gene.set, rep.labels)
#   enhancer.set <- add.tag(enhancer.set, rep.labels)
  
# }


# ##########################################################################
# ## topic gene raw weight list
# pdf(file=paste0(FIGDIRTOP,"topGeneInTopics.zscore.pdf"), width=4, height=6)
# topFeatures.raw.weight <- theta.zscore %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene", variable.name="topic", value.name="scores") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:50)
# for ( t in 1:dim(theta)[2] ) {
#     toPlot <- data.frame(Gene=topFeatures.raw.weight %>% subset(topic == t) %>% pull(Gene),
#                          Score=topFeatures.raw.weight %>% subset(topic == t) %>% pull(scores))
#     p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score*100) ) + geom_col() + theme_minimal()
#     p <- p + coord_flip() + xlab("Top 50 Genes") + ylab("Raw Weights (z-score)") + ggtitle(paste(SAMPLE[n], ", Topic ", t, sep="")) + mytheme
#     print(p)
# }
# dev.off()


# ##########################################################################
# ## top expressed genes per topic by KL specificity score list
# pdf(file=paste0(FIGDIRTOP, "topGeneInTopics.KL.pdf"), width=4, height=6)
#                                         # toSave.features <- read.delim(paste0(OUTDIR, SAMPLE[n], "/topic.KL.score_K", k, ".txt"), header=T, stringsAsFactors=F)
# topFeatures <- toSave.features %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:50)
# for ( t in 1:dim(theta)[2] ) {
#     toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
#                          Score=topFeatures %>% subset(topic == t) %>% pull(scores))
#     p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score*100) ) + geom_col() + theme_minimal()
#     p <- p + coord_flip() + xlab("Top 50 Genes") + ylab("KL score (gene specific to this topic)") + ggtitle(paste(SAMPLE[n], ", K = ", k, ", Topic ", t, sep="")) + mytheme
#     print(p)

# }
# dev.off()


# ##########################################################################
# ## Topic's top gene list, ranked by raw weight
# pdf(file=paste0(FIGDIRTOP,"topGeneInTopics.rawWeight.pdf"), width=4, height=6)
# topFeatures <- theta %>% as.data.frame() %>% mutate(genes=rownames(.)) %>% melt(id.vars="genes",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:50)
# for ( t in 1:dim(theta)[2] ) {
#     toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
#                          Score=topFeatures %>% subset(topic == t) %>% pull(scores))
#     p <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score*100) ) + geom_col() + theme_minimal()
#     p <- p + coord_flip() + xlab("Top 50 Genes") + ylab("Raw Score (gene's weight in topic)") + ggtitle(paste(SAMPLE[n], ", Topic ", t, sep="")) + mytheme
#     print(p)
# }
# dev.off()


# ##########################################################################
# ## KL score list with annotataion
# pdf(file=paste0(FIGDIRTOP,"topGeneInTopics.annotated.KL.pdf"), width=4.5, height=5)
# topFeatures <- toSave.features %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:10)
# for ( t in 1:dim(theta)[2] ) {
#     toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
#                          Score=topFeatures %>% subset(topic == t) %>% pull(scores)) %>%
#         merge(., gene.def.pathways, by="Gene", all.x=T)
#     toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
#     p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score*100, fill=Pathway) ) + geom_col(width=0.5) + theme_minimal()
#     p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("KL score (gene specific to this topic)") +
#         mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ", K = ", k, ", Topic ", t))
#     print(p4)
# }
# dev.off()


# ##########################################################################
# ## raw program TPM list with annotataion
# pdf(file=paste0(FIGDIRTOP,"topGeneInTopics.annotated.TPM.pdf"), width=4.5, height=5)
# topFeatures <- theta %>% as.data.frame() %>% mutate(genes=rownames(.)) %>% melt(id.vars="genes",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:10)
# for ( t in 1:dim(theta)[2] ) {
#     toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
#                          Score=topFeatures %>% subset(topic == t) %>% pull(scores)) %>%
#         merge(., gene.def.pathways, by="Gene", all.x=T)
#     toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
#     p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score*1000000, fill=Pathway) ) + geom_col(width=0.5) + theme_minimal()
#     p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("Raw Weight (in TPM)") +
#         mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ", K = ", k, ", Topic ", t))
#     print(p4)
# }
# dev.off()



# ##########################################################################
# ## raw program zscore list with annotataion
# pdf(file=paste0(FIGDIRTOP,"topGeneInTopics.annotated.zscore.pdf"), width=4.5, height=5)
# topFeatures <- theta.zscore %>% as.data.frame() %>% mutate(genes=rownames(.)) %>% melt(id.vars="genes",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:10)
# for ( t in 1:dim(theta)[2] ) {
#     toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
#                          Score=topFeatures %>% subset(topic == t) %>% pull(scores)) %>%
#         merge(., gene.def.pathways, by="Gene", all.x=T)
#     toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
#     p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score, fill=Pathway) ) + geom_col(width=0.5) + theme_minimal()
#     p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("z-score") +
#         mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ", K = ", k, ", Topic ", t))
#     print(p4)
# }
# dev.off()





# ##########################################################################
# ## Perturbation zscore list with annotataion
# pdf(file=paste0(FIGDIRTOP,"Perturbation_zscore.annotated.pdf"), width=4.5, height=5)
# topFeatures <- ptb.zscore %>% as.data.frame() %>% mutate(genes=rownames(.)) %>% melt(id.vars="genes",value.name="scores", variable.name="topic") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:10) %>% mutate(topic = gsub("topic_","", topic))
# for ( t in 1:dim(theta)[2] ) {
#     toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
#                          Score=topFeatures %>% subset(topic == t) %>% pull(scores)) %>%
#         merge(., gene.def.pathways, by="Gene", all.x=T)
#     toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
#     p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score, fill=Pathway) ) + geom_col(width=0.5) + theme_minimal()
#     p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("Perturbation z-score") +
#         mytheme + theme(legend.position="bottom", legend.direction="vertical") + ggtitle(paste0(SAMPLE, ", K = ", k, ", Topic ", t))
#     print(p4)
# }
# dev.off()


# if(opt$subsample.type!="ctrl") {
#     fdr.thr <- 0.1
#     pdf(file=paste0(FIGDIRTOP, "Perturbation.zscore.sig.list.pdf"), width=6, height=6)
#     ptb.zscore.long <- ptb.zscore %>% as.data.frame %>% mutate(Gene=rownames(.)) %>% melt(id.vars = "Gene", value.name = "perturbation.zscore", variable.name = "Topic")
#     for ( topic in colnames(ptb.zscore) ) {
#         t <- gsub("topic_","",topic)
#         toPlot.all.test <- all.test %>% subset(test.type=="per.cell.wilcoxon" & Topic==topic)
#         toPlot.fdr <- realPvals.df %>% subset(test.type=="per.cell.wilcoxon" & Topic == topic) %>% select(Gene,fdr)##here210809
#                                         # assemble toPlot
#         toPlot <- ptb.zscore.long %>% subset(Topic == topic) %>% merge(.,toPlot.all.test,by=c("Gene","Topic"), all.x=T) %>%
#             merge(.,toPlot.fdr,by="Gene", all.x=T) %>%
#             merge(.,gene.set.type.df,by="Gene", all.x=T) %>% ##here210809
#             ## merge(.,gene.def.pathways, by="Gene", all.x=T) %>% 
#             merge(., ref.table %>% select("Symbol", "TSS.dist.to.SNP", "GWAS.classification"), by.x="Gene", by.y="Symbol", all.x=T) %>%
#             mutate(EC_ctrl_text = ifelse(.$GWAS.classification == "EC_ctrls", "(+)", "")) %>%
#             mutate(GWAS.class.text = ifelse(grepl("CAD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb"),
#                                      ifelse(grepl("IBD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb_IBD"), ""))) %>%
#             mutate(ann.Gene = paste0(Gene, GWAS.class.text, EC_ctrl_text))

#         toPlot <- toPlot %>% mutate(significant=ifelse((adjusted.p.value >= fdr.thr | is.na(adjusted.p.value)), "", "*"))

#         toPlot.top <- toPlot %>% arrange(desc(perturbation.zscore)) %>% slice(1:25)
#         toPlot.bottom <- toPlot %>% arrange(perturbation.zscore) %>% slice(1:25)
#         toPlot.extreme <- rbind(toPlot.top, toPlot.bottom) %>%
#             mutate(color=ifelse(grepl("CAD", type), "red",
#                          ifelse(type=="non-expressed", "gray",
#                          ifelse(type=="EC_ctrls", "blue", "black"))))  %>%
#             mutate(color=ifelse(is.na(type), "black", color))
#         ## colors <- toPlot.extreme$color[order(toPlot.extreme %>% arrange(desc(perturbation.zscore)) %>% pull(color))]
#         toPlot.extreme <- toPlot.extreme %>% arrange(perturbation.zscore)
#         ## add gene distance to CAD
#         toPlot.extreme$ann.Gene <- factor(toPlot.extreme$ann.Gene, levels = toPlot.extreme$ann.Gene)
#         p <- toPlot.extreme %>%  #mutate(Gene = paste0("<span style = 'color: ", color, ";'>", Gene, "</span>")) %>%
#             ggplot(aes(x=ann.Gene, y=perturbation.zscore, fill=significant)) + geom_col() + theme_minimal() +
#             coord_flip() + xlab("Most Extreme Gene (Perturbation)") + ylab("Perturbation z-score") + ggtitle(paste(SAMPLE, " perturbations, ", topic)) +
#             scale_fill_manual(values=c("grey", "#38b4f7")) +
#             geom_text(aes(label = significant)) +
#             theme(legend.position = "none", axis.text.y = element_text(colour = toPlot.extreme$color))
#         print(p)
#     }
#     dev.off()
# } ## test for control

# ## volcano plots
# volcano.plot <- function(toplot, ep.type, ranking.type) {
#     label <- toplot %>% subset(-log10(p.adjust) > 1) %>% mutate(motif.toshow = gsub("HUMAN.H11MO.", "", motif))
#     t <- gsub("topic_", "", toplot$topic[1])
#     p <- toplot %>% ggplot(aes(x=enrichment.log2fc, y=-log10(p.adjust))) + geom_point(size=0.5) + mytheme +
#         ggtitle(paste0(SAMPLE[1], " Topic ", t, " Top 100 ", ranking.type," ", ifelse(ep.type=="promoter", "Promoter", "Enhancer"), " Motif Enrichment")) + xlab("Motif Enrichment (log2FC)") + ylab("-log10(adjusted p-value)") +
#         geom_text_repel(data=label, box.padding = 0.5,
#                         aes(label=motif.toshow), size=5,
#                         color="black")
#     print(p)
#     p <- toplot %>% ggplot(aes(x=enrichment.log2fc, y=-log10(p.value))) + geom_point(size=0.5) + mytheme +
#         ggtitle(paste0(SAMPLE[1], " Topic ", t, " Top 100 ", ranking.type," ", ifelse(ep.type=="promoter", "Promoter", "Enhancer"), " Motif Enrichment")) + xlab("Motif Enrichment (log2FC)") + ylab("-log10(p-value)") +
#         geom_text_repel(data=label, box.padding = 0.5,
#                         aes(label=motif.toshow), size=5,
#                         color="black")
#     return(p)
# }

# ## function for all volcano plots
# all.volcano.plots <- function(all.fisher.df, ep.type, ranking.type) {
#     for ( t in 1:k ){
#         toplot <- all.fisher.df %>% subset(topic==paste0("topic_",t))
#         volcano.plot(toplot, ep.type, ranking.type) %>% print()
#     }
# }


# if(opt$subsample.type!="ctrl"){
# ##########################################################################
# ## q-q plot
# pdf(file=paste0(FIGDIRTOP,"p-value.qqplot.pdf"), width=8, height=8)
# for (test.name in unique(all.test$test.type)) {
#   toPlot <- all.test %>% subset(test.type==test.name)
#   toPlot$p.value <- -1*log10(toPlot$p.value)
#   min.value <- min(toPlot %>% subset(gene.type=="expressed") %>% select(p.value) %>% unlist() %>% as.numeric(),
#                    toPlot %>% subset(gene.type=="non-expressed") %>% select(p.value) %>% unlist() %>% as.numeric())
#   max.value <- 10
#   toPlot$p.value[toPlot$p.value < 10^-10] <- 10^-10
#   exp.gene.count <- toPlot %>% subset(gene.type=="expressed") %>% select(Gene) %>% unique() %>% unlist() %>% length()
#   nonexp.gene.count <- toPlot %>% subset(gene.type=="non-expressed") %>% select(Gene) %>% unique() %>% unlist() %>% length()


#   match.df <- toPlot %>% subset(adjusted.p.value < 0.1 & gene.type=="expressed")
#   my.qqplot(y = toPlot %>% subset(gene.type=="expressed") %>% select(p.value) %>% unlist() %>% as.numeric(),
#             x = toPlot %>% subset(gene.type=="non-expressed") %>% select(p.value) %>% unlist() %>% as.numeric(),
#             xlimit = c(min.value, max.value), ylimit = c(min.value, max.value),
#             ylab = paste0("Expressed Genes (-log10 p-value, n=[", exp.gene.count, " genes])"),
#             xlab = paste0("Control - Non-Expressed Genes (-log10 p-value, n=[", nonexp.gene.count, " genes])"),
#             main = paste0(SAMPLE[n], ", K = ", k, ", ", test.name),
#             match=T, match.y=T, match.df=match.df
#   )

# }
# dev.off()


# ##########################################################################
# ## empirical fdr vs p.adjust (BH) plot
# pdf(paste0(FIGDIRTOP,"empirical.fdr.vs.p.adjust.pdf"), width=8, height=6)
# for (j in 1:length(unique(realPvals.df$test.type))) {
#   test.name<-unique(all.test$test.type)[j]
#   toPlot <- realPvals.df %>% subset(test.type==test.name)
#   p <- toPlot %>% ggplot(aes(x=adjusted.p.value, y=fdr)) + geom_point(size=0.1) + geom_abline(color="red") +
#     mytheme + xlab("Adjusted p-value") + ylab("Empirical False Discovery Rate") + ggtitle(paste0(SAMPLE[n], ", ", test.name)) + coord_fixed()
#   print(p)
# }
# dev.off()


# ##########################################################################
# ## empirical FDR heatmaps
# for(fdr.method in c("empirical.fdr", "p.adjust")) {
# pdf(file=paste0(FIGDIRTOP, fdr.method, ".sig.ptbd.gene_fill.log2fc_heatmap.pdf"), width=12, height=6)
# for (emp.fdr.thr in c(0.05, 0.1, 0.25)) {
#   for (current.test.type in realPvals.df$test.type %>% unique()) {
#       if(fdr.method == "empirical.fdr") {
#           test.list <- realPvals.df %>% subset(test.type==current.test.type) %>% mutate(adjusted.p.value=fdr)
#       } else {
#           test.list <- realPvals.df %>% subset(test.type==current.test.type)
#       }
#     genes.toInclude <- test.list %>% subset(adjusted.p.value < emp.fdr.thr) %>% pull(Gene) %>% unique()
#     toPlot <- gene.score %>% subset(., rownames(gene.score) %in% genes.toInclude)
#     toPlot <- add.snp.gene.info(toPlot, type="rownames")
#     # plot heatmap
#     cols <- rep('black', nrow(toPlot))
#     #turn red the specified rows in tf
#     cols[row.names(toPlot) %in% (gene.set.type.df %>% subset(grepl("EC_ctrls", type)) %>% pull(Gene))] <- "blue"
#     cols[row.names(toPlot) %in% (gene.set.type.df %>% subset(grepl("CAD_Loci_all", type)) %>% pull(Gene))] <- "red"
#     # cols[row.names(toPlot) %in% gene.set] <- "red"
#     rownames(toPlot)[row.names(toPlot) %in% gene.set] <- paste0("[ ", rownames(toPlot)[which(row.names(toPlot) %in% gene.set)], " ]")

#     if(nrow(toPlot) > 1) {
#         ## plotHeatmap( toPlot, cellNote=NULL, rownames(toPlot), title=title, colCol=cols)

#         toHighlight.asterisk <- merge.score.with.test(gene.score, test.list, test.col.name="p.value", p.value.thr=1, adj.p.value.thr=emp.fdr.thr, fill.all=T, fill="", overlay=T)$score.mtx
#         title=paste0(SAMPLE[n],", K = ", k, ", ", current.test.type, ", \n", ifelse(fdr.method == "empirical.fdr", "empirical fdr", "BH adjusted p-value"), " < ", emp.fdr.thr, ", number of significant genes = ", nrow(toPlot)) 
#       plotHeatmap( toPlot, cellNote=toHighlight.asterisk, rownames(toPlot), title=title, colCol=cols)

#       toHighlight.value <- merge.score.with.test(gene.score, test.list, test.col.name="p.value", p.value.thr=1, adj.p.value.thr=emp.fdr.thr, fill.all=T, fill="", overlay=F, num.thr = 1)$score.mtx
#       plotHeatmap( toPlot, cellNote=toHighlight.value, rownames(toPlot), title=title, colCol=cols)

#       toHighlight.value <- merge.score.with.test(gene.score, test.list, test.col.name="p.value", p.value.thr=1, adj.p.value.thr=emp.fdr.thr, fill.all=T, fill="", overlay=F)$score.mtx
#       plotHeatmap( toPlot, cellNote=toHighlight.value, rownames(toPlot), title=title, colCol=cols)
#     }
#   }
# }
# dev.off()
# }

    
#     ##########################################################################
#     ## full heatmap
# pdf(file = paste0(FIGDIRTOP,"Gene.full.heatmap.pdf"), width=36, height=6) 

#     current.test.type <- "per.guide.wilcoxon"
#     toPlot <- gene.score
#     toPlot <- add.snp.gene.info(toPlot, type="rownames")
#     cols <- rep('black', nrow(toPlot))
#                                         #turn red the specified rows in tf
#     cols <- rep('black', nrow(toPlot))
#                                         #turn red the specified rows in tf
#     cols[row.names(toPlot) %in% (gene.set.type.df %>% subset(grepl("EC_ctrls", type)) %>% pull(Gene))] <- "blue"
#     cols[row.names(toPlot) %in% (gene.set.type.df %>% subset(grepl("CAD_Loci_all", type)) %>% pull(Gene))] <- "red"
 
#     rownames(toPlot)[row.names(toPlot) %in% gene.set] <- paste0("[ ", rownames(toPlot)[which(row.names(toPlot) %in% gene.set)], " ]")

#     title=paste0(SAMPLE[n],", K = ", k)
#     plotHeatmap(toPlot, cellNote=NULL, rownames(toPlot), title=title, colCol=cols)

#     for(emp.fdr.thr in c(0.05, 0.1, 0.25)){
#                                         # add asterisks for significant perturbation/topic
#                                         # test.list <- realPvals.df %>% subset(test.type==current.test.type) %>% mutate(adjusted.p.value=fdr)
#                                         # toHighlight.asterisk <- merge.score.with.test(gene.score, test.list, test.col.name="p.value", p.value.thr=1, adj.p.value.thr=emp.fdr.thr, fill.all=T, fill="", overlay=T)$score.mtx
#                                         # title=paste0(SAMPLE[n],", K = ", k, ", ", current.test.type, ", \nempirical fdr < ", emp.fdr.thr)
#                                         # plotHeatmap(toPlot, cellNote=toHighlight.asterisk, rownames(toPlot), title=title, colCol=cols)
#     }
#     dev.off()


    
    
# ##########################################################################
# ## Lists of genes or perturbations to understand topics
# all.test.guide.w <- all.test %>% subset(test.type=="per.guide.wilcoxon")
# realPvals.df.guide.w <- realPvals.df %>% subset(test.type=="per.guide.wilcoxon")
# # make a plot for each topic, refer to TopFeatures list code
# # Genes with the lowest empirical FDR
# pdf(file=paste0(FIGDIRTOP, "top.perturbation.list_empirical.fdr.pdf"), width=4, height=6)
# for ( topic in realPvals.df.guide.w$Topic %>% unique() ) {
#   toPlot <- realPvals.df.guide.w %>% subset(Topic == topic) %>% arrange(fdr) %>% slice(1:50) %>% mutate(Gene = gsub("Enhancer-at-CAD-SNP-","",Gene))
#   p <- toPlot %>% ggplot(aes(x=reorder(Gene, -fdr), y=fdr) ) + geom_col() + theme_minimal() +
#     coord_flip() + xlab("Top 50 Gene (Perturbation)") + ylab("Empirical FDR") + ggtitle(paste(SAMPLE[n], " perturbations, ", topic))
#   print(p)
# }
# dev.off()


# ##########################################################################
# ## most extreme log2FC (use omega)
# # # add ABC to gene.set.type.df for this particular plot
# # gene.set.type.df$type[which(gene.set.type.df$Gene %in% gene.set)] <- "ABC"
# pdf(file=paste0(FIGDIRTOP, "top.perturbation.list_log2FC.pdf"), width=4, height=6)
# for ( topic in colnames(gene.score) ) {
#   toPlot.all.test <- all.test.guide.w %>% subset(Topic==topic)
#   toPlot.fdr <- realPvals.df.guide.w %>% subset(Topic == topic) %>% select(Gene,fdr)
#   # assemble toPlot
#   toPlot <- gene.score %>% select(all_of(topic)) %>% mutate(Gene=rownames(.)) %>% merge(.,toPlot.all.test,by="Gene", all.x=T) %>%
#     merge(.,toPlot.fdr,by="Gene", all.x=T) %>%
#     merge(.,gene.set.type.df,by="Gene") %>%
#     mutate(Gene = gsub("Enhancer-at-CAD-SNP-","",Gene))
#   colnames(toPlot)[which(colnames(toPlot)==topic)] <- "log2FC"
#   toPlot <- toPlot %>% mutate(significant=ifelse((adjusted.p.value >= fdr.thr | is.na(fdr)), "", "*"))
#   toPlot.top <- toPlot %>% arrange(desc(log2FC)) %>% slice(1:25)
#   toPlot.bottom <- toPlot %>% arrange(log2FC) %>% slice(1:25)
#   toPlot.extreme <- rbind(toPlot.top, toPlot.bottom) %>%
#     mutate(color=ifelse(grepl("CAD",type), "red",
#                         ifelse(type=="non-expressed", "grey",
#                                ifelse(type=="other", "blue", "black")))) %>%
#       mutate(Gene = ifelse(type=="ABC", paste0("[ ", Gene, " ]"), Gene))
#   colors <- toPlot.extreme$color[order(toPlot.extreme %>% arrange(desc(log2FC)) %>% pull(color))]
#   p <- toPlot.extreme %>% arrange(desc(log2FC)) %>%  #mutate(Gene = paste0("<span style = 'color: ", color, ";'>", Gene, "</span>")) %>%
#     ggplot(aes(x=reorder(Gene, log2FC), y=log2FC, fill=significant)) + geom_col() + theme_minimal() +
#     coord_flip() + xlab("Most Extreme Gene (Perturbation)") + ylab("log2 Fold Change") + ggtitle(paste(SAMPLE[n], " perturbations, ", topic)) +
#     scale_fill_manual(values=c("grey", "#38b4f7")) +
#     geom_text(aes(label = significant)) +
#     theme(legend.position = "none", axis.text.y = element_text(colour = colors))
#   print(p)#here
# }
# dev.off()


# ##########################################################################
#     ## enrichment by GWAS classification
#     pdf(file=paste0(FIGDIRTOP,"_GWAS.class.enrichment.pdf"))
#     test.condition.names <- count.by.GWAS$test.type %>% unique()
#     for( i in 1:length(test.condition.names) ){
#         test.condition <- test.condition.names[i]
#         toPlot <- count.by.GWAS %>% subset(test.type==test.condition) %>% select(GWAS.class.enrichment,GWAS.classification,gene.count.per.GWAS.category) %>% unique()
#         p <- toPlot %>% ggplot(aes(x=GWAS.classification, y=GWAS.class.enrichment)) + geom_bar(stat="identity",width=0.5,fill = "#38b4f7") + mytheme +
#             ggtitle(paste0("Fraction of genes with significant topics,", test.condition)) + xlab("GWAS Classification") + ylab("Fraction of Perturbations") 
#         print(p)
#         p <- toPlot %>% ggplot(aes(x=GWAS.classification, y=gene.count.per.GWAS.category)) + geom_bar(stat="identity",width=0.5) + mytheme +
#             ggtitle(paste0("Number of genes with significant topics,", test.condition)) + xlab("GWAS Classification") + ylab("Number of Perturbations") 
#         print(p)
#     }
#     dev.off()


#     ## GWAS gene ranking plot
#     pdf(file=paste0(FIGDIRTOP,"_GWAS.gene.rank.barplot.pdf"), width=8, height=4)
#     test.condition.names <- count.by.GWAS$test.type %>% unique()
#     GWAS.type.list <- c("CAD","IBD")
#     p.list <- vector("list", length(test.condition.names) * length(GWAS.type.list))
#     for( i in 1:length(test.condition.names) ){
#         for (GWAS.type.index in 1:length(GWAS.type.list)) {
#             test.condition <- test.condition.names[i]
#             GWAS.type <- GWAS.type.list[GWAS.type.index]

#             ## select the columns and rows we need for this plot
#             toPlot.tmp <- count.by.GWAS %>% subset(test.type==test.condition & GWAS.classification==GWAS.type) %>% select(TSS.v.SNP.ranking.in.GWAS.category, passed.filter.ranking.count, GWAS.classification,total.TSS.v.SNP.ranking.count.per.GWAS.classification) %>% unique()
#             toPlot <- toPlot.tmp %>% mutate(ranking.fraction = passed.filter.ranking.count / total.TSS.v.SNP.ranking.count.per.GWAS.classification)
#             p <- toPlot %>% ggplot(aes(x=TSS.v.SNP.ranking.in.GWAS.category, y=ranking.fraction)) + geom_bar(stat="identity", width=0.5) + mytheme +
#                 ggtitle(paste0("Fraction of ", GWAS.type, " genes per distance ranking \n with significant topics, ", test.condition)) + xlab("Rank of Distance to the Closest SNP") + ylab("Fraction of Perturbations")
#             print(p)
 
#             ## concatenate all genes ranked as 5+
#             closest.ranking.boolean <- toPlot.tmp$TSS.v.SNP.ranking.in.GWAS.category < 5
#             toPlot.closer.genes <- toPlot.tmp[which(closest.ranking.boolean),] %>% as.data.frame()
#             toPlot.farther.genes <- toPlot.tmp[which(!closest.ranking.boolean),] %>% ungroup() %>% select(-GWAS.classification, -TSS.v.SNP.ranking.in.GWAS.category) %>% apply(2,sum) %>% t() %>% as.data.frame() %>% mutate(TSS.v.SNP.ranking.in.GWAS.category = "5+", .before = passed.filter.ranking.count) %>% mutate(GWAS.classification = GWAS.type, .after = passed.filter.ranking.count)
            
#             toPlot <- rbind(toPlot.closer.genes, toPlot.farther.genes) %>% mutate(ranking.fraction = passed.filter.ranking.count / total.TSS.v.SNP.ranking.count.per.GWAS.classification)
#             p <- toPlot %>% ggplot(aes(x=TSS.v.SNP.ranking.in.GWAS.category, y=ranking.fraction)) + geom_bar(stat="identity", width=0.5) + mytheme +
#                 ggtitle(paste0("Fraction of ", GWAS.type, " genes per distance ranking \n with significant topics, ", test.condition)) + xlab("Rank of Distance to the Closest SNP") + ylab("Fraction of Perturbations")
#             print(p)
 
#             p.list[[(i-1)*2 + GWAS.type.index]] <- p
#             toPlot <- toPlot %>% select(-ranking.fraction) %>% melt(id.vars=c("TSS.v.SNP.ranking.in.GWAS.category","GWAS.classification"), variable.name="ranking.type", value.name="ranking.count")
#             toPlot$GWAS.classification <- factor(toPlot$GWAS.classification)
#             toPlot$ranking.type <- factor(toPlot$ranking.type)
#             p <- toPlot %>% ggplot(aes(x=TSS.v.SNP.ranking.in.GWAS.category, y=ranking.count, fill=ranking.type)) + geom_bar(position="dodge",stat="identity",width=0.5) + mytheme +
#                 ggtitle(paste0("Number of ", GWAS.type, " genes with significant topics, ", test.condition)) + xlab("Rank of Distance to the Closest SNP") + ylab("Number of Perturbations") + scale_fill_manual(values=wes_palette(n=length(unique(toPlot$ranking.type)), name="Darjeeling2"), name = "Type", labels = c("Significant Genes", "All Genes"))
#             print(p)

#             ## plot 4
#             toPlot <- count.by.GWAS %>% subset(test.type==test.condition & GWAS.classification==GWAS.type) %>% select(TSS.v.SNP.ranking.in.GWAS.category, GWAS.classification, passed.filter.ranking.count, total.TSS.v.SNP.ranking.count.per.GWAS.classification) %>% unique() %>% melt(id.vars=c("TSS.v.SNP.ranking.in.GWAS.category","GWAS.classification"), variable.name="ranking.type", value.name="ranking.count")
#             p <- toPlot %>% ggplot(aes(x=TSS.v.SNP.ranking.in.GWAS.category, y=ranking.count, fill=ranking.type)) + geom_bar(position="dodge",stat="identity",width=0.5) + mytheme +
#                 ggtitle(paste0("Number of ", GWAS.type, " genes with significant topics, ", test.condition)) + xlab("Rank of Distance to the Closest SNP") + ylab("Number of Perturbations") + scale_fill_manual(values=wes_palette(n=length(unique(toPlot$ranking.type)), name="Darjeeling2"), name = "Type", labels = c("Significant Genes", "All Genes"))
#             print(p)
#         }
#     }
#     dev.off()

    
#     ## GWAS gene ranking by topic ## debug
#     pdf(file=paste0(FIGDIRTOP,"_GWAS.gene.rank.barplot.by.topic.pdf"), width=8, height=4) ##todo:210812:
#     test.condition.names <- count.by.GWAS.withTopic$test.type %>% unique()
#     for (GWAS.type.index in 1:length(GWAS.type.list)) {
#         GWAS.type <- GWAS.type.list[GWAS.type.index]
#         for( i in 1:length(test.condition.names) ){
#             test.condition <- test.condition.names[i]
#             toPlot.allTopics <- count.by.GWAS.withTopic %>% subset(test.type==test.condition & GWAS.classification == GWAS.type) %>%
#                 ungroup() %>%
#                 select(Topic,
#                        TSS.v.SNP.ranking.in.GWAS.category,
#                        passed.filter.ranking.count,
#                        GWAS.classification,
#                        total.TSS.v.SNP.ranking.count.per.GWAS.classification
#                        ) %>% unique()
#             toPlot.allTopics$GWAS.classification <- factor(toPlot.allTopics$GWAS.classification)
#             for( t in sort(unique(toPlot.allTopics$Topic %>% gsub("topic_","",.))) ) {
#                 toPlot.tmp <- toPlot.allTopics %>% subset(Topic == paste0("topic_",t)) 
#                 toPlot <- toPlot.tmp %>% mutate(ranking.fraction = passed.filter.ranking.count / total.TSS.v.SNP.ranking.count.per.GWAS.classification)
#                 p1 <- toPlot %>% ggplot(aes(x=TSS.v.SNP.ranking.in.GWAS.category, y=ranking.fraction)) + geom_bar(stat="identity", width=0.5) + mytheme +
#                     ggtitle(paste0("Fraction of ", GWAS.type, " significant genes per distance ranking \n in topic ", t, ", ",  test.condition)) + xlab("Rank of Distance to the Closest SNP") + ylab("Fraction of Perturbations")
#                 ## print(p1)

#                 ## concatenate all genes ranked as 5+
#                 closest.ranking.boolean <- toPlot.tmp$TSS.v.SNP.ranking.in.GWAS.category < 5
#                 toPlot.closer.genes <- toPlot.tmp[which(closest.ranking.boolean),] %>% as.data.frame() %>% select(-Topic)
#                 toPlot.farther.genes <- toPlot.tmp[which(!closest.ranking.boolean),] %>% ungroup() %>% select(-Topic, -GWAS.classification, -TSS.v.SNP.ranking.in.GWAS.category) %>% apply(2,sum) %>% t() %>% as.data.frame() %>% mutate(TSS.v.SNP.ranking.in.GWAS.category = "5+", .before = passed.filter.ranking.count) %>% mutate(GWAS.classification = GWAS.type, .after = passed.filter.ranking.count)
                
#                 toPlot <- rbind(toPlot.closer.genes, toPlot.farther.genes) %>% mutate(ranking.fraction = passed.filter.ranking.count / total.TSS.v.SNP.ranking.count.per.GWAS.classification)
#                 p2 <- toPlot %>% ggplot(aes(x=TSS.v.SNP.ranking.in.GWAS.category, y=ranking.fraction)) + geom_bar(stat="identity", width=0.5) + mytheme +
#                 ggtitle(paste0("Fraction of ", GWAS.type, " significant genes per distance ranking \n in topic ", t, ",  ", test.condition)) + xlab("Rank of Distance to the Closest SNP") + ylab("Fraction of Perturbations")
#             ## print(p2)
 
#             toPlot <- toPlot %>% select(-ranking.fraction) %>% melt(id.vars=c("TSS.v.SNP.ranking.in.GWAS.category","GWAS.classification"), variable.name="ranking.type", value.name="ranking.count")
#             toPlot$GWAS.classification <- factor(toPlot$GWAS.classification)
#             toPlot$ranking.type <- factor(toPlot$ranking.type)
#             p3 <- toPlot %>% ggplot(aes(x=TSS.v.SNP.ranking.in.GWAS.category, y=ranking.count, fill=ranking.type)) + geom_bar(position="dodge",stat="identity",width=0.5) + mytheme +
#                 ggtitle(paste0("Number of ", GWAS.type, " genes DE in topic ", t, ", ", test.condition)) + xlab("Rank of Distance to the Closest SNP") + ylab("Number of Perturbations") + scale_fill_manual(values=wes_palette(n=length(unique(toPlot$ranking.type)), name="Darjeeling2"), name = "Type", labels = c("Significant Genes", "All Genes"))
#             ## print(p3)

#             ## plot 4
#             toPlot <- count.by.GWAS.withTopic %>% ungroup() %>% subset(Topic == paste0("topic_",t) & test.type==test.condition & GWAS.classification==GWAS.type) %>% select(TSS.v.SNP.ranking.in.GWAS.category, GWAS.classification, passed.filter.ranking.count, total.TSS.v.SNP.ranking.count.per.GWAS.classification) %>% unique() %>% melt(id.vars=c("TSS.v.SNP.ranking.in.GWAS.category","GWAS.classification"), variable.name="ranking.type", value.name="ranking.count")
#             p4 <- toPlot %>% ggplot(aes(x=TSS.v.SNP.ranking.in.GWAS.category, y=ranking.count, fill=ranking.type)) + geom_bar(position="dodge",stat="identity",width=0.5) + mytheme +
#                 ggtitle(paste0("Number of ", GWAS.type, " significant genes in topic ", t, ", ", test.condition)) + xlab("Rank of Distance to the Closest SNP") + ylab("Number of Perturbations") + scale_fill_manual(values=wes_palette(n=length(unique(toPlot$ranking.type)), name="Darjeeling2"), name = "Type", labels = c("Significant Genes", "All Genes"))
#             ## print(p4)

#                 p <- ggarrange(p4 + ggtitle("") + xlab(""), p1 + ggtitle("") + xlab(""), p3 + ggtitle("") + xlab(""), p2 + ggtitle("") + xlab(""), nrow=2, ncol=2, common.legend=T)
#                 p <- annotate_figure(p, top = text_grob(paste0("Significant ", GWAS.type, " Genes in Topic ", t, ", ", test.condition), face="bold", size=14),
#                 bottom = "Rank of Distance to the Closest SNP")                
#                 print(p)            
#                 ## toPlot <- toPlot %>% subset(Topic==paste0("topic_",t) & test.type = test.condition)
#                 ## toPlot$GWAS.classification <- factor(toPlot$GWAS.classification)
#                 ## p <- toPlot %>% ggplot(aes(x=TSS.v.SNP.ranking, y=ranking.count, fill=GWAS.classification)) + geom_bar(position="dodge",stat="identity",width=0.5) + mytheme +
#                 ##     ggtitle(paste0("Number of genes with significant topic ", t, ", ", test.condition)) + xlab("Distance Ranking to the Closest SNP") + ylab("Number of Perturbations") + scale_fill_manual(values=wes_palette(n=length(unique(toPlot$GWAS.classification)), name="Darjeeling2")) +
#                 ##     theme(legend.title = element_blank())
#                 ## print(p)
#             }
#         }
#     }
#     dev.off()
    

#     ## ## GWAS ranking by distance in kb
#     ## pdf(file=paste0(FIGDIRTOP,"_GWAS.gene.distance.by.bp.cdf.pdf"), width=8, height=4)
#     ## test.condition.names <- count.by.GWAS$test.type %>% unique()
#     ## GWAS.type.list <- c("CAD","IBD")
#     ## toPlot.list <- vector("list", length(test.condition.names) * length(GWAS.type.list))
#     ## for( i in 1:length(test.condition.names) ){
#     ##     for (GWAS.type.index in 1:length(GWAS.type.list)) {
#     ##         test.condition <- test.condition.names[i]
#     ##         GWAS.type <- GWAS.type.list[GWAS.type.index]

#     ##         toPlot <- count.by.GWAS %>% subset(test.type==test.condition & GWAS.classification==GWAS.type) %>% select(TSS.dist.to.SNP, GWAS.classification)
#     ##         toPlot$GWAS.classification <- factor(toPlot$GWAS.classification)
#     ##         p <- toPlot %>% ggplot(aes(x=TSS.dist.to.SNP)) + stat_ecdf() + mytheme +
#     ##             ggtitle(paste0(GWAS.type, " genes with significant topics, ", test.condition)) + xlab("Distance to the Closest SNP (in bp)") + ylab("Fraction of Significant Perturbed Genes") 
#     ##         print(p)
#     ##         toPlot.list[[(i-1)*2 + GWAS.type.index]] <- toPlot %>% mutate(test.type=test.condition)
#     ##     }
#     ## }
#     ## dev.off()
    
#     ## Need to double check test.type
#     ## toPlot.all <- do.call(rbind, toPlot.list) %>% mutate(GWAS.classification = paste0(GWAS.classification, ".significant")) %>% rbind(ref.table %>% ungroup() %>% subset(GWAS.classification %in% GWAS.type.list) %>% select(TSS.dist.to.SNP, GWAS.classification) %>% mutate(test.type = "None"))
#     ## for( i in 1:length(test.condition.names) ){
#     ##     test.condition <- test.condition.names[i]
#     ##     toPlot <- toPlot.all %>% subset(test.type %in% c("None", test.condition))
#     ##     p <- toPlot %>% ggplot(aes(x=TSS.dist.to.SNP, color = GWAS.classification)) + stat_ecdf() + mytheme +
#     ##         xlim(0,1000000) +
#     ##         ggtitle(paste0("Test by ", test.condition)) + xlab("Distance to the Closest SNP (in bp)") + ylab("Fraction of Perturbed Genes")
#     ##     print(p)
# }
 


# ## heatmap of ptb correlation for factor values##here210810



# ## heatmap of ptb correlation for factor values threshold with BH adjusted.p.value < 0.1

# ## heatmap of factor correlation by expressed gene raw weights

# ## heatmap of factor correlation by expressed gene zscore


# ## UMAP based on factor log2FC values (to see how perturbations cluster)
# file.name <- paste0(OUTDIRSAMPLE, "gene.score.SeuratObject.RDS")
# if(file.exists(file.name)) {
#     s.gene.score <- readRDS(file.name)
# } else {
#     s.gene.score <- CreateSeuratObject(
#         counts = gene.score %>% t(), # cells as columns
#         project = SAMPLE
#     )
#     ## do this once -- these take a long time (and save these)
#     s.gene.score <- ScaleData(s.gene.score)
#     s.gene.score <- RunPCA(s.gene.score, verbose = FALSE, features = rownames(s.gene.score)) ##debug210810:multivariate_normal_UMAP?

#     ## s.gene.score <- SCTransform(s.gene.score)
#     ## s.gene.score <- RunPCA(s.gene.score, verbose = FALSE)
#     s.gene.score <- FindNeighbors(s.gene.score, dims = 1:10)
#     s.gene.score <- FindClusters(s.gene.score, resolution = 0.06)
#     s.gene.score <- RunUMAP(s.gene.score, dims = 1:10)
#     saveRDS(s.gene.score, file=file.name)
# } 

# ## actual UMAP
# pdf(file=paste0(FIGDIRTOP,"perturbation.UMAP.based.on.factor.log2FC.pdf"))
# DimPlot(s.gene.score, reduction = "umap", label=TRUE) %>% print()
# dev.off()




# # ##########################################################################
# # ## load topic model analysis file to get expressed genes
# # cNMF.result.file <- paste0(TMDIR,SAMPLE,"/cNMF_results.k_", k, ".dt_", DENSITY.THRESHOLD, ".RData")
# # print(cNMF.result.file)
# # load(cNMF.result.file)

# file.name <- paste0(OUTDIRSAMPLE,"/cNMFAnalysis.factorMotifEnrichment.",SUBSCRIPT.SHORT,".RData")
# print(file.name)
# if(file.exists((file.name))) { 
#     load(file.name)
#     print(paste0("loading ", file.name))
# }
# motif.enrichment.variables <- c("all.enhancer.fisher.df", "all.promoter.fisher.df", 
#     "promoter.wide", "enhancer.wide", "promoter.wide.binary", "enhancer.wide.binary",
#      "enhancer.wide.10en6", "enhancer.wide.binary.10en6", "all.enhancer.fisher.df.10en6",
#      "promoter.wide.10en6", "promoter.wide.binary.10en6", "all.promoter.fisher.df.10en6",
#      "all.promoter.ttest.df", "all.promoter.ttest.df.10en6", "all.enhancer.ttest.df", "all.enhancer.ttest.df.10en6")
# motif.enrichment.variables.missing <- (!(motif.enrichment.variables %in% ls())) %>% as.numeric %>% sum 
# if ( motif.enrichment.variables.missing > 0 | opt$recompute ) {
    
#     ## load hg38 promoter region file
#     promoter.region.hg38.original <- read.table(file=paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS300bp.hg38.bed"), header = F, stringsAsFactors = F)  %>% `colnames<-`(c("chr","start","end", "gene","cell.type","strand")) %>% mutate(sequence_name = paste0(chr, ":", start, "-", end, "(", strand, ")"))

#     ## keep only the expressed gene's motifs for background
#     expressed.genes <- theta %>% rownames()
#     promoter.region.hg38 <- promoter.region.hg38.original %>% mutate(expressed = (gene %in% expressed.genes) ) %>% filter(expressed)

#     ## background stats
#         print(opt$motif.promoter.background)
#         motif.background <- read.delim(file=paste0(ifelse(opt$motif.promoter.background!="", opt$motif.promoter.background, "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv")), header=T, stringsAsFactors=F) %>% filter(!grepl("#", motif_id))  # 30 seconds
#     colnames(motif.background)[colnames(motif.background) == "strand"] <- "motif.matched.strand"

#     print("start compiling Fisher Exact Test table")

#     ## get gene regions
#     all.gene.regions <- theta.zscore %>% as.data.frame() %>% mutate(gene = rownames(.)) %>% melt(id.vars="gene", variable.name="topic", value.name="zscore")
#     top.gene.regions <- all.gene.regions %>% group_by(topic) %>% arrange(desc(zscore)) %>% slice(1:100)


#     ##### This is working! #####
#     gene.motif.list.to.tbl <- function(gene.rank.df, motif.background) {
#         ## INPUT
#         ## gene.rank.df: table with columns "gene" and "topic".
#         ## motif.background: table from fimo output, with at least "motif_id" and "sequence_name".
#         ##
#         ## OUTPUT
#         ## list contains a count table and a binary count table with all genes on the rows and all features (motifs and topics) on the columns
        
#         # try to make it faster by first merging KL gene regions and all promoter motif matches, then spread
#         ## gene.rank.df <- top.gene.regions # to be passed into the function
#         long <- merge(gene.rank.df, motif.background, by.x="gene", by.y="sequence_name", all=T) # two minutes
#         ## tbl <- long %>% group_by(motif_id, topic) %>% summarize(count=n()) # 5 seconds ish
#         tbl2 <- long %>% group_by(gene, motif_id, topic) %>% summarize(count=n()) # three minutes
#         wide <- tbl2 %>% spread(key="motif_id", value="count", fill=0) %>% mutate(topic=paste0("topic_", topic), value=1) %>% spread(key="topic", value=value) %>% select(-topic_NA)# 40 seconds
#         wide.binary <- tbl2 %>% mutate(count = ifelse(count > 0, 1, 0)) %>% spread(key="motif_id", value="count", fill=0) %>% mutate(topic=paste0("topic_", topic), value=1) %>% spread(key="topic", value=value) %>% select(-topic_NA) # <40 seconds
#         wide.binary[is.na(wide.binary)] <- 0
#         return(list(table=wide, binary.table=wide.binary))
#     }

#     ## get topic gene ranking
#     promoter.wide.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.", SUBSCRIPT.SHORT, ".txt")
#     promoter.wide.binary.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.binary.", SUBSCRIPT.SHORT, ".txt")
#     if(file.exists(promoter.wide.path) & file.exists(promoter.wide.binary.path) & !opt$recompute) { 
#         promoter.wide <- read.delim(promoter.wide.path, stringsAsFactors=F) 
#         promoter.wide.binary <- read.delim(promoter.wide.binary.path, stringsAsFactors=F) 
#         } else {
#         wide.list <- gene.motif.list.to.tbl(top.gene.regions, motif.background) 
#         promoter.wide <- wide.list[[1]]
#         promoter.wide.binary <- wide.list[[2]]
#         write.table(promoter.wide, file=paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.", SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep="\t") # < 30 seconds
#         write.table(promoter.wide.binary, file=paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.binary.", SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep="\t")
#     }


#     ## Fisher Exact Test
#     fisher.test.on.motifs <- function(wide, wide.binary) {
#         motif.ary <- colnames(wide)[!grepl("gene|topic|<NA>", colnames(wide))]
#         all.fisher.df <- do.call(rbind, lapply(motif.ary, function(motif.here) {
#             do.call(rbind,lapply(1:k, function(t) {
#                 topic <- paste0("topic_",t)
#                 tbl.here <- table(wide.binary%>%pull(all_of(motif.here)), wide.binary %>% pull(all_of(topic)))
#                 fisher.result <- fisher.test(tbl.here)
#                 fisher.df <- fisher.result %>% broom::glance() %>% mutate(motif=motif.here, topic=topic, in.promoter.global=tbl.here[2,1], global.count = sum(tbl.here[,1]), in.promoter.topic=tbl.here[1,2], topic.count=sum(tbl.here[,2]))
#                 return(fisher.df)
#             })) %>% return()
#         })) # < 10 minutes
#         all.fisher.df <- all.fisher.df %>% mutate(p.adjust = p.adjust(p.value), in.promoter.global.fraction = in.promoter.global/global.count, in.promoter.topic.fraction = in.promoter.topic / topic.count, enrichment = in.promoter.topic.fraction / in.promoter.global.fraction, enrichment.log2fc = log2(enrichment))
#         return(all.fisher.df)
#     }

#     all.promoter.fisher.df.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.fisher.enrichment.", SUBSCRIPT.SHORT,".txt")
#     if(file.exists(all.promoter.fisher.df.path) & !opt$recompute) {
#         all.promoter.fisher.df <- read.delim(all.promoter.fisher.df.path, stringsAsFactors=F)
#     } else {
#         all.promoter.fisher.df <- fisher.test.on.motifs(promoter.wide, promoter.wide.binary)
#         write.table(all.promoter.fisher.df, file=all.promoter.fisher.df.path, quote=F,row.names=F,sep="\t")
#     }
    
#     promoter.wide.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.10en6.", SUBSCRIPT.SHORT, ".txt")
#     promoter.wide.binary.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.binary.10en6.", SUBSCRIPT.SHORT, ".txt")
#     if(file.exists(promoter.wide.10en6.path) & file.exists(promoter.wide.binary.10en6.path) & !opt$recompute){
#         promoter.wide.10en6 <- read.delim(promoter.wide.10en6.path, stringsAsFactors=F)
#         promoter.wide.binary.10en6 <- read.delim(promoter.wide.binary.10en6.path, stringsAsFactors=F)
#     } else {
#         wide.list <- gene.motif.list.to.tbl(top.gene.regions, motif.background %>% subset(`p.value` < 10^-6))  ##debug:210814:why is there a `X.NA.` column?
#         promoter.wide.10en6 <- wide.list[[1]]
#         promoter.wide.binary.10en6 <- wide.list[[2]]
#         write.table(promoter.wide.10en6, file=promoter.wide.10en6.path, row.names=F, quote=F, sep="\t") # < 30 seconds
#         write.table(promoter.wide.binary.10en6, file=promoter.wide.binary.10en6.path, row.names=F, quote=F, sep="\t")
#     }
    
#     ## test for motifs with p-value < 10^-6
#     all.promoter.fisher.df.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.fisher.enrichment_motif.thr.10en6_", SUBSCRIPT.SHORT,".txt")
#     if(file.exists(all.promoter.fisher.df.10en6.path) & !opt$recompute) {
#         all.promoter.fisher.df.10en6 <- read.delim(all.promoter.fisher.df.10en6.path, stringsAsFactors=F)
#     } else {
#         all.promoter.fisher.df.10en6 <- fisher.test.on.motifs(promoter.wide.10en6, promoter.wide.binary.10en6)
#         write.table(all.promoter.fisher.df.10en6, file=all.promoter.fisher.df.10en6.path, quote=F,row.names=F,sep="\t")
#     }
    
    
#     ## t-test on motif counts
#     ttest.on.motifs <- function(wide){##todo:210812:a
#         motif.ary <- colnames(wide)[!grepl("gene|topic|<NA>", colnames(wide))]
#         all.ttest.df <- do.call(rbind, lapply(motif.ary, function(motif.here) {
#             all.topic.ttest.df <- do.call(rbind,lapply(1:k, function(t) {
#                 topic <- paste0("topic_",t)
#                 wide.here <- wide %>% select(gene, all_of(topic), all_of(motif.here))
#                 wide.here[is.na(wide.here)] <- 0
#                 top.gene.boolean <- wide.here %>% pull(all_of(topic)) == 1
#                 background.count <- wide.here[which(!top.gene.boolean),] %>% pull(all_of(motif.here))
#                 top.gene.count <- wide.here[which(top.gene.boolean),] %>% pull(all_of(motif.here))
#                 test.result <- t.test(top.gene.count, background.count)$p.value
#                 output <- data.frame(motif=motif.here,
#                                      topic=topic,
#                                      p.value = test.result) %>%
#                     mutate(global.mean=background.count %>% mean,
#                            global.var=background.count %>% var,
#                            global.count=dim(wide)[1],
#                            top.gene.mean=top.gene.count %>% mean,
#                            top.gene.var=top.gene.count %>% var,
#                            topic.count=top.gene.count %>% length,
#                            enrichment=top.gene.mean/global.mean,
#                            enrichment.log2fc=log2(enrichment)
#                         )                
#             }))
#         })) %>% mutate(p.adjust = p.adjust(p.value, method="BH"))##here210813
#         ## adjust -Inf and +Inf for log2FC
#         inf.index <- all.ttest.df %>% pull(enrichment.log2fc) %>% is.infinite %>% which
#         neg.inf.index <- ( all.ttest.df %>% pull(top.gene.mean) == 0 ) %>% which
#         pos.inf.index <- ( all.ttest.df %>% pull(global.mean) == 0 ) %>% which
#         all.ttest.df$enrichment.log2fc[inf.index] <- 0
#         min.enrichment.log2fc.value <- min(all.ttest.df$enrichment.log2fc)
#         max.enrichment.log2fc.value <- max(all.ttest.df$enrichment.log2fc)
#         all.ttest.df$enrichment.log2fc[neg.inf.index] <- min.enrichment.log2fc.value - 3
#         all.ttest.df$enrichment.log2fc[pos.inf.index] <- max.enrichment.log2fc.value + 3
#         return(all.ttest.df)
#     }

#     print("Promoter motif count t-test")
#     all.promoter.ttest.df.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.count.ttest.enrichment.", SUBSCRIPT.SHORT,".txt")
#     if(file.exists(all.promoter.ttest.df.path) & !opt$recompute) {
#         all.promoter.ttest.df <- read.delim(all.promoter.ttest.df.path, stringsAsFactors=F)
#     } else {
#         all.promoter.ttest.df <- ttest.on.motifs(promoter.wide)
#         write.table(all.promoter.ttest.df, file=all.promoter.ttest.df.path, quote=F,row.names=F,sep="\t")
#     }
    
#     print("Promoter motif (p-value < 10^-6)  count t-test")
#     all.promoter.ttest.df.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.count.ttest.enrichment_motif.thr.10en6_", SUBSCRIPT.SHORT,".txt")
#     if(file.exists(all.promoter.ttest.df.10en6.path) & !opt$recompute) {
#         all.promoter.ttest.df.10en6 <- read.delim(all.promoter.ttest.df.10en6.path, stringsAsFactors=F)
#     } else {
#         all.promoter.ttest.df.10en6 <- ttest.on.motifs(promoter.wide.10en6)
#         write.table(all.promoter.ttest.df.10en6, file=all.promoter.ttest.df.10en6.path, quote=F,row.names=F,sep="\t")
#     }  

    
#     ##  motif enrichment 
#     ## load enhancer motif matches
#         print(opt$motif.enhancer.background)
#         motif.background <- read.delim(file=paste0(ifelse(opt$motif.enhancer.background!="", opt$motif.enhancer.background, "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv")), header=F, stringsAsFactors=F) %>%
#         `colnames<-`(c("motif_id", "motif_alt_id", "enhancer_region", "enhancer_type", "gene_region","sequence_name","start","stop","strand","score","p-value","q-value","matched_sequence")) %>% filter(!grepl("#|motif_id", motif_id))  # more than 30 seconds, minutes?
#     colnames(motif.background)[colnames(motif.background) == "strand"] <- "motif.matched.strand"
#         motif.background <- motif.background %>% filter(!grepl("promoter", enhancer_type))
        
#     print("start compiling enhancer Fisher Exact Test table")

#     enhancer.wide.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.table.", SUBSCRIPT.SHORT, ".txt")
#     enhancer.wide.binary.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.table.binary.", SUBSCRIPT.SHORT, ".txt")
#     if(file.exists(enhancer.wide.path) & file.exists(enhancer.wide.binary.path) & !opt$recompute) {
#         enhancer.wide <- read.delim(enhancer.wide.path, stringsAsFactors=F)
#         enhancer.wide.binary <- read.delim(enhancer.wide.binary.path, stringsAsFactors=F)
#     } else {
#         wide.list <- gene.motif.list.to.tbl(top.gene.regions, motif.background) 
#         enhancer.wide <- wide.list[[1]]
#         enhancer.wide.binary <- wide.list[[2]]
#         write.table(enhancer.wide, file=enhancer.wide.path, row.names=F, quote=F, sep="\t") # < 30 seconds
#         write.table(enhancer.wide.binary, file=enhancer.wide.binary.path, row.names=F, quote=F, sep="\t")
#     }



#     ## Fisher Exact Test for Enhancers
#     all.enhancer.fisher.df.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.fisher.enrichment.", SUBSCRIPT.SHORT,".txt")
#     if(file.exists(all.enhancer.fisher.df.path) & !opt$recompute) {
#         all.enhancer.fisher.df <- read.delim(all.enhancer.fisher.df.path, stringsAsFactors=F)
#     } else {
#         all.enhancer.fisher.df <- fisher.test.on.motifs(enhancer.wide, enhancer.wide.binary)
#         write.table(all.enhancer.fisher.df, file=all.enhancer.fisher.df.path, quote=F,row.names=F,sep="\t")
#     }

#     print("started compiling enhancer fisher exact test table for motif with p-value < 10^-6")
    
#     enhancer.wide.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.table.10en6.", SUBSCRIPT.SHORT, ".txt")
#     enhancer.wide.binary.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.table.binary.10en6.", SUBSCRIPT.SHORT, ".txt")
#     if(file.exists(enhancer.wide.10en6.path) & file.exists(enhancer.wide.binary.10en6.path) & !opt$recompute) {
#         enhancer.wide.10en6 <- read.delim(enhancer.wide.10en6.path, stringsAsFactors=F)
#         enhancer.wide.binary.10en6 <- read.delim(enhancer.wide.binary.10en6.path, stringsAsFactors=F)
#     } else {
#         wide.list <- gene.motif.list.to.tbl(top.gene.regions, motif.background %>% subset(`p-value` < 10^-6)) 
#         enhancer.wide.10en6 <- wide.list[[1]]
#         enhancer.wide.binary.10en6 <- wide.list[[2]]
#         write.table(enhancer.wide.10en6, file=enhancer.wide.10en6.path, row.names=F, quote=F, sep="\t") # < 30 seconds
#         write.table(enhancer.wide.binary.10en6, file=enhancer.wide.binary.10en6.path, row.names=F, quote=F, sep="\t")
#     }
    
#     ## test for motifs with p-value < 10^-6
#     all.enhancer.fisher.df.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.fisher.enrichment_motif.thr.10en6_", SUBSCRIPT.SHORT,".txt")
#     if(file.exists(all.enhancer.fisher.df.10en6.path) & !opt$recompute) {
#         all.enhancer.fisher.df.10en6 <- read.delim(all.enhancer.fisher.df.10en6.path, stringsAsFactors=F)
#     } else {
#         all.enhancer.fisher.df.10en6 <- fisher.test.on.motifs(enhancer.wide.10en6, enhancer.wide.binary.10en6)
#         write.table(all.enhancer.fisher.df.10en6, file=all.enhancer.fisher.df.10en6.path, quote=F,row.names=F,sep="\t")
#     }
    
#     print("Enhancer motif count t-test")
#     all.enhancer.ttest.df.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.count.ttest.enrichment.", SUBSCRIPT.SHORT,".txt")
#     if(file.exists(all.enhancer.ttest.df.path) & !opt$recompute) {
#         all.enhancer.ttest.df <- read.delim(all.enhancer.ttest.df.path, stringsAsFactors=F)
#     } else {
#         all.enhancer.ttest.df <- ttest.on.motifs(enhancer.wide)
#         write.table(all.enhancer.ttest.df, file=all.enhancer.ttest.df.path, quote=F,row.names=F,sep="\t")
#     }
    
#     print("Enhancer motif (p-value < 10^-6)  count t-test")
#     all.enhancer.ttest.df.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.count.ttest.enrichment_motif.thr.10en6_", SUBSCRIPT.SHORT,".txt")
#     if(file.exists(all.enhancer.ttest.df.10en6.path) & !opt$recompute) {
#         all.enhancer.ttest.df.10en6 <- read.delim(all.enhancer.ttest.df.10en6, stringsAsFactors=F)
#     } else {
#         all.enhancer.ttest.df.10en6 <- ttest.on.motifs(enhancer.wide.10en6)
#         write.table(all.enhancer.ttest.df.10en6, file=all.enhancer.ttest.df.10en6.path, quote=F,row.names=F,sep="\t")
#     }
    

#     save(all.enhancer.fisher.df, all.promoter.fisher.df, promoter.wide, enhancer.wide, promoter.wide.binary, enhancer.wide.binary,
#          enhancer.wide.10en6, enhancer.wide.binary.10en6, all.enhancer.fisher.df.10en6,
#          promoter.wide.10en6, promoter.wide.binary.10en6, all.promoter.fisher.df.10en6,
#          all.promoter.ttest.df, all.promoter.ttest.df.10en6,
#          all.enhancer.ttest.df, all.enhancer.ttest.df.10en6,
#      file=paste0(file.name))


    
# }



# ##########################################################################
# ## motif enrichment plot
# ep.names <- c("enhancer", "promoter")
# for (ep.type in ep.names) {
#     pdf(file=paste0(FIGDIRTOP,"zscore.",ep.type,".motif.enrichment.pdf"))
#     all.volcano.plots(get(paste0("all.",ep.type,".fisher.df")), ep.type, ranking.type="z-score")##here210812
#     dev.off()
#     pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment_motif.thr.10e-6.pdf"))
#     all.volcano.plots(get(paste0("all.",ep.type,".fisher.df.10en6")), ep.type, ranking.type="z-score")##here210812
#     dev.off()
#     pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest.pdf"))
#     all.volcano.plots(get(paste0("all.",ep.type,".ttest.df")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score")
#     dev.off() 
#     pdf(file=paste0(FIGDIRTOP, "zscore.",ep.type,".motif.enrichment.by.count.ttest_motif.thr.10e-6.pdf"))
#     all.volcano.plots(get(paste0("all.",ep.type,".ttest.df.10en6")) %>% subset(top.gene.mean != 0 & !grepl("X.NA.",motif)), ep.type, ranking.type="z-score")
#     dev.off() 
# }



# ##########################################################################
# ## GSEA

# ## check if files already exist
# check.file <- paste0(FGSEADIR,"/fgsea_all_pathways_df_", c("raw.score", "z.score"), "_", SUBSCRIPT, ".RData") 
# if(file.exists(check.file) %>% as.numeric() %>% sum() == length(check.file) & !opt$recompute ) {
#     for (i in 1:length(check.file)) {
#         print(paste0("Loading ", check.file[i]))
#         load(check.file[i])
#     }
# } else {

# ## import gene set annotation data
# msigdb.all <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data//msigdb.v7.2.entrez.gmt")) # downloaded from https://www.gsea-msigdb.org/gsea/downloads.jsp#msigdb
# msigdb.H <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data//h.all.v7.2.entrez.gmt"))
# msigdb.c2.cp <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c2.cp.v7.2.entrez.gmt"))
# msigdb.c2 <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c2.all.v7.2.entrez.gmt"))
# msigdb.c3 <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c3.all.v7.2.entrez.gmt"))
# msigdb.c5 <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c5.all.v7.2.entrez.gmt"))

# msigdb.names <- paste0("msigdb.", c("all", "H", "c2.cp", "c2", "c3", "c5"))

# ## get conversion list (gene name to entrez id)
# convertion.list <- as.list(org.Hs.egALIAS2EG) %>% lapply(., `[[`, 1)
# gene.to.genbank.id <- do.call(rbind, convertion.list) %>% as.data.frame() %>% `colnames<-`("genbank_id")
# gene.to.genbank.id$Gene <- rownames(gene.to.genbank.id)
# gene.to.genbank.id$genbank_id <- gene.to.genbank.id$genbank_id %>% as.character() %>% as.numeric()

# ## loop over all pathways
# for (type in c("raw.score", "z.score")){ # "raw.score", "z.score"
#     if (type == "z.score") input.data <- theta.zscore else input.data <- theta
#     input.data <- input.data %>% as.data.frame() %>% mutate(Gene = rownames(.)) %>% melt(id.vars="Gene", value.name="weight", variable.name="topic") ## convert to long table
#     fgsea.df.list <- vector("list", length(msigdb.names))
#     for (msigdb.index in 1:length(msigdb.names)) {
#         msigdb.pathway <- msigdb.names[msigdb.index]
#         ## specify data set to use
#         gsea.pathways=get(msigdb.pathway)
#         ## initialize output storage variables
#         ## out.gsea <- data.frame(ID=character(), Description=character(), setSize=double(), enrichmentScore=double(), NES=character(), pvalue=double(), p.adjust=double(), qvalues=double())
#         fgseaRes <- vector("list", k)
#         fgseaRes.df <- vector("list", k) # store fgseaRes that's converted to df format (remove leading Edge)
#         for(t in c(1:k)){
#             input.data.id <- input.data %>% subset(topic==t) %>% merge(., gene.to.genbank.id, by="Gene") %>% arrange(desc(weight))

#             ## fgsea
#             gene.list <- input.data.id$weight
#             names(gene.list) <- input.data.id$genbank_id
#             ## if (short.list) {
#             ##   names(gene.list)[(number+1):length(gene.list)] <- 0
#             ## }
#             fgseaRes[[t]] <- fgsea(pathways = gsea.pathways, 
#                                    stats    = gene.list,
#                                    nperm    = 10000,
#                                    ## eps      = 0.0,
#                                    minSize  = 3,
#                                    maxSize  = 800) %>% mutate(topic=t) %>% arrange(pval) 
#             fgseaRes[[t]]$leadingEdge <- sapply(fgseaRes[[t]]$leadingEdge, paste, collapse=",") # convert list in leadingEdge to string
#             fgseaRes.df[[t]] <- data.frame(pathway = fgseaRes[[t]]$pathway,
#                                            pval = fgseaRes[[t]]$pval,
#                                            padj = fgseaRes[[t]]$padj,
#                                            ES = fgseaRes[[t]]$ES,
#                                            NES = fgseaRes[[t]]$NES,
#                                            nMoreExtreme = fgseaRes[[t]]$nMoreExtreme,
#                                            numGenesInList = fgseaRes[[t]]$size,
#                                            topic = fgseaRes[[t]]$topic)
#         }
#         fgsea.df.list[[msigdb.index]] <- do.call(rbind, fgseaRes.df) %>% mutate(database = msigdb.pathway, padj.over.topics = p.adjust(pval))
#         ## save fgseaRes list
#         save(fgseaRes, file=paste0(FGSEADIR,"/fgsea_",msigdb.pathway,"_",type,"_", SUBSCRIPT, ".RData"))
#     }
#     fgsea.df <- do.call(rbind, fgsea.df.list)
#     save(fgsea.df, file=paste0(FGSEADIR,"/fgsea_all_pathways_df_",type, "_", SUBSCRIPT, ".RData"))
#     write.table(fgsea.df, file=paste0(FGSEADIR, "/fgsea_", type, "_", SUBSCRIPT, ".txt"), row.names=F, quote=F, sep="\t")
#     write.table(fgsea.df %>% subset(padj < 0.1), file=paste0(FGSEADIR, "/fgsea_", type, "_p.adj.0.1_", SUBSCRIPT, ".txt"), row.names=F, quote=F, sep="\t")
# }

# ## scratch for plotting GSEA results
# toplot <- fgsea.df %>% subset(database == "msigdb.all" & padj < 0.1) %>% select(topic, padj, pathway, ES, NES) %>% group_by(topic) %>% summarize(sig.pathway.count = n())

# } ## bracket for checking files


# pdf(paste0(FGSEAFIG, "/msigdb.all.sig.pathways.pdf"))
# p <- toplot %>% ggplot(aes(x=topic, y=sig.pathway.count)) + geom_bar(stat="identity", fill ="#38b4f7") + mytheme +
#     ggtitle(paste0(SAMPLE, " msigdb all, GSEA, number of significant pathway per topic")) + xlab("Topic") + ylab("Number of Significant Pathways")
# print(p)
# dev.off()


# ## Pairwise Pearson correlation for perturbation's topic expression
# d <- cor(gene.score %>% t(), method="pearson")
# m <- as.matrix(d)
# pdf(file=paste0(FIGDIRTOP,"cluster.ptb.by.Pearson.corr.pdf"), width=75, height=75)
# plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, K=", k, ", ", SAMPLE, " topic clustering by Euclidean distance"))
# dev.off()
# png(file=paste0(FIGDIRTOP, "cluster.ptb.by.Pearson.corr.png"), width=1500, height=1500)
# plotHeatmap(m, labCol=rownames(m), margins=c(12,12), title=paste0("cNMF, K=", k, ", ", SAMPLE, " topic clustering by Euclidean distance"))
# dev.off()
    
# ##########################################################################
# ##### SUMMARY PLOTS #####

# file.name=paste0(opt$raw.mtx.RDS.dir, "_modified.multiTargetGuide.cell.names.RDS")
# file.txt.name=paste0(opt$raw.mtx.RDS.dir, "_modified.multiTargetGuide.cell.names.txt")
# fc.file.name=paste0(opt$raw.mtx.RDS.dir, "_FC_modified.multiTargetGuide.cell.names.RDS")
# log2fc.file.name=paste0(opt$raw.mtx.RDS.dir, "_log2FC_modified.multiTargetGuide.cell.names.RDS")
#     if(file.exists(file.name) & file.exists(fc.file.name) & file.exists(log2fc.file.name)) {
#         print(paste0("Loading ", file.name))
#         X.full = readRDS(file.name)
#         print(paste0("Loading ", fc.file.name))
#         fc.X.full = readRDS(fc.file.name)
#         print(paste0("Loading ", log2fc.file.name))
#         log2fc.X.full = readRDS(log2fc.file.name)##here210813
#     } else {
        
#         if(!("X" %in% ls())) { # if X, omega, and theta are not loaded, load them!
#             raw.input.path <- opt$raw.mtx.dir
#             raw.RDS.path <- opt$raw.mtx.RDS.dir
#             if(file.exists(raw.input.path)) {
#                 X <- read.delim(raw.input.path)
#             } else if (file.exists(raw.RDS.path)) {
#                 X <- readRDS(raw.RDS.path)
#             } else {
#                 message("Full matrix is not provided")
#             }
            
#         }
#         X.full <- X %>% as.matrix() # %>% adjust.multiTargetGuide.rownames()
        
#         ## duplicate cells with guide that targets multiple genes (e.g. "A-and-B", and change the cell names to "A", "B", and "A-and-B" accordingly ## need 512 GB memory
#         rownames(X.full) <- rownames(X.full) %>% gsub("CSNK2B-and-CSNK2B", "CSNK2B",.) ## added 210930
#         X.full <- adjust.multiTargetGuide.rownames(X.full) # this might take a lot of memory
#         ## fold change compared to control cells
#         ctrl.vec <- X.full %>% subset(grepl("negative|safe",rownames(.))) %>% as.matrix %>% apply(2, mean)
#         fc.X.full <- sweep(X.full, 2, ctrl.vec, "/")
#         log2fc.X.full <- log2(fc.X.full)
#         saveRDS(X.full, file=paste0(opt$raw.mtx.RDS.dir, "_modified.multiTargetGuide.cell.names.RDS"))
#         write.table(X.full, file=file.txt.name, sep="\t", quote=F)
#         saveRDS(fc.X.full, file=paste0(opt$raw.mtx.RDS.dir, "_FC_modified.multiTargetGuide.cell.names.RDS"))
#         saveRDS(log2fc.X.full, file=paste0(opt$raw.mtx.RDS.dir, "_log2FC_modified.multiTargetGuide.cell.names.RDS"))
#     }


# ######################################################################
# ## Seurat UMAPs
# file.name <- paste0(opt$datadir,"/", SAMPLE, "_SeuratObjectUMAP.rds")
# subset.file.name <- paste0(opt$datadir,"/", SAMPLE, "_subset_SeuratObjectUMAP.rds")
# options(future.globals.maxSize=1000*1024^2)
# if(!file.exists(file.name)){
#     s <- CreateSeuratObject(
#       counts = X.full %>% t(), # cells as columns
#       project = SAMPLE
#     )
#     # do this once -- these take a long time (and save these)
#     s <- SCTransform(s)
#     s <- RunPCA(s, verbose = FALSE)
#     s <- FindNeighbors(s, dims = 1:10)
#     s <- FindClusters(s, resolution = 0.06)
#     s <- RunUMAP(s, dims = 1:10)
#     saveRDS(s, file=file.name)
#     s.subset <- subset(x = s, downsample = 8000)
#     saveRDS(s.subset, file=subset.file.name)
# } else {
#     s <- readRDS(file.name)
# }
# s@meta.data <- s@meta.data[,-which(grepl("topic",colnames(s@meta.data)))]
# s <- AddMetaData(s, metadata = omega, col.name = paste0("K",k,"_",colnames(omega)))
# ## s.subset@meta.data <- s.subset@meta.data[,-which(grepl("topic",colnames(s.subset@meta.data)))]
# ## subset.index <- which(colnames(s.subset) %in% rownames(omega))
# ## subset.omega <- omega[subset.index, ]
# ## s.subset <- AddMetaData(s.subset, metadata = subset.omega, col.name = paste0("K",k,"_",colnames(subset.omega)))


# pdf(paste0(FIGDIRTOP, "Factor.Expression.UMAP.pdf"))
# DimPlot(s, reduction = "umap", label=TRUE) %>% print()
# meta.data.names <- colnames(s[[]])
# plot.features <- paste0("K",k,"_",colnames(omega))
# for (feature.name in plot.features) {
#     feature.vec <- s@meta.data %>% select(all_of(feature.name))
#     if(feature.vec[1,1] %>% is.numeric()) { # numeric
#         FeaturePlot(s, reduction = "umap", features=feature.name) %>% print()
#     } else { # discrete values / categories
#         Idents(s) <- s@meta.data %>% select(feature.name)
#         DimPlot(s, reduction = "umap", label=TRUE) %>% print()
#     }
# }
# dev.off()

# ## ## 10x lane sample label UMAP ## checking for batch effects
# ## CBC.sample <- colnames(s) %>% gsub("^.*-", "", .) %>% gsub("scRNAseq_2kG_","",.)
# ## CBC.sample.short <- CBC.sample %>% gsub("_.*$","",.)
# ## s <- AddMetaData(s, CBC.sample, col.name="sample.label.10X")
# ## s <- AddMetaData(s, CBC.sample.short, col.name="sample.label.short")
# ## pdf(file=paste0(FIGDIRTOP, "10X.lane.sample.label.UMAP.pdf"), width=10, height=6)
# ## Idents(s) <- s$sample.label.10X
# ## DimPlot(s, reduction = "umap", label=F, group.by="sample.label.10X") %>% print()
# ## Idents(s) <- s$sample.label.short
# ## DimPlot(s, reduction = "umap", label=F) %>% print()
# ## DimPlot(s, reduction = "umap", label=F, split.by="sample.label.short") %>% print()
# ## dev.off()






# ##HERE

# # need ~36G (16G doens't work)
# ## ann.X.full <- merge(barcode.names, X.full, by.x="CBC", by.y=0, all.y=TRUE)  # annotated omega
# ## ann.X.full <- cbind(X.full, barcode.names) # this step requires too much memory
# ## ann.X.full.filtered <- ann.X.full %>% group_by(Gene,Guide) %>% mutate(nCount=n()) %>% subset(nCount > cell.count.thr) %>% select(-nCount) %>% ungroup() %>%
# ## group_by(Gene) %>% mutate(nCount = length(unique(Guide))) %>% subset(nCount > guide.count.thr) %>% select(-nCount) %>% as.data.frame()
#     ## ann.X.full <- X.full %>% as.data.frame() %>% mutate(long.CBC = row.names(.)) %>% separate(., col=long.CBC, into=c("CBC", "Guide", "CBC"), remove=F)  # still takes too much memory
#     tokeep.index <- rownames(X.full) %in% ann.omega.filtered$long.CBC ## need to expand X.full to match ann.omega.filtered due to guides that target two promoters
#     ann.X.full.filtered <- X.full[tokeep.index,]
#     # adjust colnames, remove ENSG number
#     tmp <- colnames(ann.X.full.filtered) %>% strsplit(., split=":") %>% sapply("[[",1)
#     tmpp <- data.frame(table(tmp)) %>% subset(Freq > 1)  # keep row names that have duplicated gene names but different ENSG names
#     tmp.copy <- tmp
#     tmp.copy[grepl(paste0(tmpp$tmp,collapse="|"),tmp)] <- colnames(ann.X.full.filtered)[grepl(paste0(tmpp$tmp,collapse="|"), colnames(ann.X.full.filtered))]
#     colnames(ann.X.full.filtered) <- tmp.copy # the above section takes a while
#     ## get ctrl log2 transformed expression
#     tokeep.index <- grepl("control|targeting",rownames(ann.X.full.filtered))
#     ## ctrl.X <- ann.X.full.filtered %>% subset(grepl("control|targeting",Gene))
#     ctrl.X <- ann.X.full.filtered[tokeep.index,]
#     ## get ctrl topic weights
#     ctrl.ann.omega <- ann.omega.filtered %>% subset(grepl("control|targeting",Gene)) %>% `rownames<-`(.$long.CBC)
# X.gene.names <- rownames(ann.X.full.filtered) %>% strsplit(., split=":") %>% sapply("[[",1) %>% gsub("_multiTarget|-TSS2","",.)


# ##########################################################################
# ## functions
# get.average.ptb.gene.expression.based.on.ctrl <- function(ptb, ptb.df, ctrl.df, mode="per.guide") {
#   # ptb: column that has gene expression or topic weight, (e.g. "SWAP70" or "topic_4")
#   if (mode=="per.guide") { # first average by guide then average by perturbation
#       e.name <- ptb.df %>% rownames() %>% gsub("_multiTarget|-TSS2","",.) %>% strsplit(., split=":") %>% sapply("[[", 1) %>% unique() 
#       ptb.gene.column.index <- which(grepl(paste0("^",ptb,"$"), colnames(ptb.df)))
#       expression.of.ptb <- ptb.df[,ptb.gene.column.index] %>% as.data.frame() %>% `rownames<-`(rownames(ptb.df)) %>% `colnames<-`(ptb) %>% mutate(long.CBC=rownames(.)) %>% separate(., col=long.CBC, into=c("Gene", "Guide", "CBC"), remove=F, sep=":") %>% group_by(Guide) %>% summarise(Gene.expression=mean(get(ptb))) %>% mutate(Gene=e.name)
#     ## expression.of.ptb <- ptb.df %>% select(c(all_of(ptb), colnames(guideCounts))) %>%
#     ##   group_by(Guide) %>% summarise(Gene.expression=mean(get(ptb))) %>% mutate(Gene=e.name)
      
#       expression.ctrl.ptb <- ctrl.df[,ptb.gene.column.index] %>% as.data.frame() %>% `rownames<-`(rownames(ctrl.df)) %>% `colnames<-`(ptb) %>% mutate(long.CBC=rownames(.)) %>% separate(., col=long.CBC, into=c("Gene", "Guide", "CBC"), remove=F, sep=":") %>% group_by(Guide) %>% summarise(Gene.expression=mean(get(ptb))) %>% mutate(Gene="Control")# %>% select(c(all_of(ptb), colnames(guideCounts))) %>% group_by(Guide) %>% summarise(Gene.expression=mean(get(ptb))) %>% mutate(Gene="Control")
#       toCalculate <- rbind(expression.of.ptb, expression.ctrl.ptb)
#       toCalculate$Gene.expression <- toCalculate$Gene.expression / (toCalculate %>% subset(Gene=="Control") %>% pull(Gene.expression) %>% mean()) * 100
#     toCalculate <- toCalculate %>% group_by(Gene) %>% mutate(Gene = paste0(gsub("topic_","Topic ", Gene), "\n(n=", n(), ")"))
#     toPlot <- toCalculate %>% group_by(Gene) %>% summarise(mean.expression=mean(Gene.expression), error.bar=1.96 * sd(Gene.expression)/sqrt(n()), count=n()) # 1.96 * sd(vals)/sqrt(length(vals))
#   } else { # directly average by cell
#       e.name <- ptb.df %>% rownames() %>% gsub("_multiTarget|-TSS2","",.) %>% strsplit(., split=":") %>% sapply("[[", 1) %>% unique() 
#     # e.name <- ptb.df$Gene %>% unique()
#     ptb.gene.column.index <- which(grepl(paste0("^",ptb,"$"), colnames(ptb.df)))
#     expression.of.ptb <- ptb.df[,ptb.gene.column.index] %>% as.data.frame() %>% `rownames<-`(rownames(ptb.df)) %>% `colnames<-`(ptb) %>% mutate(long.CBC=rownames(.)) %>% separate(., col=long.CBC, into=c("Gene", "Guide", "CBC"), remove=F, sep=":") %>% mutate(Gene=e.name)
#     expression.ctrl.ptb <- ctrl.df[,ptb.gene.column.index] %>% as.data.frame() %>% `rownames<-`(rownames(ctrl.df)) %>% `colnames<-`(ptb) %>% mutate(long.CBC=rownames(.)) %>% separate(., col=long.CBC, into=c("Gene", "Guide", "CBC"), remove=F, sep=":") %>% mutate(Gene="Control")
#     ## expression.of.ptb <- ptb.df %>% select(c(all_of(ptb), colnames(guideCounts))) %>% mutate(Gene=e.name)
#     ## expression.ctrl.ptb <- ctrl.df %>% select(c(all_of(ptb), colnames(guideCounts))) %>% mutate(Gene="Control")
#     toCalculate <- rbind(expression.of.ptb, expression.ctrl.ptb) %>% select(all_of(ptb),Guide,Gene) %>% mutate(Gene.expression = get(ptb))
#     toCalculate$Gene.expression <- toCalculate$Gene.expression / (toCalculate %>% subset(Gene=="Control") %>% pull(Gene.expression) %>% mean()) * 100
#     toCalculate <- toCalculate %>% group_by(Gene) %>% mutate(Gene = paste0(gsub("topic_","Topic ", Gene), "\n(n=", n(), ")"))
#     toPlot <- toCalculate %>% group_by(Gene) %>% summarise(mean.expression=mean(Gene.expression), error.bar=1.96 * sd(Gene.expression)/sqrt(n()), count=n())
#     # 1.96 * sd(vals)/sqrt(length(vals))
#   }
#   return(list(toCalculate=toCalculate, toPlot=toPlot))
# }

# lm_eqn <- function(df){
#   m <- lm(y ~ x, df);
#   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
#                    list(a = format(unname(coef(m)[1]), digits = 2),
#                         b = format(unname(coef(m)[2]), digits = 2),
#                         r2 = format(summary(m)$r.squared, digits = 3)))
#   as.character(as.expression(eq));
# }
# lm_eqn_manual <- function(a, b){
#   eq <- substitute(italic(y) == p + q %.% italic(x),
#                    list(p = format(unname(a), digits = 2),
#                         q = format(unname(b), digits = 2)))
#   as.character(as.expression(eq)) %>% return()
# }


# normalize.by.ctrl.avg <- function(ptb, ptb.df, ctrl.df, mode="per.guide"){
#   e.name <- ptb.df$Gene %>% gsub("_multiTarget|-TSS2","",.) %>% strsplit(., split=":") %>% sapply("[[", 1) %>% unique()
#   ptb.col.index <- which(grepl(paste0("^",ptb,"$"), colnames(ctrl.df)))
#   ctrl.df.tmp <- ctrl.df[,ptb.col.index] %>% as.data.frame() %>% `rownames<-`(rownames(ctrl.df)) %>% `colnames<-`(ptb)
#   ctrl.df.tmp <- ctrl.df.tmp %>% mutate(long.CBC=rownames(.)) %>% separate(., col=long.CBC, into=c("Gene", "Guide", "CBC"), sep=":", remove=F)

#   if (mode=="per.guide") { # first average by guide then average by perturbation
#       expression.ctrl.ptb.tmp <- ctrl.df.tmp %>% group_by(Guide) %>% summarise(Gene.expression=mean(get(ptb))) %>% mutate(Gene="Control")
#       ## expression.ctrl.ptb.tmp <- ctrl.df %>% select(c(all_of(ptb), colnames(guideCounts))) %>%
#       ## group_by(Guide) %>% summarise(Gene.expression=mean(get(ptb))) %>% mutate(Gene="Control")
#     # ptb.df with mean as % control and errorbar
#       expression.of.ptb <- ptb.df %>%
#       mutate(Gene.expression = get(ptb) / (expression.ctrl.ptb.tmp$Gene.expression %>% mean()) * 100) %>%
#       group_by(Guide) %>%
#       summarise(Gene.error.bar = 1.96 * sd(Gene.expression)/sqrt(n()),
#                 Gene.expression=mean(Gene.expression)) %>% mutate(Gene=e.name) # summarise or mutate? # need to keep CBC
#     # control with error bar
#       expression.ctrl.ptb <- ctrl.df.tmp %>%
#       mutate(Gene.expression = get(ptb) / (expression.ctrl.ptb.tmp$Gene.expression %>% mean()) * 100) %>%
#       group_by(Guide) %>%
#       summarise(Gene.error.bar = 1.96 * sd(Gene.expression)/sqrt(n()),
#                 Gene.expression=mean(Gene.expression)) %>% mutate(Gene="Control")
#       toCalculate <- rbind(expression.of.ptb, expression.ctrl.ptb)

#   } else { # directly average by cell
#     expression.ctrl.ptb <- ctrl.df.tmp %>% mutate(Gene="Control")
#     expression.of.ptb <- ptb.df %>% mutate(Gene=e.name)
#     toCalculate <- rbind(expression.of.ptb %>% select(all_of(ptb),Guide,Gene), expression.ctrl.ptb %>% select(all_of(ptb),Guide,Gene)) %>% mutate(Gene.expression = get(ptb))
#   }
#   toCalculate <- toCalculate %>% group_by(Gene) %>% mutate(Gene.count = paste0(gsub("topic_","Topic ", Gene), " (n=", n(), ")"))
#   return(toCalculate)
# }



# giant.summary.plot <- function(ptb, ptb.expressed.name, expressed.gene, mode.selection, enhancer=F) { 
#   # toPlot.per.guide.DE.test.wilcox <- per.guide.DE.test.wilcox %>% subset(Topic==ptb)
#   # if (enhancer) {
#   #   array <- ann.X.full.filtered$Gene[grep("E_at_", ann.X.full.filtered$Gene)] %>% unique()
#   #   e.name <- array[grep(ptb,array)]
#   # } else {
#   #   e.name <- ptb
#   # }
  
  
#   # tmp <- ptb # "E_at_" "-no"
#   # ptb <- expressed.gene # "GOSR2"
#   # e.name <- tmp

#     all.test.guide.w <- all.test %>% subset(test.type==paste0(mode.selection,".wilcoxon"))
#     realPvals.df.guide.w <- realPvals.df %>% subset(test.type==paste0(mode.selection,"per.guide.wilcoxon"))

    
#   #TODO: add rep1rep2 to ctrl.X %>% subset()
#   # for SEP=T, subset ctrl.X to the right sample's control
#   if(SEP) {
#     label.here <- strsplit(ptb, split="-") %>% unlist() %>% nth(2) %>% paste0("-",.)
#     ctrl.X.here <- ctrl.X %>% subset(grepl(label.here,Gene))
#     ctrl.ann.omega.here <- ctrl.ann.omega %>% subset(grepl(label.here,Gene))
#   } else  {
#     ctrl.X.here <- ctrl.X
#     ctrl.ann.omega.here <- ctrl.ann.omega
#   }

#     subset.index <- which(X.gene.names == ptb) 
#     ## ann.X.full.filtered.df <- 
#     toPlot.list <- get.average.ptb.gene.expression.based.on.ctrl(expressed.gene, ann.X.full.filtered[subset.index,], ctrl.X.here, mode=mode.selection)
#     toPlot1 <- toPlot.list[["toPlot"]]
#     order.toPlot.Gene <- function(df) {
#         df$Gene <- factor(x=df$Gene, levels=df$Gene[c(which(grepl("^Control", df$Gene)) %>% min(), which(!grepl("^Control", df$Gene))%>% min())])
#         return(df)
#     }
#     order.toPlot <- function(df, column) {
#         array <- df %>% pull(all_of(column))
#         df <- df %>% mutate(!!column:=factor(x=array, levels=array[c(which(grepl("^Control", array)) %>% min(), which(!grepl("^Control", array))%>% min())]))
#         return(df)
#     }
#     toPlot1 <- order.toPlot(toPlot1, column="Gene")

#   KD.per.guide <- toPlot.list[["toCalculate"]]
#   colnames(KD.per.guide)[colnames(KD.per.guide)=="Gene.expression"] <- "mean.expression"

#   # plot 5 data
#     subset.index <- which(X.gene.names == ptb)
#     ## ptb.gene.col.index <- which(grepl(ptb.expressed.name, colnames(ann.X.full.filtered)))
#     ptb.gene.col.index <- which(colnames(ann.X.full.filtered) == ptb.expressed.name)
#     ann.X.ptb <- ann.X.full.filtered[subset.index,ptb.gene.col.index] %>% as.data.frame() %>% `colnames<-`(ptb.expressed.name) %>% mutate(long.CBC = rownames(.)) %>% separate(., col=long.CBC, into=c("Gene", "Guide", "CBC"), sep=":", remove=F)# select ptb expression column
#   ## ann.X.ptb <- ann.X.full.filtered %>% subset(Gene==ptb) %>% select(colnames(guideCounts), all_of(expressed.gene)) # select ptb expression column
#   normalized.ann.X.ptb.by.guide <- normalize.by.ctrl.avg(expressed.gene, ann.X.ptb, ctrl.X.here, mode.selection) ##210823:debug expressed.gene or ptb
#   normalized.ann.X.ptb.by.guide[is.na(normalized.ann.X.ptb.by.guide)] <- 0
#   # topic
#   ptb.omega.filtered <- ann.omega.filtered %>% subset(Gene == ptb)

 
#   if (dim(ptb.omega.filtered)[1] > 0){  # if the guide didn't cause a fitness effect
#     if(!dir.exists(paste0(FIGDIRSAMPLE,"/gene.by.topic/"))) dir.create(paste0(FIGDIRSAMPLE,"/gene.by.topic/"), recursive=T)
#     if(!enhancer) pdf(file=paste0(FIGDIRSAMPLE,"/gene.by.topic/",ptb, ".", mode.selection, ".dt_", DENSITY.THRESHOLD, ".pdf"), width=16, height=8)
#     # Do gRNAs at the promoter reduce gene expression?
#     # old plot
#     p1 <- toPlot1 %>% ggplot(aes(x=Gene,y=mean.expression)) + geom_bar(stat='identity', width=0.5, fill="#38b4f7") +
#       geom_errorbar(data = toPlot1, aes(x=Gene, ymin=mean.expression-error.bar, ymax=mean.expression+error.bar), width=.15) +
#       # ylim(min(KD.per.guide$mean.expression - 25, 0), max(KD.per.guide$mean.expression + 50, 150)) +
#       ylab(paste0(expressed.gene, " RNA Expression\n(% vs control)")) + xlab("Guides") + mytheme +
#       # scale_y_continuous(limits = c(0, max(KD.per.guide$mean.expression + 50, 150)), breaks = round(seq(0, max(KD.per.guide$mean.expression + 50, 150), 20),20) ) +
#       geom_hline(yintercept = 100, color="grey", linetype="dashed", size=0.5) +
#       geom_jitter(data=KD.per.guide, size=0.25, width=0.15, color="red") +
#               ggdist::stat_halfeye(
#                     data=KD.per.guide,
#                     ## custom bandwidth
#                     adjust = .5, 
#                     ## adjust height
#                     width = .3, 
#                     ## move geom to the right
#                     justification = -1,
#                     ## remove slab interval
#                     .width = 0, 
#                     point_colour = NA,
#                     ## change violin color 
#                     fill = "#a0dafa"
#                 ) 

 
#     ## ## referenced https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-with-ggplot2/
#     ## p1 <- KD.per.guide %>% ggplot(aes(x=Gene, y=mean.expression)) +
#     ##     ## add half-violin from {ggdist} package 
#     ##     ggdist::stat_halfeye(
#     ##                 ## custom bandwidth
#     ##                 adjust = .5, 
#     ##                 ## adjust height
#     ##                 width = .3, 
#     ##                 ## move geom to the right
#     ##                 justification = -.4, 
#     ##                 ## remove slab interval
#     ##                 .width = 0, 
#     ##                 point_colour = NA
#     ##             ) +
#     ##     geom_boxplot(
#     ##         width = .1, 
#     ##         ## remove outliers
#     ##         outlier.color = NA ## `outlier.shape = NA` works as well
#     ##     ) +
#     ##     ## add justified jitter from the {gghalves} package
#     ##     gghalves::geom_half_point(
#     ##                   ## control point size
#     ##                   size = 0.5,
#     ##                   ## draw jitter on the left
#     ##                   side = "l", 
#     ##                   ## control range of jitter
#     ##                   range_scale = .4, 
#     ##                   ## add some transparency
#     ##                   alpha = .3
#     ##               ) +
#     ##     coord_cartesian(xlim = c(1.2, NA), clip = "off") +
#     ##     mytheme +
#     ##     geom_hline(
#     ##         yintercept = 100,
#     ##         linetype = "dashed",
#     ##         color = "#38b4f7",
#     ##         size = 0.5
#     ##     ) +
#     ##     ylab(paste0(expressed.gene, " RNA Expression\n(% vs control)")) +
#     ##     xlab("Perturbation") 


    
#     for (t in 1:dim(omega)[2]) { #:dim(omega)[2]
#       topic <- paste0("topic_", t)
#       toPlot.all.test <- all.test.guide.w %>% subset(Topic==topic)
#       toPlot.fdr <- realPvals.df.guide.w %>% subset(Topic == topic) %>% select(Gene,fdr)

#       # Which topics does the gene regulate?
#       toPlot.list <- get.average.ptb.gene.expression.based.on.ctrl(topic, ptb.omega.filtered %>% `rownames<-`(ptb.omega.filtered$long.CBC), ctrl.ann.omega.here %>% `rownames<-`(ctrl.ann.omega.here$long.CBC), mode=mode.selection)
#       toPlot2 <- toPlot.list[["toPlot"]]
#       # levels(toPlot2$Gene)[grep(topic,levels(toPlot2$Gene))] <- gsub("_", " ", toPlot2$Gene[grep(topic,levels(toPlot2$Gene))]) # change topic_x to "topic x"
#       toPlot2 <- order.toPlot.Gene(toPlot2)

#       toCalculate <- toPlot.list[["toCalculate"]]
#       colnames(toCalculate)[colnames(toCalculate)=="Gene.expression"] <- "mean.expression"

#       p2.y.step <- (max(toCalculate$mean.expression/5) %/% 20) * 20
#       p2.y.step <- ifelse(p2.y.step==0, 20, p2.y.step)
#       p2.y.step <- ifelse(p2.y.step > 100, 100, p2.y.step)

#       ## ## old plot
#       ## p2 <- toPlot2 %>% ggplot(aes(x=Gene,y=mean.expression)) + geom_bar(stat='identity', width=0.5, fill="#38b4f7") +
#       ##   geom_errorbar(aes(ymin=mean.expression-error.bar, ymax=mean.expression+error.bar), width=.15) +
#       ##   # ylim(0, max(toPlot2$mean.expression + toPlot2$error.bar + 50, 150)) +
#       ##   ylab(paste0(gsub("topic_","Topic ", topic), " Expression\n(% vs control)")) + xlab("Guides") + mytheme +
#       ##   # scale_y_continuous(limits = c(0, max(toCalculate$mean.expression + 50, 150)), breaks = seq(0, max(toCalculate$mean.expression + 50, 150), p2.y.step)) +
#       ##   geom_hline(yintercept = 100, color="grey", linetype="dashed", size=0.5) +
#       ##   geom_jitter(data=toCalculate, size=0.25, width=0.15, color="red")
 
#        p2 <- toPlot2 %>% ggplot(aes(x=Gene,y=mean.expression)) + geom_bar(stat='identity', width=0.5, fill="#38b4f7") +
#           ## geom_jitter(data=toCalculate, aes(x=Gene, y=mean.expression), size=0.25, width=0.15, color="red") +
#           geom_errorbar(aes(ymin=mean.expression-error.bar, ymax=mean.expression+error.bar), width=.15) +
#         # ylim(0, max(toPlot2$mean.expression + toPlot2$error.bar + 50, 150)) +
#         ylab(paste0(gsub("topic_","Topic ", topic), " Expression\n(% vs control)")) + xlab("Guides") + mytheme +
#         # scale_y_continuous(limits = c(0, max(toCalculate$mean.expression + 50, 150)), breaks = seq(0, max(toCalculate$mean.expression + 50, 150), p2.y.step)) +
#         geom_hline(yintercept = 100, color="grey", linetype="dashed", size=0.5) +
#         geom_jitter(data=toCalculate, size=0.25, width=0.15, color="red") +
#         ggdist::stat_halfeye(
#                     data=toCalculate,
#                     ## custom bandwidth
#                     adjust = .5, 
#                     ## adjust height
#                     width = .3, 
#                     ## move geom to the right
#                     justification = -1,
#                     ## remove slab interval
#                     .width = 0, 
#                     point_colour = NA,
#                     ## change violin color 
#                     fill = "#a0dafa"
#                 ) 

#       ## ## referenced https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-with-ggplot2/
#       ## p2 <- toCalculate %>% ggplot(aes(x=Gene, y = mean.expression)) +
#       ##     ggdist::stat_halfeye(
#       ##                 ## custom bandwidth
#       ##                 adjust = .5, 
#       ##                 ## adjust height
#       ##                 width = .3, 
#       ##                 ## move geom to the right
#       ##                 justification = -.4, 
#       ##                 ## remove slab interval
#       ##                 .width = 0, 
#       ##                 point_colour = NA
#       ##             ) +
#       ##     geom_boxplot(
#       ##         width = .1, 
#       ##         ## remove outliers
#       ##         outlier.color = NA ## `outlier.shape = NA` works as well
#       ##     ) +
#       ##     ## add justified jitter from the {gghalves} package
#       ##     gghalves::geom_half_point(
#       ##                   ## control point size
#       ##                   size = 0.5,
#       ##                   ## draw jitter on the left
#       ##                   side = "l", 
#       ##                   ## control range of jitter
#       ##                   range_scale = .4, 
#       ##                   ## add some transparency
#       ##                   alpha = .3
#       ##               ) +
#       ##     coord_cartesian(xlim = c(1.2, NA), clip = "off") +
#       ##     mytheme +
#       ##     geom_hline(
#       ##         yintercept = 100,
#       ##         linetype = "dashed",
#       ##         color = "#38b4f7",
#       ##         size = 0.5
#       ##     ) +
#       ##     ylab(paste0(gsub("topic_","Topic ", topic), " Expression\n(% vs control)")) +
#       ##     xlab("Perturbation")
          




      
#       # ## alternative to p2, but basically the same thing
#       # ## stat_summary removes some values before averaging and causes control to be not at 100%
#       # p3 <- toPlot.list[["toCalculate"]] %>% ggplot(aes(x=Gene,y=Gene.expression)) +
#       #   geom_bar(stat = "summary", fun = "mean", width=0.5, fill="#38b4f7") +
#       #   stat_summary(fun.data = mean_cl_normal,
#       #                geom = "errorbar", width=0.15) +
#       #   # ylim(0, max(toPlot2$mean.expression + toPlot2$error.bar + 50, 150)) +
#       #   ylab(paste0(gsub("_"," ", topic), " Expression\n(% vs control)")) + xlab("") + mytheme +
#       #   scale_y_continuous(limits = c(0, max(toPlot2$mean.expression + toPlot2$error.bar + 50, 150)), breaks = seq(0, max(toPlot2$mean.expression + toPlot2$error.bar + 50, 150), p2.y.step)) +
#       #   geom_hline(yintercept = 100, color="grey", linetype="dashed") +
#       #   geom_jitter(size=0.25, width=0.15, color="red")

#       ## #TODO: fix KL specificity score
#       ## # What cellular program does this topic represent?
#       ## toPlot <- data.frame(Gene=topFeatures %>% subset(topic == t) %>% pull(genes),
#       ##                      Score=topFeatures %>% subset(topic == t) %>% pull(scores)) %>%
#       ##   merge(., gene.def.pathways, by="Gene", all.x=T)
#       ## toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
#       ## p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score*100, fill=Pathway) ) + geom_col(width=0.5) + theme_minimal()
#       ## p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("KL score (gene specific to this topic)") +
#       ##   mytheme + theme(legend.position="bottom", legend.direction="vertical")
      
#       # raw weight version
#       ## hand annotated files by Helen
#       toPlot <- data.frame(Gene=topFeatures.raw.weight %>% subset(topic == t) %>% pull(Gene),
#                            Score=topFeatures.raw.weight %>% subset(topic == t) %>% pull(scores)) %>%
#         merge(., gene.def.pathways, by="Gene", all.x=T) %>% arrange(desc(Score)) %>% slice(1:10)
#       ## 210804 use Gavin's new summary annotations
#       ## toPlot <- data.frame(Gene=topFeatures.raw.weight %>% subset(topic == t) %>% pull(Gene),
#       ##                      Score=topFeatures.raw.weight %>% subset(topic == t) %>% pull(scores)) %>%
#       ##   merge(., summaries %>% select(Symbol, top_class), by.x="Gene", by.y="Symbol", all.x=T) %>% unique %>% `colnames<-`(c("Gene", "Score", "Pathway")) %>% arrange(desc(Score)) %>% slice(1:10)
#       toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
#       toPlot$Pathway[toPlot$Pathway == "unclassified"] <- "Other/Unclassified"
#       p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score*100, fill=Pathway) ) + geom_col(width=0.5) + theme_minimal()
#       p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("Specificity Score (z-score)") +
#         mytheme + theme(legend.position="bottom", legend.direction="vertical")

      

#       # plot 4
#       # add ABC to gene.set.type.df for this particular plot
#       ## gene.set.type.df$type[which(gene.set.type.df$Gene %in% gene.set)] <- "ABC"
#       # assemble toPlot
#         toPlot <- gene.score %>% select(all_of(topic)) %>% mutate(Gene=rownames(.)) %>%
#             merge(.,toPlot.all.test,by="Gene", all.x=T) %>%
#             merge(.,toPlot.fdr,by="Gene", all.x=T) %>%
#             ## merge(.,gene.set.type.df,by="Gene") %>%
#             mutate(Gene = gsub("Enhancer-at-CAD-SNP-","",Gene)) %>%
#             merge(., ref.table %>% select("Symbol", "TSS.dist.to.SNP", "GWAS.classification"), by.x="Gene", by.y="Symbol", all.x=T) %>%
#             mutate(EC_ctrl_text = ifelse(.$GWAS.classification == "EC_ctrls", "(+)", "")) %>%
#             mutate(GWAS.class.text = ifelse(grepl("CAD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb"),
#                                             ifelse(grepl("IBD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb_IBD"), ""))) %>%
#             mutate(ann.Gene = paste0(Gene, GWAS.class.text, EC_ctrl_text))
#         colnames(toPlot)[which(colnames(toPlot)==topic)] <- "log2FC"
#         toPlot <- toPlot %>% mutate(significant=ifelse((adjusted.p.value >= fdr.thr | is.na(adjusted.p.value)), "", "*")) %>%
#             arrange(log2FC) %>%
#             mutate(x = seq(nrow(.)))
#         label <- toPlot %>% subset(x <= 3 | x > (nrow(toPlot)-3) | Gene == ptb | adjusted.p.value < fdr.thr)
#       ## toPlot <- gene.score %>% select(all_of(topic)) %>% mutate(Gene=rownames(.)) %>% merge(.,toPlot.all.test,by="Gene", all.x=T) %>%
#       ##   merge(.,toPlot.fdr,by="Gene", all.x=T) %>%
#       ##   merge(.,gene.set.type.df,by="Gene") %>%
#       ##   mutate(Gene = gsub("Enhancer-at-CAD-SNP-","",Gene))
#       ## colnames(toPlot)[which(colnames(toPlot)==topic)] <- "log2FC"
#       ## toPlot <- toPlot %>% mutate(significant=ifelse((fdr >= fdr.thr | is.na(fdr)), "", "*")) %>%
#       ##   arrange(log2FC) %>%
#       ##   mutate(x = seq(nrow(.)))
#       ## label <- toPlot %>% subset(fdr < fdr.thr | x <= 3 | x > (nrow(toPlot)-3) | Gene == ptb | adjusted.p.value < fdr.thr)
#       mytheme <- theme_classic() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, face = "bold"))
#       p5 <- toPlot %>% ggplot(aes(x=reorder(Gene, log2FC), y=log2FC, color=significant)) + geom_point(size=0.75) + mytheme +
#         theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + scale_color_manual(values = c("#E0E0E0", "#38b4f7")) +
#         xlab("Perturbed Genes") + ylab(paste0("Topic ", t, " Expression log2 Fold Change")) +
#         geom_text_repel(data=label, box.padding = 0.5,
#                         aes(label=ann.Gene), size=5,
#                         color="black") +
#         theme(legend.position = "none") # legend at the bottom?


#       # plot 5: topic expression vs KD efficacy
#       # data
#       ptb.omega.t <- ptb.omega.filtered %>% select(all_of(c("Gene","Gene.full.name","long.CBC","CBC","Guide")), all_of(topic))
#       normalized.ptb.omega.t.by.guide <- normalize.by.ctrl.avg(topic, ptb.omega.t, ctrl.ann.omega.here, mode=mode.selection) %>%
#         select(-Gene, -Gene.count) %>%  # remove redundant columns before merging
#         `colnames<-`(gsub("Gene", "Topic", colnames(.)))
#       merged.ptb.normalized.X.omega.t.by.guide <- merge(normalized.ann.X.ptb.by.guide, normalized.ptb.omega.t.by.guide, by="Guide") # GOSR2-plus should have 10 guides, but only 8 left
#       ##debug: rbind(deparse.level, ...):   numbers of columns of arguments do not match
      
#       toPlot <- merged.ptb.normalized.X.omega.t.by.guide %>% order.toPlot(., column="Gene.count")
#       toPlot.ptb <- toPlot %>% subset(Gene == ptb) # for linear regression
#       fit <- tryCatch(york(toPlot.ptb %>% select(Gene.expression, Gene.error.bar, Topic.expression, Topic.error.bar)),
#                       error=function(cond){
#                         return(data.frame(a = 0, b = 0))
#                       })
#       p6a <- toPlot %>% ggplot(aes(x=Gene.expression, y=Topic.expression, color = Gene.count)) +
#         geom_point(size=1) +
#         ylab(paste0("Topic ", t, " Expression\n(% vs control)")) + xlab(paste0(expressed.gene, " RNA Expression\n(% vs control)")) + mytheme +
#         scale_color_manual(values=c("grey", "red"), name="Guide")  + theme(legend.position="bottom") +
#         geom_smooth(data = toPlot.ptb, method="lm", se=F, fullrange=T, size=0.5) +
#         annotate("text", size = 5, x = min(toPlot$Gene.expression ) + 30, y = min(toPlot$Topic.expression - 15), hjust=0.2,
#                  label = lm_eqn(data.frame(x=toPlot.ptb$Gene.expression,
#                                            y=toPlot.ptb$Topic.expression)), parse = TRUE) # , parse = TRUE

#       p6 <- toPlot %>% ggplot(aes(x=Gene.expression, y=Topic.expression, color = Gene.count)) +
#         geom_point(size=0.5) +
#         ylab(paste0("Topic ", t, " Expression\n(% vs control)")) + xlab(paste0(expressed.gene, " RNA Expression\n(% vs control)")) + mytheme +
#         scale_color_manual(values=c("grey", "red"), name="Guide")  + theme(legend.position="bottom") +
#         geom_abline(intercept = fit$a[1], slope = fit$b[1], col="red", size=0.5) +
#         annotate("text", size = 5, x = min(toPlot$Gene.expression ) + 30, y = min(toPlot$Topic.expression - 15), hjust=0.2,
#                  label = lm_eqn_manual(fit$a[1], fit$b[1]), parse = TRUE)


#       # set xlim and ylim to fit annotation!

#       ## TF motif enrichment volcano plot
#       toplot <- all.promoter.ttest.df %>% subset(topic==paste0("topic_",t) & top.gene.mean != 0)
#       volcano.plot <- function(toplot, EP.string) {
#           label <- toplot %>% subset(-log10(p.adjust) > 1 & enrichment.log2fc > 0) %>% mutate(motif.toshow = gsub("HUMAN.H11MO.", "", motif))
#           t <- gsub("topic_", "", toplot$topic[1])
#           p <- toplot %>% ggplot(aes(x=enrichment.log2fc, y=-log10(p.adjust))) + geom_point(size=0.5) + mytheme +
#               ggtitle(paste0(EP.string, " Motif Enrichment")) + xlab("Motif Enrichment (log2FC)") + ylab("-log10(adjusted p-value)") +
#               geom_text_repel(data=label, box.padding = 0.3,
#                               aes(label=motif.toshow), size=6,
#                               color="black")
#           return(p)
#       }
#       p.promoter.motif <- volcano.plot(toplot, "Promoter")
#       toplot <- all.enhancer.ttest.df.10en6 %>% subset(topic==paste0("topic_",t) & top.gene.mean != 0)
#       p.enhancer.motif <- volcano.plot(toplot, "Enhancer")


#       ## edgeR p-value vs ptb effect RNA expression for top 100 specific gene in the topic
#       ##use.edgeR.results
#       edgeR.expr.gene.names <- rownames(log2fc.edgeR) %>% strsplit(split=":") %>% sapply("[[",1)
#       ptb.colindex <- which(grepl(ptb, colnames(log2fc.edgeR)))
#       theta.zscore.t <- theta.zscore[,t] %>% as.data.frame %>% `colnames<-`("topic.zscore.weight") %>% mutate(genes=rownames(.)) %>% arrange(desc(topic.zscore.weight)) %>% mutate(gene.rank = 1:n(), top100=ifelse(gene.rank <= 100, paste0("Top 100 in Topic ", t), paste0("Not in Top 100")))
#       ptb.pval.log2fc <- merge(p.value.edgeR[,c(1,ptb.colindex)] %>% `colnames<-`(c("genes", "pval")) %>% mutate(genes = strsplit(genes, split=":") %>% sapply("[[",1)), log2fc.edgeR[,c(1,ptb.colindex)] %>% `colnames<-`(c("genes", "log2fc")) %>% mutate(genes = strsplit(genes, split=":") %>% sapply("[[",1)), by="genes") %>% merge(theta.zscore.t, by="genes") %>% mutate(nlog10pval = -log10(pval))
#       ptb.pval.log2fc$nlog10pval[ptb.pval.log2fc$nlog10pval > 15] <- 15
      
#       ptb.pval.log2fc$top100 <- factor(ptb.pval.log2fc$top100) %>% ordered(levels=c(paste0("Top 100 in Topic ", t), "Not in Top 100"))
      
#       label <- ptb.pval.log2fc %>% subset(top100!="Not in Top 100")
#       ptb.10X <- ptb.10X.name.conversion$`Name used by CellRanger`[which(grepl(ptb, ptb.10X.name.conversion$Symbol))]
#       label.self <- ptb.pval.log2fc %>% subset(genes == ptb.10X)
#       p.density <- ptb.pval.log2fc %>% ggplot(aes(x=log2fc, fill=top100)) + geom_density(alpha=0.4) + mytheme + xlab("RNA Expression\n(log2FC vs Control") + scale_fill_manual(values=c("red", "gray"), name = "") + theme(legend.position = "none") # + guides(fill=guide_legend(nrow=2, byrow=T))
#       p.edgeR <- ggplot(ptb.pval.log2fc, aes(x=log2fc, y=nlog10pval)) + geom_point(size=0.1, color = "gray") +
#           geom_point(data = label, size=0.1, color = "red") +
#           mytheme +
#           geom_text_repel(data=label.self, box.padding = 0.5,
#                           max.overlaps=30,
#                           aes(label=genes), size=4,
#                           color="blue") +
#           geom_text_repel(data=label, box.padding = 0.5,
#                           max.overlaps=30,
#                           aes(label=genes), size=4,
#                           color="black") +
#           scale_color_manual(values=c("gray", "red")) +
#           geom_vline(xintercept=0, col = "#38b4f7", lty=3) +
#           theme(legend.position="bottom") + 
#           ggtitle(paste0("Topic ", t, " Perturbation ", ptb)) +
#           xlab("Average RNA Expression (log2FC vs control)") + ylab("p-value (-log10)") # +
#       ## inset_element(p.density, left = 0.1, bottom = 0.75, right=0.9, top = 0.95)
#       p.edgeR.combined <- cowplot::plot_grid(p.density, p.edgeR, aign="v", ncol=1, rel_heights=c(0.15, 0.85))

#       ## ## old
#       ## p.top.left <- ggarrange(p1, p2, p6, nrow=1, widths=c(1.5,1.5,2))
#       ## ## p.bottom.left <- ggarrange(p6a, p6, nrow=1)
#       ## p.bottom.left <- ggarrange(p.promoter.motif, p.enhancer.motif, nrow=1)
#       ## p.left <- ggarrange(p.top.left, p.bottom.left, nrow=2)
#       ## p <- ggarrange(p.left, p4, p5, ncol=3, nrow=1, widths=c(2,1,1))
#       ## # p <- ggarrange(p1, p2, p4, p5, p6, ncol=4, nrow=1, widths=c(1,1,1.5,1.5))

#       p.top.left <- ggarrange(p1, p2, p6, nrow=1, widths=c(1.5,1.5,2))
#       p.bottom.left <- ggarrange(p.promoter.motif, p.enhancer.motif, nrow=1)
#       p.left <- ggarrange(p.top.left, p.bottom.left, nrow=2)
#       p.mid <- ggarrange(p4, p.edgeR, nrow=2, heights=c(1.5,1))
#       p <- ggarrange(p.left, p.mid, p5, ncol=3, nrow=1, widths=c(2,1,1))
      
#       print(p)
#     }
#     if(!enhancer) dev.off()
#   }
# }



# ##########################################################################
# ## slide 26 for all perturbation x topic pairs
# # function for plot 5

# ## # plot 3 data
# ## toSave.features <- read.delim(paste0(OUTDIRSAMPLE, "/topic.KL.score_K", k,ifelse(SEP, ".sep", ""), ".txt"),header=T, stringsAsFactors=F) ### commented out because it's not necessary in this pipeline
# ## topFeatures <- toSave.features %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:10)
# # plot 4 data
# all.test <- read.delim(file=paste0(OUTDIRSAMPLE, "/all.test.", SUBSCRIPT, ".txt"), header=T, stringsAsFactors=F) ##here
# realPvals.df <- read.delim(file=paste0(OUTDIRSAMPLE, "/all.expressed.genes.pval.fdr.",SUBSCRIPT,".txt"), header=T, stringsAsFactors=F)
# all.test.guide.w <- all.test %>% subset(test.type==opt$test.type)
# realPvals.df.guide.w <- realPvals.df %>% subset(test.type==opt$test.type)
# fdr.thr <- 0.1
# topFeatures.raw.weight <- theta.zscore %>% as.data.frame() %>% mutate(Gene=rownames(.)) %>% melt(id.vars="Gene", variable.name="topic", value.name="scores") %>% group_by(topic) %>% arrange(desc(scores)) %>% slice(1:50)


# ## topic plots
# ## load data
# file.name <- ifelse(SEP,
#                     paste0(OUTDIRSAMPLE,"/cNMFAnalysis.factorMotifEnrichment.",SUBSCRIPT.SHORT,".sep.RData"),
#                     paste0(OUTDIRSAMPLE,"/cNMFAnalysis.factorMotifEnrichment.",SUBSCRIPT.SHORT,".RData"))
# print(file.name)
# if(file.exists((file.name)) & !("all.promoter.fisher.df" %in% ls())) { 
#     load(file.name)
#     print(paste0("loading ", file.name))
# }
# type <- "z.score"
# fgsea.df <- read.delim(file=paste0(FGSEADIR, "/fgsea_", type, "_", SUBSCRIPT, ".txt"), header=T, stringsAsFactors=F) ##HERE
# fgsea.df.GO <- fgsea.df %>% subset(database == "msigdb.c3")
# fgsea.df.all <- fgsea.df %>% subset(database == "msigdb.all")


# # ptb.array <- c("GOSR2", "PRKCE", "PHACTR1", "EIF2B2")
#     ptb.array <- c("GOSR2","TP53","MAT2A","SKI","EDN1","SMAD3","PHACTR1","EDN1","GGCX","CDKN1A","EGFL7","ELOF1","GPANK1","YLPM1","TLNRD1","ITGA5","SKIV2L","LST1","R3HCC1L","UPF2","MEAF6", "CCM2", "CCM1", "ITGB1BP","HEG1")
#     ## ptb.array <- c("TP53", "SWAP70", "GOSR2", "PRKCE", "PHACTR1", "EIF2B2", "PPIF", "DMRTA1", "ADAMTS7", "VEZT", "MEAF6", "CDH5")
# # ptb.array <- guideCounts$Gene %>% unique()
# # ptb.array <- enhancer.set
# ## ptb.array <- append(ptb.array, all.test.guide.w$Gene %>% unique()) %>% append(CAD.focus.gene.set) %>% append(gene.set.type.df %>% subset(grepl("CAD_Loci", type))  %>% pull(Gene) %>% sort())
# ptb.array <- append(ptb.array, all.test.guide.w$Gene %>% unique())
# # ptb.array <- CAD.focus.gene.set[grepl("E_at_", CAD.focus.gene.set)]

# mode.selection="per.guide"
# for (mode.selection in c("per.guide", "per.cell")){

# ## new!
# for (ptb in ptb.array) { # make a separate one for enhancers, clear up cases like E_at_AAGAB,IQCH,LOC102723493&SMAD3
#   print(paste0(ptb, "\n\n"))
#   if(grepl("E_at_", ptb)) {
#     target.gene <- strsplit(gsub("E_at_","",ptb %>% strsplit(split="-") %>% unlist() %>% nth(1)), split=",|&") %>% unlist() %>% as.character()
#     expressed.gene <- 
#     if((target.gene %in% colnames(ann.X.full.filtered)) %>% as.numeric() %>% sum() > 0 &  # target gene must express
#        (ptb %in% ann.X.full.filtered$Gene) %>% as.numeric() %>% sum() > 0) { # enhancer must have enough cells and guides
#       pdf(file=paste0(FIGDIRSAMPLE,"/gene.by.topic/",ptb, ".", mode.selection, ".pdf"), width=15, height=8)
#       for(expressed.gene in target.gene) {
#         if(expressed.gene %in% colnames(ann.X.full.filtered)){
#           print(paste0("Enhancer of ", expressed.gene, "\n\n"))
#           giant.summary.plot(ptb, expressed.gene, mode.selection, enhancer=T)
#         }
#       }
#       dev.off()
#     }
#   } else{
#       expressed.gene <- ptb %>% strsplit(split="-") %>% unlist() %>% nth(1)
#       gene.index <- which(grepl(ptb, ptb.10X.name.conversion$Symbol))
#       if(length(gene.index) == 1) {
#           ptb.expressed.name <- ptb.10X.name.conversion$`Name used by CellRanger`[gene.index]
#           expressed.gene <- ptb.expressed.name
#       } else {
#           ptb.expressed.name <- expressed.gene
#       }
#       message(paste0("Expressed gene: ", expressed.gene, ", Perturbation: ", ptb, ", with expressed name: ", ptb.expressed.name))
#       ##   message(paste0("Expressed gene: ", expressed.gene, ", Perturbation: ", ptb))
#       if((colnames(X.full) == ptb.expressed.name) %>% as.numeric() %>% sum() > 0) giant.summary.plot(ptb, ptb.expressed.name, expressed.gene, mode.selection, enhancer=F)
#       ##if((colnames(ann.X.full.filtered) == expressed.gene) %>% as.numeric() %>% sum() > 0) giant.summary.plot(ptb, ptb.expressed.name, mode.selection, enhancer=F)
#   }
# }
# #HERE
# } # if(opt$subsample.type!="ctrl")

 

# ## ##scratch:210818
# ## ## topic 29 expression log2FC versus CDH5 RNA expression log2FC ##todo:210812
# ## for (t in 1:k) {
# ## ## t <- 29 ##scratch
# ## topic <- paste0("topic_",t)
# ##     ## ptb.name <- "CDH5" ##scratch top perturbations
# ##     top.ptb.name <- gene.score %>% as.data.frame %>% mutate(Gene = rownames(.)) %>% select(all_of(topic), Gene) %>% arrange(desc(get(topic))) %>% slice(1:10) %>% pull(Gene)
# ##     for( ptb.name in top.ptb.name) {
# ##         topic.log2fc.here <- log2fc.omega %>% select(all_of(topic)) %>% mutate(Gene=rownames(.)) %>% subset(grepl(ptb.name, Gene))  ## get Topic expression for Gene G
# ##         topic.fc.here <- fc.omega %>% select(all_of(topic)) %>% mutate(Gene=rownames(.)) %>% subset(grepl(ptb.name, Gene))

# ## ## inf.index <- topic.log2fc.here %>% pull(all_of(topic)) %>% is.infinite %>% which  ## get cell index for zero Topic expression entries
# ## ## neg.inf.index <- (topic.fc.here %>% pull(all_of(topic)) == 0) %>% which
# ## ## pos.inf.index <- which(!(inf.index %in% neg.inf.index))
# ## ## topic.log2fc.here[[topic]][inf.index] <- 0
# ## ## min.log2fc.here <- topic.log2fc.here %>% pull(all_of(topic)) %>% min
# ## ## max.log2fc.here <- topic.log2fc.here %>% pull(all_of(topic)) %>% max
# ## ## topic.log2fc.here[[topic]][neg.inf.index] <- floor(min.log2fc.here - 3)
# ## ## topic.log2fc.here[[topic]][pos.inf.index] <- ceil(max.log2fc.here + 3)

# ## remove.inf <- function(topic.log2fc.here, topic.fc.here, topic) {
# ##     topic.index <- which(grepl(topic, colnames(topic.log2fc.here)))
# ##     inf.index <- topic.log2fc.here[,topic.index] %>% is.infinite %>% which
# ##     neg.inf.index <- (topic.fc.here[,topic.index] == 0) %>% which

# ##     ## inf.index <- topic.log2fc.here %>% pull(all_of(topic)) %>% is.infinite %>% which  ## get cell index for zero Topic expression entries
# ##     ## neg.inf.index <- (topic.fc.here %>% pull(all_of(topic)) == 0) %>% which

# ##     pos.inf.index <- which(!(inf.index %in% neg.inf.index))
# ##     topic.log2fc.here[inf.index,topic.index] <- 0
# ##     min.log2fc.here <- topic.log2fc.here[,topic.index] %>% min
# ##     max.log2fc.here <- topic.log2fc.here[,topic.index] %>% max
# ##     ## min.log2fc.here <- topic.log2fc.here %>% pull(all_of(topic)) %>% min
# ##     ## max.log2fc.here <- topic.log2fc.here %>% pull(all_of(topic)) %>% max
# ##     topic.log2fc.here[neg.inf.index,topic.index] <- floor(min.log2fc.here - 3)
# ##     topic.log2fc.here[pos.inf.index,topic.index] <- ceil(max.log2fc.here + 3)
# ##     return(topic.log2fc.here)
# ## }

# ## topic.log2fc.here <- remove.inf(topic.log2fc.here, topic.fc.here, topic)


# ## ## select RNA expresion
# ## ## todo: function to extract df
# ## X.full.here <- X.full %>% subset(grepl(ptb.name, rownames(.))) %>% `colnames<-`(gsub(":.*$", "", colnames(.)))
# ## fc.X.full.here <- fc.X.full %>% subset(grepl(ptb.name, rownames(.))) %>% `colnames<-`(gsub(":.*$", "", colnames(.)))
# ## log2fc.X.full.here <- log2fc.X.full %>% subset(grepl(ptb.name, rownames(.))) %>% `colnames<-`(gsub(":.*$", "", colnames(.)))

# ## get.topicFC.vs.RNAexpFC.toPlot <- function(ptb.log2fc.df, expressed.ptb.name, X.full.here) {
# ##     ## ##### function to combine per cell topic data and gene expression data
# ##     ## INPUT
# ##     ## ptb.log2fc.df: topic log2FC df subset to one perturbation only
# ##     ## expressed.ptb.name: the expressed gene
# ##     ## X.full.here: matrix subset to one perturbation only, the values could be RNA expression count, FC, or log2FC.
# ##     ## 
# ##     ## OUTPUT
# ##     ## toPlot: a dataframe that has ptb only cells and one topic column and one expressed gene column
# ##     ##
# ##     expressed.gene.index <- which(grepl(expressed.ptb.name, colnames(X.full.here)))
# ##     expressed.gene.RNA.here <- X.full.here[,expressed.gene.index] %>% as.data.frame %>% `colnames<-`(expressed.ptb.name) 
# ##     toPlot <- merge(ptb.log2fc.df, expressed.gene.RNA.here, by.x="Gene", by.y=0)
# ##     return(toPlot)
# ## }


# ## ## ptb.name.index <- which(grepl(ptb.name, colnames(log2fc.X.full.here)))
# ## ## log2fc.RNA.expression.here <- log2fc.X.full.here[,ptb.name.index] %>% as.data.frame %>% `colnames<-`(ptb.name)
# ## ## RNA.expression.here <- X.full.here[,ptb.name.index] %>% as.data.frame %>% `colnames<-`(ptb.name)
# ## ## topic.RNA.expression.here <- merge(log2fc.here, RNA.expression.here, by.x="Gene", by.y=0)

# ##         if(ptb.name=="MESDC1") ptb.name <- "TLNRD1"
# ##         if(ptb.name %in% colnames(X.full.here)) {
# ##             log2fc.X.full.here <- remove.inf(log2fc.X.full.here, fc.X.full.here, ptb.name)
# ##             topic.RNA.expression.here <- get.topicFC.vs.RNAexpFC.toPlot(topic.log2fc.here, ptb.name, X.full.here)
# ##         } else {
# ##             warning(paste0(ptb.name, " is not in the 10X gene list"))
# ##         }
        
# ##         ## plot topic FC vs RNA expression
# ##         perFACTORFIGDIR <- paste0(FIGDIRSAMPLE, "factor.summary/factor", t, "/")
# ##         if(!dir.exists(perFACTORFIGDIR)) dir.create(perFACTORFIGDIR)
# ##         perFACTORFIGDIR <- paste0(FIGDIRSAMPLE, "factor.summary/factor", t, "/", SAMPLE,"_K",k, "_dt_", DENSITY.THRESHOLD,"_factor", t, "_")
# ##         pdf(file=paste0(perFACTORFIGDIR, "ptb.", ptb.name, "_FC.vs.RNA.expression.pdf"))##here210818
# ##         p <- topic.RNA.expression.here %>% ggplot(aes(x=get(ptb.name), y=get(topic))) + geom_point() + mytheme +
# ##             ggtitle(paste0(ptb.name, " Perturbation")) + 
# ##             xlab(paste0(ptb.name, " RNA Expression")) + ylab(paste0("Topic ", gsub("topic_","",topic), " Expression (log2FC)"))
# ##         print(p)
# ##         dev.off()

# ##         ## EDN1 (or top genes in topic29) expression log2FC in CDH5 KD cells versus topic 29 log2FC in CDH5 KD cells ##todo:210812
# ##         ## get list of top perturbations
# ##         top.expressed.gene.in.topic <- theta.zscore[,t] %>% sort(decreasing=T) %>% head(10) %>% names
# ##         pdf(file=paste0(perFACTORFIGDIR, "ptb.", ptb.name, "_log2FC.vs.top.gene.RNA.expression.pdf"))
# ##         for (expressed.gene.here in top.expressed.gene.in.topic) {
# ##             print(paste0("plotting ", expressed.gene.here))
# ##             ## expressed.gene.here <- "EDN1" # loop over top genes in topic
# ##             ## old (converted to get.topicFC.vs.RNAexpFC.toPlot())
# ##             ## expressed.ptb.name.index <- which(grepl(expressed.gene.here, colnames(log2fc.X.full.here)))
# ##             ## expressed.gene.RNA.here <- X.full.here[,expressed.ptb.name.index] %>% as.data.frame %>% `colnames<-`(expressed.gene.here)
# ##             ## topic.expressed.gene.RNA.here <- merge(topic.log2fc.here, expressed.gene.RNA.here, by.x="Gene", by.y=0)
# ##             topic.expressed.gene.RNA.here <- get.topicFC.vs.RNAexpFC.toPlot(topic.log2fc.here, expressed.gene.here, X.full.here)
# ##             topic.expressed.gene.RNA.fc.here <- get.topicFC.vs.RNAexpFC.toPlot(topic.log2fc.here, expressed.gene.here, fc.X.full.here)
# ##             topic.expressed.gene.RNA.log2fc.here <- get.topicFC.vs.RNAexpFC.toPlot(topic.log2fc.here, expressed.gene.here, log2fc.X.full.here)

# ##             p.raw.x <- topic.expressed.gene.RNA.here %>% ggplot(aes(y=get(expressed.gene.here), x=get(topic))) + geom_point() + mytheme +
# ##                 ggtitle(paste0(ptb.name, " Perturbation")) + 
# ##                 ylab(paste0(expressed.gene.here, " RNA Expression")) + xlab(paste0("Topic ", gsub("topic_","",topic), " Expression (log2FC)"))

# ##             p.fc.x <- topic.expressed.gene.RNA.fc.here %>% ggplot(aes(y=get(expressed.gene.here), x=get(topic))) + geom_point() + mytheme +
# ##                 ggtitle(paste0(ptb.name, " Perturbation")) + 
# ##                 ylab(paste0(expressed.gene.here, " RNA Expression (FC)")) + xlab(paste0("Topic ", gsub("topic_","",topic), " Expression (log2FC)"))

# ##             p.log2fc.x <- topic.expressed.gene.RNA.log2fc.here %>% ggplot(aes(y=get(expressed.gene.here), x=get(topic))) + geom_point() + mytheme +
# ##                 ggtitle(paste0(ptb.name, " Perturbation")) + 
# ##                 ylab(paste0(expressed.gene.here, " RNA Expression (log2FC)")) + xlab(paste0("Topic ", gsub("topic_","",topic), " Expression (log2FC)"))

# ##             ## output plots

# ##             print(p.raw.x) ## fit a line
# ##             print(p.fc.x)
# ##             print(p.log2fc.x)

# ##         }
# ##         dev.off()
# ##     }
# ## }

# ## ##end of scratch:210818



# ## plot by test
# FACTORFIGDIR <- paste0(FIGDIRSAMPLE, "factor.summary/", SAMPLE,"_K",k,"_dt_", DENSITY.THRESHOLD,"_")
# if(!dir.exists(FACTORFIGDIR)) dir.create(FACTORFIGDIR)
# ptb.zscore.long <- ptb.zscore %>% as.data.frame %>% mutate(Gene=rownames(.)) %>% melt(id.vars = "Gene", value.name = "perturbation.zscore", variable.name = "Topic") ## Perturbation z-score plot
# for(test.type.here in c("per.cell.wilcoxon", "per.guide.wilcoxon")) {
#     all.test.w <- all.test %>% subset(test.type==test.type.here)
#     realPvals.df.w <- realPvals.df %>% subset(test.type==test.type.here)
#     for (t in 1:dim(omega)[2]) { #:dim(omega)[2]
#         figure.path <- paste0(FACTORFIGDIR, "factor", t, "_with.sig.", test.type.here, ".0.1.perturbation")
#         ## raw weight version
#         topic <- paste0("topic_",t)
#         toPlot.all.test <- all.test.w %>% subset(Topic==topic)
#         toPlot.fdr <- realPvals.df.w %>% subset(Topic == topic) %>% select(Gene,fdr)
#         ## toPlot <- data.frame(Gene=topFeatures.raw.weight %>% subset(topic == t) %>% pull(Gene),
#         ##                      Score=topFeatures.raw.weight %>% subset(topic == t) %>% pull(scores)) %>%
#         ##     merge(., gene.def.pathways, by="Gene", all.x=T) %>% arrange(desc(Score)) %>% slice(1:10)
#         ## toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
#         toPlot <- data.frame(Gene=topFeatures.raw.weight %>% subset(topic == t) %>% pull(Gene),
#                              Score=topFeatures.raw.weight %>% subset(topic == t) %>% pull(scores)) %>%
#             merge(., summaries %>% select(Symbol, top_class), by.x="Gene", by.y="Symbol", all.x=T) %>% `colnames<-`(c("Gene", "Score", "Pathway")) %>% arrange(desc(Score)) %>% slice(1:10)
#         toPlot$Pathway[is.na(toPlot$Pathway)] <- "Other/Unclassified"
#         toPlot$Pathway[toPlot$Pathway == "unclassified"] <- "Other/Unclassified"
#         p4 <- toPlot %>% ggplot(aes(x=reorder(Gene, Score), y=Score*100, fill=Pathway) ) + geom_col(width=0.5) + theme_minimal()
#         p4 <- p4 + coord_flip() + xlab("Top 10 Genes") + ylab("Raw Weights (z-score)") +
#             mytheme + theme(legend.position="bottom", legend.direction="vertical")

        
#                                         # plot 4
#                                         # add ABC to gene.set.type.df for this particular plot
#         ## gene.set.type.df$type[which(gene.set.type.df$Gene %in% gene.set)] <- "ABC"
#                                         # assemble toPlot
#         toPlot <- gene.score %>% select(all_of(topic)) %>% mutate(Gene=rownames(.)) %>%
#             merge(.,toPlot.all.test,by="Gene", all.x=T) %>%
#             merge(.,toPlot.fdr,by="Gene", all.x=T) %>%
#             ## merge(.,gene.set.type.df,by="Gene") %>%
#             mutate(Gene = gsub("Enhancer-at-CAD-SNP-","",Gene)) %>%
#             merge(., ref.table %>% select("Symbol", "TSS.dist.to.SNP", "GWAS.classification"), by.x="Gene", by.y="Symbol", all.x=T) %>%
#             mutate(EC_ctrl_text = ifelse(.$GWAS.classification == "EC_ctrls", "(+)", "")) %>%
#             mutate(GWAS.class.text = ifelse(grepl("CAD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb"),
#                                             ifelse(grepl("IBD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb_IBD"), ""))) %>%
#             mutate(ann.Gene = paste0(Gene, GWAS.class.text, EC_ctrl_text))
#         colnames(toPlot)[which(colnames(toPlot)==topic)] <- "log2FC"
#         toPlot <- toPlot %>% mutate(significant=ifelse((adjusted.p.value >= fdr.thr | is.na(adjusted.p.value)), "", "*")) %>%
#             arrange(log2FC) %>%
#             mutate(x = seq(nrow(.)))
#         label <- toPlot %>% subset(x <= 3 | x > (nrow(toPlot)-3) | adjusted.p.value < fdr.thr)
#         mytheme <- theme_classic() + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, face = "bold"))
#         p5 <- toPlot %>% ggplot(aes(x=reorder(Gene, log2FC), y=log2FC, color=significant)) + geom_point(size=0.75) + mytheme +
#             theme(axis.ticks.x=element_blank(), axis.text.x=element_blank()) + scale_color_manual(values = c("#E0E0E0", "#38b4f7")) +
#             xlab("Perturbed Genes") + ylab(paste0("Factor ", t, " Expression log2 Fold Change")) +
#             geom_text_repel(data=label, box.padding = 0.5,
#                             aes(label=ann.Gene), size=4,
#                             color="black") +
#             theme(legend.position = "none") # legend at the bottom?

#         ## top perturbation list
#         ## toPlot.all.test <- all.test.w 
#         ## toPlot.fdr <- realPvals.df.w %>% subset(Topic == topic) %>% select(Gene,adjusted.p.value)
#                                         # assemble toPlot


#         toPlot <- gene.score %>% select(all_of(topic)) %>% mutate(Gene=rownames(.)) %>% merge(.,toPlot.all.test,by="Gene", all.x=T) %>%
#             merge(.,toPlot.fdr,by="Gene", all.x=T) %>%
#             merge(.,gene.set.type.df,by="Gene") %>%
#             mutate(Gene = gsub("Enhancer-at-CAD-SNP-","",Gene)) %>%
#             merge(., ref.table %>% select("Symbol", "TSS.dist.to.SNP"), by.x="Gene", by.y="Symbol", all.x=T) %>%
#             mutate(EC_ctrl_text = ifelse(.$type == "EC_ctrls", "(+)", "")) %>%
#             mutate(GWAS.class.text = ifelse(grepl("CAD", type), paste0("_", floor(TSS.dist.to.SNP/1000),"kb"),
#                                             ifelse(grepl("IBD", type), paste0("_", floor(TSS.dist.to.SNP/1000),"kb_IBD"), ""))) %>%
#             mutate(ann.Gene = paste0(Gene, GWAS.class.text, EC_ctrl_text))
#          colnames(toPlot)[which(colnames(toPlot)==topic)] <- "log2FC"

#         ## toPlot <- gene.score %>% select(all_of(topic)) %>% mutate(Gene=rownames(.)) %>% merge(.,toPlot.all.test,by="Gene", all.x=T) %>%
#         ##     merge(.,toPlot.fdr,by="Gene", all.x=T) %>%
#         ##     merge(.,gene.set.type.df,by="Gene") %>%
#         ##     mutate(Gene = gsub("Enhancer-at-CAD-SNP-","",Gene))
#         ## colnames(toPlot)[which(colnames(toPlot)==topic)] <- "log2FC"

#         toPlot <- toPlot %>% mutate(significant=ifelse((adjusted.p.value >= fdr.thr | is.na(adjusted.p.value)), "", "*"))
#         toPlot.top <- toPlot %>% arrange(desc(log2FC)) %>% slice(1:25)
#         toPlot.bottom <- toPlot %>% arrange(log2FC) %>% slice(1:25)
#         toPlot.extreme <- rbind(toPlot.top, toPlot.bottom) %>%
#             mutate(color=ifelse(type %in% c("ABC","CAD focus"), "red",
#                          ifelse(type=="non-expressed", "grey",
#                          ifelse(type=="other", "blue", "black")))) %>%
#              mutate(significant=ifelse((adjusted.p.value >= fdr.thr | is.na(adjusted.p.value)), "", "*"))
#      #       mutate(Gene = ifelse(type=="ABC", paste0("[ ", Gene, " ]"), Gene))
#         p.ptb.list <- toPlot.extreme %>% arrange(desc(log2FC)) %>%  #mutate(Gene = paste0("<span style = 'color: ", color, ";'>", Gene, "</span>")) %>%
#             ggplot(aes(x=reorder(Gene, log2FC), y=log2FC, fill=significant)) + geom_col() + theme_minimal() +
#             coord_flip() + xlab("Most Extreme Gene (Perturbation)") + ylab("log2 Fold Change") +
#             scale_fill_manual(values=c("grey", "#38b4f7")) +
#             geom_text(aes(label = significant)) +
#             theme(legend.position = "none")#, axis.text.y = element_text(colour = toPlot.extreme$color))

        
#         ## ## enriched gene sets by fgsea
#         ## fgsea.here <- fgsea.df.all %>% subset(topic == t) %>% arrange(padj) 
        
#         ## Perturbation z-score
#         toPlot.all.test <- all.test %>% subset(test.type=="per.cell.wilcoxon" & Topic==topic)
#         toPlot.fdr <- realPvals.df %>% subset(test.type=="per.cell.wilcoxon" & Topic == topic) %>% select(Gene,fdr)
#                                         # assemble toPlot
#         toPlot <- ptb.zscore.long %>% subset(Topic == topic) %>% merge(.,toPlot.all.test,by=c("Gene","Topic"), all.x=T) %>%
#             merge(.,toPlot.fdr,by="Gene", all.x=T) %>%
#             merge(.,gene.set.type.df,by="Gene", all.x=T) %>% ##here210809
#             ## merge(.,gene.def.pathways, by="Gene", all.x=T) %>% 
#             merge(., ref.table %>% select("Symbol", "TSS.dist.to.SNP", "GWAS.classification"), by.x="Gene", by.y="Symbol", all.x=T) %>%
#             mutate(EC_ctrl_text = ifelse(.$GWAS.classification == "EC_ctrls", "(+)", "")) %>%
#             mutate(GWAS.class.text = ifelse(grepl("CAD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb"),
#                                      ifelse(grepl("IBD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb_IBD"), ""))) %>%
#             mutate(ann.Gene = paste0(Gene, GWAS.class.text, EC_ctrl_text))

#         toPlot <- toPlot %>% mutate(significant=ifelse((adjusted.p.value >= fdr.thr | is.na(adjusted.p.value)), "", "*"))

#         toPlot.top <- toPlot %>% arrange(desc(perturbation.zscore)) %>% slice(1:25)
#         toPlot.bottom <- toPlot %>% arrange(perturbation.zscore) %>% slice(1:25)
#         toPlot.extreme <- rbind(toPlot.top, toPlot.bottom) %>%
#             mutate(color=ifelse(grepl("CAD", type), "red",
#                          ifelse(type=="non-expressed", "gray",
#                          ifelse(type=="EC_ctrls", "blue", "black"))))  %>%
#             mutate(color=ifelse(is.na(type), "black", color))
#         ## colors <- toPlot.extreme$color[order(toPlot.extreme %>% arrange(desc(perturbation.zscore)) %>% pull(color))]
#         toPlot.extreme <- toPlot.extreme %>% arrange(perturbation.zscore)
#         ## add gene distance to CAD
#         toPlot.extreme$ann.Gene <- factor(toPlot.extreme$ann.Gene, levels = toPlot.extreme$ann.Gene)
#         ## color y.axis.label
#         p.ptb.zscore <- toPlot.extreme %>%  #mutate(Gene = paste0("<span style = 'color: ", color, ";'>", Gene, "</span>")) %>%
#             ggplot(aes(x=ann.Gene, y=perturbation.zscore, fill=significant)) + geom_col() + theme_minimal() +
#             coord_flip() + xlab("Most Extreme Gene (Perturbation)") + ylab("Perturbation z-score") +
#             scale_fill_manual(values=c("grey", "#38b4f7")) +
#             geom_text(aes(label = significant)) +
#             theme(legend.position = "none", axis.text.y = element_text(colour = toPlot.extreme$color))

        
#         ## enriched TF motifs
#         toplot <- all.promoter.ttest.df %>% subset(topic==paste0("topic_",t) & top.gene.mean != 0 & !grepl("X.NA.",motif)) ##here:210816
#         p.promoter.motif <- volcano.plot(toplot, ep.type="promoter", ranking.type="z-score")
#         toplot <- all.enhancer.ttest.df.10en6 %>% subset(topic==paste0("topic_",t) & top.gene.mean != 0 & !grepl("X.NA.",motif))
#         p.enhancer.motif <- volcano.plot(toplot, ep.type="enhancer", ranking.type="z-score")
        
#         ## old as of 210816
#         ## volcano.plot <- function(toplot) {
#         ##     label <- toplot %>% subset(-log10(p.adjust) > 1) %>% mutate(motif.toshow = gsub("HUMAN.H11MO.", "", motif))
#         ##     t <- gsub("topic_", "", toplot$topic[1])
#         ##     p <- toplot %>% ggplot(aes(x=enrichment.log2fc, y=-log10(p.adjust))) + geom_point(size=0.5) + mytheme +
#         ##         ggtitle(paste0("Top 100 z-score Specific Promoter Motif Enrichment")) + xlab("Motif Enrichment (log2FC)") + ylab("-log10(adjusted p-value)") +
#         ##         geom_text_repel(data=label, box.padding = 0.5,
#         ##                         aes(label=motif.toshow), size=5,
#         ##                         color="black")
#         ##     return(p)
#         ## }
#         ## p.motif <- volcano.plot(toplot)

#         ## Factor expression on UMAP
#         plot.features <- paste0("K",k,"_",colnames(omega))
#         feature.name <- plot.features[grepl(paste0("_",t, "$"), plot.features)] ## make sure the seruat object has this feature
#         if ( grepl(feature.name, colnames(s@meta.data)) %>% as.numeric() %>% sum() > 0 ) {
#             p.umap <- FeaturePlot(s, reduction = "umap", features=feature.name)
#         } else {
#             p.umap <- DimPlot(s, reduction = "umap")
#         }

#         p <- ggarrange(ggarrange(p4, p5, p.ptb.list, p.ptb.zscore, nrow=1), ggarrange(p.promoter.motif, p.enhancer.motif, p.umap, nrow=1), nrow=2)
#         ## p.left <- ggarrange(ggarrange(p4, p5, nrow=1), ggarrange(p.motif, p.umap, nrow=1), nrow=2)
#         ## p <- ggarrange(p.left, p.ptb.list, nrow=1, width=c(2.5,1))
#         p <- annotate_figure(p, top = text_grob(paste0("K = ", k, ", Factor ", t), face = "bold", size = 16))

#         ## pdf(paste0(figure.path, ".pdf"), width=10, height=12)
#         ## print(p)
#         ## dev.off()

#         png(paste0(figure.path, ".png"), width=1600, height=1200)
#         print(p)
#         dev.off()
#     }


#     ## ## convert the original pdf to png
#     ## im.convert(pdf.path, output = paste0(pdf.path, ".png"), extra.opts="-density 100") ## takes > 10 minutes
    
# }



# ## java options for allocating more memory to write Excel sheet
# options(java.parameters = "-Xmx16000m") ## 16 GB

# ## make an Excel file that has (topic x top gene x GeneCard summaries) and (topic x top perturbations x GeneCard summaries)
# ## load gene summary files
# theta.zscore.long <- theta.zscore %>% as.data.frame %>% mutate(Gene = rownames(.)) %>% melt(value.name = "score", id.vars="Gene", variable.name = "Factor")
# ## theta.zscore.annotation <- merge(theta.zscore.long %>% group_by(Factor) %>% arrange(desc(score)) %>% slice(1:100), gene.summary)
# ann_omega_long <- gene.score %>% as.data.frame %>% mutate(Gene=rownames(.)) %>% melt(value.name = "log2FC", id.vars = "Gene", variable.name = "Factor")
# ## ann.omega.long <- sqldf("select ann_omega_long.*, summaries.* from ann_omega_long left join summaries on instr(ann_omega_long.Gene, summaries.Symbol)") ## takes a while
# ann.omega.long <- merge(ann_omega_long, summaries, by.x="Gene", by.y="Symbol", all.x=T) %>% merge(., all.test %>% subset(test.type=="per.cell.wilcoxon"), by.x=c("Gene","Factor"), by.y=c("Gene","Topic"), all.x=T)
# ann.omega.long.top <- ann.omega.long %>% group_by(Factor) %>% arrange(desc(log2FC)) %>% slice(1:50) %>% as.data.frame
# ann.omega.long.bottom <- ann.omega.long %>% group_by(Factor) %>% arrange(log2FC) %>% slice(1:50) %>% as.data.frame
# ann.omega.long.sig <- ann.omega.long %>% subset(adjusted.p.value < 0.1) %>% as.data.frame
# ann.omega.long.output <- rbind(ann.omega.long.top, ann.omega.long.bottom, ann.omega.long.sig) %>% unique %>% arrange(Factor, desc(log2FC)) %>% relocate(c(adjusted.p.value,p.value,top_class,classes), .after="log2FC")

# ## load Gavin's top genes in factors with GeneCards information
# theta.zscore.long.output <- read_xlsx(paste0(opt$datadir, "210730_cNMF_topic_model_anal.xlsx"), sheet="top100_annotated")

# ## combine topic definition and perturbation information
# omega.theta.zscore.topic.analysis <- rbind(ann.omega.long.output %>% select(Factor, classes, top_class) %>% mutate(datasource = "omega"),
#                                            theta.zscore.long.output %>% select(Topic, Classes, Top_Class) %>% `colnames<-`(c("Factor", "classes", "top_class")) %>% mutate(datasource = "theta.zscore")) %>% 
#     group_by(Factor, top_class) %>%
#     mutate(top_class_count = n()) %>% ungroup() %>%
#     filter(!is.na(top_class)) %>% ## remove rows with NA in column `top_class`
#     separate_rows(classes, sep=";", convert=T) %>% ## separate classes by ";" and expand the rows
#     filter(!is.na(classes)) %>% subset(classes != "") %>%
#     group_by(Factor, classes) %>%
#     mutate(class_count = n()) %>% ungroup() %>% unique() %>%
#     arrange(Factor, desc(top_class_count), desc(class_count)) 

# ## separate top_class and class into two dfs
# classes.output <- omega.theta.zscore.topic.analysis %>% select(Factor, classes, class_count) %>% as.data.frame %>% unique %>% arrange(Factor,desc(class_count))
# top_class.output <- omega.theta.zscore.topic.analysis %>% select(Factor, top_class, top_class_count) %>% as.data.frame %>% unique %>% arrange(Factor,desc(top_class_count))
# ## write to xlsx
# topic.ptb.summary.xlsx.path <- paste0(OUTDIRSAMPLE, "Significant.or.Top50Ptb_per.cell.wilcoxon.adj.0.1_Summary_", SUBSCRIPT, ".xlsx") 
# write_xlsx(list(Annotations = ann.omega.long.output,
#                 Factor_pathway_summary = omega.theta.zscore.topic.analysis,
#                 Classes_summary = classes.output,
#                 Top_class_summary = top_class.output),
#            path=topic.ptb.summary.xlsx.path, col_names=T)

# ## heatmap for all classes
# classes.heatmap <- classes.output
# classes.heatmap$class_count <- as.numeric(classes.heatmap$class_count)
# classes.heatmap <- classes.heatmap %>% spread(key = classes, value = class_count, fill = as.numeric(0)) %>% `rownames<-`(gsub("topic", "factor",.$Factor)) %>% select(-Factor) %>% as.matrix

# ## plot heatmap
# plotHeatmap <- function(mtx){
#     heatmap.2(
#     mtx, 
#     Rowv=T, 
#     Colv=T,
#     trace='none',
#     key=T,
#     col=palette,
#     labCol=colnames(mtx),
#     margins=c(15,5), 
#     cex.main=0.1, 
#     cexCol=2.5/(nrow(mtx)^(1/3)), cexRow=1.7/(ncol(mtx)^(1/3)),
#     main=paste0(SAMPLE, ", K=", k, ", Factor Pathway Enrichment")
# )
# }

# pdf(paste0(FIGDIRTOP, "factor.pathway.enrichment.pdf"), width= 12, height = 9)
# mtx <- classes.heatmap
# plotHeatmap(mtx)
# mtx <- classes.heatmap %>% apply(2, function(x) x/sum(x)) 
# plotHeatmap(mtx)
# dev.off()




# ptb.array <- c("GOSR2", "TP53", "CDKN1A", "EDN1", "NOS3", "FGD6", "ELN")
# ptb.array <- c("TP53", "ELN", "PHB", "LRPPRC", "MESDC1")

# ## manual QC
# pdf(paste0(FIGDIRSAMPLE, "/manual.QC.pdf"))
# toPlot <- cell.per.ptb <- ann.omega %>% subset(!grepl("^neg|^safe",Gene)) %>% group_by(Gene) %>% summarize(cell.count = n()) # number of cell per perturbation
# p <- toPlot %>% ggplot(aes(x=cell.count)) + stat_ecdf() + mytheme + ggtitle("Number of Cells per Perturbation") + xlab("Number of Cells") + ylab("Fraction of Perturbed Genes")
# print(p)
# p <- toPlot %>% ggplot(aes(x=cell.count)) + geom_histogram() + mytheme + ggtitle("Number of Cells per Perturbation") + xlab("Number of Cells") + ylab("Number of Perturbed Genes")
# print(p)

# toPlot <- guide.per.ptb <- ann.omega %>% subset(!grepl("^neg|^safe",Gene)) %>% select(Gene,Guide) %>% unique() %>% group_by(Gene) %>% summarize(guide.count = n())
# p <- toPlot %>% ggplot(aes(x=guide.count)) + stat_ecdf() + mytheme + ggtitle("Number of Guides per Perturbation") + xlab("Number of Guides") + ylab("Fraction of Perturbed Genes")
# print(p)
# p <- toPlot %>% ggplot(aes(x=guide.count)) + geom_histogram() + mytheme + ggtitle("Number of Guides per Perturbation") + xlab("Number of Guides") + ylab("Number of Perturbed Genes")
# print(p)

# toPlot <- ann.omega %>% group_by(Guide) %>% summarize(count = n())
# p <- toPlot %>% ggplot(aes(x=count)) + stat_ecdf() + mytheme + ggtitle("Number of Cells per Guide") + xlab("Number of Cells") + ylab("Fraction of Perturbed Guides")
# print(p)
# p <- toPlot %>% ggplot(aes(x=count)) + stat_ecdf() + mytheme + ggtitle("Number of Cells per Guide") + xlim(0,20) + xlab("Number of Cells") + ylab("Fraction of Perturbed Guides")
# print(p)
# dev.off()

# ## cell.QC.data <- merge(cell.per.ptb, guide.per.ptb, by="Gene")




# # ##########################################################################
# # ## dotplot of log2FC for perturbation and control
# # ## CDF of log2FC for perturbation and control
# # ## per gene and per cell
# # # create a directory for this type of files (one perturbation per file)
# # FIGDIR.HERE=paste0(FIGDIRTOP,"log2FC_dotplot_CDF_each.perturbation/")
# # if(!dir.exists(FIGDIR.HERE)) dir.create(FIGDIR.HERE, recursive=T)
#  # 
# # gene.fc.ann.omega.ctrl <- fc.ann.omega %>% subset(grepl("^safe-targeting|^negative-control",Gene)) # required for plot.CDF.dotplot()
# # 
# # ptbd.gene <- fc.ann.omega$Gene %>% unique() %>% sort() %>% as.character()
# # for (index in 1:length(ptbd.gene)) {
# #   f <- ptbd.gene[index]
# #   print(paste(index, f, sep=" "))
# #   df <- fc.ann.omega %>% subset(Gene %in% f)
# #   pdf(file=paste0(FIGDIR.HERE, f, ".pdf"), width=16, height=18)
# #   plot.CDF.dotplot(df)
# #   dev.off()
# # }
# # 
# # 
# # 
# # 
# # 
# # 
# # # ## scratch, to delete when comfortable with deleting
# # # # target.gene <- c("Enhancer-at-CAD-SNP-rs17608766","Enhancer-at-CAD-SNP-rs4072980","GOSR2","TP53","ITGB1BP1","KRIT1","TGFBR1","SNAI1","HBB","LBX1")
# # # target.gene <- guideCounts %>% subset(grepl("^Enhancer",Gene.SNP)) %>% pull(Gene.marked) %>% unique()
# # # target.gene <- gsub("E_at_","",target.gene) %>% strsplit(split=",|&") %>% unlist()
# # # # target.gene <- c("PRKCE", "HIF1A", "CDH5", "CTNNA1")
# # # for (index in 1:length(target.gene)) {
# # #   f <- target.gene[index]
# # #   print(paste(index, f, sep=" "))
# # #   df <- fc.ann.omega %>% subset(Gene %in% f)
# # #   pdf(file=paste0(FIGDIR.HERE, f, ".pdf"), width=16, height=18)
# # #   plot.CDF.dotplot(df)
# # #   dev.off()
# # # }
# # # ### end of scratch

