
## Helen Kang
## Topic Model Analysis Only (no plot output)
## 210503
.libPaths("/home/groups/engreitz/Software/R_3.6.1")


packages <- c("optparse","dplyr", "cowplot", "ggplot2", "gplots", "data.table", "reshape2",
              "CountClust", "Hmisc", "tidyr", "grid", "gtable", "gridExtra","ggrepel","ramify",
              "GGally","RNOmni","usedist","ggpubr","gridExtra","GSEA",
              "org.Hs.eg.db","limma","clusterProfiler","fgsea", "conflicted",
              "cluster","textshape","readxl", "IsoplotR", "wesanderson", 
              "ggdist", "gghalves", "Seurat", "writexl")
library(textshape)
library(readxl)
# library(Seurat)
xfun::pkg_attach(packages)
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")


source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")

option.list <- list(
  make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/figures/all_genes/", help="Figure directory"),
  make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/all_genes/", help="Output directory"),
  # make_option("--olddatadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/", help="Input 10x data directory"),
  make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/", help="Input 10x data directory"),
  # make_option("--topic.model.result.dir", type="character", default="/scratch/groups/engreitz/Users/kangh/Perturb-seq_CAD/210707_snakemake_maxParallel/all_genes_acrossK/2kG.library/", help="Topic model results directory"),
  make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
  # make_option("--sep", type="logical", default=F, help="Whether to separate replicates or samples"),
  make_option("--K.list", type="character", default="2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25", help="K values available for analysis"),
  make_option("--K.val", type="numeric", default=60, help="K value to analyze"),
  make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
  make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
  make_option("--ABCdir",type="character", default="/oak/stanford/groups/engreitz/Projects/ABC/200220_CAD/ABC_out/TeloHAEC_Ctrl/Neighborhoods/", help="Path to ABC enhancer directory"),
  make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
  # make_option("--raw.mtx.dir",type="character",default="stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/data/no_IL1B_filtered.normalized.ptb.by.gene.mtx.filtered.txt", help="input matrix to cNMF pipeline"),
  # make_option("--raw.mtx.RDS.dir",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.RDS", help="input matrix to cNMF pipeline"), # the first lane: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210623_aggregate_samples/outputs/aggregated.2kG.library.mtx.cell_x_gene.expandedMultiTargetGuide.RDS"
  # make_option("--subsample.type", type="character", default="", help="Type of cells to keep. Currently only support ctrl"),
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


mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))

SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
DATADIR=opt$olddatadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
# TMDIR=opt$topic.model.result.dir
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
SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr


 
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
if ( motif.enrichment.variables.missing > 0 | opt$recompute ) {
    
    ## load hg38 promoter region file
    promoter.region.hg38.original <- read.table(file=paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS300bp.hg38.bed"), header = F, stringsAsFactors = F)  %>% `colnames<-`(c("chr","start","end", "gene","cell.type","strand")) %>% mutate(sequence_name = paste0(chr, ":", start, "-", end, "(", strand, ")"))

    ## keep only the expressed gene's motifs for background
    expressed.genes <- theta %>% rownames()
    promoter.region.hg38 <- promoter.region.hg38.original %>% mutate(expressed = (gene %in% expressed.genes) ) %>% filter(expressed)

    ## background stats
        print(opt$motif.promoter.background)
        motif.background <- read.delim(file=paste0(ifelse(opt$motif.promoter.background!="", opt$motif.promoter.background, "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv")), header=T, stringsAsFactors=F) %>% filter(!grepl("#", motif_id))  # 30 seconds
    colnames(motif.background)[colnames(motif.background) == "strand"] <- "motif.matched.strand"

    print("start compiling Fisher Exact Test table")

    ## get gene regions
    all.gene.regions <- theta.zscore %>% as.data.frame() %>% mutate(gene = rownames(.)) %>% melt(id.vars="gene", variable.name="topic", value.name="zscore")
    top.gene.regions <- all.gene.regions %>% group_by(topic) %>% arrange(desc(zscore)) %>% slice(1:100)


    ##### This is working! #####
    gene.motif.list.to.tbl <- function(gene.rank.df, motif.background) {
        ## INPUT
        ## gene.rank.df: table with columns "gene" and "topic".
        ## motif.background: table from fimo output, with at least "motif_id" and "sequence_name".
        ##
        ## OUTPUT
        ## list contains a count table and a binary count table with all genes on the rows and all features (motifs and topics) on the columns
        
        # try to make it faster by first merging KL gene regions and all promoter motif matches, then spread
        ## gene.rank.df <- top.gene.regions # to be passed into the function
        long <- merge(gene.rank.df, motif.background, by.x="gene", by.y="sequence_name", all=T) # two minutes
        ## tbl <- long %>% group_by(motif_id, topic) %>% summarize(count=n()) # 5 seconds ish
        tbl2 <- long %>% group_by(gene, motif_id, topic) %>% summarize(count=n()) # three minutes
        wide <- tbl2 %>% spread(key="motif_id", value="count", fill=0) %>% mutate(topic=paste0("topic_", topic), value=1) %>% spread(key="topic", value=value) %>% select(-topic_NA)# 40 seconds
        wide.binary <- tbl2 %>% mutate(count = ifelse(count > 0, 1, 0)) %>% spread(key="motif_id", value="count", fill=0) %>% mutate(topic=paste0("topic_", topic), value=1) %>% spread(key="topic", value=value) %>% select(-topic_NA) # <40 seconds
        wide.binary[is.na(wide.binary)] <- 0
        return(list(table=wide, binary.table=wide.binary))
    }

    ## get topic gene ranking
    promoter.wide.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.", SUBSCRIPT.SHORT, ".txt")
    promoter.wide.binary.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.binary.", SUBSCRIPT.SHORT, ".txt")
    if(file.exists(promoter.wide.path) & file.exists(promoter.wide.binary.path) & !opt$recompute) { 
        promoter.wide <- read.delim(promoter.wide.path, stringsAsFactors=F) 
        promoter.wide.binary <- read.delim(promoter.wide.binary.path, stringsAsFactors=F) 
        } else {
        wide.list <- gene.motif.list.to.tbl(top.gene.regions, motif.background) 
        promoter.wide <- wide.list[[1]]
        promoter.wide.binary <- wide.list[[2]]
        write.table(promoter.wide, file=paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.", SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep="\t") # < 30 seconds
        write.table(promoter.wide.binary, file=paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.binary.", SUBSCRIPT.SHORT, ".txt"), row.names=F, quote=F, sep="\t")
    }


    ## Fisher Exact Test
    fisher.test.on.motifs <- function(wide, wide.binary) {
        motif.ary <- colnames(wide)[!grepl("gene|topic|<NA>", colnames(wide))]
        all.fisher.df <- do.call(rbind, lapply(motif.ary, function(motif.here) {
            do.call(rbind,lapply(1:k, function(t) {
                topic <- paste0("topic_",t)
                tbl.here <- table(wide.binary%>%pull(all_of(motif.here)), wide.binary %>% pull(all_of(topic)))
                fisher.result <- fisher.test(tbl.here)
                fisher.df <- fisher.result %>% broom::glance() %>% mutate(motif=motif.here, topic=topic, in.promoter.global=tbl.here[2,1], global.count = sum(tbl.here[,1]), in.promoter.topic=tbl.here[1,2], topic.count=sum(tbl.here[,2]))
                return(fisher.df)
            })) %>% return()
        })) # < 10 minutes
        all.fisher.df <- all.fisher.df %>% mutate(p.adjust = p.adjust(p.value), in.promoter.global.fraction = in.promoter.global/global.count, in.promoter.topic.fraction = in.promoter.topic / topic.count, enrichment = in.promoter.topic.fraction / in.promoter.global.fraction, enrichment.log2fc = log2(enrichment))
        return(all.fisher.df)
    }

    all.promoter.fisher.df.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.fisher.enrichment.", SUBSCRIPT.SHORT,".txt")
    if(file.exists(all.promoter.fisher.df.path) & !opt$recompute) {
        all.promoter.fisher.df <- read.delim(all.promoter.fisher.df.path, stringsAsFactors=F)
    } else {
        all.promoter.fisher.df <- fisher.test.on.motifs(promoter.wide, promoter.wide.binary)
        write.table(all.promoter.fisher.df, file=all.promoter.fisher.df.path, quote=F,row.names=F,sep="\t")
    }
    
    promoter.wide.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.10en6.", SUBSCRIPT.SHORT, ".txt")
    promoter.wide.binary.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.promoter.motif.table.binary.10en6.", SUBSCRIPT.SHORT, ".txt")
    if(file.exists(promoter.wide.10en6.path) & file.exists(promoter.wide.binary.10en6.path) & !opt$recompute){
        promoter.wide.10en6 <- read.delim(promoter.wide.10en6.path, stringsAsFactors=F)
        promoter.wide.binary.10en6 <- read.delim(promoter.wide.binary.10en6.path, stringsAsFactors=F)
    } else {
        wide.list <- gene.motif.list.to.tbl(top.gene.regions, motif.background %>% subset(`p.value` < 10^-6))  ##debug:210814:why is there a `X.NA.` column?
        promoter.wide.10en6 <- wide.list[[1]]
        promoter.wide.binary.10en6 <- wide.list[[2]]
        write.table(promoter.wide.10en6, file=promoter.wide.10en6.path, row.names=F, quote=F, sep="\t") # < 30 seconds
        write.table(promoter.wide.binary.10en6, file=promoter.wide.binary.10en6.path, row.names=F, quote=F, sep="\t")
    }
    
    ## test for motifs with p-value < 10^-6
    all.promoter.fisher.df.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.fisher.enrichment_motif.thr.10en6_", SUBSCRIPT.SHORT,".txt")
    if(file.exists(all.promoter.fisher.df.10en6.path) & !opt$recompute) {
        all.promoter.fisher.df.10en6 <- read.delim(all.promoter.fisher.df.10en6.path, stringsAsFactors=F)
    } else {
        all.promoter.fisher.df.10en6 <- fisher.test.on.motifs(promoter.wide.10en6, promoter.wide.binary.10en6)
        write.table(all.promoter.fisher.df.10en6, file=all.promoter.fisher.df.10en6.path, quote=F,row.names=F,sep="\t")
    }
    
    
    ## t-test on motif counts
    ttest.on.motifs <- function(wide){##todo:210812:a
        motif.ary <- colnames(wide)[!grepl("gene|topic|<NA>", colnames(wide))]
        all.ttest.df <- do.call(rbind, lapply(motif.ary, function(motif.here) {
            all.topic.ttest.df <- do.call(rbind,lapply(1:k, function(t) {
                topic <- paste0("topic_",t)
                wide.here <- wide %>% select(gene, all_of(topic), all_of(motif.here))
                wide.here[is.na(wide.here)] <- 0
                top.gene.boolean <- wide.here %>% pull(all_of(topic)) == 1
                background.count <- wide.here[which(!top.gene.boolean),] %>% pull(all_of(motif.here))
                top.gene.count <- wide.here[which(top.gene.boolean),] %>% pull(all_of(motif.here))
                test.result <- t.test(top.gene.count, background.count)$p.value
                output <- data.frame(motif=motif.here,
                                     topic=topic,
                                     p.value = test.result) %>%
                    mutate(global.mean=background.count %>% mean,
                           global.var=background.count %>% var,
                           global.count=dim(wide)[1],
                           top.gene.mean=top.gene.count %>% mean,
                           top.gene.var=top.gene.count %>% var,
                           topic.count=top.gene.count %>% length,
                           enrichment=top.gene.mean/global.mean,
                           enrichment.log2fc=log2(enrichment)
                        )                
            }))
        })) %>% mutate(p.adjust = p.adjust(p.value, method="BH"))##here210813
        ## adjust -Inf and +Inf for log2FC
        inf.index <- all.ttest.df %>% pull(enrichment.log2fc) %>% is.infinite %>% which
        neg.inf.index <- ( all.ttest.df %>% pull(top.gene.mean) == 0 ) %>% which
        pos.inf.index <- ( all.ttest.df %>% pull(global.mean) == 0 ) %>% which
        all.ttest.df$enrichment.log2fc[inf.index] <- 0
        min.enrichment.log2fc.value <- min(all.ttest.df$enrichment.log2fc)
        max.enrichment.log2fc.value <- max(all.ttest.df$enrichment.log2fc)
        all.ttest.df$enrichment.log2fc[neg.inf.index] <- min.enrichment.log2fc.value - 3
        all.ttest.df$enrichment.log2fc[pos.inf.index] <- max.enrichment.log2fc.value + 3
        return(all.ttest.df)
    }

    print("Promoter motif count t-test")
    all.promoter.ttest.df.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.count.ttest.enrichment.", SUBSCRIPT.SHORT,".txt")
    if(file.exists(all.promoter.ttest.df.path) & !opt$recompute) {
        all.promoter.ttest.df <- read.delim(all.promoter.ttest.df.path, stringsAsFactors=F)
    } else {
        all.promoter.ttest.df <- ttest.on.motifs(promoter.wide)
        write.table(all.promoter.ttest.df, file=all.promoter.ttest.df.path, quote=F,row.names=F,sep="\t")
    }
    
    print("Promoter motif (p-value < 10^-6)  count t-test")
    all.promoter.ttest.df.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.count.ttest.enrichment_motif.thr.10en6_", SUBSCRIPT.SHORT,".txt")
    if(file.exists(all.promoter.ttest.df.10en6.path) & !opt$recompute) {
        all.promoter.ttest.df.10en6 <- read.delim(all.promoter.ttest.df.10en6.path, stringsAsFactors=F)
    } else {
        all.promoter.ttest.df.10en6 <- ttest.on.motifs(promoter.wide.10en6)
        write.table(all.promoter.ttest.df.10en6, file=all.promoter.ttest.df.10en6.path, quote=F,row.names=F,sep="\t")
    }  

    
    ##  motif enrichment 
    ## load enhancer motif matches
        print(opt$motif.enhancer.background)
        motif.background <- read.delim(file=paste0(ifelse(opt$motif.enhancer.background!="", opt$motif.enhancer.background, "/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv")), header=F, stringsAsFactors=F) %>%
        `colnames<-`(c("motif_id", "motif_alt_id", "enhancer_region", "enhancer_type", "gene_region","sequence_name","start","stop","strand","score","p-value","q-value","matched_sequence")) %>% filter(!grepl("#|motif_id", motif_id))  # more than 30 seconds, minutes?
    colnames(motif.background)[colnames(motif.background) == "strand"] <- "motif.matched.strand"
        motif.background <- motif.background %>% filter(!grepl("promoter", enhancer_type))
        
    print("start compiling enhancer Fisher Exact Test table")

    enhancer.wide.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.table.", SUBSCRIPT.SHORT, ".txt")
    enhancer.wide.binary.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.table.binary.", SUBSCRIPT.SHORT, ".txt")
    if(file.exists(enhancer.wide.path) & file.exists(enhancer.wide.binary.path) & !opt$recompute) {
        enhancer.wide <- read.delim(enhancer.wide.path, stringsAsFactors=F)
        enhancer.wide.binary <- read.delim(enhancer.wide.binary.path, stringsAsFactors=F)
    } else {
        wide.list <- gene.motif.list.to.tbl(top.gene.regions, motif.background) 
        enhancer.wide <- wide.list[[1]]
        enhancer.wide.binary <- wide.list[[2]]
        write.table(enhancer.wide, file=enhancer.wide.path, row.names=F, quote=F, sep="\t") # < 30 seconds
        write.table(enhancer.wide.binary, file=enhancer.wide.binary.path, row.names=F, quote=F, sep="\t")
    }



    ## Fisher Exact Test for Enhancers
    all.enhancer.fisher.df.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.motif.fisher.enrichment.", SUBSCRIPT.SHORT,".txt")
    if(file.exists(all.enhancer.fisher.df.path) & !opt$recompute) {
        all.enhancer.fisher.df <- read.delim(all.enhancer.fisher.df.path, stringsAsFactors=F)
    } else {
        all.enhancer.fisher.df <- fisher.test.on.motifs(enhancer.wide, enhancer.wide.binary)
        write.table(all.enhancer.fisher.df, file=all.enhancer.fisher.df.path, quote=F,row.names=F,sep="\t")
    }

    print("started compiling enhancer fisher exact test table for motif with p-value < 10^-6")
    
    enhancer.wide.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.table.10en6.", SUBSCRIPT.SHORT, ".txt")
    enhancer.wide.binary.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.table.binary.10en6.", SUBSCRIPT.SHORT, ".txt")
    if(file.exists(enhancer.wide.10en6.path) & file.exists(enhancer.wide.binary.10en6.path) & !opt$recompute) {
        enhancer.wide.10en6 <- read.delim(enhancer.wide.10en6.path, stringsAsFactors=F)
        enhancer.wide.binary.10en6 <- read.delim(enhancer.wide.binary.10en6.path, stringsAsFactors=F)
    } else {
        wide.list <- gene.motif.list.to.tbl(top.gene.regions, motif.background %>% subset(`p-value` < 10^-6)) 
        enhancer.wide.10en6 <- wide.list[[1]]
        enhancer.wide.binary.10en6 <- wide.list[[2]]
        write.table(enhancer.wide.10en6, file=enhancer.wide.10en6.path, row.names=F, quote=F, sep="\t") # < 30 seconds
        write.table(enhancer.wide.binary.10en6, file=enhancer.wide.binary.10en6.path, row.names=F, quote=F, sep="\t")
    }
    
    ## test for motifs with p-value < 10^-6
    all.enhancer.fisher.df.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.fisher.enrichment_motif.thr.10en6_", SUBSCRIPT.SHORT,".txt")
    if(file.exists(all.enhancer.fisher.df.10en6.path) & !opt$recompute) {
        all.enhancer.fisher.df.10en6 <- read.delim(all.enhancer.fisher.df.10en6.path, stringsAsFactors=F)
    } else {
        all.enhancer.fisher.df.10en6 <- fisher.test.on.motifs(enhancer.wide.10en6, enhancer.wide.binary.10en6)
        write.table(all.enhancer.fisher.df.10en6, file=all.enhancer.fisher.df.10en6.path, quote=F,row.names=F,sep="\t")
    }
    
    print("Enhancer motif count t-test")
    all.enhancer.ttest.df.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.count.ttest.enrichment.", SUBSCRIPT.SHORT,".txt")
    if(file.exists(all.enhancer.ttest.df.path) & !opt$recompute) {
        all.enhancer.ttest.df <- read.delim(all.enhancer.ttest.df.path, stringsAsFactors=F)
    } else {
        all.enhancer.ttest.df <- ttest.on.motifs(enhancer.wide)
        write.table(all.enhancer.ttest.df, file=all.enhancer.ttest.df.path, quote=F,row.names=F,sep="\t")
    }
    
    print("Enhancer motif (p-value < 10^-6)  count t-test")
    all.enhancer.ttest.df.10en6.path <- paste0(OUTDIRSAMPLE,"/topic.top.100.zscore.gene.ABC.enhancer.motif.count.ttest.enrichment_motif.thr.10en6_", SUBSCRIPT.SHORT,".txt")
    if(file.exists(all.enhancer.ttest.df.10en6.path) & !opt$recompute) {
        all.enhancer.ttest.df.10en6 <- read.delim(all.enhancer.ttest.df.10en6.path, stringsAsFactors=F)
    } else {
        all.enhancer.ttest.df.10en6 <- ttest.on.motifs(enhancer.wide.10en6)
        write.table(all.enhancer.ttest.df.10en6, file=all.enhancer.ttest.df.10en6.path, quote=F,row.names=F,sep="\t")
    }
    

    save(all.enhancer.fisher.df, all.promoter.fisher.df, promoter.wide, enhancer.wide, promoter.wide.binary, enhancer.wide.binary,
         enhancer.wide.10en6, enhancer.wide.binary.10en6, all.enhancer.fisher.df.10en6,
         promoter.wide.10en6, promoter.wide.binary.10en6, all.promoter.fisher.df.10en6,
         all.promoter.ttest.df, all.promoter.ttest.df.10en6,
         all.enhancer.ttest.df, all.enhancer.ttest.df.10en6,
     file=paste0(file.name))


    
}

}


