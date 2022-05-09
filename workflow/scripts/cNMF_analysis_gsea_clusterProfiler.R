## Helen Kang
## run clusterProfiler version of GSEA
## 220328

packages <- c("optparse", "data.table", "reshape2", "fgsea", "conflicted", "readxl", "writexl", "org.Hs.eg.db", "tidyr", "dplyr", "clusterProfiler", "msigdbr")
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

## get list of topic defining genes
theta.rank.list <- vector("list", ncol(theta.zscore))## initialize storage list
for(i in 1:ncol(theta.zscore)) {
    topic <- paste0("topic_", colnames(theta.zscore)[i])
    theta.rank.list[[i]] <- theta.zscore[,i] %>%
        as.data.frame %>%
        `colnames<-`("topic.zscore") %>%
        mutate(Gene = rownames(.)) %>%
        arrange(desc(topic.zscore), .before="topic.zscore") %>%
        mutate(zscore.specificity.rank = 1:n()) %>% ## add rank column
        mutate(Topic = topic) ## add topic column
}
theta.rank.df <- do.call(rbind, theta.rank.list) %>%  ## combine list to df
    `colnames<-`(c("topic.zscore", "perturbation", "zscore.specificity.rank", "ProgramID")) %>%
    mutate(ProgramID = gsub("topic_", "K60_", ProgramID)) %>%
    as.data.frame

## get list of topic genes by raw weight
theta.raw.rank.list <- vector("list", ncol(theta.raw))## initialize storage list
for(i in 1:ncol(theta.raw)) {
    topic <- paste0("topic_", colnames(theta.raw)[i])
    theta.raw.rank.list[[i]] <- theta.raw[,i] %>%
        as.data.frame %>%
        `colnames<-`("topic.raw") %>%
        mutate(Gene = rownames(.)) %>%
        arrange(desc(topic.raw), .before="topic.raw") %>%
        mutate(raw.score.rank = 1:n()) %>% ## add rank column
        mutate(Topic = topic) ## add topic column
}
theta.raw.rank.df <- do.call(rbind, theta.raw.rank.list) %>%  ## combine list to df
    `colnames<-`(c("topic.raw", "perturbation", "raw.score.rank", "ProgramID")) %>%
    mutate(ProgramID = gsub("topic_", "K60_", ProgramID)) %>%
    as.data.frame


######################################################################
## run cluster profiler GSEA on top 300 genes

## map between EntrezID and Gene Symbol
z <- org.Hs.egSYMBOL
z_mapped_genes <- mappedkeys(z)
entrez.to.symbol <- as.list(z[z_mapped_genes])

## entrez.to.symbol <- as.list(org.Hs.egSYMBOL)
symbol.to.entrez <- as.list(org.Hs.egSYMBOL2EG)

## map to entrez id (function)
symbolToEntrez <- function(df) df %>% mutate(EntrezID = symbol.to.entrez[.$gene %>% as.character] %>% sapply("[[",1) %>% as.character)

## subset to top 300 genes

ranking.type.ary <- c("zscore", "raw")
ranking.rank.colname.ary <- c("zscore.specificity.rank", "raw.score.rank")
ranking.type.varname.ary <- c("theta.rank.df", "theta.raw.rank.df")


getData <- function(t, i) {
    programID.here <- paste0("K60_", t)

    gene.df <- get(ranking.type.varname.ary[i]) %>%
        subset(ProgramID == programID.here) %>%
        mutate(gene = perturbation) %>%
        symbolToEntrez() %>%
        as.data.frame
    gene.weights <- gene.df %>% pull(get(ranking.rank.colname.ary[i])) %>% `names<-`(top.gene.df$EntrezID)

    top.gene.df <- gene.df %>%
        subset(get(ranking.rank.colname.ary[i]) <= 300) %>%
        as.data.frame
    top.genes <- unlist(mget(top.gene.df$perturbation, envir=org.Hs.egSYMBOL2EG, ifnotfound=NA))
    ## top.genes <- top.gene.df %>% pull(EntrezID) ## same as above
    geneUniverse <- unlist(mget(get(ranking.type.varname.ary[i])$perturbation %>% unique, envir=org.Hs.egSYMBOL2EG, ifnotfound=NA))

    return(list(top.genes = top.genes, geneUniverse = geneUniverse, gene.weights = gene.weights))
}


for(i in 1:length(ranking.type.ary)) {
    ranking.type.here <- ranking.type.ary[i]
    ## GO enrichment analysis. Include all of:
    ## MF: Molecular Function
    ## CC: Cellular Component
    ## BP: Biological Process
    ## ans.go <- do.call(rbind, lapply(1:60, function(t) {

    functionsToRun <- list(GO = "GO <- enrichGO(gene = top.genes, ont = 'ALL', OrgDb = 'org.Hs.eg.db', universe = geneUniverse, readable=T, pvalueCutoff=1, pAdjustMethod = 'fdr') %>% as.data.frame %>% mutate(fdr.across.ont = p.adjust, ProgramID = paste0('K60_', t))",
                           GSEA = "GSEA <- GSEA(gene.weights, TERM2GENE = m_df %>% select(gs_name, entrez_gene), pAdjustMethod = 'fdr', pvalueCutoff = 1) %>% as.data.frame %>% mutate(ProgramID = paste0('K60_', t)) ",
                           enricher.GSEA = "enricher.GSEA <- enricher(top.genes, TERM2GENE = m_df %>% select(gs_name, entrez_gene), universe = geneUniverse, pAdjustMethod = 'fdr', qvalueCutoff=1) %>% as.data.frame %>% mutate(ProgramID = paste0('K60_', t))"
                           )

    out.list <- lapply(1:length(functionsToRun), function(j) {
        out <- do.call(rbind, lapply(c(1:k), function(t) {
            data.here <- getData(t, i)        
            eval(parse(text = functionsToRun[j]))
            return(get(names(functionsToRun)[j]))
        }))
    }) %>% `names<-`(names(functionsToRun))
    
    ## out <- lapply(c(8,15), function(t) {

    ##     data.here <- getData(t, i)        
    ##     out.GO <- enrichGO(gene = top.genes, ont = "ALL", ## ont = "ALL"? # try this
    ##                        OrgDb = "org.Hs.eg.db",
    ##                        universe = geneUniverse,
    ##                        readable=T,
    ##                        pvalueCutoff=1,
    ##                        pAdjustMethod = "fdr") %>%
    ##         as.data.frame %>%
    ##         mutate(fdr.across.ont = p.adjust,
    ##                ProgramID = paste0("K60_", t))
    ## }

    ## out.GSEA.list <- GSEA(gene.weights,
    ##                       TERM2GENE = m_df %>% select(gs_name, entrez_gene),
    ##                       pAdjustMethod = "fdr") %>%
    ##     mutate(ProgramID = paste0("K60_", t))

    ## out.enricher.GSEA.list <- enricher(gene.name.list.here,
    ##                                    TERM2GENE = m_df %>% select(gs_name, entrez_gene),
    ##                                    universe = geneUniverse,
    ##                                    pAdjustMethod = "fdr") %>%
    ##     mutate(ProgramID = paste0("K60_", t))
    
    ## return(list(GO = out.GO.list, GSEA = out.GSEA.list, enricher.GSEA = out.enricher.GSEA.list))
    ## }

    out.GO <- out.list$GO
    out.GSEA <- out.list$GSEA
    out.enricher.GSEA <- out.list$enricher.GSEA
    ## }))
    write.table(out.GO, paste0(OUTDIRSAMPLE, "/clusterProfiler_top300Genes", ranking.type.here, "_GOEnrichment.txt"), sep="\t", row.names=F, quote=F)
    write.table(out.GSEA, paste0(OUTDIRSAMPLE, "/clusterProfiler_allGene", ranking.type.here, "ByWeight_GSEA.txt"), sep="\t", row.names=F, quote=F)
    write.table(out.enricher.GSEA, paste0(OUTDIRSAMPLE, "/clusterProfiler_top300Genes", ranking.type.here, "_GSEA.txt"), sep="\t", row.names=F, quote=F)
    ## write.table(ans.go, paste0(OUTDIRSAMPLE, "/clusterProfiler_top300Genes", ranking.type.here, "_GOEnrichment.txt"), sep="\t", row.names=F, quote=F)
} 

## ## try with the enricher option
## msigdb.c5 <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c5.all.v7.2.entrez.gmt"))
## msigdb.c3 <- gmtPathways(paste0("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/c3.all.v7.2.entrez.gmt"))

## ans.tf <- enricher(gene = top.genes, universe = geneUniverse, TERM2GENE = msigdb.c3, maxGSSize=800, minGSSize=3)
## tab.tf <- as.data.frame(ans.tf)


######################################################################
## ## scratch for creating a df with program 8 defining genes

## ## load gene annotations (Refseq + Uniprot)
## gene.summaries.path <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/data/combined.gene.summaries.txt"
## gene.summary <- read.delim(gene.summaries.path, stringsAsFactors=F)

## t <- 8
## programID.here <- paste0("K60_", t)

## df <- theta.rank.df %>%
##     subset(ProgramID == "K60_8" &
##            zscore.specificity.rank <= 100) %>%
##     merge(gene.summary, by.x="perturbation", by.y="Gene", all.x=T) %>%
##     arrange(zscore.specificity.rank) %>%
##     as.data.frame %>%
##     select(zscore.specificity.rank, ProgramID, perturbation, FullName, Summary)
## colnames(df)[colnames(df) == "perturbation"] <- "Gene"
## write.table(df, "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/Data/4n_K60_8_top100_programDefiningGenes.tsv", sep="\t", quote=F, row.names=F)


## ## end of scratch
