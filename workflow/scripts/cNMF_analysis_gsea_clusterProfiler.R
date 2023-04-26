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
  make_option("--recompute", type="logical", default=F, help="T for recomputing statistical tests and F for not recompute"),

  ## GSEA parameters
  make_option("--ranking.type", type="character", default="zscore", help="{zscore, raw} ranking for the top program genes"),
  make_option("--GSEA.type", type="character", default="GOEnrichment", help="{GOEnrichment, ByWeightGSEA, GSEA}")
  ## make_option("--", type="", default= , help="")
  
)
opt <- parse_args(OptionParser(option_list=option.list))

## ## all genes directories (for sdev)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/"
## opt$K.val <- 60

## ## ## K562 gwps sdev
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/"
## opt$K.val <- 35
## opt$sampleName <- "WeissmanK562gwps"
## opt$GSEA.type <- "ByWeightGSEA"
## opt$ranking.type <- "median_spectra_zscore"

## ## ENCODE mouse heart
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/figures/top2000VariableGenes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/IGVF/Cellular_Programs_Networks/230116_snakemake_mouse_ENCODE_heart/analysis/top2000VariableGenes"
## opt$K.val <- 15
## opt$sampleName <- "mouse_ENCODE_heart"

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
## FGSEADIR=paste0(OUTDIRSAMPLE,"/fgsea/")
## FGSEAFIG=paste0(FIGDIRSAMPLE,"/fgsea/")

## subscript for files
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
# SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

## adjusted p-value threshold
fdr.thr <- opt$adj.p.value.thr
p.value.thr <- opt$adj.p.value.thr

db <- ifelse(grepl("mouse", SAMPLE), "org.Mm.eg.db", "org.Hs.eg.db")
library(!!db) ## load the appropriate database

# create dir if not already
check.dir <- c(OUTDIR, FIGDIR, paste0(FIGDIR,SAMPLE,"/"), paste0(FIGDIR,SAMPLE,"/K",k,"/"), paste0(OUTDIR,SAMPLE,"/"), OUTDIRSAMPLE, FIGDIRSAMPLE)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

palette = colorRampPalette(c("#38b4f7", "white", "red"))(n = 100)


## helper function to map between ENSGID and SYMBOL
map.ENSGID.SYMBOL <- function(df) {
    ## need column `Gene` to be present in df
    ## detect gene data type (e.g. ENSGID, Entrez Symbol)
    gene.type <- ifelse(nrow(df) == sum(as.numeric(grepl("^ENS", df$Gene))),
                        "ENSGID",
                        "Gene")
    if(gene.type == "ENSGID") {
        mapped.genes <- mapIds(get(db), keys=df$Gene, keytype = "ENSEMBL", column = "SYMBOL")
        df <- df %>% mutate(ENSGID = Gene, Gene = mapped.genes)
    } else {
        mapped.genes <- mapIds(get(db), keys=df$Gene, keytype = "SYMBOL", column = "ENSEMBL")
        df <- df %>% mutate(ENSGID = mapped.genes)
    }
    return(df)
}



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
    theta.rank.list[[i]] <- theta.zscore %>%
        as.data.frame %>%
        select(all_of(i)) %>%
        `colnames<-`("topic.zscore") %>%
        mutate(Gene = rownames(.)) %>%
        arrange(desc(topic.zscore), .before="topic.zscore") %>%
        mutate(zscore.specificity.rank = 1:n()) %>% ## add rank column
        mutate(Topic = topic) ## add topic column
}
theta.rank.df <- do.call(rbind, theta.rank.list) %>%  ## combine list to df
    `colnames<-`(c("topic.zscore", "Gene", "zscore.specificity.rank", "ProgramID")) %>%
    mutate(ProgramID = gsub("topic_", paste0("K", k, "_"), ProgramID)) %>%
    as.data.frame %>% map.ENSGID.SYMBOL

## get list of topic genes by raw weight
theta.raw.rank.list <- vector("list", ncol(theta.raw))## initialize storage list
for(i in 1:ncol(theta.raw)) {
    topic <- paste0("topic_", colnames(theta.raw)[i])
    theta.raw.rank.list[[i]] <- theta.raw %>%
        as.data.frame %>%
        select(all_of(i)) %>%
        `colnames<-`("topic.raw") %>%
        mutate(Gene = rownames(.)) %>%
        arrange(desc(topic.raw), .before="topic.raw") %>%
        mutate(raw.score.rank = 1:n()) %>% ## add rank column
        mutate(Topic = topic) ## add topic column
}
theta.raw.rank.df <- do.call(rbind, theta.raw.rank.list) %>%  ## combine list to df
    `colnames<-`(c("topic.raw", "Gene", "raw.score.rank", "ProgramID")) %>%
    mutate(ProgramID = gsub("topic_", paste0("K", k, "_"), ProgramID)) %>%
    as.data.frame %>% map.ENSGID.SYMBOL


## get list of topic genes by median spectra weight
median.spectra.rank.list <- vector("list", ncol(median.spectra))## initialize storage list
for(i in 1:ncol(median.spectra)) {
    topic <- paste0("topic_", colnames(median.spectra)[i])
    median.spectra.rank.list[[i]] <- median.spectra %>%
        as.data.frame %>%
        select(all_of(i)) %>%
        `colnames<-`("median.spectra") %>%
        mutate(Gene = rownames(.)) %>%
        arrange(desc(median.spectra), .before="median.spectra") %>%
        mutate(median.spectra.rank = 1:n()) %>% ## add rank column
        mutate(Topic = topic) ## add topic column
}
median.spectra.rank.df <- do.call(rbind, median.spectra.rank.list) %>%  ## combine list to df
    `colnames<-`(c("median.spectra", "Gene", "median.spectra.rank", "ProgramID")) %>%
    mutate(ProgramID = gsub("topic_", paste0("K", k, "_"), ProgramID)) %>%
    as.data.frame %>% map.ENSGID.SYMBOL
## median.spectra.zscore.df <- median.spectra.zscore.df %>% mutate(Gene = ENSGID) ## quick fix, need to add "Gene" column to this dataframe in analysis script

######################################################################
## run cluster profiler GSEA on top 300 genes

## ## map between EntrezID and Gene Symbol
## z <- org.Hs.egSYMBOL
## z_mapped_genes <- mappedkeys(z)
## entrez.to.symbol <- as.list(z[z_mapped_genes])

## ## entrez.to.symbol <- as.list(org.Hs.egSYMBOL)
## symbol.to.entrez <- as.list(org.Hs.egSYMBOL2EG)

## ## map to entrez id (function)
## symbolToEntrez <- function(df) df %>% mutate(EntrezID = symbol.to.entrez[.$gene %>% as.character] %>% sapply("[[",1) %>% as.character)

## subset to top 300 genes

ranking.type.ary <- c("zscore", "raw", "median_spectra_zscore", "median_spectra")
score.colname.ary <- c("zscore", "raw", "median.spectra.zscore", "median.spectra")
ranking.rank.colname.ary <- c("zscore.specificity.rank", "raw.score.rank", "median.spectra.zscore.rank", "median.spectra.rank")
ranking.type.varname.ary <- c("theta.rank.df", "theta.raw.rank.df", "median.spectra.zscore.df", "median.spectra.rank.df")
db <- ifelse(grepl("mouse", SAMPLE), "org.Mm.eg.db", "org.Hs.eg.db")
library(!!db) ## load the appropriate database

getData <- function(t) {
    i <- which(ranking.type.ary == opt$ranking.type)
    programID.here <- paste0("K", k, "_", t)
    ranking.type.varname.here <- ranking.type.varname.ary[i]
    if(grepl("median.spectra", ranking.type.varname.here)) {
       ranking.score.colname.here <- score.colname.ary[i] 
    } else {
        ranking.score.colname.here <- paste0("topic.", ranking.type.ary[i])
    }
    gene.df <- get(ranking.type.varname.here) %>%
        subset(ProgramID == programID.here)
    gene.type <- ifelse(nrow(gene.df) == sum(as.numeric(grepl("^ENS", gene.df$Gene))), "ENSGID", "Gene")
    mapped.genes <- mapIds(get(db),
                           keys=gene.df$Gene,
                           keytype = ifelse(gene.type == "Gene", "SYMBOL", "ENSEMBL"),
                           column = ifelse(gene.type == "Gene", "ENSEMBL", "SYMBOL"))
    mapped.entrez.genes <- mapIds(get(db),
                           keys=gene.df$Gene,
                           keytype = ifelse(gene.type == "Gene", "SYMBOL", "ENSEMBL"),
                           column = "ENTREZID")
    gene.df <- gene.df %>%
        mutate(!!gene.type := Gene,
               !!ifelse(gene.type=="ENSGID", "Gene", "ENSGID") := mapped.genes,
               EntrezID = mapped.entrez.genes) %>%
        as.data.frame


    gene.weights <- gene.df %>% pull(get(ranking.score.colname.here)) %>% `names<-`(gene.df$EntrezID)
    gene.weights[gene.weights < 0] <- 0    

    top.gene.df <- gene.df %>%
        subset(get(ranking.rank.colname.ary[i]) <= 300) %>%
        as.data.frame
    ## top.genes <- unlist(mget(top.gene.df$Gene, envir=org.Hs.egSYMBOL2EG, ifnotfound=NA)) ## old
    top.genes <- top.gene.df %>% pull(EntrezID) ## same as above
    
    pos.gene.df <- gene.df %>%
        subset(get(ranking.score.colname.here) > 0) %>%
        as.data.frame
    pos.genes <- pos.gene.df %>% pull(EntrezID)
    ## pos.genes <- unlist(mget(pos.gene.df$Gene, envir=org.Hs.egSYMBOL2EG, ifnotfound=NA))

    ## geneUniverse <- unlist(mget(get(ranking.type.varname.ary[i])$Gene %>% unique, envir=org.Hs.egSYMBOL2EG, ifnotfound=NA))
    geneUniverse <- gene.df$EntrezID
    
    return(list(top.genes = top.genes, pos.genes = pos.genes, geneUniverse = geneUniverse, gene.weights = gene.weights))
}

m_df <- msigdbr(species = ifelse(grepl("mouse", SAMPLE), "Mus musculus", "Homo sapiens"))

## save this as a txt file and read in ## for future if needed
functionsToRun <- list(GOEnrichment = "out <- enrichGO(gene = top.genes, ont = 'ALL', OrgDb = db, universe = geneUniverse, readable=T, pvalueCutoff=1, pAdjustMethod = 'fdr') %>% as.data.frame %>% mutate(fdr.across.ont = p.adjust, ProgramID = paste0('K', k, '_', t))",
                       PosGenesGOEnrichment = "out <- enrichGO(gene = pos.genes, ont = 'ALL', OrgDb = db, universe = geneUniverse, readable=T, pvalueCutoff=1, pAdjustMethod='fdr') %>% as.data.frame %>% mutate(fdr.across.ont = p.adjust, ProgramID = paste0('K', k, '_', t))",
                       ByWeightGSEA = "out <- GSEA(gene.weights, TERM2GENE = m_df %>% select(gs_name, entrez_gene), pAdjustMethod = 'fdr', pvalueCutoff = 1) %>% as.data.frame %>% mutate(ProgramID = paste0('K', k, '_', t)) ",
                       GSEA = "out <- enricher(top.genes, TERM2GENE = m_df %>% select(gs_name, entrez_gene), universe = geneUniverse, pAdjustMethod = 'fdr', qvalueCutoff=1) %>% as.data.frame %>% mutate(ProgramID = paste0('K', k, '_', t))"
                       )

## for(i in 1:length(ranking.type.ary)) {
ranking.type.here <- opt$ranking.type
GSEA.type <- opt$GSEA.type
## ranking.type.here <- ranking.type.ary[i]
## GO enrichment analysis. Include all of:
## MF: Molecular Function
## CC: Cellular Component
## BP: Biological Process
## ans.go <- do.call(rbind, lapply(1:60, function(t) {

message("starting enrichment")

## out.list <- lapply(1:length(functionsToRun), function(j) {
out <- do.call(rbind, lapply(c(1:k) %>% rev, function(t) {
    data.here <- getData(t)
    top.genes <- data.here$top.genes
    if(sum(as.numeric(is.na(names(top.genes)))) > 0) top.genes <- top.genes[-which(is.na(names(top.genes)))] ## remove genes that doesn't have matched Entrez ID
    ## print(head(top.genes))

    pos.genes <- data.here$pos.genes
    if(sum(as.numeric(is.na(names(pos.genes)))) > 0) pos.genes <- pos.genes[-which(is.na(names(pos.genes)))]
    ## print(head(pos.genes))

    geneUniverse <- data.here$geneUniverse
    ## print(head(geneUniverse))
    
    gene.weights <- data.here$gene.weights
    if(sum(as.numeric(is.na(names(gene.weights)))) > 0) gene.weights <- gene.weights[-which(is.na(names(gene.weights)))]
    gene.weights <- gene.weights[-which(gene.weights==0)] ## can't have zero weights?
    ## print(head(gene.weights))
    
    message(paste0("Ranking type: ", ranking.type.here, ", Program ", t, ", out of ", k, ", function ", GSEA.type, ", top gene class: ", class(top.genes),
                   "\n geneUniverse class: ", class(geneUniverse), ", gene.weights class: ", class(gene.weights)))
    ## message(paste0("Function to run: \n", functionsToRun[[GSEA.type]]))
    eval(parse(text = functionsToRun[[GSEA.type]]))
    return(out)
}))
## if(j == 1) {
file.name <- paste0(OUTDIRSAMPLE, "/clusterProfiler_GeneRankingType", ranking.type.here, "_EnrichmentType", GSEA.type,".txt")
## } else if (j == 2) {
##     file.name <- paste0(OUTDIRSAMPLE, "/clusterProfiler_allGene", ranking.type.here, "_ByWeight_GSEA.txt")
## } else {
##     file.name <- paste0(OUTDIRSAMPLE, "/clusterProfiler_top300Genes", ranking.type.here, "_GSEA.txt")
## }
message(paste0("output table to ", file.name))
write.table(out, file.name, sep="\t", row.names=F, quote=F)

## }) %>% `names<-`(names(functionsToRun)) ## lapply and save in a function?

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
##                       TERM2GENE = m_dfo %>% select(gs_name, entrez_gene),
##                       pAdjustMethod = "fdr") %>%
##     mutate(ProgramID = paste0("K60_", t))

## out.enricher.GSEA.list <- enricher(gene.name.list.here,
##                                    TERM2GENE = m_df %>% select(gs_name, entrez_gene),
##                                    universe = geneUniverse,
##                                    pAdjustMethod = "fdr") %>%
##     mutate(ProgramID = paste0("K60_", t))

## return(list(GO = out.GO.list, GSEA = out.GSEA.list, enricher.GSEA = out.enricher.GSEA.list))
## }

## out.GO <- out.list$GO
## out.GSEA <- out.list$GSEA
## ## out.enricher.GSEA <- out.list$enricher.GSEA
## ## }))
## write.table(out.GO, paste0(OUTDIRSAMPLE, "/clusterProfiler_top300Genes", ranking.type.here, "_GOEnrichment.txt"), sep="\t", row.names=F, quote=F)
## write.table(out.GSEA, paste0(OUTDIRSAMPLE, "/clusterProfiler_allGene", ranking.type.here, "ByWeight_GSEA.txt"), sep="\t", row.names=F, quote=F)
## write.table(out.enricher.GSEA, paste0(OUTDIRSAMPLE, "/clusterProfiler_top300Genes", ranking.type.here, "_GSEA.txt"), sep="\t", row.names=F, quote=F)
## write.table(ans.go, paste0(OUTDIRSAMPLE, "/clusterProfiler_top300Genes", ranking.type.here, "_GOEnrichment.txt"), sep="\t", row.names=F, quote=F)
## } 

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
##     merge(gene.summary, by.x="Gene", by.y="Gene", all.x=T) %>%
##     arrange(zscore.specificity.rank) %>%
##     as.data.frame %>%
##     select(zscore.specificity.rank, ProgramID, Gene, FullName, Summary)
## write.table(df, "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/Data/4n_K60_8_top100_programDefiningGenes.tsv", sep="\t", quote=F, row.names=F)


## ## end of scratch

