## Helen Kang
## Aggregate results from topic model analysis across K
## 210630 referenced: 210514
## .libPaths("/home/groups/engreitz/Software/R_3.6.1")


packages <- c("optparse","dplyr", "data.table", "reshape2", "conflicted","ggplot2",
              "tidyr", "textshape","readxl") # , "IsoplotR"
xfun::pkg_attach(packages)
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")

## source("/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModelAnalysis.functions.R")

option.list <- list(
  make_option("--figdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/", help="Figure directory"), # "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/figures/all_genes"
  make_option("--outdir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/", help="Output directory"), # "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/analysis/all_genes/FT010_fresh_2min"
  make_option("--datadir", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/data/", help="Input 10x data directory"), # "/oak/stanford/groups/engreitz/Users/kangh/process_sequencing_data/210912_FT010_fresh_Telo_sortedEC/multiome_FT010_fresh_2min/outs/filtered_feature_bc_matrix"
  make_option("--sampleName", type="character", default="2kG.library", help="Name of Samples to be processed, separated by commas"),
  # make_option("--sep", type="logical", default=F, help="Whether to separate replicates or samples"),
  make_option("--K.list", type="character", default="14,30,60", help="K values available for analysis"),
  # make_option("--K.val", type="numeric", default=14, help="K value to analyze"),
  make_option("--K.table", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K.spectra.threshold.table.txt", help="table for defining spectra threshold"), # opt$K.table <-"/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/all_genes/2kG.library/K.spectra.threshold.table.txt"
  # make_option("--cell.count.thr", type="numeric", default=2, help="filter threshold for number of cells per guide (greater than the input number)"),
  # make_option("--guide.count.thr", type="numeric", default=1, help="filter threshold for number of guide per perturbation (greater than the input number)"),
  # make_option("--ABCdir",type="character", default="/oak/stanford/groups/engreitz/Projects/ABC/200220_CAD/ABC_out/TeloHAEC_Ctrl/Neighborhoods/", help="Path to ABC enhancer directory"),
  # make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
  # make_option("--raw.mtx.dir",type="character",default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/data/no_IL1B_filtered.normalized.ptb.by.gene.mtx.filtered.txt", help="input matrix to cNMF pipeline"),
  # make_option("--subsample.type", type="character", default="", help="Type of cells to keep. Currently only support ctrl"),

  ## fisher motif enrichment
  ## make_option("--outputTable", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2105_findK/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.k_14.df_0_2.txt", help="Output directory"),
  ## make_option("--outputTableBinary", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2105_findK/outputs/no_IL1B/topic.top.100.zscore.gene.motif.table.binary.k_14.df_0_2.txt", help="Output directory"),
  ## make_option("--outputEnrichment", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2105_findK/outputs/no_IL1B/topic.top.100.zscore.gene.motif.fisher.enrichment.k_14.df_0_2.txt", help="Output directory"),
  # make_option("--motif.promoter.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2104_remove_lincRNA/data/fimo_out_all_promoters_thresh1.0E-4/fimo.tsv", help="All promoter's motif matches"),
  # make_option("--motif.enhancer.background", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2105_findK/data/fimo_out_ABC_TeloHAEC_Ctrl_thresh1.0E-4/fimo.formatted.tsv", help="All enhancer's motif matches specific to {no,plus}_IL1B"),
  make_option("--adj.p.value.thr", type="numeric", default=0.1, help="adjusted p-value threshold"),
  make_option("--recompute", type="logical", default=F, help="T for recomputing statistical tests and F for not recompute")
  
)
opt <- parse_args(OptionParser(option_list=option.list))

## sdev for 2n1.99x singlets
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/figures/all_genes"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/"
## ## opt$datadir <- "/oak/stanford/groups/engreitz/Users/kangh/process_sequencing_data/210912_FT010_fresh_Telo_sortedEC/multiome_FT010_fresh_2min/outs/filtered_feature_bc_matrix"
## opt$sampleName <- "Perturb_2kG_dup4"

mytheme <- theme_classic() + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11), plot.title = element_text(hjust = 0.5, face = "bold"))


## ## all genes (for interactive sessions)
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/"
## opt$K.table <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K.spectra.threshold.table.txt"

## ## overdispersed Genes
## opt$figdir <- "//oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220716_snakemake_overdispersedGenes/figures/top2000VariableGenes/"
## opt$outdir <- "//oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220716_snakemake_overdispersedGenes/analysis/top2000VariableGenes/"
## opt$K.table <- ""
## opt$K.list <- "3,5,7,12,14,19,21,23,25,27,29,31,35,40,45,50,60,70,80,90,100,120"
## opt$sampleName <- "2kG.library_overdispersedGenes"

## ## sdev K562 gwps 2k most dispersed genes cNMF
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes"
## opt$outdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes"
## opt$sampleName <- "WeissmanK562gwps"
## opt$K.list <- "3,5,10,15,20,25,30,35,40,45,50,55,60,70,80,90,100,110,120"
## opt$K.table <- "/to/use/for/specifying/spectra/cut/off/threshold//default/0.2"

SAMPLE=strsplit(opt$sampleName,",") %>% unlist()
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
# STATIC.SAMPLE=c("Telo_no_IL1B_T200_1", "Telo_no_IL1B_T200_2", "Telo_plus_IL1B_T200_1", "Telo_plus_IL1B_T200_2", "no_IL1B", "plus_IL1B",  "pooled")
DATADIR=opt$datadir # "/seq/lincRNA/Gavin/200829_200g_anal/scRNAseq/"
OUTDIR=opt$outdir
OUTDIR.ACROSS.K=paste0(OUTDIR,"/",SAMPLE,"/acrossK/")
K.list <- strsplit(opt$K.list,",") %>% unlist() %>% as.numeric()
FIGDIR=opt$figdir



## adjusted p-value threshold
p.value.thr <- opt$adj.p.value.thr



## # create dir if not already
check.dir <- c(OUTDIR, FIGDIR, OUTDIR.ACROSS.K)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x) }))


# K spectra threshold table
if(file.exists(opt$K.table)) {
    K.spectra.threshold <- read.table(file=paste0(opt$K.table), header=T, stringsAsFactors=F)
} else {
    K.spectra.threshold <- data.frame(K = K.list, density.threshold=rep(0.2, length(K.list)))  ## assume 0.2 is the best threshold for filtering out outlier topics
}



## initialize storage variables
GSEA.types <- c("GOEnrichment", "PosGenesGOEnrichment", "ByWeightGSEA", "GSEA")
for (j in 1:length(GSEA.types)) {
    GSEA.type <- GSEA.types[j]
    to.eval <- paste0("clusterProfiler.", GSEA.type, ".list <- vector(\"list\",nrow(K.spectra.threshold))")
    eval(parse(text = to.eval))
}
all.MAST.df.list <- all.test.df.list <- all.fdr.df.list <- count.by.GWAS.list <- count.by.GWAS.withTopic.list <- theta.zscore.list <- theta.raw.list <-  all.promoter.ttest.df.list <- all.enhancer.ttest.df.list <- varianceExplainedByModel.list <- varianceExplainedPerProgram.list <- vector("list", nrow(K.spectra.threshold))
## loop over all values of K and aggregate results
for (n in 1:nrow(K.spectra.threshold)) {
    k <- K.spectra.threshold[n,"K"]
    DENSITY.THRESHOLD <- K.spectra.threshold[n,"density.threshold"] %>% gsub("\\.","_",.)

    ## subscript for files
    SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
    # SUBSCRIPT=paste0("k_", k,".dt_",DENSITY.THRESHOLD,".minGuidePerPtb_",opt$guide.count.thr,".minCellPerGuide_", opt$cell.count.thr)

    ## directories
    # FIGDIRSAMPLE=ifelse(SEP, paste0(FIGDIR,SAMPLE,".sep/K",k,"/"), paste0(FIGDIR, SAMPLE, "/K",k,"/"))
    # FIGDIRTOP=paste0(FIGDIRSAMPLE,"/",SAMPLE,"_K",k,"_dt_", DENSITY.THRESHOLD,"_")
    OUTDIRSAMPLE=paste0(OUTDIR, "/", SAMPLE, "/K",k, "/threshold_", DENSITY.THRESHOLD)

    ## if (SEP) {
    ##     guideCounts <- loadGuides(n, sep=T) %>% mutate(Gene=Gene.marked)
    ##     tmp.labels <- guideCounts$Gene %>% unique() %>% strsplit("-") %>% sapply("[[",2) %>% unique()
    ##     tmp.labels <- tmp.labels[!(tmp.labels %in% c("control","targeting"))]
    ##     rep1.label <- paste0("-",tmp.labels[1])
    ##     rep2.label <- paste0("-",tmp.labels[2])
    ## } else guideCounts <- loadGuides(n) %>% mutate(Gene=Gene.marked)
    
    ## cNMF direct output file (GEP)
    cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.",SUBSCRIPT.SHORT, ".RData")
    ## cNMF.result.file <- paste0(OUTDIRSAMPLE,"/cNMF_results.k_", k, ".dt_", density.threshold, ".RData")
    print(cNMF.result.file)
    if(file.exists(cNMF.result.file)) {
        print("loading cNMF result file")
        load(cNMF.result.file)
    } else {
        print(paste0("file ", cNMF.result.file, " not found"))
    }

    # ## cNMF analysis results file
    # file.name <- ifelse(SEP,
    #                     paste0(OUTDIRSAMPLE,"/cNMFAnalysis.",SUBSCRIPT,".sep.RData"),
    #                     paste0(OUTDIRSAMPLE,"/cNMFAnalysis.",SUBSCRIPT,".RData"))
    # print(file.name)
    # if(file.exists(file.name)) {
    #     print("loading the file")
    #     load(file.name)
    # } else {
    #     print("file not found")
    # }

    ## theta.zscore
    ## theta.zscore.list[[n]] <- theta.zscore %>% as.data.frame %>% `colnames<-`(paste0("factor_", colnames(.))) %>% mutate(Gene=rownames(.)) %>% melt(value.name="weight", id.vars="Gene", variable.name="Factor") %>% mutate(K=k)
    theta.zscore.list[[n]] <- theta.zscore %>% `colnames<-`(paste0("K",k,"_factor_", colnames(.)))
    ## theta.raw
    ## theta.raw.list[[n]] <- theta.raw %>% as.data.frame %>% `colnames<-`(paste0("factor_", colnames(.))) %>% mutate(Gene=rownames(.)) %>% melt(value.name="weight", id.vars="Gene", variable.name="Factor") %>% mutate(K=k)
    theta.raw.list[[n]] <- theta.raw %>% `colnames<-`(paste0("K",k,"_factor_", colnames(.)))
    ## ## theta.KL
    ## ## load KL score
    ## file.name <- paste0(OUTDIRSAMPLE, "/topic.KL.score_", SUBSCRIPT.SHORT, ".txt") %>% gsub("_k_", "_K", .)
    ## if(file.exists(file.name)) {
    ##     print(paste0("Loading ", file.name))
    ##     theta.KL.list[[n]] <- read.table(file.name, header=T, stringsAsFactors=F)
    ## } else {
    ##     print(paste0(file.name, " does not exist."))
    ## }
    
    # ## motif enrichment file (old.211025)
    # file.name <- paste0(OUTDIRSAMPLE,"/cNMFAnalysis.factorMotifEnrichment.",SUBSCRIPT.SHORT, ".RData")
    # print(file.name)
    # if(file.exists(file.name)) {
    #     print("loading the file")
    #     load(file.name)
    # } else {
    #     print("file not found")
    # }
    # if ("all.promoter.fisher.df" %in% ls()) {
    #     promoter.fisher.df.list[[n]] <- all.promoter.fisher.df %>% mutate(K = k)
    #     enhancer.fisher.df.list[[n]] <- all.enhancer.fisher.df %>% mutate(K = k)
    #     rm(list=c("all.promoter.fisher.df", "all.enhancer.fisher.df"))
    # } else {
    #     print("missing all.promoter.fisher.df and/or all.enhancer.fisher.df")
    # }
    
    ## load motif enrichment results
    for(ep.type in c("promoter", "enhancer")){
        num.top.genes <- 300
        file.name <- paste0(OUTDIRSAMPLE, "/", ep.type, ".topic.top.", num.top.genes, ".zscore.gene_motif.count.ttest.enrichment_motif.thr.pval1e-6_", SUBSCRIPT.SHORT,".txt")
        if(file.exists(file.name)) {
            eval(parse(text = paste0("all.", ep.type, ".ttest.df.list[[n]] <- read.delim(file.name, stringsAsFactors=F) %>% mutate(K = k)"))) ## store in all.{promoter, enhancer}.ttest.df.list
        } else {
            message(paste0(file.name, " does not exist"))
        }
    }


    # file.name <- paste0(OUTDIRSAMPLE,"/cNMFAnalysis.factorMotifEnrichment.",SUBSCRIPT.SHORT,".RData")
    # print(file.name)
    # if(file.exists((file.name))) { 
    #     load(file.name)
    #     print(paste0("loading ", file.name))
    # }
    # motif.enrichment.variables <- c("all.enhancer.fisher.df", "all.promoter.fisher.df", 
    #                                 "promoter.wide", "enhancer.wide", "promoter.wide.binary", "enhancer.wide.binary",
    #                                 "enhancer.wide.10en6", "enhancer.wide.binary.10en6", "all.enhancer.fisher.df.10en6",
    #                                 "promoter.wide.10en6", "promoter.wide.binary.10en6", "all.promoter.fisher.df.10en6",
    #                                 "all.promoter.ttest.df", "all.promoter.ttest.df.10en6", "all.enhancer.ttest.df", "all.enhancer.ttest.df.10en6")
    # motif.enrichment.variables.missing <- (!(motif.enrichment.variables %in% ls())) %>% as.numeric %>% sum 
    # if ( motif.enrichment.variables.missing > 0 ) {
    #     warning(paste0(motif.enrichment.variables[!(motif.enrichment.variables %in% ls())], " not available"))
    # } else {
    #     promoter.fisher.df.list[[n]] <- all.promoter.fisher.df %>% mutate(K = k)
    #     enhancer.fisher.df.list[[n]] <- all.enhancer.fisher.df %>% mutate(K = k)
    #     all.promoter.ttest.df.list[[n]] <- all.promoter.ttest.df %>% mutate(K = k)
    #     all.promoter.ttest.df.10en6.list[[n]] <- all.promoter.ttest.df.10en6 %>% mutate(K = k)
    #     all.enhancer.ttest.df.list[[n]] <- all.enhancer.ttest.df %>% mutate(K = k)
    #     all.enhancer.ttest.df.10en6.list[[n]] <- all.enhancer.ttest.df.10en6 %>% mutate(K = k)
    #     all.promoter.fisher.df.list[[n]] <- all.promoter.fisher.df
    #     all.enhancer.fisher.df.list[[n]] <- all.enhancer.fisher.df
    # }


    ## GSEA results
    ranking.types <- c("zscore", "raw", "median_spectra", "median_spectra_zscore")
    for (j in 1:length(GSEA.types)) {
        GSEA.type <- GSEA.types[j]
        to.eval <- paste0("clusterProfiler.", GSEA.type, ".list.here <- vector(\"list\",length(ranking.types))")
        eval(parse(text = to.eval))
        for (i in 1:length(ranking.types)) {
            ranking.type <- ranking.types[i]
            file.name <- paste0(OUTDIRSAMPLE,"/clusterProfiler_GeneRankingType",ranking.type,"_EnrichmentType", GSEA.type, ".txt")
            if(file.exists(file.name)) {
                message("Loading ", file.name)
                to.eval <- paste0("clusterProfiler.", GSEA.type, ".list.here[[i]] <- read.delim(file.name, header=T, stringsAsFactors = F) %>% mutate(type = ranking.type, K = k)")
                eval(parse(text = to.eval))
            } else {
                warning(paste0(file.name, " file does not exist"))
            }
        }
        to.eval <- paste0("clusterProfiler.", GSEA.type, ".list[[n]] <- do.call(rbind, clusterProfiler.", GSEA.type, ".list.here)")
        eval(parse(text = to.eval))
    }


    ## variance explained by the model
    file.name <- paste0(OUTDIRSAMPLE, "/summary.varianceExplained.df.txt")
    if(file.exists(file.name)) {
        varianceExplainedByModel.list[[n]] <- read.delim(file.name, stringsAsFactors=F) %>% mutate(K = k)
    } else {
        message(paste0(file.name, " does not exist"))
    }
    
    ## variance explained per program
    file.name <- paste0(OUTDIRSAMPLE, "/metrics.varianceExplained.df.txt")
    if(file.exists(file.name)) {
        varianceExplainedPerProgram.list[[n]] <- read.delim(file.name, stringsAsFactors=F)
    } else {
        message(paste0(file.name, " does not exist"))
    }        
    
    # ## all statistical tests
    # file.name <- paste0(OUTDIRSAMPLE, "/all.test.", SUBSCRIPT, ".txt")
    # print(file.name)
    # if(file.exists(file.name)) {
    #     print(paste0("loading ", file.name))
    #     all.test.df.list[[n]] <- read.table(file.name, header=T, stringsAsFactors=F) %>% mutate(K = k)
    # } else {
    #     print("all.test file not found")
    # }

    # file.name <- paste0(OUTDIRSAMPLE, "/all.expressed.genes.pval.fdr.", SUBSCRIPT, ".txt")
    # print(file.name)
    # if(file.exists(file.name)) {
    #     print(paste0("loading ", file.name))
    #     all.fdr.df.list[[n]] <- read.table(file.name, header=T, stringsAsFactors=F) %>% mutate(K = k)
    # } else {
    #     print("file not found")
    # }

    # # load count.by.GWAS 
    # file.name <- paste0(OUTDIRSAMPLE,"/count.by.GWAS.classes_p.adj.",p.value.thr %>% as.character,"_",SUBSCRIPT,".txt")
    # print(file.name)
    # if(file.exists(file.name)) {
    #     print(paste0("loading ", file.name))
    #     count.by.GWAS.list[[n]] <- read.delim(file=file.name, header=T, stringsAsFactors=F) %>% mutate(K = k)
    # } else {
    #     print("file not found")
    # }

    # # load count.by.GWAS.with.topic
    # file.name <- paste0(OUTDIRSAMPLE,"/count.by.GWAS.classes.withTopic_p.adj.",p.value.thr %>% as.character,"_",SUBSCRIPT,".txt")
    # print(file.name)
    # if(file.exists(file.name)) {
    #     print(paste0("loading ", file.name))
    #     count.by.GWAS.withTopic.list[[n]] <- read.delim(file=file.name, header=T, stringsAsFactors=F) %>% mutate(K = k)
    # } else {
    #     print("file not found")
    # }
    
    
}

# promoter.fisher.df <- do.call(rbind, promoter.fisher.df.list)
# enhancer.fisher.df <- do.call(rbind, enhancer.fisher.df.list)
GSEA.types <- c("GOEnrichment", "ByWeightGSEA", "GSEA") ## use external input
clusterProfiler.GO.list.here <- clusterProfiler.GSEA.list.here <- clusterProfiler.enricher.GSEA.list.here <- vector("list",length(GSEA.types))
for (j in 1:length(GSEA.types)) {
    GSEA.type <- GSEA.types[j]
    to.eval <- paste0("clusterProfiler.", GSEA.type, ".df <- do.call(rbind, clusterProfiler.", GSEA.type, ".list)")
    eval(parse(text = to.eval))
}
## clusterProfiler.GO.df <- do.call(rbind, clusterProfiler.GO.list)
## clusterProfiler.GSEA.df <- do.call(rbind, clusterProfiler.GSEA.list)
## clusterProfiler.enricher.GSEA.df <- do.call(rbind, clusterProfiler.enricher.GSEA.list)
# all.test.df <- do.call(rbind, all.test.df.list)
# all.fdr.df <- do.call(rbind, all.fdr.df.list)
# count.by.GWAS <- do.call(rbind, count.by.GWAS.list)
# count.by.GWAS.withTopic <- do.call(rbind, count.by.GWAS.withTopic.list)
theta.zscore.df <- do.call(cbind, theta.zscore.list)
theta.raw.df <- do.call(cbind, theta.raw.list)
# theta.KL.df <- do.call(rbind, theta.KL.list)
all.promoter.ttest.df <- do.call(rbind, all.promoter.ttest.df.list)
# all.promoter.ttest.df.10en6 <- do.call(rbind, all.promoter.ttest.df.10en6.list)
all.enhancer.ttest.df <- do.call(rbind, all.enhancer.ttest.df.list)
# all.enhancer.ttest.df.10en6 <- do.call(rbind, all.enhancer.ttest.df.10en6.list)
varianceExplainedByModel.df <- do.call(rbind, varianceExplainedByModel.list)
varianceExplainedPerProgram.df <- do.call(rbind, varianceExplainedPerProgram.list)

file.name <- paste0(OUTDIR.ACROSS.K, "/aggregated.outputs.findK.RData")
save(clusterProfiler.GOEnrichment.df, clusterProfiler.ByWeightGSEA.df, clusterProfiler.GSEA.df, theta.zscore.df, theta.raw.df, all.promoter.ttest.df, all.enhancer.ttest.df, varianceExplainedByModel.df, varianceExplainedPerProgram.df,
     file=file.name)

