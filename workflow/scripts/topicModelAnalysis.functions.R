######################################################################
## Topic Model Analysis Functions
## Helen Kang
## 210125

######################################################################
## Function for dropping rows with all zeros
drop_zero_rows <- function(data) {
  indices_blank <- as.numeric(which(apply(data,1,max) == 0))
  if(length(indices_blank)!=0){
    data <- as.matrix(data[-indices_blank,]);
  }
  return(data)
}

## Function for plotting heatmap  # new version (adjusted font size)
plotHeatmap <- function(mtx, labCol, title, margins=c(12,6), ...) { #original
  heatmap.2(
    mtx %>% t(), 
    Rowv=T, 
    Colv=T,
    trace='none',
    key=T,
    col=palette,
    labCol=labCol,
    margins=margins, 
    cex.main=0.8, 
    cexCol=4.8/sqrt(nrow(mtx)), cexRow=1.5/(ncol(mtx)^(1/3)),
    main=title,
    ...
  )
}



plotHeatmap <- function(mtx, cellNote=NULL, labCol, title, margins=c(12,6), SNP.annotation = T, ...) { # add annotation
  if(SNP.annotation) mtx <- add.snp.gene.info(mtx, type="rownames")
  if (!is.null(cellNote)){ 
    heatmap.2(
      mtx %>% t(),
      cellnote=cellNote %>% t(),
      notecol="black",
      notecex=0.7,
      Rowv=T,
      Colv=T,
      trace='none',
      key=T,
      col=palette,
      labCol=labCol,
      margins=margins, 
      cex.main=0.8, 
      cexCol=4.8/sqrt(nrow(mtx)), cexRow=1.5/(ncol(mtx)^(1/3)),
      main=title,
      ...
    )
  } else {
    heatmap.2(
      mtx %>% t(), 
      Rowv=T, 
      Colv=T,
      trace='none',
      key=T,
      col=palette,
      labCol=labCol,
      margins=margins, 
      cex.main=0.8, 
      cexCol=4.8/sqrt(nrow(mtx)), cexRow=1.5/(ncol(mtx)^(1/3)),
      main=title,
      ...
    )
  }
}

plotHeatmap.Vertical <- function(mtx,  title) {
  heatmap.2(
    mtx, 
    Rowv=T, 
    Colv=T,
    trace='none',
    key=T,
    col=palette,
    margins=c(6,12), 
    cex.main=0.8, 
    cexRow=3/sqrt(nrow(mtx)), cexCol=1.5/(ncol(mtx)^(1/3)),
    main=title
  )
}

# heatmap but with custom dendrogram # used for all wilcoxon heatmaps with custom dendrograms
plotHeatmapDen.old <- function(mtx, gene.den, topic.den, labCol, title, ...) { 
  heatmap.2(
    mtx %>% t(), 
    Rowv=gene.den, 
    Colv=topic.den,
    trace='none',
    key=T,
    col=palette,
    labCol=labCol,
    margins=c(12,6), 
    cex.main=0.8, 
    cexCol=4.8/sqrt(nrow(mtx)), cexRow=1.5/(ncol(mtx)^(1/3)),
    main=title,
    ...
  )
}

plotHeatmapDen <- function(mtx, gene.den, topic.den, labCol, title, ...) { 
  heatmap.2(
    mtx %>% t(), 
    Rowv=topic.den, 
    Colv=gene.den,
    trace='none',
    key=T,
    col=palette,
    labCol=labCol,
    margins=c(12,6), 
    cex.main=0.8, 
    cexCol=4.8/sqrt(nrow(mtx)), cexRow=1.5/(ncol(mtx)^(1/3)),
    main=title,
    ...
  )
}

## Function for plotting UMAP
plotUmap <- function(s, title, split=NULL, ...) {
  pdf(...)
  print(ElbowPlot(s))
  Idents(s) <- s$seurat_clusters
  DimPlot(s, reduction = "umap", label=TRUE, split.by = split) %>% print()
  DimPlot(s, reduction = "umap", label=FALSE, split.by = split) %>% print()
  
  FeaturePlot(s, reduction = "umap", features='nCount_RNA', split.by = split) %>% print()
  FeaturePlot(s, reduction = "umap", features='nFeature_RNA', split.by = split) %>% print()
  
  if ( "percent.mt" %in% colnames(s[[]]) ) FeaturePlot(s, reduction = "umap", features='percent.mt', split.by = split) %>% print()
  if ("percent.ribo" %in% colnames(s[[]])) FeaturePlot(s, reduction = "umap", features='percent.ribo', split.by = split) %>% print()

  if("plot.features" %in% ls()){
      meta.data.names <- colnames(s[[]])
      for (feature.name in plot.features) {
          feature.vec <- s@meta.data %>% select(feature.name)
          if(feature.vec[1,1] %>% is.numeric()) { # numeric
              FeaturePlot(s, reduction = "umap", features=feature.name, split.by = split) %>% print()
          } else { # discrete values / categories
              Idents(s) <- s@meta.data %>% select(feature.name)
              DimPlot(s, reduction = "umap", label=TRUE, split.by = split) %>% print()
          }
      }
  }
  
  if("omega" %in% ls()) {
      for (i in colnames(omega)) {
          FeaturePlot(s, reduction = "umap", features = i, split.by = split) %>% print()
      }
  }
  dev.off()
  
}   



## Function for plotting UMAP
plotUmap.geneSet <- function(s, title, split=NULL, file=NULL, gene.set=NULL, plot.features=NULL, ...) {
  pdf(file, ...)
  print(ElbowPlot(s))
  Idents(s) <- s$seurat_clusters
  DimPlot(s, reduction = "umap", label=TRUE, split.by = split) %>% print()
  
  FeaturePlot(s, reduction = "umap", features='nCount_RNA', split.by = split) %>% print()
  FeaturePlot(s, reduction = "umap", features='nFeature_RNA', split.by = split) %>% print()

  if ( "percent.mt" %in% colnames(s[[]]) ) FeaturePlot(s, reduction = "umap", features='percent.mt', split.by = split) %>% print()
  if ("percent.ribo" %in% colnames(s[[]])) FeaturePlot(s, reduction = "umap", features='percent.ribo', split.by = split) %>% print()
  
  if("plot.features" %in% ls() & length(plot.features) > 0){
      meta.data.names <- colnames(s[[]])
      for (feature.name in plot.features) {
          feature.vec <- s@meta.data %>% select(feature.name)
          if(feature.vec[1,1] %>% is.numeric()) { # numeric
              FeaturePlot(s, reduction = "umap", features=feature.name, split.by = split) %>% print()
          } else { # discrete values / categories
              Idents(s) <- s@meta.data %>% select(feature.name)
              DimPlot(s, reduction = "umap", label=TRUE, split.by = split) %>% print()
          }
      }
  }

  if("omega" %in% ls()) {
      for (i in colnames(omega)) {
          FeaturePlot(s, reduction = "umap", features = i, split.by = split) %>% print()
      }
  }

  if("gene.set" %in% ls() & length(gene.set)>0) {
      all.labels <- rownames(s)
      gene.set.label <- gene.set[which(gene.set %in% all.labels)]
      
      if(exists("gene.set.label") & length(gene.set.label) > 0) {
          for(gene in gene.set.label) {
            print(gene)
            FeaturePlot(s, reduction = "umap", features=gene) %>% print()
          }
      }
  }
  
  dev.off()
  
}   


rsq <- function (x, y) cor(x, y) ^ 2 # R^2


# for drawing scatter plot with overlaying heatmap [ to be used in ggpairs ]
my_bin <- function(data, mapping, ..., low = "#132B43", high = "#FDB515") {
  ggplot(data = data, mapping = mapping) +
    geom_bin2d(...) +
    theme(axis.text = element_text(size = 0.1)) + 
    scale_fill_gradient(low = low, high = high)
}



# add enhancer-cad-snp linked gene info to guideCounts or rownames
add.snp.gene.info <- function(df, type="df"){
  if(type=="df"){
    tmp <- df %>% mutate(Gene.marked = Gene)
    for (i in 1:nrow(enh.snp.to.gene)){
      tmp$Gene.marked <- gsub(enh.snp.to.gene$Enhancer_name[i], enh.snp.to.gene$alternative_name_w_gene_symbols[i],  tmp$Gene.marked)
    }
  } else {
    tmp <- df
    for (i in 1:nrow(enh.snp.to.gene)){
      rownames(tmp) <- gsub(enh.snp.to.gene$Enhancer_name[i], enh.snp.to.gene$alternative_name_w_gene_symbols[i],  rownames(tmp))
    }
  }
  return(tmp) #HERE
}



# load guides
loadGuides <- function(n, sep=T, OUTDIR="/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/topicModel/2010_Good_singlet/") {
  ### INPUT:
  # n: sample index (numeric)
  # sep: whether to separate -no and -plus in "pooled" datset (boolean)
  #. 
  # load guide info
  if ( grepl("Telo",SAMPLE[n], fixed=TRUE) ) {
    guideCounts <- read.table(paste0(DATADIR, SAMPLE[n], "/outs/filtered_feature_bc_matrix/noguidemut_calls.tsv"), header=T)
    UMICounts <- read.table(paste0(OUTDIR, SAMPLE[n], "/UMICounts.txt"), header=T)
    guideCounts <- merge(guideCounts, UMICounts, by="CBC")
    guideCounts$CBC <- paste0(guideCounts$CBC,"-1")
    # add Gene.marked
    tmp <- guideCounts %>% mutate(Gene.marked = Gene)
    for (i in 1:nrow(enh.snp.to.gene)){
      tmp$Gene.marked <- gsub(enh.snp.to.gene$Enhancer_name[i], enh.snp.to.gene$alternative_name_w_gene_symbols[i],  tmp$Gene.marked)
    }
    guideCounts <- tmp
  } else if ( grepl("IL1B",SAMPLE[n], fixed=TRUE) ) {
    guideCounts.list <- vector("list", 2)
    for ( i in 1:length(guideCounts.list) ) {
      guideCounts <- read.table(paste0(DATADIR, "Telo_", ifelse(grepl("no_IL1B", SAMPLE[n]), "no_IL1B", "plus_IL1B"), "_T200_", i, "/outs/filtered_feature_bc_matrix/noguidemut_calls.tsv"), header=T)
      UMICounts <- read.table(paste0(OUTDIR, "Telo_", ifelse(grepl("no_IL1B", SAMPLE[n]), "no_IL1B", "plus_IL1B"), "_T200_", i, "/UMICounts.txt"), header=T)
      guideCounts <- merge(guideCounts, UMICounts, by="CBC")
      guideCounts$CBC <- paste0(guideCounts$CBC,"-", i)
      guideCounts.list[[i]] <- guideCounts
    }
    guideCounts <- do.call(rbind, guideCounts.list)
    # add Gene.marked
    tmp <- guideCounts %>% mutate(Gene.marked = Gene)
    for (i in 1:nrow(enh.snp.to.gene)){
      tmp$Gene.marked <- gsub(enh.snp.to.gene$Enhancer_name[i], enh.snp.to.gene$alternative_name_w_gene_symbols[i],  tmp$Gene.marked)
    }
    guideCounts <- tmp
    
  } else { #pooled
    guideCounts.list <- vector("list", 4)
    for ( i in 1:length(guideCounts.list) ) {
      guideCounts <- read.table(paste0(DATADIR, STATIC.SAMPLE[i], "/outs/filtered_feature_bc_matrix/noguidemut_calls.tsv"), header=T)
      UMICounts <- read.table(paste0(OUTDIR, STATIC.SAMPLE[i], "/UMICounts.txt"), header=T)
      guideCounts <- merge(guideCounts, UMICounts, by="CBC")
      guideCounts$CBC <- paste0(guideCounts$CBC,"-", i)
      # add Gene.marked
      tmp <- guideCounts %>% mutate(Gene.marked = Gene)
      for (i in 1:nrow(enh.snp.to.gene)){
        tmp$Gene.marked <- gsub(enh.snp.to.gene$Enhancer_name[i], enh.snp.to.gene$alternative_name_w_gene_symbols[i],  tmp$Gene.marked)
      }
      guideCounts <- tmp
      if (sep) {
        if (i <= 2) { # no IL1B
          guideCounts$Gene <- paste0(guideCounts$Gene, "-no")
        } else {
          guideCounts$Gene <- paste0(guideCounts$Gene, "-plus")
        }
      }
      guideCounts.list[[i]] <- guideCounts
    }
    guideCounts <- do.call(rbind, guideCounts.list)
    
  }
  guideCounts$Gene.SNP <- guideCounts$Gene
  return(guideCounts)
}



#################################################
## Dendrogram functions #den.fun
# get best dendrogram
getClusteringMethod <- function(score.mtx, verb=F, default=F){
  
  topic.max.corr <- 0
  gene.max.corr <- 0
  topic.best.method <- "complete"
  gene.best.method <- "complete"
  if ( default ) {
    best.topic.den <- as.dendrogram(hclust(dist(score.mtx %>% t()), method = topic.best.method))
    best.gene.den <- as.dendrogram(hclust(dist(score.mtx), method = gene.best.method))
  } else {
    for (method in c("average", "complete", "ward.D")) {
      
      gene.cor <- cor(score.mtx %>% dist(), cophenetic(hclust(score.mtx %>% dist(), method = method))) # "ward.D" minimizes variance within cluster
      topic.cor <- cor(score.mtx %>% t() %>% dist(), cophenetic(hclust(score.mtx %>% t() %>% dist(), method = method))) # need to fix corr for K=2 (only one data point)
      
      topic.best.method <- ifelse(topic.max.corr < topic.cor, method, topic.best.method)
      topic.max.corr <- ifelse(topic.max.corr < topic.cor, topic.cor, topic.max.corr)
      gene.best.method <- ifelse(gene.max.corr < gene.cor, method, gene.best.method)
      gene.max.corr <- ifelse(gene.max.corr < gene.cor, gene.cor, gene.max.corr)
      
    }
    
    if (verb) {
      print(paste0("Gene's best clustering method: ", gene.best.method, ", with corr = ", gene.max.corr))
      print(paste0("Topic's best clustering method: ", topic.best.method, ", with corr = ", topic.max.corr))
    }
    
    best.topic.den <- as.dendrogram(hclust(dist(score.mtx %>% t()), method = topic.best.method))
    best.gene.den <- as.dendrogram(hclust(dist(score.mtx), method = gene.best.method))
  }
  # return best method name, corr, and dendrogram
  result <- list("topic.den" = best.topic.den,
                 "gene.den" = best.gene.den,
                 "topic.method" = topic.best.method,
                 "gene.method" = gene.best.method,
                 "topic.corr" = topic.max.corr,
                 "gene.corr" = gene.max.corr
  )
  return(result)
}

# function to get subgroup pairwise distances
getSubgroupDist <- function(dendrogram, sub.label=(dendrogram %>% labels), method="coph") {
  if (class(sub.label)=="matrix") { # for between group
    if (method == "coph") {
      subset.dist <- cophenetic(dendrogram)
    } else {
      subset.dist <- dist(dendrogram)
    }
    subset.pairs <- sub.label
  } else { # subset based on sub.label array
    if (method == "coph") {
      subset.dist <- dist_subset(cophenetic(dendrogram), sub.label)
    } else {
      subset.dist <- dist_subset(dist(dendrogram), sub.label) # the input should be the raw gene.score matrix
    }
    subset.pairs <- subset.dist %>% as.matrix() %>% colnames() %>% combn(2) %>% t()
  }
  subset.pairs.dist <- data.frame(subset.pairs, dist=(subset.dist %>% as.matrix())[subset.pairs])
  return(subset.pairs.dist)
}

# function to create list: for each ingroup gene, pair with all outgroup gene
getBetweenGroupLabel <- function(dendrogram.distObj, ingroup.label, complementgroup.label){
  between.pairs.label.df <- dendrogram.distObj %>% as.matrix() %>% colnames() %>% combn(2) %>% t() %>% `colnames<-`(c("X1","X2")) %>% as.data.frame() %>% 
    subset(((X1 %in% ingroup.label) & (X2 %in% complementgroup.label)) | ((X2 %in% ingroup.label) & (X1 %in% complementgroup.label))) %>% as.matrix()
  return(between.pairs.label.df)
}


# function to calculate how close a set of genes relate on a dendrogram and by Euclidean distance
test.geneSet <- function(gene.score, gene.set, verb=F, rand=F, default=F) { 
  ## Inputs:
  # gene.score: a df htat has rows of genes and columns of values (e.g. log2fc for topics)
  # gene.set: array of genes to be tested on
  # .
  ## Returns a list with fields: 
  # coph.test
  # eucl.test
  # gene.set.label: genes from gene.set that appeared in gene.score matrix
  # den.list: best dendrogram and its meta for gene.score
  # coph.test.p.value
  # eucl.test.p.value
  # rand.coph.test: select randomly a set of genes of the same size as ingroup
  # rand.eucl.test
  # rand.coph.test.p.value
  # rand.eucl.test.p.value
  # .
  
  if(length(gene.set) <= 2) { return("gene set is too small (need >2 genes)") }
  
  log2fc.den.list <- getClusteringMethod(gene.score, verb=verb, default=default)
  
  ## distances within a known set of genes
  ingroup.label <- gene.set[(gene.set %in% (log2fc.den.list[["gene.den"]] %>% labels)) %>% which()] %>% as.character()# convert these into a function
  if(ingroup.label %>% length() < 3) { return("gene set is too small (need >2 genes in the perturbed set)") }
  ingroup.coph.dist <- getSubgroupDist(log2fc.den.list[["gene.den"]], ingroup.label) # gene.pairs.dist$dist has all unique pairwise distances
  ingroup.eucl.dist <- getSubgroupDist(gene.score, ingroup.label, method="eucl") 
  
  # get complementary group label
  all.labels <- log2fc.den.list[["gene.den"]] %>% labels
  complementgroup.label <- all.labels[which(!((log2fc.den.list[["gene.den"]] %>% labels) %in% gene.set))]
  
  ## calculate {cophenetic, euclidean} distance between one in-group and each of out-group element 
  between.pairs.label.df <- getBetweenGroupLabel(log2fc.den.list[["gene.den"]] %>% cophenetic(), ingroup.label=ingroup.label, complementgroup.label=complementgroup.label)
  betweengroup.coph.dist <- getSubgroupDist(log2fc.den.list[["gene.den"]], sub.label=between.pairs.label.df, method="coph")
  
  between.pairs.label.df <- getBetweenGroupLabel(gene.score %>% dist(), ingroup.label=ingroup.label, complementgroup.label=complementgroup.label)
  betweengroup.eucl.dist <- getSubgroupDist(gene.score %>% dist(), sub.label=between.pairs.label.df, method="eucl")
  
  ## statistical test of the within-group and across-group distances
  coph.test <- t.test(ingroup.coph.dist$dist, betweengroup.coph.dist$dist) # test on cophenetic distance
  eucl.test <- t.test(ingroup.eucl.dist$dist, betweengroup.eucl.dist$dist) # test on euclidean distance
  
  ## Addon: compare to a randomly selected set of genes of the same size as ingroup
  if (!rand){
    output <- list("coph.test"=coph.test,
                   "eucl.test"=eucl.test,
                   "gene.set.label"=ingroup.label,
                   "den.list"=log2fc.den.list,
                   "coph.test.p.value"=coph.test$p.value,
                   "eucl.test.p.value"=eucl.test$p.value #,
                   # "dendrogram.list"=log2fc.den.list,
                   # "topic.den.best.method"=log2fc.den.list[["topic.method"]],
                   # "gene.den.best.method"=log2fc.den.list[["gene.method"]]
    )
  } else {
    gene.names <- rownames(gene.score) #[which(!(rownames(gene.score) %in% ingroup.label))] 
    rand.idx <- runif(length(ingroup.label), min=1, max=length(gene.names)+1) %>% floor
    rand.gene.set <- gene.names[rand.idx]
    
    between.pairs.rand.label.df <- getBetweenGroupLabel(log2fc.den.list[["gene.den"]] %>% cophenetic(), ingroup.label=ingroup.label, complementgroup.label=rand.gene.set)
    betweengroup.rand.coph.dist <- getSubgroupDist(log2fc.den.list[["gene.den"]], sub.label=between.pairs.rand.label.df, method="coph")
    
    between.pairs.rand.label.df <- getBetweenGroupLabel(gene.score %>% dist(), ingroup.label=ingroup.label, complementgroup.label=rand.gene.set)
    betweengroup.rand.eucl.dist <- getSubgroupDist(gene.score %>% dist(), sub.label=between.pairs.rand.label.df, method="eucl")
    
    rand.coph.test <- t.test(ingroup.coph.dist$dist, betweengroup.rand.coph.dist$dist) # test on cophenetic distance
    rand.eucl.test <- t.test(ingroup.eucl.dist$dist, betweengroup.rand.eucl.dist$dist) # test on euclidean distance
    
    output <- list("coph.test"=coph.test,
                   "eucl.test"=eucl.test,
                   "gene.set.label"=ingroup.label,
                   "den.list"=log2fc.den.list,
                   "coph.test.p.value"=coph.test$p.value,
                   "eucl.test.p.value"=eucl.test$p.value,
                   "rand.coph.test"=rand.coph.test,
                   "rand.eucl.test"=rand.eucl.test,
                   "rand.coph.test.p.value"=rand.coph.test$p.value,
                   "rand.eucl.test.p.value"=rand.eucl.test$p.value #,
                   # "dendrogram.list"=log2fc.den.list,
                   # "topic.den.best.method"=log2fc.den.list[["topic.method"]],
                   # "gene.den.best.method"=log2fc.den.list[["gene.method"]]
    )
  }
  return(output)
}

# truncate decimals (used in function merge.score.with.test)
trunc <- function(x, ..., prec = 0) base::trunc(x * 10^prec, ...) / 10^prec

# merge test results with log2fc so that not-significant entries are zero in the matrix

merge.score.with.test <- function(score.mtx, test.list.df, test.col.name, p.value.thr=0.05, adj.p.value.thr=0.05, fill.all=F, fill=0, overlay=F, precision=2, num.thr=0){
  ## INPUT:
  # score.mtx: log2fc matrix, (e.g. gene.score, GO.score)
  # test.list.df: df that has p.value information (e.g. gene.test)
  # test.col.name: column name for the test needed (e.g. wilcox.p, ttest.p)
  #.
  topic.names <- test.list.df$Topic %>% unique()
  category.name <- colnames(test.list.df)[1]
  tmp <- score.mtx #FIXEDFLAG
  tmp[[category.name]] <- rownames(tmp)
  tmp <- tmp %>% melt(id.vars=category.name, variable.name="Topic", value.name="log2fc")
  
  tmpp <- merge(test.list.df %>% as.data.frame() %>% select(all_of(category.name), all_of(test.col.name), all_of(paste0("adjusted.", test.col.name)), Topic), tmp, by=c(category.name, "Topic"))
  tmpp <- tmpp %>% subset(get(test.col.name) < p.value.thr & get(paste0("adjusted.", test.col.name)) < adj.p.value.thr) %>% 
    mutate(log2fc = trunc(log2fc, prec=precision)) 
  
    if(dim(tmpp)[1] > 0 ) {
    tmpp$log2fc[abs(tmpp$log2fc) < num.thr] <- "*"
  if ( overlay ) tmpp$log2fc <- "*" 
  
  if (fill.all & (length(tmpp$Topic %>% unique()) < length(topic.names))) { # include all topics whether or not it has significant perturbation
    tmppp <- tmpp %>% select(-all_of(test.col.name),-all_of(paste0("adjusted.", test.col.name))) 
    no.sig.gene.topics <- topic.names[which(!(topic.names %in% (tmpp$Topic %>% unique())))]
    tmppp <- rbind(tmppp, data.frame(Gene=tmppp$Gene[1],#HERE
                                     Topic=no.sig.gene.topics,
                                     log2fc=fill)) %>% spread(key=Topic, value=log2fc, fill=fill)
  } else tmppp <- tmpp %>% select(-all_of(test.col.name),-all_of(paste0("adjusted.", test.col.name))) %>% spread(key=Topic, value=log2fc, fill=fill)
  toUse.topic.names <- topic.names[topic.names %in% colnames(tmppp)]  # in case tmppp has less topics than those in topic.names
  tmppp <- tmppp[, c(category.name, toUse.topic.names)] # order the topic names from topic_1 to topic_k
  
  rownames(tmppp) <- tmppp[[category.name]]
  result <- list("score.mtx" = tmppp %>% select(-all_of(category.name)) %>% as.matrix(),
                 "p" = tmpp,
                 "p.value.thr" = p.value.thr,
                 "adj.p.value.thr" = adj.p.value.thr)
    } else {
        warning("No perturbation has a significant topic")
        result=0
    }
    return(result)
}

## quick way to get K.list
get.K.list <- function(sample) {
  if ( sample == "pooled" ) K.list = c(3,4,5,6,7,8,9,10,11,12,13,14,15,17,21,23,25) %>% as.numeric() else K.list=c(3,4,5,6,7,8,9,10,11,12,13,14,15,17,19,21,23,25) %>% as.numeric()
  return(K.list)
}


## function for customized qqplot
my.qqplot <- function (x, y, match=F, match.x=F, match.y=F, match.df, # only input match.df when match=T
                       plot.it = TRUE, xlab = deparse(substitute(x)),
                       ylab = deparse(substitute(y)), main, xlimit, ylimit, ...) 
{
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny > lenx) # expand the smaller data set
    sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny < lenx) 
    sy <- approx(1L:leny, sy, n = lenx)$y
  if (match) { 
    if (match.x) {
      toPlot.here <- data.frame(x=sx, y=sy, highlight=ifelse(sx %in% match.df$p.value, "p.adjust.significant", "p.adjust.not.significant")) 
    } else {
      toPlot.here <- data.frame(x=sx, y=sy, highlight=ifelse(sy %in% match.df$p.value, "p.adjust.significant", "p.adjust.not.significant")) 
    }
  } else toPlot.here <- data.frame(x=sx, y=sy)
  
  if (plot.it) # add color
    if (match){
      p <- toPlot.here %>% ggplot(aes(x,y,color=highlight)) + geom_point(shape=1) + xlab(xlab) + ylab(ylab) + mytheme + ggtitle(main) +
        xlim(xlimit) + ylim(ylimit) + coord_fixed() + scale_color_manual(values=c("black", "red")) + geom_abline()
      invisible(list(x = sx, y = sy))
      print(p)
    }
}


## make dotplot of log2FC of perturbation and ctrl side by side, for per guide and per gene
## make CDF of log2FC of perturbation and ctrl, for per guide and per gene
# requires topicModelAnalysis results to be loaded
plot.CDF.dotplot <- function(df) {
  # this function needs gene.fc.ann.omega.ctrl to be pre-computed
  this.gene <- df$Gene[1] # the input ptbd gene
  this.df <- rbind(df, gene.fc.ann.omega.ctrl)
  
  # prepare for p-values to be displayed in the title
  INT.ttest.this.gene <- wilcox.this.gene <- vector("list",2)
  INT.ttest.this.gene[[1]] <- INT.gene.test %>% subset(Gene==this.gene)
  wilcox.this.gene[[1]] <- gene.test %>% subset(Gene==this.gene)
  INT.ttest.this.gene[[2]] <- INT.avg.guide.test %>% subset(Gene==this.gene)
  wilcox.this.gene[[2]] <- perguide.gene.test %>% subset(Gene==this.gene)
  
  ## utility functions
  get.title <- function(i, per="cell") {
    topic <- paste0("topic_", i)
    if (per=="cell") index <- 1 else index <- 2
    p.value.INT <- INT.ttest.this.gene[[index]] %>% subset(Topic==topic) %>% select(Gene, ttest.p, adjusted.ttest.p) %>% 
      `colnames<-`(c("Gene", "p.value", "adjusted.p.value")) %>% mutate(test.type="INT.ttest")
    p.value.wilcox <- wilcox.this.gene[[index]] %>% subset(Topic==topic) %>% select(Gene, wilcox.p, adjusted.wilcox.p) %>% 
      `colnames<-`(c("Gene", "p.value", "adjusted.p.value")) %>% mutate(test.type="wilcoxon")
    p.value <- rbind(p.value.INT, p.value.wilcox)
    p.value.str <- paste0(p.value$test.type, " p-value: ", format(p.value$p.value, digits=4), 
                          " adjusted p-value: ", format(p.value$adjusted.p.value, digits=4)) %>% unlist() %>% paste0(collapse="\n")
    return(p.value.str)
  }
  
  
  ## section 1: Dot plot [per cell {perturbation, ctrl}]
  toPlot <- this.df %>% select(-Guide, -Row.names) %>% melt(id.vars="Gene", variable.name="topic", value.name="fc.contribution")
  toPlot$fc.contribution <- unlist(toPlot$fc.contribution) # so that geom_dotplot would work
  toPlot$Gene <- toPlot$Gene %>% as.character()
  toPlot$Gene[grepl("^safe-targeting|^negative-control",toPlot$Gene)] <- "ctrl" # combine ctrl samples
  toPlot$log2fc.contribution <- log2(toPlot$fc.contribution)
  toPlot1 <- toPlot # section 1 and 2
  
  
  ## section 3: Dot plot [per guide with error bar {perturbation, ctrl}]
  # add standard error on guide dotplots
  calc.mean.stderr <- function(df) {
    # df: 1 guide (rows = cells) x topics (columns = k topics)
    mtx <- df %>% select(grep("topic",colnames(df)))
    df.mean <- apply(mtx,2,mean)
    df.stderr <- apply(mtx,2,function(x) sd(x)/sqrt(nrow(df))*1.97)
    df.stderr[is.na(df.stderr)] <- 0
    out <- data.frame(Guide = rep(df$Guide[1],ncol(mtx)),
                      Gene = rep(df$Gene[1],ncol(mtx)),
                      topic = colnames(mtx),
                      fc.contribution = df.mean,
                      fc.contribution.stderr = df.stderr,
                      sample.size=nrow(df)
    )
  }
  toPlot <- this.df %>% group_by(Guide) %>% do(calc.mean.stderr(.)) %>% as.data.frame()
  toPlot <- toPlot %>% mutate(log2fc.contribution=log2(fc.contribution), 
                              log2fc.contribution.stderr=0.5*log2((fc.contribution + fc.contribution.stderr)/(fc.contribution - fc.contribution.stderr)))
  type <- "guide"
  toPlot$Gene <- toPlot$Gene %>% as.character()
  toPlot$Guide <- toPlot$Guide %>% as.character()
  toPlot$Gene[grepl("^safe-targeting|^negative-control",toPlot$Gene)] <- "ctrl"
  # toPlot$log2fc.contribution <- unlist(toPlot$log2fc.contribution) # so that geom_dotplot would work
  # toPlot$log2fc.contribution.stderr <- unlist(toPlot$log2fc.contribution.stderr) # so that geom_dotplot would work
  toPlot2 <- toPlot # section 3 and 4
  for (i in 1:k) {
    p.value.str <- get.title(i,"cell")
    topic <- paste0("topic_", i)
    ## section 1
    title <- paste0(SAMPLE[n], ", ", this.gene, ", K = ", k, ", topic ", i, ", per Cell\n", p.value.str)
    p1 <- toPlot1[toPlot1$topic==topic,] %>% ggplot(aes(x=Gene, y=fc.contribution)) +
      geom_point(size=0.5, position="jitter", color="red") + 
      geom_hline(aes(yintercept = 1),  linetype = "dashed") +
      ggtitle(title) + mytheme
    p1.log2 <- toPlot1[toPlot1$topic==topic,] %>% ggplot(aes(x=Gene, y=log2fc.contribution)) +
      geom_point(size=0.5, position="jitter", color="red") + 
      geom_hline(aes(yintercept = 0),  linetype = "dashed") +
      ggtitle(title) + mytheme
    ## section 2 CDF [per cell {perturbation, ctrl}]
    p2 <- toPlot1[toPlot1$topic==topic,] %>% ggplot(aes(log2fc.contribution, color=Gene)) + stat_ecdf() +
      mytheme + ylab("Fraction of cells") + ggtitle(paste0(SAMPLE[n], ", ", this.gene, ", topic ", i, ", per cell"))
    # print(p)
    ## section 3 dot plot per guide with stderr
    p.value.str <- get.title(i,"guide")
    title <- paste0(SAMPLE[n], ", ", this.gene, ", K = ", k, ", topic ", i, ", per Guide\n", p.value.str)
    p3 <- toPlot2[toPlot2$topic==topic,] %>% ggplot(aes(x=Gene, y=fc.contribution, size=sample.size)) +
      geom_pointrange(mapping = aes(x = Gene,
                                    y = fc.contribution,
                                    ymin = fc.contribution-fc.contribution.stderr,
                                    ymax = fc.contribution+fc.contribution.stderr),
                      size=0.15, color="red", position=position_jitter()) +
      geom_hline(aes(yintercept = 1),  linetype = "dashed") +
      ggtitle(title) + mytheme
    p3.log2 <- toPlot2[toPlot2$topic==topic,] %>% ggplot(aes(x=Gene, y=log2fc.contribution, size=sample.size)) +
      geom_pointrange(mapping = aes(x = Gene,
                                    y = log2fc.contribution,
                                    ymin = log2(fc.contribution-fc.contribution.stderr),
                                    ymax = log2(fc.contribution+fc.contribution.stderr)),
                      size=0.15, color="red", position=position_jitter()) +
      geom_hline(aes(yintercept = 0),  linetype = "dashed") +
      ggtitle(title) + mytheme
    # print(p)
    ## section 3: CDF [per averaged guide value {perturbation, ctrl}]
    p4 <- toPlot2[toPlot2$topic==topic,] %>% ggplot(aes(log2fc.contribution, color=Gene)) + stat_ecdf() +
      mytheme + ylab("Fraction of cells") + ggtitle(paste0(SAMPLE[n], ", ", this.gene, ", topic ", i, "\naverage log2FC for each guide"))
    # print(p)
    
    p <- ggarrange(p1,p3,p1.log2,p3.log2,p2,p4,nrow=3,ncol=2)
    print(p)
  }
  
}



## calculate empirical fdr
fdr = function(pval, realPvals, nullPvals) min(1, (sum(nullPvals<=pval)/length(nullPvals)) / (sum(realPvals<=pval)/length(realPvals)))


## initial QC plot before filtering cells before running topic model
plotSingleCellStats <- function(s.meta, mtMax=NULL, nCountMax=NULL) {
    addMaxLine <- function(maxVal=NULL) if (is.null(maxVal)) mytheme else geom_vline(xintercept=maxVal, col='red', linetype='dashed')
    print(ggplot(s.meta, aes(x=percent.mt)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,100) + xlab("% mitochondrial counts") + mytheme + addMaxLine(mtMax))
    # print(ggplot(s.meta, aes(x=percent.mt, color=MergedGuideSet)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,100) + xlab("% mitochondrial counts") + mytheme + addMaxLine(mtMax))
    # print(ggplot(s.meta, aes(x=percent.mt, color=GuideCall)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,100) + xlab("% mitochondrial counts") + mytheme + addMaxLine(mtMax))
    
    print(ggplot(s.meta, aes(x=nCount_RNA)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,NA) + xlab("# UMIs per cell") + mytheme + addMaxLine(nCountMax))
    print(ggplot(s.meta, aes(x=nCount_RNA)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,NA) + xlab("# UMIs per cell") + mytheme + addMaxLine(nCountMax)) + xlim(0, 30000)
    # print(ggplot(s.meta, aes(x=nCount_RNA, color=MergedGuideSet)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,NA) + xlab("# UMIs per cell") + mytheme + addMaxLine(nCountMax))
    # print(ggplot(s.meta, aes(x=nCount_RNA, color=GuideCall)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,NA) + xlab("# UMIs per cell") + mytheme + addMaxLine(nCountMax))
    
    print(ggplot(s.meta, aes(x=nFeature_RNA)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,NA) + xlab("# unique genes per cell") + mytheme)
    # print(ggplot(s.meta, aes(x=nFeature_RNA, color=MergedGuideSet)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,NA) + xlab("# unique genes per cell") + mytheme)
    # print(ggplot(s.meta, aes(x=nFeature_RNA, color=GuideCall)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,NA) + xlab("# unique genes per cell") + mytheme)
    
    
    print(ggplot(s.meta, aes(x=percent.ribo)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,100) + xlab("% ribosomal gene counts") + mytheme)
    # print(ggplot(s.meta, aes(x=percent.ribo, color=MergedGuideSet)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,100) + xlab("% ribosomal gene counts") + mytheme)
    # print(ggplot(s.meta, aes(x=percent.ribo, color=GuideCall)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,100) + xlab("% ribosomal gene counts") + mytheme)
  }




## Function for sweep normalize every cell's topic expression by control topic expression
sweep.wrapper <- function(mat, vec, diff=3) {
    output <- sweep(mat, 2, vec, `/`) %>% log2()
    output[!is.finite(output)] <- floor(min(output[is.finite(output)]) - diff)
    output <- output %>% as.data.frame()
    return(output)
}


## plot single cell metrics for QC
plotSingleCellStats <- function(s.meta, mtMax=NULL, nCountMax=NULL, ...) {
    pdf(...)
    addMaxLine <- function(maxVal=NULL) if (is.null(maxVal)) mytheme else geom_vline(xintercept=maxVal, col='red', linetype='dashed')
    print(ggplot(s.meta, aes(x=percent.mt)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,100) + xlab("% mitochondrial counts") + mytheme + addMaxLine(mtMax))
    print(ggplot(s.meta, aes(x=nCount_RNA)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,NA) + xlab("# UMIs per cell") + mytheme + addMaxLine(nCountMax))
    print(ggplot(s.meta, aes(x=nFeature_RNA)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,NA) + xlab("# unique genes per cell") + mytheme)
    print(ggplot(s.meta, aes(x=percent.ribo)) + stat_ecdf() + ylab("Cumulative fraction of cells") + xlim(0,100) + xlab("% ribosomal gene counts") + mytheme)
    dev.off()
}


## calculate UMAP in Seurat
calcUMAP <- function(sg) {
    sg <- SCTransform(sg)
    sg <- RunPCA(sg, verbose = FALSE)
    sg <- FindNeighbors(sg, dims = 1:10)
    sg <- FindClusters(sg, resolution = 0.06)
    sg <- RunUMAP(sg, dims = 1:10)
}

