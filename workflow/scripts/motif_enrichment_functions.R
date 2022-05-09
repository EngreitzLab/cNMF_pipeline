ttest.on.motifs <- function(wide){##todo:210812:a
    motif.ary <- colnames(wide)[!grepl("gene|topic|<NA>", colnames(wide))]
    all.ttest.df <- do.call(rbind, lapply(motif.ary, function(motif.here) {
        topic.totest <- colnames(wide)[grep("topic", colnames(wide))]
        t.totest <- gsub("topic_", "", topic.totest) %>% sort
        all.topic.ttest.df <- do.call(rbind,lapply(t.totest, function(t) {
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


## ttest.on.motifs <- function(wide, background="all"){##todo:210812:a
##     motif.ary <- colnames(wide)[!grepl("gene|topic|<NA>", colnames(wide))]
##     all.ttest.df <- do.call(rbind, lapply(motif.ary, function(motif.here) {
##         topic.totest <- colnames(wide)[grep("topic", colnames(wide))]
##         t.totest <- gsub("topic_", "", topic.totest) %>% sort
##         all.topic.ttest.df <- do.call(rbind,lapply(t.totest, function(t) {
##             topic <- paste0("topic_",t)
##             wide.here <- wide %>% select(gene, all_of(topic), all_of(motif.here))
##             wide.here[is.na(wide.here)] <- 0
##             top.gene.boolean <- wide.here %>% pull(all_of(topic)) == 1
##             top.gene.index <- top.gene.boolean %>% which
##             if(background == "all") {
##                 background.count <- wide.here[which(!top.gene.boolean),] %>% pull(all_of(motif.here))
##             } else {
##                 background.index <- wide %>%
##                     select(all_of(colnames(wide)[grep("topic", colnames(wide))])) %>% ## select all topic's boolean values
##                     apply(1, function(x) sum(x) > 0) %>%
##                     which %>%
##                     setdiff(x=., y=top.gene.index)
##                 background.count <- wide.here[background.index,] %>% pull(all_of(motif.here))
##             }
##             top.gene.count <- wide.here[which(top.gene.boolean),] %>% pull(all_of(motif.here))
##             test.result <- t.test(top.gene.count, background.count)$p.value
##             output <- data.frame(motif=motif.here,
##                                  topic=topic,
##                                  p.value = test.result) %>%
##                 mutate(global.mean=background.count %>% mean,
##                        global.var=background.count %>% var,
##                        global.count=dim(wide)[1],
##                        top.gene.mean=top.gene.count %>% mean,
##                        top.gene.var=top.gene.count %>% var,
##                        topic.count=top.gene.count %>% length,
##                        enrichment=top.gene.mean/global.mean,
##                        enrichment.log2fc=log2(enrichment)
##                        )                
##         }))
##     })) %>% mutate(p.adjust = p.adjust(p.value, method="BH"))##here210813
##     ## adjust -Inf and +Inf for log2FC
##     inf.index <- all.ttest.df %>% pull(enrichment.log2fc) %>% is.infinite %>% which
##     neg.inf.index <- ( all.ttest.df %>% pull(top.gene.mean) == 0 ) %>% which
##     pos.inf.index <- ( all.ttest.df %>% pull(global.mean) == 0 ) %>% which
##     all.ttest.df$enrichment.log2fc[inf.index] <- 0
##     min.enrichment.log2fc.value <- min(all.ttest.df$enrichment.log2fc)
##     max.enrichment.log2fc.value <- max(all.ttest.df$enrichment.log2fc)
##     all.ttest.df$enrichment.log2fc[neg.inf.index] <- min.enrichment.log2fc.value - 3
##     all.ttest.df$enrichment.log2fc[pos.inf.index] <- max.enrichment.log2fc.value + 3
##     return(all.ttest.df)
## }
