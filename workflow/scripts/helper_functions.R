## gene mapping function
map.ENSGID.SYMBOL <- function(topFeatures) {
    gene.type <- ifelse(sum(as.numeric(colnames(topFeatures) %in% "Gene")) > 0,
                        ##(median.spectra.zscore.df) == sum(as.numeric(grepl("^ENS", median.spectra.zscore %>% rownames))),
                        "Gene",
                        "ENSGID")
    db <- ifelse(grepl("mouse", SAMPLE), "org.Mm.eg.db", "org.Hs.eg.db")
    library(!!db)
    if(gene.type == "Gene") {
        if (nrow(topFeatures) == sum(as.numeric(grepl("^ENS", topFeatures$Gene)))) {
            ## put median spectra zscore into ENSGID format for PoPS
            mapped.genes <- mapIds(get(db), keys=topFeatures$Gene, keytype = "ENSEMBL", column = "SYMBOL")
            topFeatures <- topFeatures %>%
                mutate(ENSGID = Gene) %>%
                as.data.frame %>%
                mutate(Gene = mapped.genes)
            na.index <- which(is.na(topFeatures$Gene))
            if(length(na.index) > 0) topFeatures$Gene[na.index] <- topFeatures$ENSGID[na.index]
        } else {
            mapped.genes <- mapIds(get(db), keys=topFeatures$Gene, keytype = "SYMBOL", column = "ENSEMBL")
            topFeatures <- topFeatures %>%
                as.data.frame %>%
                mutate(ENSGID = mapped.genes)
        }
    } else {
        mapped.genes <- mapIds(get(db), keys=topFeatures$ENSGID, keytype = "ENSEMBL", column = "SYMBOL")
        topFeatures <- topFeatures %>%
            as.data.frame %>%
            mutate(Gene = mapped.genes)
        
    }
    return(topFeatures)
}
