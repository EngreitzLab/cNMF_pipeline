
##########################################################################################
## create table
create_topic_definition_table <- function(theta.zscore, t) {
    out <- theta.zscore[,t] %>%
        as.data.frame %>%
        `colnames<-`(c("zscore")) %>%
        mutate(Perturbation = rownames(theta.zscore), .before="zscore") %>%
        merge(gene.summary, by.x="Perturbation", by.y="Gene", all.x=T) %>%
        arrange(desc(zscore)) %>%
        mutate(Rank = 1:n(), .before="Perturbation") %>%
        mutate(ProgramID = paste0("K", k, "_", t), .before="zscore") %>%
        arrange(Rank) %>%
        mutate(My_summary = "", .after = "zscore") %>%
        select(Rank, ProgramID, Perturbation, zscore, My_summary, FullName, Summary)
}

create_topic_regulator_table <- function(all.test, program.here, fdr.thr = 0.1) {
    out <- MAST.df %>%
        subset(ProgramID == program.here &
               fdr.across.ptb < fdr.thr) %>%
        select(Perturbation, fdr.across.ptb, log2FC, log2FC.ci.hi, log2FC.ci.lo, fdr, p.value) %>%
        merge(gene.summary, by.x="Perturbation", by.y="Gene", all.x=T) %>%
        arrange(fdr.across.ptb, desc(log2FC)) %>%
        mutate(Rank = 1:n(), .before="Perturbation") %>%
        mutate(My_summary = "", .after="Perturbation") %>%
        mutate(ProgramID = program.here, .after="Rank") %>%
        ## merge(., ref.table %>% select("Symbol", "TSS.dist.to.SNP", "GWAS.classification"), by.x="Perturbation", by.y="Symbol", all.x=T) %>%
        ## mutate(EC_ctrl_text = ifelse(.$GWAS.classification == "EC_ctrls", "(+)", "")) %>%
        ## mutate(GWAS.class.text = ifelse(grepl("CAD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb"),
        ##                          ifelse(grepl("IBD", GWAS.classification), paste0("_", floor(TSS.dist.to.SNP/1000),"kb_IBD"), ""))) %>%
        ## mutate(Perturb_plus = paste0(Perturbation, GWAS.class.text, EC_ctrl_text)) %>%
    select(Rank, ProgramID, Perturbation, fdr.across.ptb, log2FC, My_summary, FullName, Summary, log2FC.ci.hi, log2FC.ci.lo, fdr, p.value) %>% ## removed Perturb_plus
        arrange(Rank)
}

