## Helen Kang
## Script to gather all program prioritization results
## 230515


## Libraries
suppressPackageStartupMessages(library(conflicted))
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("melt", "reshape2") 
conflict_prefer("slice", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("combine", "dplyr")
conflict_prefer("list", "base")
conflict_prefer("desc", "dplyr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("first", "dplyr")
conflict_prefer("select","dplyr") # multiple packages have select(), prioritize dplyr
conflict_prefer("replace_na", "dplyr")
conflict_prefer("filter", "dplyr")


suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(data.table)
    library(tidyr)
    library(reshape2)
    library(ggplot2)
    library(cowplot) ## plot_grid for arranging many figures into one panel
    library(ggpubr) ## ggarrange
    library(gplots) ## heatmap.2
    ## library(scales) ## geom_tile gradient rescale
    ## library(ggrepel)
    library(stringr)
    library(stringi)
    library(svglite)
    ## library(Seurat)
    ## library(SeuratObject)
    library(xlsx)
    library(readxl)
    library(writexl)
    library(scales) ## geom_tile rescales
})



## optparse
option.list <- list(
    make_option("--cNMF.table", type="character", default="", help="Table with cNMF program result"),
    make_option("--inputdir", type="character", default="./outputs/", help="Input directory"),
    make_option("--outdir", type="character", default="./outputs/aggregated/", help="Output directory"),
    make_option("--figdir", type="character", default="./figures/", help="Output directory"),
    make_option("--sampleName", type="character", default="2kG.library", help="Name of the sample"),
    make_option("--celltype", type="character", default="EC", help="Cell type in GWAS table"),
    make_option("--K.val", type="numeric", default=60, help="The value of K"),
    make_option("--density.thr", type="character", default="0.2", help="concensus cluster threshold, 2 for no filtering"),
    make_option("--outdirsample", type="character", default="/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/", help="path to cNMF analysis results"), ## or for 2n1.99x: "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211116_snakemake_dup4_cells/analysis/all_genes/Perturb_2kG_dup4/K60/threshold_0_2/"
    ## make_option("--num.tests", type="numeric", default=2, help="number of statistical test to do from the top of statistical.test.df.txt"),
    make_option("--trait.names", type="character", default="MAP,PP", help="name of the GWAS traits tested, separated by comma"),
    make_option("--regulator.analysis.type", type="character", default="GWASWide", help="path to statistical test recipe"),
    make_option("--perturbSeq", type="logical", default=FALSE, help="Whether this is a Perturb-seq experiment")
)
opt <- parse_args(OptionParser(option_list=option.list))


## ## Example opt inputs
## ## sdev for K562 gwps 2k overdispersed genes K=90
## opt$inputdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/program_prioritization_GenomeWide/"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/figures/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/program_prioritization_GenomeWide/aggregated/"
## opt$outdir <- paste0(opt$inputdir, "aggregated/")
## opt$sampleName <- "WeissmanK562gwps"
## opt$celltype <- "K562"
## opt$K.val <- 90
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K90/threshold_0_2/"
## opt$perturbSeq <- TRUE
## opt$regulator.analysis.type <- "GenomeWide"
## opt$trait.names <- "HbA1c,Hb,MCH,MCHC,MCV,Plt,RBC,PP,MAP,DBP,SBP"


## ## sdev for EC Perturb-seq all_genes K=60
## opt$sampleName <- "2kG.library"
## opt$inputdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/2kG.library/K60/program_prioritization_GWASWide/"
## opt$figdir <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/figures/2kG.library/all_genes/2kG.library/K60/program_prioritization_GWASWide/aggregated/"
## opt$outdir <- paste0(opt$inputdir, "aggregated/")
## opt$trait.names <- "MAP,PP,DBP,SBP"
## opt$celltype <- "EC"
## opt$K.val <- 60
## opt$outdirsample <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/"
## opt$perturbSeq <- TRUE
## opt$regulator.analysis.type <- "GWASWide"



## plotting constants
mytheme <- theme_classic() + theme(axis.text = element_text(size = 5),
                                   axis.title = element_text(size = 6),
                                   plot.title = element_text(hjust = 0.5, face = "bold", size=6),
                                   axis.line = element_line(color = "black", size = 0.25),
                                   axis.ticks = element_line(color = "black", size = 0.25))


## Directories
FIGDIR <- opt$figdir
OUTDIR <- opt$outdir
INPUTDIR <- opt$inputdir
DATADIR <- "/oak/stanford/groups/engreitz/Users/kangh/ECPerturbSeq2021-Analysis/" ## update this path
check.dir <- c(FIGDIR, OUTDIR)
invisible(lapply(check.dir, function(x) { if(!dir.exists(x)) dir.create(x, recursive=T) }))

## Constants
k <- opt$K.val
SAMPLE <- opt$sampleName
DENSITY.THRESHOLD <- gsub("\\.","_", opt$density.thr)
SUBSCRIPT.SHORT=paste0("k_", k, ".dt_", DENSITY.THRESHOLD)
OUTDIRSAMPLE <- opt$outdirsample
## OUTDIRSAMPLE <- "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210707_snakemake_maxParallel/analysis/2kG.library/all_genes/2kG.library/K60/threshold_0_2/"
celltype <- opt$celltype
trait.ary <- strsplit(opt$trait.names, split=",") %>% unlist

##########################################################################################
## Load Data
## load statistical test results
trait.df.list <- vector("list", length(trait.ary))
for(i in 1:length(trait.ary)) {
    trait.here <- trait.ary[i]
    filename <- paste0(INPUTDIR, "/", trait.here, "/", trait.here, ".program_prioritization.txt")
    trait.df.list[[i]] <- read.delim(filename, stringsAsFactors=F) %>% mutate(GWASTrait = trait.here)
}
program.prioritization.df <- do.call(rbind, trait.df.list)

fileName <- "ProgramPrioritization"
write.table(program.prioritization.df, paste0(OUTDIR, "/", fileName, ".txt"), sep="\t", quote=F, row.names=F)
write.xlsx(program.prioritization.df, file=paste0(OUTDIR, "/", fileName, ".xlsx"), sheetName="Program Prioritization Table", row.names=F, showNA=F)


## set threshold
fdr.thr <- 0.05


numSigPrograms.df <- program.prioritization.df %>% subset(LinkedGenes_fisher.p.adjust < fdr.thr) %>% group_by(GWASTrait) %>% summarize(numSigPrograms = n()) %>%
    arrange(desc(numSigPrograms))
numSigPrograms.df <- rbind(numSigPrograms.df, data.frame(GWASTrait = setdiff(trait.ary, numSigPrograms.df$GWASTrait %>% unique),
                     numSigPrograms=0))
GWASTrait.order <- numSigPrograms.df$GWASTrait %>% rev
numSigPrograms.df <- numSigPrograms.df %>% mutate(GWASTrait = factor(GWASTrait, GWASTrait.order))
    
toplot <- numSigPrograms.df

p <- toplot %>% ggplot(aes(x=GWASTrait, y=numSigPrograms)) + geom_col(fill="gray30") + coord_flip() + mytheme +
    xlab("GWAS Traits") + ylab("# V2G2P Programs") +
    ggtitle(paste0(SAMPLE, "\nK=", k))
numSignificantProgramsPerTrait.barplot <- p

pdf(paste0(FIGDIR, "/NumberOfSignificantProgramsPerTrait_barplot.pdf"), width=1.25, height=1.5)
print(p)
dev.off()


## heatmap of V2G2P significance
significance.df <- program.prioritization.df %>%
    select(ProgramID, GWASTrait, all_of(ifelse(opt$perturbSeq, "LinkedGenes_fisher.p.adjust", "ProgramGene_fisher.p.adjust"))) %>%
    `colnames<-`(c("ProgramID", "GWASTrait", "FDR")) %>%
    mutate(NegLog10FDR = -log10(FDR),
           GWASTrait = factor(GWASTrait, levels=GWASTrait.order),
           ProgramID = factor(ProgramID, levels=paste0("K", k, "_", c(1:k))))


## program ranked by the nubmer of times being a V2G2P program
Program.V2G2P.count.df <- program.prioritization.df %>%
    subset(LinkedGenes_fisher.p.adjust < fdr.thr) %>%
    group_by(ProgramID) %>%
    summarize(numV2G2PTraits = n()) %>%
    arrange(desc(numV2G2PTraits))
Program.V2G2P.count.df <- Program.V2G2P.count.df %>%
    rbind(data.frame(ProgramID = setdiff(paste0("K", k, "_", c(1:k)), Program.V2G2P.count.df$ProgramID),
                     numV2G2PTraits = 0))

toplot <- significance.df %>% mutate(significant = ifelse(NegLog10FDR > -log10(fdr.thr), "*", ""))

palette.NegLog10FDR = colorRampPalette(c("gray", "white", "red"))(n = 100)

## start pdf
pdf(paste0(FIGDIR, "/V2G2P_Significance_heatmap.pdf"), width=12, height=2)
maxNegLog10FDR <- max(toplot$NegLog10FDR)
p <- toplot %>% ggplot(aes(x=ProgramID, y=GWASTrait, fill=NegLog10FDR)) + geom_tile() + coord_equal() + mytheme +
    geom_text(aes(label = significant), nudge_y = -0.25, color = "gray") +
    xlab("Programs") + ylab("GWAS Traits") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    scale_fill_gradientn(name = "FDR (-log10)", colors=palette.NegLog10FDR, limits = c(0, maxNegLog10FDR), values=rescale(c(0, -log10(0.05), maxNegLog10FDR))) +
    ggtitle(paste0(SAMPLE, ", K=", k))
print(p)
toplot <- toplot %>% mutate(ProgramID = factor(ProgramID, levels=Program.V2G2P.count.df$ProgramID))
p <- toplot %>% ggplot(aes(x=ProgramID, y=GWASTrait, fill=NegLog10FDR)) + geom_tile() + coord_equal() + mytheme +
    geom_text(aes(label = significant), nudge_y = -0.25, color = "gray") +
    xlab("Programs") + ylab("GWAS Traits") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    scale_fill_gradientn(name = "FDR (-log10)", colors=palette.NegLog10FDR, limits = c(0, maxNegLog10FDR), values=rescale(c(0, -log10(0.05), maxNegLog10FDR))) +
    ggtitle(paste0(SAMPLE, ", K=", k))
print(p)
significance.heatmap <- p
dev.off()


## figure for program ranked by the nubmer of times being a V2G2P program
pdf(paste0(FIGDIR, "/NumV2G2PLinksToGWASTraitsPerProgram.barplot.pdf"), width=12, height=2)
toplot <- Program.V2G2P.count.df %>% mutate(ProgramID = factor(ProgramID, Program.V2G2P.count.df$ProgramID))
p <- toplot %>% ggplot(aes(x=ProgramID, y=numV2G2PTraits)) + geom_col(fill="gray30") + mytheme +
    xlab("Programs") + ylab("# V2G2P Links\nto GWAS Traits") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0(SAMPLE, ", K=", k))
print(p)
numV2G2PLinksToGWASTraitsPerProgram.barplot <- p
toplot <- Program.V2G2P.count.df %>% mutate(ProgramID = factor(ProgramID, paste0("K", k, "_", c(1:k))))
p <- toplot %>% ggplot(aes(x=ProgramID, y=numV2G2PTraits)) + geom_col(fill="gray30") + mytheme +
    xlab("Programs") + ylab("# V2G2P Links\nto GWAS Traits") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste0(SAMPLE, ", K=", k))
print(p)
dev.off()
    

## ## Aggregated plot for significance heatmap
## p.leftPlot <- plot_grid(numV2G2PLinksToGWASTraitsPerProgram.barplot, significance.heatmap,
##                         nrow=2, rel_heights=c(1,2), align="v", axis="lr")
## p.rightPlot <- plot_grid(ggplot(data = data.frame(x=c(1,2),y=c(1,2)), aes(x=x, y=y)) + geom_point() + mytheme + ggtitle("placeholder"), numSignificantProgramsPerTrait.barplot,
##                          nrow=2, rel_heights=c(1,2), align="v", axis="lr")
## p.mainPlot <- plot_grid(p.leftPlot, p.rightPlot, ncol=2, rel_widths=c(12,1.5), align="h", axis="tb")
p.mainPlot <- plot_grid(numV2G2PLinksToGWASTraitsPerProgram.barplot, ggplot(data = data.frame(x=c(1,2),y=c(1,2)), aes(x=x, y=y)) + geom_point() + mytheme + ggtitle("placeholder"), significance.heatmap, numSignificantProgramsPerTrait.barplot,
                        nrow=2, rel_heights=c(1,2), rel_widths=c(12,1.5), align = "hv", axis="tbl")
pdf(paste0(FIGDIR, "/aggregated_significance.heatmap.pdf"), width=14, height=3)
print(p.mainPlot)
dev.off()

