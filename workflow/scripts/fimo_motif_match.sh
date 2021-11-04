## extract fasta sequences from genome coordinates
## Helen Kang
## 211028


############################################################
## Directories
COORD=$1
FASTA=$2
OUTFASTA=$3


# PROJECT=$OAK/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/2104_all_genes/
# DATADIR=$PROJECT/data/
# FILEDIR=$OAK/Data/hg38
# OUTDIR=${PROJECT}/outputs/
# ABCDIR=$OAK/Projects/ABC/200220_CAD/ABC_out/TeloHAEC_Ctrl/
# TOPDATADIR=$OAK/Users/kangh/2009_endothelial_perturbseq_analysis/data/
# TOPDATADIRABC=$OAK/Users/kangh/2009_endothelial_perturbseq_analysis/data/ABC/
# SCRATCHDIR=${SCRATCH}/210509_topic_motif_enrichment/


# mkdir -p $OUTDIR
# mkdir -p $DATADIR
# mkdir -p $TOPDATADIRABC
# mkdir -p $SCRATCHDIR
# LOG=${PROJECT}/logs/
# mkdir -p $LOG
# QSUB=/home/groups/engreitz/bin/quick-sub
# hg38FASTA=/oak/stanford/groups/engreitz/Data/hg38/Sequence/hg38.fa
# hg19FASTA=/oak/stanford/groups/engreitz/Data/hg19/Sequence/hg19.fa



############################################################
## get fasta for the enhancer regions
# chr, start, end, name, class, activity_base, TargetGene, TargetGeneTSS, TargetGeneExpression, TargetGenePromoterActivityQuantile, TargetGeneIsExpressed, distance, isSelfPromoter, powerlaw_contact, powerlaw_contact_reference, hic_contact, hic_contact_pl_scaled, hic_pseudocount, hic_contact_pl_scaled_adj, ABC.Score.Numerator, ABC.Score, powerlaw.Score.Numerator, powerlaw.Score, CellType
bedtools getfasta -name -fi ${FASTA} -bed <(awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3"|"$4"|"$7}' ${COORD}) -fo ${OUTFASTA}

    # bedtools getfasta -name -fi ${hg19FASTA} -bed <(awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3"|"$4"|"$7}' ${COORD}) -fo ${TOPDATADIRABC}/${sample}_Predictions.AvgHiC.ABC0.015.minus150.fa


# ############################################################
# ## FIMO
# THRESHOLD=(1.0E-4 1.0E-3)
# THRESHOLD=(1.0E-6)
# for threshold in "${THRESHOLD[@]}"
# do
#     for sample in "${SAMPLE[@]}"
#     do
# 	JOBNAME=${sample}.ABC.FIMO.${threshold}
# 	echo "$JOBNAME"
# 	FIMO_OUTDIR=${DATADIR}/fimo_out_ABC_${sample}_thresh${threshold}/
# 	mkdir -p ${FIMO_OUTDIR}
# 	$QSUB -j $JOBNAME -s $LOG/$JOBNAME.qsh -o $LOG/$JOBNAME.qout -m 64G -t 8:00:00 "source ~/.bashrc; fimo -oc ${FIMO_OUTDIR}/ --verbosity 1 --thresh ${threshold} /oak/stanford/groups/engreitz/Data/motif/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme ${TOPDATADIRABC}/${sample}_Predictions.AvgHiC.ABC0.015.minus150.fa \
# 	echo ${FIMO_OUTDIR}/fimo.tsv | tr "|" "\t" > ${FIMO_OUTDIR}/fimo.formatted.tsv"
#     done
# done




# ## format FIMO results
# input_tsv=${FIMO_OUTDIR}/fimo.tsv 
# output_tsv=${FIMO_OUTDIR}/fimo.formatted.tsv
# Rscript format.FIMO.ABC.results.R \
# --input.tsv ${input_tsv} \ 
# --output.tsv ${output_tsv} 



