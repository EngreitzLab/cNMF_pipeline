#!/bin/bash

#################################################
## Helen Kang
## 211101
## Run cNMF on the 20 TeloHAEC full scale library samples separately (on all cells, not accounting for singlets)

conda activate cnmf_env

PROJECT=/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211101_20sample_snakemake/
SNAKEMAKEDIR=$OAK/Users/kangh/github_static/cNMF_pipeline
QSUB=$GROUP_HOME/bin/quick-sub

cd ${SNAKEMAKEDIR}

# snakemake --use-conda
 
snakemake -n --rerun-incomplete --configfile ${PROJECT}/config.json --quiet

sample=scRNAseq_2kG_11AMDox_1
datadir=/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210621_10XCloud/${sample}/outs/filtered_feature_bc_matrix
snakemake -n --configfile ${PROJECT}/config.json --config dataDir=${datadir} sampleName=${sample} --quiet


## pipeline visualization 
snakemake --rerun-incomplete --configfile ${PROJECT}/config.json --config dataDir=${datadir} sampleName=${sample} --forceall --dag | dot -Tpdf > ${PROJECT}/dag.${sample}.pdf


## run the pipeline interactively
LOGS=${PROJECT}/logs
mkdir -p ${LOGS}
LOGS_HERE=${LOGS}/${sample}
mkdir -p ${LOGS_HERE}
snakemake -k --configfile ${PROJECT}/config.json --config dataDir=${datadir} sampleName=${sample} --restart-times 1 --rerun-incomplete --nolock --jobs 35 --cluster "sbatch -n 1 -c 1 --mem {params.mem_gb}G -t {params.time} -p owners,normal -J cNMF_{rule} -o $LOGS_HERE/{rule}_{wildcards}.qout -e $LOGS_HERE/{rule}_{wildcards}.qe"

# # snakemake -k --configfile ${PROJECT}/config.json --restart-times 0 --rerun-incomplete --nolock --jobs 1000 --cluster "sbatch -n 1 -c 1 --mem {params.mem_gb}G -t {params.time} -p owners,normal,engreitz -J cNMF_{rule} -o $LOGS_HERE/{rule}_{wildcards}.qout -e $LOGS_HERE/{rule}_{wildcards}.qe"


LOGS=${PROJECT}/logs
mkdir -p ${LOGS}
while IFS= read -r sample
do
    echo ${sample}
    LOGS_HERE=${LOGS}/${sample}
    mkdir -p ${LOGS_HERE}
    shortSampleName=$( echo ${sample} | sed "s/scRNAseq_2kG_//g" )
    jobname=${shortSampleName}_qsub
    echo ${jobname}
    datadir=/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210621_10XCloud/${sample}/outs/filtered_feature_bc_matrix
    tosubmit="source ~/.bashrc; conda activate cnmf_env; cd ${SNAKEMAKEDIR}; snakemake -k --configfile ${PROJECT}/config.json --config dataDir=${datadir} sampleName=${sample} --restart-times 1 --rerun-incomplete --nolock --jobs 35 --cluster \"sbatch -n 1 -c 1 --mem {params.mem_gb}G -t {params.time} -p owners,normal -J cNMF_{rule} -o $LOGS_HERE/{rule}_{wildcards}.qout -e $LOGS_HERE/{rule}_{wildcards}.qe\" "
    echo ${tosubmit}
    $QSUB -p engreitz -j ${jobname} -s ${LOGS}/${jobname}.qsh -o ${LOGS}/${jobname}.qout -m 6G -t 168:00:00 ${tosubmit}
done < "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/210528_Broad_fastq/sampleName.txt"


## check pipeline status
tail ${LOGS}/${jobname}.qout
