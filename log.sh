#!/bin/bash

#################################################
## Helen Kang
## 211011 210625 (adapted from Perturb-seq_CAD snakemake pipeline log.sh)
## Script to test and run cNMF snakemake pipeline

# conda activate EngreitzLab

PROJECT=/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/
SNAKEMAKEDIR=$OAK/Users/kangh/cNMF_pipeline
QSUB=$GROUP_HOME/bin/quick-sub

cd ${SNAKEMAKEDIR}

# snakemake --use-conda
 
snakemake -n --rerun-incomplete --configfile config/config.json --quiet

## pipeline visualization 
snakemake --rerun-incomplete --configfile config/config.json --forceall --dag | dot -Tpdf > ${PROJECT}/dag.pdf


## running the pipeline via sbatch
LOGS=${PROJECT}/logs
mkdir -p ${LOGS}
snakemake -k --configfile config/config.json --restart-times 1 --rerun-incomplete --nolock --jobs 1000 --cluster "sbatch -n 1 -c 1 --mem {params.mem_gb}G -t {params.time} -p owners,normal -J Perturb-seq_CAD_{rule} -o $LOGS/{rule}_{wildcards}.qout -e $LOGS/{rule}_{wildcards}.qe"


## submit to computing node
LOGS=${PROJECT}/logs
mkdir -p ${LOGS}
LOGS_HERE=${LOGS}/2kG.library
mkdir -p ${LOGS_HERE}
jobname=2kG.library_qsub
tosubmit="source ~/.bashrc; conda activate cnmf_env; cd ${SNAKEMAKEDIR}; snakemake -k --configfile config/config.json --restart-times 1 --rerun-incomplete --nolock --jobs 1000 --cluster \"sbatch -n 1 -c 1 --mem {params.mem_gb}G -t {params.time} -p owners,normal -J Perturb-seq_CAD_{rule} -o $LOGS_HERE/{rule}_{wildcards}.qout -e $LOGS_HERE/{rule}_{wildcards}.qe\" "
$QSUB -p engreitz -j ${jobname} -s ${LOGS}/${jobname}.qsh -o ${LOGS}/${jobname}.qout -m 6G -t 168:00:00 ${tosubmit}
