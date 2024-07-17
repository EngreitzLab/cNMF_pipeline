#!/bin/bash

#################################################
## Helen Kang
## 220710
## Script to test and run cNMF snakemake pipeline

conda activate cnmf_env

PROJECT=/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220505_snakemake_moreK_findK/
SNAKEMAKEDIR=$OAK/Users/kangh/cNMF_pipeline
QSUB=$GROUP_HOME/bin/quick-sub

cd ${SNAKEMAKEDIR}

## dry run
snakemake -n --rerun-incomplete --configfile ${PROJECT}/config.json --quiet

## pipeline visualization 
snakemake --rerun-incomplete --configfile ${PROJECT}/config.json --forceall --dag | dot -Tsvg > ${PROJECT}/dag.pipeline.svg


## running the pipeline interactively
LOGS=${PROJECT}/logs
mkdir -p ${LOGS}
snakemake -k --configfile ${PROJECT}/config.json --restart-times 0 --rerun-incomplete --nolock --jobs 1000 --cluster "sbatch -n 1 -c 1 --mem {params.mem_gb}G -t {params.time} -p owners,normal -J Perturb-seq_CAD_{rule} -o $LOGS/{rule}_{wildcards}.qout -e $LOGS/{rule}_{wildcards}.qe"


## submit the nakemake pipeline as a job
jobname=2kG.library_qsub
echo ${jobname}
LOGS=${PROJECT}/logs
LOGS_HERE=${PROJECT}/logs/${jobname} 
mkdir -p ${LOGS} ${LOGS_HERE}
tosubmit="source ~/.bashrc; conda activate cnmf_env; cd ${SNAKEMAKEDIR}; snakemake -k --configfile ${PROJECT}/config.220505.json --restart-times 1 --rerun-incomplete --nolock --jobs 1000 --cluster \"sbatch -n 1 -c 1 --mem {params.mem_gb}G -t {params.time} -p owners,normal -J cNMF_{rule} -o $LOGS_HERE/{rule}_{wildcards}.qout -e $LOGS_HERE/{rule}_{wildcards}.qe\" "
echo ${tosubmit}
[ -f "${LOGS}/${jobname}.qout" ] && rm ${LOGS}/${jobname}.qout
$QSUB -p engreitz,normal -j ${jobname} -s ${LOGS}/${jobname}.qsh -o ${LOGS}/${jobname}.qout -m 6G -t 48:00:00 ${tosubmit}


## Broad UGER:
tosubmit="source /broad/software/scripts/useuse; source /broad/software/free/Linux/redhat_7_x86_64/pkgs/anaconda2_5.0.1/etc/profile.d/conda.sh;  use Anaconda; reuse UGER; conda activate cnmf_env; cd ${SNAKEMAKEDIR}; snakemake -k --configfile ${PROJECT}/config.json --restart-times 0 --rerun-incomplete --nolock --jobs 300 --cluster \"qsub -l h_vmem={params.mem_gb}G -l h_rt={params.time} -N V2G2P_{rule} -o $LOGS_HERE -e $LOGS_HERE\" "
