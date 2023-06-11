## Helen Kang
## Script to convert cell x gene matrix into h5ad file for IGVF format
## 230523

import numpy as np
import pandas as pd
# import scipy.sparse as sp
# import scanpy as sc
import anndata as ad
from cnmf import cNMF
import argparse
import re

## argparse
parser = argparse.ArgumentParser()
parser.add_argument('--path_to_topics', type=str, help='path to the topic (cNMF directory) to project data on', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes_acrossK')
parser.add_argument('--topic_sampleName', type=str, help='sample name for topics to project on, use the same sample name as used for the cNMF directory', default='WeissmanK562gwps')
# parser.add_argument('--tpm_counts_path', type=str, help='path to tpm input cell x gene matrix', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220716_snakemake_overdispersedGenes/analysis/top2000VariableGenes_acrossK/2kG.library_overdispersedGenes/cnmf_tmp/2kG.library_overdispersedGenes.tpm.h5ad') #/scratch/groups/engreitz/Users/kangh/cNMF_pipeline/220505_snakemake_moreK_findK/all_genes/K60/worker0/2kG.library/cnmf_tmp/2kG.library.norm_counts.h5ad')
parser.add_argument('--outdir', dest = 'outdir', type=str, help = 'path to output directory', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K10/threshold_0_2/IGVF_format/')
parser.add_argument('--k', dest = 'k', type=int, help = 'number of components', default='10')
parser.add_argument('--density_threshold', dest = 'density_threshold', type=float, help = 'component spectra clustering threshold, 2 for no filtering, recommend 0_2 (means 0.2)', default="0.2")
parser.add_argument('--barcode_dir', dest = 'barcode_dir', type=str, default='/oak/stanford/groups/engreitz/Users/kangh/collab_data/IGVF/mouse_ENCODE_heart/auxiliary_data/snrna/heart_Parse_10x_integrated_metadata.csv', help='Directory to barcodes, require columns CBC and Gene')

args = parser.parse_args()

sample = args.topic_sampleName
# output_sample = args.output_sampleName
# tpm_counts_path = args.tpm_counts_path
OUTDIR = args.outdir
selected_K = args.k
density_threshold = args.density_threshold
output_directory = args.path_to_topics
run_name = args.topic_sampleName
barcode_dir = args.barcode_dir

cnmf_obj = cNMF(output_dir=output_directory, name=run_name)
usage_norm, gep_scores, gep_tpm, topgenes = cnmf_obj.load_results(K=selected_K, density_threshold=density_threshold)


## load cell barcode inforamtion
# barcode_dir = "/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/data/K562_gwps_raw_singlecell_01_metadata.txt"
if barcode_dir.endswith('csv'):
    barcodes_df = pd.read_csv(barcode_dir, index_col="CBC")
else:
    barcodes_df = pd.read_csv(barcode_dir, sep="\t", index_col="CBC")

## organize program name
programNames = [run_name + "_K" + str(selected_K) + "_" + str(i) for i in usage_norm.columns]
programNames_df = pd.DataFrame({"ProgramNames": programNames}, index=programNames)

usage_norm.columns = programNames

## create AnnData
adata = ad.AnnData(
    X = usage_norm,
    obs = barcodes_df,
    var = programNames_df
)

## save results
fileName = OUTDIR + sample + ".k_" + str(selected_K) + ".dt_" + re.sub("[.]", "_", str(density_threshold)) + ".cellxgene.h5ad"
adata.write_h5ad(fileName)
