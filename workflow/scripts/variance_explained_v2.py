## Helen Kang
## compute variance explained for cNMF components

import numpy as np
import scipy
import pandas as pd
# import scipy.sparse as sp
import scanpy as sc
import argparse
import os
import re
from cnmf import cNMF
from sklearn.decomposition import PCA


## argparse
parser = argparse.ArgumentParser()

## add arguments
parser.add_argument('--path_to_topics', type=str, help='path to the topic (cNMF directory) to project data on', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220716_snakemake_overdispersedGenes/analysis/top2000VariableGenes_acrossK/')
parser.add_argument('--topic_sampleName', type=str, help='sample name for topics to project on, use the same sample name as used for the cNMF directory', default='2kG.library_overdispersedGenes')
parser.add_argument('--X_normalized', type=str,  help='path to normalized input cell x gene matrix from cNMF pipeline', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220716_snakemake_overdispersedGenes/analysis/top2000VariableGenes_acrossK/2kG.library_overdispersedGenes/cnmf_tmp/2kG.library_overdispersedGenes.norm_counts.h5ad')
parser.add_argument('--outdir', dest = 'outdir', type=str, help = 'path to output directory', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/220716_snakemake_overdispersedGenes/analysis/top2000VariableGenes/2kG.library_overdispersedGenes/K60/threshold_0_2/')
parser.add_argument('--k', dest = 'k', type=int, help = 'number of components', default='60')
parser.add_argument('--density_threshold', dest = 'density_threshold', type=float, help = 'component spectra clustering threshold, 2 for no filtering, recommend 0_2 (means 0.2)', default="0.2")

# ## sdev for K562 gwps
# parser.add_argument('--path_to_topics', type=str, help='path to the topic (cNMF directory) to project data on', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes_acrossK/')
# parser.add_argument('--topic_sampleName', type=str, help='sample name for topics to project on, use the same sample name as used for the cNMF directory', default='WeissmanK562gwps')
# parser.add_argument('--X_normalized', type=str,  help='path to normalized input cell x gene matrix from cNMF pipeline', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes_acrossK/WeissmanK562gwps/cnmf_tmp/WeissmanK562gwps.norm_counts.h5ad')
# parser.add_argument('--outdir', dest = 'outdir', type=str, help = 'path to output directory', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/230104_snakemake_WeissmanLabData/analysis/top2000VariableGenes/WeissmanK562gwps/K35/threshold_0_2/')
# parser.add_argument('--k', dest = 'k', type=int, help = 'number of components', default='35')
# parser.add_argument('--density_threshold', dest = 'density_threshold', type=float, help = 'component spectra clustering threshold, 2 for no filtering, recommend 0_2 (means 0.2)', default="0.2")


args = parser.parse_args()

# ## sdev debug for Disha's error
# args.X_normalized = "/oak/stanford/groups/engreitz/Users/kangh/scratch_space/230612_debug_cNMF_pipeline_variance_explained/Merge_Pauletal_subset_SMC.norm_counts.h5ad"
# args.outdir = "/oak/stanford/groups/engreitz/Users/kangh/scratch_space/230612_debug_cNMF_pipeline_variance_explained/"



sample = args.topic_sampleName
# output_sample = args.output_sampleName
# tpm_counts_path = args.tpm_counts_path
OUTDIR = args.outdir
selected_K = args.k
density_threshold = args.density_threshold
output_directory = args.path_to_topics
run_name = args.topic_sampleName

if not os.path.exists(OUTDIR):
    raise Exception("Output directory does not exist")

cnmf_obj = cNMF(output_dir=output_directory, name=run_name)
usage_norm, gep_scores, gep_tpm, topgenes = cnmf_obj.load_results(K=selected_K, density_threshold=density_threshold)


# X_original = sc.read_h5ad(args.tpm_counts_path)
X_norm = sc.read_h5ad(args.X_normalized)
# sc.pp.normalize_per_cell(X_original, counts_per_cell_after=1e6) ## normalize X to TPM
# X_original.X[0:10,:].todense().sum(axis=1) ## check normalization results


## functions
def compute_Var(X):
    if scipy.sparse.issparse(X):
        return np.sum(np.var(X.todense(), axis=0, ddof=1))
    else:
        return np.sum(np.var(X, axis=0, ddof=1))

# ## first turn X into TPM
# X_tpm_dense = X_original.X.todense()
# X_tpm_dense[0:10,].sum(axis=1) ## check TPM normalization

# X = X_norm.X.todense() ## 221203
X = X_norm.X
H_path = cnmf_obj.paths['consensus_spectra__txt'] % (selected_K, '0_2') ## median_spectra_file
H_df = pd.read_csv(H_path, sep='\t', index_col=0).T
H = H_df.to_numpy()
H = (H/H.sum(0))
W_path = cnmf_obj.paths['consensus_usages__txt'] % (selected_K, '0_2') ## median_spectra_file
W_df = pd.read_csv(W_path, sep='\t', index_col=0)
W = W_df.to_numpy()
WH = W @ H.T
# diff = X - WH

# X_col_sd = np.std(X_tpm_dense, axis=0, ddof=1) ## column normalization, shape: (1,17472)
# type(gep_tpm) ## data frame
# gep_tpm.shape ## (17472, 60)
# H_tmp = gep_tpm.to_numpy().T ## shape: (60, 17472)
# H = gep_tpm.to_numpy().T / X_col_sd ## normalize TPM spectra matrix (gene x component)
# X = X_tpm_dense / X_col_sd
# X[:,0:10].std(axis=0) ## check if standard deviation of gene expression is 1
# W = usage_norm.to_numpy()
# W[0:10,].sum(axis=1) ## to check if sum of each cell's usage is 1
# WH = W @ H
diff = X - WH
diff_sumOfSquaresError = (np.asarray(diff)**2).sum()
# X_sumOfSquares = (np.asarray(X)**2).sum()
# WH_sumOfSquares = (np.asarray(WH)**2).sum()
Var_diff = compute_Var(diff)
Var_X = compute_Var(X)
TotalVarianceExplained = 1 - Var_diff / Var_X

def computeVarianceExplained(X, H, Var_X, i):
    if not isinstance(H, (pd.DataFrame)):
        B_k = X @ H[i,:].T / np.sqrt((np.asarray(H[i,:])**2).sum())
        numerator = compute_Var(X - np.outer(B_k, H[i,:]))
    else:
        B_k = X @ H.iloc[i,:] / np.sqrt((H.iloc[i,:]**2).sum())
        numerator = compute_Var(X - np.outer(B_k, H.iloc[i,:]))
    return (1 - numerator / Var_X)

## initialize storage variable
V_k = np.empty([selected_K])

for i in range(selected_K):
    print(i)
    V_k[i] = computeVarianceExplained(X, H.T, Var_X, i)

ProgramID = ['K' + str(selected_K) + '_' + str(i+1) for i in range(selected_K)]

metrics_df = pd.DataFrame({'VarianceExplained': V_k,
                           'ProgramID': ProgramID })
metrics_summary = pd.DataFrame({'Sum' : metrics_df['VarianceExplained'].sum(),
                                'Median' : metrics_df['VarianceExplained'].median(),
                                'Max' : metrics_df['VarianceExplained'].max(),
                                'Total' : TotalVarianceExplained},
                               index = [0])
metrics_df.to_csv(os.path.join(OUTDIR, "metrics.varianceExplained.df.txt"), index = None, sep="\t")
metrics_summary.to_csv(os.path.join(OUTDIR, "summary.varianceExplained.df.txt"), index = None, sep="\t")
