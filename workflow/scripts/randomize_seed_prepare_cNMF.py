## Helen Kang
## Script to randomize starting seed for cNMF
import os
import pandas as pd
# import anndata as ad
# import scanpy as sc
import numpy as np
# import yaml
import argparse
# from sklearn.cluster import KMeans
# from sklearn.metrics import silhouette_score

# from cnmf import cNMF

## Parse in data
parser = argparse.ArgumentParser()

parser.add_argument('--nmf_params', type=str, help='path to NMF parameters to change seed', default='/scratch/groups/engreitz/Users/kangh/cNMF_pipeline/241015_V2G2P_HCASM/top2000VariableGenes/K80/worker0/HCASM.library/cnmf_tmp/HCASM.library.nmf_params.df.npz')

args = parser.parse_args()




## referenced Dylan Kotliar's cNMF code:
def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj

def save_df_to_npz(obj, filename):
    np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)

## Load parameters data frame
replicate_params = load_df_from_npz(args.nmf_params)

## Randomize nmf_seed
replicate_params['nmf_seed'] = np.random.randint(1, max(replicate_params['nmf_seed']), replicate_params.shape[0])

## save updated data frame
save_df_to_npz(replicate_params, args.nmf_params)


