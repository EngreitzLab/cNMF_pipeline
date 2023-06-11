## Convert data txt file to h5 file
## Helen Kang
## 211014

import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
# import matplotlib.pyplot as plt
#from IPython.display import Image
import scanpy as sc
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--inputPath', dest = 'inputPath', type=str, default='/oak/stanford/groups/engreitz/Users/kangh/2009_endothelial_perturbseq_analysis/cNMF/data/no_IL1B_filtered.normalized.ptb.by.gene.mtx.txt', help = 'path to count matrix file')
# parser.add_argument('--outPath', dest = 'outPath', type=str, help = 'path to output folder')
# parser.add_argument('--run_name',dest ='run_name', type=str, help = 'sample name')
parser.add_argument('--output_h5ad_mtx', dest = 'output_h5ad_mtx', type=str, help = 'path to output folder', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/211014_test_txt_to_h5ad/outputs/test.h5ad')
parser.add_argument('--output_gene_name_txt', dest = 'output_gene_name_txt', type=str, help = 'path to gene name output txt file', default='/oak/stanford/groups/engreitz/Users/kangh/TeloHAEC_Perturb-seq_2kG/211011_Perturb-seq_Analysis_Pipeline_scratch/211014_test_txt_to_h5ad/outputs/test.h5ad.all.genes.txt')
args = parser.parse_args()


# ## 200 gene library 230610, no IL1B
# args.inputPath = '/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/data/no_IL1B.raw.h5ad'
# args.output_h5ad_mtx = '/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/data/no_IL1B.h5ad'
# args.output_gene_name_txt = '/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/data/no_IL1B.h5ad.all.genes.txt'


# ## 200 gene library 230610, plus IL1B
# args.inputPath = '/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/data/plus_IL1B.raw.h5ad'
# args.output_h5ad_mtx = '/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/data/plus_IL1B.h5ad'
# args.output_gene_name_txt = '/oak/stanford/groups/engreitz/Users/kangh/tutorials/2306_V2G2P_prep/data/plus_IL1B.h5ad.all.genes.txt'


## load scanpy from txt
adata = sc.read(args.inputPath)

# remove non protein coding genes
df_columns = pd.Series(adata.var_names.values.flatten()).astype(str)
AL = df_columns[df_columns.str.contains('^AL[0-9][0-9][0-9][0-9][0-9][0-9]\\.').tolist()]
AC = df_columns[df_columns.str.contains('^AC[0-9][0-9][0-9][0-9][0-9][0-9]\\.').tolist()]
AP = df_columns[df_columns.str.contains('^AP[0-9][0-9][0-9][0-9][0-9][0-9]\\.').tolist()]
LINC = df_columns[df_columns.str.contains('LINC').tolist()]
allPattern = df_columns[df_columns.str.contains('^[A-Za-z][A-Za-z][0-9][0-9][0-9][0-9][0-9][0-9]\\.').tolist()]
toremove = allPattern.append(LINC)
tokeep = ~df_columns.isin(toremove)
adata = adata[:,tokeep]
# filter cells
sc.pp.filter_cells(adata, min_genes=200) # filter cells with fewer than 200 genes
sc.pp.filter_cells(adata, min_counts=200)  # This is a weaker threshold than above. It is just to population the n_counts column in adata
sc.pp.filter_genes(adata, min_cells=10) # filter genes detected in fewer than 3 cells

# save to h5ad file
sc.write(args.output_h5ad_mtx, adata)

# get gene names after filtering
filtered_genes = pd.DataFrame(adata.var_names.values.flatten()).astype(str)

# save gene names
filtered_genes.to_csv(args.output_gene_name_txt, sep="\t", header=False, index=False)
