#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 15:58:36 2025

Description: Testing the hvg filtering in SCLC 23169 sample.
Idea here is have a one-sample script ready to be converted to a bulk script.

@author: bellwu
"""

# %% ---- 1.0 set up environment ----
import anndata as ad
from pyxenium import xen_config as xc
import scanpy as sc
import numpy as np
# %%% ---- 1.1 Load anndata ----
xen_dir = xc.xen_bwu
adata = ad.read_h5ad(xen_dir / '64312_filtered.h5ad')
adata_all = ad.read_h5ad
# %%% ---- 1.2 Test out scanpy hvg for 3000 ----
# look at data variables
adata
adata.var.columns
adata.obs.columns
# remove hvgs
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
# take highly variable genes for 3000
hvg_3000 = adata.var['highly_variable_rank']
hvg_3000.dropna() # remove NaN from array
a_hvg = adata[:, hvg_3000.dropna()]
a_hvg.shape
# calculate how many counts removed
hvg_dif = adata.X.sum() - a_hvg.X.sum()
# percentage of counts lost
print(f"Total counts removed after hvg filter is {hvg_dif}")
print(f"Percentage of counts lost {hvg_dif / adata.X.sum() * 100}")
# %% ---- 1.2 Do the genes that are lost make sense? ----
not_hvg = hvg_3000.isna()
not_hvg = adata[:, not_hvg]
names = not_hvg.var_names
# %%%% ---- 1.2.1 testing if total_counts is good metric
adata.obs['total_counts'].sum()
a_hvg.obs['total_counts'].sum()

# test if total_counts good metric
test = adata[:, 1]
test.shape
test.obs['total_counts'].sum()
'''
See no differences in taking top 3000 hvg and total hvg with total_counts 
regardless of shape, must mean that total_counts is always going to be certain
number unless I filter out somehow.

Makes sense as total_counts calculate the counts within each cell
'''
test.X
# this will give sum of all counts of all cells for specific gene 
test.X.sum() 
# %%% ---- 1.2 Test out scanpy hvg for 2000 ----
# remove hvgs
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
# take highly variable genes for 2000
hvg_2000 = adata.var['highly_variable_rank']
hvg_2000.dropna() # remove NaN from array
a_hvg = adata[:, hvg_2000.dropna()]
a_hvg.shape
# calculate how many counts removed
hvg_dif = adata.X.sum() - a_hvg.X.sum()
# percentage of counts lost
print(f"Total counts removed after hvg filter is {hvg_dif}")
print(f"Percentage of counts lost {hvg_dif / adata.X.sum() * 100}")
# %% ---- 1.3 PC on anndata ----
'''
Below is analysis testing if different levels of hvg has an effect on the PC
outputs. What I am looking for if there are any changes to the eigenvalues
of the PCs w/o and w/ hvg filtering. More specifically are there the same number
of nontrivial PCs? 
- Obvious ratios will change since there is less total dimensions to be explained.
'''
# %%% ---- 1.3.1 normalize w/o hvg and compute PCs ----
adata.layers['counts'] = adata.X.copy() # save counts 
sc.pp.normalize_total(adata, target_sum=1e6, inplace=True)
sc.pp.log1p(adata)
sc.pp.pca(adata, use_highly_variable=False)
# %%%% ---- 1.3.1.1 investigating the PCs ----
pc_raw_ratio = adata.uns['pca']['variance_ratio']
pc_raw = adata.uns['pca']['variance']
adata.obsm['X_pca'].shape 
adata.varm['PCs'].shape
adata.varm['PCs']
adata.shape
# %%% ---- 1.3.2 normalize w/ hvg ----
sc.pp.normalize_total(adata, target_sum=1e6, inplace=True)
sc.pp.log1p(adata)
adata_hvg = adata.copy() # create copy
sc.pp.pca(adata_hvg, use_highly_variable=True)
# %%%% ---- 1.3.2.1 investigating PCs ----
adata_hvg.uns['pca']['variance_ratio'][0]
adata.uns['pca']['variance_ratio'][0]
# compute MP distribution
# %%%%% ---- 1.3.2.2 MP calculation ----
def mp_edge(p, n, sigma2=1.0):
    '''
    Compute Marchenko-Pastur lower and upper edges (λ-, λ+) for a (p x n) data
    matrix with noise variance sigma2 where n is number of features and p is
    number of samples
    '''
    q = p / n
    l_min = sigma2 * (1 - np.sqrt(q))**2
    l_max = sigma2 * (1 + np.sqrt(q))**2
    return(l_min, l_max)
# %%%% ----- 1.3.2.3 MP for hvg and raw data -----
adata_hvg.var['highly_variable'].sum()
adata.shape
hvg_min, hvg_max = mp_edge(254982, 2000)
n_min, n_max = mp_edge(254982, 5000)
## identify no. non-trivial PCs:
(adata_hvg.uns['pca']['variance'] > hvg_max).sum()
(adata.uns['pca']['variance'] > n_max).sum()




