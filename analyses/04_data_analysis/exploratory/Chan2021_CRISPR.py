"""
Created on 2026.01.16

Description: Exploratory analysis of Chan et al. 2021 data for 
CRISPR hits of own in vivo screen.

@author: bellwu
"""
# %% ---- 1.0 set up local environment ----
import scanpy as sc
from pyxenium import xen_config as xc
import pandas as pd
import numpy as np
# %% ---- 2.0 Load in preprocessed Chan et al. data ----
## directory variables
Chan_dir = xc.Chan_bwu

## reading h5ad file
adata_Chan = sc.read_h5ad(Chan_dir / "scChan_2021_CellxGene_combined.h5ad")
# %% ---- 3.0 Explore CRISPR hit genes ----
# list of hits
hits = ['NRK', 'TWNK', 'PEO1', 'SETD1B', 'PECAM1', 'SLC6A11']
# check how many hits in data
adata_Chan.var['feature_name'].isin(hits).sum()
# subset out hits
mask = adata_Chan.var['feature_name'].isin(hits)
adata_hits = adata_Chan[:, mask]
adata_hits.var['feature_name'] # check which hits present in data
# %% # ---- 3.1 Check expression levels ----
# %% # check normalization of Chan data 
X = adata_Chan.X
totals = np.array(X.sum(axis=1)).ravel() # convert to np.array and then flatten with to 1D
print(np.median(totals)) # check median total counts per cell
'''
Chan data is not normalized to any particular value...
'''
# %% # normalize Chan data to 1e4
## reset X to raw counts
adata_Chan.X = adata_Chan.raw.X 
# adata_Chan.X[:20,:20].toarray() # check first 20 cells and genes
## normalize
sc.pp.normalize_total(adata_Chan, target_sum=1e4, inplace=True)
sc.pp.log1p(adata_Chan)
# %% # check if normalized properly
X = adata_Chan.X
X_exp = X.expm1()
totals = np.array(X_exp.sum(axis=1)).ravel() # convert to np.array and then flatten with to 1D
print(np.median(totals))
# %%
