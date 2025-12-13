#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 16:06:16 2025

Description: Continuation of clustering

@author: bellwu
"""

# %% ---- 1.0 Set up environment ----
from pyxenium import xen_config as xc
import anndata as ad
import scanpy as sc
import time
# %% ---- 1.1 Load files ----
## continue after neighbour analysis
xen_dir = xc.xen_bwu
h5ad_path = xen_dir / "SCLC_PCA.h5ad"
adata = ad.read_h5ad(h5ad_path)
# %% ---- 2.0 Preprocessing ----
'''
This section includes all the preprocessing steps for the anndata:
    1) Normalization to CPM
    2) Construct PCA
    3) Build kNN from PCs (using graph construction method that is same as UMAP)
    3) Visualize 2D via UMAP projection
    4) Add leiden communities
'''
## Setting up counter for computation time
start = time.perf_counter()
end = time.perf_counter()
# %%% ---- 2.1 Normalization and transformation ----
# ## create copy of counts to new layer
# adata.layers['count'] = adata.X.copy()
# ## normalize to CPM
# sc.pp.normalize_total(adata, target_sum=1e6, inplace=True)
# ## log1p transform
# sc.pp.log1p(adata)
# %%% ---- 2.2 PCA ----
# ## start time
# start = time.perf_counter()
# ## run PCA
# sc.pp.pca(adata)
# print("PCA successful!")
# ## write h5ad
# adata.write(xen_dir / "SCLC_PCA.h5ad")
# end = time.perf_counter()
# print(f"Total run time {end - start}")
# %%% ---- 2.3 neighbour enrichment ----
## start time
start = time.perf_counter()
## NN analysis
sc.pp.neighbors(adata)
print("Neighbour enrichment successful")
adata.write(xen_dir / "SCLC_neighb.h5ad")
print("Neighbour write successful")
end = time.perf_counter()
print(f"Total run time {end - start}")
# %%% ---- 2.4 compute UMAP ----
## start time
start = time.perf_counter()
## UMPAP analysis
sc.tl.umap(adata)
print("UMAP successful!")
adata.write(xen_dir / "SCLC_umap.h5ad")
print("UMAP write successful")
end = time.perf_counter()
print(f"Total run time {end - start}")
# %%% ---- 2.5 compute leiden clusters ----
## start time
start = time.perf_counter()
## leiden analysis
sc.tl.leiden(adata)
print("Leiden successful! Pre-processing complete.")
# %% ---- 3.0 write adata ----
'''
Write the ad file containing all files and embedded clusters, etc. 
'''
adata.write(xen_dir / "SCLC_pp.h5ad")
end = time.perf_counter()
print(f"Total run time {end - start}")
