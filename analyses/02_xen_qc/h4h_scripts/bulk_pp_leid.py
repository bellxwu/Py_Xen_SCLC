#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 08:40:44 2025

Description: Perform leiden clustering on h5ad SCLC xenium dataset.

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
h5ad_path = xen_dir / "SCLC_umap.h5ad"
adata = ad.read_h5ad(h5ad_path)
# %% ---- 2.0 Leiden clustering ----
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
## start time
start = time.perf_counter()
## leiden analysis
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)
print("Leiden successful! Pre-processing complete.")
# %% ---- 3.0 write adata ----
'''
Write the ad file containing all files and embedded clusters, etc. 
'''
adata.write(xen_dir / "SCLC_pp.h5ad")
end = time.perf_counter()
print(f"Total run time {end - start}")

