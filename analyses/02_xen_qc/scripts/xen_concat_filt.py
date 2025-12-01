#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 14:48:04 2025

Description: Merging anndata files and filtering out low-quality cells.

@author: bellwu
"""
# %% ---- 1.0 Set up environment ----
from pyxenium import xen_config as xc
import anndata as ad
import scanpy as sc
# %% ---- 1.1 set directory and load h5ad files ----
xen_dir = xc.xen_bwu
# load ad files via dictionary comprehension
adatas = {
    p.stem: ad.read_h5ad(p)
    for p in xen_dir.glob("*.h5ad")}
# merge files
adata = ad.concat(adatas, label='sample')
# %% ---- 2.0 filtering anndata ---- 
'''
This section describes filtering out low-quality cells from the merged anndata
file. Strategy is taken from xenium_qc_filtering.py 
'''
# %% ---- 2.1 preparing dataset ----
adata.var_names_make_unique()
sc.pp.calculate_qc_metrics(adata, percent_top = (10, 20, 50, 150), 
                           inplace = True)
# %% ---- 2.2 cell-based parameter: via min n genes per cell  ----
s = adata.obs["n_genes_by_counts"]
s.describe() # retrieve descriptive statistics 
s.quantile(.10) # 10% quartile 
# %% ---- 2.3 cell-based parameter: via total counts per cell ----
t = adata.obs["total_counts"] # list 
t.describe()
t.quantile(.10)
# %% ---- 2.4 filtering out cells ----
print("Before:", adata.shape)
init = adata.shape[0]
# filter cells
sc.pp.filter_cells(adata, min_counts=t.quantile(.1))
sc.pp.filter_cells(adata, min_genes=s.quantile(.1))

print("After:", adata.shape)
final = adata.shape[0]

cell_rm = init - final # how many cells removed
print("Number cells removed:", cell_rm)

# %% ---- 2.5 save files ----
# create directory for storage
xc.qc_workdir.mkdir(parents=True, exist_ok=True) 
# write .txt file
with open(xc.qc_workdir / "Filtered_anndata.txt", "w") as f:
    f.write(f"Number of cells removed: {cell_rm}")

# write new adata file
adata.write(xen_dir / "SCLC_filtered.h5ad")
