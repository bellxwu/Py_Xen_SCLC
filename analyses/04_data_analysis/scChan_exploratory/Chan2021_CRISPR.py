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
import os
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from pyxenium import plot_utils as px_plt
# %% ---- 2.0 Load in Chan et al. data ----
## directory variables
Chan_dir = xc.Chan_bwu
## reading h5ad file
adata_Chan = sc.read_h5ad(Chan_dir / "scChan_2021_CellxGene_combined.h5ad")
# %% # ---- 3.0 Check expression levels ----
# # check normalization of Chan data 
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
# %% ---- 3.1 Explore CRISPR hit genes ----
# list of hits
hits = ['NRK', 'TWNK', 'PEO1', 'SETD1B', 'PECAM1', 'SLC6A11', 'ASCL1']
# check how many hits in data
adata_Chan.var['feature_name'].isin(hits).sum()
# subset out hits
mask = adata_Chan.var['feature_name'].isin(hits)
adata_hits = adata_Chan[:, mask]
adata_hits.var['feature_name'] # check which hits present in data
# %% ---- 3.2 Check expression level of hits ----
## convert ENSEMBL to geneIDs
adata_hits.var_names = adata_hits.var['feature_name'].astype(str)
## look through cell-types
adata_hits.obs['cell_type_coarse'].unique()
# %% ---- 3.3 Exploratory analyses ----
## setting up directory
os.mkdir(xc.Chan_workdir)
ChanWorkingDir = xc.Chan_workdir
# %% plot violin for coarse cell types
sc.pl.violin(adata_hits,
             keys=["NRK", "PECAM1", "SETD1B", "SLC6A11", "ASCL1"],
             groupby="cell_type_coarse",
             jitter=0.4, 
             rotation=45, 
             use_raw=False,
             show=False)
plt.savefig(ChanWorkingDir / "ViolinPlotOfHits.png", dpi = 300)
plt.close()
# %% cell counts for each cell type
adata_Chan.obs.columns
adata_Chan.obs['cell_type_coarse'].value_counts()
# %% plot violin for fine cell types
sc.pl.violin(adata_hits,
             keys=["NRK", "PECAM1", "SETD1B", "SLC6A11", "ASCL1"],
             groupby="cell_type_fine",
             jitter=0.4, 
             rotation=45, 
             use_raw=False,
             show=False)
plt.savefig(ChanWorkingDir / "ViolinPlotOfHits_fine.png", dpi = 300)
plt.close()
# %% Plot for multiple genes
plots = px_plt.OrderedViolinPlots(
    adata=adata_hits,
    GeneList=["PECAM1", "NRK", "SETD1B", "SLC6A11"],
    ColumnOfInterest='cell_type_coarse',
    AddEmpty=False
)
plots_fine = px_plt.OrderedViolinPlots(
    adata=adata_hits,
    GeneList=["PECAM1", "NRK", "SETD1B", "SLC6A11"],
    ColumnOfInterest='cell_type_fine',
    AddEmpty=False
)
# %%
