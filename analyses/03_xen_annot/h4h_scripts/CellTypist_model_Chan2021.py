"""
Created on 2026.01.15

Description: Train CellTypist model with Chan et al. 2021 CellxGene combined dataset for label
transfer to SCLC Xenium data.

@author: bellwu
"""
# %% ---- 1.0 set up local environment ----
import scanpy as sc
from pyxenium import xen_config as xc
import celltypist
import time
import pandas as pd
# %% ---- 2.0 Load in preprocessed Chan et al. data ----
Chan_dir = xc.Chan_bwu
xen_dir = xc.xen_bwu
## reading h5ad files
adata_Chan = sc.read_h5ad(Chan_dir / "scChan_2021_CellxGene_combined.h5ad")
xen_panel = pd.read_csv(xen_dir / "xen_5k_panel.txt", sep="\t")
# %% ---- 3.0 Prepare data for CellTypist ----
## set up layers for normalization
adata_Chan.layers['normalized'] = adata_Chan.raw.X.copy()
## normalize data
sc.pp.normalize_total(adata_Chan, target_sum=1e4, layer='normalized')
sc.pp.log1p(adata_Chan, layer='normalized')
## select genes only within panel
mask = adata_Chan.var['feature_name'].isin(xen_panel['gene_name'])
adata_Chan_5k = adata_Chan[:, mask]
## write anndata for future use
adata_Chan_5k.write_h5ad(xen_dir / "scChan2021_CellTypist_5k_normalized_genes.h5ad")
# %% ---- 4.0 Train CellTypist model ----
t_start = time.time()
model_Chan_5k = celltypist.train(X=adata_Chan_5k.layers['normalized'],
                                labels=adata_Chan_5k.obs['cell_type_fine'],
                                genes=adata_Chan_5k.var['feature_name'],
                                n_jobs=10,
                                max_iter=500,
                                check_expression=False)
t_end = time.time()
print(f"Time taken to train CellTypist model: {(t_end - t_start)/60} minutes")
## save trained model
model_Chan_5k.write(xen_dir / "CellTypist_Chan2021_CxG_5k_model.pkl")