#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 12:55:02 2025

Description: Check the integrity of the filter anndata file.

@author: bellwu
"""
# %% ---- 1.0 set up environment ----
import anndata as ad
from pyxenium import xen_config as xc
import scanpy as sc
# %% ---- 1.1 Load anndata ----
xen_dir = xc.xen_bwu
adata = ad.read_h5ad(xen_dir / 'SCLC_filtered.h5ad')
ad_19110 = ad.read_h5ad(xen_dir / 'SCLC_19110.h5ad')
ad_23169 = ad.read_h5ad(xen_dir / "SCLC_23169.h5ad")
# %% testing concat
adatas = dict({'ad_19110': ad_19110,
               'ad_23169': ad_23169})
adata = ad.concat(adatas, label='samples')
adata.obs['samples']
sc.pp.calculate_qc_metrics(adata, percent_top = (10, 20, 50, 150), 
                           inplace = True)
# %% --- 2.0 check labels ----
adata.obs['sample'].unique()
adata.obs.columns
adata.obsm['spatial'].shape

