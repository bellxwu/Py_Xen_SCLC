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
# %% ---- 1.1 Load anndata ----
xen_dir = xc.xen_bwu
adata = ad.read_h5ad(xen_dir / 'SCLC_filtered.h5ad')
# %% --- 2.0 check labels ----
adata.obs['sample'].unique()
adata.obs.columns
adata.obsm['spatial'].shape
