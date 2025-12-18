#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 15:57:41 2025

Description: Generate new .h5ad file from previous h5ad with goal of removing
the duplicate samples 

@author: bellwu
"""
# %% set up environment
from pyxenium import xen_config as xc
import anndata as ad
# %% load files
xen_dir = xc.xen_bwu
filtered_path = xen_dir / "SCLC_dup.h5ad"
# %% read anndata file
SCLC_filtered = ad.read_h5ad(filtered_path)
# %% remove duplicated samples
to_rm = ['SCLC_19110', 'SCLC_290442']
adata_nodup = SCLC_filtered[~SCLC_filtered.obs['sample'].isin(to_rm), :]
# %% write the h5ad
adata_nodup.write(xen_dir / "SCLC_xen.h5ad")