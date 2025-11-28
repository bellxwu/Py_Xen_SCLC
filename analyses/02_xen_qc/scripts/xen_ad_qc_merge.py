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
import spatialdata as sd
import zarr
import scanpy as sc
import pandas as pd
import seaborn as sns
# %% ---- set directory ----
xen_dir = xc.xen_bwu
zarr_23169 = xen_dir / "SCLC_23169.zarr"
# %% ---- 2.0 Validating integrety of anndata conversion ----
# load anndata
ad_23169 = ad.read_h5ad(xen_dir / "SCLC_23169.h5ad")
# load zarr
sd_23169 = sd.read_zarr(zarr_23169)
root = zarr.open_group(zarr_23169, mode="r")
print(root.tree())
print(root.attrs)

# check column metadata
ad_19110.obs.columns # check column metadata
ad_19110.var.columns
df = ad_19110.X[:5].toarray()


