#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 30 15:50:47 2025

Description: Testing indexing of .zarr files with comparison to .h5ad files.

@author: bellwu
"""

# %% ---- 1.0 Set up environment ----
import xen_config as xc
import spatialdata as sd
import anndata as ad
import zarr
import numpy as np
import pandas as pd

# %% ---- 2.0 Loading files ----
zarr_23169 = xc.xen_bwu / "SCLC_23169.zarr"
sdata = sd.read_zarr(zarr_23169)
adata = ad.read_h5ad(xc.xen_bwu / "SCLC_23169.h5ad")
# %% ---- 3.0 Validating integrety of anndata conversion ----
root = zarr.open(zarr_23169, mode="r")
root.tree()
root['tables']['table']
# compare zarr with h5ad
zarr_23169.tree()
# %% testing pd.merge / concat
a = pd.DataFrame({
    "a": np.random.choice(100,5),
    "b": np.random.choice(100,5)})
b = pd.DataFrame({
    "a": np.random.choice(100,5),
    "b": np.random.choice(100,5)})
pd.merge(a, b, how='outer')
pd.concat((a,b), axis=0) # joins rows
pd.concat((a,b), axis=1) # joins columns, duplicate names can exist

