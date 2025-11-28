#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 28 15:36:32 2025

Description: Script to read and write .h5ad and .zarr file from SCLC sample
23169. The goal having the .Zarr and .h5ad is to check non-corrupt conversion.

@author: bellwu
"""
# %% Set up environment
from pyxenium import xen_config as xc
import spatialdata_io as si
# %% Setting up directory
xen_dir = xc.xen_bwu
# %% .h5ad for 23169
dir_23169 = xc.dir_23169
sdata = si.xenium(path = dir_23169,
                 cell_boundaries=True,
                 nucleus_boundaries=True,
                 morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
# write h5ad
adata.write(xen_dir / "SCLC_23169.h5ad")
# write zarr
sdata.write(xen_dir / "SCLC_23169.zarr", overwrite=True)