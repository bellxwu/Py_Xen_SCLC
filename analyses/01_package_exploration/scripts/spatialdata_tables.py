#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 13:00:31 2025

Description: Testing out the spatialdata non-image objects in tables. 

@author: bellwu
"""

# %% Set up environment

import os
import xen_config as xc
import spatialdata as sd
import anndata
import squidpy as sq
import xarray
import spatialdata_plot

dir_64312 = xc.xen_dir / "SCLC_64312"
os.chdir(dir_64312)
zarr_64312 = dir_64312 / "SCLC_64312.zarr"

# %% Convert Xenium to SpatialData
# loading .zarr files
sdata = sd.read_zarr(zarr_64312)

# %% Exploring the transcripts within the Points SpatialElement
points = sdata['transcripts']
points.compute()
p_10 = points.head(10)
points.columns
# %% Exploring the anndata table
# selecting the anndata object
adata = sdata.tables["table"]
# index value for the cell_ids
adata.obs_names 
# the data within each cell
a_10 = adata.obs.head(10) 
# gene names
adata.var_names
# the data for each gene
adata.var
# to access unstructured metadata
adata.uns 

# taking a look at the transcript counts in the first two cells
exp_2 = adata.X[0:2, :]
arr_2 = exp_2.toarray()





