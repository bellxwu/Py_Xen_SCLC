#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 12:18:11 2025

Description: Converting xenium outputs into a .Zarr file 

@author: bellwu
"""

# %% ---- 1.0 Set up environment ----

from pyxenium import xen_config as xc
import spatialdata_io as si

# %% Writing .Zarr from xenium output
xen_dir = xc.xen_bwu
# %% .Zarr for 23169
dir_23169 = xc.dir_23169
sdata = si.xenium(path = dir_23169,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_23169.h5ad")

sdata.write(xen_dir / "SCLC_23169.zarr", overwrite=True)
# %% .Zarr for 66144
dir_66144 = xc.dir_66144
sdata = si.xenium(path = dir_66144,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_66144.h5ad")

sdata.write(dir_66144 / "SCLC_23169.zarr", overwrite=True)
# %% .Zarr 
dir_19110 = xc.dir_19110
sdata = si.xenium(path = dir_19110,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_19110.h5ad")

sdata.write(dir_19110 / "SCLC_19110.zarr", overwrite=True)

