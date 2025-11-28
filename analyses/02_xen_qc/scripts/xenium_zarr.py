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

# %% Writing .h5ad from xenium output
xen_dir = xc.xen_bwu
# %% .h5ad for 23169
# dir_23169 = xc.dir_23169
# sdata = si.xenium(path = dir_23169,
#                  cell_boundaries=True,
#                  nucleus_boundaries=True,
#                  morphology_mip=False,
#                   aligned_images=True,
#                   cells_as_circles=True)
# adata = sdata.tables['table']
# adata.write(xen_dir / "SCLC_23169.h5ad")

# sdata.write(xen_dir / "SCLC_23169.zarr", overwrite=True)
# %% .h5ad for 66144
dir_66144 = xc.dir_66144
sdata = si.xenium(path = dir_66144,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_66144.h5ad")

# sdata.write(dir_66144 / "SCLC_23169.zarr", overwrite=True)
# %% .h5ad for 19110
dir_19110 = xc.dir_19110
sdata = si.xenium(path = dir_19110,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_19110.h5ad")

<<<<<<< HEAD
# sdata.write(dir_19110 / "SCLC_19110.zarr", overwrite=True)
# %% .h5ad for 290442
dir_290442 = xc.dir_290442
sdata = si.xenium(path = dir_290442,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_dir_290442.h5ad")
# %% .h5ad for dir_4462962
dir_4462962 = xc.dir_4462962
sdata = si.xenium(path = dir_4462962,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_dir_4462962.h5ad")
# %% .h5ad for dir_64312
dir_64312 = xc.dir_64312
sdata = si.xenium(path = dir_64312,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_dir_64312.h5ad")
# %% .h5ad for dir_19110T1_2
dir_19110T1_2 = xc.dir_19110T1_2
sdata = si.xenium(path = dir_19110T1_2,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_dir_19110T1_2.h5ad")
# %% .h5ad for dir_290442_2
dir_290442_2 = xc.dir_290442_2
sdata = si.xenium(path = dir_290442_2,
                  cell_boundaries=True,
                  nucleus_boundaries=True,
                  morphology_mip=False,
                  aligned_images=True,
                  cells_as_circles=True)
adata = sdata.tables['table']
adata.write(xen_dir / "SCLC_dir_290442_2.h5ad")
=======
sdata.write(dir_19110 / "SCLC_19110.zarr", overwrite=True)

>>>>>>> 5a775120a8f2eeb17342e411c4337a7ac15a1be4
