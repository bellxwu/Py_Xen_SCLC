#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 15:33:10 2025

Description: Testing out and going through tutorials on SpatialData. 
Emphasize on the imaging aspect of the SpatialData object

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
'''
Main takeaway is to use read_zarr function and there elements within to explore
'''
# loading .zarr files
sdata = sd.read_zarr(zarr_64312)
adata = sdata.tables["table"]
adata.obs
adata.obsm['spatial']

# %% Indexing through DataTrees (xarray indexing)
sdata["morphology_focus"] # list DataTree
sdata["morphology_focus"].children # outputs all children of obj
sdata["morphology_focus"]["scale4"] # index into scale Group
sdata["morphology_focus/scale4"] # returns error
# can workaround by first selecting the DataTree
# %% Exploring the Dataset data type

# select out dataset
im_tree = sdata["morphology_focus"] # create 
scale4 = im_tree["scale4"].ds # select dataset
scale1 = im_tree["scale1"].ds
# explore dimensions
scale4.dims
scale1.dims

# exploring the coordinates
scale4.coords['x'] #  x coordinate; numerical array, float64
scale4.coords['y'] #  y coordinate; numerical array, float64
scale4.coords['c'] # list of channels
# can index through coordinates
scale4.isel(x=0) # this selects all 5 channels and all 5 y data values
scale4.isel(y=0) # all x values across all 5 channels
scale4.isel(c=0) # all x and y values from the first channel

# exploring data variables
scale4.data_vars # list data variable names
# %% Plotting with SpatialData_plot
# plot with low resolution
s4_fig = (sdata.pl.render_images(element = "morphology_focus",
                       scale = "scale4",
                       channel = "DAPI")) # plots the image from the DAPI channel
s4_fig.pl.show()
# %% plot higher resolution 
s1_fig = (sdata.pl.render_images(elements = "morphology_focus",
                        scale = "scale1",
                        channel = "DAPI"))
s1_fig.pl.show()
# %% saving image in lower scale
s4_fig.pl.show(save = xc.exp_workdir / "scale4.png",
               dpi = 600,
               figsize = (6, 8))
# %% saving image in higher scale
s1_fig.pl.show(save = xc.exp_workdir / "scale1.png",
               dpi = 600,
               figsize = (6, 8))
# %% Exploring labels and shapes variables
# labels
cell_lab = sdata["cell_labels"]
nuc_lab = sdata["nucleus_labels"]
# shapes
cell_boun = sdata['cell_boundaries']
cell_cir = sdata['cell_circles'] # lightweight cell segment
nuc_boun = sdata['nucleus_boundaries']
# %% overlaying images
s4_fig = (sdata.pl.render_images(element = "morphology_focus",
                       scale = "scale4",
                       channel = "DAPI"))
s4_lab = sdata.pl.render_labels(element = "nucleus_labels",
                                scale = "scale4",
                                channel = "DAPI")
cell_shap = sdata.pl.render_shapes(element = "cell_circles",
                                   fill_alpha = 0,
                                   outline_alpha = 0.5)
# %% show each individually
s4_fig.pl.show()
s4_lab.pl.show()
cell_shap.pl.show()
# %% overlay plots
(sdata
 .pl.render_images(
     element = "morphology_focus",
     scale = "scale4",
     channel = "DAPI")
 .pl.render_labels(
     element = "nucleus_labels",
     scale = "scale4",
     channel = "DAPI")
 .pl.render_shapes(
     element = "cell_circles",
     fill_alpha = 0,
     outline_alpha = 0.5)
 ).pl.show(
     save = xc.exp_workdir / "s4_overlay.png",
     dpi = 600,
     figsize = (6, 8))




 







