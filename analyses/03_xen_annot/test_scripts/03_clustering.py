#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 11:19:12 2025

Description: Analyze the quality of dimensionality reduction and clustering 
for preprocessed SCLC samples.

@author: bellwu
"""
# %% ---- 1.0 set up local environment ----
import scanpy as sc
from pyxenium import xen_config as xc

xen_dir = xc.xen_bwu
pp_path = xen_dir / "SCLC_pp.h5ad"
SCLC_Preprocessed_Xen = sc.read(pp_path)
# %% ---- 2.0 Look through output ----
'''
Goal here is to familiarize myself with the xenium outputs. 
'''
SCLC_Preprocessed_Xen.shape
## after PCA
SCLC_Preprocessed_Xen.uns['pca']['variance'].shape
SCLC_Preprocessed_Xen.obsm['X_pca'].shape
## after neighbor
SCLC_Preprocessed_Xen.obsp['distances'].shape
SCLC_Preprocessed_Xen.obsp['connectivities'].shape
SCLC_Preprocessed_Xen.uns['neighbors']
## after UMAP
SCLC_Preprocessed_Xen.obsm['X_umap'].shape
## after leiden
SCLC_Preprocessed_Xen.obs['leiden']
# %% ---- 3.0 Initial plots ----
SCLC_Preprocessed_Xen.obs.columns
## directory to save plots
xc.annot_workdir.mkdir()
save_directory = xc.annot_workdir
## set directory settings
sc.settings.figdir = save_directory
## plotting the UMAP
fig = sc.pl.umap(SCLC_Preprocessed_Xen, 
                 color = ['leiden','sample'],
                 show = False,
                 return_fig = True)
fig.savefig(save_directory / "UMAP_leiden_sample.pdf", bbox_inches="tight")

