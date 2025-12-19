#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 20:30:06 2025

Description: Pre-processing pipeline for individual samples from bulk file.
Purpose is to avoid batch effects seen originally with neighborhood analysis.

@author: bellwu
"""
# %% ---- 1.0 Setting up environment ----
import scanpy as sc
import anndata as ad
import time
from pyxenium import xen_config as xc

## import dir from h4h cluster
xen_dir = xc.bell_concat 
preprocessed_dir = xc.bell_pp

## reading anndata
SCLC_xen_path = xen_dir / "SCLC_xen.h5ad"
SCLC_xen = ad.read_h5ad(SCLC_xen_path)
# %% ---- 2.0 Filtering for sample ----
'''
Filter out samples within the bulk data and create a copy within a dictionary
'''
## boolean mask to filter out one sample
SCLC_sample_names = SCLC_xen.obs['sample'].unique()
## dictionary comprehension to create mask
sample_masks = {i: SCLC_xen.obs['sample'] == i for i in SCLC_sample_names}
## select samples from mask
SCLC_samples = {i: SCLC_xen[sample_masks.get(i)] for i in SCLC_sample_names}
# %% ---- 3.0 Normalization and log1p transformation ----
'''
Preprocessing done iteratively through all SCLC samples
'''
# %% ---- 3.1 Sample preprocessing ----
## write loop across all samples
for SampleName, SampleAnnData in SCLC_samples.items():
    '''
    For looop for running the SCLC xenium samples. Input is given as a 
    dictionary with SampleName as the key and SampleAnnData as its values.
    '''
    InitialStart = time.perf_counter()
    # initial normalization
    print(f"Starting analysis of {SampleName}")
    if "counts" not in SampleAnnData.layers:
        SampleAnnData.layers["counts"] = SampleAnnData.X.copy()
    sc.pp.normalize_total(SampleAnnData, target_sum=100, inplace=True)
    sc.pp.log1p(SampleAnnData)
    
    # PCA across samples
    StartTime = time.perf_counter()
    sc.pp.pca(SampleAnnData,
              n_comps=(min(SampleAnnData.shape)-1)) # compute all PCs
    EndTime = time.perf_counter()
    print(f"\nPCA of {SampleName} successful, total time {EndTime - StartTime}")
    print(f"Number of PCs computed {SampleAnnData.obsm['X_pca'].shape}")
    # write pca .h5ad
    StartTime = time.perf_counter()
    SampleAnnData.write(f"{SampleName}_pp.h5ad")
    EndTime = time.perf_counter()
    print(f"Writing .h5ad of {SampleName} successful, total time {EndTime - StartTime}")
    
    # neighbors analysis across samples
    StartTime = time.perf_counter()
    sc.pp.neighbors(SampleAnnData,
                    n_pcs=SampleAnnData.obsm['X_pca'].shape[1],
                    n_neighbors=16)
    EndTime = time.perf_counter()
    print(f"\nNeigbors of {SampleName} successful, total time {EndTime - StartTime}")
    # write neighbors .h5ad
    StartTime = time.perf_counter()
    SampleAnnData.write(f"{SampleName}_pp.h5ad")
    EndTime = time.perf_counter()
    print(f"Writing .h5ad of {SampleName} successful, total time {EndTime - StartTime}")
    
    # UMAP analysis
    StartTime = time.perf_counter()
    sc.tl.umap(SampleAnnData)
    EndTime = time.perf_counter()
    print(f"\nUMAP of {SampleName} successful, total time {EndTime - StartTime}")
    # write neighbors .h5ad
    StartTime = time.perf_counter()
    SampleAnnData.write(f"{SampleName}_pp.h5ad")
    EndTime = time.perf_counter()
    print(f"Writing .h5ad of {SampleName} successful, total time {EndTime - StartTime}")
    
    # Leiden analysis
    StartTime = time.perf_counter()
    sc.tl.leiden(SampleAnnData,
                 flavor="igraph",
                 n_iterations=2,
                 directed=False)
    EndTime = time.perf_counter()
    print(f"\nLeiden of {SampleName} successful, total time {EndTime - StartTime}")
    # write neighbors .h5ad
    StartTime = time.perf_counter()
    SampleAnnData.write(f"{SampleName}_pp.h5ad")
    EndTime = time.perf_counter()
    print(f"Writing .h5ad of {SampleName} successful, total time {EndTime - StartTime}")
    print(f"\n Total time taken to process sample {EndTime - InitialStart}")















