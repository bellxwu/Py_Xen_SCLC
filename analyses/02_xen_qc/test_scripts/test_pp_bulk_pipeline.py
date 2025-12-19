#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 13:01:39 2025

Description: Testing the preprocessing pipeline with generated small anndata
objects to see if pipeline will run.

@author: bellwu
"""
# %% ---- 1.0 Setting up environment ----
import scanpy as sc
import anndata as ad
import time
from pyxenium import xen_config as xc
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

## import Paths
xen_dir = xc.xen_bwu
SCLC_xen_path = xen_dir / "SCLC_xen.h5ad"
SCLC_xen = ad.read_h5ad(SCLC_xen_path)
test_dir = xc.test_bwu
# %% Creating a generator for random matrix for adata
## function to create random matrix
def gen_counts(times, dim=(100,1000)):
    rng = np.random.default_rng(seed=12345) # set seed for reproducibility
    i = 1
    while i <= times:
        yield csr_matrix(rng.poisson(1, dim), dtype=np.float32) 
        i += 1
## create two random matrix counts        
c1, c2 = gen_counts(2)
## 
test_groups = {'SCLC_64312': c1, 'SCLC_4462962': c2 }
test_groups = {k: ad.AnnData(v) for k, v in test_groups.items()}
# %% testing the write pipeline
## counter
start = time.perf_counter()
end = time.perf_counter()
# %% for loop
for k, v in test_groups.items():
    # for PCA
    start = time.perf_counter()
    sc.pp.pca(v, n_comps=(v.shape[0]-1))
    end = time.perf_counter()
    print(f"\nTotal PCA run time for {k} is {end - start}")
    print(f"Number of pcas computed {v.obsm['X_pca'].shape}")
    start = time.perf_counter()
    v.write(test_dir / f"{k}_pp.h5ad")
    end = time.perf_counter()
    print(f"Total PCA write time for {k} is {end - start}")
    
    # for neighbours
    start = time.perf_counter()
    sc.pp.neighbors(v, 
                    n_pcs=v.obsm['X_pca'].shape[1])
    end = time.perf_counter()
    print(f"Neighbors runtime is {end - start}")
    
    # leiden
    start = time.perf_counter()
    sc.tl.leiden(v,
                 flavor="leidenalg",
                 n_iterations=2)
    end = time.perf_counter()
    print(f"Leiden runtime is {end - start}")
    
# %% check how many pcs
SCLC_64312 = test_groups['SCLC_64312']
SCLC_64312.uns
print(f"\nTotal PCA run time for {SCLC_64312} is {end - start}")
print("test")
# %% 













