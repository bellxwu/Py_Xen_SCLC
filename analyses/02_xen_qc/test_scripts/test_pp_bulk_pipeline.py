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
import numpy as np
from scipy.sparse import csr_matrix
import psutil
import os
import gc

# ## import Paths
# xen_dir = xc.xen_bwu
# SCLC_xen_path = xen_dir / "SCLC_xen.h5ad"
# SCLC_xen = ad.read_h5ad(SCLC_xen_path)
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
c1, c2 = gen_counts(2, dim=(1000,5000))
## 
test_groups = {'SCLC_64312': c1, 'SCLC_4462962': c2 }
test_groups = {k: ad.AnnData(v) for k, v in test_groups.items()}
# %% testing the write pipeline
## counter
start = time.perf_counter()
end = time.perf_counter()
min(c1.shape)-1
# %% testing memory with libraries
proc = psutil.Process(os.getpid()) # processing object
proc.memory_info() # returns tuple of memory information
n = proc.memory_info().rss # prints in bytes'
n
## note there are 1024 bytes in 1 "kibibyte"
def bytes2human(n):
    '''
    Converts raw bytes into readable string
    '''
    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 2 ** ((i + 1) * 10)
    for s in reversed(symbols):
        if abs(n) >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.2f%s' % (value, s) # change how many decimal points
    return "%sB" % n
n
bytes2human(n)
# %% for loop
for k, v in test_groups.items():
    # for PCA
    StartMem = proc.memory_info().rss
    start = time.perf_counter()
    sc.pp.pca(v, n_comps=(min(v.shape)-1))
    end = time.perf_counter()
    print(f"\nMemory before PCA {StartMem}")
    print(f"Total PCA run time for {k} is {end - start}")
    print(f"Number of pcas computed {v.obsm['X_pca'].shape}")
    start = time.perf_counter()
    v.write(test_dir / f"{k}_pp.h5ad")
    end = time.perf_counter()
    EndMem = proc.memory_info().rss
    print(f"Total PCA write time for {k} is {end - start}")
    print(f"Total memory after PCA {EndMem}")
    print(f"Memory consumed {EndMem - StartMem}")
    
    # for neighbours
    StartMem = proc.memory_info().rss
    start = time.perf_counter()
    print(f"\nMemory before neighbors analysis {StartMem}")
    sc.pp.neighbors(v, 
                    n_pcs=min(v.obsm['X_pca'].shape))
    end = time.perf_counter()
    EndMem = proc.memory_info().rss
    print(f"Neighbors runtime is {end - start}")
    print(f"Total memory after neighbors analysis {EndMem}")
    print(f"Memory consumed {EndMem - StartMem}")
    
    # leiden
    StartMem = proc.memory_info().rss
    start = time.perf_counter()
    print(f"\nMemory before leiden analysis {StartMem}")
    sc.tl.leiden(v,
                 flavor="igraph",
                 n_iterations=2)
    end = time.perf_counter()
    EndMem = proc.memory_info().rss
    print(f"Leiden runtime is {end - start}")
    print(f"Total memory after leiden analysis {EndMem}")
    print(f"Memory consumed {EndMem - StartMem}")

    # delete v to free memory
    StartMem = proc.memory_info().rss
    print(f"\nMemory before deleting v {StartMem}")
    del v
    gc.collect()
    EndMem = proc.memory_info().rss
    print(f"Total memory after deleting v {EndMem}")
    print(f"Memory freed {StartMem - EndMem}")
# %% check how many pcs
SCLC_64312 = test_groups['SCLC_64312']
SCLC_64312.uns
print(f"\nTotal PCA run time for {SCLC_64312} is {end - start}")
print("test")
# %% final memory calculation
b = 1121026048 - 781287424
bytes2human(b)
## the entire process consumed 324M of memory

## compared to memory usage at start, memory use has increased 
MemChange = proc.memory_info().rss - 781287424
bytes2human(MemChange)

## deleting variables to see if it improves memory
MemBefore = proc.memory_info().rss
del test_groups
MemAfter = proc.memory_info().rss
## see how much memory is saved
bytes2human(MemBefore - MemAfter) # only 160Kb of memory saved 
# %% Run garbage collect
MemBefore = proc.memory_info().rss
gc.collect()
MemAfter = proc.memory_info().rss
bytes2human(MemBefore - MemAfter)
# %% removing values from dictionary
test_dict = {'SCLC_64312': c1, 'SCLC_4462962': c2 }
del test_dict['SCLC_64312']












# %%
