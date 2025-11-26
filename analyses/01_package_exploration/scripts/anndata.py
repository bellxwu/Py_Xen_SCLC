#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 19:28:43 2025

Description: Purpose of this script to learn the scverse core packages and
data structures. 

@author: bellwu
"""

# %% --- Set directory and environment ---

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

# %% --- Follow tutorial ---

# creating anndata object
rng = np.random.default_rng(seed=12345) # set seed for reproducibility
counts = csr_matrix(rng.poisson(1, size=(100, 2000)), dtype=np.float32)
print(counts)
counts[0,0]
adata = ad.AnnData(counts)
adata

# indexing
adata.obs_names
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
adata["Cell_1", "Gene_5"] # indexing gives a view 
adata[["Cell_1", "Cell_2"],["Gene_3", "Gene_6"]] # label based
adata[[1,4,6]][0:2] # fancy indexing

# adding aligned metadata
ct = rng.choice(["B", "T", "Monocyte"], size = (adata.n_obs,))
ct.shape
adata.obs["cell_type"] = pd.Categorical(ct) # Convert to categorical for efficiency
adata.obs # lists the cell with corresponding cell type
adata.obs.cell_type == "B" # can use the metadata label for boolean mask
st = rng.choice(["SCLC-A", "SCLC-B", "SCLC-C", "SCLC-D"], size = (adata.n_obs))
adata.obs["subtype"] = pd.Categorical(st) # can add multiple columns
'''
Although can add multiple columns, the addition of metadata has to be 1D array.
Alternative is using .obsm or .varm for storage of higher dimension metadata.
'''
# does anndata automatically create a copy?
test = adata[0:4, 0:4]
print(test.X) # prints values of csr
test.X[0,0] = 900 # automatically creates a copy
test.X[0,0] 
adata.X[0,0]
print(adata.X[0:4, 0:4]) 
print(test.X[0:4, 0:4])
