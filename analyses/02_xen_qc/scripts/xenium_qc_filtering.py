#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 16:21:18 2025

Description: Filtering out the problematic cells identified 

@author: bellwu
"""

# %% ---- 1.0 Set up environment ----

import os
import xen_config as xc
import spatialdata as sd
import scanpy as sc
import pandas as pd
import seaborn as sns

dir_64312 = xc.xen_dir / "SCLC_64312"
os.chdir(dir_64312)
zarr_64312 = dir_64312 / "SCLC_64312.zarr"

# %% ---- 2.0 Convert Xenium to SpatialData ----
sdata = sd.read_zarr(zarr_64312)
adata = sdata.tables["table"]

# %% ---- 3.0 Preparing dataset ----
adata.var_names_make_unique()
sc.pp.calculate_qc_metrics(adata, percent_top = (10, 20, 50, 150), 
                           inplace = True)

# %% ---- 4.0 Identifying filter parameters ----
'''
Section includes identifying the quantitative properties within qc features 
visualized with problematic qualities.
'''
adata.obs.columns # list all metrics
# %% ---- 4.1 cell-based parameter: cells with x number of genes in each ----
s = adata.obs["n_genes_by_counts"]
s.describe() # retrieve descriptive statistics 
s.quantile(.10) # 10% quartile 
# output if remove all cells with gene counts 1 std away
sd_1 = s.mean() - s.std()
(s < s.quantile(.10)).sum()
'''
Will need to remove out the low-quality cells. Thinking of anything below the 
10% quartile will filter out
'''
# %% ---- 4.2 cell-based parameter: total counts per cell ----
t = adata.obs["total_counts"] # list 
t.describe()
t.quantile(.10)
(t < 50).sum()
# %% ---- 4.3 cell-based parameter: nuclear ratio ----
n_c = adata.obs['nucleus_area'] / adata.obs['cell_area']
(n_c == 1).sum() # 2008 cells with area = nucleus

# identify how many transcripts expressed in these cells
nuc_ratio_1 = n_c == 1
nuc_ratio_1 = s[nuc_ratio_1]

# select out index
to_df = adata[nuc_ratio_1.index]

# convert to dataframe
df_nc = pd.DataFrame(
    to_df.X.toarray(),
    index = nuc_ratio_1.index,
    columns = adata.var_names)

# calculate and plot row sums
df_nc.sum(axis = 0).describe()
sns.histplot(
    df_nc.sum(axis = 0),
    kde=False
)
'''
Hard to interpret what these descriptive statistic really means for the cells
with high nuclear ratios. Will need to have annotated the cell types first but
this is good practice in indexing the anndata object.
'''
# %% ---- 4.4 gene-based parameters: general parameters ----
adata.var.columns # list all metrics
g = adata.var["n_cells_by_counts"]
g.describe()
g.quantile(.10)
# what genes appear less than 10 times?
g[g < 10]
'''
Likely won't filter out genes as can remove analysis for rare populations.
Bc is fluorescent probe-based, filtering for readouts above background was already
done for Xenium output
'''
# %% ---- 5.0 Filter out via cell-based QC metrics ----
'''
Filter out cells with properties in the bottom 10% quartile
'''
sc.pp.filter_cells(adata, min_counts=t.quantile(.1)) 
sc.pp.filter_cells(adata, min_genes=s.quantile(.1))
adata.shape
# %% ---- 6.0 Save filtered anndata file ----
adata.write(dir_64312 / "64312_filtered.h5ad")