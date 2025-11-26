#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 13:42:47 2025

Description: Loading and quality control testing of Xenium data outputs

@author: bellwu
"""
# %% ---- 1.0 Set up environment ----

import os
import xen_config as xc
import spatialdata as sd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

dir_64312 = xc.xen_dir / "SCLC_64312"
os.chdir(dir_64312)
zarr_64312 = dir_64312 / "SCLC_64312.zarr"

# %% ---- 2.0 Convert Xenium to SpatialData ----

sdata = sd.read_zarr(zarr_64312)
adata = sdata.tables["table"]
a_10 = adata.obs.head(10)
v_10 = adata.var.head(10)
# %% ---- 3.0 Preparing dataset ----

adata.var_names_make_unique()
sc.pp.calculate_qc_metrics(adata, percent_top = (10, 20, 50, 150), 
                           inplace = True)
adata.obs.columns
adata.var.columns
# %% ---- 4.0 Identifying background and negative probes ----

cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
) # control probes as percentage of total gene counts
cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
)
print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")
# %% ---- 5.0: Plotting QC metrics ----
'''
After calculating QC metrics with scanpy, need to then visualize first by plotting.
Metrics obtained can be categorized into "cell-based" or "transcript-based".
Cell-based involve all qc related to the quality of the cell
Transcript-based involve all metrics regarding the gene counts
'''
# %% ---- 5.1: Histogram plot of cell-based metrics ----
# create an empty panel of subplots 
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
# in first plot, show total transcripts per cell
axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)
# in second plot, show unique transcripts per cell
axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)
# in third plot, show approximately the size o
axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata.obs["cell_area"],
    kde=False,
    ax=axs[2],
)
# in final plot, show the ratio of the nucleus to the cell area
# this should tend to be smaller due to the nature of the SCLC cancer, 
axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,
    ax=axs[3],
)
fig.savefig(xc.qc_workdir / "qc_metric_cell.png", dpi = 300)
# %% ---- 5.2: Histogram of pct counts in genes ----
# create an empty panel of subplots 
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
# in first plot, show total transcripts per cell
axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["pct_counts_in_top_genes"],
    kde=False,
    ax=axs[0],
)
# in second plot, show unique transcripts per cell
axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)
# in third plot, show approximately the size o
axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata.obs["cell_area"],
    kde=False,
    ax=axs[2],
)
# in final plot, show the ratio of the nucleus to the cell area
# this should tend to be smaller due to the nature of the SCLC cancer, 
axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,
    ax=axs[3],
)
fig.savefig(xc.qc_workdir / "qc_metric_cell.png", dpi = 300)
# %% ---- 5.3: Histogram of gene-based metrics ----
# create empty subplots
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
# in first plot, show total transcripts per cell
axs[0].set_title("Total counts of each transcript")
sns.histplot(
    adata.var["total_counts"],
    kde=False,
    ax=axs[0],
)
# in second plot, show unique transcripts per cell
axs[1].set_title("Unique cells per transcript")
sns.histplot(
    adata.var["n_cells_by_counts"],
    kde=False,
    ax=axs[1],
)
# in third plot, show approximately the size o
axs[2].set_title("Mean expression over all cells")
sns.histplot(
    adata.var["mean_counts"],
    kde=False,
    ax=axs[2],
)
# in final plot, show the ratio of the nucleus to the cell area
# this should tend to be smaller due to the nature of the SCLC cancer, 
axs[3].set_title("Percentage cells feature is not found in")
sns.histplot(
    adata.var["pct_dropout_by_counts"],
    kde=False,
    ax=axs[3],
)
fig.savefig(xc.qc_workdir / "qc_metric_transcript.png", dpi = 300)
