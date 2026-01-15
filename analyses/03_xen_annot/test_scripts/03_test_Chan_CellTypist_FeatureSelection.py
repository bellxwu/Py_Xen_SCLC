"""
Created on 2026.01.11

Description: Probe through the Chan et al. single cell 10X dataset. Begin initial
testing of annotation approaches; specifically with CellTypist.

@author: bellwu
"""
# %% ---- 1.0 set up local environment ----
import scanpy as sc
from pyxenium import xen_config as xc
import pandas as pd
import celltypist
import time
from scipy.sparse import csr_matrix
import numpy as np
from pybiomart import Server
# %% ---- 2.0 Load in preprocessed Chan et al. data ----
## directory variables
Chan_dir = xc.Chan_bwu
xen_dir = xc.xen_bwu
## reading h5ad file
adata_Chan_SCLC = sc.read_h5ad(Chan_dir / "scChan_2021.h5ad")
adata_Chan_combined = sc.read_h5ad(Chan_dir / "scChan_2021_combined.h5ad")
adata_CxG = sc.read_h5ad(Chan_dir / "scChan_2021_CellxGene_combined.h5ad")
## reading 5k gene panel
gene_panel_5k = xen_dir / "xen_5k_panel.txt"
xen_panel = pd.read_csv(gene_panel_5k, sep="\t")
# %% ---- 3.0 Initial exploration of data ----
## of downloaded Chan SCLC data
adata_Chan_SCLC.shape
adata_Chan_SCLC.obs['clusters']
adata_Chan_SCLC.obs['clusters_fine']
adata_Chan_SCLC.var.columns
adata_Chan_SCLC.uns
adata_Chan_SCLC.obs_names

## of downloaded Chan data from CellxGene
adata_CxG.shape
adata_CxG.obs.columns
adata_CxG.var_names
adata_CxG.obs['cell_type_coarse'].unique()
adata_CxG.obs['cell_type_fine'].unique()
### checking csr_matrix values
adata_CxG.X[0:5,0:5].toarray()
adata_CxG.raw.X[0:5,0:5].toarray()
adata_CxG.layers['counts'] = adata_CxG.raw.copy()

## of downloaded Chan combined data
adata_Chan_combined.shape
adata_Chan_combined.obs.columns
adata_Chan_combined.obs['cell_type_coarse'].unique()
adata_Chan_combined.obs['cell_type_fine'].unique()
adata_Chan_combined.var_names
### checking csr_matrix values
adata_Chan_combined.X[0:5,0:5]
adata_Chan_combined.layers.keys()
'''
CellxGene combined dataset, since it has cell type annotations present but gene names are 
Ensembl IDs instead of gene symbols. The ensembl IDs are present as "feature_name" in var.

The Chan_SCLC only has cluster IDs but no reference to cell types given.

Chan_combined dataset contains the cell type annotations also and gene symbols but has no raw
'''
# %% ---- 4.0 Initial CellTypist testing ----
'''
Goal here is to start running CellTypist on Chan CellxGene data. Want to do all the testing
of the CellTypist algorithm to see how if it works. 
'''
# %% ---- 4.1 testing: merge and intersection functions ----
## merging
test_df = xen_panel.head(10)
test_df_1 = xen_panel.iloc[9:19]
test_df['gene_name']
test_df.index = test_df['gene_name']

test_intersection = test_df.index.intersection(xen_panel['gene_name'])
test_intersection = pd.merge(test_df, test_df_1, how='inner', on='gene_name')
## need to create a boolean mask for filtering
mask = test_df_1.isin(test_df['gene_name'])
test_df_1[mask] # why does this return NaN
test_df_1.loc[mask] # this doesn't work
## testing out the mask shape
len(test_df_1)
len(test_df)
len(mask)
type(mask) # mask is a DataFrame, same shape as original
mask.shape
## if I change the isin
mask2 = test_df_1['gene_name'].isin(test_df['gene_name'])
type(mask2) # this is a Series
test_df_1[mask2] # this works, does not return NaN
# %% ---- 4.2 merging xenium 5k panel with Chan geneset ----
type(adata_Chan_combined.var_names) # type of object = index
mask = adata_Chan_combined.var_names.isin(xen_panel['gene_name'])
len(mask)
adata_test = adata_Chan_combined[:,mask] # subsetted adata with 5k genes
## final adata with 5k genes only
# %% ---- 4.3 Testing feature selection CellTypist model on Chan data ----
'''
Downsample: No need to downsample here bc dataset is smaller than cells in Xenium set

HVGs: am thinking no need to do feature selection here since the gene panel is already
small. Two methods:
1) Want to see how labels looks with the 5k panel only with no feature selection
2) See how selecting for hvgs in Chan data and then filtering for these further

Will try out method 2 first. 
'''
# %% ---- 4.3.0 setting up adata for CellTypist ----
'''
Prepare the adata for CellTypist training. Need to set up layers for counts and normalized.
''' 
# %%  setting up layers for CellTypist
adata_CxG.layers['counts'] = adata_CxG.raw.X.copy()
adata_CxG.layers['normalized'] = adata_CxG.raw.X.copy()
# %% converting ensembl IDs to gene symbols
server = Server(host='http://www.ensembl.org')
mart = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
### get mapping dataframe
mapping_df = mart.query(attributes=['ensembl_gene_id', 'hgnc_symbol'])
### create dictionary for mapping to convert mappings
ens2sym_dict = dict(zip(mapping_df['Gene stable ID'], mapping_df['HGNC symbol']))
### convert adata var_names
adata_CxG.var['gene_symbol'] = adata_CxG.var_names.map(ens2sym_dict)

# %% ---- 4.3.1 normalize data ----
## normalize to 1e4 counts per cell
sc.pp.normalize_total(adata_CxG,
                      target_sum=1e4,
                      layer='normalized',
                      inplace=True)
## validate normalization
adata_CxG.layers['normalized'][0:5,0:5].toarray()
## log1p transform
sc.pp.log1p(adata_CxG, layer='normalized')
## validate log1p
adata_CxG.layers['normalized'][0:5,0:5].toarray()
# %% ---- 4.3.2 selecting HVGs in Chan data ----
## creating model for feature selection
t_start = time.time()
ModelFeatureSelection = celltypist.train(X=adata_CxG.layers['normalized'], 
                                         labels=adata_CxG.obs['cell_type_fine'],
                                         genes=adata_CxG.var['feature_name'],
                                         n_jobs=1,
                                         max_iter=3, 
                                         use_SGD=True)
t_end = time.time()
print(f"Time taken to train feature selection model: {t_end - t_start} seconds")
'''
Cannot run on desktop, need at least 10 CPUs and >32GB RAM. Will run on cluster.
'''
# %% ---- 4.3.3 extracting HVGs from model ----
## extracting gene indices from model
gene_index = np.argpartition(np.abs(ModelFeatureSelection.classifier.coef_), -100, axis=1)[:,-100:]
gene_index = np.unique(gene_index)
print(f"Number of HVGs selected: {len(gene_index)}")
# %% ---- 4.3.4 check how many of 5k genes within selected HVGs ----
## select out gene symbols from HVG indices and panel
hvg_genes = adata_CxG.var['feature_name'].values[gene_index]
panel_genes = xen_panel['gene_name'].values
## intersecting genes
intersect_genes = np.intersect1d(hvg_genes, panel_genes)
print(f"Number of genes in 5k panel that are also HVGs: {len(intersect_genes)}")
# %% ---- 4.3 Results ----
'''
Only 497 genes overlap between the 5k panel and the HVGs selected from Chan data. Very little and
might not be enough to do an accurate label transfer.

Might just use all 5000 genes. May be more computationally expensive?
'''