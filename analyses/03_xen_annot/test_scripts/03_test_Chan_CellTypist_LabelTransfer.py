"""
Created on 2026.01.14

Description: Testing how to train the model with Chan et al. data for label transfer
to SCLC Xenium data.

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
# %% ---- 2.0 Load in preprocessed Chan et al. data ----
## directory variables
Chan_dir = xc.Chan_bwu
xen_dir = xc.xen_bwu
## reading h5ad file
adata_CxG = sc.read_h5ad(Chan_dir / "scChan_2021_CellxGene_combined.h5ad")
xen_panel = pd.read_csv(xen_dir / "xen_5k_panel.txt", sep="\t")
adata_SCLC = sc.read_h5ad(xen_dir / "64312_filtered.h5ad")
adata_1e4 = sc.read_h5ad(xen_dir / "SCLC_64312_pp.h5ad") # load in recently preprocessed data
# %% ---- 3.0 Prepare data for CellTypist ----
## select genes only within panel
mask = adata_CxG.var['feature_name'].isin(xen_panel['gene_name'])
adata_CxG_5k = adata_CxG[:,mask]
# %% ---- 4.0 Train CellTypist model ----
t_start = time.time()
model_CxG_5k = celltypist.train(X=adata_CxG_5k.X,
                                labels=adata_CxG_5k.obs['cell_type_fine'],
                                genes=adata_CxG_5k.var['feature_name'],
                                check_expression=False,
                                n_jobs=1,
                                max_iter=1)
t_end = time.time()
print(f"Time taken to train CellTypist model: {(t_end - t_start)/60} minutes")
## save trained model
model_CxG_5k.write(xen_dir / "CellTypist_Chan2021_CxG_5k_model.pkl")
# %% ---- 5.0 Load in trained model and test prediction ----
'''
Testing prediction on one SCLC sample
'''
# %% # normalize adata_SCLC
adata_SCLC.layers["counts"] = adata_SCLC.X.copy()
sc.pp.normalize_total(adata_SCLC, target_sum=1e4, inplace=True)
sc.pp.log1p(adata_SCLC)
# %% # check if normalized properly
## for recently loaded adata_1e4
X = adata_1e4.X
X_exp = X.expm1()

totals = np.array(X_exp.sum(axis=1)).ravel() # convert to np.array and then flatten with to 1D
print(np.median(totals))
# %% # load and run model
model_path = xen_dir / "CellTypist_Chan2021_CxG_5k_model.pkl"
t_start = time.time()
predictions = celltypist.annotate(adata_SCLC,
                                  model = str(model_path),
                                  majority_voting=False)
t_end = time.time()
print(f"Time elapsed for CellTypist prediction: {(t_end - t_start)} seconds")
# %% -- 6.0 Probing prediction object ----
type(predictions)
predictions.predicted_labels
prediction_adata = predictions.to_adata()
# %% testing ravel
X = csr_matrix([[1, 0], [0, 2]])
X.toarray().ravel()