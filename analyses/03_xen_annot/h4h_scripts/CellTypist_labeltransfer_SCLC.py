"""
Created on 2026.01.15

Description: Label transfer testing of SCLC Xenium data using CellTypist model trained on Chan et al. 2021 data.

@author: bellwu
"""
# %% ---- 1.0 set up local environment ----
import scanpy as sc
from pyxenium import xen_config as xc
import celltypist
import time
import gc
# %% ---- 2.0 load SCLC Xenium data ----
## directory variables
xen_dir = xc.bell_pp 
label_dir = xc.bell_label

## load all h5ad files in directory
adatas = {p.stem: sc.read_h5ad(p) for p in xen_dir.glob("*.h5ad")}
## load model 
model_dir = xc.xen_bwu / "CellTypist_Chan2021_CxG_5k_model.pkl"
Chan_CellTypist_model = celltypist.models.Model.load(str(model_dir))
# %% ---- 3.0 load and run model ----
predictions_dict = {}
for SampleName, SampleAnnData in adatas.items():
    # run predictions
    print(f"Running CellTypist prediction on sample: {SampleName}")
    t_start = time.time()
    predictions = celltypist.annotate(SampleAnnData,
                                      model = str(model_dir),
                                      majority_voting=True)
    t_end = time.time()
    print(f"Time elapsed for CellTypist prediction on {SampleAnnData}: {(t_end - t_start)} seconds")
    
    # save predictions anndata with labels 
    prediction_adata = predictions.to_adata()
    prediction_adata.write_h5ad(label_dir / f"{SampleName}_CellTypist_labeled.h5ad")
    print(f"Finished CellTypist labeling for sample: {SampleName}\n")
    del SampleAnnData
    gc.collect()
