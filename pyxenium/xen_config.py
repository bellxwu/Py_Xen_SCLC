"""
Description: Configuration script for all paths for Xenium python analysis.
Analysis pipeline edited for the location of the output files on the cluster. 

Created on Thu Nov  6 13:28:30 2025
@author: bellwu
"""

from pathlib import Path
from datetime import date

# --- Base Directory ---
base_dir = Path(__file__).resolve().parents[1]

# -- Data Directories --
# Cluster data directory (within lokgroup)
# xenium raw data on cluster
xen_data_lokgroup = Path("/cluster/home/bellwu/lokgroup/Xenium_runs")
xen_1 = xen_data_lokgroup / "20241212__190528__Lok_Jalal_241210"
xen_2 = xen_data_lokgroup / "20250328__182124__Lok_Jalal_250324"
# xenium personal data on cluster
lokgroup_bellwu = Path("/cluster/home/bellwu/lokgroup/bellwu/")
bell_xen = lokgroup_bellwu / "xen_data"
bell_concat = lokgroup_bellwu / "SCLC_concat"
bell_pp = lokgroup_bellwu / "SCLC_preprocessed"

# Personal data directory
data_bwu = base_dir / "all_data" 
xen_bwu = data_bwu / "xen_data"
test_bwu = data_bwu / "test_files"


# ---- Samples ----
# from first directory
dir_23169 = xen_1 / "output-XETG00082__0051610__23169__20241212__190621"
dir_66144 = xen_1 / "output-XETG00082__0051610__66144__20241212__190621"
dir_19110 = xen_1 / "output-XETG00082__0051621__19110__20241212__190621"
dir_290442 = xen_1 / "output-XETG00082__0051621__290442__20241212__190621"
# from second directory
dir_4462962 = xen_2 / "output-XETG00082__0057746__4462962__20250328__182222"
dir_64312 = xen_2 / "output-XETG00082__0057746__64312__20250328__182222"
dir_19110T1_2 = xen_2 / "output-XETG00082__0057749__19110T1_2__20250328__182222"
dir_290442_2 = xen_2 / "output-XETG00082__0057749__290442__20250328__182222"

# -- Timestamp string --
today_str = date.today().strftime("%Y-%m-%d")

# -- Analysis directories -- 
'''
Overall project layout is split into main analyses. 
1) Package_exploration: exploring how each package needed works
2) Quality control
3) Segmentation
'''
analyses_dir = base_dir / "analyses" # base directory for all results
pack_exp_dir = analyses_dir / "01_package_exploration" # subanalysis exploring packages
qc_dir = analyses_dir / "02_xen_qc"
annot_dir = analyses_dir / "03_xen_annot"

# ------- package_exploration directory -------
'''
Directory here contains scripts and results of exploring the packages used in
scverse that is needed for xenium analyses
'''
exp_figures = pack_exp_dir / "figures"

# Working directory: package exploration
exp_workdir = exp_figures / today_str

# ------- xen_qc directory -------
'''
Directory for all qc analyses. Includes initial filtering analysis for one
sample. Also includes bulk analysis of multiple samples
'''
qc_figures = qc_dir / "figures"

# Working directory: qc
qc_workdir = qc_figures / today_str

# ------- xen_qc directory -------
'''
Directory for all clustering and annotation analyses
'''
annot_figures = annot_dir / "figures"

# Working directory: qc
annot_workdir = annot_figures / today_str



