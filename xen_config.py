"""
Description: Configuration script for all paths for Xenium python analysis.
Analysis pipeline edited for the location of the output files on the cluster. 

Created on Thu Nov  6 13:28:30 2025
@author: bellwu
"""

from pathlib import Path
from datetime import date

# --- Base Directory ---
base_dir = Path(__file__).resolve().parent

# -- Data Directories --
data_dir = Path("/cluster/home/lokgroup/Xenium_runs")
xen_1 = data_dir / "20241212__190528__Lok_Jalal_241210"
xen_2 = data_dir / "20250328__182124__Lok_Jalal_250324"
# samples
dir_23169 = xen_1 / "output-XETG00082__0051610__23169__20241212__190621"

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

# ------- package_exploration directory -------
'''
Directory here contains scripts and results of exploring the packages used in
scverse that is needed for xenium analyses
'''
exp_scripts = pack_exp_dir / "scripts" # scripts for subanalysis
exp_figures = pack_exp_dir / "figures"

# Working directory: package exploration
exp_workdir = exp_figures / today_str
exp_workdir.mkdir(parents=True, exist_ok=True)


# ------- xen_qc directory -------
'''
'''
qc_scripts = qc_dir / "scripts"
qc_figures = qc_dir / "figures"

# Working directory: qc
qc_workdir = qc_figures / today_str
qc_workdir.mkdir(parents=True, exist_ok=True)





