#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH --mem=96G
#SBATCH -J qc_pp_leid
#SBATCH -p himem
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --output=/cluster/home/bellwu/Py_Xenium/job_logs/qc_pp_leid%j.out
#SBATCH --error=/cluster/home/bellwu/Py_Xenium/job_logs/qc_pp_leid%j.err
#SBATCH -D /cluster/home/bellwu/Py_Xenium

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xenium-env

python -u analyses/02_xen_qc/h4h_scripts/qc_pp_leid.py