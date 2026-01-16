#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --mem=150G
#SBATCH -J CellTypist_SCLC_LabelTransfer
#SBATCH -p himem
#SBATCH -c 10
#SBATCH -N 1
#SBATCH -D /cluster/home/bellwu/Py_Xenium
#SBATCH --output=/cluster/home/bellwu/Py_Xenium/job_logs/CellTypist_SCLC_LabelTransfer%j.out
#SBATCH --error=/cluster/home/bellwu/Py_Xenium/job_logs/CellTypist_SCLC_LabelTransfer%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xenium-env

python -u /cluster/home/bellwu/Py_Xenium/analyses/03_xen_annot/h4h_scripts/CellTypist_labeltransfer_SCLC.py