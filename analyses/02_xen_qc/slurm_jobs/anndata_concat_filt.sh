#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH --mem=32G
#SBATCH -J anndata_concat_filter
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --output=/cluster/home/bellwu/Py_Xenium/job_logs/anndata_concat_filter%j.out
#SBATCH --error=/cluster/home/bellwu/Py_Xenium/job_logs/anndata_concat_filter%j.err
#SBATCH -D /cluster/home/bellwu/Py_Xenium

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xenium-env

python -u analyses/02_xen_qc/scripts/xen_concat_filt.py
