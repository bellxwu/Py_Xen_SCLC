#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --mem=64G
#SBATCH -J SCLC_xen_pp
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --output=/cluster/home/bellwu/Py_Xenium/job_logs/SCLC_xen_pp%j.out
#SBATCH --error=/cluster/home/bellwu/Py_Xenium/job_logs/SCLC_xen_pp%j.err
#SBATCH -D /cluster/home/bellwu/Py_Xenium

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xenium-env

python -u analyses/02_xen_qc/h4h_scripts/qc_pp_xen.py
