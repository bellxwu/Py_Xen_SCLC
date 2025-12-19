#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --mem=64G
#SBATCH -J single_pp_indiv_samples
#SBATCH -p himem
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --output=/cluster/home/bellwu/Py_Xenium/job_logs/single_pp_indiv_samples%j.out
#SBATCH --error=/cluster/home/bellwu/Py_Xenium/job_logs/single_pp_indiv_samples%j.err
#SBATCH -D /cluster/home/bellwu/Py_Xenium

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xenium-env

python -u analyses/02_xen_qc/h4h_scripts/single_pp_indiv_samples.py