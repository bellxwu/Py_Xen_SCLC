#!/bin/bash
#SBATCH -t 6:00:00
#SBATCH --mem=100G
#SBATCH -J CellTypist_Chan2021_model
#SBATCH -p veryhimem
#SBATCH -c 10
#SBATCH -N 1
#SBATCH --output=/cluster/home/bellwu/Py_Xenium/job_logs/CellTypist_Chan2021_model%j.out
#SBATCH --error=/cluster/home/bellwu/Py_Xenium/job_logs/CellTypist_Chan2021_model%j.err
#SBATCH -D /cluster/home/bellwu/Py_Xenium

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xenium-env

python -u analyses/03_xen_annot/h4h_scripts/CellTypist_model_from_Chan2021.py