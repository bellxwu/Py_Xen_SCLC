#!/bin/bash
#SBATCH -t 10
#SBATCH --mem=4G
#SBATCH -J CellTypist_Chan2021_model
#SBATCH -p himem
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -D /cluster/home/bellwu/Py_Xenium
#SBATCH --output=/cluster/home/bellwu/Py_Xenium/job_logs/test_path%j.out
#SBATCH --error=/cluster/home/bellwu/Py_Xenium/job_logs/test_path%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xenium-env

python -u /cluster/home/bellwu/Py_Xenium/analyses/03_xen_annot/h4h_scripts/test.py