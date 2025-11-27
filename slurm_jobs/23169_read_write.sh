#!/bin/bash
#SBATCH -t 1:00
#SBATCH --mem=16G
#SBATCH -J 23169_read_write
#SBATCH -p all
#SBATCH -c 4
#SBATCH -N 1
#SBATCH --output=/cluster/home/bellwu/Py_Xenium/job_logs/23169_read_write%j.out
#SBATCH --error=/cluster/home/bellwu/Py_Xenium/job_logs/23169_read_write%j.err

source ~/miniconda3/etc/profile.d/conda.sh
conda activate xen_reader

python -u /cluster/home/bellwu/Py_Xenium/analyses/02_xen_qc/scripts/xenium_zarr.py
