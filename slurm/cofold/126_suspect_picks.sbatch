#!/bin/bash
#SBATCH --job-name=mol1_suspect
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/126_suspect_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/126_suspect_%j.err
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=2 --time=00:15:00 --mem=8G
set -uo pipefail
R=/vast.mnt/project/bme-20255165-cofolding-physics
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
cd "$R" || exit 1
"$R/miniforge3/envs/sigmadock/bin/python" analysis/cyp3a4_mol1/suspect_tmpl_picks.py
