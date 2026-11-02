#!/bin/bash
#SBATCH --job-name=mol1_tfg_cmp
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/111_tfg_cmp_%x_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/111_tfg_cmp_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --mem=8G
set -uo pipefail
# Stage 111: does TFG guidance actually improve the ligand geometry it steers?
#   sbatch --dependency=afterok:<qc jobid> 111_tfg_compare.sbatch
# Pairs per LIGAND (n=20), not per pose: mol/tfg/engine.py:493 draws an extra randn per step, so
# pose (lig,seed,sample) in the guided arm is a different trajectory from the same label in the
# control and a pose-level paired test would be measuring an arbitrary alignment.
R=/vast.mnt/project/bme-20255165-cofolding-physics
ES=$R/miniforge3/envs/sigmadock
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
cd "$R" || exit 1
"$ES/bin/python" analysis/cyp3a4_mol1/compare_pool_geometry.py \
  --guided "${GUIDED:-mol1_tfg_tcryo5}" --control "${CONTROL:-mol1_tcryo5}"
