#!/bin/bash
#SBATCH --job-name=ob0gbsa_merge
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_ob0_mmgbsa/logs/45_merge_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_ob0_mmgbsa/logs/45_merge_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=4G
set -uo pipefail
# Join per-pose MM-GBSA shards back onto their pose lists (for pdb+smiles) into ONE ranker input.
#   export SCORED="glob1 glob2 ..." POSELISTS="list1 list2 ..." OUT=<name>.csv
#   sbatch --export=ALL slurm/cyp3a4_ob0_mmgbsa/45_merge_scored.sbatch
# Export in the SHELL, never as --export=ALL,SCORED=...: sbatch splits that value on the comma.
# SCORED/POSELISTS are intentionally space-separated word lists, so they are passed UNQUOTED below.
R=/vast.mnt/project/bme-20255165-cofolding-physics
ES=$R/miniforge3/envs/sigmadock
D=$R/data/cyp3a4_ob0_mmgbsa
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
cd "$R" || exit 1
# shellcheck disable=SC2086
"$ES/bin/python" analysis/cyp3a4_ob0_mmgbsa/merge_scored.py \
  --scored ${SCORED} --poselists ${POSELISTS} --out "$D/${OUT}"
RC=$?
# Exit on the python status, never on echo's. `echo "exit=$?"` as the last command makes every
# FATAL inside the merge report COMPLETED 0, so an --dependency=afterok step would run happily on
# a stale or absent merged CSV.
echo "exit=$RC"
exit $RC
