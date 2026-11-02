#!/bin/bash
#SBATCH --job-name=cyp3a4rx_45_lddtpli
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/vast.mnt/project/bme-20255165-cofolding-physics/slurm/cyp3a4_relax/logs/45_lddtpli_%j.out
#SBATCH --error=/vast.mnt/project/bme-20255165-cofolding-physics/slurm/cyp3a4_relax/logs/45_lddtpli_%j.err
#SBATCH --nodes=1 --ntasks-per-node=1 --cpus-per-task=8 --time=04:00:00 --mem=32G
set -uo pipefail
R=/vast.mnt/project/bme-20255165-cofolding-physics
# plinder is the only env with a REAL openstructure (2.8.0 + compound lib).
# openmm83/lib/python3.12/site-packages/ost holds an unrelated 2-file package of
# the same name with no __init__.py, which shadows the import as a namespace
# package and fails with "cannot import name 'conop' from 'ost' (unknown location)".
E=$R/miniforge3/envs/plinder
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
export LD_LIBRARY_PATH=$E/lib      # same GLIBCXX trap as every other openmm83 stage
cd "$R" || exit 1
"$E/bin/python" "$R/analysis/cyp3a4_relax/lddt_pli_vs_submission.py" \
    --model-zip "${MODEL:?set MODEL}" --ref-zip "${REF:?set REF}" \
    --out-csv "${OUT:?set OUT}" $([ "${PORTAL:-0}" = 1 ] && echo --portal)
rc=$?
echo "exit=$rc"
exit $rc
