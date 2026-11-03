#!/bin/bash
#SBATCH --job-name=cofold_report
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=slurm/analysis/logs/%x_%j.out
#SBATCH --error=slurm/analysis/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem=8G
set -uo pipefail
# Run one post-hoc report from analysis/cofold/ through Slurm (never inline on the login node):
#   cd cofolder && sbatch --job-name=<report> slurm/analysis/report.sbatch <script.py> [args...]
# e.g. sbatch slurm/analysis/report.sbatch compare_pool_geometry.py --guided mol1_tfg_tcryo5 --control mol1_tcryo5
# Reports read the campaign tables under $ROOT (config/cluster.env) and print to this job's log.
COFOLDER=${SLURM_SUBMIT_DIR:?submit from the cofolder repo root}
[ -f "$COFOLDER/config/cluster.env" ] || { echo "FATAL: submit from the cofolder repo root"; exit 1; }
# shellcheck source=../../config/cluster.env
source "$COFOLDER/config/cluster.env"
SCRIPT=${1:?usage: report.sbatch <script.py> [args...]}; shift
[ -f "$COFOLDER/analysis/cofold/$SCRIPT" ] || { echo "FATAL: no analysis/cofold/$SCRIPT"; exit 1; }
export COFOLD_ROOT=$ROOT PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
"$ENV_QC/bin/python" "$COFOLDER/analysis/cofold/$SCRIPT" "$@"
rc=$?
echo "exit=$rc"
exit $rc
