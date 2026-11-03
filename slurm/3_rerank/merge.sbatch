#!/bin/bash
#SBATCH --job-name=cofold_merge
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=slurm/3_rerank/logs/%x_%j.out
#SBATCH --error=slurm/3_rerank/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=4G
set -uo pipefail
# STAGE 3 / scorer A, merge -- MM-GBSA shards joined back onto the poselist (pdb + smiles) into
# scored_<SET>.csv, with every pool attributed to its checkpoint (CKPT_LABEL). Submitted by
# `slurm/run.sh rerank` with --dependency=afterany on the array: one OOMed shard must not strand
# the merge, so this job gates on CONTENT -- it names the missing shards and exits on the merge's
# own status, which the MLFF chain depends on (afterok).
source "${COFOLDER:?submit via slurm/run.sh}/slurm/lib.sh" && cofold_load || exit 1
NSH=$(cat "${SHARDS}.nshards" 2>/dev/null || echo "?")
MISSING=""
if [ "$NSH" != "?" ]; then
    for i in $(seq 0 $((NSH - 1))); do [ -s "${SHARDS}$i.csv" ] || MISSING="$MISSING $i"; done
fi
echo "set=$SET shards present: $(compgen -G "${SHARDS}[0-9]*.csv" | wc -l)/$NSH${MISSING:+  MISSING:$MISSING}"
[ -z "$MISSING" ] || echo "⚠ re-score the missing shards with ALLOW_RESUME=1 slurm/run.sh rerank $CFG"

"$ENV_QC/bin/python" "$PIPE/3_rerank/mmgbsa/merge_scored.py" \
    --scored "${SHARDS}[0-9]*.csv" --poselists "$POSELIST" \
    --pool-ckpt "$POOL=$CKPT_LABEL" --out "$SCORED"
rc=$?
# Exit on the python status, never on echo's: an `echo` as the last command would make every
# FATAL inside the merge report COMPLETED 0 and the afterok MLFF chain would run on a stale table.
echo "exit=$rc"
exit $rc
