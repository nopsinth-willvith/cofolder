#!/bin/bash
#SBATCH --job-name=cofold_gbsa
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=slurm/3_rerank/logs/%x_%A_%a.out
#SBATCH --error=slurm/3_rerank/logs/%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G
set -uo pipefail
# run.sh OVERRIDES --time from the actual poses-per-shard (14 CPU-min/pose + 30 min of build);
# the header value above is only the fallback for a hand-submitted shard.
# STAGE 3 / scorer A -- MM-GBSA, one array task = one shard of the pb_pass poselist.
# Submitted by `slurm/run.sh rerank <arm.env>`, which sizes the array and writes <SHARDS>.nshards.
#
# Per pose (~12 CPU-min, serial): build_system.py (ff14SB + Shahrokh IC6 heme + GAFF1/AM1-BCC,
# residues 25-503, 10/16 A pocket) then mmgbsa_modes.py (explicit-solvent minimisation, backbone
# restrained, then igb=2 MM-GBSA). --no-fe-bond for EVERY pose: the Fe-donor contact stays
# non-bonded, so ante-MMPBSA's :LIG split severs nothing and all poses share one physical model.
#
# The build is skipped when <WORK>/systems/<pose>/<pose>/system.prmtop exists. WORK is therefore
# tagged by PREP_PH: a pH-4.8 run against the pH-7.4 systems would skip every build and emit a
# table labelled 4.8 that is numerically 7.4. A pH variant also needs its own SET (new shards).
source "${COFOLDER:?submit via slurm/run.sh}/slurm/lib.sh" && cofold_load || exit 1
P3=$PIPE/3_rerank/mmgbsa
E=$ENV_AMBER
W=$SCRATCH/cyp3a4_ob0_mmgbsa                     # the campaign's existing pH-7.4 systems
[ "$PREP_PH" = 7.4 ] || W=${W}_ph$PREP_PH
export PREP_PH AMBERHOME=$E PATH=$E/bin:$PATH LD_LIBRARY_PATH=$E/lib:${LD_LIBRARY_PATH:-}
SH=$SLURM_ARRAY_TASK_ID
NSHARDS=$(cat "${SHARDS}.nshards") || { echo "FATAL: no ${SHARDS}.nshards -- submit via run.sh"; exit 1; }
echo "set=$SET shard=$SH/$NSHARDS pH=$PREP_PH work=$W poselist=$POSELIST"

# build this shard's systems: the same i % NSHARDS split mmgbsa_modes.py uses
awk -F, -v s="$SH" -v n="$NSHARDS" 'NR>1 && (NR-2)%n==s {print $1","$2","$3}' "$POSELIST" |
while IFS=, read -r POSE PDB SMI; do
    [ -f "$W/systems/$POSE/$POSE/system.prmtop" ] && continue
    mkdir -p "$W/systems/$POSE"
    "$E/bin/python" "$P3/build_system.py" \
        --pose "$PDB" --smiles "$SMI" --pose-id "$POSE" --out-dir "$W/systems/$POSE" \
        --first-residue 25 --no-fe-bond --pocket-free-radius 10 --pocket-buffer-radius 16 \
        > "$W/systems/$POSE/build.log" 2>&1 || echo "$POSE build FAILED (see $W/systems/$POSE/build.log)"
done

"$E/bin/python" "$P3/mmgbsa_modes.py" \
    --poselist "$POSELIST" --systems-root "$W/systems" --work-dir "$W/work" \
    --out-csv "${SHARDS}${SH}.csv" --shard "$SH" --nshards "$NSHARDS"
rc=$?
echo "exit=$rc"
exit $rc
