#!/bin/bash
#SBATCH --job-name=cofold_mlff
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=slurm/3_rerank/logs/%x_%A_%a.out
#SBATCH --error=slurm/3_rerank/logs/%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=48G
set -uo pipefail
# STAGE 3 / scorer B -- MLFF interaction energy (Orb orbmol_v2) on the MM-GBSA top-K per ligand.
#   mlff.sbatch shortlist   top-MLFF_TOPK per ligand of $SCORED (same pool + Fe gate as select)
#                           -> E_int manifest off Mol-1's native CIFs (spin per pose from Fe-lig)
#   mlff.sbatch score       array task i of MLFF_SHARDS: prep (openmm83) then Orb (orb-env, CPU)
#   mlff.sbatch merge       Orb shards -> $EINT (mlff_score = E_int_conserving, force-gated)
# Submitted by `slurm/run.sh rerank` (after the MM-GBSA merge) or `slurm/run.sh mlff`.
#
# The poses are scored AS CO-FOLDED (no relaxation): this is not the shipped relaxtop3_eint arm,
# which relaxed the whole pool first. Orb self-noise is 0.131 eV on byte-identical input, so a
# pick it changes by less than that is unresolved -- merge_eint.py counts those.
source "${COFOLDER:?submit via slurm/run.sh}/slurm/lib.sh" && cofold_load || exit 1
P=$PIPE/3_rerank/mlff
PQ=$ENV_QC/bin/python
MAN=$MLFF_DIR/manifest.csv
MODE=${1:?usage: mlff.sbatch shortlist|score|merge}
mkdir -p "$MLFF_DIR"
echo "arm=$POOL mode=$MODE scored=$SCORED topk=$MLFF_TOPK pools=$RANK_POOLS pH=$PREP_PH heme=$HEM_PROTONATION"

case "$MODE" in
shortlist)
    POPT=(); [ "$RANK_POOLS" = all ] || POPT=(--pools "$RANK_POOLS")
    "$PQ" "$PIPE/3_rerank/select/build_submission.py" --scores "$SCORED" --score-col delta_total \
        --name "mlff_$POOL" --out-dir "$MLFF_DIR" --topk "$MLFF_TOPK" ${POPT[@]+"${POPT[@]}"} \
        --fe-in-max "$FE_IN" --fe-retain-max "$FE_RETAIN" --fe-ladder "$FE_LADDER" \
        --shortlist "$MLFF_DIR/shortlist.csv" || exit 1
    # the CIF (Mol-1's own mmCIF; ligand residue l01) and Fe-ligand distance come from the QC tables
    QCS=$("$PQ" - "$MLFF_DIR/shortlist.csv" "$GEN_DIR" <<'PY'
import sys
import pandas as pd
pools = sorted(pd.read_csv(sys.argv[1]).pool.unique())
print(" ".join(f"{sys.argv[2]}/poses_qc_{p}.csv" for p in pools))
PY
) || exit 1
    # shellcheck disable=SC2086
    "$PQ" "$P/eint_manifest.py" --qc $QCS --only "$MLFF_DIR/shortlist.csv" --all-ligands \
        --lig-res l01 --out "$MAN" || exit 1
    echo "manifest: $(nrows "$MAN") poses -> $MLFF_SHARDS shards (~70 CPU-s/pose for Orb)"
    ;;
score)
    I=$SLURM_ARRAY_TASK_ID
    [ -s "$MAN" ] || { echo "FATAL: no manifest $MAN -- the shortlist job did not finish"; exit 1; }
    rm -f "$MLFF_DIR/orb_s$I.csv" "$MLFF_DIR/prep_s$I.csv"    # a stale shard must not survive a rerun
    # PREP_SEED pins the PDBFixer hydrogen placement (else ~0.03 eV of E_int is prep noise);
    # CCD_COMPONENTS is read only if the cached HEM template is missing.
    PREP_SEED=20261030 CCD_COMPONENTS=$ROOT/common/components.cif PREP_PH=$PREP_PH \
    HEM_PROTONATION=$HEM_PROTONATION LD_LIBRARY_PATH=$ENV_AMBER/lib \
        "$ENV_AMBER/bin/python" "$P/eint_prep_shard.py" --manifest "$MAN" --out-dir "$MLFF_PREP" \
        --status-csv "$MLFF_DIR/prep_s$I.csv" --shard "$I" --nshards "$MLFF_SHARDS" || exit 1
    TORCHDYNAMO_DISABLE=1 LD_LIBRARY_PATH=$ENV_MLFF/lib \
        "$ENV_MLFF/bin/python" "$P/orb_eint_score.py" --manifest "$MAN" --prep-dir "$MLFF_PREP" \
        --out-csv "$MLFF_DIR/orb_s$I.csv" --model "$MLFF_MODEL" --shard "$I" --nshards "$MLFF_SHARDS"
    rc=$?; echo "exit=$rc"; exit $rc
    ;;
merge)
    NO=$(compgen -G "$MLFF_DIR/orb_s*.csv" | wc -l)
    echo "Orb shards present: $NO/$MLFF_SHARDS"
    "$PQ" "$P/merge_eint.py" --shortlist "$MLFF_DIR/shortlist.csv" --orb "$MLFF_DIR/orb_s*.csv" \
        --out "$EINT"
    rc=$?; echo "exit=$rc"; exit $rc
    ;;
*)
    echo "FATAL: mode must be shortlist|score|merge, got '$MODE'"; exit 2 ;;
esac
