#!/bin/bash
#SBATCH --job-name=cofold_select
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=slurm/3_rerank/logs/%x_%j.out
#SBATCH --error=slurm/3_rerank/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G
set -uo pipefail
# STAGE 3 / select -- one member per ligand -> zip -> official validator -> audit from the zip's own
# bytes -> portal lDDT-PLI change vs each reference zip. Submitted by `slurm/run.sh select`.
#
# SELECT_MODE
#   rank      rank RANK_SCORES by SCORER (mmgbsa: delta_total | mlff: mlff_score), lower = better,
#             Type-II Fe gate (FE_IN / FE_RETAIN / FE_LADDER, relaxed per ligand only where forced),
#             PoseBusters re-run on the exact file down the top PB_TOPK
#   classmix  Type-II members from zip TYPEII_FROM, non-Type-II from NONTYPEII_FROM (byte copies)
#   override  every member from BASE_ZIP except the ligands FROM_ALT, taken from ALT_ZIP
#   audit     build nothing; re-run phases 2-4 on the existing $SUB/$NAME.zip
# The audit (phase 3) is the claim, not the builder's own log: it re-reads the shipped bytes.
source "${COFOLDER:?submit via slurm/run.sh}/slurm/lib.sh" && cofold_load || exit 1
S=$PIPE/3_rerank/select
PQ=$ENV_QC/bin/python
ZIP=$SUB/$NAME.zip
W=$SUB/.${NAME}_work
mkdir -p "$W"
RC=0

echo "######## 1 build ($SELECT_MODE) -> $ZIP"
case "$SELECT_MODE" in
rank)
    # shellcheck disable=SC2206
    SCORES=($RANK_SCORES)
    POOLS=$RANK_POOLS
    if [ "$ADD_REPAIRED" = 1 ]; then
        # the repaired OCNT-2313038 rows belong to pool mol1_repair -- a pool restriction must name it
        "$PQ" "$S/add_repaired.py" --table "${SCORES[@]}" --scores "$GEN_DIR/mmgbsa_mol1fix.csv" \
            --out "$W/rank_scores.csv" || { echo "BUILD FAILED (add_repaired)"; exit 1; }
        SCORES=("$W/rank_scores.csv")
        [ "$POOLS" = all ] || POOLS="$POOLS mol1_repair"
    fi
    POPT=(); [ "$POOLS" = all ] || POPT=(--pools "$POOLS")
    echo "scores=${SCORES[*]} col=$RANK_COL pools=$POOLS Fe in<=$FE_IN out<=$FE_RETAIN ladder=$FE_LADDER"
    "$PQ" "$S/build_submission.py" --scores "${SCORES[@]}" --score-col "$RANK_COL" --name "$NAME" \
        --out-dir "$SUB" --topk "$PB_TOPK" ${POPT[@]+"${POPT[@]}"} \
        --fe-in-max "$FE_IN" --fe-retain-max "$FE_RETAIN" --fe-ladder "$FE_LADDER" \
        || { echo "BUILD FAILED"; exit 1; }
    ;;
classmix)
    : "${TYPEII_FROM:?classmix needs TYPEII_FROM}" "${NONTYPEII_FROM:?classmix needs NONTYPEII_FROM}"
    "$PQ" "$S/build_class_hybrid.py" --typeii-from "$TYPEII_FROM" --nontypeii-from "$NONTYPEII_FROM" \
        --name "$NAME" --sub-dir "$SUB" --audit-dir "$SUB" || { echo "BUILD FAILED"; exit 1; }
    ;;
override)
    : "${BASE_ZIP:?override needs BASE_ZIP}" "${ALT_ZIP:?override needs ALT_ZIP}" "${FROM_ALT:?override needs FROM_ALT}"
    # shellcheck disable=SC2086
    "$PQ" "$S/build_zip_override.py" --base "$BASE_ZIP" --alt "$ALT_ZIP" --from-alt $FROM_ALT \
        --sub-dir "$SUB" --out "$W/override_scores.csv" || { echo "BUILD FAILED (override)"; exit 1; }
    # both parents are already gated and 20/20 pb_pass: a constant score, depth 1, no Fe gate --
    # a member that now fails PoseBusters is a real bug and must surface, not be swapped out
    "$PQ" "$S/build_submission.py" --scores "$W/override_scores.csv" --score-col mix --name "$NAME" \
        --out-dir "$SUB" --topk 1 || { echo "BUILD FAILED"; exit 1; }
    ;;
audit)
    [ -f "$ZIP" ] || { echo "FATAL: audit mode but no $ZIP"; exit 1; } ;;
*)
    echo "FATAL: SELECT_MODE must be rank|classmix|override|audit, got '$SELECT_MODE'"; exit 2 ;;
esac

echo "######## 2 official challenge validator"
"$ENV_VALID/bin/python" "$S/official_validate.py" --zip "$ZIP" --label "$NAME" || RC=1

echo "######## 3 audit from the zip alone (PoseBusters on the shipped bytes)"
A=$SUB/.${NAME}_audit; rm -rf "$A"; mkdir -p "$A/members"
"$PQ" -I - "$ZIP" "$A/members" <<'PY' || RC=1
import pathlib, sys, zipfile
zp, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
with zipfile.ZipFile(zp) as zf:
    names = [n for n in zf.namelist() if n.endswith(".pdb")]
    for n in names:
        (out / pathlib.Path(n).name).write_bytes(zf.read(n))
print(f"unpacked {len(names)} members")
PY
echo "base_id,pdb,run" > "$A/members.csv"
for f in "$A"/members/*.pdb; do echo "$(basename "${f%.pdb}"),$f,$NAME" >> "$A/members.csv"; done
"$PQ" "$PIPE/2_cofold/validate_frames.py" --frames-csv "$A/members.csv" \
    --out-csv "$SUB/${NAME}.zipaudit.csv" > "$A/audit.log" 2>&1 || RC=1
"$PQ" -I - "$SUB/${NAME}.zipaudit.csv" <<'PY' || RC=1
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="")))
bad = [r["base_id"] for r in rows if str(r.get("pb_pass", "")).lower() != "true"]
print(f"zip audit: pb_pass {len(rows) - len(bad)}/{len(rows)}" + (f"  FAILING {bad}" if bad else "  ALL CLEAR"))
sys.exit(1 if bad or len(rows) != 20 else 0)
PY

echo "######## 4 change vs each reference submission (portal instrument; AGREEMENT, not accuracy)"
for REFZIP in $REFS; do
    RS=$(basename "$REFZIP" .zip); RS=${RS#submission_}
    [ -f "$REFZIP" ] || { echo "---- vs $RS: reference zip missing, skipped"; continue; }
    [ "$(readlink -f "$REFZIP")" = "$(readlink -f "$ZIP")" ] && { echo "---- vs $RS: same zip, skipped"; continue; }
    echo "---- vs $RS"
    OUTC=$SUB/${NAME}_vs_${RS}.csv
    LD_LIBRARY_PATH=$ENV_OST/lib "$ENV_OST/bin/python" "$S/lddt_pli_vs_submission.py" \
        --model-zip "$ZIP" --ref-zip "$REFZIP" --out-csv "$OUTC" --portal > "$A/portal_${RS}.log" 2>&1 \
        || { echo "  portal diff FAILED (see $A/portal_${RS}.log)"; RC=1; continue; }
    "$PQ" - "$OUTC" "$TYPEII_LIST" <<'PY'
import sys
import pandas as pd
d = pd.read_csv(sys.argv[1])
t2 = set(open(sys.argv[2]).read().split())
rc = next(c for c in d.columns if "portal" in c and "bisyrmsd" in c)
d["mode"] = ["TypeII" if b in t2 else "nonTypeII" for b in d.base_id]
cols = [c for c in ["base_id", "mode", rc] + [c for c in d.columns if "com_shift" in c] if c in d]
print(d[cols].sort_values(["mode", "base_id"]).to_string(index=False))
for m, g in d.groupby("mode"):
    print(f"{m}: {int((g[rc] >= 2).sum())}/{len(g)} members move >= 2 A ({rc})")
PY
done
echo "zip: $ZIP"
echo "provenance: $SUB/$NAME.provenance.csv   audit: $SUB/${NAME}.zipaudit.csv"
echo "exit=$RC"
exit $RC
