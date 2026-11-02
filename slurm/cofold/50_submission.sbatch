#!/bin/bash
#SBATCH --job-name=ob0gbsa_submit
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_ob0_mmgbsa/logs/submit_%x_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_ob0_mmgbsa/logs/submit_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=16G
set -uo pipefail
# export NAME=<zip name> SCORES=<per-pose csv> COLS=<col[,col]>; sbatch --export=ALL 50_...sbatch
# Optional: SEEDS / POOLS restrict the candidate pool, FE_IN / FE_RETAIN apply the Type-II Fe gate
# (fe_lig_in_A / fe_lig_out_A, see build_submission.fe_gate), MAXIMIZE=1 for a score where HIGHER
# is better (a confidence, not an energy), REFS is a space-separated list of zips to compare
# against.  EXPORT THEM IN THE SHELL, never as --export=ALL,SEEDS=...: sbatch splits that
# value on the comma and only the first item arrives (the trap that silently halved a pool once).
# Every one of these changes WHICH POSES COMPETE, so a name that claims a restriction the script was
# not given is a mislabelled submission -- the printed "candidate pool"/"Fe gate" lines are the check.
# 1 select + PoseBusters on the exact files (fallback down the ranking) + zip
# 2 official challenge validator on the zip
# 3 independent audit: unpack the zip and re-run PoseBusters on its members (no intermediate CSV)
# 4 portal-instrument change vs each reference submission: how many members really move
R=/vast.mnt/project/bme-20255165-cofolding-physics
ES=$R/miniforge3/envs/sigmadock
EB=$R/miniforge3/envs/biokinema
EP=$R/miniforge3/envs/plinder
OUT=${OUTDIR:-$R/data/cyp3a4_ob0_mmgbsa/submission}
ZIP=$OUT/$NAME.zip
REFS=${REFS:-$R/data/cyp3a4_relax/submission/submission_xrpli_basefree_fe460.zip}
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
cd "$R" || exit 1
RC=0

echo "######## 1 select + PoseBusters + zip"
"$ES/bin/python" analysis/cyp3a4_ob0_mmgbsa/build_submission.py \
  --scores "$SCORES" --score-col "$COLS" --name "$NAME" --out-dir "$OUT" \
  ${SEEDS:+--seeds "$SEEDS"} ${POOLS:+--pools "$POOLS"} \
  ${FE_IN:+--fe-in-max "$FE_IN"} ${FE_RETAIN:+--fe-retain-max "$FE_RETAIN"} \
  ${MAXIMIZE:+--maximize} ${FE_LADDER:+--fe-ladder "$FE_LADDER"} \
  || { echo "BUILD FAILED"; exit 1; }

echo "######## 2 official validator"
"$EB/bin/python" analysis/cyp3a4_gamd/official_validate.py --zip "$ZIP" --label "$NAME" || RC=1

echo "######## 3 audit from the zip alone"
A=$OUT/.${NAME}_audit; rm -rf "$A"; mkdir -p "$A/members"
"$EB/bin/python" -I -c "
import sys, zipfile, pathlib
zp, out = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
with zipfile.ZipFile(zp) as zf:
    names = [n for n in zf.namelist() if n.endswith('.pdb')]
    for n in names: (out / pathlib.Path(n).name).write_bytes(zf.read(n))
print(f'unpacked {len(names)} members')
" "$ZIP" "$A/members"
echo "base_id,pdb,run" > "$A/members.csv"
for f in "$A"/members/*.pdb; do echo "$(basename "${f%.pdb}"),$f,$NAME" >> "$A/members.csv"; done
"$ES/bin/python" analysis/cyp3a4_gamd/validate_frames.py --frames-csv "$A/members.csv" \
  --out-csv "$OUT/${NAME}.zipaudit.csv" > "$A/audit.log" 2>&1
"$EB/bin/python" -I - "$OUT/${NAME}.zipaudit.csv" <<'PY' || RC=1
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="")))
bad = [r["base_id"] for r in rows if str(r.get("pb_pass", "")).lower() != "true"]
print(f"zip audit: pb_pass {len(rows) - len(bad)}/{len(rows)}" + (f"  FAILING {bad}" if bad else "  ALL CLEAR"))
sys.exit(1 if bad or len(rows) != 20 else 0)
PY

echo "######## 4 change vs each reference submission (portal instrument)"
for REF in $REFS; do
RS=$(basename "$REF" .zip); RS=${RS#submission_}
echo "---- vs $RS"
LD_LIBRARY_PATH=$EP/lib "$EP/bin/python" analysis/cyp3a4_relax/lddt_pli_vs_submission.py \
  --model-zip "$ZIP" --ref-zip "$REF" --out-csv "$OUT/${NAME}_vs_${RS}.csv" --portal > "$A/portal_${RS}.log" 2>&1 || RC=1
"$ES/bin/python" - "$OUT/${NAME}_vs_${RS}.csv" <<'PY'
import sys
import pandas as pd
t2 = set(open("/vast.mnt/project/bme-20255165-cofolding-physics/data/cyp3a4_xmodel/typeii_base_ids.txt").read().split())
d = pd.read_csv(sys.argv[1])
rc = next(c for c in d.columns if "portal" in c and "bisyrmsd" in c)
d["mode"] = ["TypeII" if b in t2 else "nonTypeII" for b in d.base_id]
cols = [c for c in ["base_id", "mode", rc] + [c for c in d.columns if "com_shift" in c] if c in d]
print(d[cols].sort_values(["mode", "base_id"]).to_string(index=False))
for m, g in d.groupby("mode"):
    print(f"{m}: {int((g[rc] >= 2).sum())}/{len(g)} members move >= 2 A ({rc})")
PY
done
echo "exit=$RC"
exit $RC
