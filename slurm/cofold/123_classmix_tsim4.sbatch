#!/bin/bash
#SBATCH --job-name=mol1_classmix
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/123_classmix_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/123_classmix_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=16G
set -uo pipefail
# Stage 123: a CLASS-SPLIT zip -- non-Type-II from mol1_tsim4, Type-II from mol1_tmplonly.
#
# Why this split and not the other way round, from this campaign's own measurements:
#   * non-Type-II: tsim4 beat tmplonly's cryo5 arm on PoseBusters validity, 7 of 9 ligands
#     favouring it, 0 against, 2 tied, p=0.016 -- the arm's one statistically supported gain.
#   * Type-II: tsim4 and the control were a dead wash on validity (4/3/4 tied, p=1.000), and
#     MM-GBSA preferred the unguided pools by a median +5.19 kcal/mol, so there is no reason to
#     move those eleven. tmplonly also draws 9 of its 11 Type-II from tcryo5/t38FX.
# So each class comes from the arm with evidence for it, and nothing is re-docked or re-ranked.
#
# ⚠ The members are BYTE COPIES of two zips that each passed their own from-the-bytes PoseBusters
# audit, and build_class_hybrid refuses a member whose source audit row is not pb_pass True. That
# is necessary but not sufficient: phases 2-4 below re-validate and RE-AUDIT the finished zip from
# its own bytes, and THAT audit is the claim.
R=/vast.mnt/project/bme-20255165-cofolding-physics
ES=$R/miniforge3/envs/sigmadock
EB=$R/miniforge3/envs/biokinema
EP=$R/miniforge3/envs/plinder
SUB=$R/data/cyp3a4_mol1/submission
NAME=${NAME:-submission_mol1_classmix_tsim4non}
T2_FROM=${T2_FROM:-submission_mol1_tmplonly_mmgbsa}
NON_FROM=${NON_FROM:-submission_mol1_tsim4_mmgbsa}
ZIP=$SUB/$NAME.zip
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
cd "$R" || exit 1
RC=0

echo "######## 0 gate"
for n in "$T2_FROM" "$NON_FROM"; do
    [ -s "$SUB/$n.zip" ] || { echo "FATAL: missing $SUB/$n.zip" >&2; exit 3; }
done
if [ -f "$ZIP" ] && [ "${ALLOW_OVERWRITE:-0}" != 1 ]; then
    echo "FATAL: $NAME.zip exists -- it and its audit would be destroyed. ALLOW_OVERWRITE=1." >&2
    exit 3
fi

echo "######## 1 assemble by class (byte copies; source audits gate each member)"
"$ES/bin/python" analysis/cyp3a4_md_eint/build_class_hybrid.py \
    --typeii-from "$T2_FROM" --nontypeii-from "$NON_FROM" --name "$NAME" \
    --sub-dir "$SUB" --audit-dir "$SUB" || exit 1
[ -s "$ZIP" ] || { echo "FATAL: $ZIP was not written" >&2; exit 1; }

echo "######## 2 official validator"
"$EB/bin/python" analysis/cyp3a4_gamd/official_validate.py --zip "$ZIP" --label "$NAME" || RC=1

echo "######## 3 audit from the zip alone"
A=$SUB/.${NAME}_audit; rm -rf "$A"; mkdir -p "$A/members"
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
  --out-csv "$SUB/${NAME}.zipaudit.csv" > "$A/audit.log" 2>&1
"$EB/bin/python" -I - "$SUB/${NAME}.zipaudit.csv" <<'PY' || RC=1
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="")))
bad = [r["base_id"] for r in rows if str(r.get("pb_pass", "")).lower() != "true"]
print(f"zip audit: pb_pass {len(rows) - len(bad)}/{len(rows)}" + (f"  FAILING {bad}" if bad else "  ALL CLEAR"))
sys.exit(1 if bad or len(rows) != 20 else 0)
PY

echo "######## 4 change vs each parent and the four-pool build"
for REF in "$SUB/$T2_FROM.zip" "$SUB/$NON_FROM.zip" "$SUB/submission_mol1_tmplsim_mmgbsa.zip"; do
    [ -s "$REF" ] || continue
    RS=$(basename "$REF" .zip); RS=${RS#submission_}
    echo "---- vs $RS"
    LD_LIBRARY_PATH=$EP/lib "$EP/bin/python" analysis/cyp3a4_relax/lddt_pli_vs_submission.py \
      --model-zip "$ZIP" --ref-zip "$REF" --out-csv "$SUB/${NAME}_vs_${RS}.csv" --portal \
      > "$A/portal_${RS}.log" 2>&1 || RC=1
    "$ES/bin/python" - "$SUB/${NAME}_vs_${RS}.csv" <<'PY'
import csv, sys
from pathlib import Path
R = Path("/vast.mnt/project/bme-20255165-cofolding-physics")
t2 = {l.strip() for l in open(R / "data/cyp3a4_xmodel/typeii_base_ids.txt") if l.strip()}
rows = list(csv.DictReader(open(sys.argv[1], newline="")))
for grp in ("TypeII", "nonTypeII"):
    sel = [r for r in rows if (r["base_id"] in t2) == (grp == "TypeII")]
    mv = [r for r in sel if r.get("portal_bisyrmsd") not in ("", None)
          and float(r["portal_bisyrmsd"]) >= 2.0]
    print(f"  {grp}: {len(mv)}/{len(sel)} members move >= 2 A")
PY
done
echo "exit=$RC"; exit $RC
