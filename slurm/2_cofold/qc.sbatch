#!/bin/bash
#SBATCH --job-name=cofold_qc
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=slurm/2_cofold/logs/qc_%x_%j.out
#SBATCH --error=slurm/2_cofold/logs/qc_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=04:00:00
#SBATCH --mem=32G
set -uo pipefail
# STAGE 2b -- QC: Mol-1 CIFs -> submission-shaped PDBs + confidence table -> PoseBusters ->
# poses_qc_<POOL>.csv, then the pb_pass poselist both re-rankers consume. ~25 min / 1000 poses.
# Submitted by `slurm/run.sh cofold` with --dependency=afterany on the generation array.
#
# Two environments on purpose: collect_poses needs the Mol-1 venv (gemmi); PoseBusters runs in the
# campaign's QC env exactly as for every other pool, so pb_pass means the same thing everywhere.
# Never change the harness between pools. pb_pass is PHYSICAL PLAUSIBILITY, not accuracy.
source "${COFOLDER:?submit via slurm/run.sh}/slurm/lib.sh" && cofold_load || exit 1
PQ=$ENV_QC/bin/python
P2=$PIPE/2_cofold
RC=0

echo "######## 1 collect + convert ($GEN)"
"$MOL1_PY" "$P2/collect_poses.py" --dump-dir "$GEN" --pool "$POOL" \
    --pdb-root "$GEN_DIR/poses" --out-csv "$GEN_DIR/poses_$POOL.csv" || RC=1
[ -s "$GEN_DIR/poses_$POOL.csv" ] || { echo "FATAL: collect wrote no pose table"; exit 1; }

echo "######## 2 PoseBusters on the converted PDBs"
# validate_frames wants base_id,pdb,run; `run` carries pose_id so the result joins back 1:1.
"$PQ" - "$GEN_DIR/poses_$POOL.csv" "$GEN_DIR/frames_$POOL.csv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="")))
with open(sys.argv[2], "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["base_id", "pdb", "run", "smiles"]); w.writeheader()
    for r in rows:
        w.writerow({"base_id": r["base_id"], "pdb": r["pdb"], "run": r["pose_id"], "smiles": r["smiles"]})
print(f"{len(rows)} frames -> {sys.argv[2]}")
PY
"$PQ" "$P2/validate_frames.py" --frames-csv "$GEN_DIR/frames_$POOL.csv" \
    --out-csv "$GEN_DIR/pb_$POOL.csv" || RC=1

echo "######## 3 join + per-target summary"
# Per-target rates, never only the pooled percentage; a target with ZERO passing poses is named
# here because it hard-fails the submission builder three stages (and ~200 CPU-h) later.
"$PQ" - "$GEN_DIR/poses_$POOL.csv" "$GEN_DIR/pb_$POOL.csv" "$QC" "$TYPEII_LIST" "$EXPECT_POSES" <<'PY' || RC=1
import sys
import pandas as pd
poses = pd.read_csv(sys.argv[1])
pb = pd.read_csv(sys.argv[2]).rename(columns={"run": "pose_id"})
keep = [c for c in ["pose_id", "pb_pass", "pb_failed", "max_cc_bond_A", "fe_lig_min_A",
                    "fe_lig_atom", "fe_lig_element"] if c in pb.columns]
out = poses.merge(pb[keep], on="pose_id", how="left", validate="1:1")
out["pb_pass"] = out.pb_pass.astype(str).str.lower().eq("true")
out.to_csv(sys.argv[3], index=False)
t2 = set(open(sys.argv[4]).read().split())
expect = int(sys.argv[5])
n, ok = len(out), int(out.pb_pass.sum())
if not n:
    # every pose skipped (a killed array, a wrong dump dir): say so instead of dying in the
    # percentage and taking the per-target report -- the whole point of this job -- with it
    sys.exit("FATAL: the pose table has no rows -- nothing was collected from the dump dir")
print(f"{sys.argv[3]}: {n} poses, pb_pass {ok} ({100*ok/n:.1f}%), {out.base_id.nunique()} ligands")
per = out.groupby("base_id").agg(n=("pose_id", "size"), ok=("pb_pass", "sum"))
for b, r in per.iterrows():
    flag = "   <-- ZERO PASSING" if r.ok == 0 else ("   <-- short" if r.n < expect else "")
    print(f"  {b} {'T2' if b in t2 else '  '} {int(r.ok):3d}/{int(r.n)}{flag}")
starved = sorted(per.index[per.ok == 0])
print(f"STARVED: {','.join(starved) or '-'}")
cc = pd.to_numeric(out.max_cc_bond_A, errors="coerce")
print(f"max_cc_bond_A median {cc.median():.3f} max {cc.max():.3f}  >1.60 A: {int((cc > 1.60).sum())}")
# The Fe band is what the Type-II gate reads; carry the closest atom's ELEMENT too -- a contact
# through carbon rather than a nitrogen donor is a defect no distance reveals.
g = out[out.pb_pass & out.base_id.isin(t2)]
near = g[g.fe_lig_min_A <= 2.6]
print(f"Type-II pb_pass poses {len(g)}; Fe-lig <= 2.6 A: {len(near)} over "
      f"{near.base_id.nunique()}/{len(t2)} ligands; elements {near.get('fe_lig_element', pd.Series()).value_counts().to_dict()}")
PY

echo "######## 4 poselist for the re-rankers (every pb_pass pose; no top-N, no score cut)"
# SET_POOLS, not just this pool: several arms may share one shard family, and a poselist rewritten
# from the arm that finished last would silently drop its siblings' poses. Every declared pool must
# have a QC table first, so the poselist is written once, complete, by the last arm to be QC'd.
QCS=(); MISS=()
for P in $SET_POOLS; do
    if [ -s "$GEN_DIR/poses_qc_$P.csv" ]; then QCS+=("$GEN_DIR/poses_qc_$P.csv"); else MISS+=("$P"); fi
done
# A poselist that already has MM-GBSA shards is a scored SET: rewriting it would desynchronise the
# shards from their pose list, so it is left alone and named.
if compgen -G "${SHARDS}[0-9]*.csv" > /dev/null; then
    echo "  ⚠ $SET already has MM-GBSA shards -- NOT rewriting $POSELIST"
elif [ ${#MISS[@]} -gt 0 ]; then
    echo "  SET=$SET is shared with ${MISS[*]}, still un-QC'd -- NOT writing $POSELIST yet."
    echo "  QC those pools; the last one writes the complete poselist."
else
    "$PQ" "$PIPE/3_rerank/mmgbsa/poselist_from_qc.py" --qc "${QCS[@]}" --out "$POSELIST" || RC=1
    [ ${#QCS[@]} -gt 1 ] && echo "  (from ${#QCS[@]} QC tables: $SET_POOLS)"
fi
echo "poselist: $POSELIST ($(nrows "$POSELIST") poses)"
# Exit status is informational: collect returns 1 for a skipped pose as well as for a real
# failure. Downstream gates on the table's content, never on this code.
echo "exit=$RC"
exit $RC
