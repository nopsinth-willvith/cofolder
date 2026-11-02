#!/bin/bash
#SBATCH --job-name=mol1_30_qc
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/30_qc_%x_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/30_qc_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=04:00:00
#SBATCH --mem=32G
set -uo pipefail
# Stage 30: Mol-1 CIFs -> submission-shaped PDBs + confidence table, then PoseBusters.
#   export POOL=mol1_full; sbatch --export=ALL --job-name=$POOL 30_collect_qc.sbatch
# Two environments on purpose: collect_poses needs the Mol-1 venv (gemmi), while validate_frames is
# the campaign's existing PoseBusters harness and runs in sigmadock, exactly as it does for every
# other pool -- so this arm's pb_pass means the same thing as the others'.
R=/vast.mnt/project/bme-20255165-cofolding-physics
REPO=/project/bme-20255165-foldx/Mol-1
ES=$R/miniforge3/envs/sigmadock
D=$R/data/cyp3a4_mol1
P=${POOL:?set POOL}
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
cd "$R" || exit 1
RC=0

echo "######## 1 collect + convert"
"$REPO/.venv/bin/python" analysis/cyp3a4_mol1/collect_poses.py \
  --dump-dir "$D/gen/$P" --pool "$P" --pdb-root "$D/poses" \
  --out-csv "$D/poses_$P.csv" || RC=1

echo "######## 2 PoseBusters on the converted PDBs"
# validate_frames wants base_id,pdb,run; `run` carries pose_id so the result joins back 1:1.
"$ES/bin/python" - "$D/poses_$P.csv" "$D/frames_$P.csv" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], newline="")))
with open(sys.argv[2], "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["base_id", "pdb", "run", "smiles"]); w.writeheader()
    for r in rows:
        w.writerow({"base_id": r["base_id"], "pdb": r["pdb"], "run": r["pose_id"],
                    "smiles": r["smiles"]})
print(f"{len(rows)} frames -> {sys.argv[2]}")
PY
"$ES/bin/python" analysis/cyp3a4_gamd/validate_frames.py \
  --frames-csv "$D/frames_$P.csv" --out-csv "$D/pb_$P.csv" || RC=1

echo "######## 3 join + summarise"
"$ES/bin/python" - "$D/poses_$P.csv" "$D/pb_$P.csv" "$D/poses_qc_$P.csv" <<'PY'
import sys
import pandas as pd
poses = pd.read_csv(sys.argv[1])
pb = pd.read_csv(sys.argv[2]).rename(columns={"run": "pose_id"})
keep = [c for c in ["pose_id", "pb_pass", "pb_failed", "max_cc_bond_A", "fe_lig_min_A",
                    "fe_lig_atom", "fe_lig_element"] if c in pb.columns]
out = poses.merge(pb[keep], on="pose_id", how="left", validate="1:1")
out["pb_pass"] = out.pb_pass.astype(str).str.lower().eq("true")
out.to_csv(sys.argv[3], index=False)
n, ok = len(out), int(out.pb_pass.sum())
print(f"{sys.argv[3]}: {n} poses, pb_pass {ok} ({100*ok/n:.1f}%), {out.base_id.nunique()} ligands")
thin = out[out.pb_pass].groupby("base_id").size()
thin = thin[thin < 5]
if len(thin):
    print("  thin after PB:", thin.to_dict())
# The Fe band is what the Type-II gate will use on BOTH ranking criteria, so report it now --
# an arm with no coordinating poses cannot support a Type-II pick however good its confidence is.
t2 = set(open("/vast.mnt/project/bme-20255165-cofolding-physics/data/cyp3a4_xmodel/"
              "typeii_base_ids.txt").read().split())
g = out[out.pb_pass & out.base_id.isin(t2)]
print(f"  Type-II pb_pass poses {len(g)}; with Fe-lig <= 2.6 A: {int((g.fe_lig_min_A <= 2.6).sum())}"
      f" over {g[g.fe_lig_min_A <= 2.6].base_id.nunique()}/{len(t2)} ligands")
PY
echo "exit=$RC"
exit $RC
