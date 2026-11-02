"""Turn the Mol-1 QC tables into the MM-GBSA stage's poselist, PoseBusters-passing poses only.

slurm/cyp3a4_ob0_mmgbsa/10_build_score.sbatch reads
`data/cyp3a4_ob0_mmgbsa/poselist_<SET>.csv` and takes columns 1,2,3 positionally as
pose_id, pdb, smiles, so the column ORDER here is load-bearing, not cosmetic.

This writes into the MM-GBSA stage's directory rather than the Mol-1 one on purpose: both arms of
this campaign score from one place, and a pool that competes in a ranking should be listed where
every other pool is listed. `srs` carries Mol-1's `ranking_score` so the column means the same
thing it does for the OB0 pools (the model's own pose score), and `lig_fe_min` carries the
PoseBusters-measured Fe-ligand distance, which is what the Type-II gate reads.

Only pb_pass poses are listed: MM-GBSA on a pose that cannot be shipped is wasted CPU.

Env: any python3 (stdlib only). CPU, instant.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

R = Path("/vast.mnt/project/bme-20255165-cofolding-physics")
FIELDS = ["pose_id", "pdb", "smiles", "pool", "query_id", "base_id", "seed", "sample",
          "srs", "lig_fe_min"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--qc", type=Path, nargs="+", required=True, help="poses_qc_<pool>.csv files")
    ap.add_argument("--out", type=Path, default=R / "data/cyp3a4_ob0_mmgbsa/poselist_mol1.csv")
    a = ap.parse_args()

    rows, seen = [], set()
    for q in a.qc:
        if not q.is_file():
            sys.exit(f"FATAL: no QC table at {q}")
        n_in, n_kept = 0, 0
        for r in csv.DictReader(q.open(newline="")):
            n_in += 1
            if str(r.get("pb_pass", "")).lower() != "true":
                continue
            if r["pose_id"] in seen:
                sys.exit(f"FATAL: duplicate pose_id {r['pose_id']} across the QC tables")
            seen.add(r["pose_id"])
            # fe_lig_min_A is PoseBusters' own measurement on the exact file to be shipped; the
            # collector's lig_fe_min is the fallback only if the PB column is absent.
            fe = r.get("fe_lig_min_A") or r.get("lig_fe_min") or ""
            rows.append({"pose_id": r["pose_id"], "pdb": r["pdb"], "smiles": r["smiles"],
                         "pool": r["pool"], "query_id": r["query_id"], "base_id": r["base_id"],
                         "seed": r["seed"], "sample": r["sample"],
                         "srs": r.get("ranking_score", ""), "lig_fe_min": fe})
            n_kept += 1
        print(f"  {q.name}: {n_in} poses in, {n_kept} pb_pass kept")
    if not rows:
        sys.exit("FATAL: no PoseBusters-passing pose in any QC table")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    per: dict[str, int] = {}
    for r in rows:
        per[r["pool"]] = per.get(r["pool"], 0) + 1
    print(f"{a.out}: {len(rows)} poses, {len({r['base_id'] for r in rows})} ligands, pools {per}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
