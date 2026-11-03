"""Add the repaired OCNT-2313038 poses to a candidate table, with honest per-pose numbers.

Mol-1 built a genuinely planar cyclohexane for this ligand in all 50 poses (ring_check.py), so
every one failed PoseBusters and the ligand had no shippable member at all. The MM-GBSA
minimisation repairs the ring; repair_ligand.py puts that ligand back into its own pose and
PoseBusters then passes it. These rows carry those repaired files into the ranking.

What comes from where, and why:

* `cp_iptm_pl`, `ranking_score`, `iptm`, `ptm`, `plddt` are INHERITED from the pose the repair was
  made from. They describe Mol-1's prediction, and the repair did not re-run Mol-1 -- it fixed one
  ring in that same prediction.
* `fe_lig_in_A` is RE-MEASURED by PoseBusters on the repaired file, never inherited. The
  minimisation carries the Fe contact by LJ and electrostatics alone (`--no-fe-bond`), so the donor
  relaxes outward: these poses sit at 2.40-3.50 A where the originals were 2.195-2.230 A. Using the
  original distance here would smuggle a coordination the shipped file does not have.
* `pool` becomes `mol1_repair` so a repaired member can never be mistaken for a raw Mol-1 pose in
  the provenance, and `repaired_from` records the parent.

Env: sigmadock. CPU, seconds.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

from build_submission import read_scores  # same directory; concatenation + duplicate-pose guard

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
INHERIT = ["smiles", "query_id", "base_id", "seed", "sample", "ranking_score", "cp_iptm_pl",
           "iptm", "ptm", "plddt", "gpde", "has_clash"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--table", type=Path, nargs="+", required=True,
                    help="candidate table(s) to extend (several are concatenated)")
    ap.add_argument("--qc", type=Path, default=R / "data/cyp3a4_mol1/poses_qc_mol1_full.csv")
    ap.add_argument("--repaired", type=Path, default=R / "data/cyp3a4_mol1/repaired_mol1fix.csv")
    ap.add_argument("--pb", type=Path, default=R / "data/cyp3a4_mol1/pb_mol1repair.csv")
    ap.add_argument("--scores", type=Path,
                    help="optional per-pose energy table keyed on the PARENT pose_id (MM-GBSA)")
    ap.add_argument("--score-cols", default="delta_total,fe_lig_in_A,fe_lig_out_A,"
                                            "min_lig_rmsd_A,vdw,eel,egb,esurf")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    base = read_scores(a.table)
    rep = pd.read_csv(a.repaired)
    pb = pd.read_csv(a.pb).rename(columns={"run": "pose_id"})
    pb["ok"] = pb.pb_pass.astype(str).str.lower().eq("true")
    pb = pb[pb.ok]
    if pb.empty:
        sys.exit("FATAL: no repaired pose passes PoseBusters -- nothing to add")
    qc = pd.read_csv(a.qc).set_index("pose_id")

    rows = []
    for _, r in rep.iterrows():
        hit = pb[pb.pose_id == r.pose_id]
        if hit.empty:
            continue                                   # repaired but still failing: not shippable
        parent = r["repaired_from"]
        if parent not in qc.index:
            sys.exit(f"FATAL: parent {parent} absent from {a.qc}")
        src = qc.loc[parent]
        row = {c: src.get(c) for c in INHERIT}
        row.update({"pose_id": r.pose_id, "pdb": r["pdb"], "pool": "mol1_repair", "ckpt": "mol1",
                    "status": "ok", "srs": src.get("ranking_score"),
                    "repaired_from": parent,
                    "fe_lig_in_A": float(hit.iloc[0].fe_lig_min_A)})
        if a.scores:
            s = pd.read_csv(a.scores)
            s = s[s.pose_id == parent]
            if s.empty:
                sys.exit(f"FATAL: --scores has no row for parent {parent}")
            for c in a.score_cols.split(","):
                if c in s.columns and c != "fe_lig_in_A":
                    row[c] = s.iloc[0][c]
        rows.append(row)
    if not rows:
        sys.exit("FATAL: no repaired pose survived both the repair and PoseBusters")

    add = pd.DataFrame(rows)
    out = pd.concat([base, add], ignore_index=True)
    out.to_csv(a.out, index=False)
    t2 = set((R / "data/cyp3a4_xmodel/typeii_base_ids.txt").read_text().split())
    lig = add.base_id.iloc[0]
    coord = int((add.fe_lig_in_A <= 2.6).sum())
    print(f"{a.out}: {len(base)} + {len(add)} repaired = {len(out)} poses, "
          f"{out.base_id.nunique()} ligands")
    print(f"  {lig} ({'Type-II' if lig in t2 else 'non-Type-II'}): {len(add)} repaired poses, "
          f"{coord} inside the 2.6 A coordination band")
    if lig in t2 and coord == 0:
        print("  ⚠ none coordinate at 2.6 A -- the Type-II gate will need a relaxed cut here")
    return 0


if __name__ == "__main__":
    sys.exit(main())
