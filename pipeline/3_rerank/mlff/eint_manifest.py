# Vendored from $COFOLD_ROOT/analysis/cyp3a4_mol1/eint_manifest.py -- unchanged apart from the COFOLD_ROOT default and --only.
"""E_int manifest for the Type-II Mol-1 poses, straight off Mol-1's own mmCIF output.

The hybrid submission ranks the 11 Type-II ligands by E_int and the 9 non-Type-II by MM-GBSA, so
only the Type-II poses need scoring and only they are listed here.

Mol-1 writes mmCIF natively, which matters: prep_cofold_system.prep_one() reads mmCIF through
gemmi.cif.read and a .pdb does not parse at all. The OB0 arms had to carry a separate cif column
through their QC tables for this; here the prediction file IS the CIF, so no second copy and no
conversion step. Verified chain layout in those files: A = 503 amino acids, B = HEM, C = `l01`,
which is why lig_res is l01 rather than LIG (LIG is the name used only in the submission PDBs).

`spin_complex` is per pose, not per ligand, and it is the whole reason this manifest exists rather
than a flat list: a ferric heme is LOW-SPIN (multiplicity 2) when a ligand N occupies the sixth
coordination site and HIGH-SPIN (6) when the site is empty. The split is taken from the
PoseBusters-measured Fe-ligand distance against the campaign's 2.6 A band -- the 122-crystal prior
puts coordinating and non-coordinating contacts in disjoint bands either side of 2.58/2.72 A -- so
a Type-II pose whose donor drifted out of the band is prepped as the 5-coordinate species it
actually is, not as the one its ligand chemotype implies.

Env: sigmadock (pandas only). CPU, seconds.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
TYPEII = R / "data/cyp3a4_xmodel/typeii_base_ids.txt"
FE_BAND = 2.6
# eint_prep_shard.py reads tag/arm/query_id/seed/sample/base_id/spin_complex/fe_coordinating by
# ATTRIBUTE off the itertuples row (line 58), so a missing one is an AttributeError mid-shard, not a
# clear error up front. `arm` is this pool's label and `query_id` equals base_id for Mol-1, whose
# queries are not stereo-expanded.
COLS = ["tag", "arm", "cif", "smiles", "lig_res", "spin_complex", "fe_coordinating",
        "pose_id", "query_id", "base_id", "pool", "seed", "sample", "fe_lig_min_A"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--qc", type=Path, nargs="+",
                    default=[R / "data/cyp3a4_mol1/poses_qc_mol1_full.csv"])
    ap.add_argument("--out", type=Path, default=R / "data/cyp3a4_mol1/eint_manifest_typeii.csv")
    ap.add_argument("--lig-res", default="l01")
    ap.add_argument("--all-ligands", action="store_true",
                    help="list non-Type-II poses too (the hybrid does not need them)")
    ap.add_argument("--only", type=Path, default=None,
                    help="CSV with a pose_id column (the MM-GBSA shortlist): list only those poses; "
                         "a listed pose absent from the QC tables is fatal")
    a = ap.parse_args()

    t2 = set(TYPEII.read_text().split())
    frames = []
    for q in a.qc:
        if not q.is_file():
            sys.exit(f"FATAL: no QC table at {q}")
        d = pd.read_csv(q)
        d["pb_pass"] = d.pb_pass.astype(str).str.lower().eq("true")
        frames.append(d)
    df = pd.concat(frames, ignore_index=True)
    n_all = len(df)
    df = df[df.pb_pass]
    if a.only:
        want = set(pd.read_csv(a.only).pose_id)
        lost = sorted(want - set(df.pose_id))
        if lost:
            sys.exit(f"FATAL: {len(lost)} shortlisted pose(s) have no pb_pass row in the QC "
                     f"tables, e.g. {lost[:3]}")
        df = df[df.pose_id.isin(want)]
    if not a.all_ligands:
        df = df[df.base_id.isin(t2)]
    if df.empty:
        sys.exit("FATAL: no PoseBusters-passing Type-II pose to score")

    missing = [c for c in ("cif", "fe_lig_min_A", "smiles") if c not in df.columns]
    if missing:
        sys.exit(f"FATAL: QC table lacks {missing}")
    bad = df[~df.cif.map(lambda p: Path(p).is_file())]
    if len(bad):
        sys.exit(f"FATAL: {len(bad)} pose(s) name a CIF that is not on disk, e.g. "
                 f"{bad.cif.iloc[0]}")

    fe = pd.to_numeric(df.fe_lig_min_A, errors="coerce")
    coord = fe <= FE_BAND
    out = pd.DataFrame({
        "tag": df.pose_id, "cif": df.cif, "smiles": df.smiles, "lig_res": a.lig_res,
        # ferric heme: N-ligated -> doublet, 5-coordinate -> sextet
        "spin_complex": coord.map({True: 2, False: 6}).astype(int),
        "fe_coordinating": coord.astype(bool),
        "arm": df.pool, "pose_id": df.pose_id, "query_id": df.base_id,
        "base_id": df.base_id, "pool": df.pool,
        "seed": df.seed, "sample": df["sample"], "fe_lig_min_A": fe,
    })[COLS]
    if out.tag.duplicated().any():
        sys.exit(f"FATAL: duplicate tag, e.g. {out.tag[out.tag.duplicated()].iloc[0]}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)

    print(f"{a.out}: {len(out)} poses of {n_all} ({out.base_id.nunique()} ligands)")
    print(f"  multiplicity split: {int((out.spin_complex == 2).sum())} low-spin (N-ligated, "
          f"Fe <= {FE_BAND} A), {int((out.spin_complex == 6).sum())} high-spin (5-coordinate)")
    thin = out.groupby("base_id").size()
    thin = thin[thin < 10]
    if len(thin):
        print(f"  thin: {thin.to_dict()}")
    absent = sorted(t2 - set(out.base_id))
    if absent:
        # A Type-II ligand with nothing here cannot be ranked by E_int at all, so the hybrid would
        # have to fall back to MM-GBSA for it -- better known now than discovered at build time.
        print(f"  ⚠ no scorable pose for {absent} -- the hybrid must fall back to MM-GBSA there")
    return 0


if __name__ == "__main__":
    sys.exit(main())
