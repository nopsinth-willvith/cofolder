"""Score one submission against another with the CHALLENGE'S OWN lDDT-PLI scorer.

NOT AN ACCURACY MEASUREMENT. lDDT-PLI needs a reference structure and the 20 blinded ligands
have none -- that is what makes them blinded, and using a PDB entry as an oracle is ruled out
for this campaign. So the reference here is ANOTHER SUBMISSION, and the number answers:

    "how differently does submission A place the ligand, in the contact space the portal
     actually scores, than submission B does?"

It cannot say which is closer to the truth. Its value is that lDDT-PLI weights protein-ligand
CONTACTS rather than coordinates, so two poses can sit 3 A apart in RMSD and still score alike
if they preserve the same contacts -- or sit 1 A apart and score very differently if they do
not. Against a reference whose portal score is KNOWN (submission_eint_templatepool scored
lDDT-PLI 0.44), the divergence is at least anchored to a measured quantity.

Uses analysis/cyp_structure_track/vendor/evaluation/evaluate_predictions.score_single_structure
-- the challenge's own code -- rather than a reimplementation, because OST re-perceives ligand
connectivity from 3D distances (1.95 A) and this project has already produced one FALSE
"broken ligand geometry" conclusion from a hand-rolled OST call.

Env: openmm83 or plinder (the only envs with openstructure). CPU.
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(os.environ.get("COFOLD_ROOT", "/project/bme-20255165-cofolding-physics"))
# the evaluator uses RELATIVE imports (from .config import ...) and evaluation/ has no
# __init__.py, so it is a PEP 420 namespace package: put its PARENT on the path and import it
# as `evaluation.evaluate_predictions`. Importing the bare file raises
# "attempted relative import with no known parent package", and adding __init__.py would be
# the wrong fix -- that is exactly what broke the OF3 overlay earlier in this campaign.
# NEITHER vendored copy of the challenge evaluator is complete:
#   vendor/evaluation/            config.py lacks BISYRMSD_NAN_PENALTY that its own
#                                 evaluate_predictions.py imports
#   vendor/blind-challenge-.../backend/   defines that, but lacks custom_scoring_functions.py
#                                 which its config.py imports
# Each is missing what the other has. Rather than splice two mismatched vintages together --
# and risk scoring with a silently different metric than the portal -- use the project's OWN
# wrapper, analysis/cyp_structure_track/ost_ligand_probe.score_pair, whose docstring states it
# IS the official call: SCRMSD + LDDT-PLI, substructure_match=True, coverage_delta=0.1. It
# imports ost directly and has no package tangle.
sys.path.insert(0, str(ROOT / "analysis/cyp_structure_track"))
# ost_ligand_probe does `import posebusters_utils` at module level; that file lives
# in the vendored evaluation tree (the one whose config.py is stale -- but this
# single module is self-contained and importable on its own).
sys.path.insert(0, str(ROOT / "analysis/cyp_structure_track/vendor/evaluation"))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--model-zip", type=Path)
    src.add_argument("--pose-csv", type=Path,
                     help="instead of a zip: score every row's `pdb` against the ref-zip member "
                          "<base_id>.pdb; the output keeps the row's own columns")
    ap.add_argument("--ref-zip", type=Path, required=True)
    ap.add_argument("--out-csv", type=Path, required=True)
    ap.add_argument("--portal", action="store_true",
                    help="also score the portal's way (SMILES-built ligands, LDDT-LP from the "
                         "SCRMSD aux) via ost_ligand_probe.score_pair_portal -> portal_* columns")
    a = ap.parse_args(argv)

    import ost
    from ost import conop
    from ost_ligand_probe import score_pair, score_pair_portal
    clib = conop.GetDefaultLib()
    if clib is None:
        raise SystemExit("OST has no default compound library; RuleBasedProcessor needs one")

    import csv as _csv
    smiles = {}
    bl = ROOT / "cyp-challenge-structure/cyp-challenge-TEST-BLINDED_structures.csv"
    if bl.is_file():
        with open(bl, newline="") as fh:
            smiles = {r["Molecule_Name"]: r["SMILES"] for r in _csv.DictReader(fh)}
    RZ = zipfile.ZipFile(a.ref_zip)
    # (row seed, model-bytes reader) per scored pose; the reference is always <base_id>.pdb
    if a.model_zip:
        MZ = zipfile.ZipFile(a.model_zip)
        items = [({"base_id": n[:-4]}, lambda n=n: MZ.read(n))
                 for n in sorted(n for n in MZ.namelist() if n.endswith(".pdb"))]
        label = a.model_zip.name
    else:
        items = [(r, lambda p=r["pdb"]: Path(p).read_bytes())
                 for r in _csv.DictReader(open(a.pose_csv, newline=""))]
        label = a.pose_csv.name
    rnames = {n for n in RZ.namelist() if n.endswith(".pdb")}
    missing = sorted({s["base_id"] for s, _ in items if f"{s['base_id']}.pdb" not in rnames})
    if missing:
        raise SystemExit(f"reference lacks {len(missing)} members: {missing[:3]}")
    print(f"model {label}\nref   {a.ref_zip.name}\n{len(items)} poses\n", flush=True)

    rows = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for i, (seed, read) in enumerate(items):
            base = seed["base_id"]
            rec = {**seed, "error": ""}
            mp, rp = td / f"m_{i}.pdb", td / f"r_{i}.pdb"
            try:
                mp.write_bytes(read()); rp.write_bytes(RZ.read(f"{base}.pdb"))
                # residue_select: score the LIG residue as OST perceives it. The SMILES path
                # needs posebusters_utils from the broken vendor tree; residue_select is what
                # the template backend itself uses.
                rec.update(score_pair(str(mp), str(rp), smiles.get(base, ""),
                                      "residue_select", clib))
                if a.portal:
                    p = score_pair_portal(str(mp), str(rp), smiles.get(base, ""), "smiles", clib)
                    if p["state"] == "no_bonded_copy":
                        p = score_pair_portal(str(mp), str(rp), "", "residue_select", clib)
                    rec.update({f"portal_{k}": v for k, v in p.items()})
            except Exception as e:  # noqa: BLE001 -- a row, never a silent gap
                rec["error"] = f"{type(e).__name__}: {e}"[:160]
            rows.append(rec)
            # `v == v` alone is a NaN test, and None PASSES it, so an unassigned pair (score
            # None, no exception) reached f"{None:.3f}" and killed the run after every row was
            # collected but before the CSV was written -- losing all results. Same isinstance
            # guard the summary block below already uses, and the OST state says why.
            v = rec.get("lddtpli_score")
            if isinstance(v, (int, float)) and v == v:
                note = f"{v:.3f}"
            else:
                note = (rec["error"][:60]
                        or f"unassigned ({rec.get('lddtpli_state_desc', '')})"[:60])
            print(f"  {rec['base_id']}  {note}", flush=True)

    import csv
    cols, seen = [], set()
    for r in rows:
        for k in r:
            if k not in seen:
                seen.add(k); cols.append(k)
    a.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow(r)

    ok = [r for r in rows if not r["error"]]
    # an unassigned pair raises nothing, so it would otherwise be counted as scored
    nun = sum(1 for r in rows if not r["error"]
              and not isinstance(r.get("lddtpli_score"), (int, float)))
    print(f"\nscored {len(ok)}/{len(rows)}"
          + (f"  ({nun} returned no assignment)" if nun else ""))
    if ok:
        import statistics as st
        for k in cols:
            vals = [r[k] for r in ok if isinstance(r.get(k), (int, float))
                    and r[k] == r[k]]
            if len(vals) > 1:
                print(f"  {k:28s} median {st.median(vals):.3f}  "
                      f"min {min(vals):.3f}  max {max(vals):.3f}")
    bad = [r for r in rows if r["error"]]
    for r in bad[:5]:
        print(f"  FAILED {r['base_id']}: {r['error']}")
    print(f"\nwrote {a.out_csv}")
    print("\nREMINDER: the reference is another PREDICTION, not ground truth. "
          "A low score means the two disagree, not that either is wrong.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
