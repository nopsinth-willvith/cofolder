"""Run THE challenge's own submission validator -- the one structure_prediction.ipynb calls.

The notebook's cell is:

    from validation.structure_validation import validate_structure_submission
    is_valid, errors = validate_structure_submission(
        zip_path, expected_ids=..., expected_ligand_smiles=...,
        require_lig_resname=True, check_numbering=True)

so that is what runs here, against the vendored copy, with the same arguments.

WORTH KNOWING BEFORE READING THE RESULT: this validator is LOOSER than this project's own
analysis/cyp3a4_relax/validate_submission_zip.py. It does NOT require a 503-residue protein
numbered 1..503 -- that check is the project's own belt-and-braces addition. What it actually
requires is: every expected id present; exactly one residue named LIG; at most 3 chains; ligand
connectivity that AssignBondOrdersFromTemplate accepts against the expected SMILES; and clean
numbering (no duplicate residue numbers per chain, no duplicate atom names in a residue, no
duplicate atom serials file-wide).

Passing it is therefore NOT the same as being a good submission: a truncated protein can satisfy
every one of these checks and still score near zero on lDDT-LP, which is scored on the receptor.
"""
from __future__ import annotations

import argparse
import csv
import sys
import warnings
from pathlib import Path

ROOT = Path("/project/bme-20255165-cofolding-physics")
VENDOR = ROOT / "analysis/cyp_structure_track/vendor"
sys.path.insert(0, str(VENDOR))

BLINDED = ROOT / "cyp-challenge-structure/cyp-challenge-TEST-BLINDED_structures.csv"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zip", dest="zip_path", type=Path, required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--smiles-csv", type=Path, default=BLINDED)
    args = ap.parse_args(argv)

    from rdkit import RDLogger
    from validation.structure_validation import validate_structure_submission

    with args.smiles_csv.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    ids = {r["Molecule_Name"] for r in rows}
    smiles = {r["Molecule_Name"]: r["SMILES"] for r in rows}

    print(f"\n=== OFFICIAL validator: {args.label or args.zip_path.name} ===")
    print(f"zip {args.zip_path} | expected ids {len(ids)}")
    RDLogger.DisableLog("rdApp.warning")
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module=r"MDAnalysis.*")
            warnings.filterwarnings("ignore", category=DeprecationWarning, module=r"MDAnalysis.*")
            ok, errors = validate_structure_submission(
                args.zip_path, expected_ids=ids, expected_ligand_smiles=smiles,
                require_lig_resname=True, check_numbering=True)
    finally:
        RDLogger.EnableLog("rdApp.warning")

    if ok:
        print("VALID -- the notebook would accept this zip")
    else:
        print(f"INVALID -- {len(errors)} error(s):")
        for m in errors[:25]:
            print(f"   - {m}")
        if len(errors) > 25:
            print(f"   ... and {len(errors) - 25} more")
    return 0
if __name__ == "__main__":
    raise SystemExit(main())
