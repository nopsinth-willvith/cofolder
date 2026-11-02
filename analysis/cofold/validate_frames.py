"""Stage G5 -- run the challenge's validators on an MD final frame.

WHY NOT validate_submission_zip.py. That one re-derives the portal's requirements and asserts a
503-residue chain numbered 1..503 with CYS442:SG. An MD frame from the pocket build is 212-239
residues by construction, so it would fail check 5 for a reason that has nothing to do with the
physics. The checks that DO transfer to a single complex are run here instead.

WHAT IS RUN, per frame:

1. POSEBUSTERS, as the official scorer calls it (vendor evaluate_predictions.py:707-708):
   PoseBusters(config="dock").bust(mol_pred, mol_true, mol_cond), where mol_pred comes from
   bonded_ligand_copies (bond orders assigned from SMILES, not guessed from 3D distance) and
   mol_cond is every non-LIG residue combined, heme and Fe included. pb_pass requires all 22
   checks in EXPECTED_PB to be exactly True -- NaN, missing or non-bool counts as a FAIL, because
   the official scorer only looks at bool columns and would silently skip the others.

2. LONGEST LIGAND C-C BOND, from the SMILES graph. PoseBusters does NOT gate this: its
   bond_lengths check is a percentage tolerance against the reference conformer, and this project
   found +20% stretches (1.73-1.83 A) that pb_pass accepted. Bonds come from the graph, never from
   a distance window -- a window flags four-membered ring diagonals at 2.08-2.19 A as stretched
   bonds, which is how an earlier version of the zip validator failed 8 of 20 good members.

3. Fe-LIGAND GEOMETRY: closest ligand heavy atom to Fe, its element, and whether a non-nitrogen
   atom sits within 2.6 A (flag_nonN_fe), mirroring qc_poses' definition.

TWO LIMITS OF RUNNING THIS ON A TRUNCATED FRAME, which the numbers do not show on their own:
  * minimum_distance_to_protein and the volume-overlap checks only see the kept shell. Every
    residue within contact range of the ligand is present, so a clash would still be caught, but
    the check is not scanning the whole protein.
  * the frames are heavy-atom only (hydrogens were dropped to match the submitted poses), which
    is the same state the portal scores.

Env: sigmadock (posebusters 0.6.5 -- the version qc_poses.py pins). CPU.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

ROOT = Path("/project/bme-20255165-cofolding-physics")
sys.path.insert(0, str(ROOT / "analysis/cyp_structure_track"))
sys.path.insert(0, str(ROOT / "analysis/cyp_structure_track/vendor/evaluation"))

EXPECTED_PB = (
    "mol_pred_loaded", "mol_cond_loaded", "sanitization", "inchi_convertible",
    "all_atoms_connected", "no_radicals", "bond_lengths", "bond_angles",
    "internal_steric_clash", "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy", "protein-ligand_maximum_distance",
    "minimum_distance_to_protein", "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors", "minimum_distance_to_waters",
    "volume_overlap_with_protein", "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors", "volume_overlap_with_waters",
)


def longest_cc_bond(mol) -> float:
    """Longest C-C bond, from the molecular GRAPH. C-S (1.80) and C-Cl (1.75) are normal and
    would pollute the statistic, so only carbon-carbon is measured."""
    conf = mol.GetConformer()
    worst = 0.0
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetSymbol() == "C" and a2.GetSymbol() == "C":
            p, q = conf.GetAtomPosition(a1.GetIdx()), conf.GetAtomPosition(a2.GetIdx())
            worst = max(worst, math.dist((p.x, p.y, p.z), (q.x, q.y, q.z)))
    return worst


def fe_geometry(pdb_path: Path):
    """Closest ligand heavy atom to Fe, and the non-nitrogen flag at 2.6 A."""
    fe, lig = None, []
    for line in pdb_path.read_text().splitlines():
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        name, resn = line[12:16].strip(), line[17:20].strip()
        el = (line[76:78].strip() or name[:1]).upper()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if name == "FE":
            fe = xyz
        elif resn == "LIG" and el != "H":
            lig.append((name, el, xyz))
    if fe is None or not lig:
        return None, None, None, None
    best = min(lig, key=lambda a: math.dist(fe, a[2]))
    d = math.dist(fe, best[2])
    non_n = any(a[1] != "N" and math.dist(fe, a[2]) <= 2.6 for a in lig)
    return round(d, 3), best[0], best[1], non_n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frames-csv", type=Path, required=True,
                    help="base_id,pdb[,run] -- the restored-numbering frames")
    ap.add_argument("--smiles-csv", type=Path,
                    default=ROOT / "cyp-challenge-structure/cyp-challenge-TEST-BLINDED_structures.csv")
    ap.add_argument("--id-col", default="Molecule_Name")
    ap.add_argument("--smiles-col", default="SMILES")
    ap.add_argument("--out-csv", type=Path, required=True)
    args = ap.parse_args(argv)

    import posebusters
    from posebusters import PoseBusters
    from rdkit import Chem, RDLogger
    # posebusters_utils imports loguru, which sigmadock does not have. qc_poses.py already
    # carries the fix (_stub_loguru, comment: "absent from sigmadock"), so reuse it rather than
    # write a second stub that could drift from it.
    from qc_poses import _stub_loguru
    _stub_loguru()
    from posebusters_utils import bonded_ligand_copies
    RDLogger.DisableLog("rdApp.*")
    print(f"posebusters {posebusters.__version__}")

    with args.smiles_csv.open(newline="") as fh:
        smiles = {r[args.id_col]: r[args.smiles_col] for r in csv.DictReader(fh)}

    buster = PoseBusters(config="dock")
    rows = []
    with args.frames_csv.open(newline="") as fh:
        frames = list(csv.DictReader(fh))

    for fr in frames:
        base, pdb = fr["base_id"], Path(fr["pdb"])
        smi = smiles.get(base)
        row = {"base_id": base, "run": fr.get("run", ""), "pdb": str(pdb)}
        if smi is None:
            row["error"] = "no SMILES"
            rows.append(row)
            continue
        try:
            complex_mol = Chem.MolFromPDBFile(str(pdb), removeHs=False, sanitize=False)
            frags = Chem.SplitMolByPDBResidues(complex_mol)
            cond = None
            for resname, frag in frags.items():
                if resname == "LIG":
                    continue
                cond = frag if cond is None else Chem.CombineMols(cond, frag)
            copies = bonded_ligand_copies(str(pdb), smi)
            # 0 copies means distance-based perception could not match the reference molecule. A
            # DISTORTED ligand fragments and lands here: OCNT-2315379's tier-2 pose had a urea C-N
            # stretched to 2.07 A, RDKit split it into 3 pieces, and `[0]` on the empty list raised
            # a bare IndexError that escalate_relax read as a clash verdict. max_cc_bond_A cannot
            # catch it either -- it measures C-C only, and every broken bond here was a C-N.
            # validate_submission_zip already fails a member with != 1 copy, so record it here too
            # rather than discovering it at zip time.
            row["n_bonded_copies"] = len(copies)
            if not copies:
                raise ValueError("0 bonded ligand copies -- perception could not match the "
                                 "reference SMILES (distorted ligand?)")
            pred = copies[0]
            true = Chem.MolFromSmiles(smi)
            res = buster.bust(mol_pred=pred, mol_true=true, mol_cond=cond)
            r0 = res.iloc[0]
            failed, nonbool = [], []
            for c in EXPECTED_PB:
                if c not in res.columns:
                    nonbool.append(c)
                elif not isinstance(r0[c], bool) and r0[c] not in (True, False):
                    nonbool.append(c)
                elif not bool(r0[c]):
                    failed.append(c)
            row["pb_pass"] = (not failed) and (not nonbool)
            row["pb_n_failed"] = len(failed)
            row["pb_failed"] = ";".join(failed)
            row["pb_nonbool"] = ";".join(nonbool)
            row["max_cc_bond_A"] = round(longest_cc_bond(pred), 3)
            d, nm, el, non_n = fe_geometry(pdb)
            row.update({"fe_lig_min_A": d, "fe_lig_atom": nm, "fe_lig_element": el,
                        "flag_nonN_fe": non_n})
        except Exception as exc:  # noqa: BLE001 -- a validator reports, it does not raise
            row["error"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)
        print(f"  {base:<16} pb_pass={row.get('pb_pass')} "
              f"failed={row.get('pb_failed') or '-'} "
              f"maxCC={row.get('max_cc_bond_A')} Fe-lig={row.get('fe_lig_min_A')} "
              f"({row.get('fe_lig_element')}) nonN_fe={row.get('flag_nonN_fe')}"
              + (f"  ERROR {row['error']}" if row.get("error") else ""))

    cols = ["base_id", "run", "pb_pass", "pb_n_failed", "pb_failed", "pb_nonbool",
            "max_cc_bond_A", "n_bonded_copies", "fe_lig_min_A", "fe_lig_atom", "fe_lig_element",
            "flag_nonN_fe", "error", "pdb"]
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    ok = sum(1 for r in rows if r.get("pb_pass"))
    print(f"\npb_pass {ok}/{len(rows)}   ->  {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
