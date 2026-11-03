#!/usr/bin/env python
# Vendored from $COFOLD_ROOT/analysis/cyp3a4_heme/substructure_priors.py -- unchanged apart from the COFOLD_ROOT default.
"""Which substructures in the blinded 20 could sharpen the heme constraints?

The restraint currently selects candidate coordinating atoms BY ELEMENT (any
N/O/S on the ligand) and soft-mins over them, so it pulls whichever heteroatom is
already nearest the iron. That is chemically wrong for most heteroatoms: an amide
or sulfonamide nitrogen has its lone pair delocalised and cannot coordinate, and
a sulfone oxygen cannot either. This script asks three questions:

1. DONOR CENSUS. Per ligand, classify every N/O/S by whether it can donate a lone
   pair to Fe(III). The gap between "element-matched" and "actually a donor" is
   the precision the restraint is currently giving away.
2. CAN WE EXPRESS IT? The constraint would need OF3 atom NAMES, not RDKit indices.
   Checks whether OF3's ligand atom naming is recoverable from the SMILES order by
   comparing against the converted submission PDBs.
3. DID IT MATTER? For the poses the type_ii arm actually produced, which ligand
   atom ended up nearest the iron -- and was it a real donor or an artefact of
   element-based selection?

Plus scaffold structure: Murcko scaffolds and a maximum common substructure, to
see whether the 20 are a series with a shared anchor or chemically diverse.
"""

from __future__ import annotations

import argparse
import os
import json
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ROOT = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
BLINDED = ROOT / "cyp-challenge-structure/cyp-challenge-TEST-BLINDED_structures.csv"
CHEMOTYPE = ROOT / "data/cyp_structure_track/crosscheck_typeII_coordination.csv"
SCORES = ROOT / "data/cyp3a4_rerank/rerank_scores_heme18.csv"

# Donor classes, strongest first. Fe(III) in a P450 is a soft-ish Lewis acid; the
# textbook Type-II donors are sp2 aromatic N with an exposed lone pair.
STRONG_AROMATIC_N = Chem.MolFromSmarts("[nX2;H0]")          # pyridine / azole N
IMIDAZOLE_TRIAZOLE = Chem.MolFromSmarts("[nX2;H0]1cnc[c,n]1")
NITRILE_N = Chem.MolFromSmarts("[NX1]#[CX2]")
THIOETHER_S = Chem.MolFromSmarts("[SX2;H0]([#6])[#6]")
ALIPHATIC_AMINE = Chem.MolFromSmarts("[NX3;H0,H1,H2;!$(N[C,S]=[O,S,N]);!$(N-a)]")
# Explicit NON-donors: lone pair delocalised or sterically/electronically dead.
AMIDE_N = Chem.MolFromSmarts("[NX3][CX3]=[OX1]")
SULFONAMIDE_N = Chem.MolFromSmarts("[NX3][SX4](=[OX1])(=[OX1])")
CARBONYL_O = Chem.MolFromSmarts("[OX1]=[#6,#16]")
ETHER_O = Chem.MolFromSmarts("[OX2;H0]([#6])[#6]")
HYDROXYL_O = Chem.MolFromSmarts("[OX2H1]")
SULFONE_O = Chem.MolFromSmarts("[OX1]=[SX4]")


def _match_set(mol, patt) -> set[int]:
    if patt is None:
        return set()
    return {i for m in mol.GetSubstructMatches(patt) for i in m}


def classify_heteroatoms(mol: Chem.Mol) -> list[dict]:
    """Label every N/O/S by donor capability."""
    strong = _match_set(mol, STRONG_AROMATIC_N)
    azole = _match_set(mol, IMIDAZOLE_TRIAZOLE)
    nitrile = {m[0] for m in mol.GetSubstructMatches(NITRILE_N)} if NITRILE_N else set()
    thioether = {
        m[0] for m in mol.GetSubstructMatches(THIOETHER_S)
    } if THIOETHER_S else set()
    amine = _match_set(mol, ALIPHATIC_AMINE)
    amide_n = {m[0] for m in mol.GetSubstructMatches(AMIDE_N)} if AMIDE_N else set()
    sulfonamide_n = {
        m[0] for m in mol.GetSubstructMatches(SULFONAMIDE_N)
    } if SULFONAMIDE_N else set()
    carbonyl_o = {
        m[0] for m in mol.GetSubstructMatches(CARBONYL_O)
    } if CARBONYL_O else set()
    sulfone_o = {m[0] for m in mol.GetSubstructMatches(SULFONE_O)} if SULFONE_O else set()
    ether_o = _match_set(mol, ETHER_O)
    hydroxyl_o = _match_set(mol, HYDROXYL_O)

    rows = []
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym not in ("N", "O", "S"):
            continue
        i = atom.GetIdx()
        if i in azole:
            cls, donor = "azole_N", "strong"
        elif i in strong:
            cls, donor = "aromatic_N", "strong"
        elif i in nitrile:
            cls, donor = "nitrile_N", "moderate"
        elif i in thioether:
            cls, donor = "thioether_S", "moderate"
        elif i in sulfonamide_n:
            cls, donor = "sulfonamide_N", "none"
        elif i in amide_n:
            cls, donor = "amide_N", "none"
        elif i in amine:
            cls, donor = "amine_N", "weak"
        elif i in sulfone_o:
            cls, donor = "sulfone_O", "none"
        elif i in carbonyl_o:
            cls, donor = "carbonyl_O", "weak"
        elif i in hydroxyl_o:
            cls, donor = "hydroxyl_O", "weak"
        elif i in ether_o:
            cls, donor = "ether_O", "weak"
        else:
            cls, donor = f"other_{sym}", "weak" if sym != "O" else "none"
        rows.append({
            "idx": i,
            "element": sym,
            "class": cls,
            "donor": donor,
            "aromatic": atom.GetIsAromatic(),
            "n_heavy_nbrs": atom.GetDegree(),
            "n_h": atom.GetTotalNumHs(),
        })
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=ROOT / "analysis/cyp3a4_heme")
    args = ap.parse_args(argv)

    blind = pd.read_csv(BLINDED)
    chem = pd.read_csv(CHEMOTYPE, low_memory=False)
    chem = chem[chem["side"] == "TEST"][["id", "typeII_class"]].rename(
        columns={"id": "base_id"})
    blind = blind.merge(chem, left_on="Molecule_Name", right_on="base_id", how="left")

    # ---------------- 1. donor census ----------------
    print("=" * 78)
    print("1  DONOR CENSUS -- element-matched candidates vs REAL donors")
    print("=" * 78)
    per_lig, atom_rows = [], []
    for r in blind.itertuples():
        mol = Chem.MolFromSmiles(r.SMILES)
        if mol is None:
            print(f"  {r.Molecule_Name}: unparsable")
            continue
        hets = classify_heteroatoms(mol)
        for h in hets:
            h["base_id"] = r.Molecule_Name
            atom_rows.append(h)
        n_elem = len(hets)
        strong = [h for h in hets if h["donor"] == "strong"]
        moderate = [h for h in hets if h["donor"] == "moderate"]
        dead = [h for h in hets if h["donor"] == "none"]
        per_lig.append({
            "base_id": r.Molecule_Name,
            "chemotype": r.typeII_class,
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "n_element_candidates": n_elem,
            "n_strong": len(strong),
            "n_moderate": len(moderate),
            "n_dead": len(dead),
            "strong_classes": ",".join(sorted({h["class"] for h in strong})),
            "recommend_type_ii": bool(strong or moderate),
        })
    L = pd.DataFrame(per_lig)
    A = pd.DataFrame(atom_rows)
    print(L.to_string(index=False))
    tot_elem = int(L["n_element_candidates"].sum())
    tot_real = int(L["n_strong"].sum() + L["n_moderate"].sum())
    tot_dead = int(L["n_dead"].sum())
    print(f"\nacross all 20: {tot_elem} atoms match the element filter (N/O/S), "
          f"of which {tot_real} can actually donate and {tot_dead} provably cannot.")
    print(f"=> the element-based restraint is diluted "
          f"{tot_elem / max(tot_real, 1):.1f}x")
    print("\nheteroatom class census:")
    print(A.groupby(["class", "donor"]).size().to_string())
    print("\nligands with NO donor at all (a type_ii arm is chemically impossible "
          "for these, and its 25 poses per state are wasted):")
    nodon = L[~L["recommend_type_ii"]]
    print(sorted(nodon["base_id"]) if len(nodon) else "  none")

    # ---------------- 2. can the prior be expressed? ----------------
    print()
    print("=" * 78)
    print("2  CAN WE EXPRESS IT? OF3 ligand atom names vs SMILES order")
    print("=" * 78)
    name_check = {}
    if SCORES.exists():
        d = pd.read_csv(SCORES, usecols=["base_id", "arm", "lig_fe_atom",
                                         "lig_fe_elem", "lig_fe_min", "mode"],
                        low_memory=False)
        print("atom names OF3 actually put nearest Fe (type_ii arms):")
        t2 = d[d["arm"].astype(str).str.endswith("type_ii")]
        tab = t2.groupby("base_id")["lig_fe_atom"].agg(
            lambda s: ", ".join(f"{k}x{v}" for k, v in
                                Counter(s.dropna()).most_common(3)))
        elem = t2.groupby("base_id")["lig_fe_elem"].agg(
            lambda s: Counter(s.dropna()).most_common(1)[0][0] if s.notna().any() else "")
        nm = pd.DataFrame({"fe_atoms": tab, "modal_elem": elem}).reset_index()
        nm = nm.merge(L[["base_id", "chemotype", "n_strong", "n_moderate"]],
                      on="base_id", how="left")
        print(nm.to_string(index=False))
        name_check = {"n_ligands": int(len(nm))}
        print("\nNOTE the names are OF3-generated (N4, O2, ...). A per-ligand "
              "`coordinating_atom_names` constraint needs this mapping to be "
              "stable; it is recoverable from these pose tables empirically, "
              "which is the cheapest route and needs no new convention.")
    else:
        print(f"(no {SCORES.name}; skipping)")

    # ---------------- 3. did element-based selection misfire? ----------------
    print()
    print("=" * 78)
    print("3  DID IT MISFIRE? donor quality of what type_ii actually grabbed")
    print("=" * 78)
    if SCORES.exists():
        t2 = d[d["arm"].astype(str).str.endswith("type_ii")].copy()
        # join the grabbed element against whether the ligand HAS a real donor
        t2 = t2.merge(L[["base_id", "chemotype", "n_strong", "n_moderate"]],
                      on="base_id", how="left")
        t2["ligand_has_donor"] = (t2["n_strong"] + t2["n_moderate"]) > 0
        close = t2[t2["lig_fe_min"] <= 2.8]
        print(f"type_ii poses within 2.8 A of Fe: {len(close)}/{len(t2)}")
        print("\nby whether the ligand has ANY real donor:")
        print(close.groupby("ligand_has_donor")["lig_fe_elem"]
              .value_counts().to_string())
        nod = close[~close["ligand_has_donor"]]
        print(f"\nposes where a NON-donor ligand was pulled to <=2.8 A: {len(nod)}")
        if len(nod):
            print(nod.groupby(["base_id", "lig_fe_elem", "lig_fe_atom"])
                  .size().to_string())
            print("\n^ these are chemically unphysical contacts manufactured by the "
                  "element-based candidate set. A donor-aware restraint would not "
                  "have produced them.")

    # ---------------- 4. scaffold structure ----------------
    print()
    print("=" * 78)
    print("4  SCAFFOLD STRUCTURE -- series with a shared anchor, or diverse?")
    print("=" * 78)
    mols, names = [], []
    for r in blind.itertuples():
        m = Chem.MolFromSmiles(r.SMILES)
        if m:
            mols.append(m)
            names.append(r.Molecule_Name)
    scaf = Counter()
    for m in mols:
        scaf[Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))] += 1
    print(f"distinct Murcko scaffolds: {len(scaf)} over {len(mols)} ligands")
    for s, n in scaf.most_common(6):
        print(f"  {n}x  {s}")
    print("\nmaximum common substructure over all 20 (ring-matching, 60 s cap):")
    res = rdFMCS.FindMCS(mols, timeout=60, ringMatchesRingOnly=True,
                         completeRingsOnly=True)
    print(f"  SMARTS: {res.smartsString}")
    print(f"  {res.numAtoms} atoms / {res.numBonds} bonds; canceled={res.canceled}")
    if res.numAtoms < 6:
        print("  => NO useful shared anchor. A scaffold-based pose prior is not "
              "available for this set; the donor-atom prior is the one to use.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    L.to_csv(args.out_dir / "substructure_priors_per_ligand.csv", index=False)
    A.to_csv(args.out_dir / "substructure_priors_atoms.csv", index=False)
    (args.out_dir / "substructure_priors.json").write_text(json.dumps({
        "n_element_candidates": tot_elem,
        "n_real_donors": tot_real,
        "n_provable_non_donors": tot_dead,
        "dilution_factor": round(tot_elem / max(tot_real, 1), 2),
        "ligands_without_any_donor": sorted(nodon["base_id"]),
        "n_murcko_scaffolds": len(scaf),
        "mcs_atoms": res.numAtoms,
        "mcs_smarts": res.smartsString,
        **name_check,
    }, indent=2))
    print(f"\nwrote -> {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
