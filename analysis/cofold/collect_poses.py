"""Turn a Mol-1 output tree into the campaign's pose-table + PDB layout.

Mol-1 writes `<dump>/<name>/seed_<seed>/predictions/<name>_sample_<k>.cif` plus a matching
`<name>_summary_confidence_sample_<k>.json`. Everything downstream in this campaign -- PoseBusters
(analysis/cyp3a4_gamd/validate_frames.py), the MM-GBSA prep, and build_submission -- consumes PDB
files whose ligand residue is named exactly `LIG`, so this converts and renames as it collects.

The target layout is taken from an existing validated submission member, not invented:
chain A protein as ATOM, chain B `HEM`, chain C `LIG`, each cofactor at residue number 1.
Mol-1 names a SMILES ligand `l01`..`l99` (mol/data/inference/json_parser.py renames by element
occurrence), and build_submission rejects a member whose ligand residue is not exactly LIG, so the
rename happens here rather than being discovered at zip time.

Confidence: `ranking_score` is kept, but it is NOT the recommended pose criterion for this system.
It is 0.8*iptm + 0.2*ptm + 0.5*disorder - 100*has_clash, where `disorder` is always 0.0 and
`has_clash` is a POLYMER-vs-POLYMER test (docs/output_format.md) -- with one protein chain it can
never fire, so ranking_score carries no protein-ligand clash information at all. The repo's own
guidance is that `chain_pair_iptm` between the chains of interest beats any whole-structure number,
so `cp_iptm_pl` (protein x ligand) is extracted alongside and is what the confidence arm ranks on.

Env: the Mol-1 venv (gemmi). CPU.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import gemmi
import numpy as np

R = Path("/vast.mnt/project/bme-20255165-cofolding-physics")
BLINDED = R / "cyp-challenge-structure/cyp-challenge-TEST-BLINDED_structures.csv"
SAMPLE_RE = re.compile(r"_sample_(\d+)\.cif$")
FIELDS = ["pose_id", "pdb", "cif", "smiles", "pool", "query_id", "base_id", "seed", "sample",
          "ranking_score", "cp_iptm_pl", "iptm", "ptm", "plddt", "gpde", "has_clash",
          "lig_fe_min", "lig_fe_elem"]


def to_pdb(cif: Path, out: Path) -> tuple[float, str, int]:
    """Write the pose as a submission-shaped PDB; return (Fe-ligand min dist, that atom, n_lig)."""
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    model = st[0]
    lig_names = set()
    for chain in model:
        for res in chain:
            # Everything that is not a standard amino acid and not the heme is the drug ligand.
            if res.name != "HEM" and not gemmi.find_tabulated_residue(res.name).is_amino_acid():
                lig_names.add(res.name)
    if len(lig_names) != 1:
        sys.exit(f"FATAL: {cif} has {sorted(lig_names)} as non-HEM heteroresidues, expected exactly one")
    lig_name = lig_names.pop()
    for chain in model:
        for res in chain:
            if res.name == lig_name:
                res.name = "LIG"
                res.seqid.num = 1
            elif res.name == "HEM":
                res.seqid.num = 1

    fe = [a.pos for ch in model for r in ch if r.name == "HEM" for a in r if a.element.name == "Fe"]
    lig = [(a.pos, a.element.name) for ch in model for r in ch if r.name == "LIG" for a in r]
    if not fe:
        sys.exit(f"FATAL: {cif} has no Fe in HEM")
    if not lig:
        sys.exit(f"FATAL: {cif} has no LIG atoms after renaming from {lig_name}")
    d = [(fe[0].dist(p), el) for p, el in lig]
    fe_min, fe_el = min(d, key=lambda t: t[0])

    st.setup_entities()
    out.parent.mkdir(parents=True, exist_ok=True)
    st.write_pdb(str(out))
    return fe_min, fe_el, len(lig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--dump-dir", type=Path, required=True, help="Mol-1 --dump_dir for one arm")
    ap.add_argument("--pool", required=True, help="pool label, e.g. mol1_full / mol1_msa64")
    ap.add_argument("--pdb-root", type=Path, required=True, help="where the converted PDBs go")
    ap.add_argument("--out-csv", type=Path, required=True)
    a = ap.parse_args()

    with BLINDED.open(newline="") as fh:
        smiles = {r["Molecule_Name"]: r["SMILES"] for r in csv.DictReader(fh)}

    cifs = sorted(a.dump_dir.glob("*/seed_*/predictions/*_sample_*.cif"))
    if not cifs:
        sys.exit(f"FATAL: no prediction CIFs under {a.dump_dir}")
    err = a.dump_dir / "ERR"
    if err.is_dir() and any(err.iterdir()):
        # A Mol-1 run can finish having produced fewer structures than asked for; the per-failure
        # reports are the only record, so they are surfaced rather than left on disk.
        print(f"WARNING: {len(list(err.iterdir()))} failure report(s) in {err}")

    rows, bad = [], []
    for cif in cifs:
        base_id = cif.parent.parent.parent.name
        seed = int(cif.parent.parent.name.split("_")[-1])
        k = int(SAMPLE_RE.search(cif.name).group(1))
        js = cif.with_name(cif.name.replace(f"_sample_{k}.cif",
                                            f"_summary_confidence_sample_{k}.json"))
        if not js.is_file():
            bad.append(f"{cif.name}: no summary confidence JSON")
            continue
        c = json.loads(js.read_text())
        cp = c.get("chain_pair_iptm")
        if cp is None:
            bad.append(f"{js.name}: no chain_pair_iptm")
            continue
        # Entity order in the input JSON fixes these indices: protein 0, HEM 1, ligand 2.
        cp_pl = float(np.asarray(cp)[0][2])
        pose_id = f"{a.pool}_{base_id}_s{seed}_k{k}"
        pdb = a.pdb_root / a.pool / f"{pose_id}.pdb"
        fe_min, fe_el, n_lig = to_pdb(cif, pdb)
        rows.append({"pose_id": pose_id, "pdb": str(pdb), "cif": str(cif),
                     "smiles": smiles[base_id], "pool": a.pool, "query_id": base_id,
                     "base_id": base_id, "seed": seed, "sample": k,
                     "ranking_score": c.get("ranking_score"), "cp_iptm_pl": cp_pl,
                     "iptm": c.get("iptm"), "ptm": c.get("ptm"), "plddt": c.get("plddt"),
                     "gpde": c.get("gpde"), "has_clash": c.get("has_clash"),
                     "lig_fe_min": round(fe_min, 3), "lig_fe_elem": fe_el})
    if bad:
        print(f"{len(bad)} pose(s) skipped, e.g.:", *bad[:5], sep="\n  ")

    a.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with a.out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    ligs = {r["base_id"] for r in rows}
    coord = sum(1 for r in rows if r["lig_fe_min"] <= 2.6)
    print(f"{a.out_csv}: {len(rows)} poses, {len(ligs)} ligands, "
          f"seeds {sorted({r['seed'] for r in rows})}")
    print(f"  Fe-ligand <= 2.6 A (the coordination band): {coord}/{len(rows)}")
    print(f"  PDBs -> {a.pdb_root / a.pool}")
    return 0 if rows and not bad else 1


if __name__ == "__main__":
    sys.exit(main())
