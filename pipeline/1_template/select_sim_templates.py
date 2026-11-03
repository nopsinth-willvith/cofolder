#!/usr/bin/env python3
"""Per-ligand CYP3A4 template selection by BOUND-LIGAND chemical similarity (ECFP4).

PRE-REGISTERED RULE, fixed before any pose was generated:
  1. Candidates are the corpus compounds that are `PUTATIVE_LIGAND` AND `in_pocket_any` --
     101 of the 113 compounds in the similarity matrix. This drops 10 buffer/cryoprotectant
     comps (CIT DMS EDO GOL HEM IMD PG0 PG4 PO4 SO4) and 2 out-of-pocket ones (STR, UC7).
  2. For each blinded ligand, rank candidates by `ecfp4_tanimoto` DESCENDING and keep the
     top-4 DISTINCT COMPOUNDS -- not the top-4 entries, which would hand the model four
     crystals of the same molecule and waste three of Mol-1's four template slots.
  3. Each compound resolves to ONE entry: of the PDB entries holding it, the best (lowest)
     `resolution_A`; ties broken by `pdb_id` ascending. Deterministic.
  4. ECFP4 ties between compounds (OCNT-2535422 is a 4-way tie at 0.125) break by the
     resolved entry's resolution ascending, then `comp_id` ascending. Without this the
     bundle order -- and therefore the selection, see below -- is nondeterministic.

EXACTLY FOUR templates per ligand, and the ORDER IS THE SELECTION. Mol-1 hard-codes
`TemplateFeatureAssemblyLine(max_templates=4)` at mol/data/template/template_featurizer.py:774
and truncates with `np.arange(num_t) < max_templates` -- FIRST-N, no sort. A fifth entry would
be silently dropped, and the .json template route logs nothing at all, so there is no runtime
witness. Hence: build 4, in rank order.

⚠ THE SIMILARITY CEILING IS THE EXPERIMENT'S MAIN LIMITATION, and it is a property of the
corpus rather than of the method. Best-neighbour ECFP4 runs ~0.125-0.340 (median ~0.20), so
"most similar" means "least dissimilar": the blinded ligands are small (MW median 329) while
the CYP3A4 co-crystal corpus is dominated by large ligands such as ritonavir (534). The same
hypothesis measured on OF3/OpenBind-0 as arm `simtmpl50` was the WORST of three template
regimes (69.4% PB-valid vs 81.8% for a fixed pocket-state-diverse bundle; the only pair
significant after Holm, median delta +0.090, p_holm 0.017). This arm re-tests it on a
DIFFERENT model, so a null is informative and a win would be surprising.

Ported from analysis/cyp3a4_finetune/build_simtmpl_queries.py with two deliberate changes:
  - reads `crosscheck_similarity_pairs.csv` (the full 20x113 matrix) instead of
    `crosscheck_nearest_neighbours.csv`, which is TRUNCATED AT RANK 3 and is why the OF3 arm
    only ever got ~3 templates per query;
  - drops that script's `n_chains == 1` filter. It existed because OF3's
    `template_cif_chain_ids` was ambiguous for multi-chain entries; Mol-1's .json route writes
    out only the named chain, so multi-chain entries are usable. This widens the pool from 96
    to 111 entries and keeps 5TE8 / 7UFE, which appear in three ligands' top-4.

Env: sigmadock (stdlib only -- no RDKit needed, the matrix is precomputed). CPU, seconds.
"""
from __future__ import annotations

import argparse
import os
import csv
import sys
from collections import defaultdict
from pathlib import Path

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
CST = R / "data/cyp_structure_track"
SIM = CST / "crosscheck_similarity_pairs.csv"
COMPS = CST / "crosscheck_pdb_ligands.csv"
ENTRIES = CST / "cyp3a4_pdb_entries.csv"
CIF_DIR = CST / "cif"
BLINDED = R / "cyp-challenge-structure/cyp-challenge-TEST-BLINDED_structures.csv"
N_TEMPLATES = 4          # Mol-1's hard cap; see the module docstring
CHAIN = "A"              # every candidate entry has an author chain A


def truthy(v: str) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes", "t")


def rows(p: Path) -> list[dict]:
    if not p.is_file():
        sys.exit(f"FATAL: missing {p}")
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-templates", type=int, default=N_TEMPLATES)
    ap.add_argument("--out-manifest", type=Path,
                    default=R / "data/cyp3a4_mol1/templates/sim4_manifest.csv")
    ap.add_argument("--out-list", type=Path,
                    default=R / "data/cyp3a4_mol1/templates/sim4_bundles.txt",
                    help="one line per ligand: '<base_id> <PDB1> <PDB2> ...' in RANK ORDER, "
                         "which make_template_json.py consumes as --templates")
    a = ap.parse_args()
    if a.n_templates > N_TEMPLATES:
        sys.exit(f"FATAL: Mol-1 keeps only the first {N_TEMPLATES} templates "
                 f"(template_featurizer.py:774); asking for {a.n_templates} would silently "
                 f"drop the tail")

    # entry metadata: resolution + chain list, for the comp -> entry resolution step
    res: dict[str, float] = {}
    nchains: dict[str, int] = {}
    for r in rows(ENTRIES):
        try:
            res[r["pdb_id"].upper()] = float(r["resolution_A"])
        except (ValueError, KeyError):
            continue
        try:
            nchains[r["pdb_id"].upper()] = int(r["n_chains"])
        except (ValueError, KeyError):
            nchains[r["pdb_id"].upper()] = 0

    # comp -> the entries holding it, and its SMILES, for the manifest
    comp_entries: dict[str, list[str]] = {}
    comp_name: dict[str, str] = {}
    for r in rows(COMPS):
        cid = r["comp_id"]
        comp_entries[cid] = [e.strip().upper() for e in (r.get("pdb_entries") or "").split(";")
                             if e.strip()]
        comp_name[cid] = (r.get("chem_name") or "").strip()

    blinded = [r["Molecule_Name"] for r in rows(BLINDED)]
    if len(blinded) != 20:
        sys.exit(f"FATAL: expected 20 blinded ligands, found {len(blinded)}")

    def resolve(cid: str) -> tuple[str, float] | None:
        """Best-resolution entry holding this comp that also has a local CIF."""
        cands = [e for e in comp_entries.get(cid, []) if (CIF_DIR / f"{e}.cif").is_file()]
        if not cands:
            return None
        cands.sort(key=lambda e: (res.get(e, 99.0), e))
        return cands[0], res.get(cands[0], float("nan"))

    # group the matrix by blinded ligand, candidates only
    by_test: dict[str, list[dict]] = defaultdict(list)
    kept = dropped = 0
    for r in rows(SIM):
        if r.get("classification") != "PUTATIVE_LIGAND" or not truthy(r.get("in_pocket_any", "")):
            dropped += 1
            continue
        kept += 1
        by_test[r["test_id"]].append(r)
    print(f"{SIM.name}: {kept} candidate pairs kept, {dropped} dropped "
          f"(not PUTATIVE_LIGAND or not in pocket)")

    missing = [l for l in blinded if l not in by_test]
    if missing:
        sys.exit(f"FATAL: {len(missing)} blinded ligand(s) absent from the matrix: {missing}")

    man: list[dict] = []
    lines: list[str] = []
    best_sims: list[float] = []
    for lig in blinded:
        cand = by_test[lig]
        # resolve first, so the tie-break can use the resolved entry's resolution
        enriched = []
        for r in cand:
            got = resolve(r["pdb_comp_id"])
            if got is None:
                continue
            pdb, rr = got
            enriched.append((float(r["ecfp4_tanimoto"]), pdb, rr, r["pdb_comp_id"], r))
        # rank: similarity DESC, then resolution ASC, then comp_id ASC -- fully deterministic
        enriched.sort(key=lambda t: (-t[0], t[2], t[3]))
        picked, seen = [], set()
        for sim, pdb, rr, cid, r in enriched:
            if cid in seen:
                continue
            seen.add(cid)
            picked.append((sim, pdb, rr, cid, r))
            if len(picked) == a.n_templates:
                break
        if len(picked) < a.n_templates:
            sys.exit(f"FATAL: only {len(picked)} candidate(s) for {lig}, need {a.n_templates}")
        best_sims.append(picked[0][0])
        lines.append(f"{lig} " + " ".join(p[1] for p in picked))
        for rank, (sim, pdb, rr, cid, r) in enumerate(picked, 1):
            man.append({"base_id": lig, "rank": rank, "template_pdb": pdb,
                        "template_chain": CHAIN, "template_ligand": cid,
                        "template_ligand_name": comp_name.get(cid, ""),
                        "ecfp4_to_test": f"{sim:.4f}",
                        "maccs_to_test": r.get("maccs_tanimoto", ""),
                        "scaffold_identical": r.get("scaffold_identical", ""),
                        "resolution_A": f"{rr:.2f}", "n_chains": nchains.get(pdb, 0),
                        "all_entries_with_comp": ";".join(comp_entries.get(cid, []))})

    a.out_manifest.parent.mkdir(parents=True, exist_ok=True)
    with a.out_manifest.open("w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(man[0]))
        wr.writeheader()
        wr.writerows(man)
    a.out_list.write_text("".join(f"{ln}\n" for ln in lines))

    print(f"\n{len(man)} rows -> {a.out_manifest}")
    print(f"{len(lines)} bundles  -> {a.out_list}")
    pdbs = {m["template_pdb"] for m in man}
    print(f"\n{len(pdbs)} distinct PDB entries across the 20 bundles")
    print(f"multi-chain entries used: "
          f"{sorted({m['template_pdb'] for m in man if int(m['n_chains']) > 1}) or 'none'}")
    bs = sorted(best_sims)
    print(f"\n⚠ BEST-NEIGHBOUR ECFP4 CEILING: min {bs[0]:.3f}  median "
          f"{0.5 * (bs[9] + bs[10]):.3f}  max {bs[-1]:.3f}")
    print("  (this is the experiment's limitation and it belongs to the corpus, not the method)")
    print(f"\n{'ligand':<14}{'rank-1..4 ligand=ECFP4[PDB]'}")
    for lig in blinded:
        got = [m for m in man if m["base_id"] == lig]
        print(f"{lig:<14}" + "  ".join(
            f"{m['template_ligand']}={m['ecfp4_to_test']}[{m['template_pdb']}]" for m in got))
    return 0


if __name__ == "__main__":
    sys.exit(main())
