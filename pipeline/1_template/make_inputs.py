"""Mol-1 input JSONs for the 20 blinded CYP3A4 ligands: one job file per ligand.

Mol-1 is a Protenix fork, so the input is the familiar list-of-jobs JSON: each job's `sequences`
holds `proteinChain` / `ligand` / `ion` entries. One file per ligand rather than one 20-job file,
because generation runs as a Slurm array with one ligand per task.

Three things here are deliberate and easy to get wrong:

* The protein sequence and the heme are LIFTED from the campaign's existing OF3 query JSON
  (data/cyp3a4_relax/of3_inputs/queries_26_*.json) rather than retyped, so this arm folds the same
  503-aa P08684 construct and the same HEM cofactor as every other pool it will be ranked against.
* `unpairedMsaPath` is an ABSOLUTE path and is checked to exist. docs/infer_json_format.md: "A path
  that does not exist counts as missing: `mol pred` warns and searches instead" -- so a typo here
  does not fail, it silently fires an MSA search at the ColabFold service from a compute node with
  no IPv6. A missing a3m is therefore fatal in this script.
* Ligand SMILES go in verbatim from the blinded CSV. They carry the correct stereochemistry
  (checked 20/20 in this campaign), and re-canonicalising through RDKit is how that gets lost, so
  nothing parses them here.
* `covalent_bonds` (--fe-bond-table) is the ONLY channel Mol-1 has for "this ligand atom sits on
  the heme iron": the v1 constraint embedders were removed at the fork
  (mol/model/modules/embedders.py:20) and the Mol-1-constraint checkpoint is webserver-only, so
  SITE/EXCLUSION/CONTACT steering is not available locally. A declared bond sets ONE binary
  token_bonds bit feeding a learned pair bias (mol/model/mol.py:266-275). It imposes no bond
  length -- no reference conformer spans the two entities, ref_space_uid is per chain+residue and
  is computed after bonding -- so it cannot crush the 2.0-2.6 A coordination distance. It is
  however off-distribution: metal coordination bonds were deliberately stripped from Mol-1's
  training mmCIFs (mol/data/core/parser.py:100-101, rewrite_biotite.py:41) and the token_bonds
  prior was learned only from bonds under 2.4 A (mol/data/core/featurizer.py:481). Whether it
  tightens or distorts the coordination is therefore an empirical question, not a code one.

MSA DEPTH IS NOT SET HERE. It is a runner flag, so the same JSONs serve every subsampling arm:
  --data.msa.min_size.test D --data.msa.sample_cutoff.test D
pins the depth at exactly D random rows (mol/model/utils.py:296-317: sample_size is drawn from
U[lower_bound, n] and then truncated to cutoff, and the rows come from torch.randperm). Left alone,
the defaults are min 1 / cutoff 16384, which means the depth is RE-DRAWN EVERY FORWARD PASS -- so
"no flags" is a stochastic-depth arm, not a full-MSA one. Full depth needs min=cutoff=16384.

Env: any python3 with the stdlib (no deps). CPU, instant.
"""

from __future__ import annotations

import argparse
import os
import csv
import json
import sys
from pathlib import Path

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
BLINDED = R / "cyp-challenge-structure/cyp-challenge-TEST-BLINDED_structures.csv"
REF_QUERY = R / "data/cyp3a4_relax/of3_inputs/queries_26_dimer.json"
MSA = R / "data/cyp_structure_track/of3_msa/CYP3A4_P08684_full/colabfold_main.a3m"


def cyp3a4_sequence() -> str:
    """The 503-aa P08684 construct, taken from the arm this one will be compared against."""
    q = json.loads(REF_QUERY.read_text())["queries"]
    seqs = {c["sequence"] for job in q.values() for c in job["chains"]
            if c.get("molecule_type") == "protein"}
    if len(seqs) != 1:
        sys.exit(f"FATAL: expected one protein sequence in {REF_QUERY}, found {len(seqs)}")
    seq = seqs.pop()
    if len(seq) != 503:
        sys.exit(f"FATAL: expected a 503-aa CYP3A4 construct, got {len(seq)}")
    return seq


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--out-dir", type=Path, default=R / "data/cyp3a4_mol1/inputs")
    ap.add_argument("--msa", type=Path, default=MSA)
    ap.add_argument("--template-json", type=Path,
                    help="a Mol-1 JSON template file (see make_template_json.py). Only the .json "
                         "route is usable for the cryo-EM entries: .a3m/.hhr hits are filtered to "
                         "structures released on or before 2021-09-30 and would be dropped.")
    ap.add_argument("--templates-dir", type=Path,
                    help="PER-LIGAND templates: <dir>/<base_id><--templates-suffix>. Every one "
                         "of the 20 must exist or the build fails -- a missing file is treated "
                         "as ABSENT, which would leave one quietly template-free system inside "
                         "a templated arm. Mutually exclusive with --template-json.")
    ap.add_argument("--templates-suffix", default=".json",
                    help="filename suffix after the base_id in --templates-dir")
    ap.add_argument("--fe-bond-table", type=Path,
                    help="CSV with base_id,donor_atom columns (see fe_donor_table.py). Each listed "
                         "ligand gets a declared covalent bond from HEM's FE to that atom; ligands "
                         "absent from the table are written unbonded, so one table drives a mixed "
                         "arm (the 11 Type-II ligands bonded, the other 9 left alone).")
    a = ap.parse_args()
    if not a.msa.is_file():
        sys.exit(f"FATAL: no MSA at {a.msa} -- mol pred would silently run its own search instead")
    if a.template_json and a.templates_dir:
        sys.exit("FATAL: --template-json and --templates-dir are mutually exclusive: one bundle "
                 "for every ligand, or one per ligand, not both")
    if a.template_json and not a.template_json.is_file():
        # Same trap as the MSA: a templatesPath that does not exist is treated as absent, and the
        # run then SEARCHES for templates instead of failing.
        sys.exit(f"FATAL: no template JSON at {a.template_json}")
    if a.templates_dir and not a.templates_dir.is_dir():
        sys.exit(f"FATAL: --templates-dir is not a directory: {a.templates_dir}")
    fe_bonds: dict[str, str] = {}
    if a.fe_bond_table:
        if not a.fe_bond_table.is_file():
            sys.exit(f"FATAL: no Fe bond table at {a.fe_bond_table}")
        with a.fe_bond_table.open(newline="") as fh:
            fe_bonds = {r["base_id"]: r["donor_atom"] for r in csv.DictReader(fh)}
        if not fe_bonds:
            sys.exit(f"FATAL: {a.fe_bond_table} lists no ligands")
    seq = cyp3a4_sequence()
    a.out_dir.mkdir(parents=True, exist_ok=True)

    with BLINDED.open(newline="") as fh:
        ligands = [(r["Molecule_Name"], r["SMILES"]) for r in csv.DictReader(fh)]
    if len(ligands) != 20:
        sys.exit(f"FATAL: expected 20 blinded ligands, found {len(ligands)}")

    # Resolve the per-ligand bundle up front, so a missing one is a build error rather than a
    # silently template-free system -- the exact failure the .json route was chosen to avoid.
    tmpl_for: dict[str, Path] = {}
    if a.templates_dir:
        missing = []
        for name, _ in ligands:
            cand = (a.templates_dir / f"{name}{a.templates_suffix}").resolve()
            if not cand.is_file():
                missing.append(cand.name)
            tmpl_for[name] = cand
        if missing:
            sys.exit(f"FATAL: --templates-dir is missing {len(missing)} bundle(s): "
                     f"{', '.join(missing)}")
        heads = len({p.read_bytes()[:4096] for p in tmpl_for.values()})
        print(f"{len(tmpl_for)} per-ligand template bundles resolved, {heads} distinct by head")

    written = []
    for name, smiles in ligands:
        chain = {"sequence": seq, "count": 1, "id": ["A"], "unpairedMsaPath": str(a.msa)}
        if tmpl_for:
            chain["templatesPath"] = str(tmpl_for[name])
        elif a.template_json:
            chain["templatesPath"] = str(a.template_json)
        job = {
            "name": name,
            "sequences": [
                {"proteinChain": chain},
                {"ligand": {"ligand": "CCD_HEM", "count": 1, "id": ["B"]}},
                {"ligand": {"ligand": smiles, "count": 1, "id": ["C"]}},
            ],
        }
        if name in fe_bonds:
            # The entity numbers are 1-based positions in `sequences` above, so proteinChain is 1,
            # CCD_HEM is 2 and the SMILES ligand is 3; `position` is 1 for every single-CCD and
            # SMILES ligand. atom2 uses Mol-1's own SMILES atom naming (element + 1-based
            # per-element occurrence, json_parser.py:456-461), which is exactly what the donor map
            # stores. A name matching no atom is a hard assert at json_to_feature.py:268-270, so a
            # typo kills the task instead of quietly folding an unbonded complex.
            job["covalent_bonds"] = [{
                "entity1": "2", "copy1": 1, "position1": "1", "atom1": "FE",
                "entity2": "3", "copy2": 1, "position2": "1", "atom2": fe_bonds[name],
            }]
        p = a.out_dir / f"{name}.json"
        p.write_text(json.dumps([job], indent=2))      # a LIST of jobs, even for one
        written.append(p)

    # ABSOLUTE paths: the generation sbatch cds into the Mol-1 repo before calling the runner, so a
    # relative entry here resolves against the wrong directory and the task dies with
    # FileNotFoundError after the GPU is already allocated (job 501431_12).
    index = a.out_dir / "index.txt"
    index.write_text("".join(f"{p.resolve()}\n" for p in written))
    print(f"{len(written)} input JSONs -> {a.out_dir}")
    print(f"  protein {len(seq)} aa + CCD_HEM + one SMILES ligand each")
    print(f"  unpairedMsaPath {a.msa} ({sum(1 for ln in a.msa.open() if ln[0] == '>')} sequences)")
    print(f"  templatesPath {a.templates_dir and f'PER-LIGAND from {a.templates_dir}' or a.template_json or 'none (template-free)'}")
    names = {n for n, _ in ligands}
    bonded = [n for n, _ in ligands if n in fe_bonds]
    if stray := sorted(set(fe_bonds) - names):
        sys.exit(f"FATAL: Fe bond table names ligands that are not in the blinded 20: {stray}")
    print(f"  covalent_bonds HEM:FE -> {len(bonded)} of {len(ligands)} ligands"
          + (": " + ", ".join(f"{n}:{fe_bonds[n]}" for n in bonded) if bonded else " (none)"))
    print(f"  array index -> {index}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
