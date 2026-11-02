"""Build Mol-1 JSON templates from CYP3A4 experimental structures, for any template family.

Two families are in use. `--cif-dir` chooses between them:
  * the Orta et al. cryo-EM entries (38FU/V/W/X/Y), the default, under data/cyp3a4_cryoem_holdout/raw
  * the four X-ray chains that submission_combined actually consumed (5VCC/6MA8/6MA6/6BCZ chain A),
    under data/cyp_structure_track/cif -- see the note on template COUNT below.

Mol-1 takes templates three ways and only one of them is usable here. `.a3m`/`.hhr` hits are
filtered to structures released on or before 2021-09-30 and silently dropped when missing from the
local release-date table -- the Orta et al. cryo-EM entries (38FU apo / 38FV ritonavir /
38FW ketoconazole / 38FX azamulin / 38FY vardenafil) are far newer, so that route would quietly
template nothing. The `.json` route has no date filter and runs no alignment: it takes the template
mmCIF text plus explicit `queryIndices`/`templateIndices` pairs, which is why it is used here.

These are the same structures the OpenBind-0 arms used (`template_cif_paths` + chain A of
data/cyp3a4_cryoem_holdout/raw/38F*.cif), so a Mol-1 template arm is comparable with those pools.
Unlike the ftconf5000+38FX pool, this one varies ONLY the template -- same stock checkpoint -- so
"the cryo template helps" is not confounded with "the fine-tune helps".

The residue pairing is built from author numbering (CYP3A4 structures are numbered in UniProt
P08684 coordinates, so query index = auth_seq_id - 1) and then VERIFIED residue by residue against
the query sequence. A position whose template residue does not match the query is dropped and
counted, so a numbering offset shows up as a collapsed pair count rather than as a template quietly
teaching the model the wrong register.

Env: the Mol-1 venv (gemmi). CPU, seconds.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import gemmi

# The query sequence comes from make_inputs, the module that WRITES the inputs these templates are
# paired against, rather than being re-derived here from data/cyp3a4_mol1/inputs/. That directory
# was the only hardcoded path left in this script and it belongs to whichever arm ran last: stage 80
# writes its inputs to inputs_tcomb4 while this read still came from inputs/. Verified identical
# (503 aa, byte-for-byte) on job 505153 before the switch, so it is not a behaviour change.
from make_inputs import cyp3a4_sequence  # noqa: E402  (same directory; sys.path[0] is the script's)

R = Path("/vast.mnt/project/bme-20255165-cofolding-physics")
RAW = R / "data/cyp3a4_cryoem_holdout/raw"


def template_entry(cif: Path, chain_id: str, query: str) -> tuple[dict, int]:
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.remove_ligands_and_waters()
    chain = st[0].find_chain(chain_id)
    if chain is None:
        sys.exit(f"FATAL: {cif.name} has no chain {chain_id}")

    q_idx, t_idx, mismatched = [], [], 0
    # t_pos counts EVERY residue of the chain that gets written, not just the amino acids. Mol-1
    # builds its template residue array over the whole written chain (template_parser.py:402-422 ->
    # template_utils.py:920-924), so a non-tabulated or non-amino-acid residue surviving
    # remove_ligands_and_waters would shift every later templateIndex by one. The per-residue
    # identity check below CANNOT see that kind of shift -- it compares sequence identity, not
    # register -- so enumerating over the full chain is what actually guarantees it. None of the
    # nine chains in use has such a residue, which is why this is a latent fix and not a bug fix.
    for t_pos, res in enumerate(chain):
        info = gemmi.find_tabulated_residue(res.name)
        if info is None or not info.is_amino_acid():
            continue
        one = info.one_letter_code.upper()
        qi = res.seqid.num - 1                 # auth numbering is UniProt numbering here
        if 0 <= qi < len(query) and query[qi] == one:
            q_idx.append(qi)
            t_idx.append(t_pos)
        else:
            mismatched += 1

    # Only the paired chain is written out, so Mol-1's "first chain is used" rule picks the right
    # one no matter what else the deposited entry contains.
    out = gemmi.Structure()
    out.add_model(gemmi.Model("1"))
    out[0].add_chain(chain.clone())
    out.setup_entities()
    return ({"mmcif": out.make_mmcif_document().as_string(),
             "queryIndices": q_idx, "templateIndices": t_idx}, mismatched)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--templates", nargs="+", default=["38FU", "38FV", "38FW", "38FX", "38FY"])
    ap.add_argument("--chain", default="A")
    ap.add_argument("--cif-dir", type=Path, default=RAW,
                    help="directory holding <ID>.cif (default: the cryo-EM raw dir)")
    ap.add_argument("--out-dir", type=Path, default=R / "data/cyp3a4_mol1/templates")
    ap.add_argument("--bundle",
                    help="also write <BUNDLE>.json holding every template at once, in the order "
                         "given. `--bundle cryo5` reproduces the old --cryo5 flag. Mol-1 keeps the "
                         "FIRST 4 templates in a bundle (InferenceTemplateFeaturizer hard-codes "
                         "max_templates=4 and does NOT sort), so a bundle of more than four "
                         "silently truncates in the order written here.")
    a = ap.parse_args()

    query = cyp3a4_sequence()
    print(f"query {len(query)} aa; templates from {a.cif_dir}")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    entries = {}
    for t in a.templates:
        cif = a.cif_dir / f"{t}.cif"
        if not cif.is_file():
            sys.exit(f"FATAL: no template CIF at {cif}")
        entry, bad = template_entry(cif, a.chain, query)
        entries[t] = entry
        n = len(entry["queryIndices"])
        p = a.out_dir / f"{t}.json"
        p.write_text(json.dumps([entry]))
        cov = 100 * n / len(query)
        print(f"  {t}: {n} paired residues ({cov:.1f}% of the query), {bad} dropped on mismatch"
              f" -> {p.name}")
        if n < 300:
            print(f"    ⚠ only {n} pairs -- check the numbering register before trusting this arm")
    if a.bundle:
        if len(a.templates) > 4:
            print(f"  ⚠ {len(a.templates)} templates in one bundle but Mol-1 features only the "
                  f"first 4: {a.templates[:4]} used, {a.templates[4:]} DROPPED")
        p = a.out_dir / f"{a.bundle}.json"
        p.write_text(json.dumps([entries[t] for t in a.templates]))
        print(f"  {a.bundle}: all {len(a.templates)} templates in one file -> {p.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
