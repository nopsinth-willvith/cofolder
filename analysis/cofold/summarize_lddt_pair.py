#!/usr/bin/env python3
"""Per-ligand + class-split lDDT-PLI / lDDT-LP / L-RMSD between two submission zips.

⚠ AGREEMENT, NOT ACCURACY. Both sides are predictions, so a low score says the two disagree,
not that either is wrong. In particular: if one zip has a PORTAL score, that score does NOT
transfer through this table -- agreeing with a 0.82-LP zip does not give you 0.82.

Reads the scored CSVs with csv.reader, never a split on ",": `scrmsd_assignment` and
`lddtpli_assignment` hold python-repr tuples containing commas, which shifts every later field
by two and silently swaps lDDT-PLI for lDDT-LP.

Unlike report_ph48_vs_ref.py this needs no .provenance.csv, so it works against a reference
built by a different stage.

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import argparse
import os
import csv
import sys
from pathlib import Path

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
T2_FILE = R / "data/cyp3a4_xmodel/typeii_base_ids.txt"


def rows(p: Path) -> dict[str, dict]:
    if not p.is_file():
        sys.exit(f"FATAL: no scored table at {p}")
    with open(p, newline="") as fh:
        return {r["base_id"]: r for r in csv.DictReader(fh)}


def num(r: dict, k: str) -> float:
    v = r.get(k, "")
    return float("nan") if v in ("", "nan", None) else float(v)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fwd", type=Path, required=True, help="A-vs-B scored csv")
    ap.add_argument("--rev", type=Path, required=True, help="B-vs-A scored csv")
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--t2-file", type=Path, default=T2_FILE)
    a = ap.parse_args()

    if not a.t2_file.is_file():
        sys.exit(f"FATAL: no Type-II id list at {a.t2_file}")
    t2 = {l.strip() for l in open(a.t2_file) if l.strip()}
    if not t2:
        sys.exit(f"FATAL: {a.t2_file} is empty")
    F, V = rows(a.fwd), rows(a.rev)
    shared = sorted(set(F) & set(V))
    if not shared:
        sys.exit("FATAL: the two directions share no ligand")
    only = (set(F) | set(V)) - set(shared)
    if only:
        print(f"⚠ {len(only)} ligand(s) in only one direction, excluded: {sorted(only)[:5]}")
    print(f"A = {a.label_a}\nB = {a.label_b}   (reference)\n{len(shared)} ligands, "
          f"{len(t2 & set(shared))} Type-II\n")

    print(f"{'ligand':<14}{'cls':<6}{'PLI A|B':>9}{'PLI B|A':>9}{'LP A|B':>9}{'LP B|A':>9}"
          f"{'RMSD A|B':>10}{'RMSD B|A':>10}{'COM':>8}")
    acc: dict[str, list] = {}
    for lig in shared:
        f, v = F[lig], V[lig]
        row = (num(f, "portal_lddt_pli"), num(v, "portal_lddt_pli"),
               num(f, "portal_lddt_lp"), num(v, "portal_lddt_lp"),
               num(f, "portal_bisyrmsd"), num(v, "portal_bisyrmsd"),
               num(f, "portal_com_shift"))
        cls = "T-II" if lig in t2 else "non"
        print(f"{lig:<14}{cls:<6}{row[0]:>9.3f}{row[1]:>9.3f}{row[2]:>9.3f}{row[3]:>9.3f}"
              f"{row[4]:>10.3f}{row[5]:>10.3f}{row[6]:>8.3f}")
        for g in ("ALL", cls):
            acc.setdefault(g, []).append(row)

    print(f"\n{'group':<14}{'n':>4}{'PLI A|B':>9}{'PLI B|A':>9}{'LP A|B':>9}{'LP B|A':>9}"
          f"{'RMSD A|B':>10}{'RMSD B|A':>10}{'COM':>8}")
    for g in ("ALL", "T-II", "non"):
        rs = acc.get(g, [])
        if not rs:
            continue
        m = [sum(r[i] for r in rs) / len(rs) for i in range(7)]
        print(f"{'mean ' + g:<14}{len(rs):>4}{m[0]:>9.3f}{m[1]:>9.3f}{m[2]:>9.3f}{m[3]:>9.3f}"
              f"{m[4]:>10.3f}{m[5]:>10.3f}{m[6]:>8.3f}")
    al = acc.get("ALL", [])
    print(f"\n{sum(1 for r in al if r[4] < 1e-3)}/{len(al)} ligands identical "
          f"(L-RMSD < 0.001 A); {sum(1 for r in al if r[4] >= 2.0)}/{len(al)} differ by >= 2 A")
    print("\n⚠ agreement, not accuracy: a PORTAL score on the reference does NOT transfer "
          "through this table.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
