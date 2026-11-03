#!/usr/bin/env python3
"""Every available metric, picked-pose for picked-pose, between two submission zips.

Compares the poses the two zips ACTUALLY SHIPPED (from each .provenance.csv), joined back to
its scoring table for the MM-GBSA decomposition. The decomposition is the point: delta_total
alone cannot say whether a difference is shape complementarity (vdw) or electrostatics (eel),
and the Fe-coordination loss measured for this pair predicts an eel effect on Type-II
specifically.

Metrics: delta_total and its four components, srs (the model's own ranking score), worst
covalent bond length, Fe-ligand distance. Also the lDDT/L-RMSD agreement, when its tables
exist -- that part is AGREEMENT, NOT ACCURACY, since the blinded 20 have no reference.

Sign tests are exact and two-sided over the paired ligands; ties are dropped, the conservative
convention. Every CSV is read with csv.DictReader, never a split on "," -- the lDDT tables'
assignment columns hold python-repr tuples containing commas.

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import argparse
import os
import csv
import math
import sys
from pathlib import Path

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
S = R / "data/cyp3a4_mol1/submission"
G = R / "data/cyp3a4_ob0_mmgbsa"
D = R / "data/cyp3a4_mol1"
T2_FILE = R / "data/cyp3a4_xmodel/typeii_base_ids.txt"
FE_COORD_PRIOR = 2.58
# Lower is better for every energy term and for max_cc; Fe is reported, not judged.
LOWER_BETTER = {"delta_total", "vdw", "eel", "egb", "esurf", "max_cc_bond_A"}


def read(p: Path, key: str = "base_id") -> dict[str, dict]:
    if not p.is_file():
        sys.exit(f"FATAL: missing {p}")
    with open(p, newline="") as fh:
        return {r[key]: r for r in csv.DictReader(fh)}


def read_by_pose(p: Path) -> dict[str, dict]:
    return read(p, "pose_id")


def num(r: dict | None, k: str) -> float:
    if not r:
        return float("nan")
    v = r.get(k, "")
    try:
        return float("nan") if v in ("", "nan", None) else float(v)
    except ValueError:
        return float("nan")


def sign_test(d: list[float]) -> tuple[int, int, int, float]:
    pos = sum(x > 0 for x in d); neg = sum(x < 0 for x in d)
    ties = sum(x == 0 for x in d); n = pos + neg
    if n == 0:
        return pos, neg, ties, 1.0
    k = min(pos, neg)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return pos, neg, ties, min(1.0, 2 * tail)


def med(v: list[float]) -> float:
    s = sorted(v); n = len(s)
    return float("nan") if not n else (s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2]))


def fmt(x: float, nd: int = 2) -> str:
    return "  -  " if not math.isfinite(x) else f"{x:.{nd}f}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a-name", default="submission_mol1_tfg_mmgbsa")
    ap.add_argument("--b-name", default="submission_mol1_tmpl_mmgbsa")
    ap.add_argument("--a-scored", type=Path, default=G / "scored_mol1_tfg.csv")
    ap.add_argument("--b-scored", type=Path, default=D / "scored_mol1_tmpl_full.csv")
    ap.add_argument("--label-a", default="tfg")
    ap.add_argument("--label-b", default="tmpl")
    a = ap.parse_args()

    t2 = {l.strip() for l in open(T2_FILE) if l.strip()}
    PA, PB = read(S / f"{a.a_name}.provenance.csv"), read(S / f"{a.b_name}.provenance.csv")
    SA, SB = read_by_pose(a.a_scored), read_by_pose(a.b_scored)
    ligs = sorted(set(PA) & set(PB))
    if not ligs:
        sys.exit("FATAL: the two zips share no ligand")
    print(f"A = {a.label_a}  ({a.a_name})")
    print(f"B = {a.label_b}  ({a.b_name})")
    print(f"{len(ligs)} ligands, {len(t2 & set(ligs))} Type-II\n")

    # Join each shipped pose to its scoring row. A miss here means the provenance and the
    # scoring table disagree, which would silently void the comparison -- so it is loud.
    rowA, rowB, miss = {}, {}, []
    for l in ligs:
        ra, rb = SA.get(PA[l]["pose_id"]), SB.get(PB[l]["pose_id"])
        if ra is None or rb is None:
            miss.append((l, PA[l]["pose_id"] if ra is None else PB[l]["pose_id"]))
        rowA[l], rowB[l] = ra, rb
    if miss:
        print(f"⚠ {len(miss)} shipped pose(s) absent from their scoring table "
              f"(components will show '-'): {miss[:3]}\n")

    # ---- MM-GBSA decomposition
    print("=== MM-GBSA decomposition of the SHIPPED poses (kcal/mol, lower = tighter) ===")
    print(f"{'ligand':<14}{'cls':<6}"
          f"{'dTOT A|B':>16}{'vdw A|B':>16}{'eel A|B':>16}{'egb A|B':>16}{'esurf A|B':>14}")
    terms = ("delta_total", "vdw", "eel", "egb", "esurf")
    deltas: dict[str, list[float]] = {t: [] for t in terms}
    d_t2: dict[str, list[float]] = {t: [] for t in terms}
    d_non: dict[str, list[float]] = {t: [] for t in terms}
    for l in ligs:
        cls = "T-II" if l in t2 else "non"
        cells = []
        for t in terms:
            va, vb = num(rowA[l], t), num(rowB[l], t)
            cells.append(f"{fmt(va):>7}|{fmt(vb):<8}" if t != "esurf"
                         else f"{fmt(va):>6}|{fmt(vb):<7}")
            d = va - vb        # A minus B; for an energy, negative = A is tighter
            deltas[t].append(d)
            (d_t2 if l in t2 else d_non)[t].append(d)
        print(f"{l:<14}{cls:<6}" + "".join(cells))

    print(f"\n=== per-term paired test (A - B; negative = {a.label_a} tighter) ===")
    print(f"{'term':<14}{'n':>4}{'mean A':>10}{'mean B':>10}{'median d':>11}"
          f"{'favours A':>11}{'favours B':>11}{'p':>8}")
    for t in terms:
        va = [num(rowA[l], t) for l in ligs]
        vb = [num(rowB[l], t) for l in ligs]
        va = [x for x in va if math.isfinite(x)]; vb = [x for x in vb if math.isfinite(x)]
        d = [x for x in deltas[t] if math.isfinite(x)]
        # favours A = A is LOWER (tighter) for these energy terms
        p, n, ties, pv = sign_test([-x for x in d])
        print(f"{t:<14}{len(d):>4}{sum(va) / len(va):>10.2f}{sum(vb) / len(vb):>10.2f}"
              f"{med(d):>+11.2f}{p:>11}{n:>11}{pv:>8.3f}")

    print(f"\n--- class split on each term (median A - B) ---")
    print(f"{'term':<14}{'Type-II':>12}{'non-Type-II':>14}{'T-II favours A':>17}{'p (T-II)':>10}")
    for t in terms:
        dt = [x for x in d_t2[t] if math.isfinite(x)]
        dn = [x for x in d_non[t] if math.isfinite(x)]
        p, n, ties, pv = sign_test([-x for x in dt])
        print(f"{t:<14}{med(dt):>+12.2f}{med(dn):>+14.2f}{f'{p}/{p + n}':>17}{pv:>10.3f}")

    # ---- non-energy metrics
    print(f"\n=== other shipped-pose metrics ===")
    print(f"{'ligand':<14}{'cls':<6}{'srs A|B':>16}{'maxCC A|B':>18}{'Fe A|B':>18}{'Fe elem A|B':>14}")
    for l in ligs:
        pa, pb = PA[l], PB[l]
        print(f"{l:<14}{'T-II' if l in t2 else 'non':<6}"
              f"{fmt(num(pa, 'srs'), 3):>7}|{fmt(num(pb, 'srs'), 3):<8}"
              f"{fmt(num(pa, 'max_cc_bond_A'), 3):>8}|{fmt(num(pb, 'max_cc_bond_A'), 3):<9}"
              f"{fmt(num(pa, 'fe_lig_min_A'), 3):>8}|{fmt(num(pb, 'fe_lig_min_A'), 3):<9}"
              f"{pa.get('fe_lig_element', '-'):>5}|{pb.get('fe_lig_element', '-'):<6}")
    for k, nd, lo in (("srs", 3, False), ("max_cc_bond_A", 3, True), ("fe_lig_min_A", 3, None)):
        va = [num(PA[l], k) for l in ligs]; vb = [num(PB[l], k) for l in ligs]
        d = [x - y for x, y in zip(va, vb) if math.isfinite(x) and math.isfinite(y)]
        orient = [-x for x in d] if lo else d        # favours A = higher, unless lower is better
        p, n, ties, pv = sign_test(orient)
        tag = "" if lo is None else ("  (lower better)" if lo else "  (higher better)")
        print(f"  {k:<16} mean A {sum(va) / len(va):>8.3f}   mean B {sum(vb) / len(vb):>8.3f}   "
              f"median d {med(d):>+8.3f}   A/B {p}/{n}  p={pv:.3f}{tag}")
    for name, P in ((a.label_a, PA), (a.label_b, PB)):
        fe = [num(P[l], "fe_lig_min_A") for l in ligs if l in t2]
        fe = [x for x in fe if math.isfinite(x)]
        if fe:
            inside = sum(1 for x in fe if x <= FE_COORD_PRIOR)
            print(f"  Type-II inside the {FE_COORD_PRIOR} A coordinating prior, {name}: "
                  f"{inside}/{len(fe)}   median {med(fe):.3f} A")

    # ---- agreement, if the lDDT tables are there
    fwd, rev = S / f"lddt_{a.label_a}_vs_{a.label_b}.csv", S / f"lddt_{a.label_b}_vs_{a.label_a}.csv"
    if fwd.is_file() and rev.is_file():
        F, V = read(fwd), read(rev)
        print(f"\n=== agreement (⚠ NOT accuracy -- both sides are predictions) ===")
        print(f"{'ligand':<14}{'cls':<6}{'PLI A|B':>9}{'PLI B|A':>9}{'LP A|B':>9}"
              f"{'LP B|A':>9}{'L-RMSD':>9}{'COM':>8}")
        acc: dict[str, list] = {}
        for l in ligs:
            r = (num(F.get(l), "portal_lddt_pli"), num(V.get(l), "portal_lddt_pli"),
                 num(F.get(l), "portal_lddt_lp"), num(V.get(l), "portal_lddt_lp"),
                 num(F.get(l), "portal_bisyrmsd"), num(F.get(l), "portal_com_shift"))
            cls = "T-II" if l in t2 else "non"
            print(f"{l:<14}{cls:<6}{r[0]:>9.3f}{r[1]:>9.3f}{r[2]:>9.3f}{r[3]:>9.3f}"
                  f"{r[4]:>9.3f}{r[5]:>8.3f}")
            for g in ("ALL", cls):
                acc.setdefault(g, []).append(r)
        print(f"\n{'group':<14}{'n':>4}{'PLI A|B':>9}{'PLI B|A':>9}{'LP A|B':>9}"
              f"{'LP B|A':>9}{'L-RMSD':>9}{'COM':>8}")
        for g in ("ALL", "T-II", "non"):
            rs = acc.get(g, [])
            if rs:
                m = [sum(x[i] for x in rs) / len(rs) for i in range(6)]
                print(f"{'mean ' + g:<14}{len(rs):>4}{m[0]:>9.3f}{m[1]:>9.3f}{m[2]:>9.3f}"
                      f"{m[3]:>9.3f}{m[4]:>9.3f}{m[5]:>8.3f}")
        al = acc.get("ALL", [])
        print(f"\n{sum(1 for x in al if x[4] < 1e-3)}/{len(al)} identical (L-RMSD < 0.001 A); "
              f"{sum(1 for x in al if x[4] >= 2.0)}/{len(al)} differ by >= 2 A")
    else:
        print(f"\n(no lDDT tables yet at {fwd.name} / {rev.name} -- agreement section skipped)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
