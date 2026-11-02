#!/usr/bin/env python3
"""Paired per-ligand comparison of the TFG pool against the unguided templated pools.

Two questions the submission log raises but does not answer:
  1. Does TFG win on the ranking column itself (delta_total), or only on PoseBusters?
  2. The Fe ladder relaxed for all 11 Type-II ligands, which no previous build did. Are the
     TFG poses systematically FURTHER from the heme Fe than the unguided arms?

Pairing is per ligand on the pool's BEST (most negative) delta_total, because that is the
quantity the selector actually acts on.
"""
import csv, math, sys
from collections import defaultdict
from pathlib import Path

G = Path("/vast.mnt/project/bme-20255165-cofolding-physics/data/cyp3a4_ob0_mmgbsa")
D = Path("/vast.mnt/project/bme-20255165-cofolding-physics/data/cyp3a4_mol1")
# The coordinating / non-coordinating Fe-ligand prior from the 122-structure CYP3A4 curation.
FE_COORD_PRIOR = 2.58


def sign_test(deltas):
    """Two-sided exact sign test. Ties dropped, the conservative convention."""
    pos = sum(d > 0 for d in deltas); neg = sum(d < 0 for d in deltas)
    ties = sum(d == 0 for d in deltas); n = pos + neg
    if n == 0:
        return pos, neg, ties, 1.0
    k = min(pos, neg)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return pos, neg, ties, min(1.0, 2 * tail)


def med(v):
    s = sorted(v)
    n = len(s)
    return float("nan") if not n else (s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2]))


def load(path, want_pools=None):
    """base_id -> list of (delta_total, fe_lig_min_A, pool). Only status-ok, scored rows."""
    out = defaultdict(list)
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            if want_pools is not None and r.get("pool") not in want_pools:
                continue
            try:
                dt = float(r["delta_total"])
            except (KeyError, ValueError, TypeError):
                continue
            if not math.isfinite(dt):
                continue
            try:
                fe = float(r.get("fe_lig_in_A") or r.get("fe_lig_min_A") or "nan")
            except ValueError:
                fe = float("nan")
            out[r["base_id"]].append((dt, fe, r.get("pool", "?")))
    return out


def main():
    tfg = load(G / "scored_mol1_tfg.csv")
    ref_pools = {"mol1_tcryo5", "mol1_t38FX", "mol1_tcomb4", "mol1_full"}
    ref = load(D / "scored_mol1_tmpl_full.csv", ref_pools)
    if not tfg or not ref:
        print("FATAL: one of the scored tables is empty", file=sys.stderr)
        return 1

    # Type-II membership straight off the zip audit so it matches the submission exactly.
    t2 = set()
    ap = D / "submission" / "submission_mol1_tfg_mmgbsa.provenance.csv"
    if ap.exists():
        with open(ap, newline="") as fh:
            for r in csv.DictReader(fh):
                if str(r.get("typeii", "")).strip().lower() in ("true", "1", "yes"):
                    t2.add(r["base_id"])

    ligs = sorted(set(tfg) & set(ref))
    print(f"TFG pool: {sum(len(v) for v in tfg.values())} scored poses, {len(tfg)} ligands")
    print(f"ref pools {sorted(ref_pools)}: {sum(len(v) for v in ref.values())} scored poses, {len(ref)} ligands")
    print(f"paired on {len(ligs)} ligands; Type-II = {len(t2)}\n")

    hdr = f"{'ligand':<14}{'cls':<7}{'best dT  TFG|ref':<24}{'delta':<10}{'bestFe TFG|ref':<20}{'<2.6A TFG|ref'}"
    print(hdr)
    d_dt, d_fe = [], []
    for l in ligs:
        bt = min(x[0] for x in tfg[l]); br = min(x[0] for x in ref[l])
        fts = [x[1] for x in tfg[l] if math.isfinite(x[1])]
        frs = [x[1] for x in ref[l] if math.isfinite(x[1])]
        ft = min(fts) if fts else float("nan"); fr = min(frs) if frs else float("nan")
        nt = sum(1 for f in fts if f <= FE_COORD_PRIOR); nr = sum(1 for f in frs if f <= FE_COORD_PRIOR)
        # delta > 0 means TFG is BETTER (more negative delta_total = tighter binding).
        d_dt.append(br - bt)
        if math.isfinite(ft) and math.isfinite(fr):
            d_fe.append(fr - ft)   # > 0 means TFG sits CLOSER to Fe
        cls = "T-II" if l in t2 else "non"
        print(f"{l:<14}{cls:<7}{bt:>9.2f}|{br:<13.2f}{br - bt:>+8.2f}  "
              f"{ft:>8.3f}|{fr:<10.3f}{nt:>5d}|{nr}")

    print(f"\n=== ranking column (delta_total, lower is better) ===")
    p, n, t, pv = sign_test(d_dt)
    print(f"  best-pose dT   {p} ligands favour TFG, {n} favour ref, {t} tied   p={pv:.3f}   "
          f"median delta {med(d_dt):+.2f} kcal/mol")

    print(f"\n=== Fe-ligand distance (coord prior {FE_COORD_PRIOR} A) ===")
    p, n, t, pv = sign_test(d_fe)
    print(f"  closest approach   {p} ligands favour TFG, {n} favour ref, {t} tied   p={pv:.3f}   "
          f"median delta {med(d_fe):+.3f} A  (positive = TFG closer to Fe)")
    for name, tbl in (("TFG", tfg), ("ref", ref)):
        t2fe = [x[1] for l in ligs if l in t2 for x in tbl[l] if math.isfinite(x[1])]
        if t2fe:
            print(f"  {name:<4} Type-II poses: n={len(t2fe)}  median {med(t2fe):.3f} A  "
                  f"min {min(t2fe):.3f}  <=2.6A {sum(1 for f in t2fe if f <= FE_COORD_PRIOR)} "
                  f"({100 * sum(1 for f in t2fe if f <= FE_COORD_PRIOR) / len(t2fe):.1f}%)")
    print("\n'better' on delta_total is agreement with MM-GBSA, not accuracy: the blinded 20 have\n"
          "no reference, and the guided arm resamples as well as steers.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
