#!/usr/bin/env python3
"""Per-pair (protein-ligand) structural confidence from Mol-1's own heads, per templated pool.

Which column answers "confidence of this PAIR":
  cp_iptm_pl  -- chain-pair ipTM for the protein<->ligand chain pair. THIS is the per-pair
                 interface confidence and the one to quote.
  ranking_score (srs) -- the model's composite pose score. ⚠ its clash term is POLYMER-ONLY,
                 so it under-penalises a ligand clash; prefer cp_iptm_pl for the pair.
  iptm / ptm  -- whole-complex interface / global TM estimates, not pair-specific.
  plddt       -- per-atom confidence averaged over the WHOLE complex, so it is dominated by
                 the 503-aa protein and barely moves with the ligand.
  gpde        -- predicted distance error; LOWER is better, unlike every column above.

⚠ Confidence is the model's self-assessment, NOT accuracy. On this target the campaign has
measured `ranking_score` ranking BELOW CHANCE on the Type-II subclass, and ipTM behaving as
random against a known oracle. Report these numbers; do not rank on them.

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

R = Path("/vast.mnt/project/bme-20255165-cofolding-physics")
D = R / "data/cyp3a4_mol1"
T2_FILE = R / "data/cyp3a4_xmodel/typeii_base_ids.txt"
COLS = ("cp_iptm_pl", "ranking_score", "iptm", "ptm", "plddt", "gpde")


def med(v: list[float]) -> float:
    s = sorted(v); n = len(s)
    return float("nan") if not n else (s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2]))


def num(r: dict, k: str) -> float:
    v = r.get(k, "")
    try:
        return float("nan") if v in ("", "nan", None) else float(v)
    except ValueError:
        return float("nan")


def load(pool: str) -> dict[str, list[dict]]:
    p = D / f"poses_qc_{pool}.csv"
    if not p.is_file():
        return {}
    by: dict[str, list[dict]] = defaultdict(list)
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            by[r["base_id"]].append(r)
    return by


def f(x: float, nd: int = 3) -> str:
    return "  -  " if not math.isfinite(x) else f"{x:.{nd}f}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pools", nargs="+",
                    default=["mol1_tcryo5", "mol1_tcomb4", "mol1_t38FX", "mol1_tfg_tcryo5"])
    a = ap.parse_args()

    t2 = {l.strip() for l in open(T2_FILE) if l.strip()}
    pools = {p: load(p) for p in a.pools}
    live = [p for p, v in pools.items() if v]
    if not live:
        sys.exit("FATAL: no QC table found for any requested pool")
    for p in a.pools:
        if not pools[p]:
            print(f"⚠ no QC table for {p} -- skipped (still generating or not QC'd)")
    ligs = sorted({l for p in live for l in pools[p]})
    print(f"\npools: {', '.join(live)}")
    print(f"{len(ligs)} ligands, {len(t2 & set(ligs))} Type-II; 50 poses per ligand per pool\n")

    # --- headline: the per-pair interface confidence
    print("=== cp_iptm_pl  (protein<->ligand chain-pair ipTM; THE per-pair number) ===")
    hdr = f"{'ligand':<14}{'cls':<6}" + "".join(
        f"{p.replace('mol1_', ''):>18}" for p in live)
    print(hdr)
    print(f"{'':<20}" + "".join(f"{'median / best':>18}" for _ in live))
    for lig in ligs:
        cells = []
        for p in live:
            rows = pools[p].get(lig, [])
            v = [num(r, "cp_iptm_pl") for r in rows]
            v = [x for x in v if math.isfinite(x)]
            cells.append(f"{f(med(v))}/{f(max(v)) if v else '  -  '}".rjust(18))
        print(f"{lig:<14}{'T-II' if lig in t2 else 'non':<6}" + "".join(cells))
    print()
    for grp, sel in (("ALL", ligs), ("Type-II", [l for l in ligs if l in t2]),
                     ("non-Type-II", [l for l in ligs if l not in t2])):
        cells = []
        for p in live:
            v = [num(r, "cp_iptm_pl") for l in sel for r in pools[p].get(l, [])]
            v = [x for x in v if math.isfinite(x)]
            cells.append(f"{f(med(v))}".rjust(18))
        print(f"{'median ' + grp:<20}" + "".join(cells))

    # --- every other head, pooled per pool, so the reader can see what moves and what does not
    print("\n=== all confidence heads, pool-level median over 1000 poses ===")
    print(f"{'column':<16}" + "".join(f"{p.replace('mol1_', ''):>18}" for p in live) + "   note")
    notes = {"cp_iptm_pl": "per-PAIR interface -- quote this one",
             "ranking_score": "⚠ clash term is POLYMER-ONLY",
             "iptm": "whole complex, not pair-specific",
             "ptm": "global fold, ~insensitive to the ligand",
             "plddt": "dominated by the 503-aa protein",
             "gpde": "predicted distance error -- LOWER is better"}
    for c in COLS:
        cells = []
        for p in live:
            v = [num(r, c) for l in ligs for r in pools[p].get(l, [])]
            v = [x for x in v if math.isfinite(x)]
            cells.append(f"{f(med(v))}".rjust(18))
        print(f"{c:<16}" + "".join(cells) + f"   {notes[c]}")

    # --- spread: is cp_iptm_pl actually discriminating between ligands?
    print("\n=== does cp_iptm_pl separate ligands at all? (per-pool, across the 20 medians) ===")
    for p in live:
        meds = []
        for lig in ligs:
            v = [num(r, "cp_iptm_pl") for r in pools[p].get(lig, [])]
            v = [x for x in v if math.isfinite(x)]
            if v:
                meds.append((med(v), lig))
        if not meds:
            continue
        meds.sort()
        spread = meds[-1][0] - meds[0][0]
        print(f"  {p:<18} range {meds[0][0]:.3f} ({meds[0][1]}) .. {meds[-1][0]:.3f} "
              f"({meds[-1][1]})   spread {spread:.3f}")
    print("\n⚠ Confidence is the model's SELF-ASSESSMENT, not accuracy. On this target "
          "ranking_score has\n  measured BELOW CHANCE on Type-II and ipTM ~random against an "
          "oracle. Report, do not rank.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
