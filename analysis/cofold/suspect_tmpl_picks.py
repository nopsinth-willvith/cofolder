#!/usr/bin/env python3
"""Which submission_mol1_tmpl_mmgbsa picks are most likely WRONG, on physical grounds?

There is no reference structure for the blinded 20, so this cannot measure error. What it can do
is score each pick against expectations that are independent of how it was selected, and rank by
how many of those it violates. Every component is printed, so the ranking is auditable and a
reader can disagree with a weight rather than with a number.

The Fe expectations come from the campaign's own curated prior, where the coordinating and
non-coordinating populations are DISJOINT: coordinating <= 2.58 A, non-coordinating >= 2.72 A.

  Type-II ligands are DEFINED by coordinating the heme iron through a nitrogen, so for them:
    * Fe-ligand distance above 2.58 A  -> not coordinating. Scored by how far past it sits.
    * closest atom to Fe not N         -> coordinating through the wrong atom. Heavy flag.
    * Fe-ligand distance below 2.00 A  -> over-contracted past any real Fe-N bond.
  Non-Type-II ligands must stay OFF the iron, so for them:
    * Fe-ligand distance below 2.72 A  -> spuriously in the coordinating band. Heavy flag.

Confidence adds a third, weaker axis: cp_iptm_pl (the protein<->ligand chain-pair ipTM) and gpde
(predicted distance error, lower better), each scored against the 20-pick distribution.

⚠ WEAKEST LINK, STATED UP FRONT: on this target `cp_iptm_pl`/ipTM has measured roughly RANDOM
against a known oracle, so the confidence axis is deliberately weighted below the Fe geometry,
which rests on 122 curated crystal structures. A pick flagged ONLY on confidence is weak evidence.

⚠ This ranks SUSPICION, not error. A pick can violate nothing and still be wrong.

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
D, S = R / "data/cyp3a4_mol1", R / "data/cyp3a4_mol1/submission"
T2 = {l.strip() for l in (R / "data/cyp3a4_xmodel/typeii_base_ids.txt").read_text().split()
      if l.strip()}
COORD_MAX, NONCOORD_MIN, OVER_CONTRACT = 2.58, 2.72, 2.00
QC_POOLS = ("mol1_full", "mol1_tcryo5", "mol1_tcomb4", "mol1_t38FX")


def num(r, k):
    v = (r or {}).get(k, "")
    try:
        return float("nan") if v in ("", "nan", None) else float(v)
    except ValueError:
        return float("nan")


def main() -> int:
    conf = {}
    for p in QC_POOLS:
        fp = D / f"poses_qc_{p}.csv"
        if fp.is_file():
            for r in csv.DictReader(open(fp, newline="")):
                conf[r["pose_id"]] = r
    pick = {r["base_id"]: r for r in csv.DictReader(
        open(S / "submission_mol1_tmpl_mmgbsa.provenance.csv", newline=""))}

    cps = [num(conf.get(r["pose_id"], {}), "cp_iptm_pl") for r in pick.values()]
    gps = [num(conf.get(r["pose_id"], {}), "gpde") for r in pick.values()]
    cps = sorted(x for x in cps if math.isfinite(x))
    gps = sorted(x for x in gps if math.isfinite(x))
    cp_lo, cp_hi = cps[0], cps[-1]
    gp_lo, gp_hi = gps[0], gps[-1]

    rows = []
    for l, r in pick.items():
        c = conf.get(r["pose_id"], {})
        fe, el = num(r, "fe_lig_min_A"), (r.get("fe_lig_element") or "?")
        cp, gp, dt = num(c, "cp_iptm_pl"), num(c, "gpde"), num(r, "delta_total")
        flags, score = [], 0.0
        if l in T2:
            if math.isfinite(fe) and fe > COORD_MAX:
                # 1 point per 0.2 A past the coordinating cut, capped at 3
                pen = min(3.0, (fe - COORD_MAX) / 0.2)
                score += pen
                flags.append(f"Fe {fe:.3f} > {COORD_MAX} coord cut (+{pen:.2f})")
            if el != "N":
                score += 2.5
                flags.append(f"donor is {el}, not N (+2.50)")
            if math.isfinite(fe) and fe < OVER_CONTRACT:
                pen = min(2.0, (OVER_CONTRACT - fe) / 0.1)
                score += pen
                flags.append(f"Fe {fe:.3f} < {OVER_CONTRACT} over-contracted (+{pen:.2f})")
        else:
            if math.isfinite(fe) and fe < NONCOORD_MIN:
                pen = min(3.0, (NONCOORD_MIN - fe) / 0.2)
                score += pen
                flags.append(f"Fe {fe:.3f} < {NONCOORD_MIN} in the coord band (+{pen:.2f})")
        # confidence axis, max 1.0 total -- deliberately below any single Fe flag
        if math.isfinite(cp) and cp_hi > cp_lo:
            pen = 0.6 * (cp_hi - cp) / (cp_hi - cp_lo)
            score += pen
            if pen > 0.3:
                flags.append(f"cp_iptm_pl {cp:.3f} low in-set (+{pen:.2f})")
        if math.isfinite(gp) and gp_hi > gp_lo:
            pen = 0.4 * (gp - gp_lo) / (gp_hi - gp_lo)
            score += pen
            if pen > 0.2:
                flags.append(f"gpde {gp:.3f} high in-set (+{pen:.2f})")
        rows.append(dict(lig=l, t2=l in T2, fe=fe, el=el, cp=cp, gp=gp, dt=dt,
                         pool=r["pool"].replace("mol1_", ""), score=score, flags=flags))

    rows.sort(key=lambda x: -x["score"])
    print("RANKED BY SUSPICION -- most likely incorrect first\n")
    print(f"{'#':<3}{'ligand':<14}{'cls':<6}{'pool':<9}{'susp':>6}"
          f"{'Fe':>8}{'el':>4}{'cp_iptm':>9}{'gpde':>7}{'mmgbsa':>9}")
    for i, r in enumerate(rows, 1):
        print(f"{i:<3}{r['lig']:<14}{'T-II' if r['t2'] else 'non':<6}{r['pool']:<9}"
              f"{r['score']:>6.2f}{r['fe']:>8.3f}{r['el']:>4}{r['cp']:>9.3f}"
              f"{r['gp']:>7.3f}{r['dt']:>9.2f}")
        for f in r["flags"]:
            print(f"{'':<32}- {f}")
    print(f"\nclean (no Fe violation): "
          f"{[r['lig'] for r in rows if not any('Fe' in f or 'donor' in f for f in r['flags'])]}")
    print("\n⚠ suspicion, not error. The Fe axis rests on 122 curated crystals; the confidence "
          "axis is\n  weighted below it because ipTM has measured ~random against an oracle here.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
