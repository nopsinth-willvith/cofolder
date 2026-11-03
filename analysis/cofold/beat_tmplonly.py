#!/usr/bin/env python3
"""Is any scored Mol-1 pose better than what submission_mol1_tmplonly_mmgbsa shipped?

tmplonly is argmin(delta_total) over only THREE pools (tcryo5, tcomb4, t38FX). Several other
pools have been MM-GBSA scored since and never competed against it. This asks, per ligand and
per pool, whether a better-scoring pose exists that would also SURVIVE the Type-II Fe gate.

⚠ "better" here means lower MM-GBSA delta_total under the same Fe gate. On a blinded target set
there is no reference structure, so this is one proxy agreeing or disagreeing with itself -- not
accuracy. Treat a margin smaller than the campaign's own near-tie band (~0.1 kcal/mol, the scale
at which the pH 4.8 re-rank flipped picks) as no difference at all.

⚠ NOT every pool is the same experiment. mol1ft_* is a DIFFERENT CHECKPOINT; mol1_febond* declare
a covalent Fe bond in the input; mol1_repair rows are grafted from another pool. A better pose
from those is a different METHOD winning, not a better template -- which is a bigger claim and is
flagged in the output rather than folded into the total.

pH 4.8 tables are excluded: a different protonation state is not comparable on delta_total.

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
S = R / "data/cyp3a4_mol1/submission"
T2 = {l.strip() for l in (R / "data/cyp3a4_xmodel/typeii_base_ids.txt").read_text().split()
      if l.strip()}
FE_IN, FE_RETAIN = 2.6, 3.5
NEAR_TIE = 0.10          # below this, call it no difference
SAME_CKPT = {"mol1_full", "mol1_tcryo5", "mol1_tcomb4", "mol1_t38FX", "mol1_tsim4"}
OTHER_METHOD = {"mol1ft_tcryo5": "fine-tuned checkpoint",
                "mol1ft_tcomb4": "fine-tuned checkpoint",
                "mol1_febond": "declared covalent Fe bond",
                "mol1_febond_tcryo5": "declared covalent Fe bond",
                "mol1_febond_tcomb4": "declared covalent Fe bond",
                "mol1_tfg_tcryo5": "TFG-guided sampler",
                "mol1_repair": "pose grafted from another pool"}
TABLES = ["data/cyp3a4_mol1/scored_mol1_tmpl_full.csv",
          "data/cyp3a4_ob0_mmgbsa/scored_mol1_tsim4.csv",
          "data/cyp3a4_ob0_mmgbsa/scored_mol1_tfg.csv",
          "data/cyp3a4_mol1/scored_mol1_expanded_full.csv"]


def num(r, k):
    v = (r or {}).get(k, "")
    try:
        return float("nan") if v in ("", "nan", None) else float(v)
    except ValueError:
        return float("nan")


def gated(r) -> bool:
    if r["base_id"] not in T2:
        return True
    fi, fo = num(r, "fe_lig_in_A"), num(r, "fe_lig_out_A")
    return math.isfinite(fi) and fi <= FE_IN and (not math.isfinite(fo) or fo <= FE_RETAIN)


def main() -> int:
    pose: dict[str, dict] = {}
    for t in TABLES:
        p = R / t
        if not p.is_file():
            print(f"⚠ missing {t}")
            continue
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                if math.isfinite(num(r, "delta_total")):
                    pose.setdefault(r["pose_id"], r)      # first wins; identical scores
    print(f"{len(pose)} distinct scored poses across {len(TABLES)} tables\n")

    pick = {r["base_id"]: r for r in csv.DictReader(
        open(S / "submission_mol1_tmplonly_mmgbsa.provenance.csv", newline=""))}

    best: dict[str, dict[str, dict]] = {}
    for r in pose.values():
        if not gated(r):
            continue
        l, pl = r["base_id"], r["pool"]
        if l not in pick:
            continue
        d = best.setdefault(l, {})
        if pl not in d or num(r, "delta_total") < num(d[pl], "delta_total"):
            d[pl] = r

    print(f"{'ligand':<14}{'cls':<6}{'tmplonly':>10}{'best elsewhere':>16}{'margin':>9}"
          f"{'from pool':<22} note")
    beaten_same, beaten_other = [], []
    for l in sorted(pick):
        dp = num(pick[l], "delta_total")
        cand = {p: v for p, v in best.get(l, {}).items()
                if p not in ("mol1_tcryo5", "mol1_tcomb4", "mol1_t38FX")}
        if not cand:
            continue
        bp = min(cand, key=lambda p: num(cand[p], "delta_total"))
        db = num(cand[bp], "delta_total")
        marg = db - dp
        note = OTHER_METHOD.get(bp, "")
        if marg < -NEAR_TIE:
            (beaten_other if bp in OTHER_METHOD else beaten_same).append((l, bp, marg))
            flag = "  <-- BEATS tmplonly"
        elif marg < 0:
            flag = "  (within the near-tie band)"
        else:
            flag = ""
        print(f"{l:<14}{'T-II' if l in T2 else 'non':<6}{dp:>10.2f}{db:>16.2f}{marg:>+9.2f}"
              f"{bp.replace('mol1_',''):<22}{note}{flag}")

    print(f"\n=== verdict ===")
    print(f"same checkpoint, ordinary co-folding (full / tsim4): beats tmplonly on "
          f"{len(beaten_same)}/20")
    for l, p, m in sorted(beaten_same, key=lambda x: x[2]):
        print(f"    {l}  {p.replace('mol1_',''):<12} {m:+.2f} kcal/mol")
    print(f"a DIFFERENT method beats tmplonly on {len(beaten_other)}/20 "
          f"(bigger claim -- see the note column)")
    for l, p, m in sorted(beaten_other, key=lambda x: x[2]):
        print(f"    {l}  {p.replace('mol1_',''):<12} {m:+.2f}   [{OTHER_METHOD[p]}]")
    tot = sum(m for _, _, m in beaten_same + beaten_other)
    print(f"\ntotal recoverable delta_total if every winner were taken: {tot:+.2f} kcal/mol")
    print("⚠ delta_total only. No reference structure exists for these ligands.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
