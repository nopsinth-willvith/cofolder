#!/usr/bin/env python3
"""Compare two Mol-1 pools on LIGAND GEOMETRY quality, paired per ligand.

Built to answer one question: does the TFG chemistry guidance produce better ligand geometry than
the same arm without it? The guidance's active terms on this target are all ligand-internal
(ExperimentalTorsion, PlanarImproper, PairwiseDistance), because geometry_featurizer skips standard
AND metal-containing residues -- HEM contributes zero RDKit terms. So pb_pass and max_cc_bond_A are
the metrics with a mechanism behind them; Fe distance is reported for context only, and a change
there would be a side effect, not the intended one.

⚠⚠ THE PAIRING UNIT IS THE LIGAND, NOT THE POSE. mol/tfg/engine.py:493 consumes an extra
torch.randn per diffusion step that it then discards, so a guided run desynchronizes from an
unguided one at the same seed. Pose (ligand, seed, sample) in one arm is therefore NOT the same
trajectory as in the other, and a pose-level paired test would be measuring an arbitrary alignment.
What IS comparable is the per-ligand aggregate: same ligand, same template, same checkpoint, same
50-pose budget. That gives n=20 paired observations, which is what the sign test below uses.

A difference here says "the guided arm is better", NOT "the guidance caused it" -- the guided arm
also resamples. Only a guidance-on/weights-zero control separates those (the pattern used in
MGBench/slurm/tfg_steer/protenix_tfg_infer.sbatch).

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import argparse
import csv
import math
import statistics as st
import sys
from collections import Counter
from pathlib import Path

R = Path("/vast.mnt/project/bme-20255165-cofolding-physics")
D = R / "data/cyp3a4_mol1"
TYPEII = set((R / "data/cyp3a4_xmodel/typeii_base_ids.txt").read_text().split())
# PoseBusters' own long-bond threshold territory; the campaign's worst C-C sits at 1.55-1.59 A.
CC_FLAG = 1.60


def load(pool: str) -> list[dict]:
    p = D / f"poses_qc_{pool}.csv"
    if not p.is_file():
        sys.exit(f"FATAL: no QC table for pool {pool} at {p}\n"
                 f"       run: export POOL={pool}; sbatch --export=ALL --job-name=$POOL "
                 f"slurm/cyp3a4_mol1/30_collect_qc.sbatch")
    rows = list(csv.DictReader(open(p, newline="")))
    if not rows:
        sys.exit(f"FATAL: {p} is empty")
    return rows


def fnum(r: dict, k: str):
    try:
        v = float(r[k])
    except (KeyError, TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def sign_test(deltas: list[float]) -> tuple[int, int, int, float]:
    """Two-sided exact sign test. Ties are dropped, which is the conservative convention."""
    pos = sum(d > 0 for d in deltas)
    neg = sum(d < 0 for d in deltas)
    ties = sum(d == 0 for d in deltas)
    n = pos + neg
    if n == 0:
        return pos, neg, ties, 1.0
    k = min(pos, neg)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return pos, neg, ties, min(1.0, 2 * tail)


def per_ligand(rows: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for r in rows:
        g = out.setdefault(r["base_id"], {"n": 0, "pass": 0, "cc": [], "fe": [], "srs": []})
        g["n"] += 1
        g["pass"] += str(r.get("pb_pass", "")).lower() == "true"
        for key, col in (("cc", "max_cc_bond_A"), ("fe", "fe_lig_min_A"), ("srs", "ranking_score")):
            v = fnum(r, col)
            if v is not None:
                g[key].append(v)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--guided", default="mol1_tfg_tcryo5", help="the TFG pool")
    ap.add_argument("--control", default="mol1_tcryo5", help="the matched unguided pool")
    # Restrict the paired tests to one ligand class. A 20-ligand class SPLIT reports medians but
    # says nothing about significance WITHIN a class, and the two classes behave differently on
    # this target (Type-II coordinate the heme iron, non-Type-II do not), so a pooled p-value can
    # hide a real effect in one class or manufacture one from the other.
    ap.add_argument("--only-class", choices=("t2", "non", "all"), default="all",
                    help="t2 = Type-II only, non = non-Type-II only (statistics recomputed at "
                         "that n, not read off a 20-ligand split)")
    a = ap.parse_args()

    RG, RC = load(a.guided), load(a.control)
    G, C = per_ligand(RG), per_ligand(RC)
    shared = sorted(set(G) & set(C))
    if not shared:
        sys.exit("FATAL: the two pools share no ligand")
    print(f"guided  = {a.guided:<22} {len(RG):>5} poses, {len(G)} ligands")
    print(f"control = {a.control:<22} {len(RC):>5} poses, {len(C)} ligands")
    if len(shared) < 20:
        print(f"⚠ only {len(shared)}/20 ligands present in BOTH pools -- the guided arm may be short")
    only = sorted(set(C) - set(G))
    if only:
        print(f"⚠ missing from the guided arm: {only}")
    if a.only_class != "all":
        keep = [l for l in shared if (l in TYPEII) == (a.only_class == "t2")]
        print(f"\nrestricted to {a.only_class.upper()}: {len(keep)} of {len(shared)} ligands "
              f"-- every statistic below is recomputed at n={len(keep)}")
        shared = keep
        if not shared:
            sys.exit(f"FATAL: no ligand in class {a.only_class}")

    print(f"\n{'ligand':<14}{'cls':<6}{'pb pass G|C':>13}{'maxCC med G|C':>18}"
          f"{'maxCC max G|C':>18}{'>1.60A G|C':>13}")
    d_pass, d_cc, d_worst = [], [], []
    for b in shared:
        g, c = G[b], C[b]
        gp, cp = g["pass"] / g["n"], c["pass"] / c["n"]
        gm, cm = st.median(g["cc"]), st.median(c["cc"])
        gx, cx = max(g["cc"]), max(c["cc"])
        gf = sum(v > CC_FLAG for v in g["cc"]) / len(g["cc"])
        cf = sum(v > CC_FLAG for v in c["cc"]) / len(c["cc"])
        d_pass.append(gp - cp)
        d_cc.append(cm - gm)        # positive = guided has SHORTER (better) bonds
        d_worst.append(cx - gx)
        print(f"{b:<14}{'T-II' if b in TYPEII else 'non':<6}"
              f"{gp:>6.0%}|{cp:<6.0%}{gm:>8.3f}|{cm:<9.3f}{gx:>8.3f}|{cx:<9.3f}"
              f"{gf:>5.0%}|{cf:<7.0%}")

    for label, rows, pool in (("guided", RG, a.guided), ("control", RC, a.control)):
        cc = [v for r in rows if (v := fnum(r, "max_cc_bond_A")) is not None]
        npass = sum(str(r.get("pb_pass", "")).lower() == "true" for r in rows)
        print(f"\n=== {label}: {pool} ===")
        print(f"  pb_pass           {npass}/{len(rows)} = {npass / len(rows):.1%}")
        print(f"  max_cc_bond_A     median {st.median(cc):.3f}  mean {st.fmean(cc):.3f}  "
              f"max {max(cc):.3f}  >{CC_FLAG}A {sum(v > CC_FLAG for v in cc) / len(cc):.1%}")
        fails = Counter(x for r in rows if (f := r.get("pb_failed", "")) and f != "-"
                        for x in f.split("|") if x)
        print(f"  top pb failures   {dict(fails.most_common(5)) or 'none'}")

    print(f"\n=== paired over {len(shared)} ligands (sign test, two-sided) ===")
    for name, d, better in (("pb_pass rate", d_pass, "guided passes more"),
                            ("median max_cc", d_cc, "guided bonds shorter"),
                            ("worst max_cc", d_worst, "guided worst-case shorter")):
        pos, neg, ties, p = sign_test(d)
        print(f"  {name:<16} {pos:>2} ligands favour guided, {neg:>2} favour control, "
              f"{ties} tied   p={p:.3f}   median delta {st.median(d):+.4f}  ({better})")
    print("\n⚠ 'better' here is agreement with a hypothesis, not accuracy: the blinded 20 have no")
    print("  reference, and the guided arm resamples as well as steers.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
