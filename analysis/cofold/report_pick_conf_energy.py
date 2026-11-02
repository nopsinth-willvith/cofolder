#!/usr/bin/env python3
"""Per-pair structural confidence AND MM-GBSA for two candidate poses per ligand.

For every blinded ligand, side by side:
  * the pose submission_mol1_tmplonly_mmgbsa actually shipped, and
  * the best Fe-GATED pose available in the mol1_tsim4 pool.

Confidence comes from the QC tables (keyed on pose_id); energies from the scored tables.
`cp_iptm_pl` is the per-PAIR number -- chain-pair ipTM for the protein<->ligand pair. plddt is
whole-complex and dominated by the 503-aa protein; gpde is a predicted error, so LOWER is better;
srs (ranking_score) has a POLYMER-ONLY clash term and has measured below chance on Type-II here.

⚠ Confidence is the model's self-assessment and `delta_total` is a force-field estimate. Neither
is accuracy: the blinded 20 have no reference structure.

⚠ The Fe gate's second test is `fe_lig_out_A` = the Fe-ligand distance AFTER the MM-GBSA
minimisation, not a coordination bond (build_submission.fe_gate's own docstring). With
--no-fe-bond the minimiser always pulls the donor off the iron, so this is a "least-escaped"
ranking cut. Both the gated and the ungated tsim4 best are reported where they differ.

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

R = Path("/vast.mnt/project/bme-20255165-cofolding-physics")
D = R / "data/cyp3a4_mol1"
G = R / "data/cyp3a4_ob0_mmgbsa"
S = D / "submission"
T2 = {l.strip() for l in (R / "data/cyp3a4_xmodel/typeii_base_ids.txt").read_text().split()
      if l.strip()}
QC_POOLS = ("mol1_tcryo5", "mol1_tcomb4", "mol1_t38FX", "mol1_tsim4")
FE_IN, FE_RETAIN = 2.6, 3.5


def num(r, k) -> float:
    if not r:
        return float("nan")
    v = r.get(k, "")
    try:
        return float("nan") if v in ("", "nan", None) else float(v)
    except ValueError:
        return float("nan")


def f(x, nd=3) -> str:
    return "   -  " if not math.isfinite(x) else f"{x:.{nd}f}"


def main() -> int:
    argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                            formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()

    # pose_id -> confidence row, across every templated pool
    conf: dict[str, dict] = {}
    for p in QC_POOLS:
        fp = D / f"poses_qc_{p}.csv"
        if not fp.is_file():
            continue
        with open(fp, newline="") as fh:
            for r in csv.DictReader(fh):
                conf[r["pose_id"]] = r

    pick = {r["base_id"]: r for r in csv.DictReader(
        open(S / "submission_mol1_tmplonly_mmgbsa.provenance.csv", newline=""))}

    best_g, best_a = {}, {}
    with open(G / "scored_mol1_tsim4.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            dt = num(r, "delta_total")
            if not math.isfinite(dt):
                continue
            l = r["base_id"]
            if l not in best_a or dt < num(best_a[l], "delta_total"):
                best_a[l] = r
            if l in T2:
                fi, fo = num(r, "fe_lig_in_A"), num(r, "fe_lig_out_A")
                ok = math.isfinite(fi) and fi <= FE_IN and (not math.isfinite(fo) or fo <= FE_RETAIN)
            else:
                ok = True
            if ok and (l not in best_g or dt < num(best_g[l], "delta_total")):
                best_g[l] = r

    ligs = sorted(pick)
    print("A = tmplonly PICK (the shipped pose)      B = tsim4 BEST, Fe-GATED")
    print(f"{len(ligs)} ligands, {len(T2 & set(ligs))} Type-II\n")
    print(f"{'ligand':<14}{'cls':<6}{'source':<9}"
          f"{'cp_iptm_pl':>11}{'srs':>8}{'plddt':>8}{'gpde':>7}"
          f"{'delta_tot':>11}{'Fe_in':>8}{'Fe_out':>8}")
    agg = {"A": [], "B": []}
    for l in ligs:
        for tag, row, src in (("A", pick[l], pick[l]["pool"].replace("mol1_", "")),
                              ("B", best_g.get(l), "tsim4")):
            if row is None:
                print(f"{l:<14}{'T-II' if l in T2 else 'non':<6}{src:<9}  -- no gated candidate --")
                continue
            c = conf.get(row["pose_id"], {})
            ci, dt = num(c, "cp_iptm_pl"), num(row, "delta_total")
            agg[tag].append((ci, dt))
            lab = l if tag == "A" else ""
            cls = ("T-II" if l in T2 else "non") if tag == "A" else ""
            print(f"{lab:<14}{cls:<6}{src:<9}"
                  f"{f(ci):>11}{f(num(row,'srs')):>8}{f(num(c,'plddt'),2):>8}"
                  f"{f(num(c,'gpde')):>7}{f(dt,2):>11}"
                  f"{f(num(row,'fe_lig_in_A') if tag=='B' else num(row,'fe_lig_min_A')):>8}"
                  f"{f(num(row,'fe_lig_out_A')):>8}")
        ba = best_a.get(l)
        if ba is not None and best_g.get(l) is not None and ba["pose_id"] != best_g[l]["pose_id"]:
            c = conf.get(ba["pose_id"], {})
            print(f"{'':<14}{'':<6}{'(tsim4 ungated — REJECTED by fe_lig_out_A)':<9}"
                  f"{f(num(c,'cp_iptm_pl')):>11}{f(num(ba,'srs')):>8}"
                  f"{f(num(c,'plddt'),2):>8}{f(num(c,'gpde')):>7}"
                  f"{f(num(ba,'delta_total'),2):>11}{f(num(ba,'fe_lig_in_A')):>8}"
                  f"{f(num(ba,'fe_lig_out_A')):>8}")
        print()

    print(f"{'':<14}{'':<6}{'':<9}{'cp_iptm_pl':>11}{'':>8}{'':>8}{'':>7}{'delta_tot':>11}")
    for tag, nm in (("A", "mean tmplonly"), ("B", "mean tsim4 gated")):
        v = [(a, b) for a, b in agg[tag] if math.isfinite(a) and math.isfinite(b)]
        if v:
            print(f"{nm:<29}{sum(a for a,_ in v)/len(v):>11.3f}{'':>23}"
                  f"{sum(b for _,b in v)/len(v):>11.2f}   (n={len(v)})")
    print("\n⚠ confidence is the model's self-assessment and delta_total a force-field estimate; "
          "neither is accuracy.\n⚠ fe_lig_out_A is POST-minimisation, a least-escaped ranking cut, "
          "not a coordination bond.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
