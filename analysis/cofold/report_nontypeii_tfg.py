#!/usr/bin/env python3
"""Focused report on the NON-Type-II blinded ligands: TFG-guided vs both unguided zips.

The Fe-coordination disqualification that rules TFG out for the 11 Type-II ligands does not
apply here: a non-Type-II ligand is not expected to coordinate the heme iron, so the 2.58 A
prior is not a criterion for it. That makes this subset the one place the guided arm can be
judged on its merits.

Every statistic is recomputed at n=9 rather than read off a 20-ligand class split, because a
class-split median says nothing about significance within the class.

Reuses the readers and the exact sign test from compare_zip_metrics rather than restating them.

⚠ AGREEMENT, NOT ACCURACY for every lDDT / L-RMSD number: the blinded 20 have no reference.

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import argparse
import os
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_zip_metrics import med, num, read, read_by_pose, sign_test  # noqa: E402

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
S = R / "data/cyp3a4_mol1/submission"
G = R / "data/cyp3a4_ob0_mmgbsa"
D = R / "data/cyp3a4_mol1"
T2_FILE = R / "data/cyp3a4_xmodel/typeii_base_ids.txt"
TERMS = ("delta_total", "vdw", "eel", "egb", "esurf")


def f(x: float, nd: int = 2) -> str:
    return "-" if not math.isfinite(x) else f"{x:.{nd}f}"


def pb_rate(pool: str) -> dict[str, tuple[int, int]]:
    """base_id -> (pb_pass, total) from a stage-30 QC table."""
    p = D / f"poses_qc_{pool}.csv"
    if not p.is_file():
        return {}
    tru = lambda v: str(v).strip().lower() in ("true", "1", "yes")
    out: dict[str, list[int]] = {}
    with open(p, newline="") as fh:
        for r in csv.DictReader(fh):
            e = out.setdefault(r["base_id"], [0, 0])
            e[1] += 1
            if tru(r.get("pb_pass")):
                e[0] += 1
    return {k: (v[0], v[1]) for k, v in out.items()}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-md", type=Path, default=S / "nontypeii_tfg_report.md")
    a = ap.parse_args()

    t2 = {l.strip() for l in open(T2_FILE) if l.strip()}
    P = {k: read(S / f"submission_mol1_{k}_mmgbsa.provenance.csv")
         for k in ("tfg", "tmplonly", "tmpl")}
    SC = {"tfg": read_by_pose(G / "scored_mol1_tfg.csv"),
          "ref": read_by_pose(D / "scored_mol1_tmpl_full.csv")}
    LD = {}
    for k in ("tfg_vs_tmpl", "tmpl_vs_tfg", "tfg_vs_tmplonly", "tmplonly_vs_tfg"):
        p = S / f"lddt_{k}.csv"
        LD[k] = read(p) if p.is_file() else {}

    ligs = sorted(l for l in P["tfg"] if l not in t2)
    if not ligs:
        sys.exit("FATAL: no non-Type-II ligand found")
    qc_g, qc_c = pb_rate("mol1_tfg_tcryo5"), pb_rate("mol1_tcryo5")

    out: list[str] = []
    w = out.append
    w("# Non-Type-II blinded ligands: TFG-guided vs unguided\n")
    w(f"**{len(ligs)} of 20 ligands.** The heme-Fe coordination requirement that disqualifies "
      "TFG for the 11 Type-II ligands does not apply to these: a non-Type-II ligand is not "
      "expected to coordinate the iron, so the 2.58 Å prior is not a criterion for it. This is "
      "the subset where the guided arm can be judged on its merits.\n")
    w("All statistics are recomputed at **n=9**, not read off a 20-ligand class split.\n")
    w("> ⚠ **Agreement, not accuracy.** The blinded 20 have no reference structure, so every "
      "lDDT / L-RMSD figure compares two predictions.\n")
    w("Zips, all built with the identical selector (`delta_total`, Fe gate 2.6/3.5, ladder "
      "3.0:4.0,3.5:4.5):\n")
    w("- **tfg** — pool `mol1_tfg_tcryo5` only (TFG ON, stock `mol-1-base.pt`, cryo5 templates)")
    w("- **tmplonly** — the 3 templated unguided pools")
    w("- **tmpl** — those 3 plus template-free `mol1_full` (the campaign's best measured zip)\n")

    # 1. pool diversity available to the selector
    w("## 1. Valid-pose supply (PoseBusters pass rate, 50 poses per ligand)\n")
    w("This is what TFG unambiguously changes: how many poses the selector may choose from.\n")
    w("| ligand | guided `mol1_tfg_tcryo5` | unguided `mol1_tcryo5` | gain |")
    w("|---|---:|---:|---:|")
    gains = []
    for l in ligs:
        g, c = qc_g.get(l), qc_c.get(l)
        if not g or not c:
            w(f"| {l} | - | - | - |")
            continue
        gp, cp = 100 * g[0] / g[1], 100 * c[0] / c[1]
        gains.append(gp - cp)
        w(f"| {l} | {g[0]}/{g[1]} = {gp:.0f}% | {c[0]}/{c[1]} = {cp:.0f}% | **+{gp - cp:.0f} pp** |")
    if gains:
        p_, n_, t_, pv = sign_test(gains)
        w(f"\nMedian gain **+{med(gains):.0f} pp**; {p_} of {len(gains)} ligands improve, "
          f"{n_} worsen, {t_} tied (sign test p={pv:.3f}).")

    # 2. MM-GBSA on the shipped poses
    w("\n## 2. MM-GBSA of the shipped poses (kcal/mol, lower = tighter)\n")
    w("| ligand | tfg dTOT | tmpl dTOT | Δ | tfg vdw | tmpl vdw | tfg eel | tmpl eel |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|")
    dd: dict[str, list[float]] = {t: [] for t in TERMS}
    for l in ligs:
        ra = SC["tfg"].get(P["tfg"][l]["pose_id"])
        rb = SC["ref"].get(P["tmpl"][l]["pose_id"])
        for t in TERMS:
            d = num(ra, t) - num(rb, t)
            if math.isfinite(d):
                dd[t].append(d)
        dtot = num(ra, "delta_total") - num(rb, "delta_total")
        w(f"| {l} | {f(num(ra, 'delta_total'))} | {f(num(rb, 'delta_total'))} | "
          f"{'-' if not math.isfinite(dtot) else format(dtot, '+.2f')} | {f(num(ra, 'vdw'))} | "
          f"{f(num(rb, 'vdw'))} | {f(num(ra, 'eel'))} | {f(num(rb, 'eel'))} |")
    w("\n| term | median Δ (tfg − tmpl) | favours tfg | favours tmpl | p |")
    w("|---|---:|---:|---:|---:|")
    for t in TERMS:
        d = dd[t]
        p_, n_, t_, pv = sign_test([-x for x in d])   # lower energy = better
        w(f"| `{t}` | {med(d):+.2f} | {p_} | {n_} | {pv:.3f} |")

    # 3. other metrics
    w("\n## 3. Other shipped-pose metrics\n")
    w("| ligand | tfg srs | tmpl srs | tfg maxCC | tmpl maxCC | tfg Fe | tmpl Fe | tfg pool seed |")
    w("|---|---:|---:|---:|---:|---:|---:|---|")
    for l in ligs:
        pa, pb = P["tfg"][l], P["tmpl"][l]
        w(f"| {l} | {f(num(pa, 'srs'), 3)} | {f(num(pb, 'srs'), 3)} | "
          f"{f(num(pa, 'max_cc_bond_A'), 3)} | {f(num(pb, 'max_cc_bond_A'), 3)} | "
          f"{f(num(pa, 'fe_lig_min_A'), 3)} | {f(num(pb, 'fe_lig_min_A'), 3)} | "
          f"s{pa.get('seed', '?')}/{pa.get('sample', '?')} |")
    for k, lower in (("srs", False), ("max_cc_bond_A", True)):
        d = [num(P["tfg"][l], k) - num(P["tmpl"][l], k) for l in ligs]
        d = [x for x in d if math.isfinite(x)]
        p_, n_, t_, pv = sign_test([-x for x in d] if lower else d)
        w(f"\n- `{k}`: median Δ {med(d):+.4f}, tfg/tmpl {p_}/{n_}, p={pv:.3f} "
          f"({'lower' if lower else 'higher'} is better)")

    # 4. agreement
    w("\n## 4. Agreement with each unguided zip\n")
    w("| ligand | PLI tfg\\|tmpl | PLI tmpl\\|tfg | LP tfg\\|tmpl | L-RMSD vs tmpl | "
      "L-RMSD vs tmplonly |")
    w("|---|---:|---:|---:|---:|---:|")
    rms = []
    for l in ligs:
        fw, rv = LD["tfg_vs_tmpl"].get(l), LD["tmpl_vs_tfg"].get(l)
        fo = LD["tfg_vs_tmplonly"].get(l)
        r = num(fw, "portal_bisyrmsd")
        if math.isfinite(r):
            rms.append(r)
        w(f"| {l} | {f(num(fw, 'portal_lddt_pli'), 3)} | {f(num(rv, 'portal_lddt_pli'), 3)} | "
          f"{f(num(fw, 'portal_lddt_lp'), 3)} | {f(r, 3)} | "
          f"{f(num(fo, 'portal_bisyrmsd'), 3)} |")
    if rms:
        w(f"\nMean L-RMSD vs tmpl **{sum(rms) / len(rms):.3f} Å**; "
          f"{sum(1 for x in rms if x >= 2.0)}/{len(rms)} differ by ≥2 Å.")

    a.out_md.write_text("\n".join(out) + "\n")
    print("\n".join(out))
    print(f"\nwrote {a.out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
