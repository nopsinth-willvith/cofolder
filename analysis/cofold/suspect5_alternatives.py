#!/usr/bin/env python3
"""For the most-suspicious tmpl_mmgbsa picks: is a BETTER-and-DIFFERENT pose actually available?

Suspicion on its own does not justify changing a pick. What matters is whether the alternatives
are a genuinely different binding mode or just a reorientation of the same one -- a low-confidence
pick whose alternatives all sit ~2 A away in the same sub-pocket is not worth touching.

For each ligand this reports, per pool, the best-scoring pose with:
  * delta_total and the margin vs the current pick,
  * cp_iptm_pl and gpde (the model's own confidence),
  * RMSD and centroid shift against the current pick.

RMSD is heavy-atom, matched by ATOM NAME, with NO symmetry correction. Mol-1 names ligand atoms
deterministically (element + 1-based per-element occurrence over the SMILES atom order), so the
same ligand carries the same names in every pool and the match is exact. For a symmetric
substituent this OVERSTATES the RMSD, so read it as an upper bound; the centroid shift is
symmetry-free and is the honest test of "same place or not".

Interpretation used below, from the OCNT-2395645 case where the two poses were 2.6 A apart with
70% shared contacts: centroid shift much smaller than RMSD = reorientation in place; centroid
shift comparable to RMSD = the ligand actually moved.

⚠ No reference structure exists for these ligands. Nothing here measures correctness.

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
D, S = R / "data/cyp3a4_mol1", R / "data/cyp3a4_mol1/submission"
POSES = D / "poses"
SUSPECT = ["OCNT-2312966", "OCNT-2395328", "OCNT-2315379", "OCNT-2395460", "OCNT-2535422"]
T2 = {l.strip() for l in (R / "data/cyp3a4_xmodel/typeii_base_ids.txt").read_text().split()
      if l.strip()}


def num(r, k):
    v = (r or {}).get(k, "")
    try:
        return float("nan") if v in ("", "nan", None) else float(v)
    except ValueError:
        return float("nan")


def lig_atoms(pdb: Path) -> dict[str, tuple[float, float, float]]:
    """Heavy ligand atoms of chain C, keyed by atom name. Hydrogens are not written."""
    out = {}
    if not pdb.is_file():
        return out
    for ln in pdb.read_text().splitlines():
        if ln[:6] not in ("ATOM  ", "HETATM") or len(ln) < 54:
            continue
        if ln[21] != "C" or ln[17:20].strip() != "LIG":
            continue
        el = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        out[ln[12:16].strip()] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return out


def rmsd_and_com(a: dict, b: dict) -> tuple[float, float]:
    keys = sorted(set(a) & set(b))
    if not keys:
        return float("nan"), float("nan")
    ss = sum((a[k][i] - b[k][i]) ** 2 for k in keys for i in range(3))
    ca = [sum(a[k][i] for k in keys) / len(keys) for i in range(3)]
    cb = [sum(b[k][i] for k in keys) / len(keys) for i in range(3)]
    com = math.sqrt(sum((ca[i] - cb[i]) ** 2 for i in range(3)))
    return math.sqrt(ss / len(keys)), com


def find_pdb(pose_id: str) -> Path:
    for d in POSES.iterdir() if POSES.is_dir() else []:
        p = d / f"{pose_id}.pdb"
        if p.is_file():
            return p
    return POSES / "missing" / f"{pose_id}.pdb"


def main() -> int:
    conf = {}
    for p in sorted(D.glob("poses_qc_*.csv")):
        for r in csv.DictReader(open(p, newline="")):
            conf[r["pose_id"]] = r
    tables = [q for d in ("data/cyp3a4_ob0_mmgbsa", "data/cyp3a4_mol1")
              for q in sorted((R / d).glob("scored_*.csv")) if "ph48" not in q.name]
    pose = {}
    for t in tables:
        for r in csv.DictReader(open(t, newline="")):
            if r.get("base_id") in SUSPECT and math.isfinite(num(r, "delta_total")):
                pose.setdefault(r["pose_id"], r)
    pick = {r["base_id"]: r for r in csv.DictReader(
        open(S / "submission_mol1_tmpl_mmgbsa.provenance.csv", newline=""))}

    verdicts = []
    for l in SUSPECT:
        cur = pick[l]
        ref = lig_atoms(find_pdb(cur["pose_id"]))
        dcur = num(cur, "delta_total")
        ccur = num(conf.get(cur["pose_id"], {}), "cp_iptm_pl")
        print(f"=== {l}  ({'Type-II' if l in T2 else 'non-Type-II'})")
        print(f"    pick {cur['pool'].replace('mol1_','')} dT {dcur:+.2f} cp_iptm {ccur:.3f} "
              f"({len(ref)} heavy ligand atoms)")
        best = {}
        for r in pose.values():
            if r["base_id"] != l:
                continue
            p = r["pool"]
            if p not in best or num(r, "delta_total") < num(best[p], "delta_total"):
                best[p] = r
        print(f"    {'pool':<20}{'dT':>8}{'vs pick':>9}{'cp_iptm':>9}{'gpde':>7}"
              f"{'RMSD':>8}{'COM':>7}  mode")
        for p, r in sorted(best.items(), key=lambda kv: num(kv[1], "delta_total")):
            if p == cur["pool"]:
                continue
            c = conf.get(r["pose_id"], {})
            rm, com = rmsd_and_com(ref, lig_atoms(find_pdb(r["pose_id"])))
            mode = "-"
            if math.isfinite(rm) and math.isfinite(com):
                mode = ("MOVED" if com > 0.5 * rm else
                        "reorientation" if rm > 1.0 else "near-identical")
            d = num(r, "delta_total")
            print(f"    {p.replace('mol1_',''):<20}{d:>8.2f}{d - dcur:>+9.2f}"
                  f"{num(c,'cp_iptm_pl'):>9.3f}{num(c,'gpde'):>7.3f}"
                  f"{rm:>8.2f}{com:>7.2f}  {mode}")
        # the best CONFIDENCE alternative, which is what suspicion would motivate taking
        alts = [(num(conf.get(r['pose_id'], {}), 'cp_iptm_pl'), p, r)
                for p, r in best.items() if p != cur["pool"]
                and math.isfinite(num(conf.get(r['pose_id'], {}), 'cp_iptm_pl'))]
        if alts:
            cb, pb, rb = max(alts)
            rm, com = rmsd_and_com(ref, lig_atoms(find_pdb(rb["pose_id"])))
            verdicts.append((l, pb, cb - ccur, num(rb, "delta_total") - dcur, rm, com))
            print(f"    -> best-CONFIDENCE alternative: {pb.replace('mol1_','')} "
                  f"cp_iptm {cb:+.3f} vs pick, energy {num(rb,'delta_total') - dcur:+.2f}, "
                  f"RMSD {rm:.2f} A, COM {com:.2f} A")
        print()

    print("=== summary: what changing each suspicious pick would actually buy ===")
    print(f"{'ligand':<14}{'from':<12}{'d cp_iptm':>11}{'d energy':>10}{'RMSD':>7}{'COM':>7}  worth it?")
    for l, p, dc, de, rm, com in verdicts:
        worth = "NO -- pays energy for a reorientation"
        if math.isfinite(com) and math.isfinite(rm) and com > 0.5 * rm and de < 1.0 and dc > 0.01:
            worth = "maybe -- genuinely different pose, small energy cost"
        elif de < -0.1 and dc > 0.01:
            worth = "YES -- better on both"
        print(f"{l:<14}{p.replace('mol1_',''):<12}{dc:>+11.3f}{de:>+10.2f}{rm:>7.2f}{com:>7.2f}  {worth}")
    print("\n⚠ RMSD is name-matched without symmetry correction -- an upper bound. COM is exact.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
