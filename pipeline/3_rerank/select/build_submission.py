"""Build a 20-member submission zip from per-pose physics scores, every member PoseBusters-checked.

Selection: per blinded ligand (the __R/__S queries pool under their base_id), poses from BOTH OB0
pools (template-free + templated) are ranked by --score-col, lower = better. With one column that
is the raw energy (e.g. MM-GBSA delta_total); with several, each is z-scored within the ligand's
candidate set (poses finite in all of them) and the z-scores are summed with equal weight -- the
combined MM-GBSA + tmol version.

--seeds / --pools restrict the CANDIDATE POOL before ranking, which is how an arm is made
comparable to an already-measured submission: `--seeds 101,102,103,104,105` is the 1300-pose pool
that produced `submission_combined.zip` (portal lDDT-PLI ~0.44), so re-ranking it swaps only the
selector, where the full seeds 101-110 pool would confound the ranker with twice the candidates.
Same flag, same meaning, as `analysis/cyp3a4_rerank/rerank.py --seeds`.

PoseBusters is RE-RUN on the exact file to be shipped (validate_frames.py: posebusters 0.6.5,
dock config, ligand bond orders from SMILES), never taken from the cached pool flag: candidates
are validated in rank order and the first that passes is the member. The coordinates shipped are
the ORIGINAL co-folded OB0 pose (the minimised MM-GBSA / tmol structures are not submission
format). A ligand with no passing pose among its top --topk is a hard failure.

Writes <out-dir>/<name>.zip (members <base_id>.pdb at the archive root, text copied unchanged)
and <name>.provenance.csv. The official validator, the from-the-zip audit and the portal
comparison run after this, in slurm/3_rerank/select.sbatch.

--scores takes one or more tables (one per scored pool, say); they are concatenated and a pose_id
appearing twice is fatal -- the same pose scored twice would compete against itself.

--shortlist OUT stops after ranking and writes the top --topk per ligand (same blinded-ligand,
pool/seed and Fe-gate filters as a real build) to OUT, without PoseBusters or a zip. It is how the
MLFF scorer picks the MM-GBSA top-K it re-scores, so both stages gate the candidate set identically.

Env: sigmadock. CPU.
"""

from __future__ import annotations

import argparse
import os
import csv
import hashlib
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))
BLINDED = R / "cyp-challenge-structure/cyp-challenge-TEST-BLINDED_structures.csv"
VALIDATE = Path(__file__).resolve().parents[2] / "2_cofold/validate_frames.py"  # the QC harness
TYPEII = set((R / "data/cyp3a4_xmodel/typeii_base_ids.txt").read_text().split())


def read_scores(paths: list[Path]) -> pd.DataFrame:
    """Concatenate per-pose score tables; a pose_id in more than one row is fatal."""
    df = pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)
    dup = df.pose_id[df.pose_id.duplicated()]
    if len(dup):
        sys.exit(f"FATAL: {dup.nunique()} pose_id(s) appear more than once across "
                 f"{[p.name for p in paths]}, e.g. {sorted(dup.unique())[:3]}")
    if len(paths) > 1:
        print(f"scores: {len(df)} rows from {len(paths)} tables")
    return df


def fe_gate(df: pd.DataFrame, in_max: float, out_max: float,
            ladder: tuple = ()) -> pd.DataFrame:
    """Keep only Type-II candidates whose Fe coordination survives the MM-GBSA minimisation.

    Non-Type-II rows pass untouched. A Type-II row is kept when it BOTH starts coordinating
    (`fe_lig_in_A` <= in_max, the campaign's 2.6 A Type-II band) and still holds the contact at the
    end (`fe_lig_out_A` <= out_max).

    The second test is the point: with --no-fe-bond the Fe-donor contact is carried by LJ and
    electrostatics alone, so the minimiser is free to pull the donor off the iron, and it always
    does -- over the 1234 scored Type-II poses the minimum `fe_lig_out_A` is 2.69 A, so a true
    coordination cut (2.6 A) keeps ZERO poses and 3.0 A keeps 41. out_max is therefore a
    "least-escaped" ranking cut, not a coordination bond, and 3.2 A is the smallest value that
    leaves every one of the 11 Type-II ligands with a candidate. An empty gate for any ligand is a
    hard failure rather than a silent fallback to the ungated pool.
    """
    if not (in_max or out_max):
        return df
    t2 = df.base_id.isin(TYPEII)

    def passes(d: pd.DataFrame, imax: float, omax: float) -> pd.Series:
        ok = pd.Series(True, index=d.index)
        if imax:
            ok &= pd.to_numeric(d.fe_lig_in_A, errors="coerce") <= imax
        if omax:
            ok &= pd.to_numeric(d.fe_lig_out_A, errors="coerce") <= omax
        return ok

    # The strict cut is applied PER LIGAND, and only a ligand that it empties is relaxed, one rung
    # at a time. Widening the cut for all eleven because one needs it would quietly buy every other
    # Type-II pick a worse pose; this way the relaxation is confined to where it is forced, and the
    # ligands that needed it are named so the provenance can carry it.
    ladder = [(in_max, out_max)] + [(i, o) for i, o in (ladder or ())]
    keep, relaxed, lost = [], {}, []
    for b, g in df[t2].groupby("base_id"):
        for imax, omax in ladder:
            ok = passes(g, imax, omax)
            if ok.any():
                keep.append(g[ok])
                if (imax, omax) != ladder[0]:
                    relaxed[b] = (imax, omax)
                break
        else:
            lost.append(b)
    if lost:
        sys.exit(f"FATAL: the Fe gate empties the candidate set for {sorted(lost)} even at the "
                 f"loosest rung {ladder[-1]}")
    out = pd.concat([df[~t2], *keep]) if keep else df[~t2]
    print(f"Fe gate in<={in_max} out<={out_max}: Type-II candidates "
          f"{int(t2.sum())} -> {int(out.base_id.isin(TYPEII).sum())}")
    if relaxed:
        print("  ⚠ relaxed, no pose at the strict cut: "
              + ", ".join(f"{b} -> in<={i} out<={o}" for b, (i, o) in sorted(relaxed.items())))
    return out


def restrict(df: pd.DataFrame, pools: str = "", seeds: str = "") -> pd.DataFrame:
    """Restrict a candidate pool to named pools and seeds, guarding the seed-namespace trap.

    Seeds are NOT globally unique across pools: the OB0 arms are seeded 101-110 and the ESMFold2
    arms 0-9, so a seed list written for one arm silently deletes every pose of another once both
    sit in the same merged table. A seed filter that empties a WHOLE pool is a mistake, not a
    selection -- name the arms with `pools` instead -- so it is fatal here rather than a quietly
    short candidate set later.

    Shared by build_submission, rank_combined_pool and fe_gated_report, which all read the same
    merged table: the filter was copied into each and the guard would otherwise protect only one.
    """
    if pools:
        df = df[df.pool.isin(pools.replace(",", " ").split())]
    if seeds:
        before = set(df.pool)
        df = df[df.seed.astype(int).isin({int(x) for x in seeds.replace(",", " ").split()})]
        gone = sorted(before - set(df.pool))
        if gone:
            sys.exit(f"FATAL: seeds {seeds} eliminate every pose of {gone} (seed namespaces differ "
                     f"per pool); name the arms with --pools instead")
    return df


def rank_table(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df[np.isfinite(df[cols]).all(axis=1)]
    if len(cols) == 1:
        df["score"] = df[cols[0]]
    else:
        # Mean/sd, NOT the median/MAD `robust_z` the cross-model selectors use (cyp3a4_xmodel/
        # select_norm_within_model.py): this is an equal-weight sum of physics energies within one
        # ligand, where the sd is the natural scale and there is no heavy tail to resist.
        # sd == 0 (a constant column, or a ligand left with one candidate by the Fe gate) would
        # make the z-score 0/0 = NaN, which np.isfinite above has already run and cannot catch:
        # the ranking would silently fall back to pose_id order with NaN in the provenance. A
        # degenerate column carries no ranking information, so it contributes 0.
        g = df.groupby("base_id")[cols]
        sd = g.transform("std", ddof=0)
        z = (df[cols] - g.transform("mean")).div(sd).where(sd > 0, 0.0)
        df["score"] = z.sum(axis=1)
    df = df.sort_values(["base_id", "score", "pose_id"])
    df["rank"] = df.groupby("base_id").cumcount() + 1
    return df


def main() -> int:
    # allow_abbrev=False: "--pool base" must not be silently accepted as --pools, and "--seed 101"
    # must not be accepted as --seeds -- either would quietly ship a different pool than intended.
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--scores", type=Path, nargs="+", required=True,
                    help="per-pose table(s): pose_id, pdb, base_id, ... (several are concatenated)")
    ap.add_argument("--score-col", required=True, help="one column, or several comma-separated (z-sum)")
    ap.add_argument("--name", required=True)
    ap.add_argument("--out-dir", type=Path, default=R / "data/cyp3a4_ob0_mmgbsa/submission")
    ap.add_argument("--topk", type=int, default=5)
    ap.add_argument("--pools", default="", help="comma list of pools to keep (default: all)")
    ap.add_argument("--seeds", default="", help="comma list of seeds to keep (default: all)")
    ap.add_argument("--fe-in-max", type=float, default=0.0,
                    help="Type-II only: require fe_lig_in_A <= this (2.6 = the Type-II band)")
    ap.add_argument("--fe-retain-max", type=float, default=0.0,
                    help="Type-II only: require fe_lig_out_A <= this, i.e. the contact survives "
                         "minimisation (3.2 is the smallest cut that leaves all 11 non-empty)")
    ap.add_argument("--fe-ladder", default="",
                    help="fallback rungs for a Type-II ligand the strict cut empties, as "
                         "'in:out,in:out' (e.g. '3.0:3.5,3.5:4.0'). Applied PER LIGAND and only "
                         "where forced; the relaxed ligands are named in the output.")
    ap.add_argument("--maximize", action="store_true",
                    help="score columns where HIGHER is better (a confidence, not an energy)")
    ap.add_argument("--shortlist", type=Path, default=None,
                    help="write the ranked top --topk per ligand here and stop (no PB, no zip)")
    args = ap.parse_args()
    cols = args.score_col.split(",")
    with BLINDED.open(newline="") as fh:
        ids = {r["Molecule_Name"] for r in csv.DictReader(fh)}

    df = read_scores(args.scores)
    if "status" in df:
        df = df[df.status == "ok"]
    df["seed"] = df.seed.astype(int)
    # Only the 20 blinded ligands may compete. `missing` below catches a ligand with no pose, but
    # not the reverse: an oracle pool merged into the same table brings 21 PDB-code base_ids, which
    # without this would be ranked, PoseBusters-checked and WRITTEN INTO THE ZIP before the member
    # count check at the end rejected it. The count is printed so the drop stays auditable.
    n_all = len(df)
    df = df[df.base_id.isin(ids)]
    df = restrict(df, args.pools, args.seeds)
    ladder = tuple(tuple(float(x) for x in rung.split(":"))
                   for rung in args.fe_ladder.replace(",", " ").split()) if args.fe_ladder else ()
    df = fe_gate(df, args.fe_in_max, args.fe_retain_max, ladder)
    if args.maximize:
        # rank_table is written for energies, where lower is better. A confidence is the other way
        # round, so the columns are negated here rather than giving rank_table a direction flag:
        # the z-sum path then still works, and `score` in the provenance keeps one consistent
        # meaning (smaller is the better pose) across every submission this script has built.
        for c in cols:
            df[c] = -pd.to_numeric(df[c], errors="coerce")
        print(f"--maximize: negated {cols}; provenance `score` stays 'lower is better'")
    print(f"candidate pool: {len(df)} poses, {df.base_id.nunique()} ligands, "
          f"pools {sorted(df.pool.unique())}, seeds {sorted(df.seed.unique())} "
          f"({n_all - len(df)} rows dropped: non-blinded ligand, pool, seed or Fe gate)")
    ranked = rank_table(df, cols)
    missing = ids - set(ranked.base_id)
    if missing:
        sys.exit(f"FATAL: no scored pose for {sorted(missing)}")
    cand = ranked[ranked["rank"] <= args.topk]
    if args.shortlist:
        args.shortlist.parent.mkdir(parents=True, exist_ok=True)
        cand.to_csv(args.shortlist, index=False)
        print(f"shortlist: top {args.topk} of {len(ranked)} ranked poses -> {len(cand)} rows, "
              f"{cand.base_id.nunique()} ligands -> {args.shortlist}")
        return 0

    work = args.out_dir / f".{args.name}_work"
    work.mkdir(parents=True, exist_ok=True)
    cand[["base_id", "pdb"]].assign(run=cand.pose_id).to_csv(work / "candidates.csv", index=False)
    subprocess.run([sys.executable, str(VALIDATE), "--frames-csv", str(work / "candidates.csv"),
                    "--out-csv", str(work / "candidates_pb.csv")], check=True)
    pb = pd.read_csv(work / "candidates_pb.csv").rename(columns={"run": "pose_id"})
    pb["pb_pass"] = pb.pb_pass.astype(str).str.lower().eq("true")
    cand = cand.merge(pb[["pose_id", "pb_pass", "pb_failed", "max_cc_bond_A", "fe_lig_min_A",
                          "fe_lig_element"]], on="pose_id", how="left", validate="1:1")

    members, fails = [], []
    for b, g in cand.groupby("base_id"):
        ok = g[g.pb_pass.fillna(False)]
        if ok.empty:
            fails.append(b)
            continue
        m = ok.iloc[0]
        members.append({**m.to_dict(), "typeii": b in TYPEII,
                        "n_failed_before": int((g["rank"] < m["rank"]).sum())})
    if fails:
        sys.exit(f"FATAL: no PoseBusters-passing pose in the top {args.topk} for {fails}")

    zp = args.out_dir / f"{args.name}.zip"
    with zipfile.ZipFile(zp, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for m in members:
            text = Path(m["pdb"]).read_text()
            if " LIG " not in text or "LIG0" in text:
                sys.exit(f"FATAL: {m['base_id']}: ligand residue is not exactly LIG in {m['pdb']}")
            zf.writestr(f"{m['base_id']}.pdb", text)
            m["sha256"] = hashlib.sha256(text.encode()).hexdigest()[:16]
    all_members = pd.DataFrame(members)
    # `ckpt` rides along when the scores table declares it (merge_scored.POOL_CKPT). It is what
    # separates two pools built from the SAME cryo-EM template under different checkpoints, so a
    # pick is never reported as "OB0" when a fine-tune produced it.
    has_ckpt = "ckpt" in all_members
    keep = ["base_id", "typeii", "pose_id", "pool", *(["ckpt"] if has_ckpt else []),
            "query_id", "seed", "sample", "rank",
            "n_failed_before", "score", *cols, "srs", "pb_pass", "max_cc_bond_A",
            "fe_lig_min_A", "fe_lig_element", "sha256", "pdb"]
    prov = all_members[keep].sort_values("base_id")
    prov.to_csv(args.out_dir / f"{args.name}.provenance.csv", index=False)
    print(prov.drop(columns=["pdb", "sha256"]).to_string(index=False))
    by_ckpt = f"; by ckpt {prov.ckpt.value_counts().to_dict()}" if has_ckpt else ""
    print(f"\n{len(prov)} members -> {zp}; picks needing a fallback past rank 1: "
          f"{int((prov['rank'] > 1).sum())}; by pool {prov.pool.value_counts().to_dict()}{by_ckpt}")
    return 0 if len(prov) == len(ids) else 1


if __name__ == "__main__":
    sys.exit(main())
