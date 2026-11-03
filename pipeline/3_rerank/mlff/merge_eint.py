"""Join the Orb E_int shards back onto the MM-GBSA shortlist -> the MLFF scorer's ranking table.

Input: the shortlist (build_submission.py --shortlist: the MM-GBSA top-K per ligand after the
blinded-ligand, pool and Fe-gate filters) and the orb_eint_score.py shard CSVs (one row per pose,
failures included). Output: the shortlist rows plus the Orb columns and `mlff_score`, which is what
`SCORER=mlff` ranks on (lower = better), so build_submission re-applies the same gate and the same
PoseBusters fallback to these K candidates as an MM-GBSA build would.

mlff_score = E_int_conserving (never E_int_crossover: that one loads Orb's untested spin gap into
the ranking). The 20 eV/A force gate is applied PER LIGAND: a pose failing it is outside the force
field's domain, so it gets no score -- unless EVERY candidate of that ligand fails, in which case
the ligand is ranked ungated and named here (the fe_gate rule: relax only where forced, and say so).
A ligand with no Orb energy at all is fatal: it cannot be ranked by this scorer.

Printed per ligand: the MM-GBSA rank-1 pose vs the MLFF pick and the E_int margin between them.
Orb's self-noise on byte-identical input is 0.131 eV, so a change of pick by less than that margin
is NOT resolved by this scorer and is counted as such.

Env: sigmadock (pandas). CPU, seconds.
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import pandas as pd

SELF_NOISE_EV = 0.131
ORB_COLS = ["status", "E_int_conserving", "E_int_crossover", "max_force", "force_gate_pass",
            "spin_complex", "n_atoms"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shortlist", type=Path, required=True)
    ap.add_argument("--orb", nargs="+", required=True, help="orb_eint_score shard CSVs / globs")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    sl = pd.read_csv(a.shortlist)
    files = sorted({f for pat in a.orb for f in glob.glob(pat)})
    if not files:
        sys.exit(f"FATAL: no Orb shard matches {a.orb}")
    orb = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    dup = orb.tag[orb.tag.duplicated()]
    if len(dup):
        sys.exit(f"FATAL: pose scored twice across the shards, e.g. {dup.iloc[0]}")
    orb = orb[["tag", *[c for c in ORB_COLS if c in orb]]].rename(
        columns={"tag": "pose_id", "status": "mlff_status"})

    # The shortlist's rank/score are MM-GBSA's; keep the rank under its own name, drop the score so
    # build_submission's re-rank cannot be confused with it.
    df = sl.rename(columns={"rank": "mmgbsa_rank"}).drop(columns=["score"], errors="ignore")
    df = df.merge(orb, on="pose_id", how="left", validate="1:1")
    df["mlff_status"] = df.mlff_status.fillna("not_scored")
    def num(c: str) -> pd.Series:      # a column every row failed before reaching is absent
        return pd.to_numeric(df[c], errors="coerce") if c in df else pd.Series(float("nan"), df.index)
    ok = df.mlff_status.eq("ok")
    gate = ok & num("force_gate_pass").eq(1)
    e = num("E_int_conserving")

    df["mlff_score"] = float("nan")
    ungated, lost = [], []
    for b, g in df.groupby("base_id"):
        if gate[g.index].any():
            idx = g.index[gate[g.index]]
        elif ok[g.index].any():
            idx = g.index[ok[g.index]]
            ungated.append(b)
        else:
            lost.append(b)
            continue
        df.loc[idx, "mlff_score"] = e[idx]

    a.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)
    print(f"{a.out}: {len(df)} shortlisted poses from {len(files)} shard file(s); Orb ok "
          f"{int(ok.sum())}, force gate pass {int(gate.sum())}; status counts "
          f"{df.mlff_status.str[:40].value_counts().to_dict()}")

    # dE_int = E_int(MM-GBSA #1) - E_int(MLFF pick), so a POSITIVE margin is how much E_int the
    # re-rank gains. It goes negative only when the MM-GBSA #1 was itself dropped (force gate or a
    # failed Orb row) and the pick is the best of what remained -- flagged `gated#1` so that sign
    # is never read as the re-rank choosing a worse pose.
    changed, unresolved = 0, 0
    print(f"\n{'ligand':<14}{'n':>3}{'ok':>4}{'gate':>5}  {'MM-GBSA #1':<34}{'MLFF pick':<34}"
          f"{'dE_int(eV)':>10}  note")
    for b, g in df.groupby("base_id"):
        top = g.sort_values("mmgbsa_rank").iloc[0]
        s = g.dropna(subset=["mlff_score"]).sort_values(["mlff_score", "pose_id"])
        if s.empty:
            print(f"{b:<14}{len(g):>3}{int(ok[g.index].sum()):>4}{int(gate[g.index].sum()):>5}  "
                  f"{top.pose_id:<34}{'-- none --':<34}")
            continue
        pick = s.iloc[0]
        margin = (e[top.name] - pick.mlff_score) if pd.notna(e[top.name]) else float("nan")
        moved = pick.pose_id != top.pose_id
        forced = moved and pd.isna(df.mlff_score[top.name])
        changed += moved
        # a forced move is the gate speaking, not a resolved energy difference
        unresolved += bool(moved and not forced and pd.notna(margin) and margin < SELF_NOISE_EV)
        note = "gated#1" if forced else ("unresolved" if moved and pd.notna(margin)
                                         and margin < SELF_NOISE_EV else "")
        print(f"{b:<14}{len(g):>3}{int(ok[g.index].sum()):>4}{int(gate[g.index].sum()):>5}  "
              f"{top.pose_id:<34}{(pick.pose_id if moved else '= same'):<34}"
              f"{margin:>10.3f}  {note}")
    print(f"\nMLFF moves {changed}/{df.base_id.nunique()} picks off the MM-GBSA rank-1; "
          f"{unresolved} of those by less than Orb's {SELF_NOISE_EV} eV self-noise (unresolved) "
          f"-- a `gated#1` row moved because the MM-GBSA #1 has no usable energy, not on margin")
    if ungated:
        print(f"⚠ force gate (20 eV/A) fails EVERY candidate of {ungated} -- ranked ungated there")
    if lost:
        print(f"FATAL: no Orb energy for any candidate of {lost} -- this scorer cannot rank them")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
