#!/usr/bin/env python3
"""Per-ligand mix of two submission zips by an EXPLICIT ligand list, not by a score.

Takes every member from --base except the ligands named in --from-alt, which come from --alt.
This is a hand-specified override, so there is no selector and nothing is ranked: the score
column exists only because build_submission.py requires one, and it is a CONSTANT. A constant
also disables the builder's fallback-down-the-ranking behaviour, which is the right thing here
-- both source zips are already pb_pass 20/20, so a member that failed would be a real bug and
should surface rather than be silently swapped out.

The override is checked, not trusted: an id that is not one of the 20, a repeated id, an id
whose two arms ship the SAME pose (a no-op override the caller probably did not intend), and a
missing pose file are each a hard FATAL or a loud warning. The printed table is the record of
which arm each of the 20 slots came from.

Env: sigmadock (stdlib only). CPU, seconds.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_plif_hybrid import prov  # noqa: E402  same provenance reader, same 20-row guard

R = Path("/vast.mnt/project/bme-20255165-cofolding-physics")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default="submission_mol1_relaxtop3_eint",
                    help="zip supplying every ligand NOT named in --from-alt")
    ap.add_argument("--alt", default="submission_mol1_tmpl_mmgbsa",
                    help="zip supplying the ligands named in --from-alt")
    ap.add_argument("--from-alt", nargs="+", required=True, metavar="BASE_ID",
                    help="ligands to take from --alt instead of --base")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    PB, PA = prov(a.base), prov(a.alt)
    want = list(dict.fromkeys(a.from_alt))
    if len(want) != len(a.from_alt):
        sys.exit(f"FATAL: --from-alt repeats an id: {a.from_alt}")
    unknown = [i for i in want if i not in PB]
    if unknown:
        sys.exit(f"FATAL: not among the 20 ligands: {unknown}")

    out, noop = [], []
    for lig, base in PB.items():
        alt = PA[lig]
        take_alt = lig in want
        src = alt if take_alt else base
        if not Path(src["pdb"]).is_file():
            sys.exit(f"FATAL: {lig}: pose file missing: {src['pdb']}")
        if take_alt and alt["pose_id"] == base["pose_id"]:
            noop.append(lig)
        out.append({"base_id": lig, "pose_id": src["pose_id"], "pool": src["pool"],
                    "ckpt": "mol1", "query_id": src["query_id"], "seed": src["seed"],
                    "sample": src["sample"], "srs": src["srs"], "pdb": src["pdb"],
                    "mix": 0.0,
                    "src_arm": "alt" if take_alt else "base",
                    "base_pose": base["pose_id"], "alt_pose": alt["pose_id"],
                    "same_pose": base["pose_id"] == alt["pose_id"]})
    if len(out) != 20:
        sys.exit(f"FATAL: {len(out)} rows, expected 20")

    print(f"base = {a.base}\nalt  = {a.alt}\noverride ({len(want)}): {want}\n")
    print(f"{'ligand':<14}{'from':<6}{'pool':<13}{'=parents?':<11}pose")
    for r in out:
        print(f"{r['base_id']:<14}{r['src_arm'].upper():<6}{r['pool']:<13}"
              f"{'identical' if r['same_pose'] else '-':<11}{r['pose_id']}")
    n_alt = sum(r["src_arm"] == "alt" for r in out)
    same = sum(r["same_pose"] for r in out)
    print(f"\n{20 - n_alt} from base, {n_alt} from alt")
    print(f"  the two parents already agreed on {same}/20 ligands, so this zip differs from")
    print(f"  base on {len(want) - len(noop)} member(s) and from alt on {20 - n_alt - (same - len(noop))}")
    if noop:
        print(f"⚠ {len(noop)} override(s) change NOTHING -- both arms ship the same pose: {noop}")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["base_id", "pose_id", "pool", "ckpt", "query_id", "seed", "sample", "srs", "pdb",
            "mix"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    audit = Path(str(a.out).replace(".csv", "_audit.csv"))
    with open(audit, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"\nwrote {a.out}\nwrote {audit}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
