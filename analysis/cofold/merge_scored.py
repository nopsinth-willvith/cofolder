"""Merge per-pose MM-GBSA scores from several pools into ONE table the ranker can read.

Each scoring run writes its own shards, and only `blinded_scored.csv` carries the `pdb`/`smiles`
columns -- the shard CSVs drop them. So a merge is a join back onto the pose lists, keyed on
pose_id, not a concatenation. Every scored row MUST find its pose file: a row without one would be
silently unrankable (build_submission needs the path to ship and to re-run PoseBusters on), so a
miss is fatal here rather than a short table later.

Pools keep their own `pool` label, which is what tells them apart downstream -- the OB0 arms
(base, tmpl), the cryo-EM-templated fine-tuned arm (fttmpl, template 38FX = the azamulin cryo-EM
structure), and the ESMFold2 arms (esmfold2, esmf2tmpl).

The label alone is not enough to say WHICH MODEL produced a pose, and that distinction became
load-bearing on 2026-10-27, when a parallel session began regenerating STOCK OB0 against the same
cryo-EM templates (pools cryo_t38F*): `fttmpl` and `cryo_t38FX` are the same template 38FX under
two different checkpoints. Reporting a pick from one as the other would attribute it to the wrong
model, so every pool must declare its checkpoint in POOL_CKPT and the merged table carries it as a
`ckpt` column. An undeclared pool is fatal rather than blank.

Env: sigmadock. CPU, seconds.
"""

from __future__ import annotations

import argparse
import csv
import glob
import sys
from collections import Counter
from pathlib import Path

# Same directory, so a plain import works when this runs as a script (the house pattern -- see
# fe_gated_report.py / rank_combined_pool.py importing from build_submission).
from make_poselist import FIELDS

R = Path("/vast.mnt/project/bme-20255165-cofolding-physics")
# What the scorer adds on top of a pose-list row.
SCORE_COLS = ["status", "delta_total", "vdw", "eel", "egb", "esurf", "min_lig_rmsd_A",
              "fe_lig_in_A", "fe_lig_out_A"]
# Built FROM make_poselist.FIELDS rather than retyped: the QC dicts there are designed to be
# extended by one entry, and a hand-copied column list would silently drop any new pose-list
# column a pool brings with it. `ckpt` sits next to `pool` because the two are read together.
_AFTER_POOL = FIELDS.index("pool") + 1
COLS = [*FIELDS[:_AFTER_POOL], "ckpt", *FIELDS[_AFTER_POOL:], *SCORE_COLS]
# Which checkpoint generated each pool. `ob0` is stock of3-ob-2025-06-30-174k; `ftconf5000` is the
# both-heads fine-tune; `esmfold2` is a different model entirely, not an OB0 arm at all.
POOL_CKPT = {"base": "ob0", "tmpl": "ob0", "fttmpl": "ftconf5000",
             "esmfold2": "esmfold2", "esmf2tmpl": "esmfold2",
             **{f"cryo_t38F{t}": "ob0" for t in "UVWXY"}, "cryo5": "ob0",
             "oracle_nopc": "ob0",
             # Mol-1 (Molence, a Protenix fork) -- one checkpoint, mol-1-base, across every arm:
             # the MSA-depth arms and the cryo-EM-templated arms differ only in their inputs.
             **{p: "mol1" for p in ("mol1_full", "mol1_default", "mol1_msa256", "mol1_msa64",
                                    "mol1_reseed", "mol1_tcryo5", "mol1_tcomb4",
                                    # The Fe-bond arms steer through the covalent_bonds JSON key,
                                    # not through weights, so they are still mol-1-base.
                                    "mol1_febond", "mol1_febond_tcomb4", "mol1_febond_tcryo5")},
             **{f"mol1_t38F{t}": "mol1" for t in "UVWXY"},
             # Per-ligand templates chosen by bound-ligand ECFP4 similarity (stage 119). Inputs
             # differ per ligand; the weights are stock mol-1-base, so the label stays "mol1".
             "mol1_tsim4": "mol1",
             # TFG chemistry guidance (stage 110/112). The WEIGHTS are stock mol-1-base, but the
             # sampler is not: the RDKit potentials steer the trajectory, and engine.py:493
             # consumes an extra randn per step, so these poses are not the unguided mol1_* arms
             # under another name and must not be attributed to them.
             **{f"mol1_tfg_{b}": "mol1+tfg" for b in ("tcryo5", "tcomb4", "t38FX", "full")},
             # Fine-tuned on the 99-structure CYP3A4 train split (job 508238, 2000 steps,
             # 1999_ema_0.999.pt). A DIFFERENT checkpoint from every mol1_* pool above, so it gets
             # its own label -- that is the whole point of this table.
             **{f"mol1ft_{b}": "mol1ft" for b in ("full", "tcomb4", "tcryo5", "t38FX")}}


def read_many(patterns: list[str]) -> list[dict]:
    rows, files = [], []
    for pat in patterns:
        hit = sorted(glob.glob(pat))
        if not hit:                          # per-pattern, so a typo'd glob is named in the error
            sys.exit(f"FATAL: no file matches {pat}")
        files += hit
    # Overlapping globs ("mmgbsa_*.csv" alongside "mmgbsa_blinded_shard*.csv") would otherwise read
    # the same shard twice and surface as a duplicate-pose_id FATAL blaming the pose.
    files = sorted(set(files))
    for f in files:
        with open(f, newline="") as fh:
            rows += list(csv.DictReader(fh))
    print(f"  {len(rows)} rows from {len(files)} file(s)")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--scored", nargs="+", required=True, help="score CSVs / shard globs")
    ap.add_argument("--poselists", nargs="+", required=True, help="pose lists supplying pdb+smiles")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    print("scored:")
    scored = read_many(a.scored)
    print("pose lists:")
    # A repeated pose_id across pose lists must not last-win: the two entries can name files from
    # DIFFERENT extractions (cyp3a4_finetune/infer_s2 was overwritten once already), so the pose
    # shipped would not be the pose scored. Same rule as the scored side, which is fatal below.
    lookup, dup = {}, []
    for r in read_many(a.poselists):
        pid = r["pose_id"]
        prev = lookup.get(pid)
        if prev and prev["pdb"] != r["pdb"]:
            dup.append(pid)
        lookup[pid] = r
    if dup:
        sys.exit(f"FATAL: {len(dup)} pose_id(s) appear in more than one pose list with different "
                 f"pose files, e.g. {dup[:5]}")

    seen, out, missing = set(), [], []
    for r in scored:
        pid = r["pose_id"]
        if pid in seen:                      # the same pose scored twice would double its weight
            sys.exit(f"FATAL: duplicate pose_id across the scored inputs: {pid}")
        seen.add(pid)
        src = lookup.get(pid)
        if src is None or not src.get("pdb"):
            missing.append(pid)
            continue
        row = {c: r.get(c) or src.get(c, "") for c in COLS}
        if row["pool"] not in POOL_CKPT:
            sys.exit(f"FATAL: pool {row['pool']!r} does not declare a checkpoint in POOL_CKPT -- "
                     f"add it there, so no pick can be attributed to the wrong model")
        row["ckpt"] = POOL_CKPT[row["pool"]]
        out.append(row)
    if missing:
        sys.exit(f"FATAL: {len(missing)} scored poses have no pose-list entry, e.g. {missing[:5]}")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS); w.writeheader(); w.writerows(out)
    per = Counter(f"{r['pool']}/{r['ckpt']}" for r in out)
    ok_per_lig = Counter(r["base_id"] for r in out if r["status"] == "ok")
    ligs = sorted({r["base_id"] for r in out})
    print(f"\n{a.out}: {len(out)} poses, {ok_per_lig.total()} status=ok, "
          f"{len(ligs)} ligands, pools/ckpt {dict(per)}")
    for b in ligs:                           # Counter gives 0 for a ligand with no scored pose
        if ok_per_lig[b] < 10:
            print(f"  thin: {b} has only {ok_per_lig[b]} scored candidates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
