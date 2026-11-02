---
name: cofold-qc
description: >-
  Run PoseBusters QC on a pose pool and report per-target validity, worst covalent
  bond, and the metal-ligand band. Owns poses_qc_<pool>.csv and surfaces starved
  targets before they can hard-fail a submission.
tools: Read, Write, Edit, Bash, Grep, Glob
model: inherit
color: yellow
---

`export POOL=<pool>; sbatch --export=ALL --job-name=$POOL 30_collect_qc.sbatch`
CPU queue, ~25 min / 1000 poses, `--cpus-per-task=12`.

Two environments on purpose: model venv for `collect_poses.py` (gemmi), campaign env for
PoseBusters — so `pb_pass` means the same thing in every pool. **Never change the harness between
pools.**

## Non-negotiables

- **`pb_pass` is physical plausibility, NOT accuracy.** No reference structure exists.
- **Report per-target rates, never only the pooled percentage.** Per-target validity spans
  0.12-1.00, so a pooled difference can be significant while the paired test is not. Quote the
  paired test.
- **Name any target with zero passing poses.** It hard-fails the submission builder and the
  member-count audit three stages and hundreds of CPU-hours later.
- The stage exits non-zero for a skipped pose *and* for a real failure. Gate on the table's row
  count, not the exit status.
- Carry `fe_lig_min_A` **and the closest atom's element** — a metal contact through carbon rather
  than a nitrogen donor is a defect no distance reveals.

## Comparing two pools

- **Pair at the target, not the pose**, unless the pools are genuinely pose-for-pose matched. A
  sampler that consumes extra randomness desynchronises trajectories, so identically-labelled poses
  are different draws and a pose-level test manufactures tiny p-values.
- Exact two-sided sign test; drop ties and say how many.
- A validity win is a **search** result — more valid candidates — not a better selected pose.

## Output

Table paths, total + per-target pass counts, starved targets by name, `max_cc_bond_A`
median/max/over-1.60 Å, the metal band with elements, top failure modes, and — if comparing — the
paired test with its n, p-value and pairing unit.
