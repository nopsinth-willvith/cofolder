---
name: cofold-arm
description: >-
  Run one co-folding arm end to end — select templates, build and gate them,
  generate the pose pool, QC it, and score it with MM-GBSA — with approval gates
  at the two irreversible spends. Use for "run a new arm", "generate poses for",
  "run the pipeline on", "score that pool". Not for building or comparing
  submissions (use /cofold-rank) or debugging a failure.
---

# Run a co-folding arm

An arm = one pose pool + its scored table: one configuration of checkpoint, templates, MSA depth
and sampler, over the whole target set, QC'd and energy-scored.

## Two approval gates — stop and ask; a subagent's report is not approval

| gate | spend | print before asking |
|---|---|---|
| **G1** before generation | GPU array, ~3 h | configuration; per-target gate coverage; pool name; how this arm differs from existing pools |
| **G2** before scoring | CPU array, ~200 CPU-h | pose count **read off the poselist**; projected CPU-hours; `SET`; clobber check |

## Sequence

0. **Frame it** (no compute). State what makes this arm different: checkpoint, templates, MSA
   depth, sampler, seeds. If it differs in **more than one**, say so — the result will not be
   attributable. Check what already exists; do not regenerate a pool.
1. **`cofold-template-selector`** → manifest + bundle list. Relay a weak selection ceiling
   **before G1**; it changes what the arm can show.
2. **`cofold-template-builder`** → bundles, inputs, per-target gate. **Do not proceed past a
   bundle the gate could not certify.** "The file exists" is not certification.
3. **GATE G1** — ask. Wait.
4. **`cofold-generator`** → array + witness audit. An unwitnessed switch is a **blocker**, not a
   caveat: the pool may not be the arm you think it is.
5. **`cofold-qc`** → per-target validity, starved targets, metal band. A starved target is
   stop-and-report: the user decides graft / re-run / partial arm.
6. **GATE G2** — ask. Wait.
7. **`cofold-physics-scorer`** → merged scored table, `POOL_CKPT` entry, shard failures.
8. **Hand back** (no compute). Report pool, validity profile, scored table. Do **not** build a
   submission and do **not** declare the arm better or worse — say what would establish that.

## Standing rules

Everything through `sbatch`. CPU → `tue.cpu1.q,tue.cpu2.q,tue.cpu3.q`. GPU →
`bme.gpumolml.q,bme.gpu2.q,tue.gpu3.q` with `--exclude=tue-gpua001,tue-gpub001` and **`tue.gpu2.q`
absent**. Never rename a stage dir with jobs queued. Never pass a comma-bearing value through
`--export`. Never drive a pipeline driver whose `SET`/`NAME` are hardcoded for a new pool — it
overwrites the arm it was written for; author a sibling. Announce absolute paths; never delete
without asking.

## Stop and ask, do not work around

A bundle the gate cannot certify · a switch with no log witness · a target with zero valid poses ·
a pose count materially different from projection · an arm differing from its control in more than
one variable.
