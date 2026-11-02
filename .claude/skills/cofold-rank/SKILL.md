---
name: cofold-rank
description: >-
  Turn scored pose pools into a ranked, validated submission and decide whether one
  candidate set beats another — MM-GBSA selection, the metal gate, confidence,
  MLFF, class splits and per-member overrides. Use for "build the submission",
  "rank these poses", "is this zip better", "which picks are suspect".
---

# Rank and build

## Build a submission

`50_submission.sbatch`: `NAME SCORES COLS OUTDIR REFS` (+ optional `POOLS` `FE_IN` `FE_RETAIN`
`FE_LADDER` `MAXIMIZE`). Four phases: select + PoseBusters + zip → official validator → audit
**from the zip's own bytes** → portal diff vs each reference. The re-audit is the claim.

- **`unset MAXIMIZE SEEDS` before an energy ranking.** `delta_total` is an energy; a leaked
  `MAXIMIZE` ships the N *worst* binders and nothing complains. `SEEDS` FATALs with several pools.
- **Never pass a comma-bearing value through `--export`** — only the first item arrives, and the
  job exits 0. This has bitten twice.
- **`POOLS` silently drops repair rows** (they belong to no named pool) — with a starved target
  that becomes a hard FATAL. Name the repair pool or omit `POOLS`.
- **A duplicate `NAME` destroys the earlier zip, provenance and audit**, no warning. Gate on it.
- Report fallback count, pool/ckpt split, and **any metal-ladder relaxation — that is a finding**.

## Ranking criteria, by measured reliability

| criterion | use |
|---|---|
| **MM-GBSA `delta_total`** | primary. Reliable for choosing *among* poses from one method; **has been misleading for choosing between methods** |
| **metal gate** (`FE_IN`/`FE_RETAIN`) | hard physical filter on the coordinating class. `fe_lig_out_A` is **post-minimisation**, a least-escaped cut, not a bond |
| **confidence** (`cp_iptm_pl`) | the per-pair number. **Measured ~random vs an oracle here** — report, do not rank. `plddt`/`ptm` are flat; `gpde` is lower-better; `ranking_score`'s clash term is polymer-only |
| **MLFF `E_int`** | `cofold-mlff-reranker` only; per-pose, never aggregate |

Decompose before explaining: `delta_total` cannot separate shape from electrostatics; `vdw`/`eel`/
`egb`/`esurf` can.

## Comparing two submissions

- **lDDT against another prediction is AGREEMENT, not accuracy.** A reference's measured score does
  **not** transfer to what agrees with it.
- **Not symmetric** — run both directions, report both.
- **Read lDDT CSVs with `csv.reader`, never `awk -F,`**: assignment columns hold tuple reprs with
  commas, which shifts fields by two and **silently swaps PLI for LP**.
- **Name local variables off the env-var namespace.** A local `REF`/`MODEL` holding a zip path is
  clobbered by the stage's own `export`, so later pairs compare the wrong zips.
- Distinguish a **different binding mode** from a **reorientation**: compare RMSD against centroid
  shift. Centroid ≪ RMSD = rotation in place; comparable = the ligand moved.
- **Stratify by class and effect size.** A pooled p-value hides a real effect in one class and
  manufactures one from the other.

## Mixing

- `build_class_hybrid.py --typeii-from A --nontypeii-from B --name N --sub-dir D --audit-dir D`
  — byte copies, each member gated on its source audit. **Re-validate and re-audit after.**
- `build_zip_override.py --base A --alt B --from-alt ID...` — per-member override by explicit list.
- Verify the split landed: 0 members moving in the class each parent supplied, movement in the other.

## Report honestly

Nulls as nulls with the p-value. Separate "suggestive but underpowered" from "established". Name
hypotheses you **refuted**, including your own. Nothing here measures correctness when the targets
have no reference structure.
