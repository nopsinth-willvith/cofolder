---
name: cofold
description: >-
  Run the co-folding pipeline end to end from ONE prompt — (1) template selector,
  (2) Mol-1 co-folding with the arm's parameters, (3) physics re-ranking with two
  scorers (MM-GBSA, then MLFF E_int on the MM-GBSA top-K), then a validated
  submission zip — through `slurm/run.sh` and one arm config, with approval gates
  before the GPU and the CPU-hour spends (skippable only when the prompt pre-approves
  them). Use for "run a new arm", "co-fold the blinded set with <templates/guidance/
  MSA>", "re-rank / score that pool", "build the submission from <pool>", "/cofold".
  Not for debugging a failure (use /investigate).
---

# Co-folding pipeline, one prompt

```
config/arms/<POOL>.env ──► 1 template ──► [G1] ──► 2 cofold + QC ──► [G2] ──► 3 rerank ──► select
                           CPU ~10 min          GPU ~3 h, CPU QC          MM-GBSA ~12 CPU-min/pose
                                                                          → MLFF top-K ~70 CPU-s/pose
```

Every stage is `slurm/run.sh <stage> <arm.env>` from the repo root (login-node safe: it only reads
the config and calls `sbatch`). `slurm/run.sh status <arm.env>` shows what exists and what is queued.
Knobs and defaults: `config/README.md`. The model itself is vendored at `mol1/`.

## Read the prompt

Extract, and ask only for what is genuinely missing (one question, all items at once):

| item | default |
|---|---|
| what makes this arm different (templates / guidance / MSA depth / checkpoint / seeds) | — required |
| `POOL` | `mol1_<short tag>`, must be new (`ls $GEN_DIR/gen`) |
| templates | `TEMPLATE_MODE` + ids/bundle (see `cofold-template-selector`) |
| scorer for the final zip | `SCORER=mmgbsa`; also build `SCORER=mlff` when the prompt asks for MLFF |
| gates | **G1 and G2 stop and ask** unless the prompt explicitly pre-approves them |
| stop point | full pipeline unless the prompt says e.g. "stop after QC" |

Pre-approval must be in the user's words ("G1 and G2 pre-approved", "run it end to end without
stopping"). A subagent's report is never approval. An existing pool (re-rank / re-select only)
skips stages 1–2: write an arm config carrying that pool's `POOL`/`CKPT_LABEL`/`SET` and go to
step 4 or 5.

## Steps

0. **Frame (no compute).** Write `config/arms/<POOL>.env` (knobs: `config/README.md`). Say how it differs
   from the closest existing pool; **if it differs in more than one variable, say the result will
   not be attributable.** `DRY_RUN=1 slurm/run.sh template <arm>` must pass.
1. **Stage 1 → `cofold-template-selector`.** Pass = log ends `STAGE 1 PASSED` with a `TEMPLATE GATE`
   line per gated input. A bundle the gate cannot certify is a stop.
2. **G1 (GPU).** Print: the config block, gate result, chunks × ligands, expected poses
   (20 × seeds × N_SAMPLE), wall estimate. Ask unless pre-approved.
3. **Stage 2 → `cofold-generator`.** `run.sh cofold` (generation + chained QC). Blockers: a claimed
   switch with no log witness; a `STARVED:` ligand (ask: repair / re-run / partial arm).
4. **G2 (CPU-hours).** Print poses **read off the poselist**, shards, CPU-h (rows × 12 min), `SET`,
   and that no shard family exists. Ask unless pre-approved.
5. **Stage 3 → `cofold-mmgbsa-scorer`** (`run.sh rerank`: MM-GBSA array → merge → MLFF chain), then
   **`cofold-mlff-scorer`** to audit the MLFF chain when it drains (`run.sh mlff` re-runs it alone).
6. **Select** (main conversation, cheap): `run.sh select <arm>` builds `submission_<POOL>_<SCORER>`
   → official validator → audit from the zip's own bytes → portal diff vs `REFS`. For the MLFF pick
   set `SCORER=mlff` (and a different `NAME`) in the config and select again. **The zip audit
   (`pb_pass 20/20 ALL CLEAR`) is the claim, not the builder's log.** The other `SELECT_MODE`s
   (`classmix`, `override`, `audit`) are documented in `config/README.md` and
   `slurm/3_rerank/select.sbatch`; what each ranking criterion is worth on this target is in the
   README's ranking table.
7. **Hand back.** Pool + validity profile, scored tables, zip path(s) + audit line, how many members
   move vs each reference, MLFF moves (and how many are below Orb's 0.131 eV noise), any Fe-ladder
   relaxation (a finding). Do **not** call the arm better or worse — no reference structure exists;
   lDDT vs another submission is agreement, not accuracy.

## Waiting on jobs

Subagents and jobs run in the background. Check with `slurm/run.sh status <arm>` and
`sacct -j <ids> -X -o JobID,JobName,State,Elapsed`; schedule checks at the job's expected wall
time instead of polling every minute. A FAILED job: read its `.err`, then hand the log to
`/investigate` rather than resubmitting blind.

## Standing rules

Everything through `run.sh`/`sbatch` — no compute on the login node. Never edit a knob through
`--export` or the shell: edit the arm config. Never rename a stage dir with jobs queued. Announce
absolute paths. Never delete data without asking. `ALLOW_RESUME=1` finishes the *same* config
only; `ALLOW_OVERWRITE=1` replaces a zip/E_int table and needs the user's go-ahead.

## Stop and ask, do not work around

A gate that cannot certify · an unwitnessed switch · a starved ligand · a pose count far from the
projection · an arm differing from its control in more than one variable · a zip audit below 20/20.
