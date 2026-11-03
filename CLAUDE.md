# cofolder — working rules

Three staged pipelines: **1 template selector → 2 co-folding → 3 physics re-ranking (MM-GBSA, then
MLFF E_int) → select**. Read [README.md](README.md) for the layout; `/cofold` runs the whole thing
from one prompt.

## How work enters this repo

- **One entry point:** `slurm/run.sh <template|cofold|rerank|mlff|select|status> <arm.env>`, from
  the repo root. Never `sbatch` a stage script by hand; never run a `pipeline/` script inline.
- **One config per arm:** `config/arms/<arm>.env`. To change a knob, **edit the config** — not the
  shell, not `--export`, not the sbatch header. `config/arms/example.env` documents every knob with
  its measured default; `config/cluster.env` holds paths, environments and partitions.
- **`DRY_RUN=1` first.** It prints the `sbatch` commands and runs every guard without submitting.
- A new arm needs a **new `POOL`**. Reusing one is refused: a re-run would mix two configurations in
  one pose directory and nothing downstream could separate them.

## Compute

Everything real goes through `sbatch`, including plotting, pytest, grep sweeps and inline audits.
`run.sh` and `lib.sh` are the only things that may run on a login node (they read config, count
lines, call `sbatch`). A one-off report goes through `slurm/analysis/report.sbatch <script.py>`.

CPU → `tue.cpu1.q,tue.cpu2.q,tue.cpu3.q`. GPU → `GPU_PARTITION` from `config/cluster.env`, which
**omits `tue.gpu2.q`** (its only node is the broken `tue-gpub001`: naming it beside the mandatory
`--exclude=tue-gpua001,tue-gpub001` leaves zero eligible nodes and the array parks at
BadConstraints forever).

Environments are invoked by absolute prefix — `conda activate` does not work inside `sbatch`.

## Adding code

Match the staged layout: durable code in `pipeline/<stage>/`, job scripts in `slurm/<stage>/`,
disposable reports in `analysis/cofold/`. Before writing anything, read the stage you are extending;
add the smallest change that fits. A `pipeline/` script must take its root from `COFOLD_ROOT`
(`lib.sh` exports it) rather than hardcoding a path, and must not import from outside the repo —
code borrowed from the parent project is **vendored** into the stage that uses it, with a first-line
comment naming its origin.

Values that may contain a comma (seeds, Fe ladders, partition lists) live in config files and are
sourced inside the job. `sbatch --export` splits on commas and exits 0, so a knob passed that way
arrives truncated. Jobs receive only `COFOLDER`, `CFG` and `ALLOW_RESUME`.

## Reporting results

- `pb_pass` is **physical plausibility, not accuracy**. The blinded targets have no reference
  structure, so nothing here measures correctness — say "agreement" when that is what was measured.
- lDDT against another prediction is agreement, and it is **not symmetric**: run both directions.
- Report **per-target** rates, never only the pooled percentage; name every starved ligand.
- Quote the **log witness** for every switch an arm claims (templates: the stage-1 gate; guidance:
  `Guidance is enabled.`). An unwitnessed switch is unverified, not fine.
- The **from-the-zip audit** (`pb_pass 20/20 ALL CLEAR`) is the claim, not the builder's own log.
- A Type-II Fe-ladder relaxation is a finding — report which ligands needed it.
- MLFF margins below Orb's **0.131 eV** self-noise are unresolved, not a winner.

## Verify before claiming done

`bash -n` and `shellcheck -x` on changed shell; `DRY_RUN=1 slurm/run.sh <stage>` on an affected
config; then the relevant `config/smoke/` job through `sbatch`. The two selection regressions
(`config/smoke/regress_*.env`) rebuild shipped submissions into a scratch directory and must
reproduce their provenance pose for pose. After a nontrivial change run `/code-review` and
`/simplify` on the diff.

## Never

Delete data without asking. Rename a stage directory while jobs are queued (array tasks resolve
paths at task start). Overwrite a zip, a scored table or an E_int table that something already
shipped from — `ALLOW_OVERWRITE=1` exists for that and needs the user's go-ahead. Commit or push
unless asked.
