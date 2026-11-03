# cofolder

Structure prediction and biomolecular co-folding as three staged pipelines: **template selector →
co-folding → physics re-ranking (two scorers)**. One config file per arm, one entry point, a
validated submission zip at the end.

Built for the CYP3A4 blinded-structure challenge: 20 ligands against a 503-aa P08684 construct with
its heme, co-folded by **Mol-1** (a Protenix fork) with the optional **TFG** chemistry-guidance port.

## The pipeline

```
config/arms/<arm>.env            one arm = one pose pool + its scores + its submission
        │
   1 ───┴─ slurm/run.sh template   select templates → bundles → 20 inputs → PRE-GPU GATE → chunks
        │                          CPU, ~10 min
   2 ────── slurm/run.sh cofold    Mol-1 GPU array (one chunk/task) → PoseBusters QC → poselist
        │                          GPU ~3 h + CPU ~25 min/1000 poses
   3 ────── slurm/run.sh rerank    scorer A MM-GBSA array → merge ─┐
        │                          scorer B MLFF E_int on the top-K ┘  CPU, ~12 CPU-min/pose
        └───── slurm/run.sh select  pick → zip → official validator → ZIP AUDIT → portal diff
                                   CPU, ~20 min
```

`slurm/run.sh status <arm.env>` shows what exists and what is queued. `DRY_RUN=1` prints the
`sbatch` commands without submitting. Run everything from the repo root; `run.sh` itself does no
compute, so it is safe on a login node.

## Run an arm

```bash
cp config/arms/example.env config/arms/myarm.env     # set POOL, CKPT_LABEL, the stage-1 block
$EDITOR config/arms/myarm.env
slurm/run.sh template config/arms/myarm.env          # wait for "STAGE 1 PASSED"
slurm/run.sh cofold   config/arms/myarm.env          # generation + QC, chained
slurm/run.sh rerank   config/arms/myarm.env          # MM-GBSA → merge → MLFF, chained
slurm/run.sh select   config/arms/myarm.env          # submission_<POOL>_<SCORER>.zip
```

Or in one prompt: **`/cofold`** — the skill drives all four stages, stopping at two approval gates
(before the GPU array and before the ~200 CPU-h scoring) unless your prompt pre-approves them.

## Layout

```
config/cluster.env        WHERE things are: roots, environments, partitions. Edit once per cluster.
config/arms/<arm>.env     WHAT an arm is: every knob, with the measured default documented.
config/arms/existing/     the pools already generated — their configs, for re-ranking
config/smoke/             regression + smoke configs (scratch outputs, nothing clobbered)

slurm/run.sh              the only entry point
slurm/lib.sh              config loading, path derivation, the checked sbatch wrapper
slurm/1_template/         template.sbatch
slurm/2_cofold/           generate.sbatch · qc.sbatch
slurm/3_rerank/           mmgbsa.sbatch · merge.sbatch · mlff.sbatch · select.sbatch
slurm/analysis/           report.sbatch — runs one analysis/cofold/ report through Slurm

pipeline/1_template/      select_sim_templates · make_template_json · make_inputs · gate
pipeline/2_cofold/        chunk_inputs · collect_poses · validate_frames
pipeline/3_rerank/mmgbsa/ poselist_from_qc · build_system · mmgbsa_modes · merge_scored
pipeline/3_rerank/mlff/   eint_manifest · eint_prep_shard · prep_* · orb_eint_score · merge_eint
pipeline/3_rerank/select/ build_submission · add_repaired · build_class_hybrid ·
                          build_zip_override · official_validate · lddt_pli_vs_submission
analysis/cofold/          10 post-hoc reports (disposable; run via slurm/analysis/report.sbatch)
mol1_tfg/                 the TFG guidance port: 4 modules + geometry_featurizer + patch
.claude/                  4 stage agents + the /cofold skill
```

Nothing in `pipeline/` reads a hardcoded path: `COFOLD_ROOT` (exported by `lib.sh` from
`config/cluster.env`) overrides the default root everywhere.

## Why one config file per arm

`sbatch --export` splits its value on commas, so `--export=ALL,POOLS=a,b` delivers `POOLS=a` and
exits 0 — that silently halved a candidate pool once. Every job here receives only `COFOLDER`,
`CFG` and `ALLOW_RESUME`, then re-sources `config/cluster.env` + the arm file, so comma-bearing
values (seeds, Fe ladders, partition lists) never cross the boundary. `lib.sh` also **unsets every
arm knob before loading the arm file**, so a stale `SEEDS` or `MAXIMIZE` left in your shell cannot
leak into a job. Change a knob by editing the config, never by exporting it.

## The guards, and what they are for

| guard | why |
|---|---|
| stage 1/2 refuse a pool that already holds poses | a re-run mixes two configurations in one directory and nothing downstream can separate them (`ALLOW_RESUME=1` finishes the *same* config) |
| the pre-GPU template gate | the `.json` template route logs nothing and discards parse errors — the gate is the only witness that templates reached the features |
| `rerank` refuses an existing shard family or a foreign `SCORED` | it would re-spend ~200 CPU-h, or overwrite a table it did not produce |
| `merge_scored --pool-ckpt` | a pool name is one checkpoint, forever: re-declaring it under another label is fatal, so no pick is attributed to the wrong model |
| `select` refuses an existing zip name | a duplicate `NAME` destroys the earlier zip, provenance and audit without warning |
| the from-the-zip audit | re-runs PoseBusters on the shipped bytes; it, not the builder's log, is the claim |
| the checked `sub()` | an unchecked `sbatch` leaves `--dependency=` empty, which Slurm accepts — the dependent job then runs on inputs that do not exist |

Stage dirs must never be renamed while jobs are queued: array tasks resolve their paths at task
start.

## What is already known about this target, so you need not rediscover it

- The template cap is **4, first-N, unsorted** — bundle order *is* the selection; a 5th entry is
  dropped silently. Regression coverage: `comb4` → `[457, 463, 466, 461]`.
- A missing `templatesPath` or MSA path is treated as **absent**, not an error.
- `MSA_MIN`/`MSA_CUT` unset is a *stochastic-depth* arm, not a full-MSA one (rows are redrawn every
  recycle), so both default to 16384.
- Guidance is verified by its **log witness** (`Guidance is enabled.`, `TFG last-step energy`), not
  by the flag; a guided run is not seed-comparable to an unguided one.
- `pb_pass` is **physical plausibility, not accuracy** — the targets have no reference structure.
- The poselist keeps **every** valid pose, which is what sets the ~200 CPU-h scoring bill.
- With `--no-fe-bond` the Fe-donor contact always relaxes outward (min 2.69 Å over 1234 Type-II
  poses), so `FE_RETAIN` is a least-escaped cut, not a coordination bond.
- lDDT between two predictions is **agreement, not accuracy**, and it is not symmetric.
- Orb's self-noise is **0.131 eV** on byte-identical input, above typical between-candidate gaps:
  per-pose only, never an aggregate claim. Egret-1 has no Fe; MACE-POLAR and UBio-MolFM are blocked
  as rankers.

## Ranking strategies, as measured here

| strategy | outcome |
|---|---|
| MM-GBSA over pooled templated arms | the strongest zip; best measured Type-II lDDT-PLI |
| dropping the template-free pool | costs 10.1 kcal/mol over 6 of 20 picks |
| per-ligand similarity templates | best metal coordination and full coverage, but loses the energy contest 19/20 (median +5.19 kcal/mol): a search device, not a ranking gain |
| TFG guidance | 100% pose validity, but a metal-distance floor that removes coordination for the whole coordinating class — usable on the non-coordinating subset only |
| confidence re-rank (`cp_iptm_pl`) | not supported: measured ~random against an oracle here |
| MLFF `E_int` re-rank | rejects net-repulsive top-ranked poses and acts as a distance floor; blind at aggregate margins |
| class split (coordinating / non-coordinating from different arms) | valid and cheap: `SELECT_MODE=classmix` |

## Porting to another cluster

1. **Apply the TFG port** to a Mol-1 checkout (only needed for `GUIDANCE=1`):
   ```bash
   cp -r mol1_tfg/tfg $MOL1/mol/tfg
   cp mol1_tfg/geometry_featurizer.py $MOL1/mol/data/core/
   git -C $MOL1 apply /path/to/mol1_tfg/tfg-port.patch    # base 6ad1a13
   ```
2. **Edit `config/cluster.env`** — roots, the Mol-1 checkout, the five conda prefixes, partitions.
   Environments are invoked by absolute prefix: `conda activate` does not work inside `sbatch`.
   `LAYERNORM_TYPE=torch` where there is no `nvcc`; `kalign` on `PATH` (a startup assertion fires
   whenever templates are on, even on the `.json` route); `TORCH_EXTENSIONS_DIR` off any
   inode-capped filesystem.
3. **Check with `DRY_RUN=1`**, then run `config/smoke/` — stage 1 in both template modes, QC on an
   existing pool, a two-pose MM-GBSA, the MLFF chain, and two selection regressions that rebuild
   shipped zips into a scratch directory and compare provenance.

## Rules that are not negotiable

Everything real goes through `sbatch` — including plotting, tests and inline audits. Announce
absolute input/output paths. Never delete without asking. Never rename a stage directory while jobs
are queued. Nothing here measures correctness while the targets have no reference structure: say
"agreement" when that is what was measured.
