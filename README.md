# cofolder

Protein–ligand co-folding and pose selection as one staged pipeline: **templates → co-folding →
physics re-ranking → a validated submission**. One config file per arm, one entry point, and a stage
agent that must witness its own stage before the next one spends anything.

An **arm** is one configuration — checkpoint, templates, MSA depth, sampler/guidance, seeds — run
over the whole target set. It yields one pose pool, one scored table, one submission.

## The pipeline

```
config/arms/<arm>.env
   1  run.sh template   select templates → bundles → input JSONs → PRE-GPU GATE → chunks
                        CPU ~10 min
   2  run.sh cofold     co-folding GPU array (one chunk/task) → PoseBusters QC → poselist
                        GPU ~3 h + CPU ~25 min/1000 poses
   3  run.sh rerank     scorer A: MM-GBSA array → merge ─┐  CPU ~12 min/pose
                        scorer B: MLFF E_int on the top-K ┘      ~70 s/pose
      run.sh select     pick → zip → official validator → ZIP AUDIT → portal diff
                        CPU ~20 min
```

```bash
cp config/arms/example.env config/arms/myarm.env            # POOL, CKPT_LABEL, the stage-1 block
slurm/run.sh template config/arms/myarm.env                 # then: cofold → rerank → select
```

Each stage chains its own follow-ups (QC behind generation, merge and MLFF behind the MM-GBSA array).
`run.sh status <arm.env>` shows what exists and what is queued; `DRY_RUN=1` prints the `sbatch`
commands and runs every guard without submitting. `run.sh` does no compute: login-node safe.

## The agents, and what they are for

A stage that silently did the wrong thing is the expensive failure here: a template that never
loaded, or a guidance flag that was never in force, yields a full pool of plausible poses answering
a different question. **So each agent owns one stage, its artifacts, and the log line proving the
stage did what the config claims — and treats an unwitnessed switch as a failed run, not a caveat.**

| agent | stage | owns | its witness |
|---|---|---|---|
| `cofold-template-selector` | 1 | the template block, bundles, input JSONs | the **pre-GPU gate**: templates reaching the model's features, per input |
| `cofold-generator` | 2 | the pose pool, `poses_qc_<POOL>.csv` | per-task config echo, guidance energy lines, per-target pose and validity counts |
| `cofold-mmgbsa-scorer` | 3A | the shards, `scored_<SET>.csv` | the cost **read off the poselist** before spending it; shards present/missing at merge |
| `cofold-mlff-scorer` | 3B | `eint_<POOL>.csv` | the force gate per pose; whether a moved pick beats the force field's self-noise |

**`/cofold`** drives all four stages from one prompt, stopping at the two irreversible spends — **G1**
before the GPU array, **G2** before the hundreds of CPU-hours of scoring — unless the prompt
pre-approves them in the user's own words. The skill holds the sequencing and the gates; the agents
hold the per-stage mechanics.

What no agent does is decide what a result *means*. Calling one arm better than another needs a
comparison the pipeline cannot make for itself, so that stays in the conversation with the user.

## Layout

```
config/cluster.env     WHERE things are: roots, environments, partitions. Edit once per cluster.
config/arms/<arm>.env  WHAT an arm is: every knob, with its measured default documented
config/arms/existing/  the pools already generated — their configs, for re-ranking
config/smoke/          smoke + regression configs (scratch outputs, nothing clobbered)
slurm/run.sh · lib.sh  the only entry point; config loading and the checked sbatch wrapper
slurm/<stage>/         the job scripts, logs/ per stage
pipeline/<stage>/      durable code: 1_template · 2_cofold · 3_rerank/{mmgbsa,mlff,select}
analysis/cofold/       post-hoc reports (disposable; run via slurm/analysis/report.sbatch)
mol1_tfg/              the TFG chemistry-guidance port: 4 modules + featurizer + patch
.claude/               4 stage agents + the /cofold skill
```

Nothing in `pipeline/` hardcodes a path: `COFOLD_ROOT` (exported by `lib.sh`) overrides the root.

## Guards, and why each exists

`sbatch --export` splits its value on commas, so `--export=ALL,POOLS=a,b` delivers `POOLS=a` and
exits 0 — that silently halved a candidate pool once. Jobs therefore receive only `COFOLDER`, `CFG`
and `ALLOW_RESUME` and re-source the configs themselves.

| guard | why |
|---|---|
| stages 1–2 refuse a pool that already holds poses | a re-run mixes two configurations in one directory, inseparably (`ALLOW_RESUME=1` finishes the *same* config) |
| the pre-GPU template gate | the `.json` template route logs nothing and discards parse errors — the gate is the only witness templates arrived |
| `rerank` refuses an existing shard family or foreign `SCORED` | it would re-spend ~200 CPU-h, or overwrite a table it did not produce |
| `merge_scored --pool-ckpt` | a pool name is one checkpoint, forever — so no pick is attributed to the wrong model |
| `select` refuses an existing zip name | a duplicate `NAME` destroys the earlier zip, provenance and audit without warning |
| the from-the-zip audit | re-runs PoseBusters on the shipped bytes; it, not the builder's log, is the claim |

Everything real goes through `sbatch`, including plotting, tests and inline audits. Never rename a
stage directory while jobs are queued (array tasks resolve paths at task start), and never delete
without asking.

## Used for: the CYP3A4 blinded-structure challenge

20 blinded ligands against a 503-aa P08684 construct with its heme, co-folded by **Mol-1** (a
Protenix fork) with the optional **TFG** chemistry guidance in `mol1_tfg/`. 11 of the 20 are
Type-II — they coordinate the heme iron — which splits almost every question below in two.

Each submission is one arm run end to end. `config/arms/existing/` holds the pools built this way
(template-free, three X-ray/cryo-EM bundles, per-ligand similarity templates, one TFG-guided arm);
selection ran MM-GBSA over them, then Orb `E_int` on its top-K.

Measured on this target, so it need not be rediscovered:

- Mol-1's template cap is **4, first-N, unsorted** — bundle order *is* the selection, and a 5th entry
  is dropped silently (regression: `comb4` → `[457, 463, 466, 461]`). A missing template or MSA path
  is **absent**, not an error.
- `MSA_MIN`/`MSA_CUT` unset is a *stochastic-depth* arm, not a full-MSA one (rows are redrawn each
  recycle), so both default to 16384.
- Guidance is verified by its log witness (`Guidance is enabled.`, `TFG last-step energy`), not by
  the flag — and a guided run is not seed-comparable to an unguided one.
- `pb_pass` is **physical plausibility, not accuracy**: the blinded targets have no reference
  structure. lDDT between two predictions is **agreement**, and it is not symmetric.
- The poselist keeps **every** valid pose — that is what sets the ~200 CPU-h scoring bill. With
  `--no-fe-bond` the Fe–donor contact always relaxes outward (min 2.69 Å over 1234 Type-II poses),
  so `FE_RETAIN` is a least-escaped cut, not a coordination bond.
- Orb's self-noise is **0.131 eV** on byte-identical input, above typical between-candidate gaps:
  per-pose only, never an aggregate claim. Egret-1 has no Fe; MACE-POLAR and UBio-MolFM are blocked.

| ranking strategy | outcome |
|---|---|
| MM-GBSA over pooled templated arms | the strongest zip; best measured Type-II lDDT-PLI |
| dropping the template-free pool | costs 10.1 kcal/mol over 6 of 20 picks |
| per-ligand similarity templates | best metal coordination and full coverage, but loses the energy contest 19/20 (median +5.19 kcal/mol): a search device, not a ranking gain |
| TFG guidance | 100% pose validity, but a metal-distance floor that removes coordination for the whole Type-II class — usable on the non-coordinating subset only |
| confidence re-rank (`cp_iptm_pl`) | not supported: measured ~random against an oracle here |
| MLFF `E_int` re-rank | rejects net-repulsive top picks and acts as a distance floor; blind at aggregate margins |
| class split (Type-II / non-Type-II from different arms) | valid and cheap: `SELECT_MODE=classmix` |

## Porting to another cluster

1. **Apply the TFG port** to a Mol-1 checkout (only for `GUIDANCE=1`):
   ```bash
   cp -r mol1_tfg/tfg $MOL1/mol/tfg
   cp mol1_tfg/geometry_featurizer.py $MOL1/mol/data/core/
   git -C $MOL1 apply /path/to/mol1_tfg/tfg-port.patch    # base 6ad1a13
   ```
2. **Edit `config/cluster.env`** — roots, the Mol-1 checkout, the conda prefixes, partitions.
   Environments are invoked by absolute prefix (`conda activate` does not work inside `sbatch`);
   `LAYERNORM_TYPE=torch` where there is no `nvcc`; `kalign` on `PATH` (asserted whenever templates
   are on, even on the `.json` route); `TORCH_EXTENSIONS_DIR` off any inode-capped filesystem.
3. **`DRY_RUN=1`, then `config/smoke/`** — stage 1 in both template modes, QC on an existing pool, a
   two-pose MM-GBSA, the MLFF chain, and two selection regressions that rebuild shipped zips into a
   scratch directory and compare provenance.
