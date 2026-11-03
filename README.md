# cofolder

Three ways to spend more compute at **inference** time on a co-folding model, and what each one was
worth on a blinded protein–ligand set:

| | contribution | where the compute goes | measured verdict |
|---|---|---|---|
| **1** | **TFG chemistry guidance ported into Mol-1** | steers each diffusion trajectory, training-free on a stock checkpoint (~2.7× wall) | physical validity 86% → **100%**; but it floors the metal distance, so Type-II (metal-coordinating) ligands lose their coordination |
| **2** | **Template selection, gated before the GPU** | conditions the sampler on 4 experimental structures, chosen per arm or per ligand | a **search and coverage** device: best metal geometry and it rescues starved ligands, but wins 1 of 20 slots on energy |
| **3** | **MM-GBSA re-ranking** | ~12 CPU-min per pose, over every valid pose, after sampling | the selection workhorse: it produced the campaign's strongest submission on the one truth measurement available, where model confidence ranked ~random |

Together they are one pipeline — `templates → co-folding → physics re-ranking → a validated
submission` — run as **arms**: one configuration (checkpoint, templates, MSA depth, guidance, seeds)
over the whole target set, yielding one pose pool, one scored table, one submission.

The target set throughout is the **CYP3A4 blinded-structure challenge**: 20 undisclosed ligands
against a 503-aa P08684 construct with its heme, co-folded by **Mol-1** (a Protenix fork). **11 of
the 20 are Type-II** — they coordinate the heme iron — and that split decides almost every number
below.

---

## 1. TFG: training-free chemistry guidance in Mol-1

`mol1_tfg/` back-ports Protenix's Training-Free Guidance into Mol-1: RDKit-derived physical-validity
potentials (bounds matrix, experimental torsions, stereo, chirality) that **steer the diffusion
sampler on a stock checkpoint** — no training, no new weights. 4 modules + `geometry_featurizer.py`
+ `tfg-port.patch` (base `6ad1a13`), applied by `APPLY.md`.

Mol-1 had no guidance machinery at all, and it is lossy about ligand chemistry on the way to the
model (bond orders discarded, the reference conformer embedded without a seed), so the potentials are
built from **the input SMILES**, not from Mol-1's featurised graph. The featuriser skips standard and
metal-containing residues, so **the heme contributes no terms**: only ligand-internal geometry is
steered. And the guided sampler draws one extra random number per step, so a guided run is **not
seed-comparable** to an unguided one, even at zero weight.

Turn it on with `GUIDANCE=1` in the arm config, and **verify it from the log, never from the flag:**
`Guidance is enabled.` once per inference call, `TFG last-step energy` (calls × terms × samples).

**What it bought** (paired per ligand against the matched unguided arm, 1000 poses each):

- PoseBusters validity **86% → 100%**, 18 ligands better / 0 worse / 2 tied, p < 0.001. The guided
  arm produced the campaign's only submission that was 20/20 valid with no repair.
- One ligand that the unguided arm could not solve at all went from 0/50 to 50/50 valid poses.

**What it cost:** a hard floor at 2.783 Å on the Fe–ligand distance. **0 of 1000** Type-II poses
reach the 2.58 Å coordination prior, where 49% of the control's passing poses do — so the selection
Fe ladder had to be relaxed for all 11 Type-II ligands (p = 0.012), and MM-GBSA `delta_total` came
out worse for 14 of 20.

**Verdict:** an inference-time *validity* device, not a ranking gain. Use it on the non-coordinating
subset; do not use it where a metal bond is the answer.

## 2. Template selection, gated before any GPU time

Stage 1 writes the structures an arm sees, in one of four modes:

| `TEMPLATE_MODE` | what it does |
|---|---|
| `none` | template-free control |
| `bundle` | reuse a bundle already built |
| `fixed` | named entries, the same for every ligand — the strongest baseline here |
| `similarity` | per-ligand: the top-4 distinct co-crystal ligands by ECFP4 similarity to the query |

Two facts shape everything: **Mol-1 keeps the first 4 templates, unsorted**, so bundle order *is* the
selection and a 5th entry is dropped silently; and the `.json` template route **logs nothing and
discards parse errors**, so a bundle can fail to load and the run will look normal. Hence the
**pre-GPU gate**: stage 1 ends in `STAGE 1 PASSED` only if the templates are present in the model's
own feature dict (one level down, in `input_feature_dict`), per input. An ungated arm is unverified,
not fine.

The similarity rule is pre-registered in `pipeline/1_template/select_sim_templates.py`: candidate
co-crystal ligands restricted to the pocket, top-4 **distinct compounds** by ECFP4, each resolved to
its best-resolution entry, ties by resolution. **Report its ceiling, not just its picks** — on this
target the best available neighbour scores 0.125–0.340 (median 0.20), so "most similar" means "least
dissimilar".

**What per-ligand similarity templates bought:** 91.0% validity at full 20/20 coverage; Type-II Fe
coordination **11/11 inside the 2.58 Å prior vs 9/11** for the fixed-bundle arm (p = 0.002), with no
ladder relaxation anywhere; validity on the non-Type-II subset p = 0.016; and the same starved ligand
rescued from 0/50 to 78%.

**What it did not buy:** in a four-pool energy contest it takes **1 of 20** slots, with a median gap
of +5.19 kcal/mol to the best pose from another pool. And the gain does **not** track match strength
(strongest-similarity stratum p = 1.000), which points at receptor context rather than ligand
similarity — a shuffled-assignment control would settle it and has not been run.

**Verdict:** templates are the search and coverage knob. Fixed X-ray/cryo-EM bundles win on energy;
per-ligand similarity wins on metal geometry and on ligands no other arm can fold. Different arms,
then, and let selection choose between them.

## 3. MM-GBSA re-ranking, and what else was tried

Stage 3 scores **every** pose that passed QC — no top-N, no score cut — which is what sets the bill:
10–14 CPU-min per pose, ~200 CPU-h per 1000 poses, sharded ≤100 ways.

Per pose: `build_system.py` (ff14SB + Shahrokh IC6 heme + GAFF1/AM1-BCC, residues 25–503, 10/16 Å
pocket, `--no-fe-bond`) → `mmgbsa_modes.py` (explicit-solvent minimisation, igb=2) → `merge_scored.py`,
which attributes the pool to its checkpoint (`--pool-ckpt`) permanently. Type-II ligands then pass an
Fe gate before selection — and because the contact always relaxes *outward* without a bond term (min
2.69 Å over 1234 Type-II poses), that gate is a **least-escaped cut, not a coordination bond**.

| ranking criterion | outcome |
|---|---|
| **MM-GBSA over pooled templated arms** | the strongest zip; best measured Type-II lDDT-PLI in the campaign |
| dropping the template-free pool from the pool set | costs 10.1 kcal/mol over 6 of 20 picks |
| one cryo-EM template pool inside the mix | takes 6 of 20 slots — pool diversity is doing real work |
| model confidence (`cp_iptm_pl`, ipTM) | **not supported**: measured ~random against an oracle here |
| MLFF `E_int` (Orb) on the MM-GBSA top-K | rejects net-repulsive top picks and acts as a distance floor, but self-noise is **0.131 eV** on byte-identical input — above typical between-candidate gaps, so per-pose only |
| relaxing an already-valid pick set | null: 0 poses repaired, worst covalent bond 1.591 → 1.599 Å |
| an Fe pre-filter before scoring | negative: saves 26 of 630 CPU-h and shifts the minimum by ±0.3 Å |
| class split (Type-II and non-Type-II from different arms) | valid and cheap: `SELECT_MODE=classmix` |

Two reporting rules that hold for all of it: `pb_pass` is **physical plausibility, not accuracy** —
the targets are blinded, so nothing here measures correctness — and lDDT between two predictions is
**agreement**, which is also not symmetric. The claim for a submission is the **audit re-run on the
shipped zip's own bytes**, never the builder's log.

---

## Running it

```
config/arms/<arm>.env
   1  run.sh template   select templates → bundles → input JSONs → PRE-GPU GATE → chunks   CPU ~10 min
   2  run.sh cofold     co-folding GPU array (one chunk/task) → PoseBusters QC → poselist   GPU ~3 h
   3  run.sh rerank     MM-GBSA array → merge → MLFF E_int on the top-K                     ~12 CPU-min/pose
      run.sh select     pick → zip → official validator → ZIP AUDIT → portal diff           CPU ~20 min
```

```bash
cp config/arms/example.env config/arms/myarm.env   # POOL, CKPT_LABEL, the stage-1 block, GUIDANCE
slurm/run.sh template config/arms/myarm.env        # then: cofold → rerank → select
```

`run.sh status <arm.env>` shows what exists and what is queued; `DRY_RUN=1` prints the `sbatch`
commands and runs every guard without submitting. `run.sh` itself does no compute: login-node safe.
Each stage chains its own follow-ups (QC behind generation; merge and MLFF behind the MM-GBSA array).

**Four agents own the four stages** — `cofold-template-selector`, `cofold-generator`,
`cofold-mmgbsa-scorer`, `cofold-mlff-scorer` — and each one's job is to produce the log line that
proves its stage did what the config claims: the gate for templates, the per-task echo and guidance
energies for generation, the pose count read off the poselist before ~200 CPU-h are spent. **An
unwitnessed switch is a failed run, not a caveat.** `/cofold` drives all four from one prompt,
stopping before the GPU array and before the scoring spend. None of them decides what a result means;
that comparison stays with the user.

```
config/cluster.env     WHERE things are: roots, environments, partitions. Edit once per cluster.
config/arms/<arm>.env  WHAT an arm is: every knob, with its measured default documented
                       existing/ = the pools already generated · smoke/ = regression configs
slurm/run.sh · lib.sh  the only entry point; config loading and the checked sbatch wrapper
pipeline/<stage>/      durable code: 1_template · 2_cofold · 3_rerank/{mmgbsa,mlff,select}
analysis/cofold/       post-hoc reports (disposable; run via slurm/analysis/report.sbatch)
mol1_tfg/ · .claude/   the TFG port · the 4 stage agents and the /cofold skill
```

Every guard has already cost something. Refused: a pool that already holds poses (a re-run mixes two
configurations inseparably), an existing shard family (~200 CPU-h), a pool re-declared under another
checkpoint label, an existing zip name (a duplicate destroys the earlier zip, its provenance and its
audit). `sbatch --export` splits on commas and exits 0, so jobs receive only `COFOLDER`, `CFG` and
`ALLOW_RESUME` and re-source the configs themselves. Never rename a stage directory while jobs are
queued: array tasks resolve paths at task start.

## Porting to another cluster

1. **Apply the TFG port** to a Mol-1 checkout (only needed for `GUIDANCE=1`):
   ```bash
   cp -r mol1_tfg/tfg $MOL1/mol/tfg
   cp mol1_tfg/geometry_featurizer.py $MOL1/mol/data/core/
   git -C $MOL1 apply /path/to/mol1_tfg/tfg-port.patch    # base 6ad1a13
   ```
2. **Edit `config/cluster.env`** — roots, the Mol-1 checkout, the conda prefixes, partitions.
   Environments are invoked by absolute prefix (`conda activate` does not work inside `sbatch`);
   `LAYERNORM_TYPE=torch` where there is no `nvcc`; `kalign` on `PATH` (asserted whenever templates
   are on, even on the `.json` route); `TORCH_EXTENSIONS_DIR` off any inode-capped filesystem.
3. **`DRY_RUN=1`, then `config/smoke/`** — both template modes, QC on an existing pool, a two-pose
   MM-GBSA, the MLFF chain, and two selection regressions that rebuild shipped zips in a scratch dir.
