# cofolder

Agentic workflow for structure prediction and biomolecular co-folding.

Migration bundle for the CYP3A4 blinded-structure challenge: the **Mol-1 + TFG** inference port,
the **co-folding and ranking agents**, and the pipeline code they drive.

## What is here

```
.claude/agents/     6 agents — selector, builder, generator, QC, physics scorer, MLFF re-ranker
.claude/skills/     cofold-arm (generate) · cofold-rank (select + compare)
mol1_tfg/           the TFG port: 4 new modules + geometry_featurizer.py + tfg-port.patch
analysis/cofold/    26 scripts — inputs, gates, poselists, selection, comparison, submission
slurm/cofold/       15 sbatch drivers — generation, QC, MM-GBSA, submission, analysis
```

Read `.claude/agents/*.md` for how to run it. Each agent is self-contained and carries its stage's
measured traps; the skills sequence them with approval gates at the two irreversible spends.

## Migration in four steps

**1. Apply the TFG port** to a Mol-1 checkout (`mol1_tfg/`):

```bash
cp -r mol1_tfg/tfg            $MOL1/mol/tfg
cp mol1_tfg/geometry_featurizer.py $MOL1/mol/data/core/
git -C $MOL1 apply /path/to/mol1_tfg/tfg-port.patch   # base 6ad1a13
```

The patch touches `mol/config/base.py`, `mol/data/inference/{infer_dataloader,json_to_feature}.py`,
`mol/model/{generator,mol}.py` (166 insertions). Guidance is off by default; enable with
`--sample_diffusion.guidance.enable true`. **Verify with the witness, not the flag**: the engine
logs `Guidance is enabled.` once per inference call and one `TFG last-step energy` line per term
per sample.

**2. Wire the environment.** Model venv with `gemmi`; `PYTHONPATH=$MOL1`; `MOL_ROOT_DIR=$MOL1/data`;
`LAYERNORM_TYPE=torch` where no `nvcc` exists; `kalign` on `PATH` (a startup assertion fires
whenever templates are enabled, even on the `.json` route); `TORCH_EXTENSIONS_DIR` on a filesystem
without an inode cap. Templates work via the **`.json` route only**.

**3. Point the scripts at the new paths.** `analysis/cofold/*.py` and `slurm/cofold/*.sbatch` carry
absolute paths from the origin cluster; change the `R=` / `ROOT =` assignments at the top of each.
Partitions and the broken-node exclusion are cluster-specific — see `cofold-generator.md`.

**4. Run one arm.** `/cofold-arm`, stopping at both gates. Then `/cofold-rank` to build and compare.

## What the agents already know, so you need not rediscover it

- The template cap is **4, first-N, unsorted** — bundle order *is* the selection.
- The `.json` template route **logs nothing and discards parse errors**; the pre-GPU gate is the
  only witness, and it reads one level down into `input_feature_dict`.
- A missing `templatesPath` or MSA path is treated as **absent**, not an error.
- The poselist keeps **every** valid pose, which sets the ~200 CPU-h scoring bill.
- `sbatch --export` **splits on commas**; a pipeline driver with a hardcoded `SET`/`NAME`
  overwrites the arm it was written for.
- lDDT between two predictions is **agreement, not accuracy**, and it is not symmetric.
- Confidence (`cp_iptm_pl`) measured **~random against an oracle** on this target; MM-GBSA is
  reliable among poses from one method and has misled between methods.
- Orb's **self-noise (0.131 eV) exceeds the margins being compared**; Egret-1 has no Fe; MACE-POLAR
  is blocked as a ranker.

## Ranking strategies, as measured on this target

| strategy | outcome |
|---|---|
| MM-GBSA over pooled templated arms | the strongest zip; best measured Type-II lDDT-PLI |
| dropping the template-free pool | costs 10.1 kcal/mol over 6 of 20 picks |
| per-ligand similarity templates | best metal coordination and full target coverage, but loses the energy contest 19/20 (median +5.19 kcal/mol) |
| TFG guidance | 100% pose validity, but a metal-distance floor that removes coordination for the whole coordinating class — usable on the non-coordinating subset only |
| confidence re-rank | not supported; the metric has no demonstrated skill here |
| MLFF `E_int` re-rank | rejects net-repulsive top-ranked poses and acts as a distance floor; blind at aggregate margins |
| class split (coordinating / non-coordinating from different arms) | valid and cheap; `build_class_hybrid.py` |

## Rules that are not negotiable

Everything real goes through `sbatch` — including plotting, tests and inline audits. Announce
absolute input/output paths. Never delete without asking. Never rename a stage directory while jobs
are queued. Nothing here measures correctness while the targets have no reference structure: say
"agreement" when that is what was measured.
