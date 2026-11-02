# Applying the TFG port to Mol-1

Base commit: **6ad1a13**. Three pieces, in this order.

```bash
cp -r tfg                      $MOL1/mol/tfg
cp geometry_featurizer.py      $MOL1/mol/data/core/
git -C $MOL1 apply $PWD/tfg-port.patch
```

`tfg-port.patch` modifies five tracked files (166 insertions):
`mol/config/base.py` · `mol/data/inference/infer_dataloader.py` ·
`mol/data/inference/json_to_feature.py` · `mol/model/generator.py` · `mol/model/mol.py`

## Enabling it

Master switch, off by default: `--sample_diffusion.guidance.enable true`.
The flags reach the parser on `mol/runner/inference.py` (it passes `parse_sys_args()`); the comment
in `base.py` claiming inference takes no CLI overrides applies only to `batch_inference.py`.

## Verify with the witness, never the flag

The engine logs to **stderr**:
- `Guidance is enabled.` — once per inference call
- `TFG last-step energy <Term>: [...]` — once per term per sample

Expect `n_calls` and `n_calls x n_terms x n_samples`. A term reading 0.0 may be zero-weight by
design or may have no features — distinguish before reporting either.

## Two behaviours to know

- **A guided run is not seed-comparable to an unguided one, even at zero weight.** The engine draws
  an extra `randn` per step, so the two RNG streams desynchronise from step 0. Treating an existing
  unguided pool as a controlled baseline confounds guidance with resampling; for clean attribution
  run a third arm with guidance on and every term weight 0.
- **The potentials cannot see a metal cofactor.** The geometry featurizer skips standard and
  metal-containing residues, so a heme contributes zero RDKit terms and the steric/interchain terms
  read 0.0. What is steered is **ligand internal geometry** — not pocket placement, not metal
  coordination.

`PocketPotential` reads undocumented env vars and has chain keys hardcoded from an unrelated
campaign. Leave its `interval` at 0.
