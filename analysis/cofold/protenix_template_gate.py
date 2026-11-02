"""Pre-GPU gate: how many templates does Protenix ACTUALLY featurise for one system?

Runs runner/inference.py's own config path (run(): two-pass parse_configs + model-specific
defaults + download_inference_cache) with ``main`` swapped for this check, then calls
InferenceDataset.__getitem__ on the first system of --input_json_path -- the exact featurisation
the GPU job performs, reading the template mask out of the returned ``input_feature_dict``.  A
template counts as present when its backbone-frame mask is non-empty.

NOTE it gates ONE system (index 0). With a per-ligand bundle per system that certifies one ligand,
so loop it over every input when the bundles differ (see slurm/cyp3a4_mol1/119_simtmpl4_arm.sbatch).

$GATE_RUNNER_MODULE selects the runner (default runner.inference; Mol-1 needs mol.runner.inference)
and the dataset module is derived from it, since Mol-1 renamed the package protenix -> mol wholesale.
$GATE_DATASET_MODULE overrides that derivation if a fork ever splits the two.

Fails unless EXACTLY $EXPECT_TEMPLATES templates carry atoms, and (if $EXPECT_COVERAGE is set,
space-separated, in input order) each carries >= 95% of its expected residue count.  This is
the check that catches the silent 5 -> 4 truncation (TemplateFeatureAssemblyLine keeps the FIRST
max_templates) and a template that parses to nothing.

STOCK Protenix and Mol-1 both hard-code max_templates=4 on the inference path, so against those
two EXPECT_TEMPLATES above 4 can never pass -- supply at most four templates and expect four. The
Protenix_maxtmpl fork adds --data.template.inference_max_templates (configs_data.py:298, honoured
at infer_dataloader.py:106); slurm/cyp3a4_cryoem_tmpl/25_protenix_multitmpl_gate.sbatch gates that
fork at EXPECT_TEMPLATES=5 and is the reason this script does not assume a cap of 4.

CPU only; the GPU-config step is skipped.  Run from the repo root with runner/inference.py's own
CLI args, e.g.
    EXPECT_TEMPLATES=4 EXPECT_COVERAGE="457 463 466 461" \
    GATE_RUNNER_MODULE=mol.runner.inference GATE_DATASET_MODULE=mol.data.inference.infer_dataloader \
    python protenix_template_gate.py --load_checkpoint_path ... --input_json_path ... \
        --dump_dir ... --seeds 101 --use_template true --device cpu
"""

import importlib
import os
import sys

import numpy as np

sys.path.insert(0, os.getcwd())
# Protenix exposes the runner as `runner.inference`; its Mol-1 fork moved it to
# `mol.runner.inference`. The gate has to import the SAME runner the GPU job will use, or it would
# certify a config path that is not the one that runs, so the module is named by the caller.
ri = importlib.import_module(os.environ.get("GATE_RUNNER_MODULE", "runner.inference"))


def gate(configs) -> None:
    # Mol-1 renamed the whole package protenix -> mol, so the dataset lives at
    # mol.data.inference.infer_dataloader there. Same class, same featurisation. DERIVED from the
    # runner's own package rather than named separately: two env vars with only two valid
    # combinations means setting one and forgetting the other is silently accepted, and the gate
    # would then certify a different featurisation than the GPU job runs -- the exact failure the
    # comment above says must not happen. GATE_DATASET_MODULE stays as an override.
    pkg = ri.__name__.split(".")[0]
    default_ds = ("mol.data.inference.infer_dataloader" if pkg == "mol"
                  else "protenix.data.inference.infer_dataloader")
    ds_mod = os.environ.get("GATE_DATASET_MODULE", default_ds)
    InferenceDataset = importlib.import_module(ds_mod).InferenceDataset

    want = int(os.environ["EXPECT_TEMPLATES"])
    cov = [int(x) for x in os.environ.get("EXPECT_COVERAGE", "").split()]
    data, _, err = InferenceDataset(configs=configs)[0]
    if err:
        raise SystemExit(f"FATAL: featurisation raised:\n{err}")
    # ONE level down, and not searched. process_one builds a flat features_dict and
    # features_dict.update(template_features) lands the template keys in it, but __getitem__ then
    # WRAPS it: infer_dataloader.py:190-191 does `data["input_feature_dict"] = feat`, so no
    # template_* key is ever at top level -- for any bundle, templated or not. Reading the top
    # level instead certified nothing and failed every bundle including the known-good comb4
    # (job 526807 vs the same invocation passing as 504467). Still an explicit path rather than a
    # recursive search, so a same-named key in some other sub-dict cannot be certified by mistake.
    fm = data.get("input_feature_dict", {}).get("template_backbone_frame_mask")
    if fm is None:
        raise SystemExit("FATAL: no template_backbone_frame_mask in the feature dict -- "
                         "templates were not featurised at all (use_template off?)")
    fm = np.asarray(fm)
    per_t = [int(np.diagonal(fm[t]).sum()) for t in range(fm.shape[0])]
    real = [n for n in per_t if n > 0]
    print(f"{data.get('sample_name')}: template slots={fm.shape[0]}, residues with a backbone "
          f"frame per slot={per_t}")
    if len(real) != want:
        raise SystemExit(f"FATAL: {len(real)} templates carry atoms, expected {want}")
    if cov and (len(cov) != len(real) or any(r < 0.95 * c for r, c in zip(real, cov))):
        raise SystemExit(f"FATAL: per-template residues {real} fall short of {cov}")
    print(f"TEMPLATE GATE PASSED: {len(real)} templates reach the features")


ri.main = gate
ri.update_gpu_compatible_configs = lambda c: c
ri.run()
