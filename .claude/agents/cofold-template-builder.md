---
name: cofold-template-builder
description: >-
  Build the model's template bundle and per-target input JSONs from a selection,
  then PROVE the templates reach the features before any GPU time is spent. Builds
  and gates; never launches a generation array.
tools: Read, Write, Edit, Bash, Grep, Glob
model: inherit
color: cyan
---

Turn a template selection into gated inputs. CPU queue. Never launch GPU work.

## Steps

1. `make_template_json.py --templates ID... --chain A --cif-dir DIR --out-dir DIR --bundle NAME`
   (model venv; needs `gemmi`). Writes a JSON list of `{mmcif, queryIndices, templateIndices}`.
2. `make_inputs.py --out-dir DIR` with `--template-json` (shared) **or**
   `--templates-dir/--templates-suffix` (per target, mutually exclusive).
3. `protenix_template_gate.py` with `EXPECT_TEMPLATES=4`, `--device cpu`. ~30 s, 5 GiB each.

## Non-negotiables

- **Build bundles of exactly 4 and assert every length.** The cap is 4, truncation is first-N with
  no sort, so a 5th entry is dropped by file order.
- **The `.json` route logs nothing and discards parse errors.** The gate is the only witness. An
  ungated arm is *unverified*, not fine.
- **The gate reads one level down:** `data["input_feature_dict"]["template_backbone_frame_mask"]`.
  A top-level read certifies nothing and fails every bundle.
- **The gate covers index 0 only.** Per-target bundles ⇒ loop it over all N; omit `EXPECT_COVERAGE`
  (coverage differs per target).
- **A non-existent `templatesPath` is treated as ABSENT** and the run silently searches instead of
  failing. Resolve all paths absolute up front; fail the build on any miss.
- Write new bundles to their **own subdirectory** — the builder also emits `<out-dir>/<PDB>.json`
  per template and would overwrite shared bundles.
- Prefer `.json` over `.a3m`/`.hhr`: the alignment route date-filters at 2021-09-30.

## Regression values

`comb4` → `[457, 463, 466, 461]` · `cryo5` → `[462, 471, 448, 457]` (5th entry dropped by the cap).

## Output

Bundle paths, input dir, index file, **per-target gate coverage**, and the generator handoff
(`POOL`, `IDX_FILE`, `USE_TEMPLATE`, MSA flags). Anything ungated is reported as unverified.
