# mol1/ — vendored Mol-1 inference code

The co-folding model lives **in this repo** so an arm runs code that is committed here, not a
checkout somewhere on a cluster. `config/cluster.env` sets `MOL1=$COFOLDER/mol1`.

## Origin and licence

| | |
|---|---|
| upstream | `https://github.com/molence/Mol-1.git` |
| copied from | `/project/bme-20255165-foldx/Mol-1` at commit **`6ad1a133f919dc79b4f829bdc24fad28b9c4d3f0`** (2026-10-27, "update readme") |
| licence | **Apache-2.0** — `LICENSE` and `THIRD_PARTY_NOTICES.md` are carried here as that licence requires, and the per-file upstream notices are untouched |
| derivation | most of Mol-1 derives from [Protenix](https://github.com/bytedance/Protenix) (ByteDance AML); each file's header names the upstream revision it came from |

Copied: the whole `mol/` package (including `mol/model/layer_norm/kernel/*`, the only non-`.py`
files in it), `pyproject.toml`, `requirements.txt`, `docs/`, `scripts/database/` +
`scripts/gen_ccd_cache.py` + `scripts/check_mol_checkpoint.py`, and the upstream README as
`README.upstream.md`. Left behind: `.venv/`, `.git/`, `data/`, `assets/`, `examples/`, `tests/`,
`scripts/msa/`, the training demo.

The whole package is copied on purpose, not the inference import closure. `mol/data/inference/
__init__.py` imports a module nothing on the path calls, and `mol/model/triangular/layers.py`
resolves two backends through `importlib.import_module`, so a hand-picked file list fails at import
inside an array task hours in, for no meaningful saving.

## What differs from upstream `6ad1a133`: the TFG port

`tfg-port.patch` is the diff against that commit — **5 files, +166/−21** (`mol/config/base.py`,
`mol/model/generator.py`, `mol/model/mol.py`, `mol/data/inference/infer_dataloader.py`,
`mol/data/inference/json_to_feature.py`) — plus two files that are new rather than edited:
`mol/tfg/` (4 modules) and `mol/data/core/geometry_featurizer.py`. All of it is **already applied
in this tree**; the patch is kept only so the delta stays auditable. In the source checkout those
changes were uncommitted, so this copy is their only durable record.

## What is NOT here, and cannot be

**Reference data and weights** — `MOL1_DATA` in `config/cluster.env` (`MOL_ROOT_DIR` at runtime).
Fetch with `scripts/database/download_mol_data.sh`; the inference preflight requires:

| path under `MOL1_DATA` | size | needed for |
|---|---|---|
| `common/components.cif` | 495 MB | always (CCD, read lazily) |
| `common/components.cif.rdkit_mol.pkl` | 143 MB | always (regenerate with `scripts/gen_ccd_cache.py`) |
| `common/clusters-by-entity-40.txt` | 22 MB | existence check only at inference |
| `common/obsolete_release_date.csv` | 0.1 MB | always |
| `common/release_date_cache.json` | 12 MB | templates on |
| `common/obsolete_to_successor.json` | 0.1 MB | templates on |
| `checkpoint/mol-1-base.pt` | 2.3 GB | the stock checkpoint (`CKPT`) |

`MOL1_DATA` must be **writable**: with templates on, template mmCIFs are downloaded into
`$MOL1_DATA/mmcif` on demand.

**The interpreter.** `MOL1_PY` is any Python carrying the dependencies below. It defaults to the
source checkout's venv, which works but keeps one external dependency; to remove it, build one
against this tree:

```bash
uv venv mol1/.venv --python 3.12
uv pip install --python mol1/.venv/bin/python -e mol1            # rewrites the editable mapping
```

⚠ **The editable-install trap.** The default venv carries `__editable__.molence-1.0.0.pth`, whose
finder maps `mol` to the *source* checkout. It is appended after `PathFinder`, so `PYTHONPATH` wins
— but if `PYTHONPATH` were ever lost, `import mol` would resolve there and the arm would silently
run code that is not in this repo. `mol1_env` (in `slurm/lib.sh`) therefore resolves `mol.__file__`
before every stage and refuses anything outside `$MOL1`, printing `mol1: <path>` when it passes.

## The environment that actually produced the campaign's pools

`requirements.txt` does **not** match it (it pins `torch==2.12.0`, `rdkit==2025.9.3`,
`numpy==2.4.4`, …). The installed versions were:

python 3.12.12 · torch 2.13.0 · triton 3.7.1 · cuequivariance{,-torch,-ops-cu13,-ops-torch-cu13}
0.10.0 · nvidia-cuda-runtime 13.0.96 · numpy 2.5.4 · scipy 1.18.1 · pandas 3.1.0 · networkx 3.7 ·
scikit-learn 1.9.1 · rdkit 2026.9.1 · biotite 1.7.1 · biopython 1.88 · gemmi 0.7.5 · modelcif 1.8 ·
pdbeccdutils 2.0.1 · ml_collections 1.1.0 · optree 0.20.0 · pydantic 2.14.1 · requests 2.34.2 ·
click 8.5.0 · tqdm 4.70.1

Everything is a **cu13** build; mixing cu12 operator wheels onto a CUDA-13 torch fails at import
with no useful message. `lmdb` and `orjson` are imported lazily and are not needed.

## Runtime notes that bite

- `LAYERNORM_TYPE=torch` wherever there is no `nvcc` — the fused LayerNorm is JIT-compiled from
  `mol/model/layer_norm/kernel/*` on first call and cannot build without it. `mol1_env` sets it.
- `kalign` must be on `PATH` whenever templates are on: `TemplateFeaturizer.__init__` asserts on it
  even for the `.json` route that never aligns. `mol1_env` prepends `ENV_KALIGN`.
- `TORCH_EXTENSIONS_DIR` must sit off any inode-capped filesystem (not `$HOME` here).
- The Triton TriAttention path writes its autotune cache **inside this package**
  (`mol/model/tri_attention/TriAttentionCache/`, gitignored), and falls back to pure torch
  *silently* if the tree is read-only. Keep `mol1/` writable, or verify the backend in the log.
- `mol/tfg/potentials.py` reads its optional pocket definition from `PROTENIX_POCKET_JSON` /
  `PROTENIX_POCKET_ID`, while `mol/config/base.py` documents `MOL_POCKET_*`. Harmless while
  `PocketPotential.interval` is 0, but the two disagree.
