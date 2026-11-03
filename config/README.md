# Configuring an arm

`config/cluster.env` says **where** things are — roots, environments, partitions. Edit it once per
cluster. The model itself is **vendored in this repo** (`MOL1=$COFOLDER/mol1`); only `MOL1_DATA`
(the CCD caches and the checkpoint) and `MOL1_PY` (an interpreter with those dependencies) point
outside it — see [mol1/VENDORED.md](../mol1/VENDORED.md).

`config/arms/<arm>.env` says **what an arm is**. One arm = one pose pool + its scores + its
submission. Write one file per arm, bash syntax, sourced inside each job (so commas are safe); the
paths from `cluster.env` (`$ROOT`, `$GEN_DIR`, `$GBSA_DIR`, `$MOL1`, …) are available. Every knob
left unset takes the default below, which is the measured reference configuration.

```bash
slurm/run.sh template config/arms/<arm>.env   # 1 template selector  CPU ~10 min
slurm/run.sh cofold   config/arms/<arm>.env   # 2 co-folding + QC    GPU ~3 h
slurm/run.sh rerank   config/arms/<arm>.env   # 3 MM-GBSA → MLFF     CPU ~200 CPU-h
slurm/run.sh select   config/arms/<arm>.env   # 3 pick + zip + audit CPU ~20 min
slurm/run.sh status   config/arms/<arm>.env   # what exists / is queued
```

## Identity — the only two required knobs

| knob | meaning |
|---|---|
| `POOL` | `[A-Za-z0-9_]`, **new per arm**: names `gen/<POOL>`, `poses_qc_<POOL>.csv`, the inputs and the templates |
| `CKPT_LABEL` | the checkpoint every pick is attributed to, permanently. Stock weights keep `mol1`; a changed sampler earns its own label (TFG → `mol1+tfg`); a fine-tune gets its own |

Optional: `SET` (the MM-GBSA shard family, default `$POOL`), `SET_POOLS` (every pool sharing that
`SET`, so QC builds `poselist_<SET>` from all their QC tables instead of the arm that finished last —
must include `POOL`), and `SCORED` (point at an existing merged table to run only `mlff` / `select`).

## Stage 1 — template selector

| knob | default | notes |
|---|---|---|
| `TEMPLATE_MODE` | `none` | `none` · `bundle` · `fixed` · `similarity` |
| `TEMPLATE_BUNDLE` | — | `bundle`: an existing bundle JSON |
| `TEMPLATE_IDS` | — | `fixed`: entries in **rank order** — Mol-1 keeps the first 4, unsorted |
| `TEMPLATE_CIF_DIR` | — | the mmCIF tree `fixed` / `similarity` draw from |
| `TEMPLATE_CHAIN` | `A` | |
| `EXPECT_TEMPLATES` | `4` | the gate **fails** unless exactly this many reach the features |
| `EXPECT_COVERAGE` | — | optional per-template residue floor; shared bundles only |
| `GATE_ALL` | `auto` | every ligand for `similarity`, the first input otherwise |
| `FE_BOND_TABLE` | — | CSV `base_id,donor_atom` declaring a HEM:FE covalent bond. Off-distribution: measured to **over-contract to 1.78–1.96 Å**, and never usable template-free |
| `MSA_A3M` | `$DEFAULT_MSA_A3M` | |

## Stage 2 — co-folding

| knob | default | notes |
|---|---|---|
| `CKPT` | `$MOL1_DATA/checkpoint/mol-1-base.pt` | the stock checkpoint; `CKPT_LABEL` must change with it |
| `SEEDS` | `101…110` | |
| `N_SAMPLE` `N_STEP` `N_CYCLE` | `5` `200` `20` | 10 seeds × 5 samples = 50 poses/ligand |
| `MSA_MIN` `MSA_CUT` | `16384` `16384` | full depth per recycle; set **both empty** for Mol-1's stochastic default, which is a different arm |
| `GUIDANCE` | `0` | `1` = TFG chemistry guidance; `TFG_MU`, `TFG_INNER` tune it |
| `CHUNK_SIZE` | `2` | ligands per array task; >2 OOMs a 22 GiB L4 with templates + full MSA |
| `GEN_TIME` `GEN_MEM` `GEN_THROTTLE` | `03:30:00` `20G` `6` | `GUIDANCE=1` defaults `GEN_TIME` to `12:00:00` (~2.7× the wall) |

## Stage 3 — physics re-ranking and selection

| knob | default | notes |
|---|---|---|
| `PREP_PH` `HEM_PROTONATION` | `7.4` `propionate` | both scorers; a pH variant needs its own `SET` |
| `POSES_PER_SHARD` `GBSA_THROTTLE` | `10` `60` | scorer A runs on **every** valid pose: ~12 CPU-min each |
| `MLFF_TOPK` `MLFF_MODEL` `MLFF_SHARDS` | `3` `orbmol_v2` `6` | scorer B, on the MM-GBSA top-K **as co-folded** (no relaxation); `0` disables. Ranks `E_int_conserving`, 20 eV/Å force-gated |
| `SELECT_MODE` | `rank` | `rank` · `classmix` (`TYPEII_FROM`, `NONTYPEII_FROM`) · `override` (`BASE_ZIP`, `ALT_ZIP`, `FROM_ALT`) · `audit` |
| `SCORER` | `mmgbsa` | `mmgbsa` → `delta_total`, `mlff` → `mlff_score` |
| `RANK_SCORES` | the arm's own table | space-separated tables competing for the pick |
| `RANK_POOLS` | `$POOL` | pools allowed to compete, or `all` |
| `FE_IN` `FE_RETAIN` `FE_LADDER` | `2.6` `3.5` `3.0:4.0,3.5:4.5` | the Type-II Fe gate, in and after minimisation |
| `PB_TOPK` | `5` | PoseBusters fallback depth down the ranking |
| `ADD_REPAIRED` | `0` | append repaired rows from a repair pool; `mmgbsa` only |
| `NAME` `SUB_DIR` `REFS` | `submission_<POOL>_<SCORER>` · `$GEN_DIR/submission` · the campaign's reference zip | `REFS=` (empty) skips the portal diff |

`slurm/lib.sh` holds these defaults in one block (`cofold_load`) and is the authority if this table
and the code ever disagree.
