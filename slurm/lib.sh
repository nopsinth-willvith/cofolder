# slurm/lib.sh -- shared by slurm/run.sh and every stage script. Sourced, never executed.
#
# ONE ARM = ONE CONFIG FILE (config/arms/<arm>.env), layered over config/cluster.env.
# Every job receives exactly COFOLDER (this repo), CFG (the arm file) and ALLOW_RESUME (see sub),
# and re-sources both files through cofold_load. No knob ever crosses `sbatch --export`: that flag
# splits its value on commas, so `--export=ALL,POOLS=a,b` delivers POOLS=a and exits 0 -- the trap
# that silently halved a candidate pool once. Comma-bearing values (SEEDS, FE_LADDER, partitions)
# therefore live in the config files.

# shellcheck disable=SC2034  # every variable set here is consumed by the scripts that source it

# Every arm-config knob. cofold_load UNSETS these before sourcing the arm file, so a value left
# exported in the submitting shell (a stale SEEDS, a MAXIMIZE from an old confidence ranking)
# cannot leak into the job through --export=ALL. Documented knob by knob in config/README.md.
COFOLD_KNOBS="POOL SET CKPT_LABEL SCORED
  TEMPLATE_MODE TEMPLATE_BUNDLE TEMPLATE_IDS TEMPLATE_CIF_DIR TEMPLATE_CHAIN EXPECT_TEMPLATES
  EXPECT_COVERAGE GATE_ALL FE_BOND_TABLE MSA_A3M
  CKPT SEEDS N_SAMPLE N_STEP N_CYCLE MSA_MIN MSA_CUT GUIDANCE TFG_MU TFG_INNER
  CHUNK_SIZE GEN_TIME GEN_MEM GEN_THROTTLE
  SET_POOLS PREP_PH HEM_PROTONATION POSES_PER_SHARD GBSA_THROTTLE MLFF_TOPK MLFF_MODEL MLFF_SHARDS
  SELECT_MODE SCORER RANK_SCORES RANK_POOLS FE_IN FE_RETAIN FE_LADDER PB_TOPK NAME REFS SUB_DIR
  ADD_REPAIRED TYPEII_FROM NONTYPEII_FROM BASE_ZIP ALT_ZIP FROM_ALT"

cofold_load() {
    : "${COFOLDER:?COFOLDER is unset -- submit through slurm/run.sh}"
    : "${CFG:?CFG is unset -- submit through slurm/run.sh <stage> <arm.env>}"
    [ -f "$CFG" ] || { echo "FATAL: no arm config at $CFG" >&2; return 1; }
    # shellcheck disable=SC2086
    unset $COFOLD_KNOBS
    # shellcheck source=../config/cluster.env
    source "$COFOLDER/config/cluster.env" || return 1
    # shellcheck disable=SC1090
    source "$CFG" || return 1

    : "${POOL:?the arm config must set POOL}"
    : "${CKPT_LABEL:?the arm config must set CKPT_LABEL -- the checkpoint every pick is attributed to}"
    case "$POOL" in *[!A-Za-z0-9_]*) echo "FATAL: POOL '$POOL' must be [A-Za-z0-9_]" >&2; return 1;; esac

    # ---- defaults: the measured reference configuration -------------------------------------
    : "${TEMPLATE_MODE:=none}" "${TEMPLATE_CHAIN:=A}" "${EXPECT_TEMPLATES:=4}" "${GATE_ALL:=auto}"
    : "${MSA_A3M:=$DEFAULT_MSA_A3M}"
    : "${CKPT:=$MOL1_DATA/checkpoint/mol-1-base.pt}"
    : "${SEEDS:=101,102,103,104,105,106,107,108,109,110}"
    : "${N_SAMPLE:=5}" "${N_STEP:=200}" "${N_CYCLE:=20}"
    # `=` not `:=`: an arm may set MSA_MIN= MSA_CUT= (empty) to get Mol-1's stochastic default
    : "${MSA_MIN=16384}" "${MSA_CUT=16384}" "${GUIDANCE:=0}"
    : "${CHUNK_SIZE:=2}" "${GEN_THROTTLE:=6}" "${GEN_MEM:=20G}"
    [ "$GUIDANCE" = 1 ] && : "${GEN_TIME:=12:00:00}"      # TFG measured 2.7x the unguided wall
    : "${GEN_TIME:=03:30:00}"
    # Pools that share this SET, i.e. whose QC tables together make poselist_<SET>. Several arms
    # may score as one shard family (the campaign's mol1tmpl = tcryo5 + t38FX + tcomb4); QC would
    # otherwise rewrite the shared poselist from the arm that happened to finish last.
    : "${SET_POOLS:=$POOL}"
    case " $SET_POOLS " in *" $POOL "*) ;;
        *) echo "FATAL: SET_POOLS '$SET_POOLS' does not include this arm's POOL '$POOL'" >&2; return 1 ;; esac
    : "${PREP_PH:=7.4}" "${HEM_PROTONATION:=propionate}"  # both scorers' protonation state
    : "${POSES_PER_SHARD:=10}" "${GBSA_THROTTLE:=60}"     # 10-14 CPU-min/pose -> ~2 h/shard
    : "${MLFF_TOPK:=3}" "${MLFF_MODEL:=orbmol_v2}" "${MLFF_SHARDS:=6}"   # ~70 CPU-s/pose (Orb)
    : "${SELECT_MODE:=rank}" "${SCORER:=mmgbsa}" "${PB_TOPK:=5}" "${ADD_REPAIRED:=0}"
    : "${FE_IN:=2.6}" "${FE_RETAIN:=3.5}" "${FE_LADDER:=3.0:4.0,3.5:4.5}"
    : "${RANK_POOLS:=$POOL}" "${NAME:=submission_${POOL}_${SCORER}}"
    : "${SUB_DIR:=$GEN_DIR/submission}"
    # `=`: an arm may set REFS= (empty) to skip the portal diff
    : "${REFS=$GEN_DIR/submission/submission_mol1_tmpl_mmgbsa.zip}"
    case "$SCORER" in mmgbsa) RANK_COL=delta_total ;; mlff) RANK_COL=mlff_score ;;
        *) echo "FATAL: SCORER must be mmgbsa|mlff, got '$SCORER'" >&2; return 1 ;; esac

    # ---- derived paths: the ONLY place an artifact name is spelled --------------------------
    SET=${SET:-$POOL}                         # MM-GBSA shard family = one pool by default
    IN=$GEN_DIR/inputs_$POOL                  # per-ligand input JSONs + index.txt
    TDIR=$GEN_DIR/templates/$POOL             # this arm's own template bundles (never shared)
    GEN=$GEN_DIR/gen/$POOL                    # Mol-1 --dump_dir
    IDX=$IN/chunks_index_$POOL.txt            # one chunk JSON per array task
    QC=$GEN_DIR/poses_qc_$POOL.csv            # THE table every downstream stage filters on
    POSELIST=$GBSA_DIR/poselist_$SET.csv      # pb_pass poses only; columns 1-3 read positionally
    SHARDS=$GBSA_DIR/mmgbsa_${SET}_shard      # + <i>.csv, and <SHARDS>.nshards (the split)
    SCORED_DEFAULT=$GBSA_DIR/scored_$SET.csv
    SCORED=${SCORED:-$SCORED_DEFAULT}         # merged MM-GBSA table, carries `ckpt`. An arm may
                                              # point it at an existing merged table instead.
    MLFF_DIR=$GBSA_DIR/mlff_$POOL             # shortlist, manifest, prep status, Orb shards
    MLFF_PREP=$SCRATCH/cofold_mlff/prep_ph${PREP_PH}_$HEM_PROTONATION   # npz per pose_id (reused)
    EINT=$GBSA_DIR/eint_$POOL.csv             # MLFF ranking table: the top-K + E_int + mlff_score
    : "${RANK_SCORES:=$([ "$SCORER" = mlff ] && echo "$EINT" || echo "$SCORED")}"
    SUB=$SUB_DIR
    PIPE=$COFOLDER/pipeline
    NSEEDS=$(tr ',' '\n' <<< "$SEEDS" | grep -c .)
    EXPECT_POSES=$((NSEEDS * N_SAMPLE))       # per ligand

    export COFOLD_ROOT=$ROOT PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
    unset PYTHONPATH
    return 0
}

# Mol-1 runtime environment (template gate on CPU, generation on GPU). Returns non-zero if the
# interpreter does not resolve `mol` inside the vendored tree -- always call it as `mol1_env ||
# exit 1`.
mol1_env() {
    export MOL_ROOT_DIR=$MOL1_DATA PYTHONPATH=$MOL1
    export TRITON_CACHE_DIR=$ROOT/.triton_cache
    # else torch builds extensions into $HOME, which sits under a 1M inode cap
    export TORCH_EXTENSIONS_DIR=$ROOT/protenix_data/torch_extensions
    mkdir -p "$TORCH_EXTENSIONS_DIR"
    # no nvcc on this cluster, so the fused LayerNorm can never build; cuEquivariance is untouched
    export LAYERNORM_TYPE=torch
    # find_kalign() asserts in TemplateFeaturizer.__init__ whenever templates are on, even for the
    # .json route that never aligns -- a startup assertion, so kalign is simply put on PATH
    export PATH=$ENV_KALIGN/bin:$PATH

    # MOL1_PY may carry an editable install of ANOTHER Mol-1 checkout: a .pth whose finder maps
    # `mol` to that tree. It is appended after PathFinder, so PYTHONPATH wins -- but if PYTHONPATH
    # is ever lost, `import mol` resolves there instead and the arm silently runs code that is not
    # in this repo, with nothing in the log to say so. So resolve the package and refuse anything
    # outside $MOL1. (mol/__init__.py imports only mol/version.py: no torch, no GPU.)
    local got
    got=$("$MOL1_PY" -c 'import mol,sys; sys.stdout.write(mol.__file__)' 2>&1) || {
        echo "FATAL: $MOL1_PY cannot import mol from $MOL1:" >&2; echo "$got" >&2; return 1; }
    case "$got" in
        "$MOL1"/*) echo "mol1: $got" ;;
        *) echo "FATAL: 'import mol' resolved to $got -- OUTSIDE the vendored tree $MOL1." >&2
           echo "       Build MOL1_PY's venv from mol1/pyproject.toml, or remove the editable" >&2
           echo "       __editable__*molence*.pth from it (mol1/VENDORED.md)." >&2; return 1 ;;
    esac
}

# sub ARGS... : sbatch --parsable with BOTH the exit status and the id checked; echoes the id.
# An unchecked sbatch leaves --dependency= empty, which Slurm accepts, so the dependent job starts
# at once on inputs that do not exist yet. `${id%%;*}` strips a federated "<jobid>;<cluster>".
# DRY_RUN=1 prints the command and returns a fake id instead of submitting.
#
# --export is EXPLICIT: this cluster's cli_filter forces --export=NONE (the job starts from the
# login environment, nothing inherited), so the job gets exactly COFOLDER, CFG and ALLOW_RESUME.
# None of them may hold a comma -- sbatch would split the value there.
sub() {
    local id exp="COFOLDER=$COFOLDER,CFG=$CFG,ALLOW_RESUME=${ALLOW_RESUME:-0}"
    case "$COFOLDER$CFG" in *,*) echo "FATAL: a comma in '$COFOLDER' or '$CFG'" >&2; return 1;; esac
    if [ "${DRY_RUN:-0}" = 1 ]; then
        echo "DRY_RUN: sbatch --export=$exp $*" >&2
        echo "DRY$RANDOM"; return 0
    fi
    id=$(sbatch --parsable --export="$exp" "$@") || { echo "SUBMIT FAILED: $*" >&2; return 1; }
    id=${id%%;*}
    [ -n "$id" ] || { echo "SUBMIT exited 0 but printed no job id: $*" >&2; return 1; }
    echo "$id"
}

# Rows of a CSV that has one record per line (poselists, index files) -- never for QC tables,
# whose pb_failed column can hold quoted commas.
nrows() { [ -s "$1" ] && echo $(($(wc -l < "$1") - 1)) || echo 0; }

# Count finished prediction CIFs in a pool.
npose() { find "$1" -name '*_sample_*.cif' 2>/dev/null | wc -l; }
