#!/bin/bash
#SBATCH --job-name=cofold_tmpl
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=slurm/1_template/logs/%x_%j.out
#SBATCH --error=slurm/1_template/logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:30:00
#SBATCH --mem=12G
set -uo pipefail
# STAGE 1 -- template selection -> bundles -> per-ligand inputs -> pre-GPU gate -> chunks.
# Submitted by `slurm/run.sh template <arm.env>`. Builds and gates only; never launches GPU work.
#
# TEMPLATE_MODE (arm config):
#   none        template-free arm
#   bundle      reuse an existing bundle JSON for every ligand   (TEMPLATE_BUNDLE=<path>)
#   fixed       build one bundle from named entries, same for all (TEMPLATE_IDS="5VCC 6MA8 ...")
#   similarity  per-ligand top-4 distinct co-crystal ligands by ECFP4 (select_sim_templates.py)
#
# Mol-1 keeps the FIRST 4 templates of a bundle, UNSORTED, so bundle order IS the selection. The
# .json route logs nothing and discards parse errors, so step 4's gate is the ONLY witness that
# templates reach the features -- an ungated templated arm is unverified, not fine.
source "${COFOLDER:?submit via slurm/run.sh}/slurm/lib.sh" && cofold_load || exit 1
mol1_env || exit 1
P1=$PIPE/1_template

echo "######## 0 gate: an unmixed pool"
NP=$(npose "$GEN")
if [ "$NP" -gt 0 ] && [ "${ALLOW_RESUME:-0}" != 1 ]; then
    echo "FATAL: $GEN already holds $NP poses -- re-building would MIX two configurations in one" >&2
    echo "       pool. Pick a new POOL, or ALLOW_RESUME=1 to re-chunk only the unfinished ligands." >&2
    exit 3
fi

echo "######## 1 select templates (TEMPLATE_MODE=$TEMPLATE_MODE)"
TOPT=()
assert_len() {   # assert_len <bundle.json> <max>: the cap is first-N, so length is the selection
    local n
    n=$("$MOL1_PY" -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" "$1") || return 1
    echo "  $(basename "$1"): $n template(s)"
    [ "$n" -ge 1 ] && [ "$n" -le "$2" ] || { echo "FATAL: $1 holds $n templates (allowed 1-$2)" >&2; return 1; }
}
case "$TEMPLATE_MODE" in
  none)
    echo "  template-free" ;;
  bundle)
    [ -f "${TEMPLATE_BUNDLE:-}" ] || { echo "FATAL: TEMPLATE_BUNDLE '${TEMPLATE_BUNDLE:-}' is not a file" >&2; exit 3; }
    # an existing bundle may hold 5 (cryo5: the 5th is dropped by the cap, knowingly) -- warn only
    assert_len "$TEMPLATE_BUNDLE" 99 || exit 1
    TOPT=(--template-json "$(readlink -f "$TEMPLATE_BUNDLE")") ;;
  fixed)
    : "${TEMPLATE_IDS:?TEMPLATE_MODE=fixed needs TEMPLATE_IDS}" "${TEMPLATE_CIF_DIR:?fixed needs TEMPLATE_CIF_DIR}"
    mkdir -p "$TDIR"   # its own dir: the builder also writes <PDB>.json per entry as a side effect
    # shellcheck disable=SC2086
    "$MOL1_PY" "$P1/make_template_json.py" --templates $TEMPLATE_IDS --chain "$TEMPLATE_CHAIN" \
        --cif-dir "$TEMPLATE_CIF_DIR" --out-dir "$TDIR" --bundle "$POOL" || exit 1
    assert_len "$TDIR/$POOL.json" "$EXPECT_TEMPLATES" || exit 1
    TOPT=(--template-json "$TDIR/$POOL.json") ;;
  similarity)
    : "${TEMPLATE_CIF_DIR:=$ROOT/data/cyp_structure_track/cif}"
    mkdir -p "$TDIR"
    "$ENV_QC/bin/python" "$P1/select_sim_templates.py" --n-templates "$EXPECT_TEMPLATES" \
        --out-manifest "$TDIR/manifest.csv" --out-list "$TDIR/bundles.txt" || exit 1
    NL=$(grep -c . "$TDIR/bundles.txt")
    [ "$NL" -eq 20 ] || { echo "FATAL: expected 20 bundle lines, got $NL" >&2; exit 1; }
    while read -r LIG PDBS; do
        [ -n "$LIG" ] || continue
        # shellcheck disable=SC2086
        "$MOL1_PY" "$P1/make_template_json.py" --templates $PDBS --chain "$TEMPLATE_CHAIN" \
            --cif-dir "$TEMPLATE_CIF_DIR" --out-dir "$TDIR" --bundle "$LIG" > /dev/null \
            || { echo "FATAL: bundle $LIG" >&2; exit 1; }
        # exactly EXPECT_TEMPLATES: a short bundle changes the selection silently
        N=$("$MOL1_PY" -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" "$TDIR/$LIG.json") || exit 1
        [ "$N" -eq "$EXPECT_TEMPLATES" ] || { echo "FATAL: $LIG.json holds $N, expected $EXPECT_TEMPLATES" >&2; exit 1; }
    done < "$TDIR/bundles.txt"
    echo "  20 per-ligand bundles of exactly $EXPECT_TEMPLATES -> $TDIR (manifest: $TDIR/manifest.csv)"
    TOPT=(--templates-dir "$TDIR" --templates-suffix .json) ;;
  *)
    echo "FATAL: TEMPLATE_MODE must be none|bundle|fixed|similarity, got '$TEMPLATE_MODE'" >&2; exit 3 ;;
esac

echo "######## 2 per-ligand inputs (same 503-aa construct, MSA, HEM and SMILES as every arm)"
"$MOL1_PY" "$P1/make_inputs.py" --out-dir "$IN" --msa "$MSA_A3M" ${TOPT[@]+"${TOPT[@]}"} \
    ${FE_BOND_TABLE:+--fe-bond-table "$FE_BOND_TABLE"} || exit 1
NJ=$(grep -c . "$IN/index.txt")
[ "$NJ" -eq 20 ] || { echo "FATAL: built $NJ inputs, expected 20" >&2; exit 1; }

echo "######## 3 pre-GPU template gate"
if [ "$TEMPLATE_MODE" = none ]; then
    echo "  skipped: template-free arm"
else
    # per-ligand bundles => gate every ligand (index 0 of each file); a shared bundle => one.
    # EXPECT_COVERAGE only makes sense for a shared bundle (coverage differs per ligand otherwise).
    ALL=$GATE_ALL; [ "$ALL" = auto ] && { [ "$TEMPLATE_MODE" = similarity ] && ALL=1 || ALL=0; }
    [ "$ALL" = 1 ] && GLIST=$(cat "$IN/index.txt") || GLIST=$(head -1 "$IN/index.txt")
    COV=""; [ "$TEMPLATE_MODE" != similarity ] && COV=${EXPECT_COVERAGE:-}
    NG=0
    for JSONP in $GLIST; do
        ( cd "$MOL1" && EXPECT_TEMPLATES=$EXPECT_TEMPLATES EXPECT_COVERAGE="$COV" \
          GATE_RUNNER_MODULE=mol.runner.inference \
          "$MOL1_PY" "$P1/protenix_template_gate.py" --load_checkpoint_path "$CKPT" \
            --input_json_path "$JSONP" --dump_dir "$IN/gate_dump" --seeds "${SEEDS%%,*}" \
            --use_msa true --use_rna_msa false --use_template true --device cpu 2>&1 \
          | grep -E "template slots|TEMPLATE GATE|FATAL" ) \
          || { echo "GATE FAILED on $(basename "$JSONP") -- do not launch this arm" >&2; exit 1; }
        NG=$((NG + 1))
    done
    echo "  gate PASSED on $NG input(s): $EXPECT_TEMPLATES templates reach the features"
fi

echo "######## 4 chunk ($CHUNK_SIZE ligands/task; >2 OOMs a 22 GiB L4 with templates + full MSA)"
"$MOL1_PY" "$PIPE/2_cofold/chunk_inputs.py" --inputs-dir "$IN" --pool "$POOL" \
    --chunk-size "$CHUNK_SIZE" --expect "$EXPECT_POSES" --gen-root "$GEN_DIR/gen" || exit 1
echo
echo "handoff -> slurm/run.sh cofold $CFG"
echo "  pool=$POOL ckpt=$CKPT_LABEL templates=$TEMPLATE_MODE chunks=$(grep -c . "$IDX") ($IDX)"
echo "STAGE 1 PASSED"
