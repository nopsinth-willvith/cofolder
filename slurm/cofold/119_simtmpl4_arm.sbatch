#!/bin/bash
#SBATCH --job-name=mol1_119_sim4
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/119_sim4_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/119_sim4_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:30:00
#SBATCH --mem=12G
set -uo pipefail
# Stage 119: PER-LIGAND templates chosen by BOUND-LIGAND chemical similarity (ECFP4).
#
# Every other template arm hands all 20 blinded ligands the SAME bundle, chosen for pocket-state
# diversity. This one gives each ligand the 4 corpus structures whose own co-crystal ligand is
# chemically closest to it. The rule is pre-registered in select_sim_templates.py's docstring.
#
# ⚠ THE HYPOTHESIS HAS A MEASURED REJECTION IN A SIBLING MODEL. On OF3/OpenBind-0 the same idea
# ran as `simtmpl50` and was the WORST of three regimes: 69.4% PB-valid vs 81.8% for a fixed
# pocket-state-diverse bundle, and tmpl-vs-simtmpl was the ONLY pair significant after Holm
# (median delta +0.090, p_holm 0.017). The stated cause is a ceiling that belongs to the CORPUS,
# not the method: best-neighbour ECFP4 is ~0.125-0.340 (median ~0.20), because the blinded
# ligands are small (MW median 329) and the corpus is dominated by large ones (ritonavir 534).
# This arm re-tests it on a DIFFERENT model -- Mol-1, whose templated arms produced the
# campaign's best measured zip -- so a null is informative and a win would be surprising.
#
# EXACTLY FOUR templates per ligand. Mol-1 hard-codes TemplateFeatureAssemblyLine(max_templates=4)
# at mol/data/template/template_featurizer.py:774 and truncates FIRST-N with no sort, so a fifth
# entry would be dropped by file order. The .json route logs NOTHING (the "Found N templates"
# line is on the .a3m/.hhr branch only), so step 3's gate is the ONLY witness that 4 bundles'
# worth of templates reached the features -- and it must run for all 20, not just the first,
# because every ligand now has a DIFFERENT bundle.
R=/vast.mnt/project/bme-20255165-cofolding-physics
REPO=/project/bme-20255165-foldx/Mol-1
PY=$REPO/.venv/bin/python
ES=$R/miniforge3/envs/sigmadock
D=$R/data/cyp3a4_mol1
CIF=$R/data/cyp_structure_track/cif
TDIR=$D/templates/sim4
LIST=$D/templates/sim4_bundles.txt
MAN=$D/templates/sim4_manifest.csv
IN=$D/inputs_tsim4
POOLNAME=mol1_tsim4
export MOL_ROOT_DIR=$REPO/data
export PYTHONPATH=$REPO
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export TRITON_CACHE_DIR=$R/.triton_cache
export TORCH_EXTENSIONS_DIR=$R/protenix_data/torch_extensions
mkdir -p "$TORCH_EXTENSIONS_DIR"   # else torch builds into $HOME, under a 1M inode cap
export LAYERNORM_TYPE=torch        # featurisation only; keeps the CUDA kernel off a CPU node
export PATH=$R/miniforge3/envs/kalign/bin:$PATH   # find_kalign() asserts in TemplateFeaturizer.__init__
# POOLS/SEEDS are a --pools/--seeds RESTRICTION downstream and survive a plain reassignment.
unset POOL POOLS SEEDS MAXIMIZE
# shellcheck source=lib.sh
source "$R/slurm/cyp3a4_mol1/lib.sh" || { echo "FATAL: cannot source lib.sh" >&2; exit 1; }
cd "$R" || exit 1

echo "######## 0 gate: an unmixed pool"
NP=$(find "$D/gen/$POOLNAME" -name '*_sample_*.cif' 2>/dev/null | wc -l)
if [ "$NP" -gt 0 ] && [ "${ALLOW_RESUME:-0}" != 1 ]; then
    echo "FATAL: $D/gen/$POOLNAME already holds $NP poses -- a re-run would MIX two template" >&2
    echo "       selections in one pool. Move it aside, or set ALLOW_RESUME=1 deliberately." >&2
    exit 3
fi

echo "######## 1 select the per-ligand templates (precomputed ECFP4 matrix, no RDKit)"
"$ES/bin/python" analysis/cyp3a4_mol1/select_sim_templates.py \
    --out-manifest "$MAN" --out-list "$LIST" || exit 1
NL=$(wc -l < "$LIST")
[ "$NL" -eq 20 ] || { echo "FATAL: expected 20 bundle lines, got $NL" >&2; exit 1; }

echo "######## 2 build one 4-template JSON bundle per ligand"
# Its own subdir: make_template_json.py also writes <out-dir>/<PDB>.json per template as a side
# effect, which in $D/templates would add ~40 stray files and needlessly rewrite the existing
# 6MA6/6MA8/6BCZ bundles. Bundle name = the base_id, so --templates-suffix is a bare .json.
mkdir -p "$TDIR"
while read -r LIG PDBS; do
    [ -n "$LIG" ] || continue
    # shellcheck disable=SC2086
    "$PY" analysis/cyp3a4_mol1/make_template_json.py \
        --templates $PDBS --chain A --cif-dir "$CIF" \
        --out-dir "$TDIR" --bundle "$LIG" > /dev/null || { echo "FATAL: bundle $LIG" >&2; exit 1; }
    [ -f "$TDIR/$LIG.json" ] || { echo "FATAL: $TDIR/$LIG.json not written" >&2; exit 1; }
    # The cap is first-N with no sort, so a bundle that is not exactly 4 long changes the
    # selection silently. Assert the length of every one.
    NE=$("$ES/bin/python" -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" \
         "$TDIR/$LIG.json") || exit 1
    [ "$NE" -eq 4 ] || { echo "FATAL: $LIG.json holds $NE templates, expected 4" >&2; exit 1; }
done < "$LIST"
echo "  20 bundles of exactly 4 -> $TDIR"

echo "######## 3 the 20 blinded inputs, each pointed at ITS OWN bundle"
"$PY" analysis/cyp3a4_mol1/make_inputs.py --out-dir "$IN" \
    --templates-dir "$TDIR" --templates-suffix .json || exit 1

echo "######## 4 pre-GPU gate on ALL 20 -- one bundle per ligand means one gate per ligand"
# EXPECT_COVERAGE is deliberately omitted: residue coverage now differs per ligand, so there is
# no single expected vector. EXPECT_TEMPLATES=4 is the assertion that matters. ~33 s each.
NG=0
while read -r JSONP; do
    [ -n "$JSONP" ] || continue
    ( cd "$REPO" \
      && EXPECT_TEMPLATES=4 GATE_RUNNER_MODULE=mol.runner.inference \
         "$PY" "$R/analysis/cyp3a4_cryoem_tmpl/protenix_template_gate.py" \
           --load_checkpoint_path "$MOL_ROOT_DIR/checkpoint/mol-1-base.pt" \
           --input_json_path "$JSONP" --dump_dir "$D/gate_dump_tsim4" --seeds 101 \
           --use_msa true --use_rna_msa false --use_template true --device cpu > /dev/null \
    ) || { echo "GATE FAILED on $(basename "$JSONP") -- not launching the arm" >&2; exit 1; }
    NG=$((NG + 1))
done < "$IN/index.txt"
[ "$NG" -eq 20 ] || { echo "FATAL: gated $NG of 20" >&2; exit 1; }
echo "  gate PASSED on all $NG ligands: 4 templates reach the features for every one"

echo "######## 5 chunk + launch"
"$PY" analysis/cyp3a4_mol1/chunk_inputs.py --inputs-dir "$IN" --pool "$POOLNAME" \
    --chunk-size 2 || exit 1
IDX=$IN/chunks_index_$POOLNAME.txt
N=$(wc -l < "$IDX")
[ "$N" -gt 0 ] || { echo "nothing to run"; exit 0; }
export POOL=$POOLNAME IDX_FILE=$IDX USE_TEMPLATE=true MSA_MIN=16384 MSA_CUT=16384
# ⚠ tue.gpu2.q is ABSENT, not merely excluded: its only node IS the broken tue-gpub001, so
# naming it alongside the mandatory --exclude leaves zero eligible nodes and Slurm parks the
# array at BadConstraints forever (job 506974). 20_generate.sbatch's own header default names
# it, which is why --partition must be overridden here.
JG=$(sub --array=0-$((N - 1))%6 --job-name=$POOLNAME --time=03:30:00 --mem=20G \
     --partition=bme.gpumolml.q,bme.gpu2.q,tue.gpu3.q,tue.gpu.q \
     --exclude=tue-gpua001,tue-gpub001 --export=ALL "$R/slurm/cyp3a4_mol1/20_generate.sbatch") \
  || exit 1
echo "$POOLNAME generation ($N chunks x 2 ligands) -> $JG"
echo "then QC with:  export POOL=$POOLNAME; sbatch --export=ALL --job-name=\$POOL slurm/cyp3a4_mol1/30_collect_qc.sbatch"
echo "STAGE 119 PASSED"
