#!/bin/bash
#SBATCH --job-name=mol1_20_gen
#SBATCH --partition=tue.gpu2.q,tue.gpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/20_gen_%x_%A_%a.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/20_gen_%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=14:00:00
#SBATCH --mem=48G
set -uo pipefail
# Stage 20: 50 Mol-1 poses for ONE blinded ligand = seeds 101-110 x N_sample 5.
#   export POOL=mol1_full MSA_MIN=16384 MSA_CUT=16384
#   sbatch --array=0-19 --export=ALL --job-name=mol1_full 20_generate.sbatch
# POOL names the arm and is the `pool` label the ranker sees. MSA_MIN/MSA_CUT set the subsampling
# depth; leave BOTH unset for the repo default, which is NOT a full-MSA run.
#
# Why one ligand per array task rather than one 20-job input (which would load the model once):
# inference.py raises after finishing if ANY sample failed, so a single bad ligand would mark the
# whole 1000-pose run failed and muddle the bookkeeping. Per-ligand tasks isolate that and run
# 8-wide on the L4 queues, which beats one serial job on a single card.
#
# MSA DEPTH, stated once so the arm is interpretable: at inference the MSA rows are re-drawn EVERY
# RECYCLE (msa_module sits inside the N_cycle loop), uniformly at random without replacement, with
# the count drawn from U[min_size.test, n] and then truncated to sample_cutoff.test. Setting both
# to D therefore means "D random rows per recycle", not one fixed D-row MSA, and torch.randperm
# does not guarantee the query row survives. The draw is seeded by --seeds, so the 10 seeds see 10
# different subsets -- that is where the conformational diversity comes from.
R=/vast.mnt/project/bme-20255165-cofolding-physics
REPO=/project/bme-20255165-foldx/Mol-1
PY=$REPO/.venv/bin/python
D=$R/data/cyp3a4_mol1
export MOL_ROOT_DIR=$REPO/data
export PYTHONPATH=$REPO
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export TRITON_CACHE_DIR=$R/.triton_cache
export TORCH_EXTENSIONS_DIR=$R/protenix_data/torch_extensions
mkdir -p "$TORCH_EXTENSIONS_DIR"   # else torch builds into $HOME, which is under a 1M inode cap
# No nvcc on this cluster, so the fused LayerNorm can never build; this is the repo's escape hatch
# and leaves the cuEquivariance triangle kernels (prebuilt wheels) untouched.
export LAYERNORM_TYPE=torch
# Templates: find_kalign() runs in TemplateFeaturizer.__init__ whenever templates are enabled, even
# for .json templates, which are paired explicitly and never aligned. It is a startup assertion
# rather than a real dependency for that route, so kalign is simply put on PATH.
export PATH=$R/miniforge3/envs/kalign/bin:$PATH

IDX=${IDX_FILE:-$D/inputs/index.txt}
LIG=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$IDX")
[ -n "$LIG" ] || { echo "FATAL: no input at index $SLURM_ARRAY_TASK_ID"; exit 1; }
[ -f "$LIG" ] || { echo "FATAL: $LIG is not a file"; exit 1; }
OUT=$D/gen/${POOL:?set POOL}
mkdir -p "$OUT"
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader
echo "pool=$POOL ligand=$(basename "$LIG" .json) msa_min=${MSA_MIN:-default} msa_cut=${MSA_CUT:-default}"
# CKPT defaults to the released weights, so every existing arm reproduces unchanged. A
# fine-tuned arm MUST also get a distinct POOL name and its own POOL_CKPT entry in
# analysis/cyp3a4_ob0_mmgbsa/merge_scored.py, or a pick gets attributed to stock Mol-1.
echo "ckpt=${CKPT:-$MOL_ROOT_DIR/checkpoint/mol-1-base.pt}"

# GUIDANCE=1 turns on the TFG chemistry guidance back-ported into mol/tfg/ -- RDKit-derived
# physical-validity potentials that steer the sampler at inference, training-free. The hook is
# generator.py:239-259, inside the denoising loop; with guidance off, parse_tfg_config(None)
# returns a fully disabled config and the sampler takes the verbatim stock Euler step.
#
# The flags DO reach the parser on this entrypoint. The comment at mol/config/base.py:271-274
# claims inference takes no config CLI overrides; that is true only of mol/runner/batch_inference.py,
# which calls parse_configs() with no arg_str. inference.py:779-806 passes parse_sys_args() twice,
# and `guidance` is registered inside model_configs, so --sample_diffusion.guidance.* is accepted.
#
# ⚠⚠ A GUIDED RUN IS NOT SEED-COMPARABLE TO AN UNGUIDED ONE, even at zero weight: mol/tfg/engine.py
# :493 computes `x_next + sigma * randn_like(x_next)` and then returns x_next, discarding the value
# but still consuming the RNG stream. The stock path draws one torch.randn per step, the TFG path
# two, so the two streams desynchronize from step 0. Treating an existing unguided pool as a
# controlled baseline therefore confounds guidance with resampling; for a clean attribution run a
# third arm with GUIDANCE=1 and every term weight 0.
#
# ⚠ Leave PocketPotential alone (interval 0): its residue set is read from PROTENIX_POCKET_JSON /
# _ID, not the MOL_* names its own config comment comments, and its chain keys are hardcoded to
# A_NS2B / B_NS3 from an unrelated campaign, so it can never match a CYP3A4 pocket.
GOPT=()
if [ "${GUIDANCE:-0}" = 1 ]; then
    # fail before the card is used, not after: an absent port would run SILENTLY unguided
    [ -f "$REPO/mol/tfg/engine.py" ] || { echo "FATAL: GUIDANCE=1 but $REPO/mol/tfg/engine.py is missing"; exit 1; }
    GOPT=(--sample_diffusion.guidance.enable true)
    [ -n "${TFG_MU:-}" ] && GOPT+=(--sample_diffusion.guidance.mu "$TFG_MU")
    [ -n "${TFG_INNER:-}" ] && GOPT+=(--sample_diffusion.guidance.steps.tfg_inner "$TFG_INNER")
    echo "guidance=ON ${GOPT[*]}"
else
    echo "guidance=off"
fi

cd "$REPO" || exit 1
"$PY" mol/runner/inference.py \
  --load_checkpoint_path "${CKPT:-$MOL_ROOT_DIR/checkpoint/mol-1-base.pt}" \
  --input_json_path "$LIG" \
  --dump_dir "$OUT" \
  --seeds "${SEEDS:-101,102,103,104,105,106,107,108,109,110}" \
  --use_seeds_in_json false \
  --model.N_cycle 20 \
  --sample_diffusion.N_step 200 \
  --sample_diffusion.N_sample 5 \
  --use_template "${USE_TEMPLATE:-false}" \
  --use_msa true \
  --use_rna_msa false \
  --need_atom_confidence false \
  --device auto \
  --triangle_attention cuequivariance \
  --triangle_multiplicative cuequivariance \
  ${MSA_MIN:+--data.msa.min_size.test "$MSA_MIN"} \
  ${MSA_CUT:+--data.msa.sample_cutoff.test "$MSA_CUT"} \
  ${GOPT[@]+"${GOPT[@]}"}
rc=$?
# The input may hold several ligands (chunk_*.json), so count the arm, not one subdirectory.
echo "arm poses now: $(find "$OUT" -name '*_sample_*.cif' 2>/dev/null | wc -l) (50 per ligand)"
# Did guidance actually run? The engine logs "Guidance is enabled." plus one last-step energy line
# per term to STDERR, which is this task's --error file. Reconstructing that path is the cheapest
# no-op assertion available, and the campaign's recurring failure is a flag that parses and does
# nothing. Non-fatal on purpose: a wrong path guess must not fail a finished 10-hour task.
if [ "${GUIDANCE:-0}" = 1 ]; then
    EF=$R/slurm/cyp3a4_mol1/logs/20_gen_${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err
    if [ -f "$EF" ]; then
        NG=$(grep -c "Guidance is enabled." "$EF" 2>/dev/null || echo 0)
        NE=$(grep -c "TFG last-step energy" "$EF" 2>/dev/null || echo 0)
        echo "guidance witness: $NG 'enabled' line(s), $NE energy line(s) in $(basename "$EF")"
        [ "$NG" -gt 0 ] || echo "⚠⚠ GUIDANCE=1 but the sampler never logged 'Guidance is enabled.'"
    else
        echo "guidance witness: could not locate $EF -- grep the .err by hand"
    fi
fi
echo "exit=$rc"
exit $rc
