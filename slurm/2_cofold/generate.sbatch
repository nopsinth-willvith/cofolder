#!/bin/bash
#SBATCH --job-name=cofold_gen
#SBATCH --partition=bme.gpumolml.q,bme.gpu2.q,tue.gpu3.q,tue.gpu.q
#SBATCH --exclude=tue-gpua001,tue-gpub001
#SBATCH --output=slurm/2_cofold/logs/gen_%x_%A_%a.out
#SBATCH --error=slurm/2_cofold/logs/gen_%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=03:30:00
#SBATCH --mem=20G
set -uo pipefail
# STAGE 2a -- co-folding. One array task = one chunk JSON (CHUNK_SIZE ligands, model loaded once).
# Submitted by `slurm/run.sh cofold <arm.env>`, which sets --array/--time/--mem/--partition from
# the config. Every knob below comes from the arm config, never from --export.
#
# MSA DEPTH: rows are re-drawn EVERY RECYCLE, uniformly without replacement, count drawn from
# U[MSA_MIN, n] then truncated to MSA_CUT, seeded by --seeds. MSA_MIN = MSA_CUT = D therefore means
# "D random rows per recycle" (and the query row is not guaranteed to survive). Leaving both unset
# in Mol-1 is a STOCHASTIC-depth arm, not a full-MSA one, so lib.sh pins both to 16384.
#
# GUIDANCE=1 turns on the TFG chemistry guidance (vendored in mol1/). A guided run is NOT
# seed-comparable to an unguided one even at zero weight (mol/tfg/engine.py:493 draws an extra
# randn per step). PocketPotential must stay at interval 0 (chain keys hardcoded from another
# campaign).
source "${COFOLDER:?submit via slurm/run.sh}/slurm/lib.sh" && cofold_load || exit 1
mol1_env || exit 1

LIG=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$IDX")
[ -n "$LIG" ] || { echo "FATAL: no input at index $SLURM_ARRAY_TASK_ID of $IDX"; exit 1; }
[ -f "$LIG" ] || { echo "FATAL: $LIG is not a file"; exit 1; }
mkdir -p "$GEN"
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader
USE_TEMPLATE=false; [ "$TEMPLATE_MODE" != none ] && USE_TEMPLATE=true
echo "pool=$POOL chunk=$(basename "$LIG") ckpt=$CKPT ckpt_label=$CKPT_LABEL"
echo "templates=$USE_TEMPLATE ($TEMPLATE_MODE) msa_min=${MSA_MIN:-default} msa_cut=${MSA_CUT:-default} seeds=$SEEDS"
echo "N_cycle=$N_CYCLE N_step=$N_STEP N_sample=$N_SAMPLE"

GOPT=()
if [ "$GUIDANCE" = 1 ]; then
    # fail before the card is used: an absent port would run SILENTLY unguided under a guided name
    for F in mol/tfg/engine.py mol/tfg/potentials.py mol/tfg/config.py mol/data/core/geometry_featurizer.py; do
        [ -f "$MOL1/$F" ] || { echo "FATAL: GUIDANCE=1 but the TFG port is incomplete: $MOL1/$F"; exit 1; }
    done
    GOPT=(--sample_diffusion.guidance.enable true)
    [ -n "${TFG_MU:-}" ] && GOPT+=(--sample_diffusion.guidance.mu "$TFG_MU")
    [ -n "${TFG_INNER:-}" ] && GOPT+=(--sample_diffusion.guidance.steps.tfg_inner "$TFG_INNER")
    echo "guidance=ON ${GOPT[*]}"
else
    echo "guidance=off"
fi

cd "$MOL1" || exit 1
"$MOL1_PY" mol/runner/inference.py \
  --load_checkpoint_path "$CKPT" \
  --input_json_path "$LIG" \
  --dump_dir "$GEN" \
  --seeds "$SEEDS" \
  --use_seeds_in_json false \
  --model.N_cycle "$N_CYCLE" \
  --sample_diffusion.N_step "$N_STEP" \
  --sample_diffusion.N_sample "$N_SAMPLE" \
  --use_template "$USE_TEMPLATE" \
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
echo "arm poses now: $(npose "$GEN") (expect $EXPECT_POSES per ligand)"

# WITNESS, not flag: the TFG engine logs "Guidance is enabled." once per inference call and one
# "TFG last-step energy" line per term per sample to STDERR. Non-fatal on purpose: a wrong path
# guess must not fail a finished multi-hour task.
if [ "$GUIDANCE" = 1 ]; then
    EF=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | grep -oP 'StdErr=\K\S+' | head -1)
    if [ -n "$EF" ] && [ -f "$EF" ]; then
        NG=$(grep -c "Guidance is enabled." "$EF" 2>/dev/null); NG=${NG:-0}
        NE=$(grep -c "TFG last-step energy" "$EF" 2>/dev/null); NE=${NE:-0}
        echo "guidance witness: $NG 'enabled' line(s), $NE energy line(s) in $(basename "$EF")"
        [ "$NG" -gt 0 ] || echo "⚠⚠ GUIDANCE=1 but the sampler never logged 'Guidance is enabled.'"
    else
        echo "guidance witness: could not locate this task's stderr -- grep the .err by hand"
    fi
fi
echo "exit=$rc"
exit $rc
