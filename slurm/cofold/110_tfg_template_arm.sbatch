#!/bin/bash
#SBATCH --job-name=mol1_tfg_drv
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/110_tfgarm_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/110_tfgarm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --mem=4G
set -uo pipefail
# Stage 110: a TEMPLATED Mol-1 arm over the 20 blinded ligands with the TFG chemistry guidance ON.
#   POOLNAME=mol1_tfg_tcryo5 TEMPLATE=cryo5 sbatch --export=ALL 110_tfg_template_arm.sbatch
#
# This is the first TFG run in this campaign on (a) the STOCK mol-1-base.pt checkpoint and (b) the
# 20 blinded CYP3A4 ligands. Every prior TFG run used a fine-tuned checkpoint on the OpenBind
# blind400 set through a different launcher, so nothing here is a reproduction -- the guidance is
# proven to move the trajectory, but not proven to help on this target.
#
# WHAT THE GUIDANCE CAN AND CANNOT TOUCH HERE. The RDKit potentials are driven ONLY by the blinded
# drug: geometry_featurizer.py skips standard residues (exclude_std_residue) and skips
# metal-containing residues wholesale, so HEM contributes ZERO RDKit terms. HEM and protein enter
# only through VinaStericPotential and InterchainBondPotential. So the testable hypothesis is about
# LIGAND INTERNAL GEOMETRY -- bond lengths, torsions, planarity, stereo -- not about pocket
# placement or Fe coordination. The campaign's standing ligand-geometry complaint is exactly that:
# worst C-C bond sits at 1.55-1.59 A across every pool. poses_qc_<pool>.csv's pb_pass and
# max_cc_bond_A are therefore the metrics this arm has to move.
#
# WHY THE MATCHED CONTROL IS NOT A CLEAN CONTROL. mol1_tcryo5 already exists at stock weights with
# the SAME inputs and the SAME seeds 101-110, so it is the natural comparison -- but
# mol/tfg/engine.py:493 consumes an extra torch.randn per step that it then discards, so a guided
# run desynchronizes from an unguided one at the same seed even at zero guidance weight. The pair
# answers "are the guided poses better"; it does NOT isolate guidance from resampling. A
# zero-weight third arm is the clean control if that attribution is wanted.
#
# ⚠⚠ THE TFG PORT IS UNCOMMITTED. mol/tfg/ and mol/data/core/geometry_featurizer.py are UNTRACKED
# on bare `main` in the Mol-1 repo, with 5 tracked hook files unstaged, no branch and no stash. The
# array tasks read that working tree at TASK START, so a `git checkout .` or `git stash` in the repo
# while this arm is queued would silently revert the hooks and leave the untracked package behind --
# an unguided run under a guided pool name. Do not touch the repo's git state until the arm drains.
R=/vast.mnt/project/bme-20255165-cofolding-physics
REPO=/project/bme-20255165-foldx/Mol-1
PY=$REPO/.venv/bin/python
D=$R/data/cyp3a4_mol1
POOLNAME=${POOLNAME:-mol1_tfg_tcryo5}
TEMPLATE=${TEMPLATE:-cryo5}
TMPL=$D/templates/$TEMPLATE.json
IN=$D/inputs_${POOLNAME#mol1_}
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH POOL SEEDS
cd "$R" || exit 1
source "$R/slurm/cyp3a4_mol1/lib.sh" || exit 1

echo "######## 0 gate: the port, the template, and an unmixed pool"
# A missing port is the one failure that would produce a plausible-looking UNGUIDED arm under a
# guided name, so it is checked here as well as in the task script.
for F in mol/tfg/engine.py mol/tfg/potentials.py mol/tfg/config.py mol/data/core/geometry_featurizer.py; do
    [ -f "$REPO/$F" ] || { echo "FATAL: TFG port incomplete, missing $REPO/$F" >&2; exit 3; }
done
echo "TFG port present: $(cd "$REPO" && git status --porcelain mol/tfg mol/data/core/geometry_featurizer.py | tr '\n' ' ')"
[ -f "$TMPL" ] || { echo "FATAL: no template bundle $TMPL" >&2; exit 3; }
# Refusing to mix guidance on/off (or checkpoints) inside one pool, same rule as 98_ft_template_arms.
N_EXIST=$(find "$D/gen/$POOLNAME" -name '*_sample_*.cif' 2>/dev/null | wc -l)
if [ "$N_EXIST" -gt 0 ] && [ "${ALLOW_RESUME:-0}" != 1 ]; then
    echo "FATAL: $D/gen/$POOLNAME already holds $N_EXIST poses; set ALLOW_RESUME=1 to pack only the" >&2
    echo "       unfinished ligands (chunk_inputs.py skips complete ones), or pick a new POOLNAME." >&2
    exit 3
fi

echo
echo "######## 1 build inputs (same protein, same MSA, same SMILES as every other arm)"
"$PY" analysis/cyp3a4_mol1/make_inputs.py --out-dir "$IN" --template-json "$TMPL" || exit 1
NJ=$(find "$IN" -maxdepth 1 -name 'OCNT-*.json' | wc -l)
[ "$NJ" -eq 20 ] || { echo "FATAL: built $NJ inputs, expected 20" >&2; exit 1; }

echo
echo "######## 2 chunk (2 ligands/task: templates + full MSA do not fit a 22GiB L4 otherwise)"
"$PY" analysis/cyp3a4_mol1/chunk_inputs.py --inputs-dir "$IN" --pool "$POOLNAME" \
    --chunk-size 2 || exit 1
IDX=$IN/chunks_index_$POOLNAME.txt
N=$(wc -l < "$IDX")
[ "$N" -gt 0 ] || { echo "FATAL: no entries in $IDX" >&2; exit 1; }

echo
echo "######## 3 launch the guided array"
# WALL: the unguided templated arms used 03:30:00 for a 2-ligand chunk. The measured TFG cost is
# 2.7x wall (271s -> 738s on the smoke pair), so 12h carries the 2.7x with headroom rather than
# losing a whole chunk to a timeout at hour four.
# PARTITIONS: tue.gpu2.q is DELIBERATELY ABSENT. Its only node IS the broken tue-gpub001, so
# requesting it alongside the mandatory --exclude leaves zero eligible nodes and Slurm parks the
# array at BadConstraints forever (job 506974). The header default of 20_generate.sbatch names it,
# which is why it is overridden here rather than merely excluded.
export POOL=$POOLNAME IDX_FILE=$IDX USE_TEMPLATE=true MSA_MIN=16384 MSA_CUT=16384 GUIDANCE=1
JG=$(sub --array=0-$((N - 1))%6 --job-name="$POOLNAME" --time=12:00:00 --mem=28G \
     --partition=bme.gpumolml.q,bme.gpu2.q,tue.gpu3.q,tue.gpu.q \
     --exclude=tue-gpua001,tue-gpub001 --export=ALL "$R/slurm/cyp3a4_mol1/20_generate.sbatch") || exit 1
echo "guided array -> $JG  ($N tasks x 2 ligands, template=$TEMPLATE, ckpt=stock mol-1-base.pt)"
echo
echo "witness to check once the first task starts:"
echo "  grep -l 'Guidance is enabled.' $R/slurm/cyp3a4_mol1/logs/20_gen_${POOLNAME}_${JG}_*.err"
echo "then QC with:  export POOL=$POOLNAME; sbatch --export=ALL --job-name=\$POOL slurm/cyp3a4_mol1/30_collect_qc.sbatch"
echo "⚠ a submission built from this pool needs its own POOL_CKPT entry in"
echo "  analysis/cyp3a4_ob0_mmgbsa/merge_scored.py, or picks get attributed to stock unguided Mol-1."
echo "STAGE 110 PASSED"
