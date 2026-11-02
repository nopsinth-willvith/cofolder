#!/bin/bash
#SBATCH --job-name=mol1_tmplsim
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/121_tmplsim_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/121_tmplsim_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mem=8G
set -uo pipefail
# Stage 121: lowest-MM-GBSA pick over ALL FOUR templated Mol-1 pools competing together.
#
#   mol1_tsim4    910 poses   per-ligand templates by bound-ligand ECFP4 similarity (stage 119)
#   mol1_tcryo5   860         fixed cryo-EM bundle (4 of its 5 entries; the cap drops the 5th)
#   mol1_tcomb4   930         fixed X-ray bundle, the four OF3 actually consumed
#   mol1_t38FX    859         single cryo-EM entry
#
# Template-free mol1_full is DELIBERATELY EXCLUDED: this zip answers "best templated pose", and
# the mol1_full comparison already exists as submission_mol1_tmpl_mmgbsa (which includes it) vs
# submission_mol1_tmplonly_mmgbsa (which does not).
#
# The selector is copied verbatim from 60_tmpl_pipeline.sbatch:197 -- same COLS, same Fe gate,
# same fallback ladder -- so the ONLY difference from submission_mol1_tmplonly_mmgbsa is that
# mol1_tsim4's 910 poses are now in the candidate pool.
#
# No add_repaired: mol1_tsim4 and mol1_tcomb4 each cover all 20 blinded ligands on their own, so
# no ligand is starved and the zip needs no grafted rows from another pool. That makes this a
# pure four-pool templated zip.
R=/vast.mnt/project/bme-20255165-cofolding-physics
D=$R/data/cyp3a4_mol1
G=$R/data/cyp3a4_ob0_mmgbsa
NAME_OUT=submission_mol1_tmplsim_mmgbsa
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
# MAXIMIZE flips the direction and delta_total is an ENERGY: a leaked value ships the 20 WORST
# binders silently. SEEDS would restrict the candidate pool and FATALs once several pools compete.
unset MAXIMIZE SEEDS
# shellcheck source=lib.sh
source "$R/slurm/cyp3a4_mol1/lib.sh" || exit 1
cd "$R" || exit 1

echo "######## 0 gate: both shard families present, and no existing zip to destroy"
for f in "$G/poselist_mol1tmpl.csv" "$G/poselist_mol1tsim4.csv"; do
    [ -s "$f" ] || { echo "FATAL: missing $f" >&2; exit 3; }
done
NT=$(ls "$G"/mmgbsa_mol1tmpl_shard*.csv 2>/dev/null | wc -l)
NS=$(ls "$G"/mmgbsa_mol1tsim4_shard*.csv 2>/dev/null | wc -l)
echo "  mol1tmpl shards $NT | mol1tsim4 shards $NS"
[ "$NT" -gt 0 ] && [ "$NS" -gt 0 ] || { echo "FATAL: a shard family is empty" >&2; exit 3; }
# 50_submission rm -rf's the audit dir and overwrites the zip/provenance in place, with no warning.
if [ -f "$D/submission/$NAME_OUT.zip" ] && [ "${ALLOW_OVERWRITE:-0}" != 1 ]; then
    echo "FATAL: $NAME_OUT.zip exists -- it would be destroyed. ALLOW_OVERWRITE=1 to proceed." >&2
    exit 3
fi

echo "######## 1 merge the two shard families into one candidate table"
# Space-separated word lists, exported in the SHELL and inherited via --export=ALL: sbatch splits
# an --export=ALL,VAR=a,b value on the comma and only the first item would arrive.
export SCORED="$G/mmgbsa_mol1tmpl_shard*.csv $G/mmgbsa_mol1tsim4_shard*.csv" \
       POSELISTS="$G/poselist_mol1tmpl.csv $G/poselist_mol1tsim4.csv" \
       OUT=scored_mol1_tmplsim.csv
JM=$(sub --export=ALL "$R/slurm/cyp3a4_ob0_mmgbsa/45_merge_scored.sbatch") || exit 1
echo "merge -> $JM"

echo "######## 2 build + validate + audit"
# POOLS names all four explicitly. It is a self-documenting check here rather than a real
# restriction (the merged table holds exactly these four), and the builder prints what it kept.
export NAME="$NAME_OUT" SCORES="$G/scored_mol1_tmplsim.csv" \
       COLS="delta_total" FE_IN=2.6 FE_RETAIN=3.5 FE_LADDER="3.0:4.0,3.5:4.5" \
       POOLS="mol1_tsim4 mol1_tcryo5 mol1_tcomb4 mol1_t38FX" OUTDIR="$D/submission" \
       REFS="$D/submission/submission_mol1_tsim4_mmgbsa.zip $D/submission/submission_mol1_tmplonly_mmgbsa.zip $D/submission/submission_mol1_tmpl_mmgbsa.zip"
JB=$(sub --dependency=afterok:"$JM" --kill-on-invalid-dep=yes --export=ALL \
     "$R/slurm/cyp3a4_ob0_mmgbsa/50_submission.sbatch") || exit 1
echo "build+validate+audit -> $JB"
echo "chain: $JM -> $JB"
echo "STAGE 121 PASSED"
