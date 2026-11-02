#!/bin/bash
#SBATCH --job-name=mol1_tfg_drv
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/112_tfg_gbsa_%x_%j.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_mol1/logs/112_tfg_gbsa_%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mem=8G
set -uo pipefail

# Stage 112: the MM-GBSA ranking pipeline over the TFG-guided templated pool.
#
# Inputs are stock mol-1-base.pt + cryo5 templates + TFG chemistry guidance ON (stage 110,
# array 524051: 1000 poses, 20 ligands x 10 seeds x 5 samples, guidance witnessed on all 200
# inference calls). The RANKING is copied verbatim from the reference build
# (60_tmpl_pipeline.sbatch:197) -- same COLS, same Fe gate, same fallback ladder -- so the only
# difference between this zip and submission_mol1_tmplonly_mmgbsa is which poses competed.
#
# Why its own SET and its own POOL_CKPT label: the WEIGHTS are stock mol-1-base, but the sampler
# is not. mol/tfg/engine.py:493 consumes an extra randn per step, so these are not the unguided
# mol1_tcryo5 poses under another name, and labelling them "mol1" would attribute a pick to
# stock unguided Mol-1. merge_scored.py:117-120 FATALs on an undeclared pool, which is the guard.
#
# Submit chained on the stage 30 QC job:
#   sbatch --dependency=afterany:<qc jobid> 112_tfg_mmgbsa.sbatch
# afterany, not afterok: stage 30 returns 1 for a skipped pose as well as for a real failure, so
# this gates on the table's own content below instead of trusting the exit status.
R=/vast.mnt/project/bme-20255165-cofolding-physics
ES=$R/miniforge3/envs/sigmadock
D=$R/data/cyp3a4_mol1
G=$R/data/cyp3a4_ob0_mmgbsa
# ARM names the variant and DERIVES every output name from it. Before this, POOL and SET were
# parameterised but NAME/OUT were hardcoded to tfg, so driving this file for a second pool would
# have silently overwritten submission_mol1_tfg_mmgbsa and its scored tables. The defaults below
# reproduce the original TFG invocation exactly.
ARM=${ARM:-tfg}
POOL=${POOL:-mol1_tfg_tcryo5}
SET=${SET:-mol1$ARM}
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
# MAXIMIZE flips build_submission's direction and delta_total is an ENERGY: a leaked value would
# ship the 20 WORST-binding poses and nothing downstream would complain. SEEDS would restrict the
# candidate pool and FATAL if it emptied one. Both cleared, as stages 60 and 100 do.
unset MAXIMIZE SEEDS
cd "$R" || exit 1
source "$R/slurm/cyp3a4_mol1/lib.sh" || exit 1

QC=$D/poses_qc_$POOL.csv
echo "######## 0 gate: QC table, pb_pass coverage, and no prior scoring of this SET"
[ -s "$QC" ] || { echo "FATAL: missing QC table $QC -- run stage 30 for $POOL first" >&2; exit 3; }
NROW=$(($(wc -l < "$QC") - 1))
[ "$NROW" -gt 0 ] || { echo "FATAL: $QC has a header and no rows" >&2; exit 3; }
echo "  $(basename "$QC"): $NROW rows"

# A poselist that already has shards means this stage ran before; re-running would double-score.
if [ -s "$G/poselist_$SET.csv" ] && compgen -G "$G/mmgbsa_${SET}_shard*.csv" > /dev/null; then
    echo "FATAL: $SET already has a poselist AND shards -- refusing to re-score. Move them aside" >&2
    echo "       or pass a different SET if you mean to build a second variant." >&2
    exit 3
fi

# Coverage decides the shape of the chain, so it is measured, not assumed. A submission needs one
# pose for each of the 20 blinded ligands: build_submission.py:204-206 FATALs on a starved ligand
# and 50_submission.sbatch:66 fails the audit on a short zip. OCNT-2313038 is 0/50 pb_pass in
# every unguided arm but mol1_tcomb4 (planar cyclohexane), and TFG's whole hypothesis is ligand
# internal geometry -- so this is also the arm's first real read-out.
COV=$("$ES/bin/python" - "$QC" <<'PY'
import csv, sys, collections
rows = list(csv.DictReader(open(sys.argv[1], newline="")))
tru = lambda v: str(v).strip().lower() in ("true", "1", "yes")
npass, tot = collections.Counter(), collections.Counter()
for r in rows:
    tot[r["base_id"]] += 1
    if tru(r.get("pb_pass")):
        npass[r["base_id"]] += 1
starved = sorted(l for l in tot if npass[l] == 0)
print(sum(npass.values()), len(tot) - len(starved), len(tot), ",".join(starved) or "-")
for l in sorted(tot):
    print(f"#  {l} {npass[l]:3d}/{tot[l]}" + ("   <-- ZERO" if npass[l] == 0 else ""), file=sys.stderr)
PY
) || { echo "FATAL: coverage check failed" >&2; exit 3; }
read -r NPASS NCOV NLIG STARVED <<< "$COV"
echo "  pb_pass $NPASS/$NROW poses; $NCOV/$NLIG ligands covered; starved: $STARVED"
[ "$NPASS" -gt 0 ] || { echo "FATAL: no pb_pass pose in the whole pool" >&2; exit 3; }

echo "######## 1 poselist (every pb_pass pose; no top-N, no srs cut)"
"$ES/bin/python" analysis/cyp3a4_mol1/poselist_from_qc.py \
    --qc "$QC" --out "$G/poselist_$SET.csv" || exit 1
N=$(($(wc -l < "$G/poselist_$SET.csv") - 1))
[ "$N" -gt 0 ] || { echo "FATAL: empty poselist" >&2; exit 3; }

# 10-14 CPU-min per pose measured over the mol1full/mol1tmpl/mol1exp arrays, serial, one core per
# shard. ~10 poses/shard lands a shard near 2 h, far inside 10_build_score's 12 h limit.
NSH=${NSHARDS:-$(( (N + 9) / 10 ))}
[ "$NSH" -lt 1 ] && NSH=1
[ "$NSH" -gt 100 ] && NSH=100
echo "######## 2 MM-GBSA: $N poses over $NSH shards (~$(( N * 12 / 60 )) CPU-h)"
export SET NSHARDS=$NSH
JG=$(sub --array=0-$((NSH - 1))%${THROTTLE:-60} --job-name=${SET}_gbsa --export=ALL \
     "$R/slurm/cyp3a4_ob0_mmgbsa/10_build_score.sbatch") || exit 1
echo "MM-GBSA array -> $JG"

echo "######## 3 merge this SET's shards only -- a TFG-only candidate pool is the point"
export SCORED="$G/mmgbsa_${SET}_shard*.csv" POSELISTS="$G/poselist_$SET.csv" \
       OUT=scored_mol1_$ARM.csv
JM=$(sub --dependency=afterany:"$JG" --kill-on-invalid-dep=yes --export=ALL \
     "$R/slurm/cyp3a4_ob0_mmgbsa/45_merge_scored.sbatch") || exit 1
echo "merge -> $JM"

# The repair rows are grafted poses from mol1_full, NOT TFG poses. They are appended only when a
# ligand would otherwise starve the zip, and then they are the single documented impurity in it.
if [ "$STARVED" = "-" ]; then
    echo "######## 4 skip add_repaired: all $NLIG ligands covered by TFG poses (pure arm)"
    SCORES_CSV=$G/scored_mol1_$ARM.csv
    SEL_POOLS="$POOL"
    JDEP=$JM
else
    echo "######## 4 append the repaired rows -- $STARVED has no pb_pass TFG pose"
    SCORES_CSV=$D/scored_mol1_${ARM}_full.csv
    SEL_POOLS="$POOL mol1_repair"
    JDEP=$(sub --dependency=afterok:"$JM" --kill-on-invalid-dep=yes --job-name=mol1_addrep_tfg \
         --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q --nodes=1 --ntasks-per-node=1 \
         --cpus-per-task=2 --time=00:15:00 --mem=8G \
         --output="$R/slurm/cyp3a4_mol1/logs/addrep_tfg_%j.out" \
         --error="$R/slurm/cyp3a4_mol1/logs/addrep_tfg_%j.err" \
         --wrap="cd $R && PYTHONNOUSERSITE=1 $ES/bin/python analysis/cyp3a4_mol1/add_repaired.py \
            --table $G/scored_mol1_$ARM.csv \
            --scores $D/mmgbsa_mol1fix.csv \
            --out $D/scored_mol1_${ARM}_full.csv") || exit 1
    echo "add repaired -> $JDEP"
fi

if [ -f "$D/submission/submission_mol1_${ARM}_mmgbsa.zip" ] && [ "${ALLOW_OVERWRITE:-0}" != 1 ]; then
    echo "FATAL: submission_mol1_${ARM}_mmgbsa.zip already exists -- 50_submission would destroy" >&2
    echo "       its zip, provenance and audit. Pass a different ARM, or ALLOW_OVERWRITE=1." >&2
    exit 3
fi

echo "######## 5 build + validate + audit the zip"
# COLS/FE_IN/FE_RETAIN/FE_LADDER verbatim from 60_tmpl_pipeline.sbatch:197. POOLS is exported in
# the SHELL and never as --export=ALL,POOLS=a,b: sbatch splits that value on the comma and only
# the first pool would arrive.
export NAME="submission_mol1_${ARM}_mmgbsa" SCORES="$SCORES_CSV" \
       COLS="delta_total" FE_IN=2.6 FE_RETAIN=3.5 FE_LADDER="3.0:4.0,3.5:4.5" \
       POOLS="$SEL_POOLS" OUTDIR="$D/submission" \
       REFS="$D/submission/submission_mol1_tmplonly_mmgbsa.zip $D/submission/submission_mol1_tmpl_mmgbsa.zip"
JB=$(sub --dependency=afterok:"$JDEP" --kill-on-invalid-dep=yes --export=ALL \
     "$R/slurm/cyp3a4_ob0_mmgbsa/50_submission.sbatch") || exit 1
echo "build+validate+audit -> $JB"
echo "chain: $JG -> $JM -> $JDEP -> $JB   (pools: $SEL_POOLS)"
echo "STAGE 112 PASSED"
