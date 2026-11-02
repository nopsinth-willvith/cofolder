#!/bin/bash
#SBATCH --job-name=ob0gbsa
#SBATCH --partition=tue.cpu1.q,tue.cpu2.q,tue.cpu3.q
#SBATCH --output=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_ob0_mmgbsa/logs/gbsa_%A_%a.out
#SBATCH --error=/project/bme-20255165-cofolding-physics/slurm/cyp3a4_ob0_mmgbsa/logs/gbsa_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G
set -uo pipefail
# sbatch --array=0-<N-1> --export=ALL,SET=oracle|blinded,NSHARDS=<N> 10_build_score.sbatch
# Per pose (~15 min, serial): build_system.py (ff14SB + Shahrokh IC6 heme + GAFF1/AM1-BCC,
# res 25-503, 10/16 A pocket) then mmgbsa_modes.py (explicit-solvent minimisation with
# backbone restrained, then igb=2 MM-GBSA). --no-fe-bond for EVERY pose: the Fe-N contact stays
# non-bonded, so ante-MMPBSA's :LIG split severs nothing (the project's fix for the Type-II
# bonded-term and cpptraj-exclusion artifacts) and all poses share one physical model.
R=/vast.mnt/project/bme-20255165-cofolding-physics
E=$R/miniforge3/envs/openmm83
# WORK_ROOT and OUT_PREFIX are overridable, and a non-default PREP_PH *requires* overriding BOTH.
# Two silent-corruption traps otherwise:
#  1. the build-skip guard below reuses any existing system.prmtop, and ~9,500 pH-7.4 systems sit
#     under the default WORK_ROOT -- a pH-4.8 run would skip every build and re-score the OLD
#     systems, emitting a CSV labelled 4.8 that is numerically 7.4;
#  2. the default OUT_PREFIX would overwrite the existing pH-7.4 shard CSVs in place.
W=${WORK_ROOT:-/scratch-shared/20255165/cyp3a4_ob0_mmgbsa}
OUT_PREFIX=${OUT_PREFIX:-$R/data/cyp3a4_ob0_mmgbsa/mmgbsa_${SET}_shard}
P=$R/data/cyp3a4_ob0_mmgbsa/poselist_${SET}.csv
if [ "${PREP_PH:-7.4}" != "7.4" ]; then
    [ "$W" != "/scratch-shared/20255165/cyp3a4_ob0_mmgbsa" ] || {
        echo "FATAL: PREP_PH=$PREP_PH with the DEFAULT WORK_ROOT -- the build-skip guard would" >&2
        echo "       re-score the existing pH-7.4 systems. Set WORK_ROOT to a pH-tagged dir." >&2
        exit 3; }
    case "$OUT_PREFIX" in
        *"/mmgbsa_${SET}_shard") echo "FATAL: PREP_PH=$PREP_PH with the DEFAULT OUT_PREFIX --" >&2
            echo "       this would overwrite the pH-7.4 results. Tag OUT_PREFIX." >&2; exit 3;;
    esac
fi
echo "pH=${PREP_PH:-7.4}  WORK_ROOT=$W  OUT_PREFIX=$OUT_PREFIX"
SH=$SLURM_ARRAY_TASK_ID
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
unset PYTHONPATH
export AMBERHOME=$E PATH=$E/bin:$PATH
export LD_LIBRARY_PATH=$E/lib:${LD_LIBRARY_PATH:-}

# build this shard's systems: the same i % NSHARDS split mmgbsa_modes.py uses
awk -F, -v s=$SH -v n=$NSHARDS 'NR>1 && (NR-2)%n==s {print $1","$2","$3}' "$P" |
while IFS=, read -r POSE PDB SMI; do
  [ -f "$W/systems/$POSE/$POSE/system.prmtop" ] && continue
  mkdir -p "$W/systems/$POSE"
  "$E/bin/python" "$R/analysis/cyp3a4_gamd/build_system.py" \
    --pose "$PDB" --smiles "$SMI" --pose-id "$POSE" --out-dir "$W/systems/$POSE" \
    --first-residue 25 --no-fe-bond --pocket-free-radius 10 --pocket-buffer-radius 16 \
    > "$W/systems/$POSE/build.log" 2>&1 || echo "$POSE build FAILED (see build.log)"
done

"$E/bin/python" "$R/analysis/cyp3a4_bioemu2/mmgbsa_modes.py" \
  --poselist "$P" --systems-root "$W/systems" --work-dir "$W/work" \
  --out-csv "${OUT_PREFIX}${SH}.csv" \
  --shard "$SH" --nshards "$NSHARDS"
