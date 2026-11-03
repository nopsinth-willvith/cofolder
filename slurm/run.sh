#!/bin/bash
# slurm/run.sh <stage> <arm.env> -- submit one pipeline stage for one arm (config/arms/<arm>.env).
#
#   template  STAGE 1  CPU job: select templates -> bundles -> inputs -> pre-GPU gate -> chunks
#   cofold    STAGE 2  GPU array (one chunk per task) -> QC job (afterany) -> pb_pass poselist
#   rerank    STAGE 3  scorer A, MM-GBSA array -> merge (afterany) -> scorer B, MLFF chain (afterok)
#   mlff      STAGE 3  scorer B only, on an existing MM-GBSA table: shortlist -> Orb array -> merge
#   select    STAGE 3  pick per ligand by SCORER -> zip -> official validator -> zip audit -> diff
#   status             what exists for this arm and what is queued (submits nothing)
#
# Login-node safe: it reads the config, counts lines and calls sbatch -- no Python, no compute.
# Prints the job ids; every log lands in slurm/<stage>/logs/. DRY_RUN=1 prints the sbatch
# commands instead of submitting. One-off overrides go in the SHELL, never in --export:
#   ALLOW_RESUME=1     finish only what is missing in an existing pool or MM-GBSA shard set
#   ALLOW_OVERWRITE=1  let `mlff` / `select` replace this arm's existing E_int table / zip
set -uo pipefail
STAGE=${1:-}; ARM=${2:-}
if [ -z "$STAGE" ] || [ -z "$ARM" ]; then
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'; exit 2
fi
COFOLDER=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CFG=$(readlink -f "$ARM")
export COFOLDER CFG
# shellcheck source=lib.sh
source "$COFOLDER/slurm/lib.sh"
cofold_load || exit 1
cd "$COFOLDER" || exit 1
L=$COFOLDER/slurm
mkdir -p "$L"/{1_template,2_cofold,3_rerank}/logs

queued() { squeue -h --me -o "%i %j %T" 2>/dev/null | awk -v p="$1" '$2 ~ p'; }
refuse_if_queued() {       # refuse_if_queued <regex> <what>: never double-submit a stage
    local q; q=$(queued "$1")
    [ -z "$q" ] || { echo "FATAL: $2 already queued/running:" >&2; echo "$q" >&2; exit 3; }
}
cpu() { echo --partition="$CPU_PARTITION"; }
log() { echo -o "$L/$1/logs/%x_$2.out" -e "$L/$1/logs/%x_$2.err"; }   # log <stage dir> <%j|%A_%a>

# scorer B: shortlist (MM-GBSA top-K) -> prep + Orb array -> merge. $1 = job to follow, or empty.
# Guarded here rather than in the `mlff` stage alone: reached through `rerank` it would otherwise
# overwrite an E_int table a submission may already have shipped from.
mlff_chain() {
    local dep=() JS JO JE
    if [ -f "$EINT" ] && [ "${ALLOW_OVERWRITE:-0}" != 1 ]; then
        echo "FATAL: $EINT exists -- the MLFF chain would overwrite it. ALLOW_OVERWRITE=1 to" >&2
        echo "       re-score (prepped npz are reused), or MLFF_TOPK=0 to skip scorer B." >&2
        return 3
    fi
    [ -n "$1" ] && dep=(--dependency="afterok:$1" --kill-on-invalid-dep=yes)
    # shellcheck disable=SC2046
    JS=$(sub $(cpu) ${dep[@]+"${dep[@]}"} --job-name="mlffsl_$POOL" --cpus-per-task=2 --mem=8G \
         --time=00:30:00 $(log 3_rerank %j) "$L/3_rerank/mlff.sbatch" shortlist) || exit 1
    # shellcheck disable=SC2046
    JO=$(sub $(cpu) --dependency="afterok:$JS" --kill-on-invalid-dep=yes --array="0-$((MLFF_SHARDS - 1))" \
         --job-name="mlff_$POOL" $(log 3_rerank %A_%a) "$L/3_rerank/mlff.sbatch" score) || exit 1
    # afterany: a failed shard must not strand the merge; merge_eint reports the poses it lacks
    # shellcheck disable=SC2046
    JE=$(sub $(cpu) --dependency="afterany:$JO" --kill-on-invalid-dep=yes --job-name="mlffm_$POOL" \
         --cpus-per-task=1 --mem=4G --time=00:15:00 $(log 3_rerank %j) "$L/3_rerank/mlff.sbatch" merge) || exit 1
    echo "        MLFF shortlist (top-$MLFF_TOPK/ligand) -> $JS -> Orb array ($MLFF_SHARDS) -> $JO -> merge -> $JE"
    echo "        table: $EINT"
}

case "$STAGE" in
template)
    refuse_if_queued "^($POOL|tmpl_$POOL)$" "stage 1/2 for $POOL (re-chunking rewrites the files tasks read)"
    NP=$(npose "$GEN")
    if [ "$NP" -gt 0 ] && [ "${ALLOW_RESUME:-0}" != 1 ]; then
        echo "FATAL: $GEN already holds $NP poses -- rebuilding would MIX two configurations in one" >&2
        echo "       pool. Pick a new POOL, or ALLOW_RESUME=1 to re-chunk only the unfinished ligands." >&2
        exit 3
    fi
    # shellcheck disable=SC2046
    J=$(sub $(cpu) --job-name="tmpl_$POOL" $(log 1_template %j) "$L/1_template/template.sbatch") || exit 1
    echo "stage 1 template -> $J   log: $L/1_template/logs/tmpl_${POOL}_$J.out (ends 'STAGE 1 PASSED')"
    ;;
cofold)
    [ -f "$IDX" ] || { echo "FATAL: no chunk index $IDX -- run stage 1 (template) first" >&2; exit 3; }
    # tmpl_ too: stage 1 rewrites $IDX, and the array below is sized from it while tasks resolve
    # their chunk at TASK START -- a concurrent re-chunk would feed them files that moved.
    refuse_if_queued "^(tmpl_$POOL|$POOL|qc_$POOL)$" "stage 1/2 for $POOL"
    N=$(grep -c . "$IDX")
    NP=$(npose "$GEN")
    if [ "$N" -eq 0 ]; then
        echo "nothing to generate: every ligand in $POOL is complete ($NP poses) -- submitting QC only"
        # shellcheck disable=SC2046
        JQ=$(sub $(cpu) --job-name="qc_$POOL" $(log 2_cofold %j) "$L/2_cofold/qc.sbatch") || exit 1
        echo "QC -> $JQ   table: $QC"; exit 0
    fi
    if [ "$NP" -gt 0 ] && [ "${ALLOW_RESUME:-0}" != 1 ]; then
        echo "FATAL: $GEN already holds $NP poses; ALLOW_RESUME=1 to finish the remaining chunks" >&2; exit 3
    fi
    # tue.gpu2.q is ABSENT from GPU_PARTITION, not merely excluded: its only node is the broken
    # tue-gpub001, so naming it beside --exclude leaves zero eligible nodes (BadConstraints forever).
    # shellcheck disable=SC2046
    JG=$(sub --array="0-$((N - 1))%$GEN_THROTTLE" --job-name="$POOL" --partition="$GPU_PARTITION" \
         --exclude="$GPU_EXCLUDE" --time="$GEN_TIME" --mem="$GEN_MEM" $(log 2_cofold %A_%a) \
         "$L/2_cofold/generate.sbatch") || exit 1
    # afterany, not afterok: one failed chunk must not strand QC; QC reports per-target counts.
    # shellcheck disable=SC2046
    JQ=$(sub $(cpu) --dependency="afterany:$JG" --kill-on-invalid-dep=yes --job-name="qc_$POOL" \
         $(log 2_cofold %j) "$L/2_cofold/qc.sbatch") || exit 1
    echo "stage 2 generate -> $JG ($N chunks x $CHUNK_SIZE ligands, expect $EXPECT_POSES poses/ligand)"
    echo "        QC       -> $JQ (afterany)   table: $QC   poselist: $POSELIST"
    # `if`, not `[ ] && echo`: a trailing failed test is the stage's exit status
    if [ "$GUIDANCE" = 1 ]; then
        echo "        witness: grep -c 'Guidance is enabled.' $L/2_cofold/logs/${POOL}_${JG}_*.err"
    fi
    ;;
rerank)
    [ "$SCORED" = "$SCORED_DEFAULT" ] || {
        echo "FATAL: this arm's SCORED points at an existing merged table ($SCORED); MM-GBSA would" >&2
        echo "       overwrite it. Use 'mlff' or 'select' for this arm." >&2; exit 3; }
    refuse_if_queued "^(gbsa|merge)_$SET$|^mlff(sl|m)?_$POOL$" "stage 3 for $SET"
    N=$(nrows "$POSELIST")
    [ "$N" -gt 0 ] || { echo "FATAL: empty or missing poselist $POSELIST -- run stage 2 (QC) first" >&2; exit 3; }
    NSH=$(( (N + POSES_PER_SHARD - 1) / POSES_PER_SHARD )); [ "$NSH" -gt 100 ] && NSH=100
    ARRAY="0-$((NSH - 1))"
    if compgen -G "${SHARDS}[0-9]*.csv" > /dev/null; then
        [ "${ALLOW_RESUME:-0}" = 1 ] || {
            echo "FATAL: $SET already has MM-GBSA shards ($(compgen -G "${SHARDS}[0-9]*.csv" | wc -l)) --" >&2
            echo "       refusing to re-score. ALLOW_RESUME=1 scores only the missing shards; a variant needs a new SET." >&2
            exit 3; }
        # resume: the i % NSHARDS split is fixed by the first launch, so reuse its shard count
        [ -s "${SHARDS}.nshards" ] || { echo "FATAL: shards exist but ${SHARDS}.nshards does not" >&2; exit 3; }
        NSH=$(cat "${SHARDS}.nshards")
        ARRAY=$(for i in $(seq 0 $((NSH - 1))); do [ -s "${SHARDS}$i.csv" ] || printf '%s,' "$i"; done)
        ARRAY=${ARRAY%,}
    elif [ -f "$SCORED" ]; then
        echo "FATAL: $SCORED exists but $SET has no shards -- it came from elsewhere; pick a new SET" >&2; exit 3
    fi
    # The shard COUNT is capped at 100, so past ~5000 poses the per-shard wall, not the count,
    # is what grows: 14 CPU-min/pose x poses-per-shard + 30 min of build overhead. A shard that
    # times out leaves no CSV, and ALLOW_RESUME would re-queue it to time out again.
    PER=$(( (N + NSH - 1) / NSH ))
    MIN=$(( PER * 14 + 30 )); [ "$MIN" -lt 120 ] && MIN=120
    if [ "$MIN" -gt 1440 ]; then
        echo "FATAL: $N poses over the 100-shard cap is $PER poses/shard = ~$((MIN / 60)) h per" >&2
        echo "       shard, past the 24 h this stage will request. Split the pool across two SETs." >&2
        exit 3
    fi
    GBSA_TIME=$(printf '%02d:%02d:00' $((MIN / 60)) $((MIN % 60)))
    echo "stage 3 MM-GBSA: $N pb_pass poses, $NSH shards x $PER poses, --time=$GBSA_TIME, ~$(( N * 12 / 60 )) CPU-h (10-14 CPU-min/pose)"
    [ "${DRY_RUN:-0}" = 1 ] || echo "$NSH" > "${SHARDS}.nshards"
    JDEP=()
    if [ -n "$ARRAY" ]; then
        # shellcheck disable=SC2046
        JG=$(sub $(cpu) --array="$ARRAY%$GBSA_THROTTLE" --job-name="gbsa_$SET" --time="$GBSA_TIME" \
             $(log 3_rerank %A_%a) "$L/3_rerank/mmgbsa.sbatch") || exit 1
        echo "        MM-GBSA array -> $JG  (shards ${SHARDS}<i>.csv)"
        JDEP=(--dependency="afterany:$JG" --kill-on-invalid-dep=yes)
    fi
    # shellcheck disable=SC2046
    JM=$(sub $(cpu) ${JDEP[@]+"${JDEP[@]}"} --job-name="merge_$SET" $(log 3_rerank %j) \
         "$L/3_rerank/merge.sbatch") || exit 1
    echo "        merge -> $JM (afterany: content-gated)   table: $SCORED"
    if [ "$MLFF_TOPK" -gt 0 ]; then mlff_chain "$JM" || exit 1; else echo "        MLFF off (MLFF_TOPK=0)"; fi
    ;;
mlff)
    [ "$MLFF_TOPK" -gt 0 ] || { echo "FATAL: MLFF_TOPK=0 in $CFG" >&2; exit 3; }
    [ -s "$SCORED" ] || { echo "FATAL: no MM-GBSA table $SCORED -- run 'rerank' first" >&2; exit 3; }
    refuse_if_queued "^mlff(sl|m)?_$POOL$" "MLFF for $POOL"
    echo "stage 3 MLFF on $SCORED (pools: $RANK_POOLS)"
    mlff_chain "" || exit 1
    ;;
select)
    refuse_if_queued "^sel_$NAME$" "select for $NAME"
    if [ -f "$SUB/$NAME.zip" ] && [ "${ALLOW_OVERWRITE:-0}" != 1 ] && [ "$SELECT_MODE" != audit ]; then
        echo "FATAL: $SUB/$NAME.zip exists -- select would destroy its zip, provenance and audit." >&2
        echo "       Set a new NAME in the arm config, or ALLOW_OVERWRITE=1." >&2; exit 3
    fi
    if [ "$SELECT_MODE" = rank ]; then
        for f in $RANK_SCORES; do [ -s "$f" ] || { echo "FATAL: no score table $f (RANK_SCORES)" >&2; exit 3; }; done
        [ "$ADD_REPAIRED" = 1 ] && [ "$SCORER" != mmgbsa ] && {
            echo "FATAL: ADD_REPAIRED=1 carries MM-GBSA energies only; it cannot compete under SCORER=$SCORER" >&2; exit 3; }
    fi
    # shellcheck disable=SC2046
    J=$(sub $(cpu) --job-name="sel_$NAME" $(log 3_rerank %j) "$L/3_rerank/select.sbatch") || exit 1
    echo "stage 3 select ($SELECT_MODE, scorer=$SCORER col=$RANK_COL) -> $J   zip: $SUB/$NAME.zip"
    ;;
status)
    echo "arm $POOL  (config $CFG)"
    echo "  ckpt_label=$CKPT_LABEL templates=$TEMPLATE_MODE guidance=$GUIDANCE msa=${MSA_MIN:-default}/${MSA_CUT:-default} seeds=$SEEDS N_sample=$N_SAMPLE"
    f() { printf '  %-10s %-4s %s\n' "$1" "$([ -e "$2" ] && echo ok || echo --)" "$2${3:+  ($3)}"; }
    f inputs "$IN/index.txt" "$( [ -f "$IN/index.txt" ] && grep -c . "$IN/index.txt") ligands"
    f chunks "$IDX" "$( [ -f "$IDX" ] && grep -c . "$IDX") to run"
    f poses "$GEN" "$(npose "$GEN") of $((20 * EXPECT_POSES))"
    f qc "$QC"
    f poselist "$POSELIST" "$(nrows "$POSELIST") pb_pass"
    f shards "${SHARDS}.nshards" "$(compgen -G "${SHARDS}[0-9]*.csv" | wc -l) of $(cat "${SHARDS}.nshards" 2>/dev/null || echo '?')"
    f mmgbsa "$SCORED" "$(nrows "$SCORED") rows"
    f mlff "$EINT" "$(nrows "$EINT") rows"
    f zip "$SUB/$NAME.zip"
    echo "  queued/running:"; queued "^($POOL|(tmpl|qc|mlffsl|mlff|mlffm)_$POOL|(gbsa|merge)_$SET|sel_$NAME)$" | sed 's/^/    /'
    ;;
*)
    echo "unknown stage '$STAGE' (template|cofold|rerank|mlff|select|status)" >&2; exit 2 ;;
esac
