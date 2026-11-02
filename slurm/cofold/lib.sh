# Shared helpers for the slurm/cyp3a4_mol1 stages. Sourced, never executed.
#
# sub(): submit a job and echo its id, with BOTH the exit status and the id checked. An unchecked
# sbatch is the specific failure this exists to prevent: it leaves --dependency= empty, which Slurm
# accepts, so the dependent job starts immediately on inputs that do not exist yet -- and a stage
# whose last command is an echo then exits 0 with nothing queued.
#
# `${id%%;*}` is load-bearing: a federated Slurm reply to --parsable is "<jobid>;<cluster>", so a
# `^[0-9]+$` test would FATAL on a submission that actually succeeded. Strip the cluster, then
# require the remainder to be non-empty. (Same handling as 9x_rescore_chain.sh:25.)
sub() {
    local id
    id=$(sbatch --parsable "$@") || { echo "SUBMIT FAILED: $*" >&2; return 1; }
    id=${id%%;*}
    [ -n "$id" ] || { echo "SUBMIT exited 0 but printed no job id: $*" >&2; return 1; }
    echo "$id"
}
