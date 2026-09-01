#!/usr/bin/env bash
#
# run_pipeline.sh -- run polkat setups in sequence, waiting for each stage's
# SLURM jobs to finish before generating the next stage's submit script.
#
# Run from the project root, inside screen or tmux -- this blocks for as long
# as the pipeline takes, and an SSH drop would otherwise stop the driver from
# advancing to the next stage (already-submitted jobs keep running).
#
#     screen -S polkat
#     mv extra/run_pipeline.sh .; chmod +x run_pipeline.sh; ./run_pipeline.sh idia 2>&1 | tee pipeline.log
#     # Ctrl-A d to detach, screen -r polkat to reattach
#
# Note: if the cluster round-robins you across login nodes, `screen -ls` will
# be empty on reconnect unless you SSH back to the same host.
#
# Start part-way through:  ./run_pipeline.sh idia 2
#
set -uo pipefail

CLUSTER=${1:-idia}
START=${2:-0}
POLL=${POLL:-300}

STAGES=(
  "setups/0_GET_INFO.py|submit_info_job.sh"
  "setups/1GC.py|submit_1GC_jobs.sh"
  "extra/setup_2GC_twostage.py|submit_2GC_jobs.sh"
  # "setups/RMSYNTH.py|submit_rmsynth_jobs.sh"
)

log() { echo "[$(date '+%F %T')] $*"; }

# Source the submit script with sbatch shadowed by a function, so we see every
# job ID while the script's own `| awk '{print $4}'` still works.
submit_and_collect() {
  local script=$1 idfile=$2
  : > "$idfile"
  (
    sbatch() {
      local out
      out=$(command sbatch "$@") || { echo "$out" >&2; return 1; }
      awk '{print $NF}' <<< "$out" >> "$idfile"
      echo "$out"
    }
    source "./$script"
  )
}

wait_for() {
  local ids=$1 q err
  while :; do
    err=$(squeue -h -j "$ids" -o "%i %T %r" 2>&1 >/tmp/squeue_out.$$)
    q=$(</tmp/squeue_out.$$); rm -f /tmp/squeue_out.$$
    if [[ -n $err ]]; then
      log "squeue query failed, retrying: $err"
      sleep "$POLL"
      continue
    fi
    [[ -z $q ]] && break
    if grep -q DependencyNeverSatisfied <<< "$q"; then
      log "dependency never satisfied -- cancelling remaining jobs"
      scancel $(awk '{print $1}' <<< "$q") 2>/dev/null
      return 1
    fi
    log "$(grep -c RUNNING <<< "$q") running, $(grep -c PENDING <<< "$q") pending"
    sleep "$POLL"
  done
  local bad
  bad=$(sacct -j "$ids" -X -n -o JobID%20,State%20 | grep -Ev 'COMPLETED' || true)
  if [[ -n $bad ]]; then
    log "non-COMPLETED jobs:"; echo "$bad"; return 1
  fi
  return 0
}

for ((i=START; i<${#STAGES[@]}; i++)); do
  setup=${STAGES[i]%%|*}
  subs=${STAGES[i]##*|}

  log "=== stage $i: $setup ==="
  python3 "$setup" "$CLUSTER" || { log "setup failed"; exit 1; }
  [[ -f ./$subs ]] || { log "$subs not generated"; exit 1; }

  idfile=.stage${i}.jobids
  submit_and_collect "$subs" "$idfile"
  mapfile -t ids < "$idfile"
  [[ ${#ids[@]} -gt 0 ]] || { log "no jobs submitted"; exit 1; }

  csv=$(IFS=,; echo "${ids[*]}")
  log "submitted ${#ids[@]} jobs: $csv"
  wait_for "$csv" || { log "stage $i failed -- stopping"; exit 1; }
  log "stage $i complete"
done

log "pipeline finished"
