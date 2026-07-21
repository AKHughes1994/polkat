#!/usr/bin/env bash
set -euo pipefail

# ---------------------------
# Auto-run inside screen
# ---------------------------
if [[ -z "${STY:-}" && "${1:-}" != "--inside-screen" ]]; then
  SESSION="global_run_$(date +%Y%m%d_%H%M%S)"
  SELF="$(readlink -f "$0" 2>/dev/null || realpath "$0")"

  if ! command -v screen >/dev/null 2>&1; then
    echo "ERROR: 'screen' not found in PATH." >&2
    exit 127
  fi

  screen -dmS "$SESSION" bash -lc "'$SELF' --inside-screen"
  echo "Started screen session: $SESSION"
  echo "Attach:  screen -r $SESSION"
  echo "Status:  tail -f \"$(dirname "$SELF")/global_benchmark.log\""
  exit 0
fi

# ---------------------------
# Run from script directory
# ---------------------------
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

LOG="$ROOT/global_benchmark.log"
: > "$LOG"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >> "$LOG"; }

run_step() {
  local name="$1"; shift
  local start end rc tmp
  tmp="$(mktemp)"

  log "START  ${name}"
  start="$(date +%s)"

  # Run command; do NOT write full output to the global log.
  # Capture only the last ~200 lines for debugging if it fails.
  if "$@" 2>&1 | tail -n 200 > "$tmp"; then
    end="$(date +%s)"
    awk -v s="$start" -v e="$end" -v n="$name" \
      'BEGIN{printf("[%s] DONE   %s   elapsed_hours=%.4f\n", strftime("%F %T"), n, (e-s)/3600)}' >>"$LOG"
    rm -f "$tmp"
  else
    rc=$?
    end="$(date +%s)"
    awk -v s="$start" -v e="$end" -v n="$name" -v r="$rc" \
      'BEGIN{printf("[%s] FAIL   %s   rc=%d   elapsed_hours=%.4f\n", strftime("%F %T"), n, r, (e-s)/3600)}' >>"$LOG"
    log "---- last 200 lines of output ----"
    cat "$tmp" >> "$LOG"
    log "---------------------------------"
    rm -f "$tmp"
    exit "$rc"
  fi
}

overall_start="$(date +%s)"
log "GLOBAL START  root=$ROOT"

run_step "INFO setup"      python3 setups/0_GET_INFO.py node
run_step "INFO submit"     bash ./submit_info_job.sh

run_step "1GC setup"       python3 setups/1GC.py node
run_step "1GC submit"      bash ./submit_1GC_jobs.sh

run_step "2GC setup"       python3 setups/2GC.py node
run_step "2GC submit"      bash ./submit_2GC_jobs.sh

run_step "RMSYNTH setup"   python3 setups/RMSYNTH.py node
run_step "RMSYNTH submit"  bash ./submit_rmsynth_jobs.sh

#run_step "SNAP setup"   python3 setups/SNAP.py node
#run_step "SNAP submit"  bash ./submit_snap_jobs.sh

overall_end="$(date +%s)"
awk -v s="$overall_start" -v e="$overall_end" \
  'BEGIN{printf("[%s] GLOBAL DONE  elapsed_hours=%.4f\n", strftime("%F %T"), (e-s)/3600)}' >>"$LOG"

log "Benchmark written to: $LOG"

