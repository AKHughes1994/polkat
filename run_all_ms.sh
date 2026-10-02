#!/usr/bin/env bash
#
# run_all_ms.sh -- run the full pipeline (0_GET_INFO -> 1GC -> 2GC -> RMSYNTH)
# on every *.ms in a directory, one MS at a time.
#
# Fill in the INPUTS globals below, then run:
#     ./run_all_ms.sh 2>&1 | tee run_all_ms.log
#
# Place this script one level above working_dir. For each MS, working_dir is
# emptied, filled with a copy of the pipeline directory and a symlink to the
# MS, and the stages are run inside it. A stage only starts once every job of
# the previous one has finished. An MS counts as done once its obsid (leading
# number of the MS name) is in the tracking file, which RMSYNTH writes as its
# very last step. An MS that fails is reported and the run moves on to the next.
#
# MSs whose obsid is already in the tracking file are skipped. Before starting,
# it lists what it found and asks for confirmation.
#
# On slurm the job queue is checked every POLL seconds (default 300):
#     POLL=600 ./run_all_ms.sh
#
# ------------------------------------------------------------------ #
# What a run actually does, step by step (for readers who don't know bash):
#
#   1. Checks the settings below are all filled in and valid, and exits
#      immediately, saying what's wrong, if not.
#   2. Reads the pipeline's own oxkat/config.py and checks it has an archive
#      directory and a tracking file configured, and that its tracking file
#      matches the one passed on the command line -- otherwise a finished
#      MS's products would be silently discarded, or never picked up as done.
#   3. Lists every *.ms in ms_dir, splits it into "already done" (its obsid
#      is in the tracking file) and "to run", prints the plan, and asks for
#      confirmation before doing anything.
#   4. For each MS still to run, in order:
#        a. Empties working_dir, copies a fresh copy of the pipeline into it,
#           and symlinks that MS in.
#        b. Runs each pipeline stage in turn (0_GET_INFO, 1GC, 2GC, RMSYNTH):
#           generates that stage's job scripts, submits them (or runs them
#           directly, in "node" mode), and waits for every job to finish
#           before starting the next stage.
#        c. Stops at the first stage that fails for that MS, records which
#           stage that was in the failures file, and moves on to the next
#           MS -- one bad MS never blocks the others.
#   5. After each MS, checks the tracking file again: if RMSYNTH's last step
#      added this obsid, that MS is counted as passed. If not, it's counted as
#      failed -- unless STAGES was deliberately shortened to stop before
#      RMSYNTH (e.g. for a quick test) and every stage that did run actually
#      succeeded, in which case it's counted as PARTIAL, not failed.
#   6. Once every MS has been tried, prints how many passed, failed, and
#      stopped early, and lists the failed and partial ones by name.
# ------------------------------------------------------------------ #
set -uo pipefail

# ------------------------------------------------------------------ #
# INPUTS -- fill these in, then run the script with no arguments
# ------------------------------------------------------------------ #

INFRA='node'      # 'idia', 'hippo', or 'node'
MS_DIR='/mnt/extraspace/tkat_reprocessing'     # directory containing the *.ms to process
TRACKING='/mnt/scratchhdd/tkat_reprocessing/tracking/mahrez_tracking.txt'   # tracking file; must match RMSYNTH_TRACKING_FILE in the pipeline's oxkat/config.py
PIPELINE='/mnt/scratchhdd/tkat_reprocessing/polkat_tkat_reprocessing'   # pipeline directory to run, e.g. 'polkat_tkat_reprocessing'
WORK='/mnt/scratchhdd/tkat_reprocessing/working_dir'       # scratch working directory, rebuilt from scratch for every MS -- must be dedicated to this script, not shared with anything else
POLL=${POLL:-300}   # slurm queue poll interval in seconds, overridable via POLL=... in the environment
FIRST_ONLY=false    # true: run only the first MS that is not already in the tracking file (debugging)

# Refuse to start unless every input above has actually been filled in.
REQUIRED=(INFRA MS_DIR TRACKING PIPELINE WORK)
MISSING=()
for name in "${REQUIRED[@]}"; do
  [[ -z ${!name} ]] && MISSING+=("$name")
done
if [[ ${#MISSING[@]} -gt 0 ]]; then
  echo "You need to specify: ${MISSING[*]} -- edit the globals at the top of $0" >&2
  exit 1
fi

# All path-like inputs must be absolute: this script cd's into $WORK partway
# through a run, so a relative path given here would silently resolve against
# the wrong directory from that point on.
PATH_INPUTS=(MS_DIR TRACKING PIPELINE WORK)
for name in "${PATH_INPUTS[@]}"; do
  [[ ${!name} = /* ]] || { echo "$name must be an absolute path (got: '${!name}')" >&2; exit 1; }
done

# MODE decides how each pipeline stage is run further down: submitted to
# slurm's job queue (idia, hippo) and waited on, or just run directly on
# this node (node) -- see run_stage_slurm/run_stage_node below.
case $INFRA in
  idia|hippo) MODE=slurm ;;
  node)       MODE=node ;;
  *)          echo "INFRA must be 'idia', 'hippo', or 'node' (got: '$INFRA')" >&2; exit 1 ;;
esac

# ROOT is fixed to this script's own location; STATE keeps each MS's own logs
# and submitted job IDs around after WORK has been wiped for the next MS. WORK
# itself is the INPUTS global above, not derived from ROOT, so it can live
# anywhere (a different disk, a different scratch area) independent of where
# this script and its state directory happen to sit.
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
STATE="$ROOT/run_all_ms_logs"

# The 4 pipeline stages, in the order they must run, and which script each
# stage's setup step is expected to generate to submit that stage's jobs.
# Comment out trailing entries for a quick partial test run (e.g. just
# 0_GET_INFO) -- see STAGES_INCLUDE_RMSYNTH below for what that changes.
STAGES=(
  "setups/0_GET_INFO.py|submit_info_job.sh"
  "setups/1GC.py|submit_1GC_jobs.sh"
  "setups/2GC.py|submit_2GC_jobs.sh"
  "setups/RMSYNTH.py|submit_rmsynth_jobs.sh"
)

# TRACKING is only ever written by RMSYNTH's own last step, so a run that
# doesn't include that stage can never end up in TRACKING even if every stage
# it did run succeeded. The main loop uses this to tell that apart from an
# actual failure, instead of reporting a deliberate partial run as FAILED.
STAGES_INCLUDE_RMSYNTH=false
[[ " ${STAGES[*]} " == *"RMSYNTH.py"* ]] && STAGES_INCLUDE_RMSYNTH=true

log() { echo "[$(date '+%F %T')] $*"; }

# Confirm the given directories actually look like what they're supposed to
# be before doing anything with them.
[[ -d $MS_DIR ]] || { echo "MS directory not found: $MS_DIR" >&2; exit 1; }
[[ -f $PIPELINE/setups/0_GET_INFO.py ]] || { echo "Not a pipeline directory: $PIPELINE" >&2; exit 1; }
PIPELINE="$(cd -- "$PIPELINE" && pwd)"
case "$PIPELINE/" in
  "$WORK"/*) echo "Pipeline directory must not be inside $WORK" >&2; exit 1 ;;
esac

# ------------------------------------------------------------------ #
# Pipeline's own move/tracking settings must line up with the args above.
# clear_work() wipes working_dir before every MS, so a finished MS's products
# only survive if RMSYNTH_ARCHIVE_DIR moves them out first; and this script's
# own resume check (in_tracking(), below) only sees what RMSYNTH_TRACKING_FILE
# is set to, since that is the file RMSYNTH actually writes completions to.
# ------------------------------------------------------------------ #

CONFIG="$PIPELINE/oxkat/config.py"
[[ -f $CONFIG ]] || { echo "Pipeline is missing oxkat/config.py: $CONFIG" >&2; exit 1; }

cfg_value() {
  sed -n "s/^$1[[:space:]]*=[[:space:]]*'\([^']*\)'.*/\1/p" "$CONFIG" | head -n1
}

# Resolves to an absolute path even if the file doesn't exist yet (TRACKING is
# created on first write), so it can be compared against config.py's value.
# 2>/dev/null: discard readlink's stderr if it fails, so the || fallback below
# can run silently instead of also printing an error.
canon() { readlink -f -- "$1" 2>/dev/null || echo "$1"; }

CFG_ARCHIVE_DIR=$(cfg_value RMSYNTH_ARCHIVE_DIR)
CFG_TRACKING_FILE=$(cfg_value RMSYNTH_TRACKING_FILE)

if [[ -z $CFG_ARCHIVE_DIR ]]; then
  echo "You need to specify RMSYNTH_ARCHIVE_DIR in $CONFIG -- without it, each MS's" >&2
  echo "finished products are never moved out of working_dir, so clear_work() deletes" >&2
  echo "them when the next MS starts." >&2
  exit 1
fi

if [[ -z $CFG_TRACKING_FILE ]]; then
  echo "You need to specify RMSYNTH_TRACKING_FILE in $CONFIG -- without it, RMSYNTH" >&2
  echo "never writes completion records, so no MS is ever picked up as done." >&2
  exit 1
fi

if [[ $(canon "$TRACKING") != $(canon "$CFG_TRACKING_FILE") ]]; then
  echo "TRACKING ($TRACKING) does not match RMSYNTH_TRACKING_FILE in $CONFIG" >&2
  echo "($CFG_TRACKING_FILE) -- RMSYNTH writes completion records to its own config" >&2
  echo "value, not to the file passed here, so this script's tracking file would" >&2
  echo "never actually be updated. Pass the same path for both." >&2
  exit 1
fi

# Sits next to TRACKING (e.g. .../mahrez_tracking.txt -> .../mahrez_tracking_failed.txt).
# Unlike TRACKING, nothing else reads this -- it's purely this script's own
# record of which stage each failed MS got stuck on, appended to like TRACKING.
FAILED_FILE="$(dirname -- "$TRACKING")/$(basename -- "$TRACKING" .txt)_failed.txt"

# The leading number of an MS's filename, used as that observation's unique
# ID everywhere below (the tracking file, the per-MS state directory, ...).
obsid() {
  local name=${1##*/}
  echo "${name%%_*}"
}

# True if obsid $1 already has a completed line in the tracking file -- i.e.
# that MS was already run successfully and can be skipped.
in_tracking() {
  [[ -f $TRACKING ]] || return 1
  awk -F'|' -v id="$1" '
    /^#/ { next }
    { gsub(/[ \t]/, "", $1); if ($1 == id) found = 1 }
    END { exit !found }' "$TRACKING"
}

# Appends one line to FAILED_FILE recording which stage a run got stuck on,
# same append-if-missing-header shape as TRACKING.
record_failure() {
  local id=$1 msname=$2 stage=$3 reason=$4
  mkdir -p -- "$(dirname -- "$FAILED_FILE")"
  [[ -f $FAILED_FILE ]] || echo '# obsid | ms_name | failed_utc | stage | reason' > "$FAILED_FILE"
  echo "$id | $msname | $(date -u '+%Y-%m-%dT%H:%M:%S') | $stage | $reason" >> "$FAILED_FILE"
}

# ------------------------------------------------------------------ #
# What needs running
# ------------------------------------------------------------------ #

mapfile -t ALL_MS < <(find "$MS_DIR" -mindepth 1 -maxdepth 1 -name '*.ms' \( -type d -o -type l \) | sort)

if [[ ${#ALL_MS[@]} -eq 0 ]]; then
  echo "No *.ms found in $MS_DIR"
  exit 0
fi

TODO=()
DONE=()
for ms in "${ALL_MS[@]}"; do
  if in_tracking "$(obsid "$ms")"; then
    DONE+=("$ms")
  else
    TODO+=("$ms")
  fi
done

[[ $FIRST_ONLY == true || $FIRST_ONLY == false ]] || { echo "FIRST_ONLY must be 'true' or 'false' (got: '$FIRST_ONLY')" >&2; exit 1; }

# With FIRST_ONLY, only the first MS to run is kept
HELD_BACK=0
if [[ $FIRST_ONLY == true && ${#TODO[@]} -gt 1 ]]; then
  HELD_BACK=$(( ${#TODO[@]} - 1 ))
  TODO=("${TODO[0]}")
fi

echo
echo "=================================================================="
echo " SUMMARY"
echo "=================================================================="
echo "  Infrastructure : $INFRA ($MODE)"
echo "  Pipeline       : $PIPELINE"
echo "  Working dir    : $WORK (emptied before each MS)"
echo "  MS directory   : $MS_DIR"
if [[ -f $TRACKING ]]; then
  echo "  Tracking file  : $TRACKING"
else
  echo "  Tracking file  : $TRACKING (does not exist yet -- nothing has been run)"
fi
echo "  Failures file  : $FAILED_FILE (which stage a failed MS got stuck on)"
echo "  Archive dir    : $CFG_ARCHIVE_DIR (from config.py -- each MS's products move here)"
stage_names=()
for entry in "${STAGES[@]}"; do
  name=${entry%%|*}; name=${name#setups/}; name=${name%.py}
  stage_names+=("$name")
done
stages_display=$(printf ' -> %s' "${stage_names[@]}"); stages_display=${stages_display# -> }
if [[ $STAGES_INCLUDE_RMSYNTH == true ]]; then
  echo "  Stages         : $stages_display"
else
  echo "  Stages         : $stages_display  (no RMSYNTH -- will never reach TRACKING, reported as PARTIAL not FAILED)"
fi
[[ $MODE == slurm ]] && echo "  Queue poll     : every ${POLL}s"
echo
echo "  Found ${#ALL_MS[@]} MS file(s): ${#DONE[@]} already in tracking, ${#TODO[@]} to run"
if [[ $FIRST_ONLY == true ]]; then
  echo "  FIRST_ONLY is true: only the first MS to run is kept ($HELD_BACK more not run)"
fi
if [[ ${#DONE[@]} -gt 0 ]]; then
  echo
  echo "  Skipping (already in tracking):"
  for ms in "${DONE[@]}"; do
    echo "    $(obsid "$ms")  ${ms##*/}"
  done
fi
if [[ ${#TODO[@]} -gt 0 ]]; then
  echo
  echo "  To run, in this order:"
  for ms in "${TODO[@]}"; do
    echo "    $(obsid "$ms")  ${ms##*/}"
  done
fi
echo "=================================================================="
echo

if [[ ${#TODO[@]} -eq 0 ]]; then
  echo "Nothing to do"
  exit 0
fi

while :; do
  read -r -p "Does this look right? Run the pipeline on these ${#TODO[@]} MS file(s)? [y/n] " answer
  case $answer in
    [Yy]|[Yy][Ee][Ss]) break ;;
    [Nn]|[Nn][Oo])     echo "Aborted"; exit 0 ;;
    *)                 echo "Please answer y or n" ;;
  esac
done

# ------------------------------------------------------------------ #
# Stage runners
# ------------------------------------------------------------------ #

# Source the submit script with sbatch shadowed by a function, so every job ID
# is recorded while the script's own `| awk '{print $4}'` still works.
submit_and_collect() {
  local script=$1 idfile=$2
  : > "$idfile"
  (
    # The submit scripts are written without set -u
    set +u
    sbatch() {
      local out
      out=$(command sbatch "$@") || { echo "$out" >&2; return 1; }
      awk '{print $NF}' <<< "$out" >> "$idfile"
      echo "$out"
    }
    source "./$script"
  )
}

# Wait until none of the job IDs in $1 (comma-separated) are queued, then
# check they all COMPLETED. One squeue call per poll, for this user's jobs only.
wait_for() {
  local ids=$1 q mine
  local pattern="^(${ids//,/|}) "
  while :; do
    if ! q=$(squeue -h -u "$USER" -o "%i %T %r" 2>&1); then
      log "squeue failed, retrying: $q"
      sleep "$POLL"
      continue
    fi
    mine=$(grep -E "$pattern" <<< "$q" || true)
    [[ -z $mine ]] && break
    if grep -q DependencyNeverSatisfied <<< "$mine"; then
      log "dependency never satisfied -- a job failed"
      return 1
    fi
    log "  $(grep -c RUNNING <<< "$mine") running, $(grep -c PENDING <<< "$mine") pending"
    sleep "$POLL"
  done
  local bad
  bad=$(sacct -j "$ids" -X -n -P -o JobID,State | grep -v '|COMPLETED$' || true)
  if [[ -n $bad ]]; then
    log "jobs not COMPLETED:"
    echo "$bad"
    return 1
  fi
  return 0
}

run_stage_slurm() {
  local subs=$1 idfile=$2 ids csv
  submit_and_collect "$subs" "$idfile"
  mapfile -t ids < "$idfile"
  [[ ${#ids[@]} -gt 0 ]] || { log "no jobs submitted"; return 1; }
  csv=$(IFS=,; echo "${ids[*]}")
  log "  submitted ${#ids[@]} job(s): $csv"
  if ! wait_for "$csv"; then
    scancel ${ids[@]} 2>/dev/null
    return 1
  fi
}

run_stage_node() {
  bash "./$1"
}

clear_work() {
  # WORK is a user-supplied global now, not a fixed path under this script, so
  # this is the actual safety gate before the recursive delete below.
  [[ $WORK != / ]] || { log "refusing to clear $WORK"; exit 1; }

  # Refuse if WORK is the same as, contains, or is inside PIPELINE or MS_DIR
  # -- wiping either would lose the pipeline checkout or the raw MS data.
  local other
  for other in "$PIPELINE" "$MS_DIR"; do
    case "$other/" in
      "$WORK"/*) log "refusing to clear $WORK: it contains $other"; exit 1 ;;
    esac
    case "$WORK/" in
      "$other"/*) log "refusing to clear $WORK: it is inside $other"; exit 1 ;;
    esac
  done

  # ROOT (this script's own directory) is only checked one way: WORK sitting
  # inside ROOT is the normal, documented layout ("place this script one
  # level above working_dir"), so only refuse the reverse -- WORK containing,
  # and so deleting, ROOT itself.
  case "$ROOT/" in
    "$WORK"/*) log "refusing to clear $WORK: it contains $ROOT"; exit 1 ;;
  esac

  mkdir -p "$WORK"
  find "$WORK" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
}

# Runs every stage for one MS inside working_dir; returns non-zero at the first
# stage that fails
run_ms() {
  local ms=$1 i setup subs stage
  local id=$(obsid "$ms")
  local state="$STATE/$id"
  # Wipe any state left over from a previous attempt at this obsid first, so
  # there's only ever one set of logs per MS -- a stale log from a stage that
  # isn't even reached this time (e.g. it now fails earlier) never lingers.
  rm -rf -- "$state"
  mkdir -p "$state"
  clear_work
  cp -a "$PIPELINE/." "$WORK/"
  ln -s "$(cd -- "$(dirname -- "$ms")" && pwd)/${ms##*/}" "$WORK/${ms##*/}"
  cd "$WORK" || return 1

  for ((i=0; i<${#STAGES[@]}; i++)); do
    setup=${STAGES[i]%%|*}
    subs=${STAGES[i]##*|}
    stage=${setup#setups/}
    stage=${stage%.py}
    log "  stage $i: $setup"
    python3 "$setup" "$INFRA" > "$state/setup_stage$i.log" 2>&1 || {
      log "  $setup failed (see $state/setup_stage$i.log)"
      record_failure "$id" "${ms##*/}" "$stage" "setup script failed, see $state/setup_stage$i.log"
      return 1
    }
    [[ -f ./$subs ]] || {
      log "  $subs not generated"
      record_failure "$id" "${ms##*/}" "$stage" "$subs not generated"
      return 1
    }
    if [[ $MODE == slurm ]]; then
      run_stage_slurm "$subs" "$state/stage$i.jobids" || {
        log "  stage $i failed"
        record_failure "$id" "${ms##*/}" "$stage" "job(s) did not complete"
        return 1
      }
    else
      run_stage_node "$subs" || log "  $subs exited non-zero, continuing"
    fi
  done
}

# ------------------------------------------------------------------ #
# Main loop
# ------------------------------------------------------------------ #

mkdir -p "$STATE"
PASSED=()
FAILED=()
PARTIAL=()
n=0

for ms in "${TODO[@]}"; do
  n=$((n+1))
  id=$(obsid "$ms")
  log "=== [$n/${#TODO[@]}] $id: ${ms##*/} ==="
  run_ms "$ms"
  ran_ok=$?
  cd "$ROOT" || exit 1
  # Pass/fail is still decided by whether TRACKING now has this obsid, since
  # only RMSYNTH's last step writes it -- but run_ms's own exit code is now
  # checked too, to tell "every stage that ran actually succeeded, RMSYNTH
  # just wasn't one of them" (STAGES was deliberately shortened, e.g. for a
  # quick test) apart from "a stage genuinely failed". Either way this loop
  # moves straight on to the next MS; there's no retry or abort here.
  if in_tracking "$id"; then
    log "=== $id complete ==="
    PASSED+=("$id")
  elif [[ $ran_ok -eq 0 && $STAGES_INCLUDE_RMSYNTH == false ]]; then
    log "=== $id stopped after the configured stage(s) -- STAGES doesn't include RMSYNTH, so it was never going to reach TRACKING. Not a failure. ==="
    PARTIAL+=("$id  ${ms##*/}")
  else
    log "=== $id FAILED (not in tracking) ==="
    FAILED+=("$id  ${ms##*/}")
  fi
done

echo
log "Finished: ${#PASSED[@]} completed, ${#FAILED[@]} failed, ${#PARTIAL[@]} stopped early (STAGES incomplete)"
for f in "${FAILED[@]}"; do
  echo "  FAILED   $f"
done
for f in "${PARTIAL[@]}"; do
  echo "  PARTIAL  $f  (ran the configured stage(s) OK; STAGES doesn't include RMSYNTH)"
done
