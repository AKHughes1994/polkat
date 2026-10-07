#!/usr/bin/env bash
#
# run_all_ms.sh -- run the full pipeline (0_GET_INFO -> 1GC -> 2GC -> RMSYNTH)
# on every *.ms in a directory, one MS at a time (or, on slurm, BATCH at a time).
#
# Fill in the INPUTS globals below, then run:
#     ./run_all_ms.sh 2>&1 | tee run_all_ms.log
#
# Place this script one level above working_dir. Each MS gets its own
# working_dir/<obsid>, emptied and filled with a copy of the pipeline directory
# and a symlink to the MS, and the stages are run inside it. A stage only
# starts once every job of the previous one has finished (for that MS). An MS counts as done once its obsid (leading
# number of the MS name) is in the tracking file, which RMSYNTH writes as its
# very last step. An MS that fails is reported and the run moves on to the next.
#
# MSs whose obsid is already in the tracking file are skipped. Before starting,
# it lists what it found and asks for confirmation.
#
# START_AT (an obsid or MS name) starts the run at that MS and leaves out every
# MS before it in the sorted list; with FIRST_ONLY=true that MS is the only one
# run, which is the way to debug one particular MS.
#
# MS_DIR is rescanned before every MS, so MSs that are still being transferred
# in when the run starts are picked up once they finish. An MS with any file
# modified in the last SETTLE seconds counts as still being written and is not
# started; if nothing else is ready the script waits for it.
#
# On slurm the job queue is checked every POLL seconds.
#
# BATCH is how many MSs run at once (slurm only). Each runs through its own
# stages independently, and as soon as one finishes (or fails) the next waiting
# MS takes its place, so there are always up to BATCH in flight.
#
# RESTART=true deletes the tracking, failures and interesting files, the
# contents of the archive and residual-MS directories and this script's logs,
# so every MS runs again. It asks twice (the second time "ARE YOU ABSOLUTELY
# SURE") and only deletes once you have also accepted the run plan.
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
#      is in the tracking file), "still being written" (modified in the last
#      SETTLE seconds) and "to run", prints the plan, and asks for
#      confirmation before doing anything.
#   4. Keeps up to BATCH MSs running at once. Whenever a slot is free, rescans
#      ms_dir and starts the next MS still to run; each MS, independently:
#        a. Empties working_dir/<obsid>, copies a fresh copy of the pipeline
#           into it, and symlinks that MS in.
#        b. Runs each pipeline stage in turn (0_GET_INFO, 1GC, 2GC, RMSYNTH):
#           generates that stage's job scripts, submits them (or runs them
#           directly, in "node" mode), and waits for every job to finish
#           before starting the next stage.
#        c. Stops at the first stage that fails for that MS, records which
#           stage that was in the failures file, and moves on to the next
#           MS -- one bad MS never blocks the others.
#   5. When an MS finishes, checks the tracking file again: if RMSYNTH's last step
#      added this obsid, that MS is counted as passed. If not, it's counted as
#      failed -- unless STAGES was deliberately shortened to stop before
#      RMSYNTH (e.g. for a quick test) and every stage that did run actually
#      succeeded, in which case it's counted as PARTIAL, not failed.
#   6. Once nothing is left to run or wait for, prints how many passed, failed, and
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
WORK='/mnt/scratchhdd/tkat_reprocessing/working_dir'       # scratch working directory, emptied at the start; each MS runs in WORK/<obsid>, rebuilt for every MS -- must be dedicated to this script, not shared with anything else
POLL=300            # slurm queue poll interval in seconds; also how often the main loop wakes
FIRST_ONLY=true     # true: run only the first MS that is not already in the tracking file (debugging)
START_AT='1549688122_2026-07-10T12-10-15_1Lf.ms'         # obsid or MS name: start the run at this MS, skipping every MS before it in the sorted list. With FIRST_ONLY=true, only this MS runs. Must not already be in the tracking file. '' = no start MS
BATCH=4             # slurm only: how many MSs run at once; a free slot is refilled as soon as an MS finishes. Ignored (1) on 'node' and when FIRST_ONLY is true
SETTLE=600          # seconds: an MS with any file modified in the last SETTLE seconds is treated as still transferring and not started
RESTART=false       # true: DELETE the tracking, failures and interesting files, the archive and residual-MS directories' contents and this script's logs, so every MS runs again. Asks twice before anything is deleted

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

[[ $BATCH =~ ^[1-9][0-9]*$ ]] || { echo "BATCH must be a positive integer (got: '$BATCH')" >&2; exit 1; }
# Batching is only for slurm (INFRA 'idia' or 'hippo'); node mode runs one MS at a time
[[ $MODE == node ]] && BATCH=1

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

# LOG_TAG is set to the obsid inside each background MS worker, so interleaved
# lines from MSs running at once can be told apart.
LOG_TAG=
log() { echo "[$(date '+%F %T')]${LOG_TAG:+ [$LOG_TAG]} $*"; }

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
  [[ $WIPE_PENDING == true ]] && return 1
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
  # Workers running at once can fail together, so create the header and append
  # under a lock.
  (
    flock 9
    [[ -f $FAILED_FILE ]] || echo '# obsid | ms_name | failed_utc | stage | reason' > "$FAILED_FILE"
    echo "$id | $msname | $(date -u '+%Y-%m-%dT%H:%M:%S') | $stage | $reason" >> "$FAILED_FILE"
  ) 9>> "$FAILED_FILE.lock"
}

# ------------------------------------------------------------------ #
# RESTART -- wipe the records and products of every earlier run
#
# With RESTART=true this script, once the plan below has been confirmed,
# deletes: the tracking file, the failures file (and its lock), the interesting
# file, everything inside the archive directory (RMSYNTH_ARCHIVE_DIR) and the
# residual-MS directory (RMSYNTH_RESIDUAL_MS_DIR) from config.py, and this
# script's own logs. It never touches MS_DIR, the pipeline or WORK's parents.
# Two confirmations are asked for before the run starts, and until the wipe
# happens every MS is treated as not done, so the plan shown is the real one.
# ------------------------------------------------------------------ #

[[ $RESTART == true || $RESTART == false ]] || { echo "RESTART must be 'true' or 'false' (got: '$RESTART')" >&2; exit 1; }

# While true, in_tracking reports nothing as done
WIPE_PENDING=false
WIPE_FILES=()
WIPE_DIRS=()

# Refuse a path that is not absolute and at least 3 levels deep, or that is, contains or sits
# inside anything this script must not lose
check_wipe_path() {
  local path=$1 other
  if [[ ! $path =~ ^/[^/]+/[^/]+/. ]]; then
    echo "RESTART: refusing to wipe '$path' (not an absolute path at least 3 levels deep)" >&2
    exit 1
  fi
  case "$(canon "$HOME")/" in
    "$path"/*) echo "RESTART: refusing to wipe $path: it contains $HOME" >&2; exit 1 ;;
  esac
  for other in "$MS_DIR" "$PIPELINE" "$WORK" "$ROOT"; do
    other=$(canon "$other")
    case "$other/" in
      "$path"/*) echo "RESTART: refusing to wipe $path: it contains $other" >&2; exit 1 ;;
    esac
    case "$path/" in
      "$other"/*) echo "RESTART: refusing to wipe $path: it is inside $other" >&2; exit 1 ;;
    esac
  done
}

# y/N, default N
confirm() {
  local reply
  read -r -p "$1 [y/N] " reply
  [[ $reply == [Yy] || $reply == [Yy][Ee][Ss] ]]
}

if [[ $RESTART == true ]]; then
  CFG_INTERESTING_FILE=$(cfg_value RMSYNTH_INTERESTING_FILE)
  CFG_RESIDUAL_MS_DIR=$(cfg_value RMSYNTH_RESIDUAL_MS_DIR)

  WIPE_FILES=("$(canon "$TRACKING")" "$(canon "$FAILED_FILE")" "$(canon "$FAILED_FILE.lock")")
  [[ -n $CFG_INTERESTING_FILE ]] && WIPE_FILES+=("$(canon "$CFG_INTERESTING_FILE")")
  WIPE_DIRS=("$(canon "$CFG_ARCHIVE_DIR")")
  [[ -n $CFG_RESIDUAL_MS_DIR ]] && WIPE_DIRS+=("$(canon "$CFG_RESIDUAL_MS_DIR")")

  for path in "${WIPE_FILES[@]}" "${WIPE_DIRS[@]}"; do
    check_wipe_path "$path"
  done

  echo
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo " RESTART IS TRUE -- THIS WILL PERMANENTLY DELETE:"
  echo
  echo "  Files (tracking / failures / interesting):"
  for path in "${WIPE_FILES[@]}"; do
    if [[ -e $path ]]; then echo "    $path"; else echo "    $path  (does not exist)"; fi
  done
  echo
  echo "  Everything inside these directories (final products):"
  for path in "${WIPE_DIRS[@]}"; do
    if [[ -d $path ]]; then
      echo "    $path  ($(find "$path" -mindepth 1 -maxdepth 1 | wc -l) entries)"
    else
      echo "    $path  (does not exist)"
    fi
  done
  echo
  echo "  This script's own logs:"
  echo "    $STATE"
  echo
  echo "  Every MS will then be run again from scratch. $MS_DIR is NOT touched."
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo

  confirm "Delete all of the above?" || { echo "Aborted, nothing was deleted"; exit 0; }
  confirm "ARE YOU ABSOLUTELY SURE? This cannot be undone." || { echo "Aborted, nothing was deleted"; exit 0; }
  WIPE_PENDING=true
fi

# Does the deletion that RESTART asked for, then lets in_tracking read the (now empty) files again
wipe_now() {
  local path
  for path in "${WIPE_FILES[@]}"; do
    if [[ -e $path ]]; then
      log "RESTART: removing $path"
      rm -f -- "$path"
    fi
  done
  for path in "${WIPE_DIRS[@]}"; do
    if [[ -d $path ]]; then
      log "RESTART: emptying $path"
      find "$path" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +
    fi
  done
  log "RESTART: removing $STATE"
  rm -rf -- "$STATE"
  WIPE_PENDING=false
}

# ------------------------------------------------------------------ #
# What needs running
# ------------------------------------------------------------------ #

# True if anything inside MS $1 was modified within the last SETTLE seconds,
# i.e. it is probably still being transferred.
ms_busy() {
  [[ -n $(find -L "$1" -newermt "$SETTLE seconds ago" -print -quit 2>/dev/null) ]]
}

# MSs started (pass or fail) in this run, keyed by obsid, so a failed MS is
# not picked up again by the next scan.
declare -A ATTEMPTED=()

# Scans MS_DIR and sorts its *.ms into ALL_MS (everything found), DONE (obsid
# in the tracking file), BUSY (still being written), and TODO (the rest, minus
# anything already attempted in this run). With START_AT set, every MS before
# it in the sorted list is left out of all of those, and START_MS is set to
# the MS it matched (empty if none did). Called before the summary and again
# before every MS, so MSs that finish transferring mid-run are picked up.
scan_ms() {
  mapfile -t ALL_MS < <(find "$MS_DIR" -mindepth 1 -maxdepth 1 -name '*.ms' \( -type d -o -type l \) | sort)
  TODO=()
  DONE=()
  BUSY=()
  START_MS=
  local ms id reached=false
  [[ -n $START_AT ]] || reached=true
  for ms in "${ALL_MS[@]}"; do
    id=$(obsid "$ms")
    if [[ $reached == false ]]; then
      if [[ $id == "$START_AT" || ${ms##*/} == "$START_AT" || ${ms##*/} == "$START_AT.ms" ]]; then
        reached=true
        START_MS=$ms
      else
        continue
      fi
    fi
    if in_tracking "$id"; then
      DONE+=("$ms")
    elif [[ -n ${ATTEMPTED[$id]:-} ]]; then
      :
    elif ms_busy "$ms"; then
      BUSY+=("$ms")
    else
      TODO+=("$ms")
    fi
  done
}

scan_ms

if [[ ${#ALL_MS[@]} -eq 0 ]]; then
  echo "No *.ms found in $MS_DIR"
  exit 0
fi

[[ $FIRST_ONLY == true || $FIRST_ONLY == false ]] || { echo "FIRST_ONLY must be 'true' or 'false' (got: '$FIRST_ONLY')" >&2; exit 1; }

# START_AT must name an MS that can actually be run
if [[ -n $START_AT ]]; then
  [[ -n $START_MS ]] || { echo "START_AT '$START_AT' does not match any MS in $MS_DIR (use an obsid or an MS name)" >&2; exit 1; }
  if in_tracking "$(obsid "$START_MS")"; then
    echo "START_AT ${START_MS##*/} is already in the tracking file ($TRACKING) -- remove its line to run it again" >&2
    exit 1
  fi
  if ms_busy "$START_MS"; then
    echo "START_AT ${START_MS##*/} is still being written (a file changed in the last ${SETTLE}s)" >&2
    exit 1
  fi
fi

# With FIRST_ONLY, only the first MS to run is kept
[[ $FIRST_ONLY == true ]] && BATCH=1
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
echo "  Working dir    : $WORK (emptied at the start; each MS runs in $WORK/<obsid>)"
if [[ $BATCH -gt 1 ]]; then
  echo "  Batch size     : up to $BATCH MSs at once (a free slot is refilled as soon as an MS finishes)"
else
  echo "  Batch size     : 1 (one MS at a time)"
fi
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
echo "  Found ${#ALL_MS[@]} MS file(s): ${#DONE[@]} already in tracking, ${#BUSY[@]} still being written, ${#TODO[@]} to run"
echo "  MS_DIR is rescanned before each MS, so MSs that finish transferring during the run are picked up"
echo "  (an MS is 'still being written' if any of its files changed in the last ${SETTLE}s)"
if [[ $RESTART == true ]]; then
  echo "  RESTART is true: the files and directories listed above WILL BE DELETED when you confirm below,"
  echo "  and every MS is treated as not done (the counts here already assume that)"
fi
if [[ -n $START_AT ]]; then
  echo "  START_AT is set: starting at ${START_MS##*/}, MSs before it in the sorted list are left out"
fi
if [[ $FIRST_ONLY == true ]]; then
  echo "  FIRST_ONLY is true: only the first MS to run is kept ($HELD_BACK more not run)"
fi
if [[ ${#BUSY[@]} -gt 0 ]]; then
  echo
  echo "  Not started yet (still being written):"
  for ms in "${BUSY[@]}"; do
    echo "    $(obsid "$ms")  ${ms##*/}"
  done
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

if [[ ${#TODO[@]} -eq 0 && ${#BUSY[@]} -eq 0 ]]; then
  echo "Nothing to do"
  exit 0
fi

while :; do
  read -r -p "Does this look right? Run the pipeline on these ${#TODO[@]} MS file(s) (waiting for ${#BUSY[@]} still being written)? [y/n] " answer
  case $answer in
    [Yy]|[Yy][Ee][Ss]) break ;;
    [Nn]|[Nn][Oo])     echo "Aborted"; exit 0 ;;
    *)                 echo "Please answer y or n" ;;
  esac
done

[[ $WIPE_PENDING == true ]] && wipe_now

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
  local ids=$1 q mine status last_status=
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
    status="$(grep -c RUNNING <<< "$mine") running, $(grep -c PENDING <<< "$mine") pending"
    [[ $status == "$last_status" ]] || log "  $status"   # only when the counts change
    last_status=$status
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

# With an obsid argument, empties just WORK/<obsid> (other MSs may be running in
# their own subdirectories); with none, empties all of WORK.
clear_work() {
  local id=${1:-}
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
  if [[ -n $id ]]; then
    rm -rf -- "$WORK/$id"
    mkdir -p "$WORK/$id"
  else
    find "$WORK" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  fi
}

# Runs every stage for one MS inside WORK/<obsid>; returns non-zero at the first
# stage that fails. Runs in a background subshell, so its cd stays its own.
run_ms() {
  local ms=$1 i setup subs stage
  local id=$(obsid "$ms")
  local state="$STATE/$id"
  # Wipe any state left over from a previous attempt at this obsid first, so
  # there's only ever one set of logs per MS -- a stale log from a stage that
  # isn't even reached this time (e.g. it now fails earlier) never lingers.
  rm -rf -- "$state"
  mkdir -p "$state"
  clear_work "$id"
  cp -a "$PIPELINE/." "$WORK/$id/"
  ln -s "$(cd -- "$(dirname -- "$ms")" && pwd)/${ms##*/}" "$WORK/$id/${ms##*/}"
  cd "$WORK/$id" || return 1

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
    # Whatever a stage reported, 2GC must have left its images and residual MS
    # before RMSYNTH runs on them
    if [[ $stage == 2GC ]]; then
      python3 tools/check_outputs.py 2GC > "$state/check_stage$i.log" 2>&1 || {
        log "  $stage outputs missing:"
        sed 's/^/    /' "$state/check_stage$i.log"
        record_failure "$id" "${ms##*/}" "$stage" "outputs missing, see $state/check_stage$i.log"
        return 1
      }
    fi
  done
  return 0
}

# ------------------------------------------------------------------ #
# Main loop
# ------------------------------------------------------------------ #

mkdir -p "$STATE"
PASSED=()
FAILED=()
PARTIAL=()
n=0

# Workers still running: obsid -> background pid, and obsid -> MS path.
declare -A PID_OF=()
declare -A MS_RUNNING=()

# On Ctrl-C / kill, stop the workers too. Submitted slurm jobs are left alone.
trap 'kill "${PID_OF[@]}" 2>/dev/null; exit 130' INT TERM

clear_work

# Judge a finished worker. Pass/fail is still decided by whether TRACKING now
# has this obsid, since only RMSYNTH's last step writes it -- but the worker's
# own exit code ($2) is checked too, to tell "every stage that ran actually
# succeeded, RMSYNTH just wasn't one of them" (STAGES was deliberately
# shortened, e.g. for a quick test) apart from "a stage genuinely failed".
finish_ms() {
  local id=$1 ran_ok=$2 ms=${MS_RUNNING[$id]}
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
  unset "PID_OF[$id]" "MS_RUNNING[$id]"
}

# Up to BATCH MSs run at once, each in a background worker going through its
# own stages independently. Each pass of this loop reaps any finished worker,
# then -- if a slot is free -- rescans MS_DIR and starts the next MS. Scanning
# is throttled to once per POLL seconds unless a worker just finished, so a
# full wait doesn't walk every MS tree often. The loop itself wakes every POLL
# seconds, so a freed slot is refilled within POLL seconds.
TICK=$POLL
last_scan=0
last_busy_msg=
just_finished=true
while :; do
  # Reap finished workers
  just_finished=false
  for id in "${!PID_OF[@]}"; do
    if ! kill -0 "${PID_OF[$id]}" 2>/dev/null; then
      wait "${PID_OF[$id]}"
      finish_ms "$id" $?
      just_finished=true
    fi
  done

  # Fill free slots
  now=$(date +%s)
  if [[ ${#PID_OF[@]} -lt $BATCH && ! ( $FIRST_ONLY == true && $n -ge 1 ) ]] \
     && [[ $just_finished == true || $last_scan -eq 0 || $((now - last_scan)) -ge $POLL ]]; then
    last_scan=$now
    scan_ms
    for ms in "${TODO[@]}"; do
      [[ ${#PID_OF[@]} -lt $BATCH ]] || break
      [[ $FIRST_ONLY == true && $n -ge 1 ]] && break
      id=$(obsid "$ms")
      n=$((n+1))
      ATTEMPTED[$id]=1
      log "=== [$n] $id: ${ms##*/} (${#PID_OF[@]} already running, ${#TODO[@]} to run including this one) ==="
      ( LOG_TAG=$id; run_ms "$ms" ) &
      PID_OF[$id]=$!
      MS_RUNNING[$id]=$ms
    done
    if [[ ${#PID_OF[@]} -lt $BATCH && ${#BUSY[@]} -gt 0 ]]; then
      busy_msg="waiting for ${#BUSY[@]} MS(s) still being written: $(for b in "${BUSY[@]}"; do obsid "$b"; done | tr '\n' ' ')"
      [[ $busy_msg == "$last_busy_msg" ]] || log "$busy_msg"   # only when the set changes
      last_busy_msg=$busy_msg
    fi
  fi

  # Done when nothing is running and nothing is left to start or wait for
  if [[ ${#PID_OF[@]} -eq 0 ]]; then
    if [[ $FIRST_ONLY == true && $n -ge 1 ]] || [[ ${#TODO[@]} -eq 0 && ${#BUSY[@]} -eq 0 ]]; then
      break
    fi
  fi
  sleep "$TICK"
done

echo
log "Finished: ${#PASSED[@]} completed, ${#FAILED[@]} failed, ${#PARTIAL[@]} stopped early (STAGES incomplete)"
for f in "${FAILED[@]}"; do
  echo "  FAILED   $f"
done
for f in "${PARTIAL[@]}"; do
  echo "  PARTIAL  $f  (ran the configured stage(s) OK; STAGES doesn't include RMSYNTH)"
done
