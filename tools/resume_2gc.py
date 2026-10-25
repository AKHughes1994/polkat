#!/usr/bin/env python3
"""
resume_2gc.py

Write submit_2GC_resume.sh: the 2GC jobs still to run after a failed 2GC
stage, so only those are resubmitted (slurm only).

Every 2GC job submitted so far is looked up in sacct, keeping the latest
attempt of each job name. Only source blocks of submit_2GC_jobs.sh with a
job that sacct reports as not COMPLETED are restarted, from that job
onwards; blocks whose jobs all completed, or that have no job in the
lookup, are left out. Dependencies on completed jobs are
dropped, since those jobs are not resubmitted. The jobs' slurm scripts in
SCRIPTS/ are reused unchanged.

Run from the MS's working directory:

    python3 tools/resume_2gc.py [--dry-run] <job id file> [<job id file> ...]

--dry-run reports what would be resubmitted without writing anything.
Exits 1, writing nothing, if every job completed.
"""

import re
import subprocess
import sys

SUBMIT = "submit_2GC_jobs.sh"
RESUME = "submit_2GC_resume.sh"

# One job: NAME=`sbatch [-d afterok:${A}:${B}] SCRIPTS/slurm_NAME.sh | awk '{print $4}'`
JOB = re.compile(r"^(\w+)=`sbatch (?:-d afterok:(\S+) )?(.*)`$")
# A source's kill-script line: echo "scancel "$A" "$B > kill_file
KILL = re.compile(r"^echo \"scancel \".* > (\S+)$")


def job_states(id_files):
    """
    Final state of each submitted job name, from sacct.

    Parameters
    ----------
    id_files : list of str
        Files of submitted slurm job IDs, one per line.

    Returns
    -------
    dict of {str: str}
        Job name to the state of its latest attempt (e.g. 'COMPLETED').
    """
    ids = [line.strip() for name in id_files for line in open(name) if line.strip()]
    out = subprocess.run(["sacct", "-j", ",".join(ids), "-X", "-n", "-P", "-o", "JobID,JobName,State"],
                         check=True, capture_output=True, text=True).stdout
    latest = {}
    for line in out.splitlines():
        jobid, name, state = line.split("|")
        if name not in latest or int(jobid) > latest[name][0]:
            latest[name] = (int(jobid), state.split()[0])  # 'CANCELLED by 123' -> 'CANCELLED'
    return {name: state for name, (_, state) in latest.items()}


def source_blocks(lines):
    """
    Split a submit script into its preamble and one block per source.

    Returns
    -------
    preamble : list of str
    blocks : list of (str, list of str)
        Source heading line and the block's lines after the heading.
    """
    preamble, blocks = [], []
    for i, line in enumerate(lines):
        current = blocks[-1][1] if blocks else preamble
        # Each source starts with a '# name (type)' line between two '#----'
        # rulers; the opening ruler and blank line belong to the heading
        if line.startswith("# ") and i > 0 and lines[i - 1].startswith("#----"):
            while current and (current[-1].startswith("#----") or not current[-1].strip()):
                current.pop()
            blocks.append((line, []))
        elif not (blocks and not current and line.startswith("#----")):
            current.append(line)
    return preamble, blocks


def main():
    dry_run = "--dry-run" in sys.argv[1:]
    states = job_states([a for a in sys.argv[1:] if a != "--dry-run"])
    preamble, blocks = source_blocks(open(SUBMIT).read().splitlines())
    out = list(preamble)

    for heading, lines in blocks:
        jobs = [(i, JOB.match(line)) for i, line in enumerate(lines) if JOB.match(line)]
        # First job sacct reports as not completed; blocks without one are left alone
        restart = next((k for k, (_, m) in enumerate(jobs) if states.get(m[1], "COMPLETED") != "COMPLETED"), None)
        if restart is None:
            if not any(m[1] in states for _, m in jobs):
                print(f"{heading[2:]}: no job in the sacct lookup, left alone")
            continue
        keep = {m[1] for _, m in jobs[restart:]}
        first = jobs[restart][1][1]
        print(f"{heading[2:]}: restarting at {first} ({states[first]}), {len(keep)} job(s)")

        out += ["", "#---------------------------------------", heading, "#---------------------------------------"]
        for i, line in enumerate(lines[jobs[restart][0] - 1:], start=jobs[restart][0] - 1):
            m, kill = JOB.match(line), KILL.match(line)
            if m:
                # Only wait on jobs that are resubmitted here
                deps = [d for d in re.findall(r"\$\{(\w+)\}", m[2] or "") if d in keep]
                dep = f"-d afterok:{':'.join('${' + d + '}' for d in deps)} " if deps else ""
                out.append(f"{m[1]}=`sbatch {dep}{m[3]}`")
            elif kill:
                names = [mm[1] for _, mm in jobs[restart:]]
                out.append('echo "scancel "$' + '" "$'.join(names) + f" > {kill[1]}")
            elif line.startswith("#") or not line.strip():
                out.append(line)

    if len(out) == len(preamble):
        sys.exit(f"Every 2GC job completed -- nothing to resubmit (see {SUBMIT})")
    if dry_run:
        return

    with open(RESUME, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"Wrote {RESUME}")


if __name__ == "__main__":
    main()
