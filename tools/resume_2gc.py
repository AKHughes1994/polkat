#!/usr/bin/env python3
"""
resume_2gc.py

Write submit_2GC_resume.sh: the 2GC jobs still to run after a failed 2GC
stage, so only those are resubmitted (slurm only).

Every 2GC job submitted so far is looked up in sacct, keeping the latest
attempt of each job name. Only source blocks of submit_2GC_jobs.sh with a
failed job are restarted, from the first failed job onwards; blocks whose
jobs all completed, or that have no job in the lookup, are left out. A job
has failed when sacct reports it as not COMPLETED, or when it is COMPLETED
but its outputs are missing: a job script carries on past a crashed command
and still exits 0, so slurm reports such a job as COMPLETED. The outputs
checked are
  - every image its wsclean calls should have written: <name>-MFS-image.fits,
    one per polarisation (<name>-MFS-Q-image.fits, ...) for more than one,
    with the -MFS- and -tNNNN parts only when -channels-out / -intervals-out
    are above 1. An image that a cleanup in the job itself or a later job of
    the same source deletes on purpose (find -delete, rm, clean_directory.py)
    is not expected, nor are channel-split (_partNNNN) runs, whose MFS images
    fix_image_naming.py removes;
  - the residual MS written by 2GC_casa_split_residuals.py.
Dependencies on completed jobs are dropped, since those jobs are not
resubmitted. The jobs' slurm scripts in SCRIPTS/ are reused unchanged.

Run from the MS's working directory:

    python3 tools/resume_2gc.py [--dry-run] <job id file> [<job id file> ...]

--dry-run reports what would be resubmitted without writing anything.
Exits 1, writing nothing, if every job completed.
"""

import fnmatch
import os
import re
import shlex
import subprocess
import sys

SUBMIT = "submit_2GC_jobs.sh"
RESUME = "submit_2GC_resume.sh"
CLEAN_SCRIPT = "tools/clean_directory.py"

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


def clean_patterns():
    """
    Image filename globs clean_directory.py deletes from a source's image
    directory, and the globs it keeps.

    Returns
    -------
    (list of str, list of str)
        Deleted and excluded globs, matched against file basenames.
    """
    if not os.path.isfile(CLEAN_SCRIPT):
        return [], []
    text = open(CLEAN_SCRIPT).read()
    deleted, _, excluded = text.partition("exclude_patterns = [")
    glob_re = re.compile(r'^\s*f"\{cfg\.IMAGES\}/\{source_name\}/([^"]+)"', re.M)
    return glob_re.findall(deleted), glob_re.findall(excluded)


def deleted_on_purpose(path, scripts, cleaned):
    """
    Whether a cleanup in `scripts` (find -delete, rm, clean_directory.py) removes `path`.
    """
    base = os.path.basename(path)
    for script in scripts:
        for line in script.splitlines():
            if "clean_directory.py" in line and "--ms-only" not in line:
                deleted, excluded = cleaned
                if (any(fnmatch.fnmatch(base, g) for g in deleted)
                        and not any(fnmatch.fnmatch(base, g) for g in excluded)):
                    return True
                continue
            try:
                tokens = shlex.split(line)
            except ValueError:
                continue
            if "find" in tokens and "-delete" in tokens:
                keep = [tokens[i + 1] for i in range(1, len(tokens) - 1)
                        if tokens[i] == "-name" and tokens[i - 1] == "!"]
                match = [tokens[i + 1] for i in range(len(tokens) - 1)
                         if tokens[i] == "-name" and tokens[i - 1] != "!"]
                if (all(fnmatch.fnmatch(base, g) for g in match)
                        and not any(fnmatch.fnmatch(base, g) for g in keep)):
                    return True
            elif tokens[:1] == ["rm"]:
                if any(fnmatch.fnmatch(path, g) or fnmatch.fnmatch(base, os.path.basename(g))
                       for g in tokens[1:] if not g.startswith("-")):
                    return True
    return False


def expected_outputs(script):
    """
    Images and the residual MS a job script should leave behind.

    Returns
    -------
    list of str
        One path per expected image (see the module docstring) or residual MS.
    """
    paths = []
    for line in script.splitlines():
        tokens = line.split()
        # ... -c <tools>/2GC_casa_split_residuals.py <ms> <residual ms>
        split = [i for i, t in enumerate(tokens) if t.endswith("2GC_casa_split_residuals.py")]
        if split and len(tokens) > split[0] + 2:
            paths.append(tokens[split[0] + 2])
            continue
        if "wsclean" not in tokens or "-predict" in tokens or "-name" not in tokens:
            continue
        opt = lambda flag, default: tokens[tokens.index(flag) + 1] if flag in tokens else default
        name = opt("-name", "")
        if re.search(r"_part\d{4}$", name):
            continue
        pol = opt("-pol", "I")
        nchan = int(opt("-channels-out", 1))
        nint = int(opt("-intervals-out", 1))
        mfs = "-MFS" if nchan > 1 else ""
        pols = [""] if len(pol) == 1 else [f"-{p}" for p in pol]
        times = [""] if nint <= 1 else [f"-t{t:04d}" for t in range(nint)]
        paths += [f"{name}{t}{mfs}{p}-image.fits" for t in times for p in pols]
    return paths


def missing_outputs(scripts, k, cleaned):
    """
    Outputs of job k (index into `scripts`, its source's job scripts in
    order) that are absent and not deleted on purpose by job k or a later job.
    """
    return [p for p in expected_outputs(scripts[k])
            if not os.path.exists(p) and not deleted_on_purpose(p, scripts[k:], cleaned)]


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
    cleaned = clean_patterns()
    out = list(preamble)

    for heading, lines in blocks:
        jobs = [(i, JOB.match(line)) for i, line in enumerate(lines) if JOB.match(line)]
        scripts = []
        for _, m in jobs:
            script = m[3].split()[0]
            scripts.append(open(script).read() if os.path.isfile(script) else "")

        # First failed job: not completed in sacct, or completed with outputs
        # missing. Blocks without one are left alone.
        restart, reason = None, ""
        for k, (_, m) in enumerate(jobs):
            state = states.get(m[1], "COMPLETED")
            if state != "COMPLETED":
                restart, reason = k, state
                break
            if m[1] in states:
                missing = missing_outputs(scripts, k, cleaned)
                if missing:
                    restart = k
                    reason = (f"COMPLETED but {len(missing)} output(s) missing: "
                              + ", ".join(os.path.basename(p) for p in missing))
                    break
        if restart is None:
            if not any(m[1] in states for _, m in jobs):
                print(f"{heading[2:]}: no job in the sacct lookup, left alone")
            continue
        keep = {m[1] for _, m in jobs[restart:]}
        first = jobs[restart][1][1]
        print(f"{heading[2:]}: restarting at {first} ({reason}), {len(keep)} job(s)")

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
