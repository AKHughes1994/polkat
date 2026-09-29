#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk
#
# Final job of the RMSYNTH stage. In order:
#
#   1. Appends every interesting target epoch of this run to
#      cfg.RMSYNTH_INTERESTING_FILE. An epoch is interesting if it is
#        BRIGHT    -- MFS Stokes I > BRIGHT_I_MJY, and/or
#        POLARISED -- peak-channel P/I > POL_FRAC_PCT with
#                     (P/I) / err(P/I) >= POL_SIGMA, and/or
#        PARANG    -- parallactic angle range across the epoch's own scan(s)
#                     > PARANG_DELTA_DEG: a lot of parallactic angle is being
#                     averaged together
#      (peak channel and its P/I +/- error as in RMSYNTH_04_summarize_target.py).
#   2. Writes <obsid>_archive_summary.txt to the working directory (date, target,
#      MFS flux density, peak-channel P/I, interesting or not, per target epoch).
#   3. Moves everything in the working directory into
#      cfg.RMSYNTH_ARCHIVE_DIR/<obsid>/, deleting that directory first if it
#      already exists. LOGS goes last, as it holds this job's own log.
#   4. Appends one line for this obsid to cfg.RMSYNTH_TRACKING_FILE: completion
#      time and the interesting sources from step 1. An obsid missing from the
#      tracking file did not finish.
#
# Any of the three config paths left blank skips its step. The obsid is the
# leading number of the master MS name, e.g. 1549688122 for
# 1549688122_2026-07-10T12-10-15_1Lf_1024ch_GX_339-4.ms.
#
# Uses only the standard library, and is run with the host python3 rather than
# inside a container, so the archive/record paths need not be bind-mounted.
#
# Input files (RESULTS = cfg.RESULTS), the same ones RMSYNTH_04 summarizes:
#
#   project_info.json
#       master_ms (falls back to working_ms) -- obsid
#       target_names, working_names          -- targets of this run
#
#   RESULTS/*<target>*_pcalmask_polarization.json (one per epoch)
#       ['MFS']:
#           I_flux_mJy, I_err_mJy   -- MFS Stokes I +/- error
#           time_ctr_isot           -- exposure-centre date/time (UTC)
#       ['timing']['middle_mjd']    -- exposure-centre MJD, used when time_ctr_isot is absent
#       ['pol_image_type'], ['CHAN'] -- per-channel P/I, see peak_chan_frac_pol()
#       parang_delta_deg            -- parallactic angle range across the epoch's own
#                                      scan(s) (deg); absent for a polarization.json
#                                      written before 1GC_08_casa_split_targets.py
#                                      recorded it

import datetime
import fcntl
import glob
import json
import os
import os.path as o
import shutil
import sys

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg
from oxkat.RMSYNTH_04_summarize_target import IDENTIFIER, fmt, peak_chan_frac_pol

BRIGHT_I_MJY = 10.0   # MFS Stokes I above this (mJy) is BRIGHT
POL_FRAC_PCT = 1.0    # Peak-channel P/I above this (%) ...
POL_SIGMA    = 10.0   # ... at at least this significance, (P/I)/err(P/I), is POLARISED

# Parallactic angle range (parang_delta_deg, from 1GC_08_casa_split_targets.py
# via the target's polarization.json) above this (deg) is PARANG: a lot of
# parallactic angle is being averaged together across this epoch.
PARANG_DELTA_DEG = 10.0

INTERESTING_HEADER = ('# obsid | source | date_centre_utc | MFS_I_mJy | '
                      'peak_chan_P/I_pct | peak_chan_freq_GHz | parang_delta_deg | reason\n')
TRACKING_HEADER    = '# obsid | completed_utc | interesting sources (reason)\n'


def msg(text):
    print(f'[RMSYNTH_05] {text}', flush=True)


def now_utc():
    return datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%S')


def mjd_to_isot(mjd):
    return (datetime.datetime(1858, 11, 17) + datetime.timedelta(days=mjd)).isoformat(timespec='seconds')


def get_obsid(project_info):
    myms = project_info.get('master_ms') or project_info['working_ms']
    return o.basename(str(myms).rstrip('/')).split('_')[0]


def epoch_record(target, pol_json_path):
    """Values reported for one target epoch, plus the reasons it is interesting (if any)."""

    with open(pol_json_path) as f:
        pol = json.load(f)
    mfs = pol.get('MFS', {})

    date = mfs.get('time_ctr_isot')
    if date is None and pol.get('timing', {}).get('middle_mjd') is not None:
        date = mjd_to_isot(pol['timing']['middle_mjd'])

    I_flux = mfs.get('I_flux_mJy')
    peak = peak_chan_frac_pol(pol)
    parang_delta = mfs.get('parang_delta_deg')

    reasons = []
    if I_flux is not None and I_flux > BRIGHT_I_MJY:
        reasons.append('BRIGHT')
    if peak is not None:
        _, _, frac, frac_err = peak
        if frac > POL_FRAC_PCT and frac_err > 0 and frac / frac_err >= POL_SIGMA:
            reasons.append('POLARISED')
    if parang_delta is not None and parang_delta > PARANG_DELTA_DEG:
        reasons.append('PARANG')

    return {
        'target':       target,
        'epoch':        o.basename(pol_json_path).replace(f'_{IDENTIFIER}_polarization.json', ''),
        'date':         date,
        'I_flux':       I_flux,
        'I_err':        mfs.get('I_err_mJy'),
        'ptype':        pol.get('pol_image_type', 'P'),
        'peak':         peak,
        'parang_delta': parang_delta,
        'reasons':      reasons,
    }


def fmt_num(value, digits):
    return 'N/A' if value is None else f'{value:.{digits}f}'


def fmt_flux(rec):
    if rec['I_flux'] is None:
        return 'N/A'
    return f"{fmt_num(rec['I_flux'], 3)} +/- {fmt_num(rec['I_err'], 3)}"


def fmt_peak(rec):
    if rec['peak'] is None:
        return 'N/A', 'N/A'
    _, freq, frac, frac_err = rec['peak']
    return f"{frac:.3f} +/- {frac_err:.3f} ({rec['ptype']})", fmt_num(freq, 4)


def append_locked(path, header, lines):
    """Append lines to path, creating it (and its directory) with header if new."""
    parent = o.dirname(o.abspath(path))
    os.makedirs(parent, exist_ok=True)
    with open(path, 'a') as f:
        try:
            fcntl.lockf(f, fcntl.LOCK_EX)
            locked = True
        except OSError as e:
            msg(f'WARNING: could not lock {path} ({e}); appending without a lock')
            locked = False
        try:
            f.seek(0, os.SEEK_END)
            if f.tell() == 0:
                f.write(header)
            f.writelines(lines)
            f.flush()
        finally:
            if locked:
                fcntl.lockf(f, fcntl.LOCK_UN)


def write_interesting(obsid, records):
    lines = []
    for rec in records:
        if not rec['reasons']:
            continue
        peak_str, freq_str = fmt_peak(rec)
        lines.append(f"{obsid} | {rec['target']} | {rec['date'] or 'N/A'} | {fmt_flux(rec)} | "
                     f"{peak_str} | {freq_str} | {fmt_num(rec['parang_delta'], 2)} | "
                     f"{', '.join(rec['reasons'])}\n")
        if 'PARANG' in rec['reasons']:
            msg(f"WARNING: {rec['target']} ({rec['epoch']}) spans "
                f"{rec['parang_delta']:.2f} deg of parallactic angle -- a lot of "
                f"parallactic angle is being averaged together.")
    if lines:
        append_locked(cfg.RMSYNTH_INTERESTING_FILE, INTERESTING_HEADER, lines)
    msg(f'{len(lines)} interesting epoch(s) appended to {cfg.RMSYNTH_INTERESTING_FILE}')


def write_tracking(obsid, records):
    interesting = [f"{rec['target']} ({', '.join(rec['reasons'])})" for rec in records if rec['reasons']]
    line = f"{obsid} | {now_utc()} | {'; '.join(interesting) if interesting else 'none'}\n"
    append_locked(cfg.RMSYNTH_TRACKING_FILE, TRACKING_HEADER, [line])
    msg(f'Marked {obsid} complete in {cfg.RMSYNTH_TRACKING_FILE}')


def write_archive_summary(path, obsid, project_info, records, missing):

    with open(path, 'w') as out:

        def line(label, value):
            out.write(f'{label.ljust(32)}: {value}\n')

        out.write(f'Archive summary for obsid {obsid}\n')
        out.write('='*60+'\n')
        line('Master MS', project_info.get('master_ms', 'N/A'))
        line('Archived (UTC)', now_utc())
        line('Archived to', o.join(cfg.RMSYNTH_ARCHIVE_DIR, obsid))

        for rec in records:
            peak_str, freq_str = fmt_peak(rec)
            out.write(f"\n--- {rec['target']} ({rec['epoch']}) ---\n")
            line('Date/time, centre (UTC)', rec['date'] or 'N/A')
            line('MFS Stokes I (mJy)', fmt_flux(rec))
            line('Peak channel P/I (%)', peak_str)
            line('Peak channel frequency (GHz)', freq_str)
            line('Parallactic angle range (deg)', fmt_num(rec['parang_delta'], 2))
            line('Interesting', ', '.join(rec['reasons']) if rec['reasons'] else 'no')

        for target in missing:
            out.write(f'\n--- {target} ---\n')
            line('Results', f'no *{target}*_{IDENTIFIER}_polarization.json in RESULTS')

    msg(f'Wrote {path}')


def archive_working_dir(obsid):
    """
    Move every entry of the working directory into RMSYNTH_ARCHIVE_DIR/<obsid>/.
    Returns True once everything has been moved.
    """

    cwd  = o.realpath(cfg.CWD)
    dest = o.realpath(o.join(cfg.RMSYNTH_ARCHIVE_DIR, obsid))

    # Refuse destinations that would delete or recurse into the working directory
    if dest == cwd or dest.startswith(cwd+os.sep) or cwd.startswith(dest+os.sep):
        msg(f'ERROR: archive destination {dest} overlaps the working directory {cwd}; not archiving')
        return False

    for record_path in (cfg.RMSYNTH_INTERESTING_FILE, cfg.RMSYNTH_TRACKING_FILE):
        if record_path and o.realpath(record_path).startswith(cwd+os.sep):
            msg(f'WARNING: {record_path} is inside the working directory and will be moved with it')

    if o.islink(dest) or o.isfile(dest):
        os.remove(dest)
    elif o.isdir(dest):
        msg(f'Deleting existing {dest}')
        shutil.rmtree(dest)
    os.makedirs(dest)

    # LOGS holds this job's own log, so it goes last
    entries = sorted(os.listdir(cwd), key=lambda name: (name == 'LOGS', name))
    msg(f'Moving {len(entries)} entries from {cwd} to {dest}')
    for name in entries:
        shutil.move(o.join(cwd, name), o.join(dest, name))

    # A scheduler that delivers its log at job end (PBS) writes it here
    os.makedirs(o.join(cwd, 'LOGS'), exist_ok=True)
    msg('Archive complete')
    return True


def main():

    with open('project_info.json') as f:
        project_info = json.load(f)

    obsid = get_obsid(project_info)
    msg(f'obsid: {obsid}')

    target_names  = project_info.get('target_names', [])
    working_names = project_info.get('working_names', [])
    targets = [t for t in target_names if t in working_names]

    records = []
    missing = []
    for target in targets:
        pol_json_paths = sorted(glob.glob(o.join(cfg.RESULTS, f'*{target}*_{IDENTIFIER}_polarization.json')))
        if not pol_json_paths:
            missing.append(target)
        for pol_json_path in pol_json_paths:
            records.append(epoch_record(target, pol_json_path))

    if cfg.RMSYNTH_INTERESTING_FILE:
        write_interesting(obsid, records)

    if cfg.RMSYNTH_ARCHIVE_DIR:
        write_archive_summary(o.join(cfg.CWD, f'{obsid}_archive_summary.txt'),
                              obsid, project_info, records, missing)
        if not archive_working_dir(obsid):
            msg('Archive not completed; obsid not marked complete in the tracking file')
            sys.exit(1)

    # Written last, so an obsid is only marked complete once its archive is done
    if cfg.RMSYNTH_TRACKING_FILE:
        write_tracking(obsid, records)


if __name__ == "__main__":

    main()
