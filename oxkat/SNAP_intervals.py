#!/usr/bin/python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import os
import pickle
import re
import time
import sys
import os.path as o
import subprocess
import numpy as np
from collections import Counter

from astropy.io import fits as pyfits
from pyrap.tables import table

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg

# Comma-separated list of scan indices (0-based, in MS scan order) to (re)image
# this run. Leave as '' to process every scan found in the MS. Set to '-1' to skip
# imaging entirely and only run the completion-check/rename pass below.
# Entries can be single indices or CASA-style ranges:
#     '2,3'   -> scans 2 and 3
#     '2~5'   -> scans 2 through 5 inclusive
#     '2~'    -> scan 2 through the last scan
#     '~3'    -> scan 0 through scan 3 inclusive
# (see the WARNING printed at the end of a run for the exact list of scans that need re-running)
#
# Imaging is the expensive part of this script -- every scan is checked against what's
# already on disk (finalized/renamed, or imaged but not yet renamed) before considering
# SNAP_SCANS at all, so a scan that's already done is never re-imaged even if listed here,
# and any leftover unrenamed scan from an earlier run gets picked up and renamed automatically.
SNAP_SCANS = ''

def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt, flush=True)


def parse_scan_spec(spec, n_scans):
    scans = set()
    for token in spec.split(','):
        token = token.strip()
        if '~' in token:
            lo_str, hi_str = token.split('~', 1)
            lo = int(lo_str) if lo_str.strip() != '' else 0
            hi = int(hi_str) if hi_str.strip() != '' else n_scans - 1
            scans.update(range(lo, hi + 1))
        else:
            scans.add(int(token))
    return sorted(scans)


def trim_mask(mask_path, imsize, output_path):
    with pyfits.open(mask_path) as hdul:
        data = hdul[0].data.copy()
        header = hdul[0].header.copy()
    ny, nx = data.shape[-2], data.shape[-1]
    if ny == imsize and nx == imsize:
        return mask_path
    cy, cx = ny // 2, nx // 2
    half = imsize // 2
    data = data[..., cy - half:cy + half, cx - half:cx + half]
    if 'CRPIX1' in header:
        header['CRPIX1'] = header['CRPIX1'] - (cx - half)
    if 'CRPIX2' in header:
        header['CRPIX2'] = header['CRPIX2'] - (cy - half)
    pyfits.writeto(output_path, data, header, overwrite=True)
    return output_path


def get_total_ints(myms):

    tt = table(myms,readonly=True)

    scan_numbers = list(set(tt.getcol('SCAN_NUMBER')))
    exposure = round(np.mean(tt.getcol('EXPOSURE')),4)

    int_arr = []

    for scan in scan_numbers:
        subtab = tt.query(query='SCAN_NUMBER=='+str(scan))
        times = subtab.getcol('TIME')
        n_int = round((times[-1] - times[0]) / exposure, 0)
        int_arr.append(int(n_int))

    return int_arr

def main():
    
    intbin    = cfg.SNAP_INTBIN
    chanout = cfg.SNAP_CHANNELSOUT
    imsize   = cfg.SNAP_IMSIZE

    # Stokes choices
    pol = 'I'
    if cfg.SNAP_POL:
        pol = 'IQUV'

    # Deconvolution choice -- don't do this unless absolutely necessary
    niter = 0
    if cfg.SNAP_DECONV:
        niter = 100_000

    if len(sys.argv) == 1:
        print('Please specify an MS file for interval imaging')
        sys.exit()
    else:
        myms = sys.argv[1]

    mask = False
    if cfg.SNAP_DECONVMASK:
        msg(f'Deconvolution mask set explicitly: {cfg.SNAP_DECONVMASK}')
        mask = cfg.SNAP_DECONVMASK
    elif cfg.SNAP_DECONV:
        candidate = None
        if cfg.WSC_MASK != '' and cfg.WSC_MASK is not False and o.exists(cfg.WSC_MASK):
            msg(f'Mask found from WSC_MASK config: {cfg.WSC_MASK}')
            candidate = cfg.WSC_MASK
        else:
            for pattern in [
                cfg.IMAGES + f'/img_{myms}_snapblind-MFS-image.mask.fits',
                cfg.IMAGES + f'/img_{myms}_datablind-MFS-image.mask.fits',
            ]:
                if o.exists(pattern):
                    msg(f'Mask found from pattern: {pattern}')
                    candidate = pattern
                    break
            if candidate is None:
                msg('No target-specific mask found, falling back to glob search')
                found = sorted(glob.glob(cfg.IMAGES + '/img_*datablind*.mask.fits'))
                if not found:
                    found = sorted(glob.glob(cfg.IMAGES + '/img_*.mask.fits'))
                if found:
                    msg(f'Mask found via glob: {found[0]}')
                    candidate = found[0]

        if candidate:
            trimmed = cfg.IMAGES + f'/img_{myms}_snap_deconv_mask.fits'
            with pyfits.open(candidate) as hdul:
                ny, nx = hdul[0].data.shape[-2], hdul[0].data.shape[-1]
            if ny != imsize or nx != imsize:
                msg(f'Mask size {nx}x{ny}px does not match imsize {imsize}px — trimming')
                mask = trim_mask(candidate, imsize, trimmed)
                msg(f'Trimmed mask written to: {mask}')
            else:
                msg(f'Mask already matches imsize {imsize}px — no trimming needed')
                mask = candidate
        else:
            msg('WARNING: SNAP_DECONV is True but no mask found in IMAGES — running without mask')

    # Get number of integrations in each scan
    int_arr = get_total_ints(myms)

    # Interval boundaries for every scan, independent of SNAP_SCANS, so t-label numbering
    # is identical whether this is a full run or a re-run of only some scans.
    int0_arr = []
    int1_arr = []
    nint_arr = []
    for scan_n, ints in enumerate(int_arr):
        int0 = int(np.sum(int_arr[:scan_n]))
        if cfg.SNAP_INTEND:
            int1 = ints - ints % intbin + int0
            nint = int((int1 - int0) / intbin)
        else:
            int1 = int0 + ints
            nint = int(np.ceil((int1 - int0) / intbin))
        int0_arr.append(int0)
        int1_arr.append(int1)
        nint_arr.append(nint)

    # Global t-label offset each scan's intervals land at after renaming (cumulative nint of
    # every earlier scan).
    t_offset_arr = []
    cum = 0
    for n in nint_arr:
        t_offset_arr.append(cum)
        cum += n
    expected_total_nint = cum

    # Predicted t-label range each scan will occupy after renaming, printed up front so
    # it's clear what the final naming should look like before any imaging happens.
    for scan_n, (nint, t0) in enumerate(zip(nint_arr, t_offset_arr)):
        t_end = t0 + nint - 1
        msg(f'Predicted: scan{scan_n:04d} -> {nint} interval(s) -> after renaming: t{t0:04d}-t{t_end:04d}')
    msg(f'Predicted total: {expected_total_nint} interval(s) across {len(int_arr)} scan(s) -> t0000-t{expected_total_nint - 1:04d}')

    SNAP_SCANS_stripped = SNAP_SCANS.strip()
    if SNAP_SCANS_stripped == '-1':
        to_image = []
        msg('SNAP_SCANS = -1 -- skipping all imaging, only checking/renaming existing output')
    elif SNAP_SCANS_stripped != '':
        to_image = parse_scan_spec(SNAP_SCANS_stripped, len(int_arr))
    else:
        to_image = list(range(len(int_arr)))

    msg(f'Imaging {len(to_image)}/{len(int_arr)} scan(s): {to_image}')

    # Regex used to count distinct per-channel (non-MFS) image files and pull out their
    # t-index, e.g. '...-t0000-0005-image.fits' or '...-t0000-0005-I-image.fits'
    # -> interval '0000', channel '0005'
    img_re = re.compile(r'-t(\d{4})-(\d{4})-(?:[IQUV]-)?image\.fits$')

    # Lists to store rename pairs and bookkeeping for scans that don't come out clean
    rename_pairs = []
    scan_images = {}
    imaging_failed_scans = []

    # Iterate through the scans explicitly requested for (re)imaging -- always (re)imaged
    # regardless of whether output already exists for them.
    for run_idx, scan_n in enumerate(to_image, start=1):

        int0 = int0_arr[scan_n]
        int1 = int1_arr[scan_n]
        nint = nint_arr[scan_n]

        msg(f'Imaging scan {scan_n:04d} ({run_idx}/{len(to_image)})')

        # Temporary name for scan image
        image_prefix = cfg.INTERVALS+f'/img_{myms}_modelsub_scan{scan_n:04d}'

        # Generate imaging call and run
        imcall = gen.generate_syscall_wsclean(mslist = [myms],
                        imgname = image_prefix,
                        datacol = 'DATA',
                        chanout = chanout,
                        imsize = imsize,
                        niter = niter,
                        nomodel = True,
                        nodirty = True,
                        makepsf = True,
                        pol = pol,
                        intervalsout = nint,
                        interval0 = int0,
                        interval1 = int1,
                        field='0',
                        mask = mask,
                        chandeconvolution = 0,  # never do reduced-channel joint deconvolution for snapshots
                        automask = 6.0,
                        autothreshold = 1.0)

        wsclean_failed = False
        for syscall in imcall:
            result = subprocess.run([syscall], shell=True)
            if result.returncode != 0:
                msg(f'ERROR: wsclean exited with code {result.returncode} for scan {scan_n:04d} '
                    f'(signal {-result.returncode} if negative)')
                wsclean_failed = True

        # wsclean leaves '*-tmp.fits' files behind (e.g. a per-channel model file) when it
        # crashes partway through a channel/Stokes -- their presence means the scan is incomplete
        # even if the process's own exit code came back 0.
        tmp_files = sorted(glob.glob(f'{image_prefix}*tmp*'))
        if tmp_files:
            msg(f'ERROR: Scan {scan_n:04d} left {len(tmp_files)} temp file(s) behind -- wsclean did not finish cleanly:')
            for f in tmp_files:
                msg(f'  {f}')
            wsclean_failed = True

        if wsclean_failed:
            imaging_failed_scans.append(scan_n)
            continue

        # Fix channel naming (necessary when breaking up frequency images due to memory constraints)
        syscall = f'python3 {cfg.TOOLS}/fix_image_naming.py {chanout} {image_prefix}'
        subprocess.run([syscall], shell=True)

        images = sorted(glob.glob(f'{image_prefix}*'))

        if not images:
            msg(f'ERROR: Scan {scan_n:04d} produced no output images -- imaging failed, skipping')
            imaging_failed_scans.append(scan_n)
            continue

        scan_images[scan_n] = images

    if imaging_failed_scans:
        failed_str = ','.join(str(s) for s in imaging_failed_scans)
        msg(f'WARNING: {len(imaging_failed_scans)} scan(s) produced no images and were skipped: '
            f'{[f"{s:04d}" for s in imaging_failed_scans]}')
        msg(f'To re-run only the failed scans, edit this file, set SNAP_SCANS = \'{failed_str}\', and re-run.')

    # Always sweep every scan number for leftover '_scanNNNN'-tagged output -- whether it was
    # just imaged above or is left over from an earlier run -- so the check/rename pass below
    # covers the whole directory every time, not just what this invocation imaged.
    for scan_n in range(len(int_arr)):
        if scan_n in scan_images:
            continue
        image_prefix = cfg.INTERVALS+f'/img_{myms}_modelsub_scan{scan_n:04d}'
        images = sorted(glob.glob(f'{image_prefix}*'))
        if images:
            scan_images[scan_n] = images

    processed_scans = sorted(scan_images)

    if not processed_scans:
        msg('No scans have any un-renamed output to check/rename.')
        return

    # A scan that silently dropped channels or intervals (e.g. a partial/OOM wsclean failure)
    # still produces some images, so check the counts match before renaming anything.
    msg('Validating channel and interval counts across scans pending rename')
    n_channels_by_scan = {}
    n_intervals_by_scan = {}
    for scan_n in processed_scans:
        intervals_seen = set()
        channels_seen = set()
        for image in scan_images[scan_n]:
            match = img_re.search(image)
            if match:
                intervals_seen.add(match.group(1))
                channels_seen.add(match.group(2))
        n_channels_by_scan[scan_n] = len(channels_seen)
        n_intervals_by_scan[scan_n] = len(intervals_seen)

    # All scans image with the same chanout, so a successful scan's channel count should match
    # every other successful scan's -- use the most common value as the reference.
    expected_channels = Counter(n_channels_by_scan.values()).most_common(1)[0][0]

    bad_scans = []
    for scan_n in processed_scans:
        n_channels = n_channels_by_scan[scan_n]
        n_intervals = n_intervals_by_scan[scan_n]
        expected_nint = nint_arr[scan_n]

        problems = []
        if n_channels != expected_channels:
            problems.append(f'{n_channels} channel(s), expected {expected_channels} (to match other scans)')
        if n_intervals != expected_nint:
            problems.append(f'{n_intervals} interval(s), expected {expected_nint}')

        if problems:
            msg(f'ERROR: Scan {scan_n:04d} failed validation -- ' + '; '.join(problems))
            bad_scans.append(scan_n)
        else:
            msg(f'Scan {scan_n:04d} OK -- {n_channels} channels x {n_intervals} intervals')

    if bad_scans:
        msg(f'Channel/interval count validation failed for scan(s): {[f"{s:04d}" for s in bad_scans]}')
        msg('Rename operations aborted.')
        sys.exit(1)

    # Build rename pairs for every scan pending rename (freshly imaged or left over from an
    # earlier run) -- 'scan' naming obeys the 't0...' wsclean naming convention -- again has to
    # be a bettter way to do this but this'll work. Basic idea is to add the scan number to the
    # t-label; e.g., if scan 0 has 100 integrations: scan0001-t0000 becomes t100 (easier
    # indexing imo and avoids scan boundaries without the need to split MS files)
    for scan_n in processed_scans:
        image_prefix = cfg.INTERVALS+f'/img_{myms}_modelsub_scan{scan_n:04d}'
        image_prefix_fix = image_prefix.replace(f'_scan{scan_n:04d}','') # remove scan string
        for image in scan_images[scan_n]:
            image_fix = image.replace(image_prefix, image_prefix_fix) # replace it in image name
            marker = image_prefix_fix + '-t'
            if marker not in image_fix: # not a normal t-labelled output (tmp files are already caught above)
                msg(f'WARNING: Skipping unexpected file for scan {scan_n:04d} (no t-label found): {image}')
                continue
            suffix = image_fix.split(marker)[-1] # get suffix with improper t-label number
            suffix_fix = suffix.split('-')
            suffix_fix[0] = '{:04d}'.format(int(suffix_fix[0]) + t_offset_arr[scan_n])
            image_fix = image_fix.replace(suffix, '-'.join(suffix_fix)) # replace t-label in image
            rename_pairs.append((image, image_fix))

    # Validate and perform renames
    msg(f'Validating {len(rename_pairs)} image rename operations')

    # Check for 1-to-1 mapping (no duplicate target names)
    target_names = [pair[1] for pair in rename_pairs]
    if len(target_names) != len(set(target_names)):
        msg('ERROR: Duplicate target filenames detected in rename operations!')
        msg('Rename operations aborted to prevent file overwrites.')
        sys.exit(1)

    # Check that all source files exist
    for old_path, new_path in rename_pairs:
        if not os.path.exists(old_path):
            msg(f'ERROR: Source file does not exist: {old_path}')
            msg('Rename operations aborted.')
            sys.exit(1)

    # All checks passed, perform the renames (source/dest are always on the same INTERVALS
    # directory, so a plain rename is safe)
    total = len(rename_pairs)
    msg(f'Removing "scan" label in image names to standardize the snapshots ({total} files)')
    report_every = max(1, total // 10)
    for i, (old_path, new_path) in enumerate(rename_pairs, start=1):
        os.rename(old_path, new_path)
        if i % report_every == 0 or i == total:
            msg(f'Renamed {i}/{total} images ({100*i//total}%)')

    if imaging_failed_scans:
        msg(f'Successfully renamed {total} images ({len(imaging_failed_scans)} scan(s) skipped due to failures)')
    else:
        msg(f'Successfully renamed {total} images')

    # Final sanity check: total number of distinct t-labels actually present after
    # renaming should match the predicted total computed up front, regardless of
    # whether this run did the imaging/renaming itself or it was already done in
    # an earlier run.
    final_prefix = cfg.INTERVALS + f'/img_{myms}_modelsub'
    final_t_labels = set()
    for image in glob.glob(f'{final_prefix}-t*'):
        match = img_re.search(image)
        if match:
            final_t_labels.add(match.group(1))
    actual_total_nint = len(final_t_labels)

    if actual_total_nint != expected_total_nint:
        msg(f'WARNING: Something has gone wrong -- found {actual_total_nint} total t-interval(s) '
            f'after renaming, expected {expected_total_nint}. Check for scans that failed to '
            f'rename or produced the wrong total number of t values.')
    else:
        msg(f'Final check OK: {actual_total_nint} total t-interval(s) match prediction.')



if __name__ == "__main__":
    main()
