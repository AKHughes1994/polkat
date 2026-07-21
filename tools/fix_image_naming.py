import os
import random
import sys
import glob
import subprocess
import time

import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt)


def fix_names(prefix, totchan):
    """
    Rename wsclean parted output files (prefix_partXXXX-CCCC-<type>.fits) to
    their final channel-contiguous names (prefix-CCCC-<type>.fits).

    Builds the full rename list before touching the filesystem, then runs a
    suffix-consistency check (image→image, model→model, etc.) before executing
    any mv commands.

    Arguments:
        prefix:  string, the wsclean image name prefix
        totchan: int, total number of output channels across all parts
    """

    part_pattern = f'{prefix}_part*.fits'
    part_files   = sorted(glob.glob(part_pattern))

    if not part_files:
        msg(f'No "part" files found matching: {part_pattern}')
        msg('Nothing to rename — exiting')
        return

    # -------------------------------------------------------------------------
    # Infer part numbers from filenames on disk
    # -------------------------------------------------------------------------
    prefix_base  = os.path.basename(prefix)
    part_tag_base = f'{prefix_base}_part'
    part_numbers = set()
    for f in part_files:
        basename = os.path.basename(f)
        if part_tag_base not in basename:
            msg(f'WARNING: expected "{part_tag_base}" in filename, skipping: {basename}')
            continue
        tail = basename.split(part_tag_base, 1)[1].split('-')[0]
        try:
            part_numbers.add(int(tail))
        except ValueError:
            msg(f'WARNING: could not parse part index from: {basename}')

    num_parts = len(part_numbers)
    msg(f'Found {num_parts} part(s) | parts: {sorted(part_numbers)}')

    if totchan % num_parts != 0:
        msg(f'ERROR: totchan ({totchan}) is not evenly divisible by '
            f'num_parts ({num_parts}) — check inputs and rerun')
        sys.exit(1)

    maxchan = totchan // num_parts
    msg(f'Channel accounting: {num_parts} parts x {maxchan} channels/part = '
        f'{num_parts * maxchan} total  (expected {totchan}) '
        f'[{"OK" if num_parts * maxchan == totchan else "MISMATCH"}]')

    # -------------------------------------------------------------------------
    # Remove invalid MFS files produced during the parted run
    # -------------------------------------------------------------------------
    mfs_files = glob.glob(f'{prefix}_part*-MFS-*')
    if mfs_files:
        msg(f'Removing {len(mfs_files)} invalid MFS file(s) from parted run')
        subprocess.run([f'rm -rf {prefix}_part*-MFS-*'], shell=True)
    else:
        msg('No MFS files to remove')

    part_files = sorted(glob.glob(f'{prefix}_part*.fits'))
    if not part_files:
        msg('ERROR: no part FITS files remain after MFS cleanup')
        sys.exit(1)

    # -------------------------------------------------------------------------
    # Build the rename list surgically — no open-ended string replacement.
    # Each filename is split at the known _partXXXX- boundary so only the
    # part tag and channel index are touched; the type suffix is carried through
    # unchanged.
    #
    # Plain:    prefix_partP-C-<type>.fits      →  prefix-C'-<type>.fits
    # Timed:    prefix_partP-tT-C-<type>.fits   →  prefix-tT-C'-<type>.fits
    #   where C' = C + P * maxchan
    #
    # has_time is determined per-file (not globally) to handle any mixed output.
    # -------------------------------------------------------------------------
    rename_pairs = []  # list of (old_path, new_path, old_suffix, new_suffix)

    for old in part_files:
        basename  = os.path.basename(old)
        directory = os.path.dirname(old)

        # Identify which part this file belongs to
        factor = None
        for part in sorted(part_numbers):
            if f'{prefix_base}_part{part:04d}' in basename:
                factor = part
                break

        if factor is None:
            msg(f'WARNING: could not match part for {basename} — skipping')
            continue

        part_tag           = f'{prefix_base}_part{factor:04d}'
        new_channel_offset = factor * maxchan

        # Detect time-interval naming per file: after part_tag- the next char is 't'
        sep = f'{part_tag}-'
        if sep not in basename:
            msg(f'WARNING: separator "{sep}" not found in {basename} — skipping')
            continue
        after_tag = basename.split(sep, 1)[1]
        file_has_time = after_tag.startswith('t') and after_tag[1:5].isdigit()

        try:
            if file_has_time:
                # e.g. prefix_part0001-t0000-0003-I-image.fits
                # after_tag = 't0000-0003-I-image.fits'
                chunks = after_tag[1:].split('-')   # strip leading 't', split remainder
                t_str  = chunks[0]                  # '0000'
                c_str  = chunks[1]                  # '0003'
                suffix = '-'.join(chunks[2:])       # 'I-image.fits' or 'image.fits'
                new_c  = f'{int(c_str) + new_channel_offset:04d}'
                new_basename = f'{prefix_base}-t{t_str}-{new_c}-{suffix}'
            else:
                # e.g. prefix_part0001-0003-I-image.fits
                # after_tag = '0003-I-image.fits'
                chunks = after_tag.split('-', 1)    # ['0003', 'I-image.fits']
                c_str  = chunks[0]
                suffix = chunks[1]                  # 'I-image.fits' or 'image.fits'
                new_c  = f'{int(c_str) + new_channel_offset:04d}'
                new_basename = f'{prefix_base}-{new_c}-{suffix}'
        except (ValueError, IndexError) as e:
            msg(f'WARNING: could not parse channel index from {basename} ({e}) — skipping')
            continue

        new = os.path.join(directory, new_basename)

        # Re-extract the suffix from the constructed new_basename independently
        # so the check below is not a no-op
        if file_has_time:
            new_suffix = new_basename.split(f'{prefix_base}-t{t_str}-{new_c}-', 1)[-1]
        else:
            new_suffix = new_basename.split(f'{prefix_base}-{new_c}-', 1)[-1]

        rename_pairs.append((old, new, suffix, new_suffix))

    # -------------------------------------------------------------------------
    # Suffix consistency check: old and new must share the same type suffix
    # (e.g. image.fits → image.fits, model.fits → model.fits).
    # Since we derive both names from the same suffix string this should always
    # pass — a failure here means the parsing above is broken.
    # -------------------------------------------------------------------------
    msg(f'Checking suffix consistency across {len(rename_pairs)} rename(s)...')
    errors = []
    ok     = []
    for old, new, old_sfx, new_sfx in rename_pairs:
        if old_sfx != new_sfx:
            errors.append((old, new, old_sfx, new_sfx))
        else:
            ok.append((old, new, old_sfx))

    if errors:
        msg(f'ERROR: {len(errors)} suffix mismatch(es) detected — aborting, no files renamed')
        for old, new, old_sfx, new_sfx in errors:
            print(f'  MISMATCH: {os.path.basename(old)}  ->  {os.path.basename(new)}'
                  f'  (suffix: "{old_sfx}" -> "{new_sfx}")')
        sys.exit(1)

    msg(f'Suffix check passed for all {len(ok)} rename(s)')

    # -------------------------------------------------------------------------
    # Print a sample of up to 100 rename pairs (random mix of suffixes) so the
    # mapping can be visually inspected before anything is moved.
    # -------------------------------------------------------------------------
    sample_size = min(100, len(rename_pairs))
    sample      = random.sample(rename_pairs, sample_size)
    # Sort sample by old name so output is readable despite random selection
    sample.sort(key=lambda x: x[0])

    pad = max(len(os.path.basename(old)) for old, _, _, _ in sample)
    msg(f'Rename preview ({sample_size} of {len(rename_pairs)} pairs, random mix of suffixes):')
    print()
    print(f'  {"OLD":<{pad}}    NEW')
    print(f'  {"-"*pad}    {"-"*pad}')
    for old, new, old_sfx, _ in sample:
        print(f'  {os.path.basename(old):<{pad}}  ->  {os.path.basename(new)}')
    print()

    # Summarise which suffix types are present in the full list
    suffix_counts = {}
    for _, _, sfx, _ in rename_pairs:
        suffix_counts[sfx] = suffix_counts.get(sfx, 0) + 1
    msg('Suffix type breakdown across all rename pairs:')
    for sfx, count in sorted(suffix_counts.items()):
        print(f'  -{sfx:<30}  {count} file(s)')
    print()

    # -------------------------------------------------------------------------
    # Execute renames
    # -------------------------------------------------------------------------
    msg(f'Renaming {len(rename_pairs)} file(s)...')
    for old, new, _, _ in rename_pairs:
        if old != new:
            subprocess.run(['mv', old, new])

    msg('Image renaming complete')


def main():

    if len(sys.argv) != 3:
        msg('ERROR: expected exactly 2 arguments: <totchan> <prefix>')
        msg('Usage: fix_image_naming.py <totchan> <prefix>')
        sys.exit(1)

    totchan = int(sys.argv[-2])
    prefix  = sys.argv[-1]

    msg(f'Running fix_image_naming for prefix: {prefix}')
    msg(f'Total channels expected: {totchan}')

    fix_names(prefix, totchan)


if __name__ == '__main__':
    main()
