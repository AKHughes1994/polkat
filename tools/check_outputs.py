# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

# Check that a pipeline stage left the outputs it should, for every field of the
# project. Run from the working directory:
#
#     python3 tools/check_outputs.py 2GC
#
# The fields are those 2GC.py processes: the targets, the primary calibrator and
# the polarization angle calibrator per scan, and the secondary calibrators per
# scan (not when CAL_SKIP_CALS is set), restricted to working_names when
# PRE_FIELDS is set. They come from project_info.json, never from what is on
# disk, so a field whose stage failed silently is still expected.
#
# 2GC: every field has IMAGES/<field>/ with a final MFS image (*pcalmask-MFS-*image.fits)
# and, unless RMSYNTH_RESIDUAL_MS_DIR is blank, every target also has its residual MS there.
# Each missing output is printed, and the exit status is 1 if any is missing.

import glob
import json
import os
import os.path as o
import sys

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg


def fields(project_info):
    """(name, ms, is_target) for every field 2GC processes."""

    filtered = cfg.PRE_FIELDS != ''
    allowed = set(project_info.get('working_names', []))
    out = []

    for name, ms in zip(project_info['target_names'], project_info['target_ms']):
        if not filtered or name in allowed:
            out.append((name, ms, True))

    def per_scan(name, ms_list):
        for ms in ms_list:
            out.append((f"{name}_scan{ms.split('_scan')[-1].replace('.ms', '')}", ms, False))

    if not filtered or project_info['primary_name'] in allowed:
        per_scan(project_info['primary_name'], project_info['primary_ms'])

    pacal_name = project_info.get('polang_name', '')
    if pacal_name != '' and project_info.get('polang_ms'):
        if not filtered or pacal_name in allowed:
            per_scan(pacal_name, project_info['polang_ms'])

    if not cfg.CAL_SKIP_CALS:
        for name, ms_list in zip(project_info['secondary_names'], project_info['secondary_ms']):
            if not filtered or name in allowed:
                per_scan(name, ms_list)

    return out


def residual_ms(ms):
    """Path of the residual MS 2GC makes from a target's MS."""

    return o.join(cfg.RMSYNTH_RESIDUAL_MS_DIR, o.basename(ms))


def check_2gc(project_info):
    problems = []
    checked = fields(project_info)

    for name, ms, is_target in checked:
        img_dir = o.join(cfg.IMAGES, name.replace(' ', '_'))
        if not o.isdir(img_dir) or not os.listdir(img_dir):
            problems.append(f'{name}: {img_dir} is missing or empty')
        elif not glob.glob(o.join(img_dir, '*pcalmask-MFS-*image.fits')):
            problems.append(f'{name}: no final MFS image (*pcalmask-MFS-*image.fits) in {img_dir}')

        if is_target and cfg.RMSYNTH_RESIDUAL_MS_DIR != '' and not o.isdir(residual_ms(ms)):
            problems.append(f'{name}: no residual MS {residual_ms(ms)}')

    print(f'2GC: {len(checked)} field(s) checked, {len(problems)} problem(s)')
    for problem in problems:
        print(f'  MISSING  {problem}')
    return problems


def main():
    if len(sys.argv) != 2:
        sys.exit('Usage: python3 tools/check_outputs.py 2GC')

    with open('project_info.json') as f:
        project_info = json.load(f)

    if sys.argv[1] == '2GC':
        problems = check_2gc(project_info)
    else:
        sys.exit(f'No checks for stage {sys.argv[1]}')

    sys.exit(1 if problems else 0)


if __name__ == '__main__':
    main()
