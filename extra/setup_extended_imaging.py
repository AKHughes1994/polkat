#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.a.uk
#
# Joint, full-UV-range imaging of an arbitrary set of measurement sets, for
# recovering extended/diffuse structure that the standard short-baseline
# uv-cut (WSC_BASELINE_CUT / WSC_TAPERMASK) is designed to suppress in 2GC.
# Unlike 2GC/setup_2GC_twostage.py this does no calibration at all: it is
# just blind image -> mask -> final image (or final image only, if a mask
# is supplied), run once across whichever MS are given, jointly.

import argparse
import glob
import os.path as o
import sys
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg


EPILOG = '''\
The MS pointer is the only required argument, and must end in:
  .ms   a glob pattern for the measurement sets to combine, e.g. "field*.ms"
        or "*obs123*.ms" (quote it so the shell doesn't expand it first)
  .txt  a text file listing MS paths, one per line (blank lines and lines
        starting with '#' are ignored)

Examples:
  python3 extra/setup_extended_imaging.py "field*.ms"
  python3 extra/setup_extended_imaging.py my_ms_list.txt -i chpc -d data
  python3 extra/setup_extended_imaging.py "field*.ms" -m img_field_mask.fits

Regardless of config.py, this script always forces mf-weighting on and the
inner Tukey taper / minuv-l / maxuv-l cuts off, so every image uses the full
baseline range.
'''


def build_parser():
    parser = argparse.ArgumentParser(
        prog='setup_extended_imaging.py',
        description='Jointly image a set of MS with the full UV range, for extended/diffuse structure.',
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('ms_pointer',
        help="Either a glob pattern ending in '.ms', or a '.txt' file listing MS paths, one per line.")
    parser.add_argument('-i', '--infra', dest='infra', default='idia', type=str.lower,
        choices=['idia', 'chpc', 'hippo', 'node', 'glam'],
        help="Compute infrastructure (idia / chpc / hippo / node / glam). Default: idia.")
    parser.add_argument('-m', '--mask', dest='mask', default=None,
        help='Path to an existing FITS mask. If given, the blind image and mask-generation steps are skipped.')
    parser.add_argument('-d', '--datacol', dest='datacol', default='corrected_data',
        help="Column to image, case-insensitive (e.g. 'corrected', 'CORRECTED_DATA', 'data'). Default: corrected_data.")
    return parser


def normalize_datacol(raw):
    col = raw.strip().upper()
    if col == 'DATA':
        return 'DATA'
    if not col.endswith('_DATA'):
        col += '_DATA'
    return col


def make_label(ms_pointer):
    base = o.basename(ms_pointer)
    for ext in ('.txt', '.ms'):
        if base.lower().endswith(ext):
            base = base[:-len(ext)]
    base = base.replace('*', '').replace('?', '').strip('_-.')
    return base if base else 'extended'


def resolve_ms_list(ms_pointer, parser):
    lowered = ms_pointer.lower()
    if lowered.endswith('.txt'):
        mode = 'list file'
        if not o.isfile(ms_pointer):
            parser.error(f"'{ms_pointer}' ends in '.txt' but is not a file that exists.")
        with open(ms_pointer) as f:
            lines = [line.strip() for line in f if line.strip() and not line.strip().startswith('#')]
        bad_lines = [line for line in lines if not line.lower().endswith('.ms')]
        if bad_lines:
            parser.error(
                f"'{ms_pointer}' must list only measurement sets (paths ending in '.ms'); "
                f"found {len(bad_lines)} invalid line(s), e.g. '{bad_lines[0]}'.")
        ms_list = lines
    elif lowered.endswith('.ms'):
        mode = 'glob pattern'
        ms_list = sorted(glob.glob(ms_pointer))
    else:
        parser.error(f"'{ms_pointer}' must end in '.ms' (a glob pattern) or '.txt' (a list file) -- see --help.")

    if not ms_list:
        parser.error(f"No measurement sets resolved from '{ms_pointer}'.")

    missing = [m for m in ms_list if not o.isdir(m)]
    ms_list = [m for m in ms_list if o.isdir(m)]
    if missing:
        print(gen.col('WARNING') + f'{len(missing)} listed MS not found on disk, skipping:')
        for m in missing:
            print(gen.col() + m)
    if not ms_list:
        parser.error('None of the resolved measurement sets exist on disk.')

    return mode, ms_list


def main():

    if len(sys.argv) == 1 or [a.lower() for a in sys.argv[1:]].count('help'):
        build_parser().print_help()
        sys.exit(0)

    parser = build_parser()
    args = parser.parse_args()

    mode, ms_list = resolve_ms_list(args.ms_pointer, parser)
    datacol = normalize_datacol(args.datacol)
    label = make_label(args.ms_pointer)
    code = gen.get_target_code(label)

    gen.preamble()
    print(gen.col() + 'Extended imaging (joint, full UV range) setup')
    gen.print_spacer()

    IMAGES = cfg.IMAGES
    SCRIPTS = cfg.SCRIPTS
    LOGS = cfg.LOGS
    TOOLS = cfg.TOOLS

    gen.setup_dir(IMAGES)
    gen.setup_dir(SCRIPTS)
    gen.setup_dir(LOGS)

    INFRASTRUCTURE, CONTAINER_PATH = gen.set_infrastructure(['', args.infra])
    USE_SINGULARITY = cfg.USE_SINGULARITY
    CONTAINER_RUNNER = 'singularity exec ' if CONTAINER_PATH is not None else ''

    PYTHON3_CONTAINER = gen.get_container(CONTAINER_PATH, cfg.PYTHON3_PATTERN, USE_SINGULARITY)
    WSCLEAN_CONTAINER = gen.get_container(CONTAINER_PATH, cfg.WSCLEAN_PATTERN, USE_SINGULARITY)

    gen.print_spacer()
    print(gen.col('MS pointer') + args.ms_pointer)
    print(gen.col('Input mode') + mode)
    print(gen.col('Matched MS') + str(len(ms_list)))
    for myms in ms_list:
        print(gen.col() + myms)
    print(gen.col('Data column') + datacol)
    print(gen.col('MF weighting') + 'Forced on')
    print(gen.col('Tukey taper / UV cuts') + 'Forced off (full baseline range)')
    if args.mask:
        print(gen.col('Mask') + args.mask + ' (skipping blind image + mask generation)')
    else:
        print(gen.col('Mask') + 'None (will be generated from a blind image)')
    print(gen.col('Output label') + label)
    gen.print_spacer()

    blind_img_prefix = IMAGES + f'/img_{label}_extblind'
    final_img_prefix = IMAGES + f'/img_{label}_extimage'
    kill_file = SCRIPTS + '/kill_extended_imaging_jobs_' + label + '.sh'

    steps = []
    n = 0
    mask = args.mask

    if not mask:

        step = {}
        step['step'] = n
        step['comment'] = f'Shallow blind wsclean (Stokes I) jointly across {len(ms_list)} MS for extended-structure masking'
        step['dependency'] = None
        step['id'] = 'WSEBL' + code
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step, INFRASTRUCTURE, cfg.WSC_ABSMEM)
        prefix = CONTAINER_RUNNER + WSCLEAN_CONTAINER + ' ' if USE_SINGULARITY else ''
        syscall = ''
        imcall = gen.generate_syscall_wsclean(mslist=ms_list,
                imgname=blind_img_prefix,
                datacol=datacol,
                chanout=cfg.WSC_BLIND_CHANNELSOUT,
                nomodel=True,
                pol='I',
                intervalsout=False,
                mfweight=True,
                localrms=cfg.WSC_LOCALRMS_BLIND,
                automask=cfg.WSC_AUTOMASK_BLIND,
                autothreshold=cfg.WSC_AUTOTHRESHOLD_BLIND,
                threshold=cfg.WSC_THRESHOLD_BLIND,
                tukeytaper=False,
                minuvl='',
                maxuvl='',
                sourcelist=False,
                absmem=absmem)
        for call in imcall:
            syscall += prefix + call + '\n\n'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Make extended-structure cleaning mask'
        step['dependency'] = n - 1
        step['id'] = 'MASKX' + code
        syscall = CONTAINER_RUNNER + WSCLEAN_CONTAINER + ' ' if USE_SINGULARITY else ''
        breizorro_call, mask = gen.generate_syscall_breizorro(
                restoredimage=f'{blind_img_prefix}-MFS-image.fits',
                outfile=f'{blind_img_prefix}-MFS-image.mask.fits')
        syscall += breizorro_call
        step['syscall'] = syscall
        steps.append(step)
        n += 1

    step = {}
    step['step'] = n
    step['comment'] = f'Run wsclean, masked deconvolution jointly across {len(ms_list)} MS with full UV range ({datacol})'
    step['dependency'] = n - 1 if n > 0 else None
    step['id'] = 'WSEFN' + code
    step['slurm_config'] = cfg.SLURM_WSCLEAN
    step['pbs_config'] = cfg.PBS_WSCLEAN
    absmem = gen.absmem_helper(step, INFRASTRUCTURE, cfg.WSC_ABSMEM)
    prefix = CONTAINER_RUNNER + WSCLEAN_CONTAINER + ' ' if USE_SINGULARITY else ''
    syscall = ''
    imcall = gen.generate_syscall_wsclean(mslist=ms_list,
            imgname=final_img_prefix,
            datacol=datacol,
            mask=mask,
            chanout=cfg.WSC_PCAL_CHANNELSOUT,
            nomodel=True,
            sourcelist=False,
            mfweight=True,
            tukeytaper=False,
            minuvl='',
            maxuvl='',
            absmem=absmem)
    for call in imcall:
        syscall += prefix + call + '\n\n'
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    if cfg.WSC_MAX_CHANNELS < cfg.WSC_PCAL_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
        step = {}
        step['step'] = n
        step['comment'] = 'Homogenize the extended-image resolution across frequency channels'
        step['dependency'] = n - 1
        step['id'] = 'HOEFN' + code
        prefix = CONTAINER_RUNNER + PYTHON3_CONTAINER + ' ' if USE_SINGULARITY else ''
        syscall = prefix + f'python3 {TOOLS}/fix_image_naming.py {cfg.WSC_PCAL_CHANNELSOUT} {final_img_prefix}\n\n'
        syscall += prefix + f'python3 {TOOLS}/homogenize_beams.py {final_img_prefix}'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

    if not cfg.SKIP_PB:
        if cfg.BAND == 'not yet determined':
            print(gen.col('WARNING') + 'project_info.json not found / band unknown, skipping primary beam correction')
        else:
            step = {}
            step['step'] = n
            step['comment'] = 'Apply primary beam correction to the extended image'
            step['dependency'] = n - 1
            step['id'] = 'PBEFN' + code
            syscall = ''
            images = [f'{final_img_prefix}-MFS-image.fits']
            if cfg.WSC_MAX_CHANNELS < cfg.WSC_PCAL_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                images.append(f'{final_img_prefix}-MFS-image.homogenized.fits')
            for image in images:
                syscall += CONTAINER_RUNNER + PYTHON3_CONTAINER + ' ' if USE_SINGULARITY else ''
                syscall += f'python3 {TOOLS}/pbcor_katbeam.py --band {cfg.BAND[0]} {image}\n\n'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------

    submit_file = 'submit_extended_imaging_jobs.sh'

    f = open(submit_file, 'w')
    f.write('#!/usr/bin/env bash\n')
    f.write('export SINGULARITY_BINDPATH=' + cfg.BINDPATH + '\n')
    f.write('\n#---------------------------------------\n')
    f.write('# ' + label + '\n')
    f.write('#---------------------------------------\n')

    id_list = []
    for step in steps:

        nd = step['id']
        id_list.append(nd)
        if step['dependency'] is not None:
            dependency = steps[step['dependency']]['id']
        else:
            dependency = None
        syscall = step['syscall']
        slurm_config = step['slurm_config'] if 'slurm_config' in step.keys() else cfg.SLURM_DEFAULTS
        pbs_config = step['pbs_config'] if 'pbs_config' in step.keys() else cfg.PBS_DEFAULTS
        comment = step['comment']

        run_command = gen.job_handler(syscall=syscall,
                        jobname=nd,
                        infrastructure=INFRASTRUCTURE,
                        dependency=dependency,
                        slurm_config=slurm_config,
                        pbs_config=pbs_config)

        f.write('\n# ' + comment + '\n')
        f.write(run_command)

    if INFRASTRUCTURE != 'node':
        f.write('\n# Generate kill script\n')
    if INFRASTRUCTURE == 'idia' or INFRASTRUCTURE == 'hippo':
        f.write('echo "scancel "$' + '" "$'.join(id_list) + ' > ' + kill_file + '\n')
    elif INFRASTRUCTURE == 'chpc':
        f.write('echo "qdel "$' + '" "$'.join(id_list) + ' > ' + kill_file + '\n')

    f.close()
    gen.make_executable(submit_file)

    gen.print_spacer()
    print(gen.col('Run file') + submit_file)
    gen.print_spacer()

    # ------------------------------------------------------------------------------


if __name__ == "__main__":

    main()
