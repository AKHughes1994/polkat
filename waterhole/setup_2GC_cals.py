#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.a.uk

import json
import math
import os
import os.path as o
import glob
import re
import sys
import yaml
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg


SET_REFANT        = True
ADAPTIVE_CHANNELS = True


def get_adaptive_freq_intervals(yaml_path, n_model_channels):
    """Parse a QuartiCal YAML and return (ratio, overrides, zeros) where:
      overrides — list of (term, original_fi, new_fi) for terms that need coarsening
      zeros     — list of term names with freq_interval=0 (whole-band, left alone)
    ratio = ceil(PRE_NCHANS / n_model_channels)."""
    ratio = math.ceil(cfg.PRE_NCHANS / n_model_channels)
    overrides = []
    zeros = []
    with open(yaml_path) as f:
        qc = yaml.safe_load(f)
    for term in qc.get('solver', {}).get('terms', []):
        fi = int(str(qc.get(term, {}).get('freq_interval', 0)).strip())
        if fi == 0:
            zeros.append(term)
        elif fi < ratio:
            overrides.append((term, fi, ratio))
    return ratio, overrides, zeros


def main():

    USE_SINGULARITY = cfg.USE_SINGULARITY

    gen.preamble()
    print(gen.col()+'2GC Calibrator Imaging & Self-Calibration setup')
    gen.print_spacer()

    # ------------------------------------------------------------------------------
    # Setup paths, required containers, infrastructure
    # ------------------------------------------------------------------------------

    OXKAT      = cfg.OXKAT
    DATA       = cfg.DATA
    IMAGES     = cfg.IMAGES
    SCRIPTS    = cfg.SCRIPTS
    TOOLS      = cfg.TOOLS
    GAINTABLES = cfg.GAINTABLES
    LOGS       = cfg.LOGS

    gen.setup_dir(GAINTABLES)
    gen.setup_dir(IMAGES)
    gen.setup_dir(cfg.LOGS)
    gen.setup_dir(cfg.SCRIPTS)

    INFRASTRUCTURE, CONTAINER_PATH = gen.set_infrastructure(sys.argv)
    CONTAINER_RUNNER = 'singularity exec ' if CONTAINER_PATH is not None else ''

    PYTHON3_CONTAINER   = gen.get_container(CONTAINER_PATH, cfg.PYTHON3_PATTERN,   USE_SINGULARITY)
    CASA_CONTAINER      = gen.get_container(CONTAINER_PATH, cfg.CASA_PATTERN,      USE_SINGULARITY)
    WSCLEAN_CONTAINER   = gen.get_container(CONTAINER_PATH, cfg.WSCLEAN_PATTERN,   USE_SINGULARITY)
    QUARTICAL_CONTAINER = gen.get_container(CONTAINER_PATH, cfg.QUARTICAL_PATTERN, USE_SINGULARITY)

    with open('project_info.json') as f:
        project_info = json.load(f)

    band  = project_info['band']
    myms  = project_info['working_ms']

    ref_ant_arg = None
    if SET_REFANT:
        raw_ref_ant = project_info.get('ref_ant', None)
        if raw_ref_ant is None:
            print(gen.col('Reference Antenna')+'WARNING: ref_ant not found in project_info.json — defaulting to antenna index 0')
            ref_ant_arg = 0
        elif not isinstance(raw_ref_ant, str):
            print(gen.col('Reference Antenna')+f'WARNING: ref_ant value "{raw_ref_ant}" is not a string — defaulting to antenna index 0')
            ref_ant_arg = 0
        else:
            try:
                ref_ant_arg = int(raw_ref_ant.split(',')[0].strip())
                print(gen.col('Reference Antenna')+str(ref_ant_arg))
            except (ValueError, IndexError):
                print(gen.col('Reference Antenna')+f'WARNING: ref_ant value "{raw_ref_ant}" could not be parsed as a comma-separated list of ints — defaulting to antenna index 0')
                ref_ant_arg = 0

    # UV taper settings
    if not cfg.WSC_TAPERMASK:
        tukeytaper = False
        minuvl     = ''
    else:
        tukeytaper = cfg.WSC_TUKEYTAPER
        minuvl     = cfg.WSC_MINUVL

    # ------------------------------------------------------------------------------
    # Validate config
    # ------------------------------------------------------------------------------

    cals_to_image = cfg.CALS_TO_IMAGE.strip()
    if not cals_to_image:
        print(gen.col('ERROR')+'CALS_TO_IMAGE is empty — nothing to do.')
        sys.exit()

    field_list   = [f.strip() for f in cals_to_image.split(',') if f.strip()]
    combine_scan = cfg.CAL_COMBINESCAN

    # Build the calibrator MS names from scantimes log written during INFO setup.
    scantimes_patterns = [
        f'scantimes*{myms}*log',
        f'scantimes*{o.basename(myms)}*log',
    ]
    scantimes_matches = []
    for pattern in scantimes_patterns:
        scantimes_matches.extend(glob.glob(pattern))
    scantimes_matches = sorted(set(scantimes_matches))

    if not scantimes_matches:
        print(gen.col('ERROR')+'Missing scantimes log in CWD.')
        print(gen.col('CWD')+os.getcwd())
        print(gen.col('Tried patterns')+', '.join(scantimes_patterns))
        print(gen.col('Hint')+'Re-run INFO/1GC setup to regenerate scantimes logs.')
        sys.exit()

    # If there are multiple matches, pick the latest by modification time.
    scantimes_log = max(scantimes_matches, key=o.getmtime)
    if len(scantimes_matches) > 1:
        print(gen.col('WARNING')+f'Multiple scantimes logs matched ({len(scantimes_matches)}); using latest: {scantimes_log}')

    field_scan_map = {}
    with open(scantimes_log, 'r') as fh:
        for line in fh:
            if '|' in line:
                msg = line.split('|', 1)[1].strip()
            else:
                msg = line.strip()
            if not msg:
                continue
            if re.match(r'^\d+\s+', msg):
                parts = msg.split()
                if len(parts) >= 2:
                    scan = int(parts[0])
                    field = parts[1]
                    field_scan_map.setdefault(field, set()).add(scan)

    if not field_scan_map:
        print(gen.col('ERROR')+f'No scan rows parsed from {scantimes_log}. Re-run INFO/scan_times.py.')
        sys.exit()

    print(gen.col('Scantimes log')+scantimes_log)
    print(gen.col('Cal scan summary')+f'combine_scan={combine_scan}')
    for field in field_list:
        scans = sorted(field_scan_map.get(field, []))
        if scans:
            scan_list = ','.join(str(s) for s in scans)
            print(gen.col('  '+field)+f'{len(scans)} scan(s): {scan_list}')
        else:
            print(gen.col('  '+field)+'0 scan(s): not present in scantimes log')

    # Build the list of split MS names expected from casa_split_cals.py
    ms_list = []
    for field in field_list:
        scans = sorted(field_scan_map.get(field, []))
        if not scans:
            print(gen.col('WARNING')+f'Field {field} not found in {scantimes_log}')
            continue

        if combine_scan:
            ms_list.append((field, myms.replace('.ms', f'_{field}.ms')))
        else:
            for scan in scans:
                ms_list.append((field, myms.replace('.ms', f'_{field}_scan{scan:03d}.ms')))

    if not ms_list:
        print(gen.col('ERROR')+'No calibrator MS targets built from scantimes and CALS_TO_IMAGE.')
        sys.exit()

    # ------------------------------------------------------------------------------
    # Build steps for each calibrator MS
    # ------------------------------------------------------------------------------

    freq_int_overrides = []
    if ADAPTIVE_CHANNELS:
        ratio, freq_int_overrides, zeros = get_adaptive_freq_intervals(cfg.CAL_2GC_FILE, cfg.WSC_CAL_CHANNELSOUT)
        yaml_name = o.basename(cfg.CAL_2GC_FILE)
        parts = []
        if zeros:
            parts.append(', '.join(f'{t}: 0 (all channels, skipping)' for t in zeros))
        if freq_int_overrides:
            parts.append(', '.join(f'{t}: {old}→{new}' for t, old, new in freq_int_overrides))
        summary = ' | '.join(parts) if parts else 'no overrides needed'
        print(gen.col('Adaptive Channels')+f'{yaml_name} | ratio {ratio} — {summary}')

    stamp      = gen.timenow()
    all_steps  = []   # list of (steps, kill_file, label)
    codes_seen = []
    ii         = 1

    for (fieldname, calms) in ms_list:

        filename_fieldname = gen.scrub_target_name(o.basename(calms).replace('.ms', ''))
        code = gen.get_target_code(filename_fieldname)
        if code in codes_seen:
            code += '_' + str(ii)
            ii += 1
        codes_seen.append(code)

        gain_outdir = GAINTABLES + '/2GC_' + filename_fieldname + f'_{stamp}.qc/'
        log_outdir  = LOGS       + '/2GC_' + filename_fieldname + f'_{stamp}.qc/'

        mask_img_prefix  = IMAGES + f'/img_{o.basename(calms)}_mask'
        diag_img_prefix  = IMAGES + f'/img_{o.basename(calms)}_diagnostic'
        kill_file        = SCRIPTS + f'/kill_2GC_cal_{filename_fieldname}.sh'

        gen.print_spacer()
        print(gen.col('Calibrator') + fieldname)
        print(gen.col('MS')         + calms)
        print(gen.col('Code')       + code)

        steps = []
        n     = 0

        # Split is scheduled as the first submit step; this tree starts after it.

        # ---- Blind image + mask --------------------------------------------------
        step = {}
        step['step']       = n
        step['comment']    = f'Shallow blind wsclean on DATA for {fieldname}'
        step['dependency'] = None
        step['id']         = 'WSMSK' + code
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config']   = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step, INFRASTRUCTURE, cfg.WSC_ABSMEM)
        prefix = CONTAINER_RUNNER + WSCLEAN_CONTAINER + ' ' if USE_SINGULARITY else ''
        syscall = ''
        imcall = gen.generate_syscall_wsclean(
            mslist        = [calms],
            imgname       = mask_img_prefix,
            datacol       = 'DATA',
            chanout       = cfg.WSC_BLIND_CHANNELSOUT,
            imsize        = cfg.WSC_CAL_IMSIZE,
            nomodel       = True,
            pol           = 'I',
            mask          = False,
            intervalsout  = False,
            mfweight      = True,
            localrms      = False,
            automask      = cfg.WSC_CALS_BLINDMASK,
            autothreshold = cfg.WSC_AUTOTHRESHOLD_BLIND,
            tukeytaper    = tukeytaper,
            minuvl        = minuvl,
            absmem        = absmem)
        for call in imcall:
            syscall += prefix + call + '\n\n'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        step = {}
        step['step']       = n
        step['comment']    = f'Make cleaning mask for {fieldname}'
        step['dependency'] = n - 1
        step['id']         = 'MKMSK' + code
        syscall  = CONTAINER_RUNNER + WSCLEAN_CONTAINER + ' ' if USE_SINGULARITY else ''
        syscall += gen.generate_syscall_breizorro(
            restoredimage = f'{mask_img_prefix}-MFS-image.fits',
            outfile       = f'{mask_img_prefix}-MFS-image.mask.fits')[0]
        step['syscall'] = syscall
        steps.append(step)
        n += 1
        mask = f'{mask_img_prefix}-MFS-image.mask.fits'

        # ---- WSDMA: data-masked image using WSC_SHALLOWMASK ---------------------
        step = {}
        step['step']       = n
        step['comment']    = f'Run wsclean masked deconvolution of DATA for {fieldname}'
        step['dependency'] = n - 1
        step['id']         = 'WSDMA' + code
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config']   = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step, INFRASTRUCTURE, cfg.WSC_ABSMEM)
        prefix  = CONTAINER_RUNNER + WSCLEAN_CONTAINER + ' ' if USE_SINGULARITY else ''
        syscall = ''
        imcall  = gen.generate_syscall_wsclean(
            mslist            = [calms],
            imgname           = diag_img_prefix,
            datacol           = 'DATA',
            mask              = mask,
            weight            = cfg.WSC_WEIGHT_CAL,
            chanout           = cfg.WSC_CAL_CHANNELSOUT,
            imsize            = cfg.WSC_CAL_IMSIZE,
            tukeytaper        = tukeytaper,
            minuvl            = minuvl,
            automask          = cfg.WSC_SHALLOWMASK,
            localrms          = cfg.WSC_INTER_LOCALRMS,
            autothreshold     = cfg.WSC_INTER_AUTOTHRESHOLD,
            joinpolarizations = False,
            multiscale        = False,
            splitpol          = True,
            mfweight          = False,
            intervalsout      = False,
            threshold         = False,
            nomodel           = True,
            sourcelist        = False,
            absmem            = absmem)
        for call in imcall:
            syscall += prefix + call + '\n\n'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        # ---- fix_image_naming if needed ------------------------------------------
        if cfg.WSC_MAX_CHANNELS < cfg.WSC_CAL_CHANNELSOUT:
            step = {}
            step['step']       = n
            step['comment']    = f'Fix image naming for {fieldname}'
            step['dependency'] = n - 1
            step['id']         = 'FXDMA' + code
            prefix  = CONTAINER_RUNNER + PYTHON3_CONTAINER + ' ' if USE_SINGULARITY else ''
            syscall = prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_CAL_CHANNELSOUT} {diag_img_prefix}'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        # ---- Predict + optional mod_model_selfcal (V-only, mask as spatial) ------
        step = {}
        step['step']       = n
        step['comment']    = f'Run wsclean-predict on data model for {fieldname}'
        step['dependency'] = n - 1
        step['id']         = 'PRDMA' + code
        step['slurm_config'] = cfg.SLURM_PREDICT
        step['pbs_config']   = cfg.PBS_WSCLEAN
        absmem  = gen.absmem_helper(step, INFRASTRUCTURE, cfg.WSC_ABSMEM)
        prefix  = CONTAINER_RUNNER + WSCLEAN_CONTAINER + ' ' if USE_SINGULARITY else ''
        syscall = prefix + 'python3 ' + TOOLS + '/fix_nan_models.py ' + diag_img_prefix + '\n\n'
        if cfg.MOD_MODEL_SELFCAL and cfg.WSC_POL != 'I':
            syscall += prefix + (f'python3 {TOOLS}/mod_model_selfcal.py '
                                 f'--identifier {diag_img_prefix} --stokes V '
                                 f'--spatial {mask}\n\n')
        syscall += prefix + gen.generate_syscall_predict(
            msname = calms,
            imgname = diag_img_prefix,
            chanout = cfg.WSC_CAL_CHANNELSOUT,
            absmem  = absmem)
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        # ---- Single round selfcal with CAL_2GC_YAML_COMPLEX ----------------------
        step = {}
        step['step']       = n
        step['comment']    = f'Run Quartical self-calibration on {fieldname}'
        step['dependency'] = n - 1
        step['id']         = 'CL2GC' + code
        syscall     = CONTAINER_RUNNER + QUARTICAL_CONTAINER + ' ' if USE_SINGULARITY else ''
        extra_args  = f'output.gain_directory={gain_outdir} output.log_directory={log_outdir}'
        if ref_ant_arg is not None:
            extra_args += f' solver.reference_antenna={ref_ant_arg}'
        for term, _, new_fi in freq_int_overrides:
            extra_args += f' {term}.freq_interval={new_fi}'
        if not cfg.CAL_1GC_APPLYPARANG:
            extra_args += ' output.apply_p_jones_inv=true'
        if minuvl != '':
            extra_args += f' input_ms.select_uv_range=[{minuvl},0]'
        syscall += gen.generate_syscall_quartical(
            yaml       = cfg.CAL_2GC_FILE,
            myms       = calms,
            extra_args = extra_args)
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        # ---- Parang safety check (only when parang not applied in 1GC) -----------
        if not cfg.CAL_1GC_APPLYPARANG:
            fallback_extra = f'output.gain_directory={gain_outdir} output.log_directory={log_outdir}'
            if minuvl != '':
                fallback_extra += f' input_ms.select_uv_range=[{minuvl},0]'
            fallback_qc_cmd = (CONTAINER_RUNNER + QUARTICAL_CONTAINER + ' ' if USE_SINGULARITY else '') + \
                gen.generate_syscall_quartical(yaml=cfg.CAL_2GC_YAML_COMPLEX, myms=calms, extra_args=fallback_extra)

            step = {}
            step['step']       = n
            step['comment']    = f'Check parang selfcal succeeded for {fieldname}; fall back if not'
            step['dependency'] = n - 1
            step['id']         = 'CHKPJ' + code
            step['slurm_config'] = cfg.SLURM_DEFAULTS
            step['pbs_config']   = cfg.PBS_DEFAULTS
            prefix  = CONTAINER_RUNNER + PYTHON3_CONTAINER + ' ' if USE_SINGULARITY else ''
            syscall = prefix + (
                f'python3 {TOOLS}/check_and_fix_parang_selfcal.py '
                f'{log_outdir} {calms} "{fallback_qc_cmd}"')
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        # ---- Final image with CORRECTED_DATA ------------------------------------
        step = {}
        step['step']       = n
        step['comment']    = f'Final image of {fieldname} with CORRECTED_DATA'
        step['dependency'] = n - 1
        step['id']         = 'WSCMA' + code
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config']   = cfg.PBS_WSCLEAN
        absmem  = gen.absmem_helper(step, INFRASTRUCTURE, cfg.WSC_ABSMEM)
        prefix  = CONTAINER_RUNNER + WSCLEAN_CONTAINER + ' ' if USE_SINGULARITY else ''
        syscall = ''
        imcall  = gen.generate_syscall_wsclean(
            mslist            = [calms],
            imgname           = diag_img_prefix,
            datacol           = 'CORRECTED_DATA',
            mask              = mask,
            weight            = cfg.WSC_WEIGHT_CAL,
            chanout           = cfg.WSC_CAL_CHANNELSOUT,
            imsize            = cfg.WSC_CAL_IMSIZE,
            tukeytaper        = tukeytaper,
            minuvl            = minuvl,
            automask          = cfg.WSC_AUTOMASK,
            localrms          = cfg.WSC_INTER_LOCALRMS,
            autothreshold     = cfg.WSC_INTER_AUTOTHRESHOLD,
            joinpolarizations = False,
            multiscale        = False,
            splitpol          = True,
            mfweight          = False,
            intervalsout      = False,
            threshold         = False,
            nomodel           = True,
            sourcelist        = False,
            absmem            = absmem)
        for call in imcall:
            syscall += prefix + call + '\n\n'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        # ---- fix_image_naming + homogenize for final image ----------------------
        if cfg.WSC_MAX_CHANNELS < cfg.WSC_CAL_CHANNELSOUT:
            step = {}
            step['step']       = n
            step['comment']    = f'Fix naming and homogenize final image for {fieldname}'
            step['dependency'] = n - 1
            step['id']         = 'HODSC' + code
            prefix  = CONTAINER_RUNNER + PYTHON3_CONTAINER + ' ' if USE_SINGULARITY else ''
            syscall  = prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_CAL_CHANNELSOUT} {diag_img_prefix}\n\n'
            syscall += prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {diag_img_prefix}'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        if cfg.WSC_POL != 'I':
            step = {}
            step['step']       = n
            step['comment']    = f'Make polarization intensity images for {fieldname}'
            step['dependency'] = n - 1
            step['id']         = 'MKLPI' + code
            syscall  = CONTAINER_RUNNER + PYTHON3_CONTAINER + ' ' if USE_SINGULARITY else ''
            syscall += f'python3 {cfg.TOOLS}/make_pol_images.py {cfg.IMAGES} {fieldname} True'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        all_steps.append((steps, kill_file, fieldname))

    # ------------------------------------------------------------------------------
    # Add split step at the very beginning (single CASA call for all fields)
    # ------------------------------------------------------------------------------

    split_step = {}
    split_step['step']       = 0
    split_step['comment']    = 'Split calibrator fields from main MS'
    split_step['dependency'] = None
    split_step['id']         = 'SPCAL'
    split_step['slurm_config'] = cfg.SLURM_DEFAULTS
    split_step['pbs_config']   = cfg.PBS_DEFAULTS
    syscall  = CONTAINER_RUNNER + CASA_CONTAINER + ' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_casa(casascript=cfg.TOOLS + '/casa_split_cals.py')
    split_step['syscall'] = syscall

    # ------------------------------------------------------------------------------
    # Write submit file
    # ------------------------------------------------------------------------------

    submit_file = 'submit_2GC_cal_jobs.sh'
    f = open(submit_file, 'w')
    f.write('#!/usr/bin/env bash\n')
    f.write('export SINGULARITY_BINDPATH=' + cfg.BINDPATH + '\n')

    f.write('\n#---------------------------------------\n')
    f.write('# Split calibrator fields\n')
    f.write('#---------------------------------------\n')
    run_command = gen.job_handler(
        syscall        = split_step['syscall'],
        jobname        = split_step['id'],
        infrastructure = INFRASTRUCTURE,
        dependency     = None,
        slurm_config   = split_step['slurm_config'],
        pbs_config     = split_step['pbs_config'])
    f.write('\n# ' + split_step['comment'] + '\n')
    f.write(run_command)

    for content in all_steps:
        steps     = content[0]
        kill_file = content[1]
        fieldname = content[2]
        id_list   = [split_step['id']]

        f.write('\n#---------------------------------------\n')
        f.write('# ' + fieldname)
        f.write('\n#---------------------------------------\n')

        for step in steps:
            nd = step['id']
            id_list.append(nd)

            # First imaging step depends on split completion.
            if step['dependency'] is None:
                dependency = split_step['id']
            else:
                dependency = steps[step['dependency']]['id']

            syscall = step['syscall']
            slurm_config = step.get('slurm_config', cfg.SLURM_DEFAULTS)
            pbs_config   = step.get('pbs_config',   cfg.PBS_DEFAULTS)

            run_command = gen.job_handler(
                syscall        = syscall,
                jobname        = nd,
                infrastructure = INFRASTRUCTURE,
                dependency     = dependency,
                slurm_config   = slurm_config,
                pbs_config     = pbs_config)

            f.write('\n# ' + step['comment'] + '\n')
            f.write(run_command)

        if INFRASTRUCTURE != 'node':
            f.write('\n# Generate kill script for ' + fieldname + '\n')
        if INFRASTRUCTURE in ('idia', 'hippo'):
            kill = 'echo "scancel "$' + '" "$'.join(id_list) + ' > ' + kill_file + '\n'
            f.write(kill)
        elif INFRASTRUCTURE == 'chpc':
            kill = 'echo "qdel "$' + '" "$'.join(id_list) + ' > ' + kill_file + '\n'
            f.write(kill)

    f.close()
    gen.make_executable(submit_file)

    gen.print_spacer()
    print(gen.col('Run file') + submit_file)
    gen.print_spacer()


if __name__ == "__main__":
    main()
