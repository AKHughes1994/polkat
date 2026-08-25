#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.a.uk

import glob
import json
import math
import os
import os.path as o
import sys
import yaml
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg


def get_adaptive_freq_intervals(yaml_path, n_model_channels):
    """Parse a QuartiCal YAML and return (ratio, overrides, zeros, ok) where:
      overrides — list of (term, original_fi, new_fi) for terms that need coarsening
      zeros     — list of term names with freq_interval=0 (whole-band, left alone)
      ok        — list of (term, fi) already at or above ratio (no change needed)
    ratio = ceil(PRE_NCHANS / n_model_channels)."""
    ratio = math.ceil(cfg.PRE_NCHANS / n_model_channels)
    overrides = []
    zeros = []
    ok = []
    with open(yaml_path) as f:
        qc = yaml.safe_load(f)
    for term in qc.get('solver', {}).get('terms', []):
        fi = int(str(qc.get(term, {}).get('freq_interval', 0)).strip())
        if fi == 0:
            zeros.append(term)
        elif fi < ratio:
            overrides.append((term, fi, ratio))
        else:
            ok.append((term, fi))
    return ratio, overrides, zeros, okmeasurements (of unique sources) from a single instrument to date. I will highlight some key results from the catalog, including a 


def main():

    USE_SINGULARITY = cfg.USE_SINGULARITY

    gen.preamble()
    print(gen.col()+'2GC (TRICOLOR flagging, imaging & DI phase self-calibration) setup')
    gen.print_spacer()


    # ------------------------------------------------------------------------------
    #
    # Setup paths, required containers, infrastructure
    #
    # ------------------------------------------------------------------------------


    OXKAT = cfg.OXKAT
    DATA = cfg.DATA
    IMAGES = cfg.IMAGES
    SCRIPTS = cfg.SCRIPTS
    TOOLS = cfg.TOOLS
    GAINTABLES = cfg.GAINTABLES
    LOGS = cfg.LOGS

    # ------------------------------------------------------------------------------
    #
    # Global flags and parameters
    #
    # ------------------------------------------------------------------------------

    # Set this flag to True to enable batch submission (4 at a time), False for no batching
    BATCH_FIELD_SUBMISSION = False

    # Set this flag to True to force all jobs to use 8 cores (useful for debugging/testing)
    FORCE_CORES_8 = True

    # Set to a directory path to use a single shared temp directory, or '' to use per-field temp directories
    OVERRIDE_TEMP_DIR = '/mnt/scratchssd'

    # Set this flag to False to disable all cleanup operations (find/delete/rm commands)
    ENABLE_CLEANUP = True

    # If True, read 'ref_ant' from project_info.json and pass the first antenna index
    # to QuartiCal via solver.reference_antenna. Helps stabilise phase solutions by
    # anchoring the gain phases to a single well-behaved antenna.
    SET_REFANT = True

    # If True, automatically scale each QuartiCal term's freq_interval so that it
    # covers at least as many channels as the ratio (PRE_NCHANS / model_channelsout).
    # This prevents QuartiCal being asked to solve on finer frequency bins than the
    # number of model channels written by wsclean, which would be under-constrained.
    ADAPTIVE_CHANNELS = True

    # If True, replace the standard full-sky predict + amplitude selfcal (PRCMI/CL2GC)
    # with a direction-dependent selfcal block (mirroring the 3GC_peel workflow).
    # Looks for DIR.reg in the CWD to define the calibration direction; extracts that
    # direction from the WSCMI inter model, populates a DIR1_DATA column in stage2_ms,
    # re-predicts the full sky model into MODEL_DATA, then solves with CAL_DDECAL_YAML.
    DO_DD_SELFCAL = False

    # Path to the QuartiCal YAML to use for the direction-dependent selfcal step above.
    # Should define G+dE (or equivalent DD) terms — distinct from CAL_2GC_YAML_COMPLEX.
    CAL_DDECAL_YAML = DATA + '/quartical/2GC_complex_2dir.yaml'

    # Set this flag to True to use custom parameters for the final pcal (WSCMA) image of targets
    PCAL_TARGET_IMAGE = False
    PCAL_TARGET_IMSIZE       = 8500    # Image size (pixels)
    PCAL_TARGET_CHANOUT      = 512     # Number of output channels
    PCAL_TARGET_MAXCHAN      = 128      # Max deconvolution channels (set equal to chanout = all channels)

    SELFCAL_MOD_DIR = DATA + '/selfcal_mod'

    # If True, skip imaging/self-calibrating the primary/polarization-angle/secondary calibrators
    # entirely and process targets only. Now set in oxkat/config.py (CAL_SKIP_CALS) so 1GC's
    # splitting step and RMSYNTH's extraction step can share the same setting.
    CAL_2GC_SKIP_CALS = cfg.CAL_SKIP_CALS

    gen.setup_dir(GAINTABLES)
    gen.setup_dir(IMAGES)
    gen.setup_dir(cfg.LOGS)
    gen.setup_dir(cfg.SCRIPTS)


    INFRASTRUCTURE, CONTAINER_PATH = gen.set_infrastructure(sys.argv)
    if CONTAINER_PATH is not None:
        CONTAINER_RUNNER='singularity exec '
    else:
        CONTAINER_RUNNER=''


    PYTHON3_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.PYTHON3_PATTERN,USE_SINGULARITY)
    CASA_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.CASA_PATTERN,USE_SINGULARITY)
    TRICOLOUR_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.TRICOLOUR_PATTERN,USE_SINGULARITY)
    WSCLEAN_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.WSCLEAN_PATTERN,USE_SINGULARITY)
    QUARTICAL_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.QUARTICAL_PATTERN,USE_SINGULARITY)

    # Get target and secondary information from project json

    with open('project_info.json') as f:
        project_info = json.load(f)
    
    band = project_info['band']
    target_ids = project_info['target_ids'] 
    target_names = project_info['target_names']
    target_ms = project_info['target_ms']
    pcal_names = project_info['secondary_names']
    pcal_ids = project_info['secondary_ids']
    pcal_ms = project_info['secondary_ms']
    bpcal_name = project_info['primary_name']
    bpcal_id = project_info['primary_id']
    primary_ms = project_info['primary_ms']
    pacal_name = project_info.get('polang_name', '')
    pacal_id = project_info.get('polang_id', '')
    pacal_ms = project_info.get('polang_ms', [])

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

    # Check if PRE_FIELDS was used to determine if filtering should be applied
    # working_names is always present, but only filter when PRE_FIELDS is specified
    working_names = project_info.get('working_names', [])
    apply_field_filter = cfg.PRE_FIELDS != ''
    
    # Build list of allowed field names for filtering if PRE_FIELDS is active
    if apply_field_filter:
        allowed_fields = set(working_names)
    else:
        allowed_fields = None

    # Build unified field list with MS files and flags
    field_list = []
    
    # Add targets (integrated across scans)
    for tt in range(len(target_ids)):
        field_name = target_names[tt]
        # Skip if PRE_FIELDS is active and this field is not in working_names
        if apply_field_filter and field_name not in allowed_fields:
            continue
        field_list.append({
            'name': field_name,
            'id': target_ids[tt],
            'ms': target_ms[tt],
            'is_target': True,
            'is_secondary': False
        })

    # Add primary calibrator (per scan) -- always processed, never skipped by CAL_SKIP_CALS
    # Skip if PRE_FIELDS is active and primary is not in working_names
    if not apply_field_filter or bpcal_name in allowed_fields:
        for myms in primary_ms:
            scan_str = myms.split('_scan')[-1].replace('.ms', '')
            field_list.append({
                'name': f"{bpcal_name}_scan{scan_str}",
                'id': bpcal_id,
                'ms': myms,
                'is_target': False,
                'is_secondary': False
            })

    # Add polarization angle calibrator (per scan) if it exists -- always processed,
    # never skipped by CAL_SKIP_CALS
    if pacal_name != '' and pacal_id != '' and len(pacal_ms) > 0:
        # Skip if PRE_FIELDS is active and pacal is not in working_names
        if not apply_field_filter or pacal_name in allowed_fields:
            for myms in pacal_ms:
                scan_str = myms.split('_scan')[-1].replace('.ms', '')
                field_list.append({
                    'name': f"{pacal_name}_scan{scan_str}",
                    'id': pacal_id,
                    'ms': myms,
                    'is_target': False,
                    'is_secondary': False
                })

    # Add secondaries (per scan) -- these are the only fields CAL_SKIP_CALS skips
    for ss in range(len(pcal_ids)):
        pcal_name = pcal_names[ss]
        # Skip if PRE_FIELDS is active and this secondary is not in working_names
        if apply_field_filter and pcal_name not in allowed_fields:
            continue
        for myms in pcal_ms[ss]:
            scan_str = myms.split('_scan')[-1].replace('.ms', '')
            field_list.append({
                'name': f"{pcal_name}_scan{scan_str}",
                'id': pcal_ids[ss],
                'ms': myms,
                'is_target': False,
                'is_secondary': True
            })

    # Skip secondary calibrator imaging/self-calibration -- targets, the primary, and the
    # polarization-angle calibrator are always processed regardless of this flag.
    if CAL_2GC_SKIP_CALS:
        n_cals = len([f for f in field_list if f['is_secondary']])
        field_list = [f for f in field_list if not f['is_secondary']]
        print(gen.col('2GC Cals')+f'CAL_2GC_SKIP_CALS is True -- skipping {n_cals} secondary field/scan(s)')

    # Determine if the blind/datamask images are going to be tapered
    if not cfg.WSC_TAPERMASK:
        tukeytaper = False
        minuvl = ''
        maxuvl = ''
        print(gen.col('UV Range')+f'Full UV range will be used for masking/self-calibration')
    else:
        tukeytaper = cfg.WSC_TUKEYTAPER
        minuvl = cfg.WSC_MINUVL
        maxuvl = cfg.WSC_MAXUVL
        print(gen.col('UV Range')+f'Restricted UV range will be used for masking/self-calibration: [{minuvl}, {maxuvl}]')   

    # ------------------------------------------------------------------------------
    # Calibrator (non-target) imaging parameters
    # These are used when is_target=False to reduce image size/channels for calibrators
    # ------------------------------------------------------------------------------
    CALIBRATOR_IMSIZE = 2560
    CALIBRATOR_MAXCHAN = 128
    CALIBRATOR_CHANNELSOUT = 64   

    # ------------------------------------------------------------------------------
    #
    # 2GC recipe definition
    #
    # ------------------------------------------------------------------------------


    freq_int_overrides_stage1 = []
    freq_int_overrides_stage2 = []
    freq_int_overrides_cal    = []
    if ADAPTIVE_CHANNELS:
        # Stage 2 YAML depends on whether DD selfcal is active
        stage2_yaml = CAL_DDECAL_YAML if DO_DD_SELFCAL else cfg.CAL_2GC_YAML_COMPLEX
        ratio1,   freq_int_overrides_stage1, zeros1,    ok1    = get_adaptive_freq_intervals(cfg.CAL_2GC_YAML,             cfg.WSC_DMASK_CHANNELSOUT)
        ratio2,   freq_int_overrides_stage2, zeros2,    ok2    = get_adaptive_freq_intervals(stage2_yaml,                  cfg.WSC_DMASK_CHANNELSOUT)
        ratio_cal, freq_int_overrides_cal,   zeros_cal, ok_cal = get_adaptive_freq_intervals(cfg.CAL_2GC_FILE,  CALIBRATOR_CHANNELSOUT)
        s2_label = 'Adaptive Chan (S2-DD)' if DO_DD_SELFCAL else 'Adaptive Chan (S2)'
        for label, yaml_name, ratio, overrides, zeros, ok in [
                ('Adaptive Chan (S1)',  o.basename(cfg.CAL_2GC_YAML),            ratio1,    freq_int_overrides_stage1, zeros1,    ok1),
                (s2_label,             o.basename(stage2_yaml),                  ratio2,    freq_int_overrides_stage2, zeros2,    ok2),
                ('Adaptive Chan (Cal)', o.basename(cfg.CAL_2GC_FILE), ratio_cal, freq_int_overrides_cal,   zeros_cal, ok_cal)]:
            parts = []
            parts += [f'{t}: 0 (all channels)'  for t in zeros]
            parts += [f'{t}: {old}→{new}'        for t, old, new in overrides]
            parts += [f'{t}: {fi} (no change)'   for t, fi in ok]
            print(gen.col(label)+f'{yaml_name} | ratio {ratio} — '+', '.join(parts))

    field_steps = []
    codes = []
    ii = 1
    stamp = gen.timenow()

    # Loop over all fields (targets and secondaries)
    for field in field_list:

        fieldname = field['name']
        myms = field['ms']
        is_target = field['is_target']

        steps = []        
        filename_fieldname = gen.scrub_target_name(fieldname)

        code = gen.get_target_code(fieldname)
        if code in codes:
            code += '_'+str(ii)
            ii += 1
        codes.append(code)
    
        # Define output parameters for 2GC step
        gain_outdir_2GC = GAINTABLES+'/2GC_'+str(filename_fieldname)+f'_{stamp}.qc/'
        log_outdir_2GC = LOGS+'/2GC_'+str(filename_fieldname)+f'_{stamp}.qc/'

        # Create source-specific subdirectory in IMAGES
        # All sources (targets and secondaries) get their own flat directory
        img_dir = f"{IMAGES}/{filename_fieldname}"
        gen.setup_dir(img_dir)
        
        # Create temp directory for wsclean intermediate files
        if OVERRIDE_TEMP_DIR == '':
            temp_dir = f"{IMAGES}/{filename_fieldname}_temp"
        else:
            temp_dir = f"{OVERRIDE_TEMP_DIR}/{filename_fieldname}_temp"
        gen.setup_dir(temp_dir)

        # Image prefixes
        img_prefix = img_dir+f'/img_{myms}_datablind'
        data_img_prefix = img_dir+f'/img_{myms}_datamask'
        inter_img_prefix = img_dir+f'/img_{myms}_intermask'
        final_img_prefix = img_dir+f'/img_{myms}_finalmask'
        pcal_img_prefix = img_dir+f'/img_{myms}_pcalmask'
        fulluv_img_prefix = img_dir+f'/img_{myms}_fulluv'

        # Field-specific kill file
        kill_file = SCRIPTS+'/kill_2GC_jobs_'+filename_fieldname+'.sh'

        gen.print_spacer()
        field_type = 'Target' if is_target else 'Secondary'
        print(gen.col(field_type)+fieldname)
        print(gen.col('Measurement Set')+myms)
        print(gen.col('Code')+code)

        n = 0
        
        # Run tricolour only for targets, skip for secondaries
        if is_target:
            step = {}
            step['step'] = n
            step['comment'] = 'Run Tricolour on '+myms
            step['dependency'] = None
            step['id'] = 'TRILE'+code
            step['glam_config'] = cfg.GLAM_COREHEAVY
            step['slurm_config'] = cfg.SLURM_TRICOLOUR
            step['pbs_config'] = cfg.PBS_TRICOLOUR
            syscall = CONTAINER_RUNNER+TRICOLOUR_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_tricolour(myms = myms,
                    config = DATA+'/tricolour/target_flagging_1_narrow.yaml',
                    datacol = 'DATA',
                    strategy = 'polarisation')
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        
        # Check for mask file; mask can be a string path or a list of paths.
        # If a list, do a case-insensitive search for the field name in each entry.
        mask = cfg.WSC_MASK
        user_specified_mask = False
        if mask:
            if isinstance(mask, list):
                matched = next((m for m in mask if fieldname.lower() in m.lower()), None)
                if matched:
                    mask = matched
                    user_specified_mask = True
                else:
                    mask = None  # No match found; fall back to blind imaging
            else:
                # Check if field name is in the mask path for single string masks
                if fieldname.lower() in mask.lower():
                    user_specified_mask = True
                else:
                    mask = None  # No match; fall back to blind imaging

        if not mask:

            step = {}
            step['step'] = n
            step['comment'] = 'Shallow blind wsclean on DATA column for source {}'.format(fieldname)
            step['dependency'] = n - 1 if is_target else None  # Depends on tricolour for targets only
            step['id'] = 'WSDBL'+code
            step['glam_config'] = cfg.GLAM_MEDIUM
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            # Adjust parameters based on field type
            imsize_blind = CALIBRATOR_IMSIZE if not is_target else cfg.WSC_IMSIZE
            maxchan_blind = CALIBRATOR_MAXCHAN if not is_target else cfg.WSC_MAX_CHANNELS
            weight_blind = 'uniform' if not is_target else cfg.WSC_WEIGHT
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                    tempdir = temp_dir,
                    imgname = img_prefix,
                    datacol = 'DATA',
                    chanout = cfg.WSC_BLIND_CHANNELSOUT,
                    nomodel = True,
                    pol = 'I',
                    intervalsout = False,
                    mfweight = True,
                    localrms = cfg.WSC_LOCALRMS_BLIND,
                    automask = cfg.WSC_AUTOMASK_BLIND,
                    autothreshold = cfg.WSC_AUTOTHRESHOLD_BLIND,
                    tukeytaper=tukeytaper,
                    minuvl = minuvl,
                    maxuvl = maxuvl,
                    mask = mask,
                    weight=weight_blind,
                    imsize = imsize_blind,
                    maxchan = maxchan_blind,
                    cores = cores,
                    absmem = absmem)
            for call in imcall: 
                syscall += prefix + call + '\n\n'
            # Clean up intermediate images (keep only MFS-image.fits files)
            if ENABLE_CLEANUP:
                syscall += f'find {img_dir} -type f -name "*{os.path.basename(img_prefix)}*" ! -name "*MFS-image.fits" -delete\n\n'
            step['syscall'] = syscall
            step['grouping'] = None
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Make cleaning mask for ' + fieldname
            step['dependency'] = n - 1
            step['id'] = 'MASK0'+code
            step['glam_config'] = cfg.GLAM_SMALL
            syscall  = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            # Use threshold=10.0 for calibrators (non-targets)
            breizorro_thresh = 10.0 if not is_target else cfg.BREIZORRO_THRESH
            syscall += gen.generate_syscall_breizorro(restoredimage = f"{img_prefix}-MFS-image.fits", outfile = f"{img_prefix}-MFS-image.mask.fits", thresh=breizorro_thresh)[0]
            step['syscall'] = syscall
            steps.append(step)
            n += 1
            mask = f"{img_prefix}-MFS-image.mask.fits"
            print(gen.col('Mask')+ 'None')

        else:
            print(gen.col('Mask')+mask)

        # Default finalmask to current mask (used for final WSCMA image)
        finalmask = mask

        step = {}
        step['step'] = n
        step['comment'] = 'Run wsclean, masked deconvolution of the DATA column for source {}'.format(fieldname)
        step['dependency'] = n - 1 
        step['id'] = 'WSDMA'+code
        step['glam_config'] = cfg.GLAM_WSC_SOURCE if is_target else cfg.GLAM_WSC_CAL
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
        syscall = ''
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        imsize_dmask = CALIBRATOR_IMSIZE if not is_target else cfg.WSC_IMSIZE
        maxchan_dmask = CALIBRATOR_MAXCHAN if not is_target else cfg.WSC_MAX_CHANNELS
        chanout_dmask = CALIBRATOR_CHANNELSOUT if not is_target else cfg.WSC_DMASK_CHANNELSOUT
        automask_dmask = 2.0 * cfg.WSC_SHALLOWMASK if not is_target else cfg.WSC_SHALLOWMASK
        weight_dmask = 'uniform' if not is_target else cfg.WSC_WEIGHT
        imcall = gen.generate_syscall_wsclean(mslist = [myms],
                tempdir = temp_dir,
                imgname = data_img_prefix,
                mfweight = False,
                datacol = 'DATA',
                mask = mask,
                chanout = chanout_dmask,
                intervalsout = False,
                tukeytaper=tukeytaper,
                minuvl = minuvl,
                maxuvl = maxuvl,
                nomodel = True,
                sourcelist = False,
                automask = automask_dmask,
                localrms = cfg.WSC_SHALLOWMASK_LOCALRMS,
                autothreshold = cfg.WSC_INTER_AUTOTHRESHOLD,
                weight=weight_dmask,
                imsize = imsize_dmask,
                maxchan = maxchan_dmask,
                cores = cores,
                absmem = absmem)
        for call in imcall: 
            syscall += prefix + call + '\n\n'
            if ENABLE_CLEANUP:
                syscall += f'find {img_dir} -type f -name "*{os.path.basename(data_img_prefix)}*" ! -name "*model.fits" ! -name "*MFS*image.fits" -delete\n\n'
        step['syscall'] = syscall
        step['grouping'] = 2 if ENABLE_CLEANUP else None
        steps.append(step)
        n += 1

        if is_target and cfg.WSC_MAX_CHANNELS < cfg.WSC_DMASK_CHANNELSOUT:
            step = {}
            step['step'] = n
            step['comment'] = f'Fix image naming for DATAMASK images'
            step['dependency'] = n - 1
            step['id'] = 'FXDMA'+code
            step['glam_config'] = cfg.GLAM_SMALL
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall = prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_DMASK_CHANNELSOUT} {data_img_prefix}'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Run wsclean-predict of DATAMASK model for source {}'.format(fieldname)
        step['dependency'] = n - 1
        step['id'] = 'PRDMA'+code
        step['glam_config'] = cfg.GLAM_MEDIUM
        step['slurm_config'] = cfg.SLURM_PREDICT
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + data_img_prefix + '\n\n'
        if cfg.MOD_MODEL_SELFCAL and cfg.WSC_POL != 'I':
            if is_target:
                # Look for a per-target .reg file; fall back to current mask FITS
                target_reg = f'{SELFCAL_MOD_DIR}/{filename_fieldname}.reg'
                spatial_arg = target_reg if os.path.exists(target_reg) else mask
                syscall += prefix + (f'python3 {TOOLS}/mod_model_selfcal.py '
                                     f'--identifier {data_img_prefix} --stokes V '
                                     f'--spatial {spatial_arg}\n\n')  # --zero-all-neg-I
            else:
                # Calibrators: V-only, use their mask as spatial input
                syscall += prefix + (f'python3 {TOOLS}/mod_model_selfcal.py '
                                     f'--identifier {data_img_prefix} --stokes V '
                                     f'--spatial {mask}\n\n')
        syscall += prefix + gen.generate_syscall_predict(msname = myms,
                imgname = data_img_prefix,
                cores = cores,
                absmem = absmem,
                tempdir = temp_dir,
                chanout = chanout_dmask)
        if ENABLE_CLEANUP:
            syscall += '\n\nfind ' + img_dir + ' -type f -name "*' + os.path.basename(data_img_prefix) + '*" ! -name "*MFS*image.fits" ! -name "*MFS*model.fits" -delete'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        # Stage 1: phase-only self-calibration on original MS
        step = {}
        step['step'] = n
        step['comment'] = 'Run Quartical phase-only self-calibration (stage 1) on {}'.format(fieldname)
        step['dependency'] = n - 1
        step['id'] = 'C02G2'+code
        step['glam_config'] = cfg.GLAM_COREHEAVY
        cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
        syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
        extra_args = f'output.gain_directory={gain_outdir_2GC} output.log_directory={log_outdir_2GC}'
        if not is_target and not cfg.CAL_1GC_APPLYPARANG:
            # Calibrators have no stage 2 -- this is their last (only) 2GC calibration step,
            # so rotate to the sky frame here if 1GC left the parang correction for later.
            extra_args += ' output.apply_p_jones_inv=true'
        if ref_ant_arg is not None:
            extra_args += f' solver.reference_antenna={ref_ant_arg}'
        # Apply adaptive freq_interval overrides: target        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
s use stage1 YAML, non-targets use calibrator YAML
        fi_overrides = freq_int_overrides_stage1 if is_target else freq_int_overrides_cal
        for term, _, new_fi in fi_overrides:
            extra_args += f' {term}.freq_interval={new_fi}'
        if minuvl != '' or maxuvl != '':
            min_uv = minuvl if minuvl != '' else '0'
            max_uv = maxuvl if maxuvl != '' else '0'
            extra_args += f' input_ms.select_uv_range=[{min_uv},{max_uv}]'
        cal_yaml = cfg.CAL_2GC_FILE if not is_target else cfg.CAL_2GC_YAML
        syscall += gen.generate_syscall_quartical(yaml = cal_yaml,
                myms = myms,
                extra_args = extra_args)
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        if is_target:

            imsize_inter = cfg.WSC_IMSIZE
            maxchan_inter = cfg.WSC_MAX_CHANNELS
            chanout_inter = cfg.WSC_DMASK_CHANNELSOUT
            weight_inter = cfg.WSC_WEIGHT

            # WSCMI: image CORRECTED_DATA of original MS after phase selfcal
            step = {}
            step['step'] = n
            step['comment'] = 'Run wsclean on CORRECTED_DATA (stage 1 phase selfcal) for {}'.format(fieldname)
            step['dependency'] = n - 1
            step['id'] = 'WSCMI'+code
            step['glam_config'] = cfg.GLAM_WSC_SOURCE
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                    tempdir = temp_dir,
                    imgname = inter_img_prefix,
                    datacol = 'CORRECTED_DATA',
                    mask = mask,
                    chanout = chanout_inter,
                    tukeytaper=tukeytaper,
                    minuvl = minuvl,
                    maxuvl = maxuvl,
                    cores = cores,
                    nomodel=True,
                    sourcelist = False,
                    qu_automask_scale = 1.0,
                    localrms = cfg.WSC_INTER_LOCALRMS,
                    automask = cfg.WSC_INTER_AUTOMASK,
                    autothreshold = cfg.WSC_INTER_AUTOTHRESHOLD,
                    weight = weight_inter,
                    imsize = imsize_inter,
                    maxchan = maxchan_inter,
                    absmem = absmem)
            for call in imcall:
                syscall += prefix + call + '\n\n'
                if ENABLE_CLEANUP:
                    syscall += f'find {img_dir} -type f -name "*{os.path.basename(inter_img_prefix)}*" ! -name "*model.fits" ! -name "*MFS*image.fits" -delete\n\n'
            step['syscall'] = syscall
            step['grouping'] = 2 if ENABLE_CLEANUP else None
            steps.append(step)
            n += 1

            if cfg.WSC_MAX_CHANNELS < cfg.WSC_DMASK_CHANNELSOUT:
                step = {}
                step['step'] = n
                step['comment'] = f'Fix image naming for INTER images'
                step['dependency'] = n - 1
                step['id'] = 'FXCMI'+code
                step['glam_config'] = cfg.GLAM_SMALL
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall = prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_DMASK_CHANNELSOUT} {inter_img_prefix}'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            # Split CORRECTED_DATA → stage 2 MS (DATA = phase-selfcal'd visibilities)
            stage2_ms = myms.replace('.ms', '_stage2.ms')

            step = {}
            step['step'] = n
            step['comment'] = 'Split CORRECTED_DATA into stage 2 MS for {}'.format(fieldname)
            step['dependency'] = n - 1
            step['id'] = 'SPLCD'+code
            step['glam_config'] = cfg.GLAM_CASA
            syscall = CONTAINER_RUNNER + CASA_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_casa(casascript=f'{cfg.TOOLS}/casa_split_corrected_ms.py {myms} {stage2_ms}')
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            gain_outdir_stage2 = GAINTABLES+'/2GC_'+str(filename_fieldname)+f'_{stamp}_stage2.qc/'
            log_outdir_stage2  = LOGS+'/2GC_'+str(filename_fieldname)+f'_{stamp}_stage2.qc/'

            if DO_DD_SELFCAL and is_target:

                # Direction-dependent selfcal block (mirrors 3GC_peel workflow).
                # DIR.reg in the CWD defines the calibration direction to extract.
                dir_region     = 'DIR.reg'
                dir_img_prefix = inter_img_prefix + '-' + dir_region.split('.')[0]
                dd_recipe      = 'MODEL_DATA~DIR1_DATA:DIR1_DATA'
                dd_subtract    = '[1]'

                # Single combined step: fix NaN → extract DIR → add column →
                # predict DIR into MODEL_DATA → copy to DIR1_DATA → re-predict full sky model
                step = {}
                step['step'] = n
                step['comment'] = f'DD prep for {fieldname}: fix NaN, extract DIR, add column, predict DIR, copy, re-predict full model'
                step['dependency'] = n - 1
                step['id'] = 'DDPRE'+code
                step['glam_config'] = cfg.GLAM_MEDIUM
                step['slurm_config'] = cfg.SLURM_PREDICT
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step, INFRASTRUCTURE, cfg.WSC_ABSMEM)
                cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
                py_pfx  = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                wsc_pfx = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' '  if USE_SINGULARITY else ''
                syscall  = py_pfx  + 'python3 ' + TOOLS + '/fix_nan_models.py ' + inter_img_prefix + '\n\n'
                if cfg.MOD_MODEL_SELFCAL and cfg.WSC_POL != 'I':
                    target_reg = f'{SELFCAL_MOD_DIR}/{filename_fieldname}.reg'
                    spatial_arg = target_reg if os.path.exists(target_reg) else mask
                    syscall += py_pfx + (f'python3 {TOOLS}/mod_model_selfcal.py '
                                         f'--identifier {inter_img_prefix} --stokes V '
                                         f'--spatial {spatial_arg}\n\n')
                syscall += py_pfx  + 'python3 ' + OXKAT + '/3GC_split_model_images.py '
                syscall +=           '--region ' + dir_region + ' --prefix ' + inter_img_prefix + '\n\n'
                syscall += py_pfx  + 'python3 ' + TOOLS + '/add_MS_column.py '
                syscall +=           '--colname DIR1_DATA ' + stage2_ms + '\n\n'
                syscall += wsc_pfx + gen.generate_syscall_predict(msname = stage2_ms,
                        imgname  = dir_img_prefix,
                        cores    = cores,
                        absmem   = absmem,
                        tempdir  = temp_dir,
                        chanout  = chanout_inter) + '\n\n'
                syscall += py_pfx  + 'python3 ' + TOOLS + '/copy_MS_column.py '
                syscall +=           '--fromcol MODEL_DATA --tocol DIR1_DATA ' + stage2_ms + '\n\n'
                syscall += wsc_pfx + gen.generate_syscall_predict(msname = stage2_ms,
                        imgname  = inter_img_prefix,
                        cores    = cores,
                        absmem   = absmem,
                        tempdir  = temp_dir,
                        chanout  = chanout_inter)
                step['syscall'] = syscall
                steps.append(step)
                n += 1

                step = {}
                step['step'] = n
                step['comment'] = 'Run QuartiCal DD self-calibration (stage 2) on {}'.format(fieldname)
                step['dependency'] = n - 1
                step['id'] = 'CL2GC'+code
                step['glam_config'] = cfg.GLAM_COREHEAVY
                cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
                syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
                extra_args = (f'output.gain_directory={gain_outdir_stage2} output.log_directory={log_outdir_stage2}'
                              f' input_model.recipe={dd_recipe} output.subtract_directions={dd_subtract}')
                if not cfg.CAL_1GC_APPLYPARANG:
                    # Last 2GC calibration step for targets -- rotate to the sky frame here
                    # if 1GC left the parang correction for later.
                    extra_args += ' output.apply_p_jones_inv=true'
                if ref_ant_arg is not None:
                    extra_args += f' solver.reference_antenna={ref_ant_arg}'
                for term, _, new_fi in freq_int_overrides_stage2:
                    extra_args += f' {term}.freq_interval={new_fi}'
                if minuvl != '' or maxuvl != '':
                    extra_args += f' input_ms.select_uv_range=[{min_uv},{max_uv}]'
                syscall += gen.generate_syscall_quartical(yaml = CAL_DDECAL_YAML,
                        myms = stage2_ms,
                        extra_args = extra_args)
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            else:

                step = {}
                step['step'] = n
                step['comment'] = 'Run wsclean-predict of INTER model into stage 2 MS for {}'.format(fieldname)
                step['dependency'] = n - 1
                step['id'] = 'PRCMI'+code
                step['glam_config'] = cfg.GLAM_MEDIUM
                step['slurm_config'] = cfg.SLURM_PREDICT
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
                prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + inter_img_prefix + '\n\n'
                if cfg.MOD_MODEL_SELFCAL and cfg.WSC_POL != 'I':
                    target_reg = f'{SELFCAL_MOD_DIR}/{filename_fieldname}.reg'
                    spatial_arg = target_reg if os.path.exists(target_reg) else mask
                    syscall += prefix + (f'python3 {TOOLS}/mod_model_selfcal.py '
                                         f'--identifier {inter_img_prefix} --stokes V '
                                         f'--spatial {spatial_arg}\n\n')
                syscall += prefix + gen.generate_syscall_predict(msname = stage2_ms,
                        imgname = inter_img_prefix,
                        cores = cores,
                        absmem = absmem,
                        tempdir = temp_dir,
                        chanout = chanout_inter)
                if ENABLE_CLEANUP:
                    syscall += '\n\nfind ' + img_dir + ' -type f -name "*' + os.path.basename(inter_img_prefix) + '*" ! -name "*MFS*image.fits" ! -name "*MFS*model.fits" -delete'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

                step = {}
                step['step'] = n
                step['comment'] = 'Run Quartical amplitude self-calibration (stage 2) on {}'.format(fieldname)
                step['dependency'] = n - 1
                step['id'] = 'CL2GC'+code
                step['glam_config'] = cfg.GLAM_COREHEAVY
                cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
                syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
                extra_args = f'output.gain_directory={gain_outdir_stage2} output.log_directory={log_outdir_stage2}'
                if not cfg.CAL_1GC_APPLYPARANG:
                    # Last 2GC calibration step for targets -- rotate to the sky frame here
                    # if 1GC left the parang correction for later.
                    extra_args += ' output.apply_p_jones_inv=true'
                if ref_ant_arg is not None:
                    extra_args += f' solver.reference_antenna={ref_ant_arg}'
                for term, _, new_fi in freq_int_overrides_stage2:
                    extra_args += f' {term}.freq_interval={new_fi}'
                if minuvl != '' or maxuvl != '':
                    extra_args += f' input_ms.select_uv_range=[{min_uv},{max_uv}]'
                syscall += gen.generate_syscall_quartical(yaml = cfg.CAL_2GC_YAML_COMPLEX,
                        myms = stage2_ms,
                        extra_args = extra_args)
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            myms = stage2_ms

        step = {}
        step['step'] = n
        step['comment'] = f'Run wsclean, final masked deconvolution of the CORRECTED_DATA (self-calibrated) for {fieldname}'
        step['dependency'] = n - 1
        step['id'] = 'WSCMA'+code
        step['glam_config'] = cfg.GLAM_WSC_SOURCE if is_target else cfg.GLAM_WSC_CAL
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
        syscall = ''
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        imsize_pcal  = cfg.WSC_IMSIZE if is_target else CALIBRATOR_IMSIZE
        maxchan_pcal = cfg.WSC_MAX_CHANNELS if is_target else CALIBRATOR_MAXCHAN
        chanout_pcal = cfg.WSC_PCAL_CHANNELSOUT if is_target else CALIBRATOR_CHANNELSOUT
        weight_pcal  = cfg.WSC_WEIGHT if is_target else 'uniform'
        if PCAL_TARGET_IMAGE and is_target:
            imsize_pcal  = PCAL_TARGET_IMSIZE
            chanout_pcal = PCAL_TARGET_CHANOUT
            maxchan_pcal = PCAL_TARGET_MAXCHAN
            print(gen.col() + f"PCAL image will be made for target {fieldname} with imsize={imsize_pcal}, chanout={chanout_pcal}, maxchan={maxchan_pcal}")
        imcall = gen.generate_syscall_wsclean(mslist = [myms],
                tempdir = temp_dir,
                imgname = pcal_img_prefix,
                datacol = 'CORRECTED_DATA',
                mask = finalmask,
                chanout = chanout_pcal,
                nomodel=True,
                sourcelist = False,
                tukeytaper=tukeytaper,
                minuvl = minuvl,
                maxuvl = maxuvl,
                weight = weight_pcal,
                imsize = imsize_pcal,
                maxchan = maxchan_pcal,
                cores = cores,
                absmem = absmem)
        for call in imcall: 
            syscall += prefix + call + '\n\n'
            if ENABLE_CLEANUP:
                syscall += f'rm -f {pcal_img_prefix}*dirty*\n'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        # Homogenize beams if maxchan_pcal is greater than WSC_MAX_CHANNELS
        if maxchan_pcal < chanout_pcal:
            step = {}
            step['step'] = n
            step['comment'] = f'Homogenize the PCAL resolution across frequency channels'
            step['dependency'] = n - 1
            step['id'] = 'HOCMA' + code
            step['glam_config'] = cfg.GLAM_LARGE
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {chanout_pcal} {pcal_img_prefix}\n\n'
            syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {pcal_img_prefix}\n\n'
            if ENABLE_CLEANUP:
                syscall += f'find . -name "{os.path.basename(pcal_img_prefix)}*homogenized.fits" ! -name "*MFS*" -delete'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        step = {}
        step['step'] = n
        step['comment'] = f'Run wsclean, full UV range (no tukey taper / minuv-l cut), 16-channel Stokes I deconvolution of the CORRECTED_DATA (self-calibrated) for {fieldname}'
        step['dependency'] = n - 1
        step['id'] = 'WSFUV'+code
        step['glam_config'] = cfg.GLAM_WSC_SOURCE if is_target else cfg.GLAM_WSC_CAL
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        cores = '8' if FORCE_CORES_8 else step['glam_config']['CPUS']
        syscall = ''
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        imcall = gen.generate_syscall_wsclean(mslist = [myms],
                tempdir = temp_dir,
                imgname = fulluv_img_prefix,
                datacol = 'CORRECTED_DATA',
                mask = finalmask,
                chanout = 16,
                nomodel = True,
                sourcelist = False,
                pol = 'I',
                tukeytaper = False,
                minuvl = '0',
                maxuvl = '',
                weight = weight_pcal,
                imsize = imsize_pcal,
                maxchan = maxchan_pcal,
                cores = cores,
                absmem = absmem)
        for call in imcall:
            syscall += prefix + call + '\n\n'
            if ENABLE_CLEANUP:
                syscall += f'rm -f {fulluv_img_prefix}*dirty*\n'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        if cfg.WSC_POL != 'I':
            step = {}
            step['step'] = n
            step['comment'] = 'Make Polarization Intensity Images'
            step['dependency'] = n - 1
            step['id'] = 'MKLPI'+code
            step['glam_config'] = cfg.GLAM_SMALL
            syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += f"python3 {cfg.TOOLS}/make_pol_images.py {img_dir} {filename_fieldname} false"
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Create zoom cutouts for '+fieldname
        step['dependency'] = n - 1
        step['id'] = 'CTOUT'+code
        step['glam_config'] = cfg.GLAM_SMALL
        syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        # Pass the relative path from IMAGES directory for make_cutouts.py
        cutout_path = img_dir.replace(IMAGES+'/', '')
        syscall += f"python3 {cfg.TOOLS}/make_cutouts.py {cutout_path}"
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Clean directory for '+fieldname
        step['dependency'] = n - 1
        step['id'] = 'CLEAN'+code
        step['glam_config'] = cfg.GLAM_SMALL
        syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += f"python3 {cfg.TOOLS}/clean_directory.py {filename_fieldname}\n\n"
        if ENABLE_CLEANUP:
            syscall += f"rm -rf {temp_dir}"
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        field_type_str = 'target' if is_target else 'secondary'
        field_steps.append((steps,kill_file,fieldname,field_type_str))


    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    # Batch submission control variables
    ENABLE_BATCH_FIELD_SUBMISSION = False  # Set to True to enable batching
    FIELD_BATCH_SIZE = 5  # Number of sources to process in parallel if batching is enabled

    submit_file = 'submit_2GC_jobs.sh'

    f = open(submit_file,'w')
    f.write('#!/usr/bin/env bash\n')
    f.write('export SINGULARITY_BINDPATH='+cfg.BINDPATH+'\n')

    if INFRASTRUCTURE == 'glam':
        # Use Glamdring-specific submission script
        f.write('set -euo pipefail\n\n')

        # First generate the individual job scripts for all fields
        for content in field_steps:
            steps = content[0]
            fieldname = content[2]
            for step in steps:
                nd = step['id']
                syscall = step['syscall']
                if 'glam_config' in step.keys():
                    glam_config = step['glam_config']
                else:
                    glam_config = cfg.GLAM_STANDARD
                grouping = step.get('grouping', None)
                gen.job_handler(syscall = syscall,
                                jobname = nd,
                                infrastructure = INFRASTRUCTURE,
                                glam_config = glam_config,
                                grouping = grouping)

        # Batch logic for field submission
        if ENABLE_BATCH_FIELD_SUBMISSION:
            batch_size = FIELD_BATCH_SIZE
            num_fields = len(field_steps)
            prev_batch_last_job = None
            for batch_start in range(0, num_fields, batch_size):
                batch = field_steps[batch_start:batch_start+batch_size]
                for i, content in enumerate(batch):
                    steps = content[0]
                    fieldname = content[2]
                    field_type = content[3]

                    f.write('\n#---------------------------------------\n')
                    f.write('# '+fieldname+f' ({field_type})')
                    f.write('\n#---------------------------------------\n')

                    # If this is the first job in the batch and not the first batch, set dependency
                    if i == 0 and prev_batch_last_job is not None:
                        steps[0]['dependency'] = prev_batch_last_job
                    gen.write_glam_submission_script(f, steps)

                # Remember the last job id in this batch for dependency
                prev_batch_last_job = batch[-1][0][-1]['id']
        else:
            # Now write the submission commands for each field (no batching)
            for content in field_steps:
                steps = content[0]
                fieldname = content[2]
                field_type = content[3]

                f.write('\n#---------------------------------------\n')
                f.write('# '+fieldname+f' ({field_type})')
                f.write('\n#---------------------------------------\n')
                gen.write_glam_submission_script(f, steps)
        
    else:
        # Original behavior for other infrastructures
        for content in field_steps:  
            steps = content[0]
            kill_file = content[1]
            fieldname = content[2]
            field_type = content[3]
            id_list = []

            f.write('\n#---------------------------------------\n')
            f.write('# '+fieldname+f' ({field_type})')
            f.write('\n#---------------------------------------\n')

            for step in steps:

                nd = step['id']
                id_list.append(nd)
                if step['dependency'] is not None:
                    dependency = steps[step['dependency']]['id']
                else:
                    dependency = None
                syscall = step['syscall']
                if 'slurm_config' in step.keys():
                    slurm_config = step['slurm_config']
                else:
                    slurm_config = cfg.SLURM_DEFAULTS
                if 'pbs_config' in step.keys():
                    pbs_config = step['pbs_config']
                else:
                    pbs_config = cfg.PBS_DEFAULTS
                comment = step['comment']

                run_command = gen.job_handler(syscall = syscall,
                                jobname = nd,
                                infrastructure = INFRASTRUCTURE,
                                dependency = dependency,
                                slurm_config = slurm_config,
                                pbs_config = pbs_config)


                f.write('\n# '+comment+'\n')
                f.write(run_command)

            if INFRASTRUCTURE != 'node':
                f.write('\n# Generate kill script for '+fieldname+'\n')
            if INFRASTRUCTURE == 'idia' or INFRASTRUCTURE == 'hippo':
                kill = 'echo "scancel "$'+'" "$'.join(id_list)+' > '+kill_file+'\n'
                f.write(kill)
            elif INFRASTRUCTURE == 'chpc':
                kill = 'echo "qdel "$'+'" "$'.join(id_list)+' > '+kill_file+'\n'
                f.write(kill)

        
    f.close()

    gen.make_executable(submit_file)

    gen.print_spacer()
    print(gen.col('Run file')+submit_file)
    gen.print_spacer()

    # ------------------------------------------------------------------------------



if __name__ == "__main__":


    main()
