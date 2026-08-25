#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.a.uk

import glob
import json
import math
import os.path as o
import sys
import yaml
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg


SET_REFANT       = True
ADAPTIVE_CHANNELS = True
DO_DD_SELFCAL    = False

# Never True when CAL_1GC_APPLYPARANG is True — 1GC already puts CORRECTED_DATA
# in the sky frame, so re-applying parang here would double-correct it.
PARANGMODEL = cfg.CAL_2GC_PARANGMODEL and not cfg.CAL_1GC_APPLYPARANG


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
    return ratio, overrides, zeros, ok


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

    SELFCAL_MOD_DIR = DATA + '/selfcal_mod'

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

    # Get target information from project json

    with open('project_info.json') as f:
        project_info = json.load(f)
    
    band = project_info['band']
    target_ids = project_info['target_ids']
    target_names = project_info['target_names']
    target_ms = project_info['target_ms']

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
    #
    # 2GC recipe definition
    #
    # ------------------------------------------------------------------------------


    freq_int_overrides_stage1 = []
    freq_int_overrides_stage2 = []
    if ADAPTIVE_CHANNELS:
        # Stage 2 YAML depends on whether DD selfcal is active
        stage2_yaml = CAL_DDECAL_YAML if DO_DD_SELFCAL else cfg.CAL_2GC_YAML_COMPLEX
        ratio1, freq_int_overrides_stage1, zeros1, ok1 = get_adaptive_freq_intervals(cfg.CAL_2GC_YAML, cfg.WSC_DMASK_CHANNELSOUT)
        ratio2, freq_int_overrides_stage2, zeros2, ok2 = get_adaptive_freq_intervals(stage2_yaml,       cfg.WSC_DMASK_CHANNELSOUT)
        s2_label = 'Adaptive Chan (S2-DD)' if DO_DD_SELFCAL else 'Adaptive Chan (S2)'
        for label, yaml_name, ratio, overrides, zeros, ok in [
                ('Adaptive Chan (S1)', o.basename(cfg.CAL_2GC_YAML), ratio1, freq_int_overrides_stage1, zeros1, ok1),
                (s2_label,            o.basename(stage2_yaml),       ratio2, freq_int_overrides_stage2, zeros2, ok2)]:
            parts = []
            parts += [f'{t}: 0 (all channels)'   for t in zeros]
            parts += [f'{t}: {old}→{new}'         for t, old, new in overrides]
            parts += [f'{t}: {fi} (no change)'    for t, fi in ok]
            print(gen.col(label)+f'{yaml_name} | ratio {ratio} — '+', '.join(parts))

    target_steps = []
    codes = []
    ii = 1
    stamp = gen.timenow()

    # Loop over targets
    for tt in range(0,len(target_ids)):

        targetname   = target_names[tt]
        myms = target_ms[tt]

        if targetname not in project_info['working_names']:

            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('MS')+'not found, skipping')

        else:

            targetindex   = str(project_info['working_ids'][project_info['working_names'].index(targetname)])
            save_flag = f'after_pcal_{targetname}'
            
            steps = []        
            filename_targetname = gen.scrub_target_name(targetname)

            code = gen.get_target_code(targetname)
            if code in codes:
                code += '_'+str(ii)
                ii += 1
            codes.append(code)
        
            # Define output parameters for 2GC step
            gain_outdir_2GC    = GAINTABLES+'/2GC_'+str(filename_targetname)+f'_{stamp}.qc/'
            log_outdir_2GC     = LOGS+'/2GC_'+str(filename_targetname)+f'_{stamp}.qc/'
            gain_outdir_stage2 = GAINTABLES+'/2GC_'+str(filename_targetname)+f'_{stamp}_stage2.qc/'
            log_outdir_stage2  = LOGS+'/2GC_'+str(filename_targetname)+f'_{stamp}_stage2.qc/'

            # Stage 2 MS (split from CORRECTED_DATA after stage 1 phase selfcal)
            stage2_ms = myms.replace('.ms', '_stage2.ms')

            # Image prefixes
            img_prefix = IMAGES+f'/img_{myms}_datablind'
            data_img_prefix = IMAGES+f'/img_{myms}_datamask'
            inter_img_prefix = IMAGES+f'/img_{myms}_intermask'
            pcal_img_prefix = IMAGES+f'/img_{myms}_pcalmask'
            uniform_img_prefix = IMAGES+f'/img_{myms}_uniform'

            # Target-specific kill file
            kill_file = SCRIPTS+'/kill_2GC_jobs_'+filename_targetname+'.sh'

            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('Measurement Set')+myms)
            print(gen.col('Code')+code)
            mod_regions = sorted(glob.glob(f'{SELFCAL_MOD_DIR}/*{filename_targetname}*.reg'))
            target_mod_region = mod_regions[0] if mod_regions else None
            if cfg.MOD_MODEL_SELFCAL and cfg.WSC_POL != 'I':
                if target_mod_region is not None:
                    print(gen.col('Selfcal MOD Region')+f'{target_mod_region} (using first of {len(mod_regions)} match(es))')
                else:
                    print(gen.col('Selfcal MOD Region')+'None found (will fallback to imaging mask)')

            n = 0
            step = {}
            step['step'] = n
            step['comment'] = 'Run Tricolour on '+myms
            step['dependency'] = None
            step['id'] = 'TRILE'+code
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

            
            # Check for mask file, if doesn't exist make one
            mask = cfg.WSC_MASK
            if not mask:

                step = {}
                step['step'] = n
                step['comment'] = 'Shallow blind wsclean on CORRECTED_DATA column for source {}'.format(targetname)
                step['dependency'] = n - 1
                step['id'] = 'WSDBL'+code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                syscall = ''
                prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                imcall = gen.generate_syscall_wsclean(mslist = [myms],
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
                        threshold = cfg.WSC_THRESHOLD_BLIND,
                        tukeytaper=tukeytaper,
                        minuvl = minuvl,
                        maxuvl = maxuvl,
                        absmem = absmem)
                for call in imcall: 
                    syscall += prefix + call + '\n\n'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

                step = {}
                step['step'] = n
                step['comment'] = 'Make cleaning mask for ' + targetname
                step['dependency'] = n - 1
                step['id'] = 'MASK0'+code
                syscall  = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += gen.generate_syscall_breizorro(restoredimage = f"{img_prefix}-MFS-image.fits", outfile = f"{img_prefix}-MFS-image.mask.fits")[0]
                step['syscall'] = syscall
                steps.append(step)
                n += 1
                mask = f"{img_prefix}-MFS-image.mask.fits"
                print(gen.col('Mask')+ 'None')

            else:
                print(gen.col('Mask')+mask)

            step = {}
            step['step'] = n
            if PARANGMODEL:
                step['comment'] = f'Parang-correct then run wsclean, masked deconvolution for source {targetname} (sky-frame model; DATA itself stays feed-frame)'
            else:
                step['comment'] = 'Run wsclean, masked deconvolution of the DATA column for source {}'.format(targetname)
            step['dependency'] = n - 1
            step['id'] = 'WSDMA'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            syscall = ''
            if PARANGMODEL:
                # DATA is still feed-frame here; parang-correct into CORRECTED_DATA
                # so the datamask model is built in the sky frame (avoids smearing
                # polarised flux when this image is time-averaged), while leaving
                # DATA itself untouched for the feed-frame QuartiCal solve below.
                datamask_datacol = 'CORRECTED_DATA'
                prefix_py = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += prefix_py + f'python3 {TOOLS}/casa_correct_parang.py {myms}\n\n'
            else:
                datamask_datacol = 'DATA'
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                    imgname = data_img_prefix,
                    mfweight = False,
                    datacol = datamask_datacol,
                    mask = mask,
                    chanout = cfg.WSC_DMASK_CHANNELSOUT,
                    # qu_autoscale = False,
                    #intervalsout = False,
                    tukeytaper=tukeytaper,
                    minuvl = minuvl,
                    maxuvl = maxuvl,
                    automask = cfg.WSC_SHALLOWMASK,
                    localrms = cfg.WSC_SHALLOWMASK_LOCALRMS,
                    autothreshold = cfg.WSC_INTER_AUTOTHRESHOLD,
                    nomodel = True,
                    sourcelist = False,
                    absmem = absmem)
            for call in imcall: 
                syscall += prefix + call + '\n\n'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            if cfg.WSC_MAX_CHANNELS < cfg.WSC_DMASK_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                step = {}
                step['step'] = n
                step['comment'] = f'Fix image naming for DATA images'
                step['dependency'] = n - 1
                step['id'] = 'FXDMA' +code
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_DMASK_CHANNELSOUT} {data_img_prefix}\n\n'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Run wsclean-predict on masked deconvolution of the model for source {}'.format(targetname)
            step['dependency'] = n - 1 
            step['id'] = 'PRDMA'+code
            step['slurm_config'] = cfg.SLURM_PREDICT
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + data_img_prefix + '\n\n'
            if cfg.MOD_MODEL_SELFCAL and cfg.WSC_POL != 'I':
                spatial_arg  = target_mod_region if target_mod_region is not None else mask
                syscall += prefix + (f'python3 {TOOLS}/mod_model_selfcal.py '
                                     f'--identifier {data_img_prefix} --stokes V '
                                     f'--spatial {spatial_arg}\n\n')  # --zero-all-neg-I
            syscall += prefix + gen.generate_syscall_predict(msname = myms,
                    imgname = data_img_prefix,
                    chanout = cfg.WSC_DMASK_CHANNELSOUT,
                    absmem = absmem)
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Run Quartical phase self-calibration (stage 1) on the target {}'.format(targetname)
            step['dependency'] = n - 1
            step['id'] = 'C02G2'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
            extra_args = f'output.gain_directory={gain_outdir_2GC} output.log_directory={log_outdir_2GC}'
            if PARANGMODEL:
                # DATA is feed-frame; MODEL_DATA is sky-frame (built from the
                # parang-corrected WSDMA image) — forward-rotate the model to
                # match. Do not derotate the output here; this isn't the final solve.
                extra_args += ' input_model.apply_p_jones=true'
            if ref_ant_arg is not None:
                extra_args += f' solver.reference_antenna={ref_ant_arg}'
            for term, _, new_fi in freq_int_overrides_stage1:
                extra_args += f' {term}.freq_interval={new_fi}'
            if maxuvl != '' or minuvl != '':
                minuv_val = minuvl if minuvl != '' else '0'
                maxuv_val = maxuvl if maxuvl != '' else '0'
                extra_args += f' input_ms.select_uv_range=[{minuv_val},{maxuv_val}]'
            syscall += gen.generate_syscall_quartical(yaml = cfg.CAL_2GC_YAML,
                    myms = myms,
                    extra_args = extra_args)
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            # Split CORRECTED_DATA into stage 2 MS for amplitude self-calibration
            step = {}
            step['step'] = n
            step['comment'] = 'Split CORRECTED_DATA into stage 2 MS for {}'.format(targetname)
            step['dependency'] = n - 1
            step['id'] = 'SPLCD'+code
            step['slurm_config'] = cfg.SLURM_DEFAULTS
            step['pbs_config'] = cfg.PBS_DEFAULTS
            syscall = CONTAINER_RUNNER + CASA_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_casa(casascript=f'{cfg.TOOLS}/casa_split_corrected_ms.py {myms} {stage2_ms}')
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            if PARANGMODEL:
                step['comment'] = f'Parang-correct then run wsclean, masked deconvolution of the stage-1 self-calibrated data for {targetname} (sky-frame model; stage2_ms DATA stays feed-frame)'
            else:
                step['comment'] = f'Run wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated, stage 1) for {targetname}'
            step['dependency'] = n - 1
            step['id'] = 'WSCMI'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            syscall = ''
            if PARANGMODEL:
                # stage2_ms DATA (split from the stage-1 solve) is still feed-frame;
                # parang-correct into CORRECTED_DATA so the inter model is built in
                # the sky frame, leaving stage2_ms DATA untouched for the stage-2
                # amplitude solve below.
                inter_datacol = 'CORRECTED_DATA'
                prefix_py = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += prefix_py + f'python3 {TOOLS}/casa_correct_parang.py {stage2_ms}\n\n'
            else:
                inter_datacol = 'DATA'
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [stage2_ms],
                    imgname = inter_img_prefix,
                    datacol = inter_datacol,
                    mask = mask,
                    chanout = cfg.WSC_DMASK_CHANNELSOUT,
                    tukeytaper=tukeytaper,
                    minuvl = minuvl,
                    maxuvl = maxuvl,
                    # qu_automask_scale = 1.0,
                    localrms = cfg.WSC_INTER_LOCALRMS,
                    automask = cfg.WSC_INTER_AUTOMASK,
                    autothreshold = cfg.WSC_INTER_AUTOTHRESHOLD,
                    nomodel=True,
                    sourcelist = False,
                    absmem = absmem)
            for call in imcall: 
                syscall += prefix + call + '\n\n'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            if cfg.WSC_MAX_CHANNELS < cfg.WSC_DMASK_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                step = {}
                step['step'] = n
                step['comment'] = f'Fix image naming and homogenize beams for INTER images'
                step['dependency'] = n - 1
                step['id'] = 'HOCMI' + code
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_DMASK_CHANNELSOUT} {inter_img_prefix}\n\n'
                syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {inter_img_prefix}'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            if DO_DD_SELFCAL:

                # Path to the QuartiCal YAML for direction-dependent selfcal (used when DO_DD_SELFCAL=True).
                # Should define G+dE (or equivalent DD) terms — distinct from CAL_2GC_YAML_COMPLEX.
                CAL_DDECAL_YAML = DATA + '/quartical/2GC_complex_2dir.yaml'

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
                step['comment'] = f'DD prep for {targetname}: fix NaN, extract DIR, add column, predict DIR, copy, re-predict full model'
                step['dependency'] = n - 1
                step['id'] = 'DDPRE'+code
                step['slurm_config'] = cfg.SLURM_PREDICT
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step, INFRASTRUCTURE, cfg.WSC_ABSMEM)
                py_pfx  = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                wsc_pfx = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' '  if USE_SINGULARITY else ''
                syscall  = py_pfx + 'python3 ' + TOOLS + '/fix_nan_models.py ' + inter_img_prefix + '\n\n'
                if cfg.MOD_MODEL_SELFCAL and cfg.WSC_POL != 'I':
                    spatial_arg = target_mod_region if target_mod_region is not None else mask
                    syscall += py_pfx + (f'python3 {TOOLS}/mod_model_selfcal.py '
                                         f'--identifier {inter_img_prefix} --stokes V '
                                         f'--spatial {spatial_arg}\n\n')
                syscall += py_pfx + 'python3 ' + OXKAT + '/3GC_split_model_images.py '
                syscall +=          '--region ' + dir_region + ' --prefix ' + inter_img_prefix + '\n\n'
                syscall += py_pfx + 'python3 ' + TOOLS + '/add_MS_column.py '
                syscall +=          '--colname DIR1_DATA ' + stage2_ms + '\n\n'
                syscall += wsc_pfx + gen.generate_syscall_predict(msname = stage2_ms,
                        imgname = dir_img_prefix,
                        chanout = cfg.WSC_DMASK_CHANNELSOUT,
                        absmem  = absmem) + '\n\n'
                syscall += py_pfx + 'python3 ' + TOOLS + '/copy_MS_column.py '
                syscall +=          '--fromcol MODEL_DATA --tocol DIR1_DATA ' + stage2_ms + '\n\n'
                syscall += wsc_pfx + gen.generate_syscall_predict(msname = stage2_ms,
                        imgname = inter_img_prefix,
                        chanout = cfg.WSC_DMASK_CHANNELSOUT,
                        absmem  = absmem)
                step['syscall'] = syscall
                steps.append(step)
                n += 1

                step = {}
                step['step'] = n
                step['comment'] = 'Run QuartiCal DD self-calibration (stage 2) on {}'.format(targetname)
                step['dependency'] = n - 1
                step['id'] = 'CL2GC'+code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
                extra_args = (f'output.gain_directory={gain_outdir_stage2} output.log_directory={log_outdir_stage2}'
                              f' input_model.recipe={dd_recipe} output.subtract_directions={dd_subtract}')
                if PARANGMODEL:
                    # Forward-rotate the model to match feed-frame DATA.
                    extra_args += ' input_model.apply_p_jones=true'
                if not cfg.CAL_1GC_APPLYPARANG:
                    # Last solve: derotate CORRECTED_DATA back to the sky frame.
                    extra_args += ' output.apply_p_jones_inv=true'
                if ref_ant_arg is not None:
                    extra_args += f' solver.reference_antenna={ref_ant_arg}'
                for term, _, new_fi in freq_int_overrides_stage2:
                    extra_args += f' {term}.freq_interval={new_fi}'
                if maxuvl != '' or minuvl != '':
                    extra_args += f' input_ms.select_uv_range=[{minuv_val},{maxuv_val}]'
                syscall += gen.generate_syscall_quartical(yaml = CAL_DDECAL_YAML,
                        myms = stage2_ms,
                        extra_args = extra_args)
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            else:

                step = {}
                step['step'] = n
                step['comment'] = 'Run wsclean-predict on intermask model for source {}'.format(targetname)
                step['dependency'] = n - 1
                step['id'] = 'PRCMI'+code
                step['slurm_config'] = cfg.SLURM_PREDICT
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + inter_img_prefix + '\n\n'
                if cfg.MOD_MODEL_SELFCAL and cfg.WSC_POL != 'I':
                    spatial_arg = target_mod_region if target_mod_region is not None else mask
                    syscall += prefix + (f'python3 {TOOLS}/mod_model_selfcal.py '
                                         f'--identifier {inter_img_prefix} --stokes V '
                                         f'--spatial {spatial_arg}\n\n')
                syscall += prefix + gen.generate_syscall_predict(msname = stage2_ms,
                        imgname = inter_img_prefix,
                        chanout = cfg.WSC_DMASK_CHANNELSOUT,
                        absmem = absmem)
                step['syscall'] = syscall
                steps.append(step)
                n += 1

                step = {}
                step['step'] = n
                step['comment'] = 'Run Quartical refined self-calibration (stage 2) on the target {}'.format(targetname)
                step['dependency'] = n - 1
                step['id'] = 'CL2GC'+code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
                extra_args = f'output.gain_directory={gain_outdir_stage2} output.log_directory={log_outdir_stage2}'
                if ref_ant_arg is not None:
                    extra_args += f' solver.reference_antenna={ref_ant_arg}'
                for term, _, new_fi in freq_int_overrides_stage2:
                    extra_args += f' {term}.freq_interval={new_fi}'
                if PARANGMODEL:
                    # Forward-rotate the model to match feed-frame DATA.
                    extra_args += ' input_model.apply_p_jones=true'
                if not cfg.CAL_1GC_APPLYPARANG:
                    # Last solve: derotate CORRECTED_DATA back to the sky frame.
                    extra_args += ' output.apply_p_jones_inv=true'
                if maxuvl != '' or minuvl != '':
                    extra_args += f' input_ms.select_uv_range=[{minuv_val},{maxuv_val}]'
                syscall += gen.generate_syscall_quartical(yaml = cfg.CAL_2GC_YAML_COMPLEX,
                        myms = stage2_ms,
                        extra_args = extra_args)
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            # Safety check: verify parang solve succeeded; fall back to feed-frame + CASA applycal if not
            if not cfg.CAL_1GC_APPLYPARANG:
                fallback_yaml = CAL_DDECAL_YAML if DO_DD_SELFCAL else cfg.CAL_2GC_YAML_COMPLEX
                fallback_extra_args = f'output.gain_directory={gain_outdir_stage2} output.log_directory={log_outdir_stage2}'
                if DO_DD_SELFCAL:
                    fallback_extra_args += f' input_model.recipe={dd_recipe} output.subtract_directions={dd_subtract}'
                if maxuvl != '' or minuvl != '':
                    fallback_extra_args += f' input_ms.select_uv_range=[{minuv_val},{maxuv_val}]'
                fallback_qc_cmd = (CONTAINER_RUNNER + QUARTICAL_CONTAINER + ' ' if USE_SINGULARITY else '') + \
                    gen.generate_syscall_quartical(yaml=fallback_yaml, myms=stage2_ms, extra_args=fallback_extra_args)

                step = {}
                step['step'] = n
                step['comment'] = f'Check parang selfcal succeeded for {targetname}; fall back to feed-frame + CASA applycal if not'
                step['dependency'] = n - 1
                step['id'] = 'CHKPJ'+code
                step['slurm_config'] = cfg.SLURM_DEFAULTS
                step['pbs_config'] = cfg.PBS_DEFAULTS
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall = prefix + (
                    f'python3 {TOOLS}/check_and_fix_parang_selfcal.py '
                    f'{log_outdir_stage2} {stage2_ms} "{fallback_qc_cmd}"'
                    + (' --parangmodel' if PARANGMODEL else '')
                )
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = f'Run wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated, working) for {targetname}'
            step['dependency'] = n - 1
            step['id'] = 'WSCMA'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [stage2_ms],
                    imgname = pcal_img_prefix,
                    datacol = 'CORRECTED_DATA',
                    mask = mask,
                    chanout = cfg.WSC_PCAL_CHANNELSOUT,
                    #intervalsout = False,
                    nomodel=True,
                    sourcelist = False,
                    absmem = absmem)
            for call in imcall:
                syscall += prefix + call + '\n\n'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            if cfg.WSC_MAX_CHANNELS < cfg.WSC_PCAL_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                step = {}
                step['step'] = n
                step['comment'] = f'Homogenize the PCAL resolution across frequency channels'
                step['dependency'] = n - 1
                step['id'] = 'HOCMA' + code
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_PCAL_CHANNELSOUT} {pcal_img_prefix}\n\n'
                syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {pcal_img_prefix}'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            if cfg.WSC_UNIFORM_IMAGE:
                step = {}
                step['step'] = n
                step['comment'] = 'Run high angular resolution, wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated) for {}'.format(targetname)
                step['dependency'] = n - 1
                step['id'] = 'WSUNI'+code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                syscall = ''
                prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                imcall = gen.generate_syscall_wsclean(mslist = [stage2_ms],
                    imgname = uniform_img_prefix,
                    datacol = 'CORRECTED_DATA',
                    mask = mask,
                    chanout = cfg.WSC_BLIND_CHANNELSOUT,
                    nomodel=True,
                    intervalsout=False,
                    weight=cfg.WSC_WEIGHT_HIGHRES,
                    mfweight=True,
                    tukeytaper=False,
                    minuvl = '',
                    maxuvl = '',
                    pol='I',
                    sourcelist = False,
                    absmem = absmem)
                for call in imcall: 
                    syscall += prefix + call + '\n\n'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

                if not cfg.SKIP_PB:
                    step = {}
                    step['step'] = n
                    step['comment'] = 'Apply primary beam correction to '+targetname+'(UNIFORM) image'
                    step['dependency'] = n - 1
                    step['id'] = 'PBUNI'+code
                    syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                    syscall += 'python3 '+TOOLS+'/pbcor_katbeam.py --band '+band[0]+' '+uniform_img_prefix+'-MFS-image.fits'
                    step['syscall'] = syscall
                    steps.append(step)
                    n += 1

            if cfg.WSC_POL != 'I':

                step = {}
                step['step'] = n
                step['comment'] = 'Make Polarization Intensity Images for '+targetname
                step['dependency'] = n - 1
                step['id'] = 'MKLPI'+code
                syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += f"python3 {cfg.TOOLS}/make_pol_images.py {cfg.IMAGES} {targetname} True"
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            target_steps.append((steps,kill_file,targetname))


    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    submit_file = 'submit_2GC_jobs.sh'

    f = open(submit_file,'w')
    f.write('#!/usr/bin/env bash\n')
    f.write('export SINGULARITY_BINDPATH='+cfg.BINDPATH+'\n')

    for content in target_steps:  
        steps = content[0]
        kill_file = content[1]
        targetname = content[2]
        id_list = []

        f.write('\n#---------------------------------------\n')
        f.write('# '+targetname)
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
            f.write('\n# Generate kill script for '+targetname+'\n')
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
