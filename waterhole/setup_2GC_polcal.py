#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.a.uk

import glob
import json
import os.path as o
import sys
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg


def main():

    USE_SINGULARITY = cfg.USE_SINGULARITY

    gen.preamble()
    print(gen.col()+'2GC Polarization Angle Calibrator (Imaging & DI phase self-calibration) setup')
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
    myms = project_info['working_ms']

    # Determine if the blind/datamask images are going to be tapered
    if not cfg.WSC_TAPERMASK:
        tukeytaper = False
        minuvl = ''
    else:
        tukeytaper = cfg.WSC_TUKEYTAPER
        minuvl = cfg.WSC_MINUVL   

    # ------------------------------------------------------------------------------
    #
    # 2GC Polarization Angle Calibrator recipe definition
    #
    # ------------------------------------------------------------------------------


    stamp = gen.timenow()

    # Get polarization angle calibrator info from project_info
    pacal_name = project_info.get('polang_name', '')
    
    if not pacal_name:
        gen.print_spacer()
        print(gen.col('ERROR')+'No polarization angle calibrator specified in project_info')
        sys.exit()

    # Set up the MS name - will be created by the split step from myms
    polms = myms.replace('.ms', f'_{pacal_name}.ms')
    targetname = pacal_name
    filename_targetname = gen.scrub_target_name(targetname)

    steps = []        
    code = gen.get_target_code(targetname)

    # Define output parameters for 2GC step
    gain_outdir_2GC = GAINTABLES+'/2GC_'+str(filename_targetname)+f'_{stamp}.qc/'
    log_outdir_2GC = LOGS+'/2GC_'+str(filename_targetname)+f'_{stamp}.qc/'

    # Image prefixes - using diagnostic for all imaging steps
    diagnostic_img_prefix = IMAGES+f'/img_{polms}_diagnostic'
    
    # Check if mask already exists from previous run
    existing_mask = IMAGES+f'/img_{polms.split("/")[-1]}_mask-MFS-image.mask.fits'
    mask_exists = o.isfile(existing_mask)

    # Target-specific kill file
    kill_file = SCRIPTS+'/kill_2GC_polcal_jobs.sh'

    gen.print_spacer()
    print(gen.col('Polarization Angle Calibrator')+targetname)
    print(gen.col('Measurement Set')+polms)
    print(gen.col('Code')+code)
    if mask_exists:
        print(gen.col('Existing Mask')+existing_mask)

    n = 0
    step = {}
    step['step'] = n
    step['comment'] = 'Split polarization angle calibrator from main MS'
    step['dependency'] = None
    step['id'] = 'SPPAC'+code
    syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_casa(casascript=cfg.TOOLS+'/casa_split_polcal.py')
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = 'Run Tricolour on polarization angle calibrator'
    step['dependency'] = n - 1
    step['id'] = 'TRPAC'+code
    step['slurm_config'] = cfg.SLURM_TRICOLOUR
    step['pbs_config'] = cfg.PBS_TRICOLOUR
    syscall = CONTAINER_RUNNER+TRICOLOUR_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_tricolour(myms = polms,
            config = DATA+'/tricolour/target_flagging_1_narrow.yaml',
            datacol = 'DATA',
            strategy = 'polarisation')
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    
    # Check for existing mask file, if doesn't exist make one
    # Define mask_img_prefix for potential re-masking later
    mask_img_prefix = f"{cfg.IMAGES}/img_{polms}_mask"
    
    if mask_exists:
        mask = existing_mask
        print(gen.col('Initially using existing mask')+mask)
    else:
        # Create new mask
        step = {}
        step['step'] = n
        step['comment'] = 'Shallow blind wsclean on DATA column for polarization angle calibrator'
        step['dependency'] = n - 1
        step['id'] = 'WSMSK'+code
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        syscall = ''
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        imcall = gen.generate_syscall_wsclean(mslist = [polms],
                imgname = mask_img_prefix,
                datacol = 'DATA',
                chanout = cfg.WSC_BLIND_CHANNELSOUT,
                imsize = cfg.WSC_CAL_IMSIZE,
                nomodel = True,
                pol = 'I',
                intervalsout = False,
                mfweight = True,
                localrms = False,
                automask = 10.0,
                autothreshold = 3.0,
                tukeytaper=tukeytaper,
                minuvl = minuvl,
                absmem = absmem)
        for call in imcall: 
            syscall += prefix + call + '\n\n'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Make cleaning mask for polarization angle calibrator'
        step['dependency'] = n - 1
        step['id'] = 'MKMSK'+code
        syscall  = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += gen.generate_syscall_breizorro(restoredimage = f"{mask_img_prefix}-MFS-image.fits", outfile = f"{mask_img_prefix}-MFS-image.mask.fits")[0]
        step['syscall'] = syscall
        steps.append(step)
        n += 1
        mask = f"{mask_img_prefix}-MFS-image.mask.fits"
        print(gen.col('Created new mask')+mask)

    step = {}
    step['step'] = n
    step['comment'] = 'Run wsclean, masked deconvolution of DATA column for polarization angle calibrator'
    step['dependency'] = n - 1 
    step['id'] = 'WSD01'+code
    step['slurm_config'] = cfg.SLURM_WSCLEAN
    step['pbs_config'] = cfg.PBS_WSCLEAN
    absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
    syscall = ''
    prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
    imcall = gen.generate_syscall_wsclean(mslist = [polms],
            imgname = diagnostic_img_prefix,
            datacol = 'DATA',
            mask = mask,
            weight = cfg.WSC_WEIGHT_CAL,
            chanout = cfg.WSC_CAL_CHANNELSOUT,
            imsize = cfg.WSC_CAL_IMSIZE,
            tukeytaper = tukeytaper,
            minuvl = minuvl,
            automask = cfg.WSC_SHALLOWMASK,
            joinpolarizations = True,
            multiscale = False,
            splitpol = False,
            mfweight = False,
            intervalsout = False,
            localrms = False,
            threshold = False,
            nomodel = True,
            sourcelist = False,
            absmem = absmem)
    for call in imcall: 
        syscall += prefix + call + '\n\n'
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    if cfg.WSC_MAX_CHANNELS < cfg.WSC_CAL_CHANNELSOUT:
        step = {}
        step['step'] = n
        step['comment'] = f'Fix naming for diagnostic images'
        step['dependency'] = n - 1
        step['id'] = 'FXD01' +code
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_CAL_CHANNELSOUT} {diagnostic_img_prefix}'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

    step = {}
    step['step'] = n
    step['comment'] = 'Run wsclean-predict on diagnostic model for polarization angle calibrator'
    step['dependency'] = n - 1 
    step['id'] = 'PRD01'+code
    step['slurm_config'] = cfg.SLURM_PREDICT
    step['pbs_config'] = cfg.PBS_WSCLEAN
    absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
    prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + diagnostic_img_prefix + '\n\n'
    syscall += prefix + gen.generate_syscall_predict(msname = polms,
            imgname = diagnostic_img_prefix,
            chanout = cfg.WSC_CAL_CHANNELSOUT,
            absmem = absmem)
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = 'Run Quartical self-calibration on polarization angle calibrator (round 1)'
    step['dependency'] = n - 1
    step['id'] = 'QC2G1'+code
    syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_quartical(yaml = cfg.CAL_2GC_YAML,
            myms = polms,
            extra_args = f'output.gain_directory={gain_outdir_2GC} output.log_directory={log_outdir_2GC}')
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = f'Re-image polarization angle calibrator with CORRECTED_DATA (self-calibrated, round 1)'
    step['dependency'] = n - 1
    step['id'] = 'WSD02'+code
    step['slurm_config'] = cfg.SLURM_WSCLEAN
    step['pbs_config'] = cfg.PBS_WSCLEAN
    absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
    syscall = ''
    prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
    imcall = gen.generate_syscall_wsclean(mslist = [polms],
            imgname = mask_img_prefix,
            datacol = 'CORRECTED_DATA',
            mask = mask,
            weight = cfg.WSC_WEIGHT_CAL,
            chanout = cfg.WSC_BLIND_CHANNELSOUT,
            pol='I',
            imsize = cfg.WSC_CAL_IMSIZE,
            tukeytaper = tukeytaper,
            minuvl = minuvl,
            automask = cfg.WSC_AUTOMASK,
            joinpolarizations = True,
            multiscale = False,
            splitpol = False,
            mfweight = False,
            intervalsout = False,
            localrms = False,
            threshold = False,
            nomodel = True,
            sourcelist = False,
            absmem = absmem)
    for call in imcall:
        syscall += prefix + call + '\n\n'
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    # Make new mask after round 1 selfcal
    step = {}
    step['step'] = n
    step['comment'] = 'Make cleaning mask from round 1 self-calibrated image'
    step['dependency'] = n - 1
    step['id'] = 'MKM02'+code
    syscall  = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_breizorro(restoredimage = f"{mask_img_prefix}-MFS-image.fits", outfile = mask, thresh = 8.0)[0]
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = f'Re-image with new mask after round 1 self-calibration'
    step['dependency'] = n - 1
    step['id'] = 'WSD03'+code
    step['slurm_config'] = cfg.SLURM_WSCLEAN
    step['pbs_config'] = cfg.PBS_WSCLEAN
    absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
    syscall = ''
    prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
    imcall = gen.generate_syscall_wsclean(mslist = [polms],
            imgname = diagnostic_img_prefix,
            datacol = 'CORRECTED_DATA',
            mask = mask,
            weight = cfg.WSC_WEIGHT_CAL,
            chanout = cfg.WSC_CAL_CHANNELSOUT,
            imsize = cfg.WSC_CAL_IMSIZE,
            tukeytaper = tukeytaper,
            minuvl = minuvl,
            automask = cfg.WSC_AUTOMASK,
            joinpolarizations = True,
            multiscale = False,
            splitpol = False,
            mfweight = False,
            intervalsout = False,
            localrms = False,
            threshold = False,
            nomodel = True,
            sourcelist = False,
            absmem = absmem)
    for call in imcall:
        syscall += prefix + call + '\n\n'
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    if cfg.WSC_MAX_CHANNELS < cfg.WSC_CAL_CHANNELSOUT:
        step = {}
        step['step'] = n
        step['comment'] = f'Fix naming for diagnostic images'
        step['dependency'] = n - 1
        step['id'] = 'FXD03' + code
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_CAL_CHANNELSOUT} {diagnostic_img_prefix}'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

    step = {}
    step['step'] = n
    step['comment'] = 'Run wsclean-predict on diagnostic model for polarization angle calibrator'
    step['dependency'] = n - 1 
    step['id'] = 'PRD03'+code
    step['slurm_config'] = cfg.SLURM_PREDICT
    step['pbs_config'] = cfg.PBS_WSCLEAN
    absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
    prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + diagnostic_img_prefix + '\n\n'
    syscall += prefix + gen.generate_syscall_predict(msname = polms,
            imgname = diagnostic_img_prefix,
            chanout = cfg.WSC_CAL_CHANNELSOUT,
            absmem = absmem)
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = f'Run Quartical selfcal (round 2) on polarization angle calibrator'
    step['dependency'] = n - 1
    step['id'] = 'QC2G2'+code
    syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_quartical(yaml = cfg.CAL_2GC_YAML,
            myms = polms,
            extra_args = f'output.gain_directory={gain_outdir_2GC} output.log_directory={log_outdir_2GC}')
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = f'Final image of polarization angle calibrator with CORRECTED_DATA (self-calibrated, round 2)'
    step['dependency'] = n - 1
    step['id'] = 'WSD04'+code
    step['slurm_config'] = cfg.SLURM_WSCLEAN
    step['pbs_config'] = cfg.PBS_WSCLEAN
    absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
    syscall = ''
    prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
    imcall = gen.generate_syscall_wsclean(mslist = [polms],
            imgname = diagnostic_img_prefix,
            datacol = 'CORRECTED_DATA',
            mask = mask,
            weight = cfg.WSC_WEIGHT_CAL,
            chanout = cfg.WSC_CAL_CHANNELSOUT,
            imsize = cfg.WSC_CAL_IMSIZE,
            tukeytaper = tukeytaper,
            minuvl = minuvl,
            automask = cfg.WSC_AUTOMASK,
            joinpolarizations = True,
            multiscale = False,
            splitpol = False,
            mfweight = False,
            intervalsout = False,
            localrms = False,
            threshold = False,
            nomodel = True,
            sourcelist = False,
            absmem = absmem)
    for call in imcall:
        syscall += prefix + call + '\n\n'
    step['syscall'] = syscall
    steps.append(step)
    n += 1

    if cfg.WSC_MAX_CHANNELS < cfg.WSC_CAL_CHANNELSOUT:
        step = {}
        step['step'] = n
        step['comment'] = f'Fix naming and homogenize the diagnostic resolution across frequency channels'
        step['dependency'] = n - 1
        step['id'] = 'HODSC' + code
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_CAL_CHANNELSOUT} {diagnostic_img_prefix}\n\n'
        syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {diagnostic_img_prefix}'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

    if cfg.WSC_POL != 'I':
        step = {}
        step['step'] = n
        step['comment'] = 'Make Polarization Intensity Images for polarization angle calibrator'
        step['dependency'] = n - 1
        step['id'] = 'MKLPI'+code
        syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += f"python3 {cfg.TOOLS}/make_pol_images.py {cfg.IMAGES}"
        step['syscall'] = syscall
        steps.append(step)
        n += 1


    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    submit_file = 'submit_2GC_polcal_jobs.sh'

    f = open(submit_file,'w')
    f.write('#!/usr/bin/env bash\n')
    f.write('export SINGULARITY_BINDPATH='+cfg.BINDPATH+'\n')

    id_list = []

    f.write('\n#---------------------------------------\n')
    f.write('# Polarization Angle Calibrator: '+targetname)
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
        f.write('\n# Generate kill script for polarization angle calibrator\n')
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
