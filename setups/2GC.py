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


# Never True when CAL_1GC_APPLYPARANG is True — 1GC already puts CORRECTED_DATA
# in the sky frame, so re-applying parang here would double-correct it.
PARANGMODEL = cfg.CAL_2GC_PARANGMODEL and not cfg.CAL_1GC_APPLYPARANG


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
        
            # Define output parameters for 2GC steo
            gain_outdir_2GC = GAINTABLES+'/2GC_'+str(filename_targetname)+f'_{stamp}.qc/'
            log_outdir_2GC = LOGS+'/2GC_'+str(filename_targetname)+f'_{stamp}.qc/'

            # Image prefixes
            img_prefix = IMAGES+f'/img_{myms}_datablind'
            data_img_prefix = IMAGES+f'/img_{myms}_datamask'
            pcal_img_prefix = IMAGES+f'/img_{myms}_pcalmask'
            uniform_img_prefix = IMAGES+f'/img_{myms}_uniform'

            # Target-specific kill file
            kill_file = SCRIPTS+'/kill_2GC_jobs_'+filename_targetname+'.sh'

            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('Measurement Set')+myms)
            print(gen.col('Code')+code)

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
                step['comment'] = 'Shallow blind wsclean on DATA column for source {}'.format(targetname)
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
                    intervalsout = False,
                    tukeytaper=tukeytaper,
                    minuvl = minuvl,
                    maxuvl = maxuvl,
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
                step['comment'] = f'Homogenize the MASK resolution across frequency channels'
                step['dependency'] = n - 1
                step['id'] = 'HODMA' +code
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_DMASK_CHANNELSOUT} {data_img_prefix}'
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
            syscall += prefix + gen.generate_syscall_predict(msname = myms,
                    imgname = data_img_prefix,
                    absmem = absmem)
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            if not cfg.SKIP_PB:
                step = {}
                step['step'] = n
                step['comment'] = f'Apply primary beam correction to {targetname} (MASK) image'
                step['dependency'] = n - 1
                step['id'] = 'PBDMA'+code
                syscall = ''
                images = []
                for stoke_i in cfg.WSC_POL:
                    stoke = f'-{stoke_i}'
                    if len(cfg.WSC_POL) == 1:
                        stoke = ''
                    if cfg.WSC_INTERVALSOUT:
                        for t in range(cfg.WSC_INTERVALSOUT):
                            images.append(f'{data_img_prefix}-t{t:04d}-MFS{stoke}-image.fits')
                            if cfg.WSC_MAX_CHANNELS < cfg.WSC_DMASK_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                                images.append(f'{data_img_prefix}-t{t:04d}-MFS{stoke}-image.homogenized.fits')
                    else:
                        images.append(f'{data_img_prefix}-MFS{stoke}-image.fits')
                        if cfg.WSC_MAX_CHANNELS < cfg.WSC_DMASK_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                            images.append(f'{data_img_prefix}-MFS{stoke}-image.homogenized.fits')
                for image in images:
                    syscall += CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                    syscall += 'python3 '+TOOLS+'/pbcor_katbeam.py --band '+band[0]+f' {image}\n\n'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Run Quartical self-calibration on the target {}'.format(targetname)
            step['dependency'] = n - 1
            step['id'] = 'CL2GC'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
            extra_args = f'output.gain_directory={gain_outdir_2GC} output.log_directory={log_outdir_2GC}'
            if PARANGMODEL:
                # DATA is feed-frame; MODEL_DATA is sky-frame (built from the
                # parang-corrected WSDMA image) — forward-rotate the model to match.
                extra_args += ' input_model.apply_p_jones=true'
            if not cfg.CAL_1GC_APPLYPARANG:
                # This is the only/final solve: derotate CORRECTED_DATA back to the sky frame.
                extra_args += ' output.apply_p_jones_inv=true'
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

            # Safety check: if parang was handled by QuartiCal, verify the solve
            # succeeded and fall back to a feed-frame solve + CASA applycal if not
            if not cfg.CAL_1GC_APPLYPARANG:
                # Build the fallback QC command (identical but without apply_p_jones_inv)
                fallback_extra_args = f'output.gain_directory={gain_outdir_2GC} output.log_directory={log_outdir_2GC}'
                if maxuvl != '' or minuvl != '':
                    fallback_extra_args += f' input_ms.select_uv_range=[{minuv_val},{maxuv_val}]'
                fallback_qc_cmd = (CONTAINER_RUNNER + QUARTICAL_CONTAINER + ' ' if USE_SINGULARITY else '') + \
                    gen.generate_syscall_quartical(yaml=cfg.CAL_2GC_YAML, myms=myms, extra_args=fallback_extra_args)

                step = {}
                step['step'] = n
                step['comment'] = f'Check parang selfcal succeeded for {targetname}; fall back to feed-frame solve + CASA applycal if not'
                step['dependency'] = n - 1
                step['id'] = 'CHKPJ'+code
                step['slurm_config'] = cfg.SLURM_DEFAULTS
                step['pbs_config'] = cfg.PBS_DEFAULTS
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall = prefix + (
                    f'python3 {TOOLS}/check_and_fix_parang_selfcal.py '
                    f'{log_outdir_2GC} {myms} "{fallback_qc_cmd}"'
                    + (' --parangmodel' if PARANGMODEL else '')
                )
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = f'Run wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated) for {targetname}'
            step['dependency'] = n - 1
            step['id'] = 'WSCMA'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                    imgname = pcal_img_prefix,
                    datacol = 'CORRECTED_DATA',
                    mask = mask,
                    chanout = cfg.WSC_PCAL_CHANNELSOUT,
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

            if not cfg.SKIP_PB:
                step = {}
                step['step'] = n
                step['comment'] = 'Apply primary beam correction to '+targetname+'(PCAL) image'
                step['dependency'] = n - 1
                step['id'] = 'PBPCL'+code
                syscall = ''
                images = []
                for stoke_i in cfg.WSC_POL:
                    stoke = f'-{stoke_i}'
                    if len(cfg.WSC_POL) == 1:
                        stoke = ''
                    if cfg.WSC_INTERVALSOUT:
                        for t in range(cfg.WSC_INTERVALSOUT):
                            images.append(f'{pcal_img_prefix}-t{t:04d}-MFS{stoke}-image.fits')
                            if cfg.WSC_MAX_CHANNELS < cfg.WSC_PCAL_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                                images.append(f'{pcal_img_prefix}-t{t:04d}-MFS{stoke}-image.homogenized.fits')
                    else:
                        images.append(f'{pcal_img_prefix}-MFS{stoke}-image.fits')
                        if cfg.WSC_MAX_CHANNELS < cfg.WSC_PCAL_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                            images.append(f'{pcal_img_prefix}-MFS{stoke}-image.homogenized.fits')
                for image in images:
                    syscall += CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                    syscall += 'python3 '+TOOLS+'/pbcor_katbeam.py --band '+band[0]+f' {image}\n\n'
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
                imcall = gen.generate_syscall_wsclean(mslist = [myms],
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
                
                # Only make Plin images
                only_Plin = True
                if project_info['polang_name'] == '':
                    only_Plin = False
    
                step = {}
                step['step'] = n
                step['comment'] = 'Make Polarization Intensity Images for '+targetname
                step['dependency'] = n - 1
                step['id'] = 'MKLPI'+code
                syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += f"python3 {cfg.TOOLS}/make_pol_images.py {cfg.IMAGES} {targetname} {only_Plin}"
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
