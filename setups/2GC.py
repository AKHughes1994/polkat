#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


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


    # Get target information from project json

    with open('project_info.json') as f:
        project_info = json.load(f)
    
    band = project_info['band']
    target_ids = project_info['target_ids'] 
    target_names = project_info['target_names']
    target_ms = project_info['target_ms']
    myms = project_info['working_ms']
    pcals = project_info['target_cal_map']

    # ------------------------------------------------------------------------------
    #
    # 2GC recipe definition
    #
    # ------------------------------------------------------------------------------


    target_steps = []
    codes = []
    ii = 1

    # Loop over targets

    for tt in range(0,len(target_ids)):

        targetname   = target_names[tt]
        pcal = pcals[tt]
        restore_flag  = 'after_Xf_solutions'

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
        
            # Image prefix
            name_ms = myms.replace('.ms', f'_{targetname}.ms') # This makes the naming convention the same as oxkat
            img_prefix = IMAGES+f'/img_{name_ms}_datablind'
            data_img_prefix = IMAGES+f'/img_{name_ms}_datamask'
            pcal_img_prefix = IMAGES+f'/img_{name_ms}_pcalmask'
            uniform_img_prefix = IMAGES+f'/img_{name_ms}_uniform'

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
                    datacol = 'CORRECTED_DATA',
                    fields = targetindex,
                    strategy = 'polarisation')
            step['syscall'] = syscall
            steps.append(step)
            n += 1

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
                        datacol = 'CORRECTED_DATA',
                        chanout = cfg.WSC_MASK_CHANNELSOUT,
                        nomodel = True,
                        pol = 'I',
                        intervalsout = False,
                        mfweight = True,
                        localrms = False,
                        automask = 10.0,
                        autothreshold = 3.0,
                        tukeytaper=False,
                        field=targetindex,
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

            step = {}
            step['step'] = n
            step['comment'] = 'Apply primary beam correction to '+targetname+' (BLIND) image'
            step['dependency'] = n - 2
            step['id'] = 'PBDBL'+code
            syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += 'python3 '+TOOLS+'/pbcor_katbeam.py --band '+band[0]+' '+img_prefix+'-MFS-image.fits'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Run wsclean, masked deconvolution of the CORRECTED_DATA column for source {}'.format(targetname)
            step['dependency'] = n - 1 
            step['id'] = 'WSDMA'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                    imgname = data_img_prefix,
                    datacol = 'CORRECTED_DATA',
                    mask = img_prefix+'-MFS-image.mask.fits',
                    chanout = cfg.WSC_IMAGE_CHANNELSOUT,
                    field=targetindex,
                    tukeytaper=False,
                    nomodel = True,
                    sourcelist = False,
                    absmem = absmem)
            for call in imcall: 
                syscall += prefix + call + '\n\n'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            if cfg.WSC_MAX_CHANNELS < cfg.WSC_IMAGE_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                step = {}
                step['step'] = n
                step['comment'] = f'Homogenize the MASK resolution across frequency channels'
                step['dependency'] = n - 1
                step['id'] = 'HODMA' +code
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_IMAGE_CHANNELSOUT} {data_img_prefix}\n\n'
                syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {data_img_prefix}'
                step['syscall'] = prefix + syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Run wsclean-predict on masked deconvolution of the model for source {}'.format(targetname)
            step['dependency'] = n - 1 
            step['id'] = 'PRDMA'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + data_img_prefix + '\n\n'
            syscall += prefix + gen.generate_syscall_predict(msname = myms,
                    imgname = data_img_prefix,
                    field=targetindex,
                    absmem = absmem)
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            step['comment'] = f'Apply primary beam correction to {targetname} (MASK) image'
            step['dependency'] = n - 1
            step['id'] = 'PBDMA'+code
            syscall = ''
            images = []
            for stoke in cfg.WSC_POL:
                if cfg.WSC_INTERVALSOUT:
                    for t in range(cfg.WSC_INTERVALSOUT):
                        images.append(f'{data_img_prefix}-t{t:04d}-MFS-{stoke}-image.fits')
                        if cfg.WSC_MAX_CHANNELS < cfg.WSC_IMAGE_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                            images.append(f'{data_img_prefix}-t{t:04d}-MFS-{stoke}-image.homogenized.fits')
                else:
                    images.append(f'{data_img_prefix}-MFS-{stoke}-image.fits')
                    if cfg.WSC_MAX_CHANNELS < cfg.WSC_IMAGE_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                        images.append(f'{data_img_prefix}-MFS-{stoke}-image.homogenized.fits')
            for image in images:
                syscall += CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += 'python3 '+TOOLS+'/pbcor_katbeam.py --band '+band[0]+f' {image}\n\n'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Run phase self-calibration on the target {}'.format(targetname)
            step['dependency'] = n - 1
            step['id'] = 'CL2GC'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            syscall = CONTAINER_RUNNER + CASA_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_casa(casascript=cfg.OXKAT+f'/2GC_casa_pcal_source.py -f {save_flag} -s {cfg.CAL_2GC_PSOLINT} {targetname}')
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            restore_flag = save_flag

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
                    mask = img_prefix+'-MFS-image.mask.fits',
                    chanout = cfg.WSC_IMAGE_CHANNELSOUT,
                    nomodel=True,
                    field=targetindex,
                    sourcelist = False,
                    absmem = absmem)
            for call in imcall: 
                syscall += prefix + call + '\n\n'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            if cfg.WSC_MAX_CHANNELS < cfg.WSC_IMAGE_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                step = {}
                step['step'] = n
                step['comment'] = f'Homogenize the PCAL resolution across frequency channels'
                step['dependency'] = n - 1
                step['id'] = 'HOCMA' + code
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_IMAGE_CHANNELSOUT} {pcal_img_prefix}\n\n'
                syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {pcal_img_prefix}'
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                step['syscall'] = prefix + syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Apply primary beam correction to '+targetname+'(PCAL) image'
            step['dependency'] = n - 1
            step['id'] = 'PBPCL'+code
            syscall = ''
            images = []
            for stoke in cfg.WSC_POL:
                if cfg.WSC_INTERVALSOUT:
                    for t in range(cfg.WSC_INTERVALSOUT):
                        images.append(f'{pcal_img_prefix}-t{t:04d}-MFS-{stoke}-image.fits')
                        if cfg.WSC_MAX_CHANNELS < cfg.WSC_IMAGE_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                            images.append(f'{pcal_img_prefix}-t{t:04d}-MFS-{stoke}-image.homogenized.fits')
                else:
                    images.append(f'{pcal_img_prefix}-MFS-{stoke}-image.fits')
                    if cfg.WSC_MAX_CHANNELS < cfg.WSC_IMAGE_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                        images.append(f'{pcal_img_prefix}-MFS-{stoke}-image.homogenized.fits')
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
                    mask = img_prefix+'-MFS-image.mask.fits',
                    chanout = cfg.WSC_MASK_CHANNELSOUT,
                    nomodel=True,
                    intervalsout=False,
                    field=targetindex,
                    weight=cfg.WSC_WEIGHT_HIGHRES,
                    mfweight=True,
                    pol='I',
                    sourcelist = False,
                    absmem = absmem)
                for call in imcall: 
                    syscall += prefix + call + '\n\n'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

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
                step['comment'] = 'Make Polarization Intensity Images'
                step['dependency'] = n - 1
                step['id'] = 'MKLPI'+code
                syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += f"python3 {cfg.TOOLS}/make_pol_images.py {cfg.IMAGES}"
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
