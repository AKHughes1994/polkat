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

    # Build unified field list with MS files and flags
    field_list = []
    
    # Add targets (integrated across scans)
    for tt in range(len(target_ids)):
        field_list.append({
            'name': target_names[tt],
            'id': target_ids[tt],
            'ms': target_ms[tt],
            'is_target': True
        })
    
    # Add secondaries (per scan)
    for ss in range(len(pcal_ids)):
        for myms in pcal_ms[ss]:
            scan_str = myms.split('_scan')[-1].replace('.ms', '')
            field_list.append({
                'name': f"{pcal_names[ss]}_scan{scan_str}",
                'id': pcal_ids[ss],
                'ms': myms,
                'is_target': False
            })

    # Determine if the blind/datamask images are going to be tapered
    if not cfg.WSC_TAPERMASK:
        tukeytaper = False
        minuvl = ''
    else:
        tukeytaper = cfg.WSC_TUKEYTAPER
        minuvl = cfg.WSC_MINUVL   

    # ------------------------------------------------------------------------------
    #
    # 2GC recipe definition
    #
    # ------------------------------------------------------------------------------


    field_steps = []
    codes = []
    ii = 1
    stamp = gen.timenow()

    # Loop over all fields (targets and secondaries)
    for field in field_list:

        fieldname = field['name']
        myms = field['ms']
        is_target = field['is_target']

        if is_target and fieldname not in project_info['working_names']:

            gen.print_spacer()
            print(gen.col('Field')+fieldname)
            print(gen.col('MS')+'not found, skipping')

        else:

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

            # Image prefixes
            img_prefix = IMAGES+f'/img_{myms}_datablind'
            data_img_prefix = IMAGES+f'/img_{myms}_datamask'
            pcal_img_prefix = IMAGES+f'/img_{myms}_pcalmask'
            uniform_img_prefix = IMAGES+f'/img_{myms}_uniform'
            notaper_img_prefix = IMAGES+f'/img_{myms}_notaper'

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
                step['comment'] = 'Shallow blind wsclean on DATA column for source {}'.format(fieldname)
                step['dependency'] = n - 1 if is_target else None  # Depends on tricolour for targets only
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
                step['comment'] = 'Make cleaning mask for ' + fieldname
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
            step['comment'] = 'Run wsclean, masked deconvolution of the DATA column for source {}'.format(fieldname)
            step['dependency'] = n - 1 
            step['id'] = 'WSDMA'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                    imgname = data_img_prefix,
                    mfweight = False,
                    datacol = 'DATA',
                    mask = mask,
                    chanout = cfg.WSC_DMASK_CHANNELSOUT,
                    intervalsout = False,
                    tukeytaper=tukeytaper,
                    minuvl = minuvl,
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
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_DMASK_CHANNELSOUT} {data_img_prefix}\n\n'
                syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {data_img_prefix}'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Run wsclean-predict on masked deconvolution of the model for source {}'.format(fieldname)
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

            step = {}
            step['step'] = n
            step['comment'] = 'Run Quartical self-calibration on {}'.format(fieldname)
            step['dependency'] = n - 1
            step['id'] = 'CL2GC'+code
            syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_quartical(yaml = cfg.CAL_2GC_YAML,
                    myms = myms,
                    extra_args = f'output.gain_directory={gain_outdir_2GC} output.log_directory={log_outdir_2GC}')
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            step['comment'] = f'Run wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated) for {fieldname}'
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
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_PCAL_CHANNELSOUT} {pcal_img_prefix}\n\n'
                syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {pcal_img_prefix}'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            # High-resolution and no-taper images only for targets
            if is_target and cfg.WSC_UNIFORM_IMAGE:
                step = {}
                step['step'] = n
                step['comment'] = 'Run high angular resolution, wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated) for {}'.format(fieldname)
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
                    tukeytaper=False,  # Hardcoded no taper
                    minuvl = '',       # Hardcoded no minuvl
                    pol='I',
                    sourcelist = False,
                    absmem = absmem)
                for call in imcall: 
                    syscall += prefix + call + '\n\n'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

                # PB correction for highres image (hardcoded)
                step = {}
                step['step'] = n
                step['comment'] = 'Apply primary beam correction to '+fieldname+'(UNIFORM) image'
                step['dependency'] = n - 1
                step['id'] = 'PBUNI'+code
                syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += 'python3 '+TOOLS+'/pbcor_katbeam.py --band '+band[0]+' '+uniform_img_prefix+'-MFS-image.fits'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            # No-taper image only for targets
            if is_target:
                step = {}
                step['step'] = n
                step['comment'] = 'Run no-taper wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated) for {}'.format(fieldname)
                step['dependency'] = n - 1
                step['id'] = 'WSNOT'+code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                syscall = ''
                prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                imcall = gen.generate_syscall_wsclean(mslist = [myms],
                    imgname = notaper_img_prefix,
                    datacol = 'CORRECTED_DATA',
                    mask = mask,
                    chanout = cfg.WSC_BLIND_CHANNELSOUT,
                    nomodel=True,
                    intervalsout=False,
                    weight=cfg.WSC_WEIGHT,
                    mfweight=True,
                    tukeytaper=False,  # Hardcoded no taper
                    minuvl = '',       # Hardcoded no minuvl
                    pol=cfg.WSC_POL,   # Use default polarization
                    sourcelist = False,
                    absmem = absmem)
                for call in imcall: 
                    syscall += prefix + call + '\n\n'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

                # PB correction for notaper image (hardcoded)
                step = {}
                step['step'] = n
                step['comment'] = 'Apply primary beam correction to '+fieldname+'(NOTAPER) image'
                step['dependency'] = n - 1
                step['id'] = 'PBNOT'+code
                syscall = ''
                for stoke_i in cfg.WSC_POL:
                    stoke = f'-{stoke_i}'
                    if len(cfg.WSC_POL) == 1:
                        stoke = ''
                    syscall += CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                    syscall += 'python3 '+TOOLS+'/pbcor_katbeam.py --band '+band[0]+' '+notaper_img_prefix+f'-MFS{stoke}-image.fits\n\n'
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

            field_type_str = 'target' if is_target else 'secondary'
            field_steps.append((steps,kill_file,fieldname,field_type_str))


    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    submit_file = 'submit_2GC_jobs.sh'

    f = open(submit_file,'w')
    f.write('#!/usr/bin/env bash\n')
    f.write('export SINGULARITY_BINDPATH='+cfg.BINDPATH+'\n')

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
