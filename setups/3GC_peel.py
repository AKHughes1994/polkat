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
    print(gen.col()+'3GC (peeling) setup')
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
    else:
        tukeytaper = cfg.WSC_TUKEYTAPER
        minuvl = cfg.WSC_MINUVL   


    # ------------------------------------------------------------------------------
    #
    # 3GC peeling recipe definition
    #
    # ------------------------------------------------------------------------------


    target_steps = []
    codes = []
    ii = 1
    stamp = gen.timenow()

    # Loop over targets

    for tt in range(0,len(target_ids)):

        targetname = target_names[tt]
        myms = target_ms[tt]
        CAL_3GC_PEEL_REGION = cfg.CAL_3GC_PEEL_REGION
        skip = False

        if not o.isdir(myms):
            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('MS')+'not found, skipping')
            skip = True

        if CAL_3GC_PEEL_REGION == '':
            CAL_3GC_PEEL_REGION = sorted(glob.glob('*'+targetname+'*peel*.reg'))

        if CAL_3GC_PEEL_REGION == []:
            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('Measurement Set')+myms)
            print(gen.col()+'Please provide a DS9 region file definining the source you wish to peel.')
            print(gen.col()+'This can be specified in the config or by placing a file of the form:')
            print(gen.col()+'       *'+targetname+'*peel*.reg')
            print(gen.col()+'in this folder. Skipping.')
            skip = True

        # Make array if it is a string of regions
        if isinstance(CAL_3GC_PEEL_REGION, str):
            CAL_3GC_PEEL_REGION = CAL_3GC_PEEL_REGION.split(',')

        if not skip:

            steps = []        
            filename_targetname = gen.scrub_target_name(targetname)


            code = gen.get_target_code(targetname)
            if code in codes:
                code += '_'+str(ii)
                ii += 1
            codes.append(code)
        
            # Define output parameters for 3GC step
            gain_outdir_3GC = GAINTABLES+'/3GC_peel_'+str(filename_targetname)+f'_{stamp}.qc/'
            log_outdir_3GC = LOGS+'/3GC_peel_'+str(filename_targetname)+f'_{stamp}.qc/'

            # Define recipt and peel directions
            subtract_dirs = str([(k + 1) for k in range(len(CAL_3GC_PEEL_REGION))]).replace(' ','')
            recipe = ':'.join(['~'.join(['MODEL_DATA'] + [f'DIR{k+1}_DATA' for k in range(len(CAL_3GC_PEEL_REGION))])] + [f'DIR{k+1}_DATA' for k in range(len(CAL_3GC_PEEL_REGION))])

            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('Measurement Set')+myms)
            print(gen.col('Code')+code)
            print(gen.col('Peeling region(s):'), CAL_3GC_PEEL_REGION)


            # Image prefixes
            prepeel_img_prefix = IMAGES+'/img_'+myms+'_prepeel'
            peel_img_prefix = IMAGES+'/img_'+myms+'_peel'
            pcal_img_prefix = IMAGES+f'/img_{myms}_pcalmask'
            dir_img_prefixes = [prepeel_img_prefix+'-'+region.split('/')[-1].split('.')[0] for region in CAL_3GC_PEEL_REGION]


            # Target-specific kill file
            kill_file = SCRIPTS+'/kill_3GC_peel_jobs_'+filename_targetname+'.sh'

            # Check for mask file, if doesn't exist make one from the pcalmask MFS image if not return error and skip source 
            # We enforce that a mask is required for peeling
            n = 0
            mask = cfg.WSC_MASK
            if not mask:

                # Look for pcalmask image:
                img_prefix = sorted(glob.glob(f'{pcal_img_prefix}*MFS-I-image.fits') + glob.glob(f'{pcal_img_prefix}*MFS-image.fits'))

                # Look for pcalmask homogenized image (in the case of max channelization)
                if img_prefix == []:
                    img_prefix = sorted(glob.glob(f'{pcal_img_prefix}*MFS-I-image.homogenized.fits') + glob.glob(f'{pcal_img_prefix}*MFS-image.homogenized.fits'))

                # If not mask or masking capabilities, skip source
                if img_prefix == []:
                    print(gen.col()+targetname+' does not have manual mask specified nor 2GC image to make mask')
                    print(gen.col()+'a mask is required for peeling, specify mask or the source will be skipped')
                    continue
                else:
                    img_prefix = img_prefix[0]

                step = {}
                step['step'] = n
                step['comment'] = 'Make cleaning mask for ' + targetname
                step['dependency'] = None
                step['id'] = 'MASK1'+code
                syscall  = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += gen.generate_syscall_breizorro(restoredimage = img_prefix, outfile = img_prefix.replace('.fits','.mask.peel.fits'), thresh = 5.0)[0]
                step['syscall'] = syscall
                steps.append(step)
                n += 1
                mask = img_prefix.replace('.fits','.mask.peel.fits')

            print(gen.col('Mask')+mask)


            step = {}
            step['step'] = n
            step['comment'] = 'Run masked wsclean, high freq/angular resolution, on CORRECTED_DATA column of '+myms
            step['dependency'] = None if n == 0 else n - 1
            step['id'] = 'WSDMA'+code
            step['slurm_config'] = cfg.SLURM_EXTRALONG
            step['pbs_config'] = cfg.PBS_EXTRALONG
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                    imgname = prepeel_img_prefix,
                    datacol = 'CORRECTED_DATA',
                    mask = mask,
                    weight = cfg.CAL_3GC_PEEL_BRIGGS,
                    chanout = cfg.CAL_3GC_PEEL_NCHAN,
                    maxchan = cfg.CAL_3GC_PEEL_MAXCHAN,
                    pol = cfg.CAL_3GC_PEEL_POL,
                    intervalsout = False,
                    tukeytaper = False,
                    automask = 5.0,
                    autothreshold = 1.0,
                    minuvl = '',
                    nomodel = True,
                    sourcelist = False,
                    mfweight=False,
                    absmem = absmem)
            for call in imcall: 
                syscall += prefix + call + '\n\n'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

    
            if cfg.CAL_3GC_PEEL_MAXCHAN < cfg.CAL_3GC_PEEL_NCHAN:
                step = {}
                step['step'] = n
                step['comment'] = f'Fix naming for prepeel image'
                step['dependency'] = n - 1
                step['id'] = 'HODMA' +code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.CAL_3GC_PEEL_NCHAN} {prepeel_img_prefix}'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Fix any NaN values in blanked wsclean sub-band models'
            step['dependency'] = n - 1
            step['id'] = 'FXNAN'+code
            syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += 'python3 '+TOOLS+'/fix_nan_models.py '+prepeel_img_prefix
            step['syscall'] = syscall
            steps.append(step)
            n+= 1

            step = {}
            step['step'] = n
            step['comment'] = 'Extract problem source(s) defined by region(s) into a separate set of model images'
            step['dependency'] = n - 1
            step['id'] = 'IMSPL'+code
            prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall =''
            for k, region in enumerate(CAL_3GC_PEEL_REGION):
                syscall += prefix + 'python3 '+OXKAT+'/3GC_split_model_images.py '
                syscall += '--region '+region+' '
                syscall += '--prefix '+prepeel_img_prefix+' '
                syscall += '\n'
            step['syscall'] = syscall
            steps.append(step)
            n+= 1

            step = {}
            step['step'] = n
            step['comment'] = 'Adding additional model columns to ' + myms
            step['dependency'] = n - 1
            step['id'] = 'ADDIR'+code
            prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall = ''
            for k, region in enumerate(CAL_3GC_PEEL_REGION):
                syscall += prefix + 'python3 '+TOOLS+'/add_MS_column.py '
                syscall += f'--colname DIR{k + 1}_DATA '
                syscall += myms
                syscall += '\n'
            step['syscall'] = syscall
            steps.append(step)
            n+= 1


            # Populate new data columns with problematic direction models
            z = 0
            for k, prefix in enumerate(dir_img_prefixes):

                step = {}
                step['step'] = n + z
                step['comment'] = f'Predict problem source DIR{k+1} visibilities into MODEL_DATA column of '+myms
                step['dependency'] = n + z - 1
                step['id'] = f'WSPR{k+1}'+code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                syscall = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += gen.generate_syscall_predict(msname = myms, imgname = prefix, pol = cfg.CAL_3GC_PEEL_POL, chanout = cfg.CAL_3GC_PEEL_NCHAN, absmem = absmem)
                step['syscall'] = syscall
                steps.append(step)

                step = {}
                step['step'] = n + z + 1
                step['comment'] = f'Copy MODEL_DATA to DIR{k+1}_DATA'
                step['dependency'] = n + z
                step['id'] = f'CPMO{k+1}'+code
                syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += 'python3 '+TOOLS+'/copy_MS_column.py '
                syscall += '--fromcol MODEL_DATA '
                syscall += f'--tocol DIR{k + 1}_DATA '
                syscall += myms
                step['syscall'] = syscall
                steps.append(step)
                z += 2
            n += z

            # Repopulate MODEL_DATA with 
            step = {}
            step['step'] = n
            step['comment'] = 'Re-predict full sky model visibilities into MODEL_DATA column of '+myms
            step['dependency'] = n - 1
            step['id'] = 'WSPRF'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            syscall = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_predict(msname = myms, imgname = prepeel_img_prefix, pol = cfg.CAL_3GC_PEEL_POL, chanout = cfg.CAL_3GC_PEEL_NCHAN, absmem = absmem)
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Copy CORRECTED_DATA to DATA'
            step['dependency'] = n - 1
            step['id'] = 'CPCOR'+code
            syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += 'python3 '+TOOLS+'/copy_MS_column.py '
            syscall += '--fromcol CORRECTED_DATA '
            syscall += '--tocol DATA '
            syscall += myms
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Run QuartiCal to solve for G (full model) and dE (problem source), peel out problem source'
            step['dependency'] = n - 1
            step['id'] = 'CL3GC'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            syscall = CONTAINER_RUNNER+QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_quartical(yaml = cfg.CAL_3GC_PEEL_YAML,
                    myms = myms,
                    extra_args = f'output.gain_directory={gain_outdir_3GC} output.log_directory={log_outdir_3GC} input_model.recipe={recipe} output.subtract_directions={subtract_dirs}')
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            step = {}
            step['step'] = n
            step['comment'] = f'Run wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated/peeled) for {targetname}'
            step['dependency'] = n - 1
            step['id'] = 'WSCPE'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                    imgname = peel_img_prefix,
                    datacol = 'CORRECTED_DATA',
                    mask = mask,
                    chanout = cfg.WSC_IMAGE_CHANNELSOUT,
                    nomodel=True,
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
                step['id'] = 'HOCPE' + code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_IMAGE_CHANNELSOUT} {peel_img_prefix}\n\n'
                syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {peel_img_prefix}'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Apply primary beam correction to '+targetname+'(PEEL) image'
            step['dependency'] = n - 1
            step['id'] = 'PBPPE'+code
            syscall = ''
            images = []
            for stoke_i in cfg.WSC_POL:
                stoke = f'-{stoke_i}'
                if len(cfg.WSC_POL) == 1:
                    stoke = ''
                if cfg.WSC_INTERVALSOUT:
                    for t in range(cfg.WSC_INTERVALSOUT):
                        images.append(f'{peel_img_prefix}-t{t:04d}-MFS{stoke}-image.fits')
                        if cfg.WSC_MAX_CHANNELS < cfg.WSC_IMAGE_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                            images.append(f'{peel_img_prefix}-t{t:04d}-MFS{stoke}-image.homogenized.fits')
                else:
                    images.append(f'{peel_img_prefix}-MFS{stoke}-image.fits')
                    if cfg.WSC_MAX_CHANNELS < cfg.WSC_IMAGE_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
                        images.append(f'{peel_img_prefix}-MFS{stoke}-image.homogenized.fits')
            for image in images:
                syscall += CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += 'python3 '+TOOLS+'/pbcor_katbeam.py --band '+band[0]+f' {image}\n\n'
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


    submit_file = 'submit_3GC_peel_jobs.sh'

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

            step_id = step['id']
            id_list.append(step_id)
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
                            jobname = step_id,
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
