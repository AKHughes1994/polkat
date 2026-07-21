#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import json
import os.path as o
import sys
import re
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg


def main():
    # Set to a directory path to use a single shared temp directory, or '' to use per-field temp directories
    OVERRIDE_TEMP_DIR = '/mnt/scratchssd'

    USE_SINGULARITY = cfg.USE_SINGULARITY

    gen.preamble()
    print(gen.col()+'Snapshot (i.e., second-timescale) imaging setup')
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
    IMAGES = cfg.IMAGES
    INTERVALS = cfg.INTERVALS


    gen.setup_dir(IMAGES)
    gen.setup_dir(cfg.LOGS)
    gen.setup_dir(cfg.SCRIPTS)
    gen.setup_dir(cfg.INTERVALS, relabel=True)


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
    
    # Build unified list of fields including calibrators with scan naming
    target_names = []
    valid_field_names = set()  # Set to track all valid field names for validation
    
    # Add working targets
    for name in project_info['working_names']:
        target_names.append(name)
        valid_field_names.add(name)
    
    # Add primary calibrator scans (fieldname_scanXX format)
    bpcal_name = project_info['primary_name']
    primary_ms = project_info['primary_ms']
    for myms_cal in primary_ms:
        scan_str = myms_cal.split('_scan')[-1].replace('.ms', '')
        field_name = f"{bpcal_name}_scan{scan_str}"
        target_names.append(field_name)
        valid_field_names.add(field_name)
    
    # Add polarization angle calibrator scans if present
    pacal_name = project_info.get('polang_name', '')
    pacal_ms = project_info.get('polang_ms', [])
    if pacal_name != '' and len(pacal_ms) > 0:
        for myms_cal in pacal_ms:
            scan_str = myms_cal.split('_scan')[-1].replace('.ms', '')
            field_name = f"{pacal_name}_scan{scan_str}"
            target_names.append(field_name)
            valid_field_names.add(field_name)
    
    # Add secondary calibrator scans
    pcal_names = project_info['secondary_names']
    pcal_ms = project_info['secondary_ms']
    for ss in range(len(pcal_names)):
        for myms_cal in pcal_ms[ss]:
            scan_str = myms_cal.split('_scan')[-1].replace('.ms', '')
            field_name = f"{pcal_names[ss]}_scan{scan_str}"
            target_names.append(field_name)
            valid_field_names.add(field_name)
    
    myms = project_info['working_ms']

    if cfg.SNAP_FIELDS != '':
        target_names = cfg.SNAP_FIELDS.split(',')

    predict_pol = cfg.WSC_POL

    pol = 'I'
    if cfg.SNAP_POL == True:
        pol = 'IQUV'
    
    # ------------------------------------------------------------------------------
    # Calibrator (non-target) imaging parameters
    # These are used for calibrator fields to reduce image size/channels
    # ------------------------------------------------------------------------------
    CALIBRATOR_IMSIZE = 2560
    CALIBRATOR_MAXCHAN = 128
    CALIBRATOR_CHANNELSOUT = 64
    
    # ------------------------------------------------------------------------------
    #
    # SNAP recipe definition
    #
    # ------------------------------------------------------------------------------


    target_steps = []

    # Initialize workflow by sequentially spliting
    steps = []
    kill_file = SCRIPTS+'/kill_snap_split_jobs.sh'    
    n = 0    

    for tt in range(0,len(target_names)):
        
        targetname   = target_names[tt]

        if targetname not in valid_field_names:

            gen.print_spacer()
            print(gen.col('Snap Target')+targetname)
            print('Target not in MS file, skipping')

        else:

            # Logistics
            code = gen.get_target_code(targetname)
            dependency = n - 1
            if n == 0:
                dependency = None

            step = {}
            step['step'] = n
            step['comment'] = 'Splitting out field '+targetname
            step['dependency'] = dependency
            step['id'] = 'SNPTS'+code
            syscall = CONTAINER_RUNNER + CASA_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_casa(casascript=cfg.OXKAT+f'/SNAP_split_sources.py {targetname}')
            step['syscall'] = syscall
            steps.append(step)
            n += 1
            last_split_code = 'SNPTS' + code # Save this as a variable to set dependencies for imaging
        
    target_steps.append((steps,kill_file,'Split MS Files'))

    
    # Loop over targets for imaging -- run imaging concurrently 
    for tt in range(0,len(target_names)):

        targetname   = target_names[tt]

        if targetname not in valid_field_names:

            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('MS')+'not found, skipping')

        else:

            # Logistics
            code = gen.get_target_code(targetname)
            steps = []        
            filename_targetname = gen.scrub_target_name(targetname)

            # Determine if this is a target or calibrator
            # If it has _scanXX in the name, it's a calibrator
            is_target = not re.search(r'_scan\d+$', targetname)

            # Target-specific kill file
            kill_file = SCRIPTS+'/kill_snap_jobs_'+filename_targetname+'.sh'

            # Specify target ms name
            target_ms = myms.replace('.ms', f'_{targetname}_snapshot.ms')

            # Create source-specific subdirectory in IMAGES
            img_dir = f"{IMAGES}/{filename_targetname}"
            gen.setup_dir(img_dir)
            
            # Create temp directory for wsclean intermediate files
            if OVERRIDE_TEMP_DIR == '':
                temp_dir = f"{IMAGES}/{filename_targetname}_temp"
            else:
                temp_dir = f"{OVERRIDE_TEMP_DIR}/{filename_targetname}_temp"
            gen.setup_dir(temp_dir)

            gen.print_spacer()
            field_type = 'Target' if is_target else 'Calibrator'
            print(gen.col(field_type)+targetname)
            print(gen.col('Measurement Set')+myms)
            print(gen.col('Code')+code)

            n = 0 # Start counting

            # If model image(s) have been specified use it to predict [DEFAULT assumes 2GC pcalmask]
            model_image_prefix = img_dir + '/img_' + target_ms.replace('_snapshot','') + '_' + cfg.SNAP_MODELIDENTIFIER
            if cfg.SNAP_MODELIDENTIFIER != '' and glob.glob(model_image_prefix + '*model.fits' ) != [] and cfg.WSC_PCAL_CHANNELSOUT == cfg.SNAP_CHANNELSOUT:
                pass

            # Say (for example) you *accidentally* removed pcalmask model, then this is necessary
            else:

                model_mask = cfg.SNAP_MODELMASK
                model_image_prefix = img_dir+'/img_'+target_ms+'_snapmask'
                blind_image_prefix = img_dir+'/img_'+target_ms+'_snapblind'                

                # If you don't have a mask you need a blind clean as well to make a mask
                if cfg.SNAP_MODELMASK == '' and not o.exists(f"{blind_image_prefix}-MFS-image.mask.fits"):
                    
                    step = {}
                    step['step'] = n
                    step['comment'] = 'Shallow blind wsclean on DATA column of ' + target_ms
                    step['dependency'] = last_split_code
                    step['id'] = 'WBSNA'+code
                    step['slurm_config'] = cfg.SLURM_WSCLEAN
                    step['pbs_config'] = cfg.PBS_WSCLEAN
                    absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                    syscall = ''
                    prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                    # Adjust parameters based on field type
                    imsize_blind = CALIBRATOR_IMSIZE if not is_target else cfg.WSC_IMSIZE
                    maxchan_blind = CALIBRATOR_MAXCHAN if not is_target else cfg.WSC_MAX_CHANNELS
                    imcall = gen.generate_syscall_wsclean(mslist = [target_ms],
                        tempdir = temp_dir,
                        imgname = blind_image_prefix,
                        datacol = 'DATA',
                        chanout = cfg.WSC_BLIND_CHANNELSOUT,
                        nomodel = True,
                        pol = 'I',
                        intervalsout = False,
                        mfweight = True,
                        localrms = False,
                        automask = 10.0,
                        autothreshold = 3.0,
                        tukeytaper=False,
                        field='0',
                        imsize = imsize_blind,
                        maxchan = maxchan_blind,
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
                    step['id'] = 'MKSNA'+code
                    syscall  = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                    syscall += gen.generate_syscall_breizorro(restoredimage = f"{blind_image_prefix}-MFS-image.fits", outfile = f"{blind_image_prefix}-MFS-image.mask.fits")[0]
                    step['syscall'] = syscall
                    steps.append(step)
                    n += 1
                
                    model_mask = f"{blind_image_prefix}-MFS-image.mask.fits"

                # If you do have the mask save it's name
                if o.exists(f"{blind_image_prefix}-MFS-image.mask.fits"):
                    model_mask = f"{blind_image_prefix}-MFS-image.mask.fits"

                if n == 0:
                    dependency = last_split_code
                else:
                    dependency = n - 1 

                step = {}
                step['step'] = n
                step['comment'] = 'Run wsclean, masked deconvolution of the DATA column of ' + target_ms
                step['dependency'] = dependency
                step['id'] = 'WMSNA'+code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                syscall = ''
                prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                # Adjust parameters based on field type
                imsize_snap = CALIBRATOR_IMSIZE if not is_target else cfg.WSC_IMSIZE
                maxchan_snap = CALIBRATOR_MAXCHAN if not is_target else cfg.WSC_MAX_CHANNELS
                chanout_snap = CALIBRATOR_CHANNELSOUT if not is_target else cfg.SNAP_CHANNELSOUT
                pol_snap = 'I' if not is_target else pol
                imcall = gen.generate_syscall_wsclean(mslist = [target_ms],
                    tempdir = temp_dir,
                    imgname = model_image_prefix,
                    datacol = 'DATA',
                    mask = model_mask,
                    chanout = chanout_snap,
                    field= '0',
                    pol = pol_snap,
                    imsize = imsize_snap,
                    maxchan = maxchan_snap,
                    nomodel = True,
                    sourcelist = False,
                    absmem = absmem)
                for call in imcall: 
                    syscall += prefix + call + '\n\n'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

                predict_pol = pol_snap # If you make a new image, predict the polarisations based on the imaged polarisations

                if maxchan_snap < chanout_snap:
                    step = {}
                    step['step'] = n
                    step['comment'] = f'Fix Naming of the snap (MASK) images'
                    step['dependency'] = n - 1
                    step['id'] = 'HOSNA'+code
                    prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                    syscall =  f'python3 {cfg.TOOLS}/fix_image_naming.py {chanout_snap} {model_image_prefix}'
                    step['syscall'] = prefix + syscall 
                    steps.append(step)
                    n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Run wsclean-predict to populate snapshot MS DATA column of ' + target_ms
            if n == 0:
                step['dependency'] = last_split_code
            else:
                step['dependency'] = n - 1 
            step['id'] = 'PRSNA'+code
            step['slurm_config'] = cfg.SLURM_PREDICT
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            # Set chanout_snap if not already defined (when using existing model)
            if 'chanout_snap' not in locals():
                chanout_snap = CALIBRATOR_CHANNELSOUT if not is_target else cfg.SNAP_CHANNELSOUT
            
            syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + model_image_prefix + '\n\n'
            syscall += prefix + gen.generate_syscall_predict(msname = target_ms,
                imgname = model_image_prefix,
                field = '0',
                pol = predict_pol,
                chanout = chanout_snap,
                tempdir = temp_dir,
                absmem = absmem)
            step['syscall'] = syscall
            steps.append(step)
            n += 1
       
       
            # Now that models have been made we can perform the snapshot imaging routine     
            step = {}
            step['step'] = n
            step['comment'] = 'Perform UV-subtraction of the static sky model for ' + target_ms
            step['dependency'] = n - 1
            step['id'] = 'UVSNA'+code
            syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += 'python3 '+cfg.OXKAT+f'/SNAP_uvsub.py {target_ms}'
            step['syscall'] = syscall
            steps.append(step)
            n += 1
            
            step = {}
            step['step'] = n
            step['comment'] = 'Performing per-interval snapshot imaging of the uvsubtracted visbilities for '+targetname
            step['dependency'] = n - 1
            step['id'] = 'IISNA'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            syscall = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += 'python3 '+cfg.OXKAT+f'/SNAP_intervals.py {target_ms} {targetname}'
            step['syscall'] = syscall
            steps.append(step)
            n += 1
        

            step = {}
            step['step'] = n
            step['comment'] = 'Performing static model restoration for '+targetname
            step['dependency'] = n - 1
            step['id'] = 'RESNA'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += 'python3 '+cfg.OXKAT+f'/SNAP_restore.py {model_image_prefix} {target_ms} {targetname}'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            maxchan_snap = CALIBRATOR_MAXCHAN if not is_target else cfg.WSC_MAX_CHANNELS
            chanout_snap = CALIBRATOR_CHANNELSOUT if not is_target else cfg.SNAP_CHANNELSOUT
            if chanout_snap > maxchan_snap:
                # Use target-specific subdirectory under INTERVALS
                intervals_subdir = cfg.INTERVALS+f'/{filename_targetname}'
                restored_image_prefix = intervals_subdir+f'/img_{target_ms}_restored'
                step = {}
                step['step'] = n
                step['comment'] = 'Max channels is less than total channels, making a (really rough) MFS image for '+targetname
                step['dependency'] = n - 1
                step['id'] = 'MFSNA'+code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += 'python3 '+cfg.TOOLS+f'/make_rough_mfs.py {restored_image_prefix}'
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            # More thought needs to go into the usefulness of Homogenizing + SNAPSHOT imaging
            #if cfg.WSC_MAX_CHANNELS < cfg.SNAP_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
            #
            #   snap_img_prefix = cfg.INTERVALS + f'/img_{target_ms}_restored'

            #    step = {}
            #    step['step'] = n
            #    step['comment'] = f'Homogenize the SNAPSHOT imaging beams'
            #    step['dependency'] = n - 1
            #    step['id'] = 'HOSNA' +code
            #    syscall =  f'python3 {cfg.TOOLS}/homogenize_beams.py {snap_img_prefix}'
            #    prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            #    step['syscall'] = prefix + syscall
            #    steps.append(step)
            #    n += 1

            if cfg.SNAP_POL:
                step = {}
                step['step'] = n
                step['comment'] = 'Make Polarization Intensity Images'
                step['dependency'] = n - 1
                step['id'] = 'PISNA'+code
                syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += f"python3 {cfg.TOOLS}/make_pol_images.py {cfg.INTERVALS}"
                step['syscall'] = syscall
                steps.append(step)
                n += 1

            step = {}
            step['step'] = n
            step['comment'] = 'Clean temp directory for '+targetname
            step['dependency'] = n - 1
            step['id'] = 'CLSNA'+code
            syscall = f"rm -rf {temp_dir}"
            step['syscall'] = syscall
            steps.append(step)
            n += 1

            target_steps.append((steps,kill_file,targetname))
            

    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    submit_file = 'submit_snap_jobs.sh'

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
            if type(step['dependency']) == int:
                dependency = steps[step['dependency']]['id']
            elif type(step['dependency']) == str:
                dependency = step['dependency']
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
