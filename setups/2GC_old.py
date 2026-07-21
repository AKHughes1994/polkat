#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.a.uk

import glob
import json
import os
import os.path as o
import sys
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg


def main():
    # Set this flag to True to enable batch submission (4 at a time), False for no batching
    BATCH_FIELD_SUBMISSION = True

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
    bpcal_name = project_info['primary_name']
    bpcal_id = project_info['primary_id']
    primary_ms = project_info['primary_ms']
    pacal_name = project_info.get('polang_name', '')
    pacal_id = project_info.get('polang_id', '')
    pacal_ms = project_info.get('polang_ms', [])
    
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
            'is_target': True
        })
    
    # Add primary calibrator (per scan)
    # Skip if PRE_FIELDS is active and primary is not in working_names
    if not apply_field_filter or bpcal_name in allowed_fields:
        for myms in primary_ms:
            scan_str = myms.split('_scan')[-1].replace('.ms', '')
            field_list.append({
                'name': f"{bpcal_name}_scan{scan_str}",
                'id': bpcal_id,
                'ms': myms,
                'is_target': False
            })
    
    # Add polarization angle calibrator (per scan) if it exists
    if pacal_name != '' and pacal_id != '' and len(pacal_ms) > 0:
        # Skip if PRE_FIELDS is active and pacal is not in working_names
        if not apply_field_filter or pacal_name in allowed_fields:
            for myms in pacal_ms:
                scan_str = myms.split('_scan')[-1].replace('.ms', '')
                field_list.append({
                    'name': f"{pacal_name}_scan{scan_str}",
                    'id': pacal_id,
                    'ms': myms,
                    'is_target': False
                })
    
    # Add secondaries (per scan)
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
                'is_target': False
            })

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
        temp_dir = f"{IMAGES}/{filename_fieldname}_temp"
        gen.setup_dir(temp_dir)

        # Image prefixes
        img_prefix = img_dir+f'/img_{myms}_datablind'
        data_img_prefix = img_dir+f'/img_{myms}_datamask'
        inter_img_prefix = img_dir+f'/img_{myms}_intermask'
        pcal_img_prefix = img_dir+f'/img_{myms}_pcalmask'
        uniform_img_prefix = img_dir+f'/img_{myms}_uniform'
        notaper_img_prefix = img_dir+f'/img_{myms}_notaper'

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

        
        # Check for mask file, if doesn't exist make one
        mask = cfg.WSC_MASK
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
            cores = step['glam_config']['CPUS']
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
                    weight=weight_blind,
                    imsize = imsize_blind,
                    maxchan = maxchan_blind,
                    cores = cores,
                    absmem = absmem)
            for call in imcall: 
                syscall += prefix + call + '\n\n'
            # Clean up intermediate images (keep only MFS-image.fits files)
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

        step = {}
        step['step'] = n
        step['comment'] = 'Run wsclean, masked deconvolution of the DATA column for source {}'.format(fieldname)
        step['dependency'] = n - 1 
        step['id'] = 'WSDMA'+code
        step['glam_config'] = cfg.GLAM_WSC_SOURCE if is_target else cfg.GLAM_WSC_CAL
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        cores = step['glam_config']['CPUS']
        syscall = ''
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        # Adjust parameters based on field type
        imsize_dmask = CALIBRATOR_IMSIZE if not is_target else cfg.WSC_IMSIZE
        maxchan_dmask = CALIBRATOR_MAXCHAN if not is_target else cfg.WSC_MAX_CHANNELS
        chanout_dmask = CALIBRATOR_CHANNELSOUT if not is_target else cfg.WSC_DMASK_CHANNELSOUT
        # Use 1.5x WSC_SHALLOWMASK for automask for calibrators (non-targets)
        automask_dmask = 1.5 * cfg.WSC_SHALLOWMASK if not is_target else cfg.WSC_SHALLOWMASK
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
                fitspectralpol = 0,
                automask = automask_dmask,
                weight=weight_dmask,
                imsize = imsize_dmask,
                maxchan = maxchan_dmask,
                cores = cores,
                absmem = absmem)
        for call in imcall: 
            syscall += prefix + call + '\n\n'
            # Clean up intermediate images (keep only model.fits and MFS*image.fits files)
            syscall += f'find {img_dir} -type f -name "*{os.path.basename(data_img_prefix)}*" ! -name "*model.fits" ! -name "*MFS*image.fits" -delete\n\n' # DEBUG
            #syscall += f'find {img_dir} -type f -name "*{os.path.basename(data_img_prefix)}*dirty*" -delete\n'
            #syscall += f'find {img_dir} -type f -name "*{os.path.basename(data_img_prefix)}*psf*" -delete\n\n'
        step['syscall'] = syscall
        step['grouping'] = 2
        steps.append(step)
        n += 1

        # Skip homogenization for non-targets
        if is_target and cfg.WSC_MAX_CHANNELS < cfg.WSC_DMASK_CHANNELSOUT:
            step = {}
            step['step'] = n
            step['comment'] = f'Fix image naming for MASK images'
            step['dependency'] = n - 1
            step['id'] = 'FXDMA' +code
            step['glam_config'] = cfg.GLAM_SMALL
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_DMASK_CHANNELSOUT} {data_img_prefix}'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Run wsclean-predict on masked deconvolution of the model for source {}'.format(fieldname)
        step['dependency'] = n - 1 
        step['id'] = 'PRDMA'+code
        step['glam_config'] = cfg.GLAM_MEDIUM
        step['slurm_config'] = cfg.SLURM_PREDICT
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        cores = step['glam_config']['CPUS']
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        # Clean up intermediate images (keep only model.fits and MFS*image.fits files)
        syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + data_img_prefix + '\n\n'
        syscall += prefix + gen.generate_syscall_predict(msname = myms,
                imgname = data_img_prefix,
                cores = cores,
                absmem = absmem,
                tempdir = temp_dir)
        syscall += '\n\nfind ' + img_dir + ' -type f -name "*' + os.path.basename(data_img_prefix) + '*" ! -name "*MFS-image.fits" -delete'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Run Quartical self-calibration (stage 1) on {}'.format(fieldname)
        step['dependency'] = n - 1
        step['id'] = 'C02G2'+code
        step['glam_config'] = cfg.GLAM_COREHEAVY
        cores = step['glam_config']['CPUS']
        syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
        # Build extra args with UV range selection if needed
        extra_args = f'output.gain_directory={gain_outdir_2GC} output.log_directory={log_outdir_2GC}'
        if minuvl != '' or maxuvl != '':
            min_uv = minuvl if minuvl != '' else '0'
            max_uv = maxuvl if maxuvl != '' else '0'
            extra_args += f' input_ms.select_uv_range=[{min_uv},{max_uv}]'
        syscall += gen.generate_syscall_quartical(yaml = cfg.CAL_2GC_YAML,
                myms = myms,
                extra_args = extra_args)
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        step = {}
        step['step'] = n
        step['comment'] = f'Run wsclean, masked deconvolution of the CORRECTED_DATA (stage 1 self-calibrated) for {fieldname}'
        step['dependency'] = n - 1
        step['id'] = 'WSCMI'+code
        step['glam_config'] = cfg.GLAM_WSC_SOURCE if is_target else cfg.GLAM_WSC_CAL
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        cores = step['glam_config']['CPUS']
        syscall = ''
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        # Adjust parameters based on field type
        imsize_inter = CALIBRATOR_IMSIZE if not is_target else cfg.WSC_IMSIZE
        maxchan_inter = CALIBRATOR_MAXCHAN if not is_target else cfg.WSC_MAX_CHANNELS
        chanout_inter = CALIBRATOR_CHANNELSOUT if not is_target else cfg.WSC_DMASK_CHANNELSOUT
        weight_inter = 'uniform' if not is_target else cfg.WSC_WEIGHT
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
                fitspectralpol = 0,
                weight=weight_inter,
                imsize = imsize_inter,
                maxchan = maxchan_inter,
                absmem = absmem)
        for call in imcall: 
            syscall += prefix + call + '\n\n'
            # Clean up intermediate images (keep only model.fits and MFS*image.fits files)
            syscall += f'find {img_dir} -type f -name "*{os.path.basename(inter_img_prefix)}*" ! -name "*model.fits" ! -name "*MFS*image.fits" -delete\n\n'
        step['syscall'] = syscall
        step['grouping'] = 2
        steps.append(step)
        n += 1

        # Skip homogenization for non-targets
        if is_target and cfg.WSC_MAX_CHANNELS < cfg.WSC_DMASK_CHANNELSOUT:
            step = {}
            step['step'] = n
            step['comment'] = f'Fix image naming for INTER images'
            step['dependency'] = n - 1
            step['id'] = 'FXCMI' + code
            step['glam_config'] = cfg.GLAM_SMALL
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_DMASK_CHANNELSOUT} {inter_img_prefix}'
            step['syscall'] = syscall
            steps.append(step)
            n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Run wsclean-predict on intermask model for source {}'.format(fieldname)
        step['dependency'] = n - 1 
        step['id'] = 'PRCMI'+code
        step['glam_config'] = cfg.GLAM_MEDIUM
        step['slurm_config'] = cfg.SLURM_PREDICT
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        cores = step['glam_config']['CPUS']
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall = prefix + 'python3 '+TOOLS+'/fix_nan_models.py ' + inter_img_prefix + '\n\n'
        syscall += prefix + gen.generate_syscall_predict(msname = myms,
                imgname = inter_img_prefix,
                cores = cores,
                absmem = absmem,
                tempdir = temp_dir)
        #syscall += '\n\nfind ' + img_dir + ' -type f -name "*' + os.path.basename(inter_img_prefix) + '*" ! -name "*MFS-image.fits" -delete'
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Run Quartical self-calibration (stage 2) on {}'.format(fieldname)
        step['dependency'] = n - 1
        step['id'] = 'CL2GC'+code
        step['glam_config'] = cfg.GLAM_COREHEAVY
        cores = step['glam_config']['CPUS']
        syscall = CONTAINER_RUNNER + QUARTICAL_CONTAINER+' ' if USE_SINGULARITY else ''
        # Build extra args with UV range selection if needed
        extra_args = f'output.gain_directory={gain_outdir_2GC} output.log_directory={log_outdir_2GC}'
        if minuvl != '' or maxuvl != '':
            min_uv = minuvl if minuvl != '' else '0'
            max_uv = maxuvl if maxuvl != '' else '0'
            extra_args += f' input_ms.select_uv_range=[{min_uv},{max_uv}]'
        syscall += gen.generate_syscall_quartical(yaml = cfg.CAL_2GC_YAML,
                myms = myms,
                extra_args = extra_args)
        step['syscall'] = syscall
        steps.append(step)
        n += 1

        step = {}
        step['step'] = n
        step['comment'] = f'Run wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated) for {fieldname}'
        step['dependency'] = n - 1
        step['id'] = 'WSCMA'+code
        step['glam_config'] = cfg.GLAM_WSC_SOURCE if is_target else cfg.GLAM_WSC_CAL
        step['slurm_config'] = cfg.SLURM_WSCLEAN
        step['pbs_config'] = cfg.PBS_WSCLEAN
        absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
        cores = step['glam_config']['CPUS']
        syscall = ''
        prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
        # Adjust parameters based on field type
        imsize_pcal = CALIBRATOR_IMSIZE if not is_target else cfg.WSC_IMSIZE
        maxchan_pcal = CALIBRATOR_MAXCHAN if not is_target else cfg.WSC_MAX_CHANNELS
        chanout_pcal = CALIBRATOR_CHANNELSOUT if not is_target else cfg.WSC_PCAL_CHANNELSOUT
        weight_pcal = 'uniform' if not is_target else cfg.WSC_WEIGHT
        imcall = gen.generate_syscall_wsclean(mslist = [myms],
                tempdir = temp_dir,
                imgname = pcal_img_prefix,
                datacol = 'CORRECTED_DATA',
                mask = mask,
                chanout = chanout_pcal,
                tukeytaper=tukeytaper,
                minuvl = minuvl,
                maxuvl = maxuvl,
                nomodel=True,
                sourcelist = False,
                fitspectralpol = 0,
                weight = weight_pcal,
                imsize = imsize_pcal,
                maxchan = maxchan_pcal,
                cores = cores,
                absmem = absmem)
        for call in imcall: 
            syscall += prefix + call + '\n\n'
            syscall += f'rm -f {pcal_img_prefix}*dirty*\n\n'
        step['syscall'] = syscall
        step['grouping'] = 2
        steps.append(step)
        n += 1

        # Skip homogenization for non-targets
        if is_target and (cfg.WSC_MAX_CHANNELS < cfg.WSC_PCAL_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM):
            step = {}
            step['step'] = n
            step['comment'] = f'Homogenize the PCAL resolution across frequency channels'
            step['dependency'] = n - 1
            step['id'] = 'HOCMA' + code
            step['glam_config'] = cfg.GLAM_MEDIUM
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall =  prefix + f'python3 {cfg.TOOLS}/fix_image_naming.py {cfg.WSC_PCAL_CHANNELSOUT} {pcal_img_prefix}\n\n'
            syscall +=  prefix + f'python3 {cfg.TOOLS}/homogenize_beams.py {pcal_img_prefix}\n\n'
            syscall += f'find . -name "{os.path.basename(pcal_img_prefix)}*homogenized.fits" ! -name "*MFS*" -delete'
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
            step['glam_config'] = cfg.GLAM_MEDIUM
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            cores = step['glam_config']['CPUS']
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            # Calculate smaller cellsize for higher resolution (0.8 x default)
            cellsize_val = float(cfg.WSC_CELLSIZE.replace('asec','').replace('arcsec',''))
            cellsize_highres = f'{cellsize_val * 0.8:.2f}asec'
            
            # First do blind imaging to create mask with correct cell size
            imcall_blind = gen.generate_syscall_wsclean(mslist = [myms],
                tempdir = temp_dir,
                imgname = uniform_img_prefix,
                datacol = 'CORRECTED_DATA',
                chanout = cfg.WSC_BLIND_CHANNELSOUT,
                nomodel = True,
                pol = 'I',
                intervalsout = False,
                mfweight = True,
                localrms = cfg.WSC_LOCALRMS_BLIND,
                automask = cfg.WSC_AUTOMASK_BLIND,
                autothreshold = cfg.WSC_AUTOTHRESHOLD_BLIND,
                tukeytaper=False,
                minuvl = '',
                maxuvl = '',
                weight=cfg.WSC_WEIGHT_HIGHRES,
                cellsize = cellsize_highres,
                cores = cores,
                absmem = absmem)
            for call in imcall_blind:
                syscall += prefix + call + '\n\n'
                syscall += f'rm -f {uniform_img_prefix}*dirty*\n\n'
            
            # Generate mask from blind image
            syscall += prefix + gen.generate_syscall_breizorro(restoredimage = f"{uniform_img_prefix}-MFS-image.fits", 
                                                      outfile = f"{uniform_img_prefix}-MFS-image.mask.fits")[0] + '\n\n'
            uniform_mask = f"{uniform_img_prefix}-MFS-image.mask.fits"
            
            # Now do masked imaging with the new mask
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                tempdir = temp_dir,
                imgname = uniform_img_prefix,
                datacol = 'CORRECTED_DATA',
                mask = uniform_mask,
                chanout = cfg.WSC_BLIND_CHANNELSOUT,
                nomodel=True,
                intervalsout=False,
                weight=cfg.WSC_WEIGHT_HIGHRES,
                mfweight=True,
                tukeytaper=False,  # Hardcoded no taper
                minuvl = '',       # Hardcoded no minuvl
                maxuvl = '',       # Hardcoded no maxuvl
                pol='I',
                cellsize = cellsize_highres,
                sourcelist = False,
                cores = cores,
                absmem = absmem)
            for call in imcall: 
                syscall += prefix + call + '\n\n'
                syscall += f'rm -f {uniform_img_prefix}*dirty*\n\n'
            # Clean up intermediate images (keep only MFS image.fits files)
            syscall += f'find {img_dir} -type f -name "*{os.path.basename(uniform_img_prefix)}*" ! -name "*MFS*image.fits" -delete\n\n'
            step['syscall'] = syscall
            step['grouping'] = None
            steps.append(step)
            n += 1

        # No-taper image only for targets
        if is_target:
            step = {}
            step['step'] = n
            step['comment'] = 'Run no-taper wsclean, masked deconvolution of the CORRECTED_DATA (self-calibrated) for {}'.format(fieldname)
            step['dependency'] = n - 1
            step['id'] = 'WSNOT'+code
            step['glam_config'] = cfg.GLAM_MEDIUM
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
            cores = step['glam_config']['CPUS']
            syscall = ''
            prefix = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            imcall = gen.generate_syscall_wsclean(mslist = [myms],
                tempdir = temp_dir,
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
                maxuvl = '',       # Hardcoded no maxuvl
                pol='IQUV',        
                sourcelist = False,
                cores = cores,
                absmem = absmem)
            for call in imcall: 
                syscall += prefix + call + '\n\n'
                syscall += f'rm -f {notaper_img_prefix}*dirty*\n\n'
            # Clean up intermediate images (keep only MFS image.fits files)
            syscall += f'find {img_dir} -type f -name "*{os.path.basename(notaper_img_prefix)}*" ! -name "*MFS*image.fits" -delete\n\n'
            step['syscall'] = syscall
            step['grouping'] = None
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
            syscall += f"python3 {cfg.TOOLS}/make_pol_images.py {img_dir} {filename_fieldname}"
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
