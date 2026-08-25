#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

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
    print(gen.col()+'Extraction of the polarization flux densities')
    if cfg.CAL_1GC_DIAGNOSTICS:
        print(gen.col() + 'Including systematic calculations using Primary and Polarization Angle Calibrator')
    gen.print_spacer()

    # ------------------------------------------------------------------------------
    #
    # Setup paths, required containers, infrastructure
    #
    # ------------------------------------------------------------------------------


    gen.setup_dir(cfg.LOGS)
    gen.setup_dir(cfg.SCRIPTS)
    gen.setup_dir(cfg.GAINTABLES)
    gen.setup_dir(cfg.IMAGES)
    gen.setup_dir(cfg.RESULTS)
    gen.setup_dir(cfg.GAINPLOTS)
    gen.setup_dir(cfg.VISPLOTS)

    INFRASTRUCTURE, CONTAINER_PATH = gen.set_infrastructure(sys.argv)
    if CONTAINER_PATH is not None:
        CONTAINER_RUNNER='singularity exec '
    else:
        CONTAINER_RUNNER=''

    PYTHON3_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.PYTHON3_PATTERN,USE_SINGULARITY)
    # ALBUS_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.ALBUS_PATTERN,USE_SINGULARITY)  # DEPRECATED: ALBUS step below is commented out; polkat-albus image no longer exists
    CASA_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.CASA_PATTERN,USE_SINGULARITY)
    TRICOLOUR_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.TRICOLOUR_PATTERN,USE_SINGULARITY)
    SPINIFEX_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.PYTHON3_PATTERN,USE_SINGULARITY)


    # ------------------------------------------------------------------------------
    #
    # Build field list from project_info.json
    #
    # ------------------------------------------------------------------------------

    with open('project_info.json') as f:
        project_info = json.load(f)
    
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

    # Build unified field list
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

    # Skip secondary calibrators -- targets, the primary, and the polarization-angle
    # calibrator are always processed regardless of this flag.
    if cfg.CAL_SKIP_CALS:
        n_cals = len([f for f in field_list if f['is_secondary']])
        field_list = [f for f in field_list if not f['is_secondary']]
        print(gen.col('RMSYNTH Cals')+f'CAL_SKIP_CALS is True -- skipping {n_cals} secondary field/scan(s)')

    # ------------------------------------------------------------------------------
    #
    # Check that images actually exist for each remaining field before building the
    # job graph. A field can still be listed in project_info.json with no images on
    # disk (e.g. 2GC hasn't been run yet) -- extraction would just fail at runtime.
    # Drop those fields here instead, and say so.
    #
    # ------------------------------------------------------------------------------

    fields_with_images = []
    for field in field_list:
        img_dir = cfg.IMAGES + '/' + gen.scrub_target_name(field['name'])
        has_images = o.isdir(img_dir) and len(glob.glob(img_dir + '/*-I-image.fits')) > 0
        if has_images:
            fields_with_images.append(field)
        else:
            print(gen.col('NOTE')+f"No images for {field['name']} -- skipping")

    n_missing = len(field_list) - len(fields_with_images)
    if n_missing:
        print(gen.col('RMSYNTH Images')+f'skipping {n_missing}, {len(fields_with_images)} remaining')
    field_list = fields_with_images

    if not field_list:
        print(gen.col('ERROR')+'No fields with images found -- nothing to extract. Run 2GC.py first.')
        sys.exit()

    # Print field information
    gen.print_spacer()
    print(gen.col('Fields to process:'))
    for field in field_list:
        field_type = 'TARGET' if field['is_target'] else 'PCAL'
        print(gen.col(f"  {field['name']}") + f" (ID: {field['id']}, Type: {field_type})")

    # ------------------------------------------------------------------------------
    #
    # Verify every target has a known position before building the job graph.
    # Matched by fuzzy/normalized similarity (gen.match_field_name()), the same
    # rule RMSYNTH_01_extract_fluxes.py's read_position_file() uses at runtime --
    # a target that passes here is guaranteed to resolve the same way there.
    # read_position_file() hard-fails if it doesn't find a match, so this exists
    # purely to catch that before a job is ever submitted, not after.
    #
    # ------------------------------------------------------------------------------

    position_file = cfg.DATA + '/positions/XRB_pos_list.txt'
    known_names = []
    if o.exists(position_file):
        with open(position_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t', 1)
                if len(parts) == 2 and parts[0].strip():
                    known_names.append(parts[0].strip())
    else:
        print(gen.col('ERROR')+f'Position file {position_file} not found')
        sys.exit()

    print(gen.col('Position matches:'))
    missing_targets = []
    for field in field_list:
        if not field['is_target']:
            continue
        best_name, best_ratio, matched = gen.match_field_name(field['name'], known_names)
        if matched:
            print(gen.col(f"  {field['name']}") + f"-> '{best_name}' ({best_ratio:.0%})")
        else:
            closest = f"closest was '{best_name}' at {best_ratio:.0%}" if best_name else 'position file is empty'
            missing_targets.append((field['name'], closest))

    if missing_targets:
        print(gen.col('ERROR')+f"Target(s) with no position-file match >= 80% similarity in {position_file}:")
        for name, closest in missing_targets:
            print(gen.col('ERROR')+f"  {name}: {closest}")
        print(gen.col('ERROR')+'Add or fix the target(s) in the position file before re-running RMSYNTH.py')
        sys.exit()

    # ------------------------------------------------------------------------------
    #
    # RMSynth recipe definition
    #
    # ------------------------------------------------------------------------------


    myms = project_info['working_ms']
    code = gen.get_code(myms)
    steps = []
    extraction_job_ids = []

    step_i = 0
    
    # Create one extraction job per field (all with no dependencies)
    for field in field_list:
        fieldname = field['name']
        field_code = str(field['name'])[-3:]  # Last 3 digits of field ID
        
        step = {}
        step['step'] = step_i
        step['comment'] = f'Extract polarization properties for {fieldname}'
        step['dependency'] = None
        step['glam_config'] = cfg.GLAM_SMALL  # Add glam_config
        step['id'] = f'EXT{field_code}_{step_i}'
        extraction_job_ids.append(step['id'])
        
        syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += 'python-pycasa '+cfg.OXKAT+f'/RMSYNTH_01_extract_fluxes.py {fieldname}'
        step['syscall'] = syscall
        steps.append(step)
        step_i += 1

    # SPINIFEX step depends on ALL sources in working list
    step = {}
    step['step'] = step_i
    step['comment'] = 'Run SPINIFEX on all sources'
    step['dependency'] = extraction_job_ids  # List of all extraction job IDs
    step['glam_config'] = cfg.GLAM_SMALL  # Add glam_config
    step['id'] = 'SPFEX'+code
    syscall = CONTAINER_RUNNER+SPINIFEX_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += 'python-spinifex '+cfg.OXKAT+'/RMSYNTH_03_run_SPINIFEX.py'
    step['syscall'] = syscall
    steps.append(step)
    step_i += 1

    # Run RMSynth step depends on ALL extraction jobs
    step = {}
    step['step'] = step_i
    step['comment'] = 'Run RMSynth on all rmsynth.txt files'
    step['dependency'] = extraction_job_ids  # List of all extraction job IDs
    step['glam_config'] = cfg.GLAM_SMALL  # Add glam_config
    step['id'] = 'RMSYN'+code
    syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += 'python3 '+cfg.OXKAT+'/RMSYNTH_02_run_rmsynth.py'
    step['syscall'] = syscall
    steps.append(step)
    rmsyn_step_index = step_i
    step_i += 1

    # Diagnostic step depends on ALL extraction jobs (primary/polarization-angle calibrator only)
    if cfg.CAL_1GC_DIAGNOSTICS:
        step = {}
        step['step'] = step_i
        step['comment'] = 'Calculate systematic effects by performing image plane analysis on polarization/primary calibrator'
        step['dependency'] = extraction_job_ids  # List of all extraction job IDs
        step['glam_config'] = cfg.GLAM_SMALL  # Add glam_config
        step['id'] = 'POSYS'+code
        syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += gen.generate_syscall_casa_short(casascript=cfg.OXKAT+f'/RMSYNTH_01B_systematics.py')
        step['syscall'] = syscall
        steps.append(step)
        step_i += 1

    # Summarize MFS/RM synthesis outputs for every target; depends on RMSynth (rmsynth1d/rmclean1d outputs)
    step = {}
    step['step'] = step_i
    step['comment'] = 'Summarize MFS/RM synthesis outputs for all targets'
    step['dependency'] = rmsyn_step_index
    step['glam_config'] = cfg.GLAM_SMALL  # Add glam_config
    step['id'] = 'RMSUM'+code
    syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += 'python3 '+cfg.OXKAT+'/RMSYNTH_04_summarize_target.py'
    step['syscall'] = syscall
    steps.append(step)
    step_i += 1

    # DEPRECATED:
    # ALBUS step depends on ALL extraction jobs
    # step = {}
    # step['step'] = step_i
    # step['comment'] = 'Run ALBUS on all extracted sources'
    # step['dependency'] = extraction_job_ids  # List of all extraction job IDs
    # step['glam_config'] = cfg.GLAM_MEDIUM  # Add glam_config
    # step['id'] = 'ALBUS'+code
    # syscall = CONTAINER_RUNNER+ALBUS_CONTAINER+' ' if USE_SINGULARITY else ''
    # syscall += 'python3 '+cfg.OXKAT+'/RMSYNTH_03_run_ALBUS.py'
    # step['syscall'] = syscall
    # steps.append(step)


    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    submit_file = 'submit_rmsynth_jobs.sh'
    kill_file = cfg.SCRIPTS+'/kill_rmsynth_jobs.sh'

    f = open(submit_file,'w')
    f.write('#!/usr/bin/env bash\n')
    f.write('export SINGULARITY_BINDPATH='+cfg.BINDPATH+'\n')

    if INFRASTRUCTURE == 'glam':
        # Use Glamdring-specific submission script
        f.write('set -euo pipefail\n\n')
        
        # First generate the individual job scripts
        id_list = []
        for step in steps:
            step_id = step['id']
            id_list.append(step_id)
            if step['dependency'] is not None:
                # Handle both single dependency (int) and multiple dependencies (list)
                if isinstance(step['dependency'], list):
                    dependency = step['dependency']  # Already a list of job IDs
                else:
                    dependency = steps[step['dependency']]['id']
            else:
                dependency = None
            syscall = step['syscall']
            if 'glam_config' in step.keys():
                glam_config = step['glam_config']
            else:
                glam_config = cfg.GLAM_STANDARD
            
            # Extract grouping parameter if present
            grouping = step.get('grouping', None)
            
            # Generate the job scripts using job_handler
            gen.job_handler(syscall = syscall,
                            jobname = step_id,
                            infrastructure = INFRASTRUCTURE,
                            glam_config = glam_config,
                            grouping = grouping)
        
        # Now write the submission commands
        gen.write_glam_submission_script(f, steps)
        
    else:
        # Original behavior for other infrastructures
        id_list = []

        for step in steps:

            step_id = step['id']
            id_list.append(step_id)
            if step['dependency'] is not None:
                # Handle both single dependency (int) and multiple dependencies (list)
                if isinstance(step['dependency'], list):
                    # Multiple dependencies - join with colons
                    dependency = ':'.join(step['dependency'])
                else:
                    # Single dependency - get the job ID from the step index
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


        if INFRASTRUCTURE == 'idia' or INFRASTRUCTURE == 'hippo':
            kill = '\necho "scancel "$'+'" "$'.join(id_list)+' > '+kill_file+'\n'
            f.write(kill)
        elif INFRASTRUCTURE == 'chpc':
            kill = '\necho "qdel "$'+'" "$'.join(id_list)+' > '+kill_file+'\n'
            f.write(kill)
    
    f.close()

    gen.make_executable(submit_file)

    gen.print_spacer()
    print(gen.col('Run file')+submit_file)
    gen.print_spacer()

    # ------------------------------------------------------------------------------



if __name__ == "__main__":


    main()
