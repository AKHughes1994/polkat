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
    print(gen.col()+'1GC (referenced calibration) setup')
    gen.print_spacer()

    if cfg.PRE_FIELDS != '':
        print(gen.col('Field selection')+cfg.PRE_FIELDS)
    if cfg.PRE_SCANS != '':
        print(gen.col('Scan selection')+cfg.PRE_SCANS)
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
    gen.setup_dir(cfg.GAINPLOTS)
    gen.setup_dir(cfg.VISPLOTS)

    INFRASTRUCTURE, CONTAINER_PATH = gen.set_infrastructure(sys.argv)
    if CONTAINER_PATH is not None:
        CONTAINER_RUNNER='singularity exec '
    else:
        CONTAINER_RUNNER=''


    CASA_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.CASA_PATTERN,USE_SINGULARITY)
    WSCLEAN_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.WSCLEAN_PATTERN,USE_SINGULARITY)
    SHADEMS_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.SHADEMS_PATTERN,USE_SINGULARITY)
    PYTHON3_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.PYTHON3_PATTERN,USE_SINGULARITY)


    # ------------------------------------------------------------------------------
    #
    # 1GC recipe definition
    #
    # ------------------------------------------------------------------------------


    with open('project_info.json') as f:
        project_info = json.load(f)

    myms  = project_info['working_ms']
    code = gen.get_code(myms)

    n = 0
    steps = []
    if cfg.CAL_1GC_RENAME_FEEDS:
        step = {}
        step['step'] = 0
        step['comment'] = 'Swap H and V polarization labels'
        step['dependency'] = None
        step['id'] = 'SWPHV'+code
        syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += f"python3 {cfg.OXKAT}/1GC_01_correct_parang.py"
        step['syscall'] = syscall
        step['glam_config'] = cfg.GLAM_CASA
        steps.append(step)
        n += 1

        step = {}
        step['step'] = n
        step['comment'] = 'Set the feed offset angle to zero'
        step['dependency'] = n - 1
        step['id'] = 'ZEROF'+code
        syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += gen.generate_syscall_casa(casascript=cfg.OXKAT+'/1GC_02_zero_feed.py ' + myms)
        step['syscall'] = syscall
        step['glam_config'] = cfg.GLAM_CASA
        steps.append(step)
        n += 1    

    step = {}
    step['step'] = n
    step['comment'] = 'Apply basic flagging steps to all fields'
    step['dependency'] = None if n == 0 else n - 1  
    step['id'] = 'FGBAS'+code
    syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_casa(casascript=cfg.OXKAT+'/1GC_03_casa_basic_flags.py')
    step['syscall'] = syscall
    step['glam_config'] = cfg.GLAM_CASA
    steps.append(step)
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = 'Run auto-flaggers on calibrators'
    step['dependency'] = n - 1
    step['id'] = 'FGCAL'+code
    syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_casa(casascript=cfg.OXKAT+'/1GC_04_casa_autoflag_cals_DATA.py')
    step['syscall'] = syscall
    step['glam_config'] = cfg.GLAM_CASA
    steps.append(step)
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = 'Using reference calibrators perform full polarization calibration'
    step['dependency'] = n - 1
    step['id'] = 'CL1GC'+code
    syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_casa(casascript=cfg.OXKAT+'/1GC_05_casa_refcal.py')
    step['syscall'] = syscall
    step['glam_config'] = cfg.GLAM_COREHEAVY
    steps.append(step)
    cl1gc_step = n
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = 'Plot the final gain tables'
    step['dependency'] = cl1gc_step
    step['id'] = 'PLTAB'+code
    syscall = CONTAINER_RUNNER+SHADEMS_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += 'python3 '+cfg.OXKAT+'/1GC_06_plot_gaintables.py cal_1GC_*'
    step['syscall'] = syscall
    step['glam_config'] = cfg.GLAM_CASA
    steps.append(step)
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = 'Plot the corrected calibrator visibilities'
    step['dependency'] = cl1gc_step
    step['id'] = 'PLVIS'+code
    syscall = CONTAINER_RUNNER+SHADEMS_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += 'python3 '+cfg.OXKAT+'/1GC_07_plot_visibilities.py'
    step['syscall'] = syscall
    step['glam_config'] = cfg.GLAM_MEDIUM
    steps.append(step)
    plvis_step = n
    n += 1

    step = {}
    step['step'] = n
    step['comment'] = 'Split out the target MS files'
    step['dependency'] = plvis_step
    step['id'] = 'SPTRG'+code
    syscall = CONTAINER_RUNNER+SHADEMS_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_casa(casascript=cfg.OXKAT+'/1GC_08_casa_split_targets.py')
    step['syscall'] = syscall
    step['glam_config'] = cfg.GLAM_CASA
    steps.append(step)

   
    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    submit_file = 'submit_1GC_jobs.sh'
    kill_file = cfg.SCRIPTS+'/kill_1GC_jobs.sh'

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
