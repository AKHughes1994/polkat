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
    print(gen.col()+'RM Synthesis analysis of the polarization properties')
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

    with open(cfg.RMSYN_INFO_FILE, 'r') as j:
        rmsynth_info = json.load(j)
    for k, name in enumerate(rmsynth_info['image_directory']):
        print(gen.col() + 'Trying to fit {} in directory {}'.format(rmsynth_info['source_name'][k], rmsynth_info['image_directory'][k]))
    gen.print_spacer()

    INFRASTRUCTURE, CONTAINER_PATH = gen.set_infrastructure(sys.argv)
    if CONTAINER_PATH is not None:
        CONTAINER_RUNNER='singularity exec '
    else:
        CONTAINER_RUNNER=''

    PYTHON3_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.PYTHON3_PATTERN,USE_SINGULARITY)
    SPINIFEX_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.SPINIFEX_PATTERN,USE_SINGULARITY)
    CASA_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.CASA_PATTERN,USE_SINGULARITY)
    TRICOLOUR_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.TRICOLOUR_PATTERN,USE_SINGULARITY)


    # ------------------------------------------------------------------------------
    #
    # RMSynth recipe definition
    #
    # ------------------------------------------------------------------------------


    with open('project_info.json') as f:
        project_info = json.load(f)

    myms  = project_info['working_ms']
    code = gen.get_code(myms)
    steps = []

    # Single-target run: exactly one science target survived the PRE
    # averaging step (i.e. is present in working_names). In that case we
    # add a final step that summarizes the key MFS/RM synthesis outputs
    # into one text file, instead of leaving them scattered across RESULTS.
    target_names   = project_info.get('target_names', [])
    working_names  = project_info.get('working_names', [])
    single_targets = [t for t in target_names if t in working_names]
    if len(single_targets) == 1:
        print(gen.col('Single Target')+f'Detected single target {single_targets[0]} -- will write a combined output summary')
    else:
        print(gen.col('Single Target')+f'Not applicable ({len(single_targets)} targets present)')


    step_i = 0
    step = {}
    step['step'] = step_i
    step['comment'] = 'Extract polarization properties using rmsynth_info.txt (found in data/rmsynth directory)'
    step['dependency'] = None
    step['slurm_config'] = cfg.SLURM_RM
    step['id'] = 'POEXT'+code
    syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
    # syscall += gen.generate_syscall_casa_short(casascript=cfg.OXKAT+f'/RMSYNTH_01_extract_fluxes.py')
    syscall += 'python-pycasa '+cfg.OXKAT+'/RMSYNTH_01_extract_fluxes.py'
    step['syscall'] = syscall
    steps.append(step)
    step_i += 1

    if cfg.CAL_1GC_DIAGNOSTICS:
        step = {}
        step['step'] = step_i
        step['comment'] = 'Calculate systematic effects by performing image plane analysis on polarization/primary calibrator'
        step['dependency'] = step_i - 1
        step['slurm_config'] = cfg.SLURM_RM
        step['id'] = 'POSYS'+code
        syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += gen.generate_syscall_casa_short(casascript=cfg.OXKAT+f'/RMSYNTH_01B_systematics.py')
        step['syscall'] = syscall
        steps.append(step)
        step_i += 1      

    if cfg.POLANG_NAME != '':
        step = {}
        step['step'] = step_i
        step['comment'] = 'Run RM Synthesis on all of the extracted IQUV curves'
        step['dependency'] = step_i - 1
        step['slurm_config'] = cfg.SLURM_RM
        step['id'] = 'RMSYN'+code
        syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += 'python3 '+cfg.OXKAT+'/RMSYNTH_02_run_rmsynth.py'
        step['syscall'] = syscall
        steps.append(step)
        step_i += 1

    step = {}
    step['step'] = step_i
    step['comment'] = 'Run SPINIFEX on all sources'
    step['dependency'] = step_i - 1
    step['slurm_config'] = cfg.SLURM_RM
    step['id'] = 'SPFEX'+code
    syscall = CONTAINER_RUNNER+SPINIFEX_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += 'python-spinifex '+cfg.OXKAT+'/RMSYNTH_03_run_SPINIFEX.py'
    step['syscall'] = syscall
    steps.append(step)
    step_i += 1

    if len(single_targets) == 1:
        step = {}
        step['step'] = step_i
        step['comment'] = f'Summarize MFS/RM synthesis outputs for single target {single_targets[0]}'
        step['dependency'] = step_i - 1
        step['slurm_config'] = cfg.SLURM_RM
        step['id'] = 'RMSUM'+code
        syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += 'python3 '+cfg.OXKAT+'/RMSYNTH_04_summarize_target.py'
        step['syscall'] = syscall
        steps.append(step)
        step_i += 1

    # step = {}
    # step['step'] = step_i
    # step['comment'] = 'Run ALBUS on all extracted sources'
    # step['dependency'] = step_i - 1
    # step['slurm_config'] = cfg.SLURM_RM
    # step['id'] = 'ALBUS'+code
    # syscall = CONTAINER_RUNNER+ALBUS_CONTAINER+' ' if USE_SINGULARITY else ''
    # syscall += 'python3 '+cfg.OXKAT+'/RMSYNTH_03_run_ALBUS.py'
    # step['syscall'] = syscall
    # steps.append(step)
    # step_i += 1


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
