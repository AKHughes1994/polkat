#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import json
import os.path as o
import sys
import random
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg


def main():

    USE_SINGULARITY = cfg.USE_SINGULARITY

    gen.preamble()
    print(gen.col()+'0GC (fixing project info) setup')
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

    # Determine MS to pass to scan_times, preferring working_ms from project_info.
    scantimes_ms = None
    if o.isfile('project_info.json'):
        with open('project_info.json') as f:
            project_info = json.load(f)
        scantimes_ms = project_info.get('working_ms')
    if not scantimes_ms:
        myms = sorted(glob.glob('*.ms'))
        if len(myms) > 0:
            scantimes_ms = myms[0]
    if not scantimes_ms:
        print(gen.col('ERROR')+'Could not determine MS for scan_times (project_info working_ms or *.ms).')
        gen.print_spacer()
        sys.exit()


    # ------------------------------------------------------------------------------
    #
    # 0GC recipe definition
    #
    # ------------------------------------------------------------------------------


    n = 0
    steps = []
    step = {}
    step['step'] = n
    step['comment'] = 'Fix project info'
    step['dependency'] = None if n == 0 else n - 1  
    step['id'] = 'FIXPI'+ str(random.randint(100,999))
    syscall = CONTAINER_RUNNER+CASA_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += gen.generate_syscall_casa(casascript=cfg.TOOLS+'/fix_project_info.py')
    syscall += '\n'
    syscall += CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
    syscall += ' python3 '+cfg.TOOLS+'/scan_times.py '+scantimes_ms
    step['syscall'] = syscall
    steps.append(step)
    n += 1

   
    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    submit_file = 'submit_0GC_jobs.sh'
    kill_file = cfg.SCRIPTS+'/kill_0GC_jobs.sh'

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
