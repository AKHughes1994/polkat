#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import os
import os.path as o
import sys
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg


def main():

    gen.print_spacer()
    print(gen.col('MAKE SOME MOVIES'))
    gen.print_spacer()

    USE_SINGULARITY = cfg.USE_SINGULARITY
    INFRASTRUCTURE, CONTAINER_PATH = gen.set_infrastructure(sys.argv)
    POLKAT_CONTAINER = gen.get_container(CONTAINER_PATH,cfg.MOVIE_PATTERN,True)
    if CONTAINER_PATH is not None:
        CONTAINER_RUNNER='singularity exec '
    else:
        CONTAINER_RUNNER=''

    intervals = sorted(glob.glob(cfg.INTERVALS+'*/'))
    rootdir = os.getcwd()

    gen.print_spacer()
    print(gen.col('Found INTERVALs '), intervals)
    gen.print_spacer()    


    steps = []
    for k, interval in enumerate(intervals):

        step = {}
        step['step'] = k
        step['comment'] = f'Making movie for images in {interval}'
        step['dependency'] = None if k == 0 else k-1
        step['id'] = 'MKMOV'+ str(k)
        syscall = CONTAINER_RUNNER+POLKAT_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += f"python3 {cfg.TOOLS}/make_movie.py {interval}"
        step['syscall'] = syscall
        steps.append(step)

    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    submit_file = 'submit_movie_jobs.sh'
    kill_file = cfg.SCRIPTS+'/kill_movie_jobs.sh'

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

if __name__ == "__main__":

    main()
