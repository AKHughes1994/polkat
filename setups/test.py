#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk


import glob
import json
import os.path as o
import sys
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import generate_imaging as gi
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


    with open('project_info.json') as f:
        project_info = json.load(f)

    band = project_info['band']

    OXKAT = cfg.OXKAT
    DATA = cfg.DATA
    IMAGES = cfg.IMAGES
    SCRIPTS = cfg.SCRIPTS
    TOOLS = cfg.TOOLS

    INFRASTRUCTURE, CONTAINER_PATH = gen.set_infrastructure(sys.argv)
    if CONTAINER_PATH is not None:
        CONTAINER_RUNNER='singularity exec '
    else:
        CONTAINER_RUNNER=''


    CASA_CONTAINER       = gen.get_container(CONTAINER_PATH,cfg.CASA_PATTERN,USE_SINGULARITY)
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

    img_prefix = f"{cfg.IMAGES}/img_test_diagnostic"

    steps = []
    step = {}
    step_n = 0
    step['step'] = step_n
    step['comment'] = f'Image calibrator as a diagnostic for calibration systematics'
    step['dependency'] = None if step_n == 0 else step_n - 1
    step['id'] = 'DIAGN'
    syscall = gi.generate_syscall_wsclean(mslist = [myms],
                imgname = img_prefix,
                datacol = 'CORRECTED_DATA',
                field = '1',
                weight=cfg.WSC_WEIGHT_CAL,
                imsize = cfg.WSC_CAL_IMSIZE,
                chanout = cfg.WSC_IMAGE_CHANNELSOUT,
                pol='IQUV',
                mask = False,
                automask = 5.0,
                autothreshold = 1.0,
                splitpol = True,
                localrms= True,
                threshold = False,
                nomodel  = True,
                sourcelist = False)
    prefix = '\n\n' + CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
    step['syscall'] = prefix + prefix.join(syscall)
    steps.append(step)
    step_n += 1


    if cfg.WSC_MAX_CHANNELS < cfg.WSC_IMAGE_CHANNELSOUT or cfg.WSC_HOMOGENIZEBEAM:
        step = {}
        step['step'] = step_n
        step['comment'] = f'Homogenize the resolution across frequency channels'
        step['dependency'] = None if step_n == 0 else step_n - 1
        step['id'] = 'HOMOG'
        syscall =  f'python3 {cfg.TOOLS}/homogenize_beams.py {img_prefix}'
        prefix = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        step['syscall'] = prefix + syscall
        steps.append(step)
        step_n += 1

    step = {}
    step['step'] = step_n
    step['comment'] = 'Apply primary beam correction to test (BLIND) image'
    step['dependency'] = None if step_n == 0 else step_n - 1
    step['id'] = 'PBBLD'+code
    syscall = ''
    images = glob.glob(f'{img_prefix}*-MFS*-image.fits')
    if glob.glob(f'{img_prefix}*-MFS*-image.homogenized.fits') != []:
        images.extend(glob.glob(f'{img_prefix}*-MFS*-image.homogenized.fits'))
    for image in images:
        syscall += CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
        syscall += 'python3 '+TOOLS+'/pbcor_katbeam.py --band '+band[0]+f' {image}\n'
    step['syscall'] = syscall
    steps.append(step)
    step_n += 1

    

    # ------------------------------------------------------------------------------
    #
    # Write the run file and kill file based on the recipe
    #
    # ------------------------------------------------------------------------------


    submit_file = 'submit_1GC_jobs.sh'
    kill_file = cfg.    SCRIPTS+'/kill_1GC_jobs.sh'

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
