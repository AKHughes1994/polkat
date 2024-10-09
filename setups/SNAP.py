#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


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
    gen.setup_dir(cfg.INTERVALS)


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

    with open('prefields_info.json') as f:
        prefields_info = json.load(f)
    
    target_names = project_info['target_names']
    myms = project_info['working_ms']

    if cfg.SNAP_FIELDS != '':
        target_names = cfg.SNAP_FIELDS.split(',')
    
    # ------------------------------------------------------------------------------
    #
    # SNAP recipe definition
    #
    # ------------------------------------------------------------------------------


    target_steps = []
    ii = 1

    # Initialize workflow by sequentially spliting
    steps = []
    kill_file = SCRIPTS+'/kill_snap_split_jobs_.sh'    
    step_i = 0    
    for tt in range(0,len(target_names)):
        
        targetname   = target_names[tt]

        if targetname not in prefields_info['field_names']:

            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('MS')+'not found, skipping')

        else:

            # Logistics
            code = gen.get_target_code(targetname)
            dependency = step_i - 1
            if step_i == 0:
                dependency = None

            step = {}
            step['step'] = step_i
            step['comment'] = 'Splitting out field '+targetname
            step['dependency'] = dependency
            step['id'] = 'SNPTS'+code
            syscall = CONTAINER_RUNNER + CASA_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_casa(casascript=cfg.OXKAT+f'/SNAP_split_sources.py {targetname}')
            step['syscall'] = syscall
            steps.append(step)
            step_i += 1
            last_split_code = 'SNPTS'+code # Save this as a variable to set dependencies for imaging
        
    target_steps.append((steps,kill_file,'Split MS Files'))


    # Loop over targets for imaging -- run imaging concurrently 
    for tt in range(0,len(target_names)):

        targetname   = target_names[tt]

        if targetname not in prefields_info['field_names']:

            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('MS')+'not found, skipping')

        else:

            # Logistics
            code = gen.get_target_code(targetname)
            steps = []        
            filename_targetname = gen.scrub_target_name(targetname)

            # Target-specific kill file
            kill_file = SCRIPTS+'/kill_snap_jobs_'+filename_targetname+'.sh'

            # Define target ms 
            target_ms = myms.replace('.ms', f'_{targetname}.ms')

            gen.print_spacer()
            print(gen.col('Target')+targetname)
            print(gen.col('Measurement Set')+myms)
            print(gen.col('Code')+code)
            
            chans_out = cfg.SNAP_CHANS
            num_frequency_splits = cfg.SNAP_FREQ_SPLITS
            best_mask_path = cfg.SNAP_MODEL_MASK_PATH
            
            snap_imaging_mask = cfg.SNAP_MASK_PATH
            snap_deconv = cfg.SNAP_DECONV
            
            if snap_deconv and (snap_imaging_mask == ''):
                print('WARNING WILL DECONVOLVE SNAPSHOT IMAGES WITH NO MASK')
            
            #always make an 8 channel model, no matter what the desired channels out of snapshot actually is
            if chans_out < 8:
                chans_out = 8
            
            if chans_out == 1:
                fitspectralpol = 0
                join_chans = False
            elif chans_out == 2:
                join_chans = True
                fitspectralpol = 1
            elif chans_out <= 4:
                join_chans = True
                fitspectralpol = 2
            elif chans_out <= 6:
                join_chans = True
                fitspectralpol = 3
            else:
                join_chans = True
                fitspectralpol = 4
            
            if best_mask_path == '':
                print('BEST MASK NOT SPECIFIED, DEFAULT MASK USED')
                best_mask_path = 'fits'
            
            # Image prefix
            #blind_img_prefix = IMAGES+'/img_'+target_ms+'_snapblind'+str(chans_out)+'chan'
            mask_img_prefix = IMAGES+'/img_'+target_ms+'_snapmask'+str(chans_out)+'chan'
            
            
            
            if cfg.SNAP_POL == True:
                pols = 'IQUV'
                join_pols = True
            else:
                pols = 'I'
                join_pols = False
                
            
            
            step_i = 0
            
            if num_frequency_splits == 1:
                
                step = {}
                step['step'] = step_i
                step['comment'] = 'Masked wsclean on DATA column of '+target_ms
                step['dependency'] = last_split_code
                step['id'] = 'SNPMA'+code
                step['slurm_config'] = cfg.SLURM_WSCLEAN
                step['pbs_config'] = cfg.PBS_WSCLEAN
                absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                syscall = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                syscall += gen.generate_syscall_wsclean(mslist = [target_ms],
                            imgname = mask_img_prefix,
                            mask = best_mask_path,
                            datacol = 'DATA',
                            localrms = False,
                            field='0',
                            absmem = absmem,
                            chanout = chans_out,
                            joinchannels = join_chans,
                            fitspectralpol=fitspectralpol,
                            pol = pols,
                            joinpolarizations = join_pols)
                step['syscall'] = syscall
                steps.append(step)
                step_i += 1
            
            else:
            
                #at the moment everything must be divisible!!!! ALso assuming averaged to 1024 channels
                pre_nchans = cfg.PRE_NCHANS #number of channels the MS has been averaged to
                chans_per_part = (pre_nchans/num_frequency_splits)
                
                if not chans_per_part.is_integer():
                    sys.exit('Channels in MS not divisble by the number of frequency splits')
                  
                
                for ei in range(num_frequency_splits):
                    step = {}
                    step['step'] = step_i
                    step['comment'] = 'Masked wsclean on DATA column of '+target_ms
                    
                    if ei == 0:
                        step['dependency'] = last_split_code
                    else:
                        step['dependency'] = step_i - 1
                    
                    step['id'] = 'SNPMA'+code+'part'+str(ei)
                    step['slurm_config'] = cfg.SLURM_WSCLEAN
                    step['pbs_config'] = cfg.PBS_WSCLEAN
                    absmem = gen.absmem_helper(step,INFRASTRUCTURE,cfg.WSC_ABSMEM)
                    syscall = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
                    syscall += gen.generate_syscall_wsclean(mslist = [target_ms],
                                imgname = mask_img_prefix+'part'+str(ei),
                                mask = best_mask_path,
                                datacol = 'DATA',
                                startchan = 0 + (ei*chans_per_part),
                                endchan = chans_per_part + (ei*chans_per_part), #not minus 1 here because the endchan parameter is exclusionary of last channel
                                localrms = False,
                                field='0',
                                absmem = absmem,
                                chanout = int(chans_out/num_frequency_splits),
                                joinchannels = True,
                                pol = pols,
                                joinpolarizations = join_pols)
                    step['syscall'] = syscall
                    steps.append(step)
                    step_i += 1

            step = {}
            step['step'] = step_i
            step['comment'] = 'Making individual ms files for the scans of '+targetname
            step['dependency'] = step_i - 1
            step['id'] = 'SNPSS'+code
            syscall = CONTAINER_RUNNER + CASA_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += gen.generate_syscall_casa(casascript=cfg.OXKAT+f'/SNAP_split_scans.py {target_ms}')
            step['syscall'] = syscall
            steps.append(step)
            step_i += 1

            step = {}
            step['step'] = step_i
            step['comment'] = 'Performing UV-subtraction of the static sky model for '+targetname
            step['dependency'] = step_i - 1
            step['id'] = 'SNPUV'+code
            syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += 'python3 '+cfg.OXKAT+f'/SNAP_uvsub.py {targetname}'
            step['syscall'] = syscall
            steps.append(step)
            step_i += 1

            step = {}
            step['step'] = step_i
            step['comment'] = 'Performing per-interval dirty imaging of the uvsubtracted visbilities for '+targetname
            step['dependency'] = step_i - 1
            step['id'] = 'SNPIT'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            syscall = CONTAINER_RUNNER+WSCLEAN_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += 'python3 '+cfg.OXKAT+f'/SNAP_intervals.py {targetname}'
            step['syscall'] = syscall
            steps.append(step)
            step_i += 1

            step = {}
            step['step'] = step_i
            step['comment'] = 'Performing static model restoration for '+targetname
            step['dependency'] = step_i - 1
            step['id'] = 'SNPRM'+code
            step['slurm_config'] = cfg.SLURM_WSCLEAN
            step['pbs_config'] = cfg.PBS_WSCLEAN
            syscall = CONTAINER_RUNNER+PYTHON3_CONTAINER+' ' if USE_SINGULARITY else ''
            syscall += 'python3 '+cfg.OXKAT+f'/SNAP_restore.py {targetname}'
            step['syscall'] = syscall
            steps.append(step)
            step_i += 1

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

            step_id = step['id']
            id_list.append(step_id)
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
