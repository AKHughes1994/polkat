#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import glob
import datetime
import time
import os
import os.path as o
import sys
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg


# ------------------------------------------------------------------------

def preamble():
    print('---------------------+----------------------------------------------------------')
    print('                     |')
    print('                     | v2.0')  
    print('   p o l  k  a  t    | The poorly coded, younger brother of oxkat:')
    print('   C A S A / Q C     | Feel free to email questions/concerns to:')
    print('                     | hughesakh@gmail.com')
    print('                     | fraser.cowie@physics.ox.ac.uk')
    print('                     |')
    print('---------------------+----------------------------------------------------------')
    print(now()+'Observing band is '+cfg.BAND)
    print(col()+'Intermediate flag tables will be backed up')


def now():
    # stamp = time.strftime('[%H:%M:%S] ')
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    # msg = '\033[92m'+stamp+'\033[0m' # time in green
    msg = stamp
    return msg


def col(txt=''):
    colstr = ' '+txt.ljust(20)+'| '
    return colstr


def print_spacer():
    print('---------------------+----------------------------------------------------------')


def set_infrastructure(args):

    if len(args) == 1:
        print(col()+'Please specify infrastructure (idia / chpc / hippo / node / glam)')
        print_spacer()
        sys.exit()

    if args[1].lower() not in ['idia','chpc','hippo','node','glam']:
        print(col()+'Please specify infrastructure (idia / chpc / hippo / node / glam)')
        print_spacer()
        sys.exit()

    if args[1].lower() == 'idia':
        infrastructure = 'idia'
        CONTAINER_PATH = cfg.IDIA_CONTAINER_PATH
    elif args[1].lower() == 'chpc':
        infrastructure = 'chpc'
        CONTAINER_PATH = cfg.CHPC_CONTAINER_PATH
    elif args[1].lower() == 'node':
        infrastructure = 'node'
        CONTAINER_PATH = cfg.NODE_CONTAINER_PATH
    elif args[1].lower() == 'glam':
        infrastructure = 'glam'
        CONTAINER_PATH = cfg.GLAM_CONTAINER_PATH
    elif args[1].lower() == 'hippo':
        infrastructure = 'hippo'
        CONTAINER_PATH = None

        
    print(col('Infrastructure')+infrastructure.upper())
    if cfg.USE_SINGULARITY:
        print(col('Singularity')+'Enabled')
        print(col('Searching')+str(CONTAINER_PATH))
    else:
        print(col('Singularity')+'Not enabled')


    return infrastructure,CONTAINER_PATH


container_list = []

def get_container(pathlist,pattern,use_singularity):
    
    # Pass a list, if container isn't in it then append?]
    # Call another function to print it

    # For running without containers
    if pathlist is None: # Retain backwards compatibility with hippo fix
        return ''
    if not use_singularity:
        return ''

    ll = []
    for path in pathlist:
        # Search for a file matching pattern in path
        path = path.rstrip('/')+'/'
        ll.extend(sorted(glob.glob(path+'*'+pattern+'*img')))
        ll.extend(sorted(glob.glob(path+'*'+pattern+'*sif')))

    # Exclude stimela's casa4.7 and casarest containers
    if 'casa' in pattern.lower():
        for ii in ll:
            if 'casa47' in ii or 'casarest' in ii:
                ll.remove(ii)


    if len(ll) == 0:
        print(col(pattern)+'not found!')
        print_spacer()
        sys.exit()
    container = ll[-1]
    opstr = container.split('/')[-1]
    if len(ll) > 1:
        opstr += ' (multiple matches)'
    if opstr not in container_list:
        print(col('Found container')+opstr)
        container_list.append(opstr)
    return container


def setup_dir(DIR,relabel=False):

    # Make scripts folder if it doesn't exist

    if not o.isdir(DIR):
        os.mkdir(DIR)
    elif o.isdir(DIR) and relabel and (len(os.listdir(DIR)) != 0):
        dirs = glob.glob(DIR+'*')
        n = len(dirs)
        os.rename(DIR,DIR+str(n))
        os.mkdir(DIR)


def timenow():

    # Return a date and time string suitable for being part of a filename

    now = str(datetime.datetime.now()).replace(' ','-').replace(':','-').split('.')[0]
    return now


def get_code(myms):

    # Last three digits of the data set ID

    myms = myms.split('/')[-1]
    code = myms.split('_')[0][-3:]
    code = code.replace('-','_')
    code = code.replace('.','X')
    return code


def get_mms_code(myms):

    # Last three digits of the data set ID

    myms = myms.split('/')[-1]
    code = myms.split('.')[-2][-3:]
    return code


def get_target_code(targetname):

    # Last three digits of the target name

    code = targetname.replace('-','_').replace('.','p').replace(' ','')[-3:]
    return code


def scrub_target_name(targetname):
#    scrubbed = targetname.replace('+','p').replace(' ','_')
    scrubbed = targetname.replace(' ','_')
    return scrubbed


def make_executable(infile):

    # https://stackoverflow.com/questions/12791997/how-do-you-do-a-simple-chmod-x-from-within-python

    mode = os.stat(infile).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(infile, mode)


def is_odd(xx):
    if (xx % 2) == 0:
        odd = False
    else:
        odd = True
    return odd


def job_handler(syscall,
                jobname,
                infrastructure,
                dependency = None,
                slurm_config = cfg.SLURM_DEFAULTS,
                slurm_account = cfg.SLURM_ACCOUNT,
                slurm_reservation = cfg.SLURM_RESERVATION,
                slurm_nodelist = cfg.SLURM_NODELIST,
                slurm_exclude = cfg.SLURM_EXCLUDE,
                pbs_config = cfg.PBS_DEFAULTS,
                glam_config = cfg.GLAM_STANDARD,
                grouping = None,
                bind = cfg.BIND):
                # slurm_time=cfg.SLURM_TIME,
                # slurm_partition=cfg.SLURM_PARTITION,
                # slurm_ntasks=cfg.SLURM_NTASKS,
                # slurm_nodes=cfg.SLURM_NODES,
                # slurm_cpus=cfg.SLURM_CPUS,
                # slurm_mem=cfg.SLURM_MEM,
                # pbs_program=cfg.PBS_PROGRAM,
                # pbs_walltime=cfg.PBS_WALLTIME,
                # pbs_queue=cfg.PBS_QUEUE,
                # pbs_nodes=cfg.PBS_NODES,
                # pbs_ppn=cfg.PBS_PPN,
                # pbs_mem=cfg.PBS_MEM):


    if infrastructure == 'idia' or infrastructure == 'hippo':


        slurm_time = slurm_config['TIME']
        slurm_partition = slurm_config['PARTITION']
        slurm_ntasks = slurm_config['NTASKS']
        slurm_nodes = slurm_config['NODES']
        slurm_cpus = slurm_config['CPUS']
        slurm_mem = slurm_config['MEM']
        
        # HACK: Override idia settings if hippo here
        # (really this should be broken down as slurm vs. non-slurm scheduler
        if infrastructure == 'hippo':
            if int(slurm_cpus) > 20:
                slurm_cpus='20'
            if int(slurm_cpus) < 20:
                slurm_mem = '60000'
            else:
                slurm_mem = '64000'
            slurm_partition = 'debug'

        slurm_runfile = cfg.SCRIPTS+'/slurm_'+jobname+'.sh'
        slurm_logfile = cfg.LOGS+'/slurm_'+jobname+'.log'

        run_command = jobname+"=`sbatch "
        if dependency:
            #run_command += "-d afterok:${"+dependency+"} "
            run_command += '-d afterok:'+'${'+dependency.replace(':','}:${')+'} '
        run_command += slurm_runfile+" | awk '{print $4}'`"

        if cfg.SLURM_NODELIST != '':
            slurm_nodelist = '#SBATCH --nodelist='+cfg.SLURM_NODELIST+'\n'
        else:
            slurm_nodelist = ''

        if cfg.SLURM_EXCLUDE != '':
            slurm_exclude = '#SBATCH --exclude='+cfg.SLURM_EXCLUDE+'\n'
        else:
            slurm_exclude = ''

        if cfg.SLURM_ACCOUNT != '':
            slurm_account = '#SBATCH --account='+cfg.SLURM_ACCOUNT+'\n'
        else:
            slurm_account = ''

        if cfg.SLURM_RESERVATION != '':
            slurm_reservation = '#SBATCH --reservation='+cfg.SLURM_RESERVATION+'\n'
        else:
            slurm_reservation = ''

        f = open(slurm_runfile,'w')
        f.writelines(['#!/bin/bash\n',
            '#file: '+slurm_runfile+':\n',
            '#SBATCH --job-name='+jobname+'\n',
            '#SBATCH --time='+slurm_time+'\n',
            '#SBATCH --partition='+slurm_partition+'\n'
            '#SBATCH --ntasks='+slurm_ntasks+'\n',
            '#SBATCH --nodes='+slurm_nodes+'\n',
            '#SBATCH --cpus-per-task='+slurm_cpus+'\n',
            '#SBATCH --mem='+slurm_mem+'\n',
            '#SBATCH --output='+slurm_logfile+'\n',
            slurm_nodelist,
            slurm_exclude,
            slurm_account,
            slurm_reservation,
            'SECONDS=0\n',
            syscall+'\n',
            'echo "****ELAPSED "$SECONDS" '+jobname+'"\n'])
#            'sleep 10\n'])
        f.close()

        make_executable(slurm_runfile)

    elif infrastructure == 'chpc':

        pbs_program = pbs_config['PROGRAM']
        pbs_walltime = pbs_config['WALLTIME']
        pbs_queue = pbs_config['QUEUE']
        pbs_nodes = pbs_config['NODES']
        pbs_ppn = pbs_config['PPN']
        pbs_mem = pbs_config['MEM']

        pbs_runfile = cfg.SCRIPTS+'/pbs_'+jobname+'.sh'
        pbs_logfile = cfg.LOGS+'/pbs_'+jobname+'.log'
        pbs_errfile = cfg.LOGS+'/pbs_'+jobname+'.err'

        run_command = jobname+"=`qsub "
        if dependency:
          run_command += "-W depend=afterok:${"+dependency+"} "
        run_command += pbs_runfile+" | awk '{print $1}'`"

        f = open(pbs_runfile,'w')
        f.writelines(['#!/bin/bash\n',
            '#PBS -N '+jobname+'\n',
            '#PBS -P '+pbs_program+'\n'
            '#PBS -l walltime='+pbs_walltime+'\n',
            '#PBS -l nodes='+pbs_nodes+':ppn='+pbs_ppn+',mem='+pbs_mem+'\n',
            '#PBS -q '+pbs_queue+'\n'
            '#PBS -o '+pbs_logfile+'\n'
            '#PBS -e '+pbs_errfile+'\n'
            'SECONDS=0\n'
            'module load chpc/singularity\n'
            'cd '+cfg.CWD+'\n',
            syscall+'\n',
            'echo "****ELAPSED "$SECONDS" "'+jobname+'"\n',
            'sleep 10\n'])
        f.close()

        make_executable(pbs_runfile)

    elif infrastructure == 'glam':

        # GLAM uses custom addqueue scripts
        glam_cpus = glam_config['CPUS']
        glam_mem = glam_config['MEM']

        # Split syscall into lines
        syscall_lines = [line for line in syscall.splitlines() if line.strip()]
        
        # Group commands based on grouping parameter
        if grouping is None:
            # Put all commands in a single file
            command_groups = [syscall_lines]
        elif grouping == 1:
            # One command per file (original behavior)
            command_groups = [[line] for line in syscall_lines]
        else:
            # Group N commands per file
            command_groups = []
            for i in range(0, len(syscall_lines), grouping):
                command_groups.append(syscall_lines[i:i+grouping])
        
        # Create runfiles for each group
        for k, group in enumerate(command_groups):
            glam_runfile = cfg.SCRIPTS+'/glam_'+jobname+f'_part{k:03d}.sh'
            
            f = open(glam_runfile,'w')
            f.write('#!/bin/bash\n')
            for line in group:
                f.write(line+'\n')
            f.close()
            
            make_executable(glam_runfile)

        # Return empty run_command as GLAM uses custom addqueue scripts
        run_command = ''

    elif infrastructure == 'node':

        node_logfile = cfg.LOGS + '/oxk_' + jobname + '.log'

        lines = []
        first = True

        for line in syscall.splitlines():
            if line.strip():  # non-blank
                if first:
                    lines.append(f"{{ time {line} ; }} 2>&1 | tee {node_logfile}")
                    first = False
                else:
                    lines.append(f"{{ time {line} ; }} 2>&1 | tee -a {node_logfile}")
            else:
                lines.append("")

        run_command = "\n".join(lines)

    run_command += '\n'

    return run_command


def mem_string_to_gb(mem):
    headroom = 0.98 # fraction of memory specified in IDIA/CHPC config to convert to absmem (hippo a special case)
    mem = mem.upper().replace('B','')
    if 'M' in mem:
        factor = 1e-3
    if 'G' in mem:
        factor = 1
    if 'T' in mem:
        factor = 1e3
    absmem = float(''.join(x for x in mem if x.isdigit()))
    absmem = int(absmem * factor * headroom)
    return absmem


def absmem_helper(step,infrastructure,absmem):
    if infrastructure == 'chpc':
        config_mem = step['pbs_config']['MEM']
    elif infrastructure == 'idia':
        config_mem = step['slurm_config']['MEM']
    elif infrastructure == 'glam':
        config_mem = step['glam_config']['MEM']
    elif infrastructure == 'hippo':
        slurm_cpus = step['slurm_config']['CPUS']
        if int(slurm_cpus) > 20:
            slurm_cpus='20'
        if int(slurm_cpus) < 20:
            config_mem = '60gb'
        else:
            config_mem = '64gb'
    if infrastructure != 'node':
        absmem = mem_string_to_gb(config_mem)
    return absmem


def get_scan_times(scanpickle):
    scan_times = []
    ss = pickle.load(open(scanpickle,'rb'))
    fields = []
    for ii in ss:
        fields.append(ii[1])
    fields = numpy.unique(fields).tolist()
    for field in fields:
        scans = []
        intervals = []
        for ii in ss:
            if ii[1] == field:
                scans.append(ii[0])
                intervals.append(ii[5])
        scan_times.append((field,scans,intervals))
    return scan_times



def generate_syscall_casa(casascript,casalogfile='',extra_args=''):

    if casalogfile != '':
        syscall = 'casa --logfile {} '.format(casalogfile)
    else:
        syscall = 'casa --log2term '
    syscall += '--nogui '
    if extra_args != '':
        syscall += extra_args
    syscall += '-c ' + casascript

#    syscall += '\n'

    return syscall

def generate_syscall_casa_short(casascript):

    syscall = 'casa --nogui -c ' + casascript
    return syscall


def generate_syscall_cubical(parset,myms,extra_args=''):

    # now = timenow()
    # outname = 'cube_'+prefix+'_'+myms.split('/')[-1]+'_'+now

    # # Debugging stuff
    # syscall = 'bash -c "/sbin/sysctl vm.max_map_count ; '
    # syscall += 'df -h /dev/shm ; '

    # syscall += 'gocubical '+parset+' '
    # syscall += '--data-ms='+myms+' '
    # syscall += '--out-name='+outname

    # # Move output to logs...
    # syscall += ' ; mv '+outname+'* '+LOGS+'"'

    syscall = 'gocubical '+parset+' '
    syscall += '--data-ms='+myms+' '
    if extra_args != '':
        syscall += extra_args
#        syscall += '--out-name '+outname+' '

    return syscall

def generate_syscall_quartical(yaml,myms,extra_args=''):

    syscall = 'goquartical '+yaml+' '
    syscall += 'input_ms.path='+myms+' '
    if extra_args != '':
        syscall += extra_args

    return syscall

def generate_syscall_tricolour(myms = '',
                          config = '',
                          datacol = 'DATA',
                          subtractcol = '',
                          fields = 'all',
                          strategy = 'polarisation'):

    syscall = 'tricolour '
    syscall += '--config '+config+' '
    syscall += '--data-column '+datacol+' '
    if subtractcol != '':
        syscall += '--subtract-model-colum '+subtractcol+' '
    if fields != 'all':
        syscall += '--field-names '+fields+' '
    syscall += '--flagging-strategy '+strategy+' '
    syscall += myms

    return syscall


def generate_syscall_predict(msname,
                            imgname,
                            field = cfg.WSC_FIELD,
                            pol = cfg.WSC_POL,
                            intervalsout = cfg.WSC_INTERVALSOUT,
                            nwlayersfactor = cfg.WSC_NWLAYERSFACTOR,
                            chanout = cfg.WSC_DMASK_CHANNELSOUT,
                            usewgridder = cfg.WSC_USEWGRIDDER,
#                            imsize = cfg.WSC_IMSIZE,
#                            cellsize = cfg.WSC_CELLSIZE,
#                            predictchannels = cfg.WSC_PREDICTCHANNELS,
                            mem = cfg.WSC_MEM,
                            absmem = cfg.WSC_ABSMEM,
                            parallelreordering = cfg.WSC_PARALLELREORDERING,
                            parallelgridding = cfg.WSC_PARALLELGRIDDING,
                            cores = False,
                            tempdir = None):

    # Generate system call to run wsclean in predict mode
    syscall = 'wsclean '
    syscall += '-predict '
    if cores:
        syscall += '-j '+str(cores)+' '
    syscall += '-field '+str(field)+' '
    syscall += '-pol ' + str(pol) + ' '

    if parallelreordering != 0:
        syscall += '-parallel-reordering '+str(parallelreordering)+' '
    if parallelgridding != 0:
        syscall += '-parallel-gridding '+str(parallelgridding)+' '
    if tempdir is not None:
        syscall += '-temp-dir '+tempdir+' '
    if intervalsout:
        syscall += '-intervals-out ' + str(intervalsout) + ' '
    syscall += '-channels-out '+str(chanout)+' '
    syscall += '-name '+imgname+' '
    if absmem < 0:
        syscall += '-mem '+str(mem)+' '
    else:
        syscall += '-abs-mem '+str(absmem)+' '
#    syscall += '-predict-channels '+str(predictchannels)+' '
    syscall += msname + ' '

    return syscall

def generate_syscall_wsclean(mslist,
                          imgname,
                          datacol,
                          continueclean = cfg.WSC_CONTINUE,
                          field = cfg.WSC_FIELD,
                          makepsf = cfg.WSC_MAKEPSF,
                          nodirty = cfg.WSC_NODIRTY,
                          startchan = cfg.WSC_STARTCHAN,
                          endchan = cfg.WSC_ENDCHAN,
                          minuvl = cfg.WSC_MINUVL,
                          maxuvl = cfg.WSC_MAXUVL,
                          even = cfg.WSC_EVEN,
                          odd = cfg.WSC_ODD,
                          chanout = cfg.WSC_PCAL_CHANNELSOUT,
                          maxchan = cfg.WSC_MAX_CHANNELS,
                          chandeconvolution = cfg.WSC_CHANDECONV,
                          interval0 = cfg.WSC_INTERVAL0,
                          interval1 = cfg.WSC_INTERVAL1,
                          intervalsout = cfg.WSC_INTERVALSOUT,
                          imsize = cfg.WSC_IMSIZE,
                          cellsize = cfg.WSC_CELLSIZE,
                          weight = cfg.WSC_WEIGHT,
                          tapergaussian = cfg.WSC_TAPERGAUSSIAN,
                          niter = cfg.WSC_NITER,
                          gain = cfg.WSC_GAIN,
                          mgain = cfg.WSC_MGAIN,
                          multiscale = cfg.WSC_MULTISCALE,
                          multiscale_bias = cfg.WSC_MULTISCALE_BIAS,
                          scales = cfg.WSC_SCALES,
                          nonegative = cfg.WSC_NONEGATIVE,
                          sourcelist = cfg.WSC_SOURCELIST,
                          bda = cfg.WSC_BDA,
                          bdafactor = cfg.WSC_BDAFACTOR,
                          nwlayersfactor = cfg.WSC_NWLAYERSFACTOR,
                          joinchannels = cfg.WSC_JOINCHANNELS,
                          joinpolarizations = cfg.WSC_JOINPOLARIZATIONS,
                          squarechans = cfg.WSC_SQUARECHANS,
                          qu_automask_scale = cfg.WSC_QU_AUTOMASK_SCALE,
                          pol = cfg.WSC_POL,
                          splitpol = cfg.WSC_SPLITPOL,
                          padding = cfg.WSC_PADDING,
                          nomodel = cfg.WSC_NOMODEL,
                          mask = cfg.WSC_MASK,
                          mfweight = cfg.WSC_MFWEIGHT,
                          threshold = cfg.WSC_THRESHOLD,
                          autothreshold = cfg.WSC_AUTOTHRESHOLD,
                          automask = cfg.WSC_AUTOMASK,
                          localrms = cfg.WSC_LOCALRMS,
                          localrms_strength = cfg.WSC_LOCALRMS_STRENGTH,
                          stopnegative = cfg.WSC_STOPNEGATIVE,
                          fitspectralpol = cfg.WSC_FITSPECTRALPOL,
                          circularbeam = cfg.WSC_CIRCULARBEAM,
                          mem = cfg.WSC_MEM,
                          absmem = cfg.WSC_ABSMEM,
                          usewgridder = cfg.WSC_USEWGRIDDER,
                          wgridderaccuracy = cfg.WSC_WGRIDDERACCURACY,
                          useidg = cfg.WSC_USEIDG,
                          idgmode = cfg.WSC_IDGMODE,
                          tukeytaper=cfg.WSC_TUKEYTAPER,
                          paralleldeconvolution = cfg.WSC_PARALLELDECONVOLUTION,
                          parallelreordering = cfg.WSC_PARALLELREORDERING,
                          parallelgridding = cfg.WSC_PARALLELGRIDDING,
                          cores = False,
                          tempdir = None):

    # Generate system call to run wsclean based imaging for 2GC (and beyond)
    if  imsize % 2 != 0:
        print(col('wsclean')+'Do not use odd image sizes')
        sys.exit()

    if continueclean and bda:
        print(col('wsclean')+'Cannot continue deconvolution if BDA is enabled')
        sys.exit()

    if even and odd:
        print(col('wsclean')+'Even and odd timeslots selections are both enabled, defaulting to all.')
        even = False
        odd = False

    # -----------
    syscall = 'wsclean '
    syscall += '-log-time '
    if cores:
        syscall += '-j '+str(cores)+' '
    if absmem < 0:
        syscall += '-mem '+str(mem)+' '
    else:
        syscall += '-abs-mem '+str(absmem)+' '
    if continueclean:
        syscall += '-continue '
    if parallelreordering != 0:
        syscall += '-parallel-reordering '+str(parallelreordering)+' '
    if tempdir is not None:
        syscall += '-temp-dir '+tempdir+' '

    # Outputs  
    if makepsf:
        syscall += '-make-psf '
    if nodirty:
        syscall += '-no-dirty '
    if sourcelist: # and fitspectralpol != 0:
        syscall += '-save-source-list '


    # Data selection
    syscall += '-data-column '+datacol+' '
    syscall += '-field '+str(field)+' '
    if startchan != -1 and endchan != -1:
        syscall += '-channel-range '+str(startchan)+' '+str(endchan)+' '
    if minuvl != '':
        syscall += '-minuv-l '+str(minuvl)+' '
    if maxuvl != '':
        syscall += '-maxuv-l '+str(maxuvl)+' '
    if even:
        syscall += '-even-timesteps '
    if odd:
        syscall += '-odd-timesteps '
    if interval0 or interval1:
        syscall += '-interval '+str(interval0)+' '+str(interval1)+' '
    if intervalsout:
        syscall += '-intervals-out '+str(intervalsout)+' '

    # Image dimensions
    syscall += '-size '+str(imsize)+' '+str(imsize)+' '
    syscall += '-scale '+cellsize+' '

    # Gridding
    if usewgridder:
        syscall += '-gridder wgridder '
#        syscall += '-wgridder-accuracy '+str(wgridderaccuracy)+' '
    if parallelgridding != 0:
        syscall += '-parallel-gridding '+str(parallelgridding)+' '
    if bda and not useidg:
        syscall += '-baseline-averaging '+str(bdafactor)+' '
        syscall += '-no-update-model-required '
    elif not bda and nomodel:
        syscall += '-no-update-model-required '
    if not usewgridder and not useidg:
        syscall += '-padding '+str(padding)+' '
        syscall += '-nwlayers-factor '+str(nwlayersfactor)+' '
    if useidg:
        syscall += '-use-idg '
        syscall += '-idg-mode '+idgmode+' '

    # Weighting
    syscall += '-weight '+str(weight)+' '
    if tapergaussian != '':
        syscall += '-taper-gaussian '+str(tapergaussian)+' '
    if mfweight:
        syscall  += '-mf-weighting '
    else:
        syscall += '-no-mf-weighting '
    if tukeytaper:
        if minuvl == '':
            syscall += f'-minuv-l 0.0 -taper-inner-tukey {tukeytaper} '
        else:
            syscall += f'-taper-inner-tukey {tukeytaper} '

    # Deconvolution
    if paralleldeconvolution != 0:
        syscall += '-parallel-deconvolution '+str(paralleldeconvolution)+' '    
    if multiscale:
        syscall += '-multiscale '
        syscall += f'-multiscale-scale-bias {multiscale_bias} '
        # syscall += '-multiscale-scales '+scales+' ' # WSCLEAN docs seem to not favour this option for multiscale cleaning anymore
    syscall += '-niter '+str(niter)+' '
    syscall += '-gain '+str(gain)+' '
    syscall += '-mgain '+str(mgain)+' '
    if chandeconvolution and fitspectralpol != 0:
        syscall += '--deconvolution-channels '+str(chandeconvolution)+' '
    if joinchannels:
        syscall += '-join-channels '
        if squarechans and pol != 'I':
            syscall += '-squared-channel-joining '

    if nonegative:
        syscall += '-no-negative '
    if stopnegative:
        syscall += '-stop-negative '
    if circularbeam:
        syscall += '-circular-beam '

    # Masking
    if mask:
        if mask.lower() == 'fits':
            mymask = glob.glob(cfg.IMAGES + '/*mask.fits')[0]
            syscall += '-fits-mask '+mymask+' '
        else:
            syscall += '-fits-mask '+mask+' '
    if automask:
        syscall += '-auto-mask '+str(automask)+' '
    if autothreshold:
        syscall += '-auto-threshold '+str(autothreshold)+' '
        if localrms:
            syscall += '-local-rms '
            syscall += '-local-rms-strength '+str(localrms_strength)+' '
    if threshold:
        syscall += '-threshold '+str(threshold)+' '


    # Andrew additions for more complex imaging strategies
    wsclean_syscall_base = syscall # Copy syscall options to apply more advanced techniques
    syscall_arr = []
    _syscall_names = []   # tracks the -name value for each entry in syscall_arr

    # Add intervals out option to wsclean call

    # Fraser suggestions to break up chanels for high-mem costly imaging
    if maxchan < chanout:
        nchans = cfg.PRE_NCHANS
        intchans = int(nchans / chanout * maxchan) # Each image will be composed of this many channels
        nint = int(chanout / maxchan)

        if len(mslist) > 1:
            print(col('WARNING') + f'You are using multiple MS files; make sure they all have {nchans} channels')

        # Make sure all the integers are dividable
        if nchans % maxchan != 0 or nchans % chanout:
            print(col('ERROR') + f'MS file has {nchans} channels, which is not divisble by {maxchan} or {chanout}; please choose better numbers!')
            sys.exit()

        if chanout % maxchan != 0:
            print(col('ERROR') + 'Maxchan is not an integer multiple of the total output channels; please choose better numbers!')
            sys.exit()

        for k in range(nint):
            _part_name = f'{imgname}_part{k:04d}'
            syscall_arr.append(wsclean_syscall_base + f'-channels-out {maxchan} -channel-range {k * intchans} {(k + 1) * intchans} -name {_part_name} ')
            _syscall_names.append(_part_name)

    else:
        syscall_arr.append(wsclean_syscall_base + f'-channels-out {chanout} -name {imgname} ')
        _syscall_names.append(imgname)

    # Add option to split the deconvolution into IQU & IV steps
    # Sources with large rotation measures may not want polynomial fitting to the QU channels
    if splitpol and pol != 'I':
        pol_QU  = pol.replace('V','').replace('I','')
        pol_IV  = pol.replace('Q', '').replace('U','')

        spectralpol_IV = ''
        if fitspectralpol != 0:
            spectralpol_IV = '-fit-spectral-pol '+str(fitspectralpol) + ' '

        joinpol_QU = ''

        if joinpolarizations and len(pol_QU) >= 2:
            joinpol_QU = '-join-polarizations '

        # When splitpol is active, join-polarizations applies to QU only — not IV
        joinpol_IV = ''

        k = len(syscall_arr)
        syscall_arr    += syscall_arr
        _syscall_names += _syscall_names

        # Resolve the QU/IV split deconvolution behaviour from qu_automask_scale:
        #   0       : disabled — no changes to QU or IV, local-rms not added
        #   'auto'  : QU gets local-rms; IV auto-mask floored at 3.0, auto-threshold floored at 1.0
        #   float   : QU gets local-rms; IV auto-mask multiplied by scale, auto-threshold floored at 1.0
        try:
            _scale = float(qu_automask_scale)
            _qu_disabled = (_scale == 0)
            _qu_auto_mode = False
        except (TypeError, ValueError):
            _scale = None
            _qu_disabled = False
            _qu_auto_mode = True  # 'auto' or unrecognised string

        for _k in range(k):
            # QU: strip --deconvolution-channels (only appropriate for IV)
            qu_syscall = syscall_arr[_k]
            if chandeconvolution:
                qu_syscall = qu_syscall.replace(
                    f'--deconvolution-channels {chandeconvolution} ', '')
            # Enable local-rms for QU unless disabled
            if not _qu_disabled and not localrms:
                qu_syscall = qu_syscall.replace(
                    '-auto-threshold ',
                    '-local-rms -local-rms-strength '+str(localrms_strength)+' -auto-threshold ')
            syscall_arr[_k] = qu_syscall + f'-pol {pol_QU} {joinpol_QU} '
            # IV: strip squared-channel-joining unconditionally (not appropriate for Stokes V)
            iv_syscall = syscall_arr[_k + k].replace('-squared-channel-joining ', '')
            # Apply IV mask/threshold overrides unless disabled
            if not _qu_disabled:
                if _qu_auto_mode:
                    _iv_automask = max(automask, 3.0) if automask else 3.0
                else:
                    _iv_automask = round(automask * _scale, 1) if automask else automask
                if automask and _iv_automask != automask:
                    iv_syscall = iv_syscall.replace(
                        f'-auto-mask {automask} ',
                        f'-auto-mask {_iv_automask} ')
                _iv_autothreshold = max(autothreshold, 1.0) if autothreshold else 1.0
                if autothreshold and _iv_autothreshold != autothreshold:
                    iv_syscall = iv_syscall.replace(
                        f'-auto-threshold {autothreshold} ',
                        f'-auto-threshold {_iv_autothreshold} ')
            syscall_arr[_k + k] = iv_syscall + f'-pol {pol_IV} {spectralpol_IV} {joinpol_IV} '

    else:
        joinpol = ''
        spectralpol = ''

        if joinpolarizations and len(pol) > 1:
            joinpol     = '-join-polarizations '

        if fitspectralpol != 0:
            spectralpol = '-fit-spectral-pol '+str(fitspectralpol) + ' '   

        k = len(syscall_arr)
        for _k in range(k):     
            syscall_arr[_k] += f'-pol {pol} {spectralpol} {joinpol}'   

    # End by appending ms files to the wsclean calls
    for k in range(len(syscall_arr)):
        for myms in mslist:
            syscall_arr[k] += myms + ' '

    # In splitpol mode, if a pol group has only one Stokes parameter WSClean drops
    # the Stokes label from output filenames.  Insert mv commands right after each
    # affected wsclean call to reinstate it.
    if splitpol and pol != 'I':
        _qu_single = len(pol_QU) == 1
        _iv_single = len(pol_IV) == 1

        if _qu_single or _iv_single:
            def _rename_cmd(name, stokes):
                parts = []
                for suf in ('dirty', 'image', 'model', 'residual'):
                    parts.append(
                        f'for f in {name}-*-{suf}.fits; do '
                        f'[ -f "$f" ] && [[ ! "$f" =~ -[IQUV]-{suf}\\.fits$ ]] && mv "$f" "${{f%-{suf}.fits}}-{stokes}-{suf}.fits"; '
                        f'done'
                    )
                inner = ' && '.join(parts)
                return f"bash -c '{inner}'"

            _orig_k   = len(syscall_arr) // 2
            _renames  = {}   # index in syscall_arr -> rename shell command
            for _j in range(_orig_k):
                if _qu_single:
                    _renames[_j] = _rename_cmd(_syscall_names[_j], pol_QU)
                if _iv_single:
                    _renames[_j + _orig_k] = _rename_cmd(_syscall_names[_j + _orig_k], pol_IV)

            # Interleave: insert rename command immediately after each wsclean call
            final_arr = []
            for _idx, _call in enumerate(syscall_arr):
                final_arr.append(_call)
                if _idx in _renames:
                    final_arr.append(_renames[_idx])
            syscall_arr = final_arr

    return syscall_arr

def generate_syscall_makemask(restoredimage,
                            outfile = '',
                            thresh = cfg.MAKEMASK_THRESH,
                            boxsize = cfg.MAKEMASK_BOXSIZE,
                            smallbox = cfg.MAKEMASK_SMALLBOX,
                            islandsize = cfg.MAKEMASK_ISLANDSIZE,
                            dilation = cfg.MAKEMASK_DILATION,
                            zoompix = cfg.DDF_NPIX):

    # Generate call to MakeMask.py and dilate the result
  
    if outfile == '':
        outfile = restoredimage.replace('.fits','.mask.fits')

    syscall = 'bash -c "'
    syscall += 'python3 '+cfg.TOOLS+'/pyMakeMask.py '
    syscall += '--dilate='+str(dilation)+' '
    syscall += '--boxsize='+str(boxsize)+' '
    syscall += '--smallbox='+str(smallbox)+' '
    syscall += '--islandsize='+str(islandsize)+' '
    syscall += '--threshold='+str(thresh)+' '
    syscall += '--outfile='+str(outfile)+' '
    syscall += restoredimage

    if zoompix != '':
        zoomfits = outfile.replace('.fits','.zoom'+str(zoompix)+'.fits')
        syscall += ' && fitstool.py -z '+str(zoompix)+' -o '+zoomfits+' '
        syscall += outfile
    syscall += '"'

    return syscall,outfile

def generate_syscall_breizorro(restoredimage,
                            outfile = '',
                            thresh = cfg.BREIZORRO_THRESH,
                            boxsize = cfg.BREIZORRO_BOXSIZE,
                            fillholes = cfg.BREIZORRO_FILLHOLES,
                            dilation = cfg.BREIZORRO_DILATION):

    # Generate call to MakeMask.py and dilate the result
  
    if outfile == '':
        outfile = restoredimage.replace('.fits','.mask.fits')

    syscall = 'breizorro '
    if dilation:
        syscall += '--dilate='+str(dilation)+' '
    syscall += '--boxsize='+str(boxsize)+' '
    syscall += '--threshold='+str(thresh)+' '
    syscall += '--outfile='+str(outfile)+' '
    if fillholes:
        syscall += '--fill-holes '
    syscall += '--restored-image='+str(restoredimage)

    return syscall,outfile

def generate_syscall_ddfacet(mspattern,
                          imgname,
                          ddid = cfg.DDF_DDID,
                          field = cfg.DDF_FIELD,
                          colname = cfg.DDF_COLNAME,
                          chunkhours = cfg.DDF_CHUNKHOURS,
                          datasort = cfg.DDF_DATASORT,
                          predictcolname = cfg.DDF_PREDICTCOLNAME,
                          initdicomodel = cfg.DDF_INITDICOMODEL,
                          outputalso = cfg.DDF_OUTPUTALSO,
                          outputimages = cfg.DDF_OUTPUTIMAGES,
                          outputcubes = cfg.DDF_OUTPUTCUBES,
                          npix = cfg.DDF_NPIX,
                          cell = cfg.DDF_CELL,
                          diammax = cfg.DDF_DIAMMAX,
                          diammin = cfg.DDF_DIAMMIN,
                          nfacets =cfg.DDF_NFACETS,
                          psfoversize = cfg.DDF_PSFOVERSIZE,
                          padding = cfg.DDF_PADDING,
                          robust = cfg.DDF_ROBUST,
                          sparsification = cfg.DDF_SPARSIFICATION,
                          ncpu = cfg.DDF_NCPU,
                          cachereset = cfg.DDF_CACHERESET,
                          cachedir = cfg.DDF_CACHEDIR,
                          cachehmp = cfg.DDF_CACHEHMP,
                          beam = cfg.DDF_BEAM,
                          beamnband = cfg.DDF_BEAMNBAND,
                          dtbeammin = cfg.DDF_DTBEAMMIN,
                          fitsparangleincdeg = cfg.DDF_FITSPARANGLEINCDEG,
                          beamcentrenorm = cfg.DDF_BEAMCENTRENORM,
                          beamsmooth = cfg.DDF_BEAMSMOOTH,
                          feedswap = cfg.DDF_FEEDSWAP,
                          nband = cfg.DDF_NBAND,
                          ndegridband = cfg.DDF_NDEGRIDBAND,
                          ddsols = cfg.DDF_DDSOLS,
                          ddmodegrid = cfg.DDF_DDMODEGRID,
                          ddmodedegrid = cfg.DDF_DDMODEDEGRID,
                          gain = cfg.DDF_GAIN,
                          fluxthreshold = cfg.DDF_FLUXTHRESHOLD,
                          cyclefactor = cfg.DDF_CYCLEFACTOR,
                          rmsfactor = cfg.DDF_RMSFACTOR,
                          deconvmode = cfg.DDF_DECONVMODE,
                          ssd_deconvpeakfactor = cfg.DDF_SSD_DECONVPEAKFACTOR,
                          ssd_maxminoriter = cfg.DDF_SSD_MAXMINORITER,
                          ssd_maxmajoriter = cfg.DDF_SSD_MAXMAJORITER,
                          ssd_enlargedata = cfg.DDF_SSD_ENLARGEDATA,
                          hogbom_deconvpeakfactor = cfg.DDF_HOGBOM_DECONVPEAKFACTOR,
                          hogbom_maxminoriter = cfg.DDF_HOGBOM_MAXMINORITER,
                          hogbom_maxmajoriter = cfg.DDF_HOGBOM_MAXMAJORITER,
                          hogbom_polyfitorder = cfg.DDF_HOGBOM_POLYFITORDER,
                          mask = cfg.DDF_MASK,
                          masksigma = cfg.DDF_MASKSIGMA,
                          conservememory = cfg.DDF_CONSERVEMEMORY):

    syscall = 'DDF.py '
    # [Data]
    syscall += '--Data-MS '+mspattern+'//'+ddid+'//'+field+' '
    syscall += '--Data-ColName '+colname+' '
    syscall += '--Data-ChunkHours '+str(chunkhours)+' '
    syscall += '--Data-Sort '+str(datasort)+' '
    # [Predict]
    if predictcolname:
        syscall += '--Predict-ColName '+predictcolname+' '
    if initdicomodel != '':
        syscall += '--Predict-InitDicoModel '+initdicomodel+' '
    # [Output]
    syscall += '--Output-Name '+imgname+' '
    syscall += '--Output-Mode Clean '
    syscall += '--Output-Also '+outputalso+' '
    syscall += '--Output-Images '+outputimages+' '
    if outputcubes != '':
        syscall += '--Output-Cubes '+outputcubes+' '
    # [Image]
    syscall += '--Image-NPix '+str(npix)+' '
    syscall += '--Image-Cell '+str(cell)+' '
    # [Facets]
    syscall += '--Facets-DiamMax '+str(diammax)+' '
    syscall += '--Facets-DiamMin '+str(diammin)+' '
    syscall += '--Facets-NFacets '+str(nfacets)+' '
    syscall += '--Facets-PSFOversize '+str(psfoversize)+' '
    syscall += '--Facets-Padding '+str(padding)+' '
    # [Weight]
    syscall += '--Weight-Mode Briggs '
    syscall += '--Weight-Robust '+str(robust)+' '
    # [CF]
    # syscall += '--CF-wmax 0 '
    # syscall += '--CF-Nw 100 '
    # [Comp]
    syscall += '--Comp-GridDecorr 0.01 '
    syscall += '--Comp-DegridDecorr 0.01 '
    syscall += '--Comp-Sparsification '+str(sparsification)+' '
    # [Parallel]
    syscall += '--Parallel-NCPU '+str(ncpu)+' '
    # [Cache]    
    syscall += '--Cache-Reset '+str(cachereset)+' '
    syscall += '--Cache-Dir '+str(cachedir)+' '
    syscall += '--Cache-HMP '+str(cachehmp)+' '
    # [Beam]
    if beam == '':
        syscall += '--Beam-Model None '
    else:
        syscall += '--Beam-Model FITS '
        syscall += "--Beam-FITSFile \'"+str(beam)+"\' "
        syscall += '--Beam-NBand '+str(beamnband)+' '
        syscall += '--Beam-DtBeamMin '+str(dtbeammin)+' '
        syscall += '--Beam-FITSParAngleIncDeg '+str(fitsparangleincdeg)+' '
        syscall += '--Beam-CenterNorm '+str(beamcentrenorm)+' '
        syscall += '--Beam-FITSFeedSwap '+str(feedswap)+' '
        syscall += '--Beam-Smooth '+str(beamsmooth)+' '
    # [Freq]
    syscall += '--Freq-NBand '+str(nband)+' '
    syscall += '--Freq-NDegridBand '+str(ndegridband)+' '
    # [DDESolutions]
    if ddsols != '':
        syscall += '--DDESolutions-DDSols '+ddsols+' '
        syscall += '--DDESolutions-DDModeGrid '+ddmodegrid+' '
        syscall += '--DDESolutions-DDModeDeGrid '+ddmodedegrid+' '
    # [Deconv]
    syscall += '--Deconv-Gain '+str(gain)+' '
    syscall += '--Deconv-FluxThreshold '+str(fluxthreshold)+' '
    syscall += '--Deconv-CycleFactor '+str(cyclefactor)+' '
    syscall += '--Deconv-RMSFactor '+str(rmsfactor)+' '
    if deconvmode.lower() == 'ssd':
        syscall += '--Deconv-Mode SSD '
        syscall += '--Deconv-PeakFactor '+str(ssd_deconvpeakfactor)+' '
        syscall += '--Deconv-MaxMajorIter '+str(ssd_maxmajoriter)+' '
        syscall += '--Deconv-MaxMinorIter '+str(ssd_maxminoriter)+' '
        syscall += '--SSDClean-NEnlargeData '+str(ssd_enlargedata)+' '
    elif deconvmode.lower() == 'hogbom':
        syscall += '--Deconv-Mode Hogbom '
        syscall += '--Deconv-PeakFactor '+str(hogbom_deconvpeakfactor)+' '
        syscall += '--Deconv-MaxMajorIter '+str(hogbom_maxmajoriter)+' '
        syscall += '--Deconv-MaxMinorIter '+str(hogbom_maxminoriter)+' '    
        syscall += '--Hogbom-PolyFitOrder '+str(hogbom_polyfitorder)+' '
    # [Mask]
    if mask.lower() == 'fits':
        mymask = sorted(glob.glob('*mask.fits')[0])
        syscall += '--Mask-Auto 0 '
        syscall += '--Mask-External '+mymask+' '
    elif mask.lower() == 'auto':
        syscall += '--Mask-Auto 1 '
        syscall += '--Mask-SigTh '+str(masksigma)+' '
    else:
        syscall += '--Mask-Auto 0 '
        syscall += '--Mask-External '+mask+' '
    # [Misc]
    syscall += '--Misc-ConserveMemory '+str(conservememory)+' '
    syscall += '--Log-Memory 1 '
    syscall += '--Log-Boring 1 '

    return syscall


def generate_syscall_killms(myms,
                        baseimg,
                        outsols,
                        nodesfile,
                        dicomodel = cfg.KMS_DICOMODEL,
                        tchunk = cfg.KMS_TCHUNK,
                        incol = cfg.KMS_INCOL,
                        outcol = cfg.KMS_OUTCOL,
                        beam = cfg.KMS_BEAM,
                        beamat = cfg.KMS_BEAMAT,
                        dtbeammin = cfg.KMS_DTBEAMMIN,
                        centrenorm = cfg.KMS_CENTRENORM,
                        nchanbeamperms = cfg.KMS_NCHANBEAMPERMS,
                        fitsparangleincdeg = cfg.KMS_FITSPARANGLEINCDEG,
                        fitsfeedswap = cfg.KMS_FITSFEEDSWAP,
                        maxfacetsize = cfg.KMS_MAXFACETSIZE,
                        uvminmax = cfg.KMS_UVMINMAX,
                        fieldid = cfg.KMS_FIELDID,
                        ddid = cfg.KMS_DDID,
                        ncpu = cfg.KMS_NCPU,
                        dobar = cfg.KMS_DOBAR,
                        debugpdb = cfg.KMS_DEBUGPDB,
                        solvertype= cfg.KMS_SOLVERTYPE,
                        dt = cfg.KMS_DT,
                        nchansols = cfg.KMS_NCHANSOLS,
                        niterkf = cfg.KMS_NITERKF,
                        covq = cfg.KMS_COVQ):

    # Generate system call to run killMS

    syscall = 'kMS.py '
    # [VisData]
    syscall+= '--MSName '+myms+' '
    syscall+= '--TChunk '+str(tchunk)+' '
    syscall+= '--InCol '+incol+' '
    syscall+= '--OutCol '+outcol+' '
    # [Beam]
    if beam == '':
        syscall+= '--BeamModel None '
    else:
        syscall+= '--BeamModel FITS '
        syscall+= '--BeamAt '+beamat+' '
        syscall+= '--DtBeamMin '+str(dtbeammin)+' '
        syscall+= '--CenterNorm '+str(centrenorm)+' '
        syscall+= '--NChanBeamPerMS '+str(nchanbeamperms)+' '
        syscall+= "--FITSFile \'"+str(beam)+"\' "
        syscall+= '--FITSParAngleIncDeg '+str(fitsparangleincdeg)+' '
        syscall+= '--FITSFeedSwap '+str(fitsfeedswap)+' '
    # [ImageSkyModel]
    syscall+= '--BaseImageName '+baseimg+' '
    if dicomodel != '':
        syscall+= '--DicoModel '+dicomodel+' '
    syscall+= '--NodesFile '+nodesfile+' '
    syscall+= '--MaxFacetSize '+str(maxfacetsize)+' '
    # [DataSelection]
    syscall+= '--UVMinMax '+uvminmax+' '
    syscall+= '--FieldID '+str(fieldid)+' '
    syscall+= '--DDID '+str(ddid)+' '
    # [Weighting]
    syscall+= '--Weighting Natural '
    # [Actions]
    syscall+= '--NCPU '+str(ncpu)+' '
    syscall+= '--DoBar '+str(dobar)+' '
#    syscall+= '--DebugPdb '+str(debugpdb)+' '
    # [Solutions]
    syscall+= '--OutSolsName '+outsols+' '
    # [Solvers]
    syscall+= '--SolverType '+solvertype+' '
    syscall+= '--PolMode Scalar '
    syscall+= '--dt '+str(dt)+' '
    syscall+= '--NChanSols '+str(nchansols)+' '
    # [KAFCA]
    syscall+= '--NIterKF '+str(niterkf)+' '
    syscall+= '--CovQ '+str(covq)+' '

    return syscall


def generate_syscall_pybdsf(fitsfile,
                        thresh_pix = cfg.PYBDSF_THRESH_PIX,
                        thresh_isl = cfg.PYBDSF_THRESH_ISL,
                        catalogtype = cfg.PYBDSF_CATALOGTYPE,
                        catalogformat = cfg.PYBDSF_CATALOGFORMAT):

    if catalogtype == 'srl':
        opfile = fitsfile+'.srl'
    elif catalogtype == 'gaul':
        opfile = fitsfile+'.gaul'

    if catalogformat == 'fits':
        opfile += '.fits'

    syscall = "python3 -c '"
    syscall += "import bdsf; "
    syscall += "img = bdsf.process_image(\""+fitsfile+"\","
    syscall += "thresh_pix="+str(thresh_pix)+","
    syscall += "thresh_isl="+str(thresh_isl)+","
    syscall += "adaptive_rms_box=True) ; "
    syscall += "img.write_catalog(outfile=\""+opfile+"\","
    syscall += "format=\""+catalogformat+"\","
    syscall += "catalog_type=\""+catalogtype+"\","
    syscall += "clobber=True,incl_empty=True)'"

    return syscall,opfile


def generate_syscall_clustercat(srl,
                        ndir = cfg.CLUSTERCAT_NDIR,
                        centralradius = cfg.CLUSTERCAT_CENTRALRADIUS,
                        ngen = cfg.CLUSTERCAT_NGEN,
                        fluxmin = cfg.CLUSTERCAT_FLUXMIN,
                        ncpu = cfg.CLUSTERCAT_NCPU):

    opfile = srl.replace('.srl.fits','.srl.fits.'+str(ndir)+'.dirs.ClusterCat.npy')
    syscall = 'ClusterCat.py --SourceCat '+srl+' '
    syscall += '--NGen '+str(ngen)+' '
    syscall += '--NCluster '+str(ndir)+' '
    syscall += '--FluxMin='+str(fluxmin)+' '
    syscall += '--CentralRadius='+str(centralradius)+' '
    syscall += '--NCPU='+str(ncpu)+' '
    syscall += '--DoPlot=0 '
    syscall += '--OutClusterCat='+opfile

    return syscall, opfile


def generate_syscall_crystalball(myms,
                        model,
                        outcol,
                        region,
                        num_workers=32,
                        mem_fraction=90):

    syscall = 'crystalball '
    syscall += '-sm '+model+' '
    syscall += '-o '+outcol+' '
    syscall += '-w '+region+' '
    syscall += '--spectra '
    syscall += '-j '+num_workers+' '
    syscall += '-mf '+mem_fraction+' '
    syscall += myms

    return syscall




# ------------------------------------------------------------------------


# ========================================================================
# GLAMDRING-SPECIFIC QUEUE SELECTION FUNCTIONS
# ========================================================================
# Functions for automatically selecting the best Slurm queue based on
# resource availability on the Glamdring cluster.
# Based on get_queue_summary.py functionality.
# ========================================================================

import re
import subprocess
from typing import Dict, List, Any, Tuple, Set


_NUM_PAIR_RE = re.compile(
    r"^\s*([0-9]+(?:\.[0-9]+)?)\s*/\s*([0-9]+(?:\.[0-9]+)?)\s*$"
)


def run_n() -> str:
    """
    Execute the 'n' command and return its output.
    
    Returns:
        str: stdout from the 'n' command
    """
    cp = subprocess.run(["n"], check=True, text=True, capture_output=True)
    return cp.stdout


def parse_free_value(cell: str) -> float:
    """
    Parse a 'free/total' cell value and return the free amount.
    
    Args:
        cell: String in format "123/456" or "12.5/24.0"
    
    Returns:
        float: The free value (first number)
    
    Raises:
        ValueError: If the cell doesn't match the expected format
    """
    m = _NUM_PAIR_RE.match(cell)
    if not m:
        raise ValueError
    return float(m.group(1))


def parse_n_table(text: str) -> Dict[str, Dict[str, List[Any]]]:
    """
    Parse the output table from 'n' command into structured data.
    
    Memory units from 'n' are assumed to be in GB.
    
    Args:
        text: Raw output from 'n' command
    
    Returns:
        Dict mapping queue names to:
            {"states": [...], "cores": [free...], "mem": [free_gb...]}
    """
    data: Dict[str, Dict[str, List[Any]]] = {}

    for raw in text.splitlines():
        line = raw.strip()
        if not line.startswith("│"):
            continue

        cells = [c.strip() for c in line.split("│")[1:-1]]
        if len(cells) < 6:
            continue

        queue, node, state, cores_cell, mem_cell, _os = cells[:6]

        if node in ("QUEUE TOTAL", "ALL TOTAL"):
            continue
        if not queue:
            continue

        try:
            free_cores = parse_free_value(cores_cell)
            free_mem = parse_free_value(mem_cell)
        except ValueError:
            continue

        bucket = data.setdefault(queue, {"states": [], "cores": [], "mem": []})
        bucket["states"].append(state)
        bucket["cores"].append(free_cores)
        bucket["mem"].append(free_mem)

    return data


def normalise(q: str) -> str:
    """Normalize queue name to uppercase."""
    return q.strip().upper()


def allowed_states(include_mixed: bool) -> Set[str]:
    """
    Return the set of allowed node states.
    
    Args:
        include_mixed: If True, include MIXED state; otherwise only IDLE and OFFtillNEEDED
    
    Returns:
        Set of allowed state strings
    """
    s = {"IDLE", "OFFtillNEEDED"}
    if include_mixed:
        s.add("MIXED")
    return s


def count_nodes(
    data: Dict[str, Dict[str, List[Any]]],
    min_cores: float,
    min_mem: float,
    states_ok: Set[str],
    excluded: Set[str],
) -> Dict[str, int]:
    """
    Count qualifying nodes per queue based on resource requirements.
    
    Args:
        data: Parsed queue data from parse_n_table
        min_cores: Minimum free cores required per node
        min_mem: Minimum free memory required per node (GB)
        states_ok: Set of allowed node states
        excluded: Set of excluded queue names (normalized)
    
    Returns:
        Dict mapping queue names to count of qualifying nodes
    """
    counts: Dict[str, int] = {}

    for queue, info in data.items():
        if normalise(queue) in excluded:
            continue

        count = 0
        for st, c, m in zip(info["states"], info["cores"], info["mem"]):
            if st in states_ok and c >= min_cores and m >= min_mem:
                count += 1

        counts[queue] = count

    return counts


def pick_best(counts: Dict[str, int]) -> Tuple[List[str], int]:
    """
    Select the queue(s) with the most qualifying nodes.
    
    Args:
        counts: Dict mapping queue names to qualifying node counts
    
    Returns:
        Tuple of (list of best queue names, best count)
    """
    best_count = -1
    best: List[str] = []
    for q, c in counts.items():
        if c > best_count:
            best_count = c
            best = [q]
        elif c == best_count:
            best.append(q)
    return best, best_count


def select_best_queue(
    min_cores: float,
    min_mem: float,
    cpu_only_states: bool = False,
    verbose: bool = False,
) -> Tuple[str, int]:
    """
    Select the best Slurm queue based on resource availability.
    
    This function runs the 'n' command, parses its output, and selects
    the queue with the most nodes meeting the specified requirements.
    
    Memory units:
    - min_mem is interpreted as GB of free memory per node
    - The 'Mem Free/TOT' column from 'n' is assumed to be in GB
    
    Args:
        min_cores: Minimum free cores required per node
        min_mem: Minimum free memory required per node (GB)
        cpu_only_states: If True, only allow IDLE and OFFtillNEEDED states
                        (exclude MIXED). Default False.
        verbose: If True, print detailed information about the selection.
                Default False.
    
    Returns:
        Tuple of (best queue name, number of qualifying nodes)
        If no queue qualifies, returns ("redwood", 0) as default
    
    Example:
        >>> queue, count = select_best_queue(min_cores=32, min_mem=100)
        >>> print(f"Selected queue: {queue} with {count} available nodes")
    """
    # GPU queues to exclude
    excluded = {normalise(x) for x in ("GPULONG", "OPTGPU", "CMBGPU")}
    
    # Determine allowed states
    states_ok = allowed_states(include_mixed=not cpu_only_states)

    try:
        text = run_n()
    except Exception as e:
        if verbose:
            print(f"Error running 'n' command: {e}")
        return "redwood", 0

    # Parse and analyze
    data = parse_n_table(text)
    counts = count_nodes(
        data,
        min_cores=min_cores,
        min_mem=min_mem,
        states_ok=states_ok,
        excluded=excluded,
    )

    best_queues, best_count = pick_best(counts)

    # Default to redwood if nothing qualifies
    if best_count <= 0:
        best_queues = ["redwood"]
        best_count = 0

    if verbose:
        print(f"Requirements: free_cores >= {min_cores}, free_mem >= {min_mem} GB")
        print(f"Allowed states: {', '.join(sorted(states_ok))}")
        print(f"Excluded queues: {', '.join(sorted(excluded))}")
        print("\nPer-queue qualifying node counts:")
        for q in sorted(counts):
            print(f"  {q}: {counts[q]}")

        if best_count == 0:
            print("\nNo queues meet requirements --- defaulting to: redwood")
        elif len(best_queues) == 1:
            print(f"\nBest queue: {best_queues[0]} ({best_count} qualifying nodes)")
        else:
            print(f"\nBest queues (tie at {best_count} nodes): {', '.join(sorted(best_queues))}")

    # Return first queue if there's a tie
    return best_queues[0], best_count


def write_glam_submission_script(file_obj, steps, field_name=None):
    """
    Write addqueue job submission commands to an opened file object.
    
    This function generates a bash script that submits jobs using addqueue,
    with automatic queue selection and dependency handling for Glamdring.
    
    Args:
        file_obj: An opened file object to write commands to
        steps: List of step dictionaries containing job configuration.
               Each step should have:
               - 'jobname': Name of the job
               - 'glam_config': Dict with 'CPUS', 'MEM', 'QUEUE' keys
               - 'dependency': Optional job name this depends on
        field_name: Optional field name to include in job names
    
    The function will:
    - Automatically determine the best queue using select_best_queue()
    - Handle multi-part scripts (part000, part001, etc.) as dependency chains
    - Extract Slurm job IDs and use them for --runafter dependencies
    - Create log files in cfg.LOGS with standard naming
    
    Example usage:
        with open('submit_jobs.sh', 'w') as f:
            f.write('#!/usr/bin/env bash\\n')
            f.write('set -euo pipefail\\n\\n')
            write_glam_submission_script(f, steps)
    """
    import glob
    
    file_obj.write('\n# Auto-generated job submission commands for Glamdring\n\n')
    
    # Track job IDs for dependencies
    job_ids = {}  # Maps jobname to variable holding its JOBID
    
    # Run queue selection once for all jobs and store results
    file_obj.write('# Queue Selection Summary:\n')
    
    # Get current queue status once
    try:
        text = run_n()
        data = parse_n_table(text)
        excluded = {normalise(x) for x in ("GPULONG", "OPTGPU", "CMBGPU")}
        states_ok = allowed_states(include_mixed=True)
        
        # Print overall queue status
        file_obj.write('#\n')
        file_obj.write('# Current cluster status:\n')
        for queue_name in sorted(data.keys()):
            if normalise(queue_name) not in excluded:
                total_nodes = len(data[queue_name]['cores'])
                avg_free_cores = sum(data[queue_name]['cores']) / max(1, total_nodes)
                avg_free_mem = sum(data[queue_name]['mem']) / max(1, total_nodes)
                file_obj.write(f'#   {queue_name}: {total_nodes} nodes, avg {avg_free_cores:.1f} free cores, avg {avg_free_mem:.1f} GB free mem\n')
        file_obj.write('#\n')
    except Exception as e:
        file_obj.write(f'# Warning: Could not query cluster status: {e}\n')
        data = {}
    
    # Pre-compute queue selections for each step
    step_info = []
    
    for step in steps:
        jobname = step['id']
        glam_config = step.get('glam_config', cfg.GLAM_STANDARD)
        dependency = step.get('dependency', None)
        
        # Resolve dependency index to jobname
        dependency_jobname = None
        if dependency is not None:
            if isinstance(dependency, list):
                # Multiple dependencies - keep as list of jobnames
                dependency_jobname = []
                for dep in dependency:
                    if isinstance(dep, int):
                        # dependency is a step index, look up the jobname
                        if 0 <= dep < len(steps):
                            dependency_jobname.append(steps[dep]['id'])
                    else:
                        # dependency is already a jobname string
                        dependency_jobname.append(dep)
            elif isinstance(dependency, int):
                # Single dependency as step index, look up the jobname
                if 0 <= dependency < len(steps):
                    dependency_jobname = steps[dependency]['id']
            else:
                # dependency is already a jobname string
                dependency_jobname = dependency
        
        # Extract resource requirements
        cpus = glam_config.get('CPUS', '32')
        mem_gb = glam_config.get('MEM', '100GB')
        
        # Parse memory string to GB number (with headroom for queue selection)
        mem_value = mem_string_to_gb(mem_gb)
        
        # Parse memory value without headroom for addqueue command
        mem_gb_clean = mem_gb.upper().replace('B', '').replace('G', '').replace('T', '').replace('M', '')
        mem_numeric = int(float(''.join(x for x in mem_gb_clean if x.isdigit() or x == '.')))
        
        # Select best queue based on requirements
        queue, node_count = select_best_queue(
            min_cores=float(cpus),
            min_mem=float(mem_value),
            cpu_only_states=False,
            verbose=False
        )
        
        # Store all info for this step
        step_info.append({
            'step': step,
            'jobname': jobname,
            'cpus': cpus,
            'mem_gb': mem_gb,
            'mem_value': mem_value,
            'mem_numeric': mem_numeric,
            'queue': queue,
            'node_count': node_count,
            'dependency': dependency_jobname
        })
        
        file_obj.write(f'#   {jobname}: Queue={queue}, Cores={cpus}, Mem={mem_gb} ({node_count} nodes available)\n')
    
    file_obj.write('\n')
    
    # Now process each step with pre-computed queue info
    for info in step_info:
        jobname = info['jobname']
        cpus = info['cpus']
        mem_gb = info['mem_gb']
        mem_numeric = info['mem_numeric']
        queue = info['queue']
        node_count = info['node_count']
        dependency = info['dependency']
        
        # Find all script files for this job (may be multi-part)
        # Use exact match to avoid matching POEXT077_1 with POEXT077_10, etc.
        script_pattern = cfg.SCRIPTS + '/glam_' + jobname + '_part*.sh'
        script_files = sorted(glob.glob(script_pattern))
        
        # If no multi-part files, check for single script
        if not script_files:
            single_script = cfg.SCRIPTS + '/glam_' + jobname + '.sh'
            if glob.glob(single_script):
                script_files = [single_script]
        
        if not script_files:
            file_obj.write(f'# WARNING: No script files found for {jobname}\n')
            continue
        
        # Log file for this job
        logfile = cfg.LOGS + '/slurm_' + jobname + '.log'
        
        file_obj.write(f'\n# ---- Job: {jobname} ----\n')
        file_obj.write(f'# Queue: {queue}, Cores: {cpus}, Mem: {mem_gb}\n')
        
        # Previous job ID for chaining multi-part scripts
        prev_jobid_var = None
        
        # If this job depends on another job, start with that dependency
        if dependency:
            if isinstance(dependency, list):
                # Multiple dependencies - join only first and last (range notation)
                dep_vars = []
                for dep_name in dependency:
                    if dep_name in job_ids:
                        dep_vars.append(job_ids[dep_name])
                if dep_vars:
                    if len(dep_vars) == 1:
                        # Single dependency
                        prev_jobid_var = dep_vars[0]
                    else:
                        # Multiple dependencies: use first:last notation (includes all in between)
                        prev_jobid_var = f'${{{dep_vars[0]}}}:${{{dep_vars[-1]}}}'
            else:
                # Single dependency
                if dependency in job_ids:
                    prev_jobid_var = job_ids[dependency]
        
        # Submit each part
        for idx, script_file in enumerate(script_files):
            script_basename = script_file.split('/')[-1]
            
            # Variable name and job name for this part's job ID
            if len(script_files) > 1:
                jobid_var = f"{jobname}_part{idx:03d}"
                job_submit_name = f"{jobname}_part{idx:03d}"
                part_logfile = cfg.LOGS + f'/slurm_{jobname}_part{idx:03d}.log'
            else:
                jobid_var = jobname
                job_submit_name = jobname
                part_logfile = cfg.LOGS + f'/slurm_{jobname}.log'
            
            # Build addqueue command (memory is per-cpu)
            mem_per_cpu = int(mem_numeric / int(cpus))
            file_obj.write(f'\n# Submit {script_basename}\n')
            file_obj.write(f'out="$(addqueue -c "{job_submit_name}" -n 1x{cpus} -m {mem_per_cpu} ')
            file_obj.write(f'-o "{part_logfile}" -q "{queue}" -s ')
            
            # Add dependency if exists
            if prev_jobid_var:
                # Check if it's a formatted multi-dependency string (contains ${ and :)
                if '${' in str(prev_jobid_var) and ':' in str(prev_jobid_var):
                    # Already formatted as ${VAR1}:${VAR2}, use directly
                    file_obj.write(f'--runafter "{prev_jobid_var}" ')
                else:
                    # Single variable, wrap with ${}
                    file_obj.write(f'--runafter "${{{prev_jobid_var}}}" ')
            
            file_obj.write(f'--sbatch "{script_file}")"\n')
            
            # Extract job ID
            file_obj.write(f'{jobid_var}="$(printf \'%s\\n\' "$out" | ')
            file_obj.write(f'awk \'/Submitted batch job/ {{print $NF; exit}}\')"\n')
            
            # Error check
            file_obj.write(f'if [[ -z "${{{jobid_var}:-}}" ]]; then\n')
            file_obj.write(f'  echo "ERROR: Could not extract JOBID for {script_basename}" >&2\n')
            file_obj.write(f'  echo "---- addqueue output ----" >&2\n')
            file_obj.write(f'  echo "$out" >&2\n')
            file_obj.write(f'  exit 1\n')
            file_obj.write(f'fi\n')
            file_obj.write(f'echo "Submitted {script_basename} with JOBID=${{{jobid_var}}}"\n')
            
            # This part becomes the dependency for the next part
            prev_jobid_var = jobid_var
        
        # Store the last job ID for this step (for other steps that depend on it)
        job_ids[jobname] = prev_jobid_var
    
    file_obj.write('\n# All jobs submitted successfully\n')


# ========================================================================
