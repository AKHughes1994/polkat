#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk

import glob
import datetime
import time
import os
import os.path as o
import sys

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg

def col(txt=''):
    colstr = ' '+txt.ljust(20)+'| '
    return colstr

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
                          chanout = cfg.WSC_IMAGE_CHANNELSOUT,
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
                          squarepolarizations = cfg.WSC_SQUAREPOLARIZATIONS,
                          pol = cfg.WSC_POL,
                          splitpol = cfg.WSC_SPITPOL,
                          padding = cfg.WSC_PADDING,
                          nomodel = cfg.WSC_NOMODEL,
                          mask = cfg.WSC_MASK,
                          mfweight = cfg.WSC_MFWEIGHT,
                          threshold = cfg.WSC_THRESHOLD,
                          autothreshold = cfg.WSC_AUTOTHRESHOLD,
                          automask = cfg.WSC_AUTOMASK,
                          localrms = cfg.WSC_LOCALRMS,
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
                          parallelgridding = cfg.WSC_PARALLELGRIDDING):

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
    if absmem < 0:
        syscall += '-mem '+str(mem)+' '
    else:
        syscall += '-abs-mem '+str(absmem)+' '
    if continueclean:
        syscall += '-continue '
    if parallelreordering != 0:
        syscall += '-parallel-reordering '+str(parallelreordering)+' '

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
    if interval0 and interval1:
        syscall += '-intervals '+str(interval0)+' '+str(interval1)+' '
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
        syscall += f'-minuv-l 0.0 -taper-inner-tukey {tukeytaper} '


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
    if chandeconvolution:
        syscall += '--deconvolution-channels '+str(chandeconvolution)+' '
    if joinchannels:
        syscall += '-join-channels '

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
    if threshold:
        syscall += '-threshold '+str(threshold)+' '


    # Andrew additions for more complex imaging strategies 
    wsclean_syscall_base = syscall # Copy syscall options then append additional options
    syscall_arr = []

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
            syscall_arr.append(wsclean_syscall_base + f'-channels-out {maxchan} -channel-range {k * intchans} {(k + 1) * intchans - 1} -name {imgname}_part{k:04d} ') 

    else:
        syscall_arr.append(wsclean_syscall_base + f'-channels-out {chanout} -name {imgname} ')

    # Add option to split the deconvolution into IV and QU steps
    # Sources with large rotation measures may not want polynomial fitting to the QU channels
    if splitpol and pol != 'I':
        spectralpol = ''

        pol_QU = pol.replace('I', '').replace('V','')
        pol_IV  = pol.replace('Q', '').replace('U','')

        if fitspectralpol != 0:
            spectralpol = '-fit-spectral-pol '+str(fitspectralpol) + ' '

        joinpol_QU = ''
        squarepol = ''

        if joinpolarizations and len(pol.replace('I', '').replace('V','')) >= 2:
            joinpol_QU     = '-join-polarizations '
            if squarepolarizations:
                squarepol = '-squared-channel-joining '    

        joinpol_IV = ''
        if joinpolarizations and len(pol.replace('Q', '').replace('U','')) >= 2: # this makes sure to not join polarization for 'I' only imaging
            joinpol_IV     = '-join-polarizations '

        k = len(syscall_arr)
        syscall_arr += syscall_arr

        for _k in range(k):
            syscall_arr[_k] += f'-pol {pol_IV} {spectralpol} {joinpol_IV} '
            syscall_arr[_k + k] += f'-pol {pol_QU} {joinpol_QU} {squarepol} '

    else:
        joinpol = ''
        squarepol = ''
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


    return syscall_arr

