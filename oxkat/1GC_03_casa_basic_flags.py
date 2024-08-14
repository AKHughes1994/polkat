# ian.heywood@physics.ox.ac.uk


import numpy

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

args = sys.argv
for item in sys.argv:
    parts = item.split('=')
    if parts[0] == 'myms':
        myms = parts[1]


clearstat()
clearstat()


# ------------------------------------------------------------------------
# Frequency ranges to flag over all baselines

if CAL_1GC_BAD_FREQS != []:

    myspw = ','.join(CAL_1GC_BAD_FREQS)

    # myspw = ''
    # for badfreq in CAL_1GC_BAD_FREQS:
    #     myspw += '*:'+badfreq+','
    # myspw = myspw.rstrip(',')

    flagdata(vis = myms, 
        mode = 'manual', 
        spw = myspw, flagbackup=False)


# ------------------------------------------------------------------------
# Frequency ranges to flag over a subset of baselines

if CAL_1GC_BL_FREQS != []:

    myspw = ','.join(CAL_1GC_BL_FREQS)

    # myspw = ''
    # for badfreq in CAL_1GC_BL_FREQS:
    #     myspw += '*:'+badfreq+','
    # myspw = myspw.rstrip(',')

    if CAL_1GC_AGGRESSIVE_FLAG:
        flagdata(vis = myms,
            mode = 'manual',
            spw = myspw,
            #uvrange = CAL_1GC_BL_FLAG_UVRANGE,
            flagbackup=False)
    else:
        flagdata(vis = myms,
            mode = 'manual',
            spw = myspw,
            uvrange = CAL_1GC_BL_FLAG_UVRANGE,
            flagbackup=False)

# ------------------------------------------------------------------------
# Clipping, quacking, zeros, autos
# Note that clip will always flag NaN/Inf values even with a range 

#flagdata(vis = myms,
#    mode = 'quack',
#    quackinterval = 8.0,
#    quackmode = 'beg')

flagdata(vis = myms,
    mode = 'manual',
    autocorr = True, flagbackup=False)

flagdata(vis = myms,
    mode = 'clip',
    clipzeros = True, flagbackup=False)

flagdata(vis = myms,
    mode = 'clip',
    clipminmax = [0.0,100.0], flagbackup=False)


# ------------------------------------------------------------------------
# Manual Flagging
#flagdata(vis=myms, mode='manual', antenna='39,41', flagbackup=False)

# ------------------------------------------------------------------------
# Save the flags
flagmanager(vis = myms,
        mode = 'save',
        versionname = 'basic')

clearstat()
clearstat()
