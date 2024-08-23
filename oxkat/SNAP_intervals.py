import glob
import numpy
import os
import pickle
import sys
import os.path as o
import subprocess

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg


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

def main():
    
    intervals_averaging = cfg.SNAP_AVERAGE
    intervals_channels = cfg.SNAP_CHANS
    pol = cfg.SNAP_POL
    deconv = cfg.SNAP_DECONV
    img_size = cfg.SNAP_SIZE
    snap_mask = cfg.SNAP_MASK_PATH
    
    if len(sys.argv) == 1:
        print('Please specify a field for interval imaging')
        sys.exit()
    else:
        targetname = sys.argv[1]

    # Load in the pickle file containing scan information
    scan_pickle = glob.glob(f'scantimes*.ms.p')[0]
    
    scan_times = get_scan_times(scan_pickle)
    for ss in scan_times:
        scans = ss[1]
        intervals = ss[2] #this I think is a list of the number of intervals for each scan
        
        # Only execut on target of interest
        if ss[0] == targetname:
            
            # Iterate through scans
            for i in range(0,len(scans)):
                
                number_intervals = int(intervals[i]/intervals_averaging)
                real_averaging_factor = intervals[i]/number_intervals 
                print('The real averaging factor for '+ '*'+targetname+'*scan'+str(scans[i])+'.ms' + ' is '+str(real_averaging_factor))
                
                print('*'+targetname+'*scan'+str(scans[i])+'.ms')
                myms = glob.glob('*'+targetname+'*scan'+str(scans[i])+'.ms')[0]
                img_prefix = cfg.INTERVALS+'/img_'+myms+'_modelsub'+str(intervals_channels)+'chan'
                
                #TO ADD OR NOT TO ADD JOIN CHANNELS?
                # Build command line syscall
                syscall = 'wsclean -intervals-out '+str(number_intervals)+' -interval 0 '+str(intervals[i])+' '
                syscall += f'-log-time -field 0 -no-dirty -make-psf -size {img_size} {img_size} -scale {cfg.WSC_CELLSIZE} '
                syscall += '-baseline-averaging 10 -no-update-model-required '
                
                if deconv:
                    syscall += '-gridder wgridder -niter 100 -name '+img_prefix+' '
                else:
                    syscall += '-gridder wgridder -niter 0 -name '+img_prefix+' '
                if pol:
                    syscall += '-pol IQUV -join-polarizations '
                if snap_mask != '':
                    syscall += '-fits-mask '+snap_mask+' '
                syscall += f'-weight {cfg.WSC_WEIGHT} -data-column DATA -padding 1.2 -no-mf-weighting -channels-out ' + str(intervals_channels) + ' ' + myms

                # Run syscall
                subprocess.run([syscall], shell=True)

if __name__ == "__main__":
    main()
