import os
import sys

if len(sys.argv) == 1:
    print('Please specify a Measurement Set to split')
    sys.exit()
else:
    target_ms = str(sys.argv[1])

# Open MS and get scan numbers
tb.open(target_ms)    
scans = sorted(set(tb.getcol("SCAN_NUMBER")))
tb.close()

# Split out each scan into its own MS file
for scan in scans:
    scan_ms = target_ms.replace('.ms', f'_scan{scan}.ms')
    mstransform(vis = target_ms,    
                                outputvis = scan_ms,
                                scan=scan,
                                datacolumn='all')
