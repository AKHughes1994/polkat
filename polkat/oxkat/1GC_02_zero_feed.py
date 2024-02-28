#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


import json
import os
import sys
import numpy as np
    
# Get ms file information
with open('project_info.json') as f:
    project_info = json.load(f)

myms = project_info['working_ms']

# Open the ms file table and set receptor angle to 0.0 for all entries
tb.open(f"{myms}::FEED", nomodify=False)
fa = tb.getcol("RECEPTOR_ANGLE")
fa[...] = 0.0
tb.putcol("RECEPTOR_ANGLE", fa)
tb.flush()
tb.close()
   
# This code will re-open and output the receptor angles, to double check it worked
tb.open(f"{myms}::FEED", nomodify=True)
fa = tb.getcol("RECEPTOR_ANGLE")
print('MAX/MIN angles are (should both be 0.0): ', np.amax(fa), np.amin(fa))
tb.close()

# Save flags
flagmanager(vis=myms, versionname='after_swap_and_zero', mode='save')
