#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk


import glob
import os
import os.path as o
import subprocess
import sys
import time
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg


# If the first ragavi-gains call takes longer than this, assume the interactive
# bokeh HTML output is the bottleneck and fall back to PNG-only plotting with a
# non-interactive matplotlib backend for every remaining table.
HTML_TIMEOUT_SECONDS = 60.0


def main():


    GAINPLOTS = cfg.GAINPLOTS
    GAINTABLES = cfg.GAINTABLES
    gen.setup_dir(GAINPLOTS)


    include = sys.argv[1]
    if len(sys.argv) > 2:
        exclude = sys.argv[2]
    else:
        exclude = ''

    caltabs = sorted([item for item in glob.glob(GAINTABLES+'/'+include) if not os.path.basename(item).endswith('flagversions')])
    if exclude != '':
        exclude = glob.glob(GAINTABLES+'/'+exclude)

    caltabs = [caltab for caltab in caltabs if caltab not in exclude]

    png_only = False    # set once the first call trips HTML_TIMEOUT_SECONDS
    timed_first = False # the first call that actually ran (skips don't count)

    for idx, caltab in enumerate(caltabs):
        gaintype = caltab.split('.')[-1]
        gaintype = ''.join([i for i in gaintype if not i.isdigit()])
        if gaintype == 'Ga' or gaintype == 'Gp':
            gaintype = 'G'
        htmlname = GAINPLOTS+'/'+caltab.split('/')[-1]+'.html'
        plotname = GAINPLOTS+'/'+caltab.split('/')[-1]+'.png'

        # In PNG-only mode no HTML is written, so the PNG becomes the marker
        # for work that has already been done.
        existing = plotname if png_only else htmlname
        if os.path.isfile(existing):
            print(existing+' exists, skipping')
            continue

        syscall = 'ragavi-gains -g '+gaintype+' -t '+caltab+' --plotname='+plotname
        if not png_only:
            syscall += ' --htmlname='+htmlname

        t0 = time.time()
        subprocess.run([syscall],shell=True)
        elapsed = time.time() - t0

        if not timed_first:
            timed_first = True
            print('First gain plot ('+os.path.basename(caltab)+') took '+format(elapsed,'.1f')+' s')
            if elapsed > HTML_TIMEOUT_SECONDS:
                # Inherited by every subsequent ragavi-gains subprocess.
                os.environ['MPLBACKEND'] = 'Agg'
                png_only = True
                n_remaining = len(caltabs) - idx - 1
                print('')
                print('  That is over the '+format(HTML_TIMEOUT_SECONDS,'.0f')+' s threshold -- the interactive HTML')
                print('  output is the likely bottleneck. Falling back to PNG-only plotting:')
                print('    MPLBACKEND=Agg exported for all remaining ragavi-gains calls')
                print('    --htmlname dropped -- '+str(n_remaining)+' remaining table(s) will produce .png only')
                print('    no .html will be written for them; delete the .png to re-plot')
                print('')

if __name__ == "__main__":


    main()
