# ian.heywood@physics.ox.ac.uk
import sys

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

suffix = sys.argv[-1]

if PRE_FIELDS != '':
    target_names = user_targets

for target in target_names:

    opms = ''

    for mm in target_ms:
        if target in mm:
            opms = mm.replace('.ms', suffix + '.ms')

    if opms != '':

        mstransform(vis=myms,
            outputvis=opms,
            field=target,
            usewtspectrum=True,
            realmodelcol=True,
            datacolumn='corrected')

        if SAVE_FLAGS:
            flagmanager(vis=opms,
                mode='save',
                versionname='post-1GC')

    else:

        print('Target/MS mismatch in project info for '+target+', please check.')
