import sys, os, json

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

ms    = str(sys.argv[-1])
flags = str(sys.argv[-1])

# Flag the data
flagdata(vis=ms, mode='list', inpfile=flags, flagbackup=False)

if SAVE_FLAGS:
    flagmanager(vis=ms, versionname='after_manual', mode='save')
