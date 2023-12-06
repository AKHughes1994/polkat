import sys, os, json
    
exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

ms = str(sys.argv[-1])
tb.open("{}::FEED".format(ms), nomodify=False)
fa = tb.getcol("RECEPTOR_ANGLE")
fa[...] = 0.0
tb.putcol("RECEPTOR_ANGLE", fa)
tb.flush()
tb.close()
   
tb.open("{}::FEED".format(ms), nomodify=False)
fa = tb.getcol("RECEPTOR_ANGLE")
print(fa)
tb.close()

if SAVE_FLAGS:
    flagmanager(vis=ms, versionname='before_cal', mode='save')

