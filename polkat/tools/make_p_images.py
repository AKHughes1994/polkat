import glob, os, sys, subprocess

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

subprocess.run(['rm -rf {}/*-P-*'.format(IMAGES)], shell=True)

for qname in glob.glob('{}/*-Q-image.fits'.format(IMAGES)):
    pname = qname.replace('-Q-','-P-')
    uname = qname.replace('-Q-','-U-')
    
    immath(imagename = [qname,uname], mode='poli', outfile = pname)
