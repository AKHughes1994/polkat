#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


import json
import os.path as o
import subprocess
import sys
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))


from oxkat import generate_jobs as gen
from oxkat import config as cfg


def main():


    VISPLOTS = cfg.VISPLOTS
    gen.setup_dir(VISPLOTS)


    with open('project_info.json') as f:
        project_info = json.load(f)

    with open('prefields_info.json') as f:
        prefield_info = json.load(f)

    field_ids = prefield_info['field_ids']
    myms = project_info['working_ms']
    
    for field_id in field_ids:
        syscall = f'python3 {cfg.TOOLS}/correct_parang.py --field {field_id} -sc "DATA" -rc "DATA" -npa -ad -cs 20000 {myms}'
        subprocess.run([syscall],shell=True)


if __name__ == "__main__":
    main()
