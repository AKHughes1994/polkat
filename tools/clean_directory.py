import os.path as o
import sys
import subprocess
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))
from oxkat import config as cfg


# Clean up the crap
syscall   = 'rm -rf '
syscall += cfg.CWD + '/*.log '
syscall += cfg.CWD + '/*.txt '
syscall += cfg.CWD + '/*.last '
syscall += cfg.CWD + '/*scan*.ms '
syscall += cfg.CWD + '/*1024ch_*.ms '
syscall += cfg.CWD + '/*.parmdb '
syscall += cfg.CWD + '/*.skel '
syscall += cfg.CWD + '/*.html '
syscall += cfg.IMAGES + '/*datablind*-00* '
syscall += cfg.IMAGES + '/*residual* '
syscall += cfg.IMAGES + '/*datamask*-00* '
syscall += cfg.IMAGES + '/*model* '
syscall += cfg.IMAGES + '/*uniform*-00* '
syscall += cfg.IMAGES + '/*snapmask* '
syscall += cfg.IMAGES + '/*snapblind* '
syscall += cfg.IMAGES + '/*psf* '
syscall += cfg.IMAGES + '/*dirty* '
syscall += cfg.IMAGES + '/*dirty* '
syscall += cfg.IMAGES + '/*.pb.* '
syscall += cfg.IMAGES + '/*.wt.* '
syscall += cfg.IMAGES + '/*J1939*-00* '
syscall += cfg.INTERVALS + '/*psf* '
syscall += cfg.IMAGES + '/*J1331*-00*-V-* '
syscall += cfg.INTERVALS + '/*image.fits '

subprocess.run([syscall], shell=True)
