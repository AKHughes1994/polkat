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
# syscall += cfg.CWD + '/*1024ch_*.ms '
syscall += cfg.CWD + '/*.parmdb '
syscall += cfg.CWD + '/*.skel '
syscall += cfg.CWD + '/*.html '
syscall += cfg.IMAGES + '/*datablind-0* '
syscall += cfg.IMAGES + '/*datablind-MFS* '
syscall += cfg.IMAGES + '/*datamask-0* '
syscall += cfg.IMAGES + '/*datamask-MFS* '
syscall += cfg.IMAGES + '/*uniform-0* '
syscall += cfg.IMAGES + '/*uniform-MFS* '
syscall += cfg.IMAGES + '/*pcalmask-0* '
syscall += cfg.IMAGES + '/*pcalmask-MFS* '
syscall += cfg.IMAGES + '/*snapmask* '
syscall += cfg.IMAGES + '/*diagnostic-MFS* '
syscall += cfg.IMAGES + '/*diagnostic-0* '
syscall += cfg.IMAGES + '/*snapblind* '
syscall += cfg.IMAGES + '/*psf* '
syscall += cfg.IMAGES + '/*dirty* '
syscall += cfg.IMAGES + '/*.pb.* '
syscall += cfg.IMAGES + '/*.wt.* '
syscall += cfg.IMAGES + '/*kernel* '
subprocess.run([syscall], shell=True)

syscall   = 'rm -rf '
syscall += cfg.INTERVALS + '/*psf* '
syscall += cfg.INTERVALS + '/*modelsub*.fits '
syscall += cfg.INTERVALS + '/*-psf*.fits '
subprocess.run([syscall], shell=True)
