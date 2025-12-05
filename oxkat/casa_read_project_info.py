#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


import json
import sys
exec(open('oxkat/config.py').read())


def str_iterator(inlist):
	xx = []
	for yy in inlist:
		xx.append(str(yy))
	return xx

with open('project_info.json') as f:
	project_info = json.load(f)

# Load basic MS and observation info
master_ms = str(project_info['master_ms'])
myms = str(project_info['working_ms'])
band = str(project_info['band'])
nchan = int(project_info['nchan'])
ref_ant = str(project_info['ref_ant'])

# Primary calibrator info
bpcal = str(project_info['primary_id'])
bpcal_name = str(project_info['primary_name'])
primary_tag = str(project_info['primary_tag'])

# Secondary calibrator info
pcal_names = str_iterator(project_info['secondary_names'])
pcals = str_iterator(project_info['secondary_ids'])
pcal_dirs = project_info['secondary_dirs']
pcal_ms = project_info['secondary_ms']  # List of lists: [[sec1_scan1.ms, sec1_scan2.ms], [sec2_scan1.ms, ...]]

# Target info
target_names = str_iterator(project_info['target_names'])
targets = str_iterator(project_info['target_ids'])
target_dirs = project_info['target_dirs']
target_cal_map = str_iterator(project_info['target_cal_map'])
target_ms = str_iterator(project_info['target_ms'])

# Polarization angle calibrator info
pacal = str(project_info['polang_id']) 
pacal_name = str(project_info['polang_name'])

# Filter fields if PRE_FIELDS is specified in config
if PRE_FIELDS != '':

	pre_field_list = PRE_FIELDS.split(',')

	user_targets = []
	user_pcals = []
	user_pcal_ms = []
	user_cal_map = []

	# Determine if pre_field_list uses field names or IDs
	names = False
	for src in [bpcal_name]+target_names+pcal_names:
		if src in pre_field_list:
			names = True

	# Filter by field names
	if names:
		if bpcal_name not in pre_field_list:
			print('Pre-field selection does not include a primary calibrator')
			sys.exit()
		for src in target_names:
			if src in pre_field_list:
				user_targets.append(src)
		for src in pcal_names:
			if src in pre_field_list:
				idx = pcal_names.index(src)
				user_pcals.append(src)
				user_pcal_ms.append(pcal_ms[idx])  # Keep scan-numbered MS list for this secondary

	# Filter by field IDs
	if not names:
		if bpcal not in pre_field_list:
			print('Pre-field selection does not include a primary calibrator')
			sys.exit()
		for src in targets:
			idx = targets.index(src)
			if src in pre_field_list:
				user_targets.append(target_names[idx])
		for src in pcals:
			idx = pcals.index(src)
			if src in pre_field_list:
				user_pcals.append(pcal_names[idx])
				user_pcal_ms.append(pcal_ms[idx])  # Keep scan-numbered MS list for this secondary

	# Build target-calibrator mapping for selected targets
	for src in user_targets:
		idx = target_names.index(src)
		user_cal_map.append(target_cal_map[idx])

	# Update lists to filtered versions
	pcal_ms = user_pcal_ms








