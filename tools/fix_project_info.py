# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk
import numpy as np
import glob, json

exec(open('oxkat/config.py').read())
exec(open('oxkat/casa_read_project_info.py').read())

myfields = PRE_FIELDS
myscans = PRE_SCANS
myoutputchans = int(PRE_NCHANS)
mytimebins = PRE_TIMEBIN

# Load in the master ms
opms = master_ms.replace('.ms','_'+str(myoutputchans)+'ch.ms')

flagmanager(opms, versionname='after_swap_and_zero', mode='restore')

# Remove old flag versions
for fname in ['1GC_flags', 'after_pcal_SwiftJ1727', 'applycal_1', 'autoflag_cals_data', 'basic', 'bpcal_residual_flags']:
    try:
        flagmanager(opms, versionname=fname, mode='delete')
    except:
        pass

# Get names and field IDs for sources
tb.open(opms+'/FIELD')
names = tb.getcol('NAME')
ids   = tb.getcol('SOURCE_ID')
tb.done()

# Append the working names and IDs to project info
with open('project_info.json','r') as j:
    project_info = json.load(j)

project_info['working_names'] = names.tolist()
project_info['working_ids'] = ids.tolist()


# -----------------------------------------------------------------------
# Re-select reference antenna based on current flag state
# -----------------------------------------------------------------------

if CAL_1GC_REF_ANT != 'auto':
    # User has specified a fixed refant — respect it
    ref_ant = CAL_1GC_REF_ANT
    print(f'fix_project_info: Using user-specified reference antenna: {ref_ant}')
else:
    # Resolve primary field ID in the working MS via working_names/working_ids
    primary_name   = project_info['primary_name']
    working_names  = project_info['working_names']
    working_ids    = project_info['working_ids']
    if primary_name in working_names:
        primary_id = int(working_ids[working_names.index(primary_name)])
        print(f'fix_project_info: Primary "{primary_name}" -> working field ID {primary_id}')
    else:
        print(f'fix_project_info: WARNING primary "{primary_name}" not found in working_names; falling back to project_info primary_id')
        primary_id = int(project_info['primary_id'])
    ref_pool    = CAL_1GC_REF_POOL

    # Get all antenna names from the MS
    tb.open(opms + '/ANTENNA')
    ant_names = [a.lower() for a in tb.getcol('NAME').tolist()]
    tb.close()

    # Query all unique antennas present for the primary field
    _field_query = f'FIELD_ID=={primary_id}'
    print(f'fix_project_info: Querying antennas for primary field (FIELD_ID={primary_id})')
    print(f'fix_project_info: Query: {_field_query}')
    tb.open(opms)
    _field_tab = tb.query(query=_field_query)
    _ant1 = set(_field_tab.getcol('ANTENNA1').tolist())
    _ant2 = set(_field_tab.getcol('ANTENNA2').tolist())
    _field_tab.close()
    _field_ants = _ant1 | _ant2
    print(f'fix_project_info: Antenna indices present for primary field: {sorted(_field_ants)}')
    print(f'fix_project_info: Antenna names present for primary field: {[ant_names[_i] for _i in sorted(_field_ants) if _i < len(ant_names)]}')
    tb.close()

    pc_list      = []
    idx_list     = []
    ant_name_list = []

    tb.open(opms)
    for ant in ref_pool:
        if ant not in ant_names:
            continue
        idx = ant_names.index(ant)
        # Query baselines involving this antenna for the primary field
        _query = f'(ANTENNA1=={idx} || ANTENNA2=={idx}) && FIELD_ID=={primary_id}'
        print(f'fix_project_info: Query for {ant} (idx={idx}): {_query}')
        subtab = tb.query(query=_query)
        flags = subtab.getcol('FLAG')
        subtab.close()
        print(f'fix_project_info:   Rows returned: {flags.shape}, dtype: {flags.dtype}')
        if flags.size == 0:
            print(f'fix_project_info: Antenna {idx}:{ant} has no data for primary field — skipping')
            continue
        vals, counts = np.unique(flags, return_counts=True)
        print(f'fix_project_info:   FLAG unique values: {dict(zip(vals.tolist(), counts.tolist()))}')
        if len(vals) == 1:
            flag_pc = 100.0 if bool(vals[0]) else 0.0
        else:
            flag_pc = 100.0 * float(counts[1]) / float(np.sum(counts))
        print(f'fix_project_info: Antenna {idx}:{ant} is {flag_pc:.4f}% flagged')
        if flag_pc < 80.0:
            pc_list.append(flag_pc)
            idx_list.append(str(idx))
            ant_name_list.append(ant)
    tb.close()

    pc_list       = np.array(pc_list)
    idx_list      = np.array(idx_list)
    ant_name_list = np.array(ant_name_list)

    # Sort antennas by flag percentage (ascending)
    ranked_list = [x for _, x in sorted(zip(pc_list, idx_list))]

    # Prioritise m060 if its flag fraction is within 1-sigma of the median
    if 'm060' in ant_name_list:
        m060_pos     = np.where(ant_name_list == 'm060')[0][0]
        m060_flag_pc = pc_list[m060_pos]
        m060_idx     = idx_list[m060_pos]
        median_flag  = np.median(pc_list)
        std_flag     = np.std(pc_list)
        if abs(m060_flag_pc - median_flag) <= std_flag:
            if m060_idx in ranked_list:
                ranked_list.remove(m060_idx)
            ranked_list.insert(0, m060_idx)
            print(f'fix_project_info: m060 ({m060_flag_pc:.4f}%) within 1-sigma of median '
                  f'({median_flag:.4f}%), placing at front of list')

    ref_ant = ','.join(ranked_list)
    print(f'fix_project_info: New ranked reference antenna ordering: {ref_ant}')

project_info['ref_ant'] = ref_ant

with open('project_info.json','w') as j:
    json.dump(project_info, j, indent=4, sort_keys=True)
