#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk
#
# Per-target summary of MFS flux/polarization and RM synthesis results.
# Runs once for every science target present in this run (see
# setups/RMSYNTH.py), i.e. every entry of target_names that is also in
# working_names. Writes one plain-text file per target to the CWD.
#
# Matches the *actual* output of oxkat/RMSYNTH_01_extract_fluxes.py for this
# branch: single component (no per-component nesting), scalar MFS values
# (not lists -- each json is already one epoch), and a timestamped filename
# ({timestamp_prefix}{src_name}_pcalmask_polarization.json) rather than
# dev's img_<ms>_pcalmask_polarization.json. A target may have more than one
# matching json when RMSYNTH_01 ran in time-interval mode (one file per
# interval, src_name suffixed with _t0000 etc.) -- all matches are summarized.
#
# Input files and the keys/columns read from each (RESULTS = cfg.RESULTS):
#
#   RESULTS/*<target>*_pcalmask_polarization.json
#       ['MFS']:
#           I_flux_mJy, I_err_mJy   -- MFS Stokes I +/- error
#           alpha, alpha_err        -- spectral index (absent if not fit)
#           chi2, ndof, chi2_red    -- spectral index fit quality (absent if not fit)
#           time_ctr_mjd            -- MJD at the centre of the exposure
#           time_ctr_isot           -- ISOT date/time at the centre of the exposure
#           time_dt                 -- exposure time (s)
#           freq_GHz                -- MFS central frequency
#           I_RA_deg, I_DEC_deg     -- fitted position (degrees)
#
#   RESULTS/*<target>*_pcalmask_rmsynth_RMclean.json
#       phiPeakPIfit_rm2       -- RM (rad/m^2), post RM-CLEAN
#       dPhiObserved_rm2       -- "effective" RM error (from the corrected-MAD FDF noise,
#                                  rather than the purely theoretical fit error)
#       polAngleFit_deg        -- polarization angle (deg)
#       dPolAngleFitObserved_deg -- "effective" angle error (as above)
#       snrPIfit                -- SNR of the RM synthesis peak
#       ampPeakPIfitEff         -- debiased ("effective") peak polarized amplitude, used
#                                  together with Ifreq0 (below) for the representative
#                                  fractional polarization
#
#   RESULTS/*<target>*_pcalmask_rmsynth_RMsynth.json
#   (sibling of the RMclean.json above, before RM-CLEAN; carries keys that
#    RMclean.json does not)
#       freq0_Hz                               -- central frequency used for the RM fit
#       Ifreq0                                  -- Stokes I at freq0
#
# Position (hmsdms) is derived from I_RA_deg/I_DEC_deg above via
# astropy.coordinates.SkyCoord, not read from a file.
#
# NOTE: dev's polarization-angle calibrator diagnostic cross-check (written by
# RMSYNTH_01B_systematics.py, glob 'img_*_diagnostic_rmsynth_RMclean.json' in
# dev) is intentionally omitted here -- RMSYNTH_01B_systematics.py on this
# branch names its output per-calibrator (see cal_json_path() there), not with
# a 'diagnostic' substring, so dev's glob pattern would never match. Wiring
# this up needs that naming traced separately; left out rather than guessed.

import glob
import json
import os.path as o
import sys

from astropy.coordinates import SkyCoord
from astropy import units as u

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg

# Must match IDENTIFIER in RMSYNTH_01_extract_fluxes.py
IDENTIFIER = 'pcalmask'


def fmt(value, err=None, unit=''):
    if value is None:
        return 'N/A'
    if err is None:
        return f'{value}{unit}'
    return f'{value} +/- {err}{unit}'


def hz_to_ghz(hz):
    return None if hz is None else hz / 1.0e9


def representative_frac_pol(rmclean, rmsynth):
    """Pol fraction (%) = debiased peak amplitude (post RM-CLEAN) / Stokes I at freq0."""
    amp    = rmclean.get('ampPeakPIfitEff')
    ifreq0 = rmsynth.get('Ifreq0')
    if amp is None or not ifreq0:
        return None
    return amp / ifreq0 * 100.0


def write_rmsynth_block(line, rmclean, rmsynth):
    line('RM (rad/m^2)', fmt(rmclean.get('phiPeakPIfit_rm2'), rmclean.get('dPhiObserved_rm2')))
    line('Polarization angle (deg)', fmt(rmclean.get('polAngleFit_deg'), rmclean.get('dPolAngleFitObserved_deg')))
    line('SNR of rmsynth', fmt(rmclean.get('snrPIfit')))
    line('Central frequency for rmsynth (GHz)', fmt(hz_to_ghz(rmsynth.get('freq0_Hz'))))
    line('Stokes I at freq0', fmt(rmsynth.get('Ifreq0')))
    line('Pol. fraction, representative (AMP_eff/I0, %)', fmt(representative_frac_pol(rmclean, rmsynth)))


def summarize_epoch(out, line, pol_json_path):

    with open(pol_json_path) as f:
        mfs = json.load(f).get('MFS', {})

    rmclean_path = pol_json_path.replace(f'_{IDENTIFIER}_polarization.json',
                                          f'_{IDENTIFIER}_rmsynth_RMclean.json')
    rmsynth_path = rmclean_path.replace('_RMclean.json', '_RMsynth.json')
    rmclean = {}
    if o.isfile(rmclean_path):
        with open(rmclean_path) as f:
            rmclean = json.load(f)
    rmsynth = {}
    if o.isfile(rmsynth_path):
        with open(rmsynth_path) as f:
            rmsynth = json.load(f)

    ra_deg = mfs.get('I_RA_deg')
    if ra_deg is not None and ra_deg < 0:
        ra_deg += 360.0
    dec_deg = mfs.get('I_DEC_deg')

    hmsdms = None
    if ra_deg is not None and dec_deg is not None:
        hmsdms = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg).to_string('hmsdms')

    epoch_label = o.basename(pol_json_path).replace(f'_{IDENTIFIER}_polarization.json', '')
    out.write(f'\n--- {epoch_label} ---\n')
    line('MFS Stokes I (mJy)', fmt(mfs.get('I_flux_mJy'), mfs.get('I_err_mJy')))
    if mfs.get('alpha') is not None:
        line('Spectral index', fmt(mfs.get('alpha'), mfs.get('alpha_err')))
        line('Spectral fit chi2 / ndof / chi2_red',
             f"{fmt(mfs.get('chi2'))} / {fmt(mfs.get('ndof'))} / {fmt(mfs.get('chi2_red'))}")
    else:
        line('Spectral index', 'N/A')
    line('Date/time, centre (ISOT, UTC)', fmt(mfs.get('time_ctr_isot')))
    line('MJD, centre', fmt(mfs.get('time_ctr_mjd')))
    line('Exposure time (s)', fmt(mfs.get('time_dt')))
    line('Central frequency (GHz)', fmt(mfs.get('freq_GHz')))
    if rmclean or rmsynth:
        write_rmsynth_block(line, rmclean, rmsynth)
    line('Position (deg)', f'{fmt(ra_deg)}, {fmt(dec_deg)}')
    line('Position (hmsdms)', hmsdms if hmsdms is not None else 'N/A')


def summarize_target(target, myms):

    pol_json_paths = sorted(glob.glob(o.join(cfg.RESULTS, f'*{target}*_{IDENTIFIER}_polarization.json')))
    if not pol_json_paths:
        print(gen.col('Target Summary')+f'Skipping {target}: no *{target}*_{IDENTIFIER}_polarization.json found in {cfg.RESULTS}')
        return

    out_path = f'{myms}_summary.txt'
    with open(out_path, 'w') as out:
        out.write(f'Summary for {target} ({myms})\n')
        out.write('='*60+'\n')

        def line(label, value):
            out.write(f'{label.ljust(40)}: {value}\n')

        for pol_json_path in pol_json_paths:
            summarize_epoch(out, line, pol_json_path)

    print(gen.col('Target Summary')+f'Wrote {out_path}')


def main():

    with open('project_info.json') as f:
        project_info = json.load(f)

    target_names  = project_info.get('target_names', [])
    target_ms     = project_info.get('target_ms', [])
    working_names = project_info.get('working_names', [])
    targets = [t for t in target_names if t in working_names]

    if not targets:
        print(gen.col('Target Summary')+f'Skipping: no targets found in working_names')
        return

    for target in targets:
        myms = target_ms[target_names.index(target)]
        summarize_target(target, myms)


if __name__ == "__main__":

    main()
