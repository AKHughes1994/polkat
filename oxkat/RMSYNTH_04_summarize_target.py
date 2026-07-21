#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk
#
# Single-target summary of MFS flux/polarization and RM synthesis results.
# Only runs (see dev/setups/RMSYNTH.py) when exactly one science target is
# present in this run, i.e. len([t for t in target_names if t in
# working_names]) == 1. Writes one plain-text file to the CWD.
#
# Input files and the keys/columns read from each (RESULTS = cfg.RESULTS):
#
#   RESULTS/img_<target_ms>_pcalmask_polarization.json
#       ['MFS']['component<N>']:
#           I_flux_mJy[0], I_err_mJy[0]        -- MFS Stokes I +/- error
#           alpha[0], alpha_err[0]             -- spectral index (None if not fit)
#           chi2_red[0], ndof[0]               -- reduced chi2 of the alpha fit and its dof
#                                                  (None if not fit)
#           time_ctr_mjd[0]                    -- MJD at the centre of the exposure
#           time_ctr_isot[0]                   -- ISOT date/time at the centre of the exposure
#           time_dt[0]                         -- exposure time (s)
#           freq_GHz[0]                        -- MFS central frequency
#           I_RA_deg[0], I_DEC_deg[0]          -- fitted position (degrees)
#           LP_frac[0], LP_frac_err[0]         -- image-plane ("calculated") pol fraction
#
#   RESULTS/img_<target_ms>_pcalmask_component<N>_rmsynth_RMclean.json
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
#   RESULTS/img_<target_ms>_pcalmask_component<N>_rmsynth_RMsynth.json
#   (sibling of the RMclean.json above, before RM-CLEAN; carries keys that
#    RMclean.json does not)
#       freq0_Hz                               -- central frequency used for the RM fit
#       Ifreq0                                  -- Stokes I at freq0
#
#   RESULTS/img_*_diagnostic_rmsynth_RMclean.json / _RMsynth.json
#       Same keys as above, but for the polarization-angle calibrator's own
#       diagnostic RM synthesis (written by RMSYNTH_01B_systematics.py when
#       cfg.CAL_1GC_DIAGNOSTICS is set). Used as a calibrator sanity check,
#       printed if found.
#
# Position (hmsdms) is derived from I_RA_deg/I_DEC_deg above via
# astropy.coordinates.SkyCoord, not read from a file.

import glob
import json
import os.path as o
import sys

from astropy.coordinates import SkyCoord
from astropy import units as u

sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import generate_jobs as gen
from oxkat import config as cfg


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


def write_rmsynth_block(line, rmclean, rmsynth, include_snr=True):
    line('RM (rad/m^2)', fmt(rmclean.get('phiPeakPIfit_rm2'), rmclean.get('dPhiObserved_rm2')))
    line('Polarization angle (deg)', fmt(rmclean.get('polAngleFit_deg'), rmclean.get('dPolAngleFitObserved_deg')))
    if include_snr:
        line('SNR of rmsynth', fmt(rmclean.get('snrPIfit')))
    line('Central frequency for rmsynth (GHz)', fmt(hz_to_ghz(rmsynth.get('freq0_Hz'))))
    line('Stokes I at freq0', fmt(rmsynth.get('Ifreq0')))
    line('Pol. fraction, representative (AMP_eff/I0, %)', fmt(representative_frac_pol(rmclean, rmsynth)))


def main():

    with open('project_info.json') as f:
        project_info = json.load(f)

    target_names  = project_info.get('target_names', [])
    working_names = project_info.get('working_names', [])
    single_targets = [t for t in target_names if t in working_names]

    if len(single_targets) != 1:
        print(gen.col('Single Target Summary')+f'Skipping: found {len(single_targets)} target(s) in working_names, expected 1')
        return

    target  = single_targets[0]
    myms    = project_info['target_ms'][target_names.index(target)]

    pol_json_path = o.join(cfg.RESULTS, f'img_{myms}_pcalmask_polarization.json')
    with open(pol_json_path) as f:
        pol_data = json.load(f)
    mfs = pol_data.get('MFS', {})

    rmclean_paths = sorted(glob.glob(o.join(cfg.RESULTS, f'img_{myms}_pcalmask_component*_rmsynth_RMclean.json')))

    out_path = f'{myms}_summary.txt'
    with open(out_path, 'w') as out:
        out.write(f'Single-target summary for {target} ({myms})\n')
        out.write('='*60+'\n')

        def line(label, value):
            out.write(f'{label.ljust(40)}: {value}\n')

        for rmclean_path in rmclean_paths:
            component = o.basename(rmclean_path).split(f'img_{myms}_pcalmask_')[1].split('_rmsynth_RMclean.json')[0]

            with open(rmclean_path) as f:
                rmclean = json.load(f)
            rmsynth_path = rmclean_path.replace('_RMclean.json', '_RMsynth.json')
            rmsynth = {}
            if o.isfile(rmsynth_path):
                with open(rmsynth_path) as f:
                    rmsynth = json.load(f)

            comp_mfs = mfs.get(component, {})

            ra_deg = comp_mfs.get('I_RA_deg', [None])[0]
            if ra_deg is not None and ra_deg < 0:
                ra_deg += 360.0
            dec_deg = comp_mfs.get('I_DEC_deg', [None])[0]

            hmsdms = None
            if ra_deg is not None and dec_deg is not None:
                hmsdms = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg).to_string('hmsdms')

            time_ctr_mjd  = comp_mfs.get('time_ctr_mjd', [None])[0]
            time_ctr_isot = comp_mfs.get('time_ctr_isot', [None])[0]

            out.write(f'\n--- {component} ---\n')
            line('MFS Stokes I (mJy)', fmt(comp_mfs.get("I_flux_mJy", [None])[0], comp_mfs.get("I_err_mJy", [None])[0]))
            if 'alpha' in comp_mfs:
                line('Spectral index', fmt(comp_mfs.get("alpha", [None])[0], comp_mfs.get("alpha_err", [None])[0]))
                chi2_red = comp_mfs.get("chi2_red", [None])[0]
                ndof     = comp_mfs.get("ndof", [None])[0]
                if chi2_red is not None:
                    line('Spectral index fit, chi2_red (ndof)', f'{fmt(chi2_red)} ({fmt(ndof)})')
            else:
                line('Spectral index', 'N/A')
            line('Date/time, centre (ISOT, UTC)', fmt(time_ctr_isot))
            line('MJD, centre', fmt(time_ctr_mjd))
            line('Exposure time (s)', fmt(comp_mfs.get("time_dt", [None])[0]))
            line('Central frequency (GHz)', fmt(comp_mfs.get("freq_GHz", [None])[0]))
            line('Pol. fraction, calculated (%)', fmt(comp_mfs.get("LP_frac", [None])[0], comp_mfs.get("LP_frac_err", [None])[0]))
            write_rmsynth_block(line, rmclean, rmsynth)
            line('Position (deg)', f'{fmt(ra_deg)}, {fmt(dec_deg)}')
            line('Position (hmsdms)', hmsdms if hmsdms is not None else 'N/A')

        # Polarization-angle calibrator diagnostic cross-check, if present
        diag_rmclean_paths = sorted(glob.glob(o.join(cfg.RESULTS, 'img_*_diagnostic_rmsynth_RMclean.json')))
        if diag_rmclean_paths:
            diag_rmclean_path = diag_rmclean_paths[0]
            with open(diag_rmclean_path) as f:
                diag_rmclean = json.load(f)
            diag_rmsynth_path = diag_rmclean_path.replace('_RMclean.json', '_RMsynth.json')
            diag_rmsynth = {}
            if o.isfile(diag_rmsynth_path):
                with open(diag_rmsynth_path) as f:
                    diag_rmsynth = json.load(f)

            out.write('\n--- polarization-angle calibrator diagnostic ---\n')
            write_rmsynth_block(line, diag_rmclean, diag_rmsynth, include_snr=False)

    print(gen.col('Single Target Summary')+f'Wrote {out_path}')


if __name__ == "__main__":

    main()
