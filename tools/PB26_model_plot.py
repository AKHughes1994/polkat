#!/usr/bin/env python3
"""
Temporary, standalone preview of the two EVPA models in tools/manual_XF_solver.py
(compute_3c286_evpa, compute_3c138_evpa), duplicated here so they can be plotted
without a CASA session. Saves the same two PNGs manual_XF_solver.py would, into
the current working directory. Delete this file once you're done checking them.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def compute_3c286_evpa(freq_ghz):
    freq_ghz = np.asarray(freq_ghz, dtype=float)
    x = np.log10(freq_ghz)
    evpa_deg = np.full_like(freq_ghz, np.nan)

    low_mask = (freq_ghz >= 0.5) & (freq_ghz <= 1.0)
    high_mask = (freq_ghz > 1.0) & (freq_ghz <= 50.0)

    evpa_deg[low_mask] = (26.0 + 57.0 * x[low_mask]
                           + 615.0 * x[low_mask] ** 2
                           + 3790.0 * x[low_mask] ** 3)
    evpa_deg[high_mask] = (26.1 + 17.1 * x[high_mask]
                            - 16.1 * x[high_mask] ** 2
                            + 5.75 * x[high_mask] ** 3)
    return evpa_deg


def compute_3c138_evpa(freq_ghz):
    freq_ghz = np.asarray(freq_ghz, dtype=float)
    x = np.log10(freq_ghz)
    evpa_deg = np.full_like(freq_ghz, np.nan)

    low_mask = (freq_ghz >= 0.5) & (freq_ghz <= 1.0)
    high_mask = (freq_ghz > 1.0) & (freq_ghz <= 4.0)

    evpa_deg[low_mask] = (-21.9 + 71.1 * x[low_mask]
                           - 1435.0 * x[low_mask] ** 2
                           - 11110.0 * x[low_mask] ** 3
                           - 23190.0 * x[low_mask] ** 4)
    evpa_deg[high_mask] = (-22.0 + 95.4 * x[high_mask]
                            - 340.0 * x[high_mask] ** 2
                            + 534.0 * x[high_mask] ** 3
                            - 308.0 * x[high_mask] ** 4)
    return evpa_deg


XF_TARGET_POLANG_3C286 = 33.0   # config default for 3C286/J1331
XF_TARGET_POLANG_3C138 = -11.0  # config default for 3C138/J0521

PLOT_PAD_GHZ = 0.2  # x-axis padding either side of each model's own calc range

for name, func, freq_range, config_val in [
    ('3c286', compute_3c286_evpa, (0.56, 4.0), XF_TARGET_POLANG_3C286),
    ('3c138', compute_3c138_evpa, (0.56, 4.0), XF_TARGET_POLANG_3C138),
]:
    freq_ghz = np.linspace(freq_range[0], freq_range[1], 2000)
    target_polang_array = func(freq_ghz)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(freq_ghz, target_polang_array, 'b-', linewidth=2, label=f'{name.upper()} EVPA model')
    ax.axhline(config_val, color='red', linestyle='--', linewidth=1.5, alpha=0.7,
               label=f'Config value: {config_val}°')
    ax.set_xlim(freq_range[0] - PLOT_PAD_GHZ, freq_range[1] + PLOT_PAD_GHZ)
    ax.set_xlabel('Frequency [GHz]')
    ax.set_ylabel('EVPA [deg]')
    ax.set_title(f'{name.upper()} Frequency-Dependent EVPA Model (Perley, Butler et al. 2026)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{name}_evpa_model.png', dpi=150)
    plt.close(fig)
    print(f'{name}: EVPA range {target_polang_array.min():.2f} to {target_polang_array.max():.2f} deg '
          f'-> saved {name}_evpa_model.png')
