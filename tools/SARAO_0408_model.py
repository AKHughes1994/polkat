import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def casa_flux_model(lnunu0, iref, *args):
    """
    Compute CASA-style flux density model.
    
    The model is of the form:
    S = S_ref * (nu/nu_ref)^(a_0 + a_1*log10(nu/nu_ref) + a_2*log10(nu/nu_ref)^2 + ...)
    
    Parameters
    ----------
    lnunu0 : np.ndarray
        Log10 of the frequency ratio (nu/nu_ref)
    iref : float
        Reference flux density in Jy
    *args : float
        Polynomial coefficients for the spectral index expansion
        
    Returns
    -------
    np.ndarray
        Model flux densities in Jy
    """
    # Compute polynomial exponent as sum of terms
    exponent = np.sum([arg * (lnunu0 ** power) 
                       for power, arg in enumerate(args)], axis=0)
    
    return iref * (10**lnunu0) ** exponent


def fit_flux_model(nu, s, nu0, sigma, sref, order=5):
    """
    Fit a CASA-style flux model to observed flux densities.
    
    Fits a model of the form:
    S = S_ref * (nu/nu_ref)^(a_0 + a_1*log10(nu/nu_ref) + a_2*log10(nu/nu_ref)^2 + ...)
    
    If the requested fit fails, iteratively falls back to lower orders
    until a successful fit is achieved or order 0 is reached.
    
    Parameters
    ----------
    nu : np.ndarray
        Frequencies in Hz
    s : np.ndarray
        Flux densities in Jy
    nu0 : float
        Reference frequency in Hz
    sigma : np.ndarray
        Uncertainties in flux densities (Jy)
    sref : float
        Initial guess for flux density at reference frequency
    order : int, optional
        Desired polynomial order (default: 5)
        order=1: spectral index only
        order=2: spectral index + curvature
        
    Returns
    -------
    list
        [reference_freq, flux_ref, a_0, a_1, a_2, ...]
        CASA-style flux model parameters
    """
    # Initial parameter guess: [S_ref, spectral_index, higher_order_terms...]
    init = [sref, -0.7] + [0] * (order - 1)
    lnunu0 = np.log10(nu / nu0)
    
    # Try fitting at requested order, fall back to lower orders if needed
    for fitorder in range(order, -1, -1):
        try:
            popt, _ = curve_fit(casa_flux_model, lnunu0, s, 
                                p0=init[:fitorder + 1], sigma=sigma)
        except RuntimeError:
            print(f"Warning: Fitting flux model of order {fitorder} failed. "
                  f"Trying lower order fit.")
        else:
            # Pad coefficients to match requested order
            coeffs = np.pad(popt, ((0, order - fitorder),), "constant")
            return [nu0] + coeffs.tolist()
    
    # Fallback: return weighted mean as zeroth-order model
    coeffs = [np.average(s, weights=1./(sigma**2))] + [0] * order
    return [nu0] + coeffs.tolist()


def convert_flux_model(nu=np.linspace(0.9, 4, 200)*1e9, 
                       a=1, b=0, c=0, d=0, 
                       Reffreq=2.7e9):
    """
    Convert from katpoint log-flux model to CASA-style flux model.
    
    Converts from the form:
        log10(S) = a + b*log10(nu/MHz) + c*log10(nu/MHz)^2 + d*log10(nu/MHz)^3
    
    To CASA form:
        S = S_ref * (nu/nu_ref)^(a_0 + a_1*log10(nu/nu_ref) + ...)
    
    Parameters
    ----------
    nu : np.ndarray
        Frequencies in Hz
    a, b, c, d : float
        Log-flux model polynomial coefficients
    Reffreq : float
        Reference frequency for CASA model (Hz)
        
    Returns
    -------
    list
        [reference_freq, flux_ref, a_0, a_1, a_2]
        CASA-style flux model parameters
    """
    MHz = 1e6
    
    # Compute flux densities from log model
    S = 10**(a + b*np.log10(nu/MHz) + 
             c*np.log10(nu/MHz)**2 + 
             d*np.log10(nu/MHz)**3)
    
    # Fit CASA model to these flux densities
    return fit_flux_model(nu, S, Reffreq, 
                          np.ones_like(nu), sref=1, order=3)


# Example: Convert flux model for calibrator J0408-6545
# Parameters from katpoint model (epoch 2016)
# name=0408-65 epoch=2016 ra=04h08m20.4s dec=-65d45m09s
a = -0.9790
b = 3.3662
c = -1.1216
d = 0.0861

# Convert to CASA format
freq_range = np.linspace(0.9, 4, 200) * 1e9  # 0.9-4 GHz (L to S-band)
reffreq, fluxdensity, spix0, spix1, spix2 = convert_flux_model(
    freq_range, a, b, c, d
)

f_cal_alt = 'J0408-6545'

# Generate comparison plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Compute input model (katpoint style)
MHz = 1e6
S_input = 10**(a + b*np.log10(freq_range/MHz) + 
               c*np.log10(freq_range/MHz)**2 + 
               d*np.log10(freq_range/MHz)**3)

# Compute fitted CASA model
lnunu0 = np.log10(freq_range / reffreq)
S_casa = casa_flux_model(lnunu0, fluxdensity, spix0, spix1, spix2)

# Plot 1: Flux density comparison
ax1.plot(freq_range/1e9, S_input, 'b-', linewidth=2, label='Input (katpoint)')
ax1.plot(freq_range/1e9, S_casa, 'r--', linewidth=2, label='Fitted (CASA)')
ax1.set_ylabel('Flux Density (Jy)', fontsize=12)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_title(f'Flux Model Comparison: {f_cal_alt}', fontsize=14, fontweight='bold')

# Plot 2: Residuals (percentage)
residual_pct = 100 * (S_casa - S_input) / S_input
ax2.plot(freq_range/1e9, residual_pct, 'k-', linewidth=1.5)
ax2.axhline(0, color='gray', linestyle='--', linewidth=1)
ax2.set_xlabel('Frequency (GHz)', fontsize=12)
ax2.set_ylabel('Residual (%)', fontsize=12)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('flux_model_comparison.png', dpi=300, bbox_inches='tight')
plt.clf()
plt.close()
print(f"Plot saved as 'flux_model_comparison.png'")

# Print CASA model parameters
print(f"\nCASA Model Parameters for {f_cal_alt}:")
print(f"  Reference frequency: {reffreq/1e9:.3f} GHz")
print(f"  Reference flux:      {fluxdensity:.4f} Jy")
print(f"  Spectral index:      {spix0:.4f}")
print(f"  Curvature:           {spix1:.4f}")
print(f"  2nd order:           {spix2:.4f}")
plt.show()
