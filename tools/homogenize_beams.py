#!/usr/bin/env python3
"""
Advanced beam homogenization tool with intelligent resource management and optimal parallelization.

This script provides:
- Automatic system resource detection and optimization
- Intelligent memory management for large image datasets
- Parallel processing with both multiprocessing and threading approaches
- Enhanced Welzl algorithm implementation
- MFS image creation with chunked median computation
- Comprehensive logging and progress reporting

Author: andrew.hughes@physics.ox.ac.uk, fraser.cowie@physics.ox.ac.uk
Enhanced implementation with advanced parallelization and resource management.
"""

import os
import sys
import glob
import subprocess
import gc
import time
import uuid
import scipy
import matplotlib.pyplot as plt
import numpy as np
import bottleneck as bn
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from scipy.spatial import ConvexHull
import psutil
import multiprocessing as mp
from multiprocessing import shared_memory
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
import threading
from typing import List, Tuple, Optional, Dict, Any, Union

import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def msg(txt: str) -> None:
    """Print timestamped message with guaranteed flush"""
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp + txt, flush=True)


def get_available_memory_gb() -> float:
    """
    Return the available memory in GB, preferring SLURM allocation over
    node-wide psutil figures.

    SLURM environment variables checked (in order):
      SLURM_MEM_PER_NODE  — total MB allocated to this job
      SLURM_MEM_PER_CPU   — MB per CPU; multiplied by allocated CPUs
                            (SLURM_CPUS_PER_TASK x SLURM_NTASKS, not SLURM_CPUS_ON_NODE)
    Falls back to psutil.virtual_memory().available if neither is set.
    """
    slurm_mem_node = os.environ.get('SLURM_MEM_PER_NODE')
    if slurm_mem_node:
        try:
            allocated_gb = int(slurm_mem_node) / 1024.0
            msg(f'SLURM allocation detected: SLURM_MEM_PER_NODE={slurm_mem_node} MB '
                f'-> {allocated_gb:.1f} GB')
            return allocated_gb
        except ValueError:
            pass

    slurm_mem_cpu  = os.environ.get('SLURM_MEM_PER_CPU')
    # Use allocated CPUs (not SLURM_CPUS_ON_NODE which reflects the full node)
    _cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    _ntasks        = os.environ.get('SLURM_NTASKS')
    _alloc_cpus    = (int(_cpus_per_task) * int(_ntasks) if _cpus_per_task and _ntasks
                      else int(_cpus_per_task) if _cpus_per_task
                      else int(_ntasks) if _ntasks
                      else None)
    if slurm_mem_cpu and _alloc_cpus:
        try:
            allocated_gb = int(slurm_mem_cpu) * _alloc_cpus / 1024.0
            msg(f'SLURM allocation detected: SLURM_MEM_PER_CPU={slurm_mem_cpu} MB x '
                f'{_alloc_cpus} allocated CPUs -> {allocated_gb:.1f} GB')
            return allocated_gb
        except ValueError:
            pass

    # No SLURM allocation found — use node-wide available memory
    available_gb = psutil.virtual_memory().available / 1024**3
    msg(f'No SLURM allocation detected — using node-wide available memory: {available_gb:.1f} GB')
    return available_gb

def get_divisors(n: int) -> List[int]:
    """
    Get all divisors of n in descending order for optimal chunk sizing.
    
    Parameters
    ----------
    n : int
        Number to find divisors for
        
    Returns
    -------
    List[int]
        Sorted list of divisors in descending order
    """
    divisors = []
    sqrt_n = int(np.sqrt(n))
    for i in range(1, sqrt_n + 1):
        if n % i == 0:
            divisors.append(n // i)  # Larger divisor first
            if i != n // i:
                divisors.append(i)
    return sorted(divisors, reverse=True)

def get_system_resources(memory_threshold: float = 0.8) -> Dict[str, Any]:
    """
    Comprehensive system resource detection and analysis.
    
    Parameters
    ----------
    memory_threshold : float
        Fraction of available memory to consider as usable (default: 0.8 = 80%)
    
    Returns
    -------
    Dict[str, Any]
        Dictionary containing detailed system resource information
    """
    # CPU information
    cpu_physical = psutil.cpu_count(logical=False)
    cpu_logical = psutil.cpu_count(logical=True)

    # Respect SLURM CPU allocation if present — same priority order as compute_max_workers
    # SLURM_CPUS_ON_NODE intentionally excluded (reflects total node, not allocation)
    import re as _re
    _cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    _ntasks        = os.environ.get('SLURM_NTASKS')
    _job_cpus      = os.environ.get('SLURM_JOB_CPUS_PER_NODE', '')
    _job_cpus_m    = _re.match(r'(\d+)', _job_cpus)

    if _cpus_per_task and _ntasks:
        slurm_cpu_cap = int(_cpus_per_task) * int(_ntasks)
        _cpu_src = f'SLURM_CPUS_PER_TASK={_cpus_per_task} x SLURM_NTASKS={_ntasks}'
    elif _cpus_per_task:
        slurm_cpu_cap = int(_cpus_per_task)
        _cpu_src = f'SLURM_CPUS_PER_TASK={_cpus_per_task}'
    elif _job_cpus_m:
        slurm_cpu_cap = int(_job_cpus_m.group(1))
        _cpu_src = f'SLURM_JOB_CPUS_PER_NODE={_job_cpus} -> {slurm_cpu_cap}'
    elif _ntasks:
        slurm_cpu_cap = int(_ntasks)
        _cpu_src = f'SLURM_NTASKS={_ntasks}'
    else:
        slurm_cpu_cap = None
        _cpu_src = 'os.cpu_count() (no SLURM allocation detected)'

    # Cap physical/logical counts to SLURM allocation if present
    if slurm_cpu_cap is not None:
        cpu_physical = min(cpu_physical, slurm_cpu_cap)
        cpu_logical  = min(cpu_logical,  slurm_cpu_cap)
        msg(f'SLURM CPU cap applied: {slurm_cpu_cap} CPUs [{_cpu_src}]')
    else:
        msg(f'No SLURM CPU cap [{_cpu_src}]')
    
    # Memory information — respect SLURM allocation if present
    memory = psutil.virtual_memory()
    total_memory_gb = memory.total / (1024**3)
    available_memory_gb = get_available_memory_gb()
    
    # Calculate usable memory with single threshold
    usable_memory_gb = available_memory_gb * memory_threshold
    usable_memory_pixels = int(usable_memory_gb * 1024**3 / 4)  # 4 bytes per float32
    
    # Disk information (for I/O optimization)
    disk_usage = psutil.disk_usage('.')
    free_disk_gb = disk_usage.free / (1024**3)
    
    # Calculate optimal threading parameters
    optimal_threads_io = min(12, cpu_logical)  # I/O bound tasks
    optimal_processes_cpu = max(1, min(cpu_physical, int(available_memory_gb / 2)))  # CPU bound tasks
    
    msg("System Resource Analysis:")
    msg(f"  Physical CPU cores: {cpu_physical}")
    msg(f"  Logical CPU cores: {cpu_logical}")
    msg(f"  Total memory: {total_memory_gb:.1f} GB")
    msg(f"  Available memory: {available_memory_gb:.1f} GB")
    msg(f"  Usable memory ({int(memory_threshold*100)}%): {usable_memory_gb:.1f} GB")
    msg(f"  Usable memory pixels: {usable_memory_pixels:,}")
    msg(f"  Free disk space: {free_disk_gb:.1f} GB")
    msg(f"  Optimal I/O threads: {optimal_threads_io}")
    msg(f"  Optimal CPU processes: {optimal_processes_cpu}")
    
    return {
        'cpu_physical': cpu_physical,
        'cpu_logical': cpu_logical,
        'memory_total_gb': total_memory_gb,
        'memory_available_gb': available_memory_gb,
        'usable_memory_gb': usable_memory_gb,
        'usable_memory_pixels': usable_memory_pixels,
        'free_disk_gb': free_disk_gb,
        'optimal_threads_io': optimal_threads_io,
        'optimal_processes_cpu': optimal_processes_cpu
    }

def filter_beam_outliers(images: List[str], sigma_threshold: float = 3.0) -> Tuple[List[str], List[str]]:
    """
    Filter out PSF files with abnormal beam sizes using iterative sigma clipping until convergence.
    Applied to both BMAJ and BMIN with the same threshold.
    
    Parameters
    ----------
    images : List[str]
        List of PSF file paths
    sigma_threshold : float
        Sigma threshold for outlier rejection (default: 3.0)
        
    Returns
    -------
    Tuple[List[str], List[str]]
        (good_files, rejected_files)
    """
    beam_info = []
    
    msg(f"Analyzing {len(images)} PSF files for beam outliers")
    
    for im in sorted(images):
        try:
            header = fits.getheader(im)
            bmaj = header.get('BMAJ', 0)
            bmin = header.get('BMIN', 0)
            bpa = header.get('BPA', 0)
            freq = header.get('CRVAL3', 0)
            
            if bmaj > 1e-12 and bmin > 1e-12 and freq > 0:
                beam_info.append({
                    'filename': im,
                    'bmaj': bmaj,
                    'bmin': bmin,
                    'bpa': bpa,
                    'freq': freq,
                    'bmaj_normalized': bmaj * freq / 1e9,  # Normalize by GHz
                    'bmin_normalized': bmin * freq / 1e9   # Normalize by GHz
                })
        except Exception as e:
            msg(f"Warning: Could not read beam from {im}: {e}")
            continue
    
    if not beam_info:
        msg("Warning: No valid beam parameters found")
        return [], [f for f in images]
    
    # Iterative sigma clipping until convergence
    current_beam_info = beam_info.copy()
    iteration = 0
    max_iterations = 20  # Safety limit
    
    msg(f"Starting iterative sigma clipping with {sigma_threshold}σ threshold for both BMAJ and BMIN")
    
    while iteration < max_iterations:
        iteration += 1
        
        # Calculate statistics on current set for both BMAJ and BMIN
        normalized_bmaj = np.array([info['bmaj_normalized'] for info in current_beam_info])
        normalized_bmin = np.array([info['bmin_normalized'] for info in current_beam_info])
        
        median_bmaj = np.nanmedian(normalized_bmaj)
        median_bmin = np.nanmedian(normalized_bmin)
        std_bmaj = np.std(normalized_bmaj)
        std_bmin = np.std(normalized_bmin)
        
        if std_bmaj == 0 and std_bmin == 0:
            msg(f"Iteration {iteration}: No variation in beam sizes, stopping")
            break
        
        # Identify outliers in current iteration - reject if EITHER BMAJ OR BMIN is an outlier
        outliers = []
        remaining = []
        
        for info in current_beam_info:
            bmaj_deviation = abs(info['bmaj_normalized'] - median_bmaj) / std_bmaj if std_bmaj > 0 else 0
            bmin_deviation = abs(info['bmin_normalized'] - median_bmin) / std_bmin if std_bmin > 0 else 0
            
            # Reject if either BMAJ or BMIN exceeds threshold
            if bmaj_deviation > sigma_threshold or bmin_deviation > sigma_threshold:
                outliers.append(info)
                info['rejection_reason'] = f"BMAJ: {bmaj_deviation:.1f}σ, BMIN: {bmin_deviation:.1f}σ"
            else:
                remaining.append(info)
        
        msg(f"Iteration {iteration}: BMAJ median={median_bmaj:.3f}, std={std_bmaj:.3f}")
        msg(f"                BMIN median={median_bmin:.3f}, std={std_bmin:.3f}")
        msg(f"                rejected {len(outliers)} outliers, {len(remaining)} remaining")
        
        # Check for convergence (no outliers found)
        if len(outliers) == 0:
            msg(f"Sigma clipping converged after {iteration} iterations")
            break
            
        # Check if we'd reject everything
        if len(remaining) == 0:
            msg(f"Warning: All files would be rejected, stopping at iteration {iteration-1}")
            break
            
        # Update for next iteration
        current_beam_info = remaining
    
    if iteration >= max_iterations:
        msg(f"Warning: Sigma clipping stopped at maximum iterations ({max_iterations})")
    
    # Final statistics
    if current_beam_info:
        final_bmaj = np.array([info['bmaj_normalized'] for info in current_beam_info])
        final_bmin = np.array([info['bmin_normalized'] for info in current_beam_info])
        final_median_bmaj = np.nanmedian(final_bmaj)
        final_median_bmin = np.nanmedian(final_bmin)
        final_std_bmaj = np.std(final_bmaj)
        final_std_bmin = np.std(final_bmin)
        
        msg(f"Final beam statistics:")
        msg(f"  Final median normalized BMAJ: {final_median_bmaj:.3f} arcsec·GHz")
        msg(f"  Final std normalized BMAJ: {final_std_bmaj:.3f} arcsec·GHz")
        msg(f"  Final median normalized BMIN: {final_median_bmin:.3f} arcsec·GHz")
        msg(f"  Final std normalized BMIN: {final_std_bmin:.3f} arcsec·GHz")
        msg(f"  Convergence achieved: {iteration < max_iterations}")
    
    # Separate good and rejected files
    good_filenames = {info['filename'] for info in current_beam_info}
    good_files = [info['filename'] for info in current_beam_info]
    rejected_files = [info['filename'] for info in beam_info if info['filename'] not in good_filenames]
    
    # Create detailed rejection info
    rejected_info = []
    if rejected_files:
        final_median_bmaj = np.nanmedian([info['bmaj_normalized'] for info in current_beam_info]) if current_beam_info else 0
        final_median_bmin = np.nanmedian([info['bmin_normalized'] for info in current_beam_info]) if current_beam_info else 0
        final_std_bmaj = np.std([info['bmaj_normalized'] for info in current_beam_info]) if current_beam_info else 1
        final_std_bmin = np.std([info['bmin_normalized'] for info in current_beam_info]) if current_beam_info else 1
        
        for info in beam_info:
            if info['filename'] in rejected_files:
                bmaj_deviation = abs(info['bmaj_normalized'] - final_median_bmaj) / final_std_bmaj if final_std_bmaj > 0 else 0
                bmin_deviation = abs(info['bmin_normalized'] - final_median_bmin) / final_std_bmin if final_std_bmin > 0 else 0
                rejected_info.append({
                    'filename': info['filename'],
                    'bmaj_normalized': info['bmaj_normalized'],
                    'bmin_normalized': info['bmin_normalized'],
                    'bmaj_deviation_sigma': bmaj_deviation,
                    'bmin_deviation_sigma': bmin_deviation,
                    'freq_ghz': info['freq'] / 1e9,
                    'bmaj': info['bmaj'],
                    'bmin': info['bmin'],
                    'bpa': info['bpa']
                })
    
    # Report results
    msg(f"Sigma clipping results:")
    msg(f"  Original files: {len(beam_info)}")
    msg(f"  Good files: {len(good_files)}")
    msg(f"  Rejected files: {len(rejected_files)}")
    msg(f"  Rejection rate: {len(rejected_files)/len(beam_info)*100:.1f}%")
    
    if rejected_info:
        msg(f"  Final rejected outliers (>{sigma_threshold}σ in BMAJ or BMIN):")
        for rej in rejected_info:
            msg(f"    {os.path.basename(rej['filename'])}: "
                f"BMAJ={rej['bmaj_deviation_sigma']:.1f}σ, "
                f"BMIN={rej['bmin_deviation_sigma']:.1f}σ "
                f"({rej['freq_ghz']:.2f} GHz)")
        
        # Save detailed rejected files list to text file
        rejected_filename = "rejected_psf_files.txt"
        try:
            with open(rejected_filename, 'w') as f:
                f.write("# Rejected PSF files from iterative sigma clipping (BMAJ and BMIN)\n")
                f.write(f"# Sigma threshold: {sigma_threshold}\n")
                f.write(f"# Converged after {iteration} iterations\n")
                if current_beam_info:
                    f.write(f"# Final median normalized BMAJ: {final_median_bmaj:.3f} arcsec·GHz\n")
                    f.write(f"# Final std normalized BMAJ: {final_std_bmaj:.3f} arcsec·GHz\n")
                    f.write(f"# Final median normalized BMIN: {final_median_bmin:.3f} arcsec·GHz\n")
                    f.write(f"# Final std normalized BMIN: {final_std_bmin:.3f} arcsec·GHz\n")
                f.write("#\n")
                f.write("# Columns: filename, bmaj(arcsec), bmin(arcsec), bpa(deg), freq(GHz), bmaj_norm(arcsec·GHz), bmin_norm(arcsec·GHz), bmaj_dev(sigma), bmin_dev(sigma)\n")
                f.write("#\n")
                
                for rej in rejected_info:
                    f.write(f"{rej['filename']:<80} "
                           f"{rej['bmaj']:<10.6f} "
                           f"{rej['bmin']:<10.6f} "
                           f"{rej['bpa']:<8.2f} "
                           f"{rej['freq_ghz']:<8.3f} "
                           f"{rej['bmaj_normalized']:<12.3f} "
                           f"{rej['bmin_normalized']:<12.3f} "
                           f"{rej['bmaj_deviation_sigma']:<8.2f} "
                           f"{rej['bmin_deviation_sigma']:<8.2f}\n")
            
            msg(f"  Saved detailed rejection list to: {rejected_filename}")
            
        except Exception as e:
            msg(f"Warning: Could not save rejected files list: {e}")
        
        # Also save a simple list of just the filenames for easy use
        simple_rejected_filename = "rejected_psf_filenames.txt"
        try:
            with open(simple_rejected_filename, 'w') as f:
                f.write("# List of rejected PSF filenames from iterative sigma clipping (BMAJ and BMIN)\n")
                f.write(f"# Converged after {iteration} iterations\n")
                for rej in rejected_info:
                    f.write(f"{rej['filename']}\n")
            
            msg(f"  Saved simple filename list to: {simple_rejected_filename}")
            
        except Exception as e:
            msg(f"Warning: Could not save simple rejected files list: {e}")
    
    return good_files, rejected_files

# ============================================================================
# FITS FILE MANIPULATION
# ============================================================================

def get_image(fitsfile: str) -> Optional[np.ndarray]:
    """
    Load FITS image data with proper degenerate axis handling.
    
    Parameters
    ----------
    fitsfile : str
        Path to FITS file
        
    Returns
    -------
    Optional[np.ndarray]
        Image data as 2D array, or None if loading fails
    """
    try:
        with fits.open(fitsfile) as hdul:
            data = hdul[0].data
            # Handle degenerate axes - always use last 2 dimensions
            while len(data.shape) > 2:
                data = data[0]
            return np.array(data, dtype=np.float32)
    except Exception as e:
        msg(f"Error loading {fitsfile}: {e}")
        return None

def create_fits(newimage: np.ndarray, newheader: fits.Header, fitsfile: str) -> None:
    """
    Create a new FITS file with proper 4D structure for radio astronomy.
    
    Parameters
    ----------
    newimage : np.ndarray
        Image data to write
    newheader : fits.Header
        FITS header for the new file
    fitsfile : str
        Output file path
    """
    # Ensure 4D structure for radio astronomy FITS
    if len(newimage.shape) == 2:
        data_to_write = newimage.reshape(1, 1, newimage.shape[0], newimage.shape[1])
    elif len(newimage.shape) == 3:
        data_to_write = newimage.reshape(1, newimage.shape[0], newimage.shape[1], newimage.shape[2])
    else:
        data_to_write = newimage
    
    hdu = fits.PrimaryHDU(data=data_to_write, header=newheader)
    hdu.writeto(fitsfile, overwrite=True)
    msg(f'Created FITS file: {fitsfile}')

def flush_fits(newimage: np.ndarray, newheader: fits.Header, fitsfile: str) -> None:
    """
    Update existing FITS file with new data and header.
    
    Parameters
    ----------
    newimage : np.ndarray
        New image data
    newheader : fits.Header
        New header
    fitsfile : str
        FITS file to update
    """
    with fits.open(fitsfile, mode='update') as f:
        input_hdu = f[0]
        input_hdu.header = newheader
        
        # Handle different dimensionalities
        if len(input_hdu.data.shape) == 2:
            input_hdu.data[:, :] = newimage
        elif len(input_hdu.data.shape) == 3:
            input_hdu.data[0, :, :] = newimage
        else:
            input_hdu.data[0, 0, :, :] = newimage
        f.flush()

# ============================================================================
# WELZL ALGORITHM IMPLEMENTATION
# ============================================================================

def is_singular(A: np.ndarray, tolerance: float = 1e-12) -> bool:
    """
    Check if matrix is close to singular with improved precision.
    
    Parameters
    ----------
    A : np.ndarray
        Matrix to check
    tolerance : float
        Numerical tolerance for singularity detection
        
    Returns
    -------
    bool
        True if matrix is singular or near-singular
    """
    try:
        # Use SVD for singularity detection
        _, s, _ = np.linalg.svd(A)
        return np.min(s) < tolerance * np.max(s)
    except np.linalg.LinAlgError:
        return True

def center_form_to_geometric(F: np.ndarray, c: np.ndarray) -> Optional[Tuple[np.ndarray, float, float, float]]:
    """
    Convert ellipse from center form to geometric parameters with enhanced precision.
    
    Parameters
    ----------
    F : np.ndarray
        Ellipse matrix in center form
    c : np.ndarray
        Ellipse center coordinates
        
    Returns
    -------
    Optional[Tuple[np.ndarray, float, float, float]]
        (center, major_radius, minor_radius, rotation_angle) or None if degenerate
    """
    try:
        # Compute eigenvalues and eigenvectors
        w, V = np.linalg.eigh(F)
        
        # Check for positive definiteness (valid ellipse)
        if np.any(w <= 0):
            return None
        
        # Ensure consistent orientation (sin(t) >= 0)
        if V[1, 0] < 0:
            V[:, 0] = -V[:, 0]
        
        # Compute rotation angle
        rotation_angle = np.arctan2(V[1, 0], V[0, 0])
        
        # Compute semi-axes with safeguards
        semi_major = 1.0 / np.sqrt(np.min(w))
        semi_minor = 1.0 / np.sqrt(np.max(w))
        
        return c, semi_major, semi_minor, rotation_angle
        
    except (np.linalg.LinAlgError, ValueError, ZeroDivisionError):
        return None

def ellipse_from_boundary5(S: np.ndarray) -> Optional[Tuple[np.ndarray, float, float, float]]:
    """
    Compute unique ellipse through 5 boundary points with enhanced numerical stability.
    
    Parameters
    ----------
    S : np.ndarray
        Array of shape (5, 2) containing boundary points
        
    Returns
    -------
    Optional[Tuple[np.ndarray, float, float, float]]
        Ellipse parameters or None if computation fails
    """
    if S.shape != (5, 2):
        return None
    
    try:
        x, y = S[:, 0], S[:, 1]
        A = np.column_stack((x**2, y**2, 2 * x * y, x, y))
        
        if is_singular(A):
            return None
        
        # Solve for ellipse coefficients
        sol = np.linalg.solve(A, -np.ones(S.shape[0]))
        
        # Compute ellipse center
        center_matrix = -2 * np.array([[sol[0], sol[2]], [sol[2], sol[1]]])
        if is_singular(center_matrix):
            return None
            
        c = np.linalg.solve(center_matrix, sol[3:5])
        
        # Construct ellipse matrix in center form
        A_ext = np.vstack([
            np.hstack([np.eye(3), -np.array([[sol[0], sol[2], sol[1]]]).T]),
            np.array([c[0]**2, 2 * c[0] * c[1], c[1]**2, -1])
        ])
        
        s = np.linalg.solve(A_ext, np.array([0, 0, 0, 1]))
        F = np.array([[s[0], s[1]], [s[1], s[2]]])
        
        return center_form_to_geometric(F, c)
        
    except (np.linalg.LinAlgError, ValueError):
        return None

def ellipse_from_boundary4(S: np.ndarray) -> Optional[Tuple[np.ndarray, float, float, float]]:
    """
    Compute smallest ellipse through 4 boundary points with enhanced precision.
    
    Parameters
    ----------
    S : np.ndarray
        Array of shape (4, 2) containing boundary points
        
    Returns
    -------
    Optional[Tuple[np.ndarray, float, float, float]]
        Ellipse parameters or None if computation fails
    """
    if S.shape != (4, 2):
        return None
    
    try:
        # Sort points in clockwise order
        Sc = S - np.mean(S, axis=0)
        angles = np.arctan2(Sc[:, 1], Sc[:, 0])
        S = S[np.argsort(-angles), :]
        
        # Find diagonal intersection
        A = np.column_stack([S[2, :] - S[0, :], S[1, :] - S[3, :]])
        if is_singular(A):
            return None
        
        b = S[1, :] - S[0, :]
        s = np.linalg.solve(A, b)
        diag_intersect = S[0, :] + s[0] * (S[2, :] - S[0, :])
        
        # Transform to origin
        S = S - diag_intersect
        
        # Rotation transformation
        AC = S[2, :] - S[0, :]
        theta = np.arctan2(AC[1], AC[0])
        cos_theta, sin_theta = np.cos(theta), np.sin(theta)
        rot_mat = np.array([[cos_theta, sin_theta], [-sin_theta, cos_theta]])
        S = rot_mat @ S.T
        S = S.T
        
        # Shear transformation with numerical safeguards
        denominator = S[3, 1] - S[1, 1]
        if abs(denominator) < 1e-12:
            return None
        m = (S[1, 0] - S[3, 0]) / denominator
        shear_mat = np.array([[1.0, m], [0.0, 1.0]])
        S = shear_mat @ S.T
        S = S.T
        
        # Make quadrilateral cyclic
        b = np.linalg.norm(S, axis=1)
        if np.any(b == 0):
            return None
        d = (b[1] * b[3]) / (b[2] * b[0])
        if d <= 0:
            return None
        
        d_quarter = d ** 0.25
        stretch_mat = np.diag([d_quarter, 1.0 / d_quarter])
        S = stretch_mat @ S.T
        S = S.T
        
        # Solve cubic equation for optimal swing angle
        a = np.linalg.norm(S, axis=1)
        coeff = np.zeros(4)
        coeff[0] = -4 * a[1]**2 * a[2] * a[0]
        coeff[1] = -4 * a[1] * (a[2] - a[0]) * (a[1]**2 - a[2] * a[0])
        coeff[2] = (3 * a[1]**2 * (a[1]**2 + a[2]**2) - 
                   8 * a[1]**2 * a[2] * a[0] + 
                   3 * (a[1]**2 + a[2]**2) * a[0]**2)
        coeff[3] = coeff[1] / 2.0
        
        roots = np.roots(coeff)
        valid_roots = roots[(-1 < np.real(roots)) & (np.real(roots) < 1) & (np.abs(np.imag(roots)) < 1e-10)]
        
        if len(valid_roots) == 0:
            return None
        
        theta = np.arcsin(np.real(valid_roots[0]))
        
        # Apply transformation D_theta
        cos_theta_sqrt = np.cos(theta) ** 0.5
        if cos_theta_sqrt == 0:
            return None
        D_mat = np.array([[1.0 / cos_theta_sqrt, np.sin(theta) / cos_theta_sqrt],
                         [0.0, cos_theta_sqrt]])
        S = D_mat @ S.T
        S = S.T
        
        # Find enclosing circle
        boundary = S[:-1, :]
        A = np.vstack([-2 * boundary.T, np.ones(boundary.shape[0])]).T
        if is_singular(A):
            return None
        
        b = -np.sum(boundary**2, axis=1)
        s = np.linalg.solve(A, b)
        
        circle_c = s[:2]
        circle_r_sq = np.sum(circle_c**2) - s[2]
        if circle_r_sq <= 0:
            return None
        circle_r = np.sqrt(circle_r_sq)
        
        # Compute total transformation
        T_mat = D_mat @ stretch_mat @ shear_mat @ rot_mat
        
        # Find original ellipse parameters
        ellipse_c = np.linalg.solve(T_mat, circle_c) + diag_intersect
        ellipse_F = T_mat.T @ T_mat / (circle_r**2)
        
        return center_form_to_geometric(ellipse_F, ellipse_c)
        
    except (np.linalg.LinAlgError, ValueError, ZeroDivisionError):
        return None

def ellipse_from_boundary3(S: np.ndarray) -> Optional[Tuple[np.ndarray, float, float, float]]:
    """
    Compute smallest ellipse through 3 boundary points with enhanced stability.
    
    Parameters
    ----------
    S : np.ndarray
        Array of shape (3, 2) containing boundary points
        
    Returns
    -------
    Optional[Tuple[np.ndarray, float, float, float]]
        Ellipse parameters or None if computation fails
    """
    if S.shape != (3, 2):
        return None
    
    try:
        c = np.mean(S, axis=0)
        Sc = S - c
        
        if is_singular(Sc):
            return None
        
        # Matrix inversion with condition check
        StS = Sc.T @ Sc
        if is_singular(StS):
            return None
            
        F = 1.5 * np.linalg.inv(StS)
        return center_form_to_geometric(F, c)
        
    except (np.linalg.LinAlgError, ValueError):
        return None

def welzl(interior: np.ndarray, boundary: np.ndarray = None) -> Optional[Tuple[np.ndarray, float, float, float]]:
    """
    Enhanced Welzl algorithm for minimum enclosing ellipse with improved robustness.
    
    Parameters
    ----------
    interior : np.ndarray
        Points to be contained within the ellipse
    boundary : np.ndarray, optional
        Points that must lie on the ellipse boundary
        
    Returns
    -------
    Optional[Tuple[np.ndarray, float, float, float]]
        Ellipse parameters (center, major_radius, minor_radius, rotation) or None
    """
    if boundary is None:
        boundary = np.zeros((0, 2))
    
    # Check stopping conditions
    if interior.shape[0] == 0 or boundary.shape[0] == 5:
        if boundary.shape[0] <= 2:
            return None
        elif boundary.shape[0] == 3:
            return ellipse_from_boundary3(boundary)
        elif boundary.shape[0] == 4:
            return ellipse_from_boundary4(boundary)
        else:
            return ellipse_from_boundary5(boundary)
    
    # Random point selection with deterministic fallback
    try:
        i = np.random.randint(interior.shape[0])
    except ValueError:
        return None
    
    p = interior[i, :]
    interior_wo_p = np.delete(interior, i, 0)
    
    # Recursive call without point p
    ellipse = welzl(interior_wo_p, boundary)
    
    # Check if point p is contained in the ellipse
    if is_in_ellipse(p, ellipse):
        return ellipse
    else:
        # Point p must be on the boundary
        new_boundary = np.vstack([boundary, p]) if boundary.size > 0 else p.reshape(1, -1)
        return welzl(interior_wo_p, new_boundary)

def is_in_ellipse(point: np.ndarray, ellipse: Optional[Tuple]) -> bool:
    """
    Check if point is contained within ellipse with numerical tolerance.
    
    Parameters
    ----------
    point : np.ndarray
        Point coordinates
    ellipse : Optional[Tuple]
        Ellipse parameters or None
        
    Returns
    -------
    bool
        True if point is inside ellipse
    """
    if ellipse is None or point is None:
        return False
    
    try:
        c, a, b, t = ellipse
        v = point - c
        
        cos_t, sin_t = np.cos(t), np.sin(t)
        rot_mat = np.array([[cos_t, sin_t], [-sin_t, cos_t]])
        F = rot_mat.T @ np.diag([1.0/(a**2), 1.0/(b**2)]) @ rot_mat
        
        return v.T @ F @ v <= 1.0 + 1e-10  # Small tolerance for numerical precision
        
    except (ValueError, TypeError, IndexError):
        return False

def sample_ellipse(ellipse: Tuple, num_pts: int = 100, endpoint: bool = True) -> np.ndarray:
    """
    Sample points uniformly on ellipse boundary with enhanced precision.
    
    Parameters
    ----------
    ellipse : Tuple
        Ellipse parameters (major_radius, minor_radius, rotation_angle)
    num_pts : int
        Number of points to sample
    endpoint : bool
        Whether to include endpoint for closed curves
        
    Returns
    -------
    np.ndarray
        Sampled points on ellipse
    """
    a, b, p = ellipse
    
    cos_p, sin_p = np.cos(p), np.sin(p)
    rot_mat = np.array([[cos_p, -sin_p], [sin_p, cos_p]])
    
    theta = np.linspace(0, 2 * np.pi, num_pts, endpoint=endpoint)
    z = np.column_stack((a * np.cos(theta), b * np.sin(theta)))
    
    return (rot_mat @ z.T).T

def plot_ellipse(ellipse: Optional[Tuple], num_pts: int = 100, style: str = "-") -> None:
    """
    Plot ellipse with enhanced visualization options.
    
    Parameters
    ----------
    ellipse : Optional[Tuple]
        Ellipse parameters or None
    num_pts : int
        Number of points for plotting
    style : str
        Plot style string
    """
    if ellipse is None:
        return
    
    try:
        x = sample_ellipse(ellipse, num_pts)
        plt.plot(x[:, 0], x[:, 1], style, label='Minimum enclosing ellipse', 
                zorder=100000, linewidth=2)
    except Exception as e:
        msg(f"Warning: Could not plot ellipse - {e}")

# ============================================================================
# CONVEX HULL AND BEAM ANALYSIS
# ============================================================================

def convexhull(images: List[str], sigma_threshold: float = 3.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate convex hull from beam parameters with outlier filtering.
    
    Parameters
    ----------
    images : List[str]
        List of image file paths
    sigma_threshold : float
        Sigma threshold for outlier rejection
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        (hull_points, all_points) for ellipse fitting
    """
    # Filter out beam outliers first
    good_images, rejected_images = filter_beam_outliers(images, sigma_threshold)
    
    if not good_images:
        raise ValueError("All PSF files rejected as outliers - check sigma threshold or data quality")
    
    msg(f"Using {len(good_images)} PSF files for convex hull computation")
    
    beam_params = []
    
    for im in sorted(good_images):
        try:
            header = fits.getheader(im)
            bmaj = header.get('BMAJ', 0)
            bmin = header.get('BMIN', 0)
            bpa = header.get('BPA', 0)
            
            # Convert beam parameters
            a = bmaj / (8 * np.log(2)) ** 0.5
            b = bmin / (8 * np.log(2)) ** 0.5
            p = np.radians(bpa) + np.pi * 0.5
            beam_params.append([a, b, p])
            
        except Exception as e:
            msg(f"Warning: Could not read beam from {im}: {e}")
            continue
    
    if not beam_params:
        raise ValueError("No valid beam parameters found in filtered images")
    
    # Generate ellipse points for each beam
    all_points = []
    for a, b, p in beam_params:
        try:
            points = sample_ellipse((a, b, p), endpoint=False)
            all_points.append(points)
        except Exception as e:
            msg(f"Warning: Could not sample ellipse for beam {a}, {b}, {p}: {e}")
            continue
    
    if not all_points:
        raise ValueError("Could not generate ellipse points from beam parameters")
    
    all_points = np.vstack(all_points)
    
    try:
        hull = ConvexHull(all_points)  # ← This should be the actual ConvexHull computation
        hull_points = all_points[hull.vertices]
        return hull_points, all_points
    except Exception as e:
        msg(f"Error computing convex hull: {e}")
        raise

def get_homogenized_beam(identifier: str, sigma_threshold: float = 3.0) -> Tuple[float, float, float]:    
    """
    Compute homogenized beam parameters using enhanced Welzl algorithm.
    
    Parameters
    ----------
    identifier : str
        Image identifier pattern
        
    Returns
    -------
    Tuple[float, float, float]
        (semi_major, semi_minor, position_angle) of homogenized beam
    """
    images = glob.glob(f'{identifier}-*[!MFS]-psf.fits')
    
    if not images:
        raise FileNotFoundError(f"No PSF files found for identifier: {identifier}")
    
    msg(f"Found {len(images)} PSF files for beam analysis")
    
    # Get convex hull points
    hull_points, all_points = convexhull(images, sigma_threshold)
    msg(f"Computed convex hull with {len(hull_points)} vertices from {len(all_points)} beam points")
    
    # Apply Welzl algorithm with multiple attempts
    ellipse = None
    max_attempts = 5
    
    for attempt in range(max_attempts):
        try:
            ellipse = welzl(hull_points)
            if ellipse is not None:
                break
        except Exception as e:
            msg(f"Welzl attempt {attempt + 1} failed: {e}")
            if attempt < max_attempts - 1:
                np.random.seed(attempt + 42)  # Different seed for each attempt
    
    if ellipse is None:
        raise RuntimeError("Failed to compute minimum enclosing ellipse after multiple attempts")
    
    # Create diagnostic plot
    try:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.set_ylabel('Dec offset (degrees)', fontsize=12)
        ax.set_xlabel('RA offset (degrees)', fontsize=12)
        ax.scatter(hull_points[:, 0], hull_points[:, 1], color='red', 
                  label='Convex hull points', zorder=10000, marker='o', s=16, alpha=0.8)
        ax.scatter(all_points[:, 0], all_points[:, 1], color='green', 
                  label='All beam-sampled points', marker='s', s=8, alpha=0.5)
        
        # Set appropriate axis limits
        pad_factor = 1.3
        ax.set_xlim(np.min(hull_points[:, 0]) * pad_factor, np.max(hull_points[:, 0]) * pad_factor)
        ax.set_ylim(np.min(hull_points[:, 1]) * pad_factor, np.max(hull_points[:, 1]) * pad_factor)
        
        plot_ellipse(ellipse[1:], style="k--")
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        plot_name = f'{identifier.replace("IMAGES/", "").replace("INTERVALS/", "")}_ellipse.png'
        plt.tight_layout()
        plt.savefig(plot_name, dpi=150, bbox_inches='tight')
        plt.close()
        msg(f"Saved beam ellipse plot: {plot_name}")
        
    except Exception as e:
        msg(f"Warning: Could not create ellipse plot: {e}")
    
    center, semi_major, semi_minor, rotation = ellipse
    msg(f'Ellipse center: ({center[0]:.2e}, {center[1]:.2e}) - should be near zero')
    msg(f'Semi-major axis: {semi_major:.2e}, Semi-minor axis: {semi_minor:.2e}')
    msg(f'Rotation angle: {np.degrees(rotation):.2f} degrees')
    
    return semi_major, semi_minor, rotation

# ============================================================================
# MEMORY MANAGEMENT AND CHUNK PROCESSING
# ============================================================================

def estimate_image_memory_usage(images: List[str], safety_factor: float = 1.0) -> Tuple[float, float]:
    """
    Estimate total memory required to load all images simultaneously.
    
    Parameters
    ----------
    images : List[str]
        List of image file paths
    safety_factor : float
        Safety multiplier for memory estimation (default: 1.0, can be increased for safety margin)
        
    Returns
    -------
    Tuple[float, float]
        (estimated_memory_gb, safety_factor) - memory usage in GB and the safety factor used
    """
    try:
        # Sample first few images to estimate size
        sample_size = min(3, len(images))
        total_pixels = 0
        samples_loaded = 0
        
        for i in range(sample_size):
            img = get_image(images[i])
            if img is not None:
                total_pixels += img.size
                samples_loaded += 1
        
        if samples_loaded > 0:
            avg_pixels_per_image = total_pixels / samples_loaded
            total_estimated_pixels = avg_pixels_per_image * len(images)
            # 4 bytes per float32 pixel
            estimated_gb = (total_estimated_pixels * 4 * safety_factor) / (1024**3)
            return estimated_gb, safety_factor
        else:
            return float('inf'), safety_factor  # Could not estimate
            
    except Exception as e:
        msg(f"Warning: Could not estimate memory usage: {e}")
        return float('inf'), safety_factor

def estimate_memory_from_dimensions(num_images: int, ny: int, nx: int, safety_factor: float = 1.5) -> Tuple[float, float]:
    """
    Estimate total memory required based on number of images and dimensions.
    
    Parameters
    ----------
    num_images : int
        Number of images
    ny : int
        Image height in pixels
    nx : int
        Image width in pixels
    safety_factor : float
        Safety multiplier for memory estimation (default: 1.5)
        
    Returns
    -------
    Tuple[float, float]
        (estimated_memory_gb, safety_factor) - memory usage in GB and the safety factor used
    """
    try:
        total_pixels = num_images * ny * nx
        # 4 bytes per float32 pixel
        estimated_gb = (total_pixels * 4 * safety_factor) / (1024**3)
        return estimated_gb, safety_factor
    except Exception as e:
        msg(f"Warning: Could not estimate memory usage: {e}")
        return float('inf'), safety_factor

def can_load_all_images_in_memory(images: List[str], system_info: Dict[str, Any]) -> bool:
    """
    Determine if all images can be safely loaded into memory simultaneously.
    
    Parameters
    ----------
    images : List[str]
        List of image file paths
    system_info : Dict[str, Any]
        System resource information
        
    Returns
    -------
    bool
        True if all images can fit in memory
    """
    estimated_gb, actual_safety_factor = estimate_image_memory_usage(images)
    available_gb = system_info['usable_memory_gb']
        
    msg(f"Memory analysis:")
    msg(f"  Estimated image data: {estimated_gb:.1f} GB (safety factor: {actual_safety_factor}x)")
    msg(f"  Available memory: {available_gb:.1f} GB")
    
    can_fit = estimated_gb <= available_gb
    msg(f"  Can load all images in memory: {can_fit}")
    
    return can_fit

def calculate_optimal_chunk_size(ny: int, nx: int, num_images: int, 
                                resources: Dict[str, Any], n_workers: Optional[int] = None, 
                                worker_type: str = 'auto') -> Tuple[int, int, int]:
    """
    Calculate optimal chunk size for parallel processing based on system resources and workload.
    
    Determines the best spatial chunking strategy to maximize throughput while respecting
    memory constraints. For CPU-bound multiprocessing, emphasizes load balancing across
    cores. For I/O-bound threading, optimizes memory usage and parallel I/O efficiency.
    
    Parameters
    ----------
    ny : int
        Image height in pixels
    nx : int
        Image width in pixels  
    num_images : int
        Number of images to process
    resources : Dict[str, Any]
        System resource information from get_system_resources()
    n_workers : Optional[int]
        Number of workers, auto-detected if None
    worker_type : str
        Processing strategy: 'auto', 'processes', or 'threads'
        
    Returns
    -------
    Tuple[int, int, int]
        (chunk_height, chunk_width, estimated_peak_memory) in pixels and bytes
    
    Notes
    -----
    Auto-detection chooses multiprocessing for datasets that fit in memory,
    threading for larger datasets requiring I/O-bound processing.
    """
    # Auto-detect worker type based on memory constraints if not specified
    if worker_type == 'auto':
        estimated_memory_gb, actual_safety_factor = estimate_memory_from_dimensions(num_images, ny, nx)
        can_fit_in_memory = estimated_memory_gb <= resources['usable_memory_gb']
        worker_type = 'processes' if can_fit_in_memory else 'threads'
        msg(f"Auto-detected worker type: {worker_type} (estimated memory: {estimated_memory_gb:.1f} GB, safety factor: {actual_safety_factor}x, can fit: {can_fit_in_memory})")

    # Set default number of workers based on type
    if n_workers is None:
        if worker_type == 'processes':
            n_workers = resources['optimal_processes_cpu']
        else:
            n_workers = resources['optimal_threads_io']
    
    # Calculate memory constraints based on worker type
    if worker_type == 'processes':
        memory_per_worker = resources['usable_memory_pixels'] / n_workers
        max_chunk_area_per_worker = memory_per_worker
        strategy_desc = "multiprocessing (shared memory)"
    else:
        total_memory = resources['usable_memory_pixels']
        max_chunk_area_per_worker = total_memory / (n_workers * num_images)
        memory_per_worker = total_memory
        strategy_desc = "threading (I/O bound)"    

    msg(f'Calculating optimal chunk size for {ny}x{nx} images')
    msg(f'  Strategy: {strategy_desc}')
    msg(f'  Available memory: {resources["usable_memory_pixels"]:,} pixels')
    msg(f'  Number of workers: {n_workers}')
    msg(f'  Number of images: {num_images}')
    msg(f'  Memory per worker: {memory_per_worker:,} pixels')
    msg(f'  Max chunk area per worker: {max_chunk_area_per_worker:.0f} pixels')

    # Get divisors for spatial chunking
    y_divisors = get_divisors(ny)
    x_divisors = get_divisors(nx)
    
    if worker_type == 'processes':
        # For CPU-bound multiprocessing, aim for optimal work distribution
        total_pixels = ny * nx
        min_pixels_per_chunk = 256 * 256  # Minimum chunk size for process overhead
        max_reasonable_chunks = total_pixels // min_pixels_per_chunk
        
        # Limit workers to available work
        effective_workers = min(n_workers, max(1, max_reasonable_chunks))
        
        if effective_workers < n_workers:
            msg(f'  Reducing workers from {n_workers} to {effective_workers} for efficiency')
            n_workers = effective_workers
        
        # Target chunks for load balancing
        target_chunks = n_workers * 1.5
        target_chunk_area = total_pixels / target_chunks
        
        msg(f'  Target chunks: {target_chunks:.0f}')
        msg(f'  Target chunk area: {target_chunk_area:.0f} pixels')
        
        best_chunk_y, best_chunk_x = ny, nx  # Start with whole image
        best_score = 0
        
        for chunk_y in y_divisors:
            for chunk_x in x_divisors:
                chunk_area = chunk_y * chunk_x
                
                # Must fit in memory
                if chunk_area > max_chunk_area_per_worker:
                    continue
                
                # Must have minimum work per chunk
                if chunk_area < min_pixels_per_chunk:
                    continue
                
                total_chunks = (ny // chunk_y) * (nx // chunk_x)
                chunks_per_worker = total_chunks / n_workers
                
                # Calculate work distribution score
                if chunks_per_worker < 0.5:
                    distribution_score = chunks_per_worker * 2
                elif chunks_per_worker > 3.0:
                    distribution_score = 1.0 / (chunks_per_worker - 2.0)
                else:
                    distribution_score = 1.0
                
                # Prefer reasonable chunk sizes
                size_score = min(1.0, chunk_area / target_chunk_area)
                
                # Prefer square chunks (minor factor)
                aspect_ratio = max(chunk_y, chunk_x) / min(chunk_y, chunk_x)
                square_score = 1.0 / aspect_ratio
                
                # Combined score
                # combined_score = distribution_score * 0.6 + size_score * 0.3 + square_score * 0.1
                combined_score = distribution_score * 1.0

                if combined_score > best_score:
                    best_chunk_y = chunk_y
                    best_chunk_x = chunk_x
                    best_score = combined_score
    
    else:
        # For I/O-bound threading, use existing multi-criteria optimization
        best_chunk_y, best_chunk_x = 64, 64  # Minimum fallback
        best_score = 0
        max_chunk_area = max_chunk_area_per_worker

        for chunk_y in y_divisors:
            for chunk_x in x_divisors:
                chunk_area = chunk_y * chunk_x
                
                # Check bounds
                if chunk_area > max_chunk_area:
                    continue
                
                # Multi-criteria scoring
                total_chunks = (ny // chunk_y) * (nx // chunk_x)
                
                # Size efficiency score (prefer larger chunks)
                size_efficiency = min(1.0, chunk_area / (max_chunk_area_per_worker))
                
                # Aspect ratio score (prefer square chunks)
                aspect_ratio = max(chunk_y, chunk_x) / min(chunk_y, chunk_x)
                square_score = 1.0 / aspect_ratio
                
                # Worker distribution score for I/O parallelism
                chunks_per_worker = total_chunks / n_workers
                distribution_score = min(1.0, chunks_per_worker / 4.0)
                
                # Memory efficiency score
                memory_usage = chunk_area * num_images
                memory_efficiency = 1.0 - (memory_usage / memory_per_worker)
                memory_efficiency = max(0.1, memory_efficiency)
                
                # Combined score
                combined_score = (size_efficiency * 0.2 + 
                                square_score * 0.0 + 
                                distribution_score * 0.4 + 
                                memory_efficiency * 0.4)
                
                if combined_score > best_score:
                    best_chunk_y = chunk_y
                    best_chunk_x = chunk_x
                    best_score = combined_score
    
    # Calculate final statistics
    total_chunks = (ny // best_chunk_y) * (nx // best_chunk_x)
    chunks_per_worker = total_chunks / n_workers
    
    if worker_type == 'processes':
        # For processes with shared memory: image stack in shared memory + chunk processing per worker
        total_image_memory = num_images * ny * nx * 4  # 4 bytes per float32 for entire shared image stack
        chunk_processing_memory = best_chunk_y * best_chunk_x * 4 * n_workers  # Additional memory for chunk processing per worker
        peak_memory_estimate = total_image_memory + chunk_processing_memory
    else:
        # For threads: Each thread loads one chunk worth of data from all images at a time
        # Peak memory is when all threads are processing their chunks simultaneously
        memory_per_thread = num_images * best_chunk_y * best_chunk_x * 4  # 4 bytes per float32
        peak_memory_estimate = memory_per_thread * n_workers
    
    memory_efficiency_pct = (peak_memory_estimate / (resources['usable_memory_gb'] * 1024**3)) * 100
    
    msg(f'  Optimal chunk size: {best_chunk_y}x{best_chunk_x} pixels')
    msg(f'  Chunk area: {best_chunk_y * best_chunk_x:,} pixels')
    msg(f'  Aspect ratio: {max(best_chunk_y, best_chunk_x)/min(best_chunk_y, best_chunk_x):.2f}')
    msg(f'  Total chunks: {total_chunks}')
    msg(f'  Chunks per worker (avg): {chunks_per_worker:.1f}')
    
    if worker_type == 'processes':
        total_image_memory_gb = (num_images * ny * nx * 4) / (1024**3)
        chunk_memory_gb = (best_chunk_y * best_chunk_x * 4 * n_workers) / (1024**3)
        msg(f'  Shared image stack memory: {total_image_memory_gb:.2f} GB')
        msg(f'  Chunk processing memory: {chunk_memory_gb:.2f} GB') 
    else:
        memory_per_chunk_gb = (num_images * best_chunk_y * best_chunk_x * 4) / (1024**3)
        msg(f'  Memory per chunk: {memory_per_chunk_gb:.2f} GB')
    
    msg(f'  Estimated peak memory: {peak_memory_estimate / (1024**3):.2f} GB')
    msg(f'  Memory efficiency: {memory_efficiency_pct:.1f}% of available')
    msg(f'  Optimization score: {best_score:.3f}')
    
    # Warnings for suboptimal configurations
    if memory_efficiency_pct > 90:
        msg(f'  WARNING: High memory usage ({memory_efficiency_pct:.1f}%) - consider reducing chunk size')
    elif chunks_per_worker < 1:
        msg(f'  WARNING: Fewer chunks than workers - some workers will be idle')
    elif worker_type == 'processes' and chunks_per_worker < 1.5:
        msg(f'  INFO: Consider smaller chunks for better load balancing')
    elif worker_type == 'threads' and chunks_per_worker > 10:
        msg(f'  INFO: Many chunks per thread ({chunks_per_worker:.1f}) - good for I/O parallelism')
    
    return best_chunk_y, best_chunk_x, peak_memory_estimate

# ============================================================================
# PARALLEL PROCESSING STRATEGIES
# ============================================================================

def process_chunk_cpu_worker(args):
    """
    Worker function for CPU-bound chunk processing in multiprocessing environment.
    
    Extracts a chunk from shared memory image stack and computes median across 
    the image axis. Designed for CPU-intensive median computation on data
    that is shared across processes.
    
    Parameters
    ----------
    args : Tuple
        Packed arguments containing:
        - coords : Tuple[int, int]
            Chunk coordinates (i, j) for positioning
        - shm_name : str
            Name of shared memory block containing image stack
        - shape : Tuple[int, int, int]
            Shape of the shared image stack (num_images, height, width)
        - chunk_y : int
            Chunk height in pixels
        - chunk_x : int  
            Chunk width in pixels
        - ny : int
            Total image height
        - nx : int
            Total image width
    
    Returns
    -------
    Tuple
        (median_chunk, y_start, y_end, x_start, x_end) where:
        - median_chunk : np.ndarray
            Computed median as 2D float32 array
        - y_start, y_end, x_start, x_end : int
            Chunk boundaries for result placement
    
    Notes
    -----
    Uses shared memory to avoid copying the entire image stack to each worker process.
    Each process computes median for independent spatial chunks.
    """
    coords, shm_name, shape, chunk_y, chunk_x, ny, nx = args
    i, j = coords
    y_start, y_end = i, min(i + chunk_y, ny)
    x_start, x_end = j, min(j + chunk_x, nx)
    
    try:
        # Attach to existing shared memory
        existing_shm = shared_memory.SharedMemory(name=shm_name)
        # Create numpy array view of shared memory
        image_stack = np.ndarray(shape, dtype=np.float32, buffer=existing_shm.buf)
        
        # Extract chunk from shared memory stack and compute median
        chunk_data = image_stack[:, y_start:y_end, x_start:x_end]
        median_chunk = bn.median(chunk_data, axis=0)
        
        # Don't close the shared memory here - it will be cleaned up by the main process
        return (median_chunk, y_start, y_end, x_start, x_end)
        
    except Exception as e:
        return (None, y_start, y_end, x_start, x_end)

def compute_median_preloaded_shared(shm_name: str, shape: Tuple[int, int, int],
                                   resources: Dict[str, Any],
                                   chunk_size: Optional[int] = None,
                                   n_processes: Optional[int] = None) -> np.ndarray:
    """
    Compute median from images already loaded in shared memory (zero-copy approach).
    
    This function works directly with existing shared memory created by the caller,
    eliminating the need for any data copying. Optimized for maximum memory efficiency.
    
    Parameters
    ----------
    shm_name : str
        Name of existing shared memory block containing image stack
    shape : Tuple[int, int, int]
        Shape of the image stack (num_images, height, width)
    resources : Dict[str, Any]
        System resource information for optimization
    chunk_size : Optional[int]
        Manual chunk size override, auto-detect if None
    n_processes : Optional[int]
        Number of processes for parallel computation, auto-detect if None
        
    Returns
    -------
    np.ndarray
        Median image as 2D float32 array
    
    Notes
    -----
    Memory Strategy:
    - Works directly with existing shared memory (ZERO COPY)
    - No memory duplication at any point
    - Each process works on independent spatial chunks
    - Most memory-efficient approach possible
    
    Performance Characteristics:
    - CPU-bound (median computation is the bottleneck)
    - Optimal when all images fit in memory simultaneously
    - Scales with number of CPU cores available
    """
    num_images, ny, nx = shape
    
    msg(f"Computing median from shared memory stack: {num_images} images ({ny}x{nx}) - ZERO COPY")
    
    # Calculate optimal chunk size for multiprocessing
    if chunk_size is None:
        chunk_y, chunk_x, _ = calculate_optimal_chunk_size(ny, nx, num_images, resources, n_processes, worker_type='processes')
    else:
        chunk_y = chunk_x = chunk_size

    # Use optimal CPU processes (not I/O threads)
    if n_processes is None:
        n_processes = resources['optimal_processes_cpu']

    msg(f"Using chunk size: {chunk_y}x{chunk_x} pixels")
    msg(f"Using {n_processes} CPU processes with existing shared memory")
    
    try:
        # Attach to existing shared memory to verify it exists
        shm = shared_memory.SharedMemory(name=shm_name)
        shared_array = np.ndarray(shape, dtype=np.float32, buffer=shm.buf)
        
        # Result array
        result = np.zeros((ny, nx), dtype=np.float32)
        
        # Generate chunk coordinates
        chunk_coords = [(i, j) for i in range(0, ny, chunk_y) for j in range(0, nx, chunk_x)]
        
        msg(f"Processing {len(chunk_coords)} chunks using {n_processes} CPU processes")
        
        start_time = time.time()
        
        # Prepare arguments for worker function
        worker_args = [(coords, shm_name, shape, chunk_y, chunk_x, ny, nx) for coords in chunk_coords]
        
        # Use ProcessPoolExecutor for CPU-bound median computation
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            futures = [executor.submit(process_chunk_cpu_worker, args) for args in worker_args]
            
            for i, future in enumerate(as_completed(futures)):
                median_chunk, y_start, y_end, x_start, x_end = future.result()
                if median_chunk is not None:
                    result[y_start:y_end, x_start:x_end] = median_chunk
                else:
                    msg(f"Warning: Failed to process chunk {i+1}")
                
                progress = 100 * (i + 1) / len(futures)
                elapsed = time.time() - start_time
                rate = (i + 1) / elapsed
                eta = (len(futures) - i - 1) / rate if rate > 0 else 0
                msg(f"Progress: {i+1}/{len(futures)} chunks ({progress:.1f}%) | "
                    f"Rate: {rate:.1f} chunks/s | ETA: {eta:.0f}s")
        
        elapsed = time.time() - start_time
        msg(f"Median computation completed in {elapsed:.1f}s")
        
        return result
        
    except Exception as e:
        msg(f"Shared memory processing failed: {e}")
        raise  # Let compute_median_chunked handle the fallback

def process_chunk_threaded(args: Tuple) -> Tuple:

    """
    Process a single image chunk using threaded I/O with optimized memory management.
    
    Loads images sequentially, extracts chunk regions, and computes median to minimize
    memory usage for datasets that cannot fit entirely in memory.
    
    Parameters
    ----------
    args : Tuple
        (chunk_id, images, y_start, y_end, x_start, x_end)
        
    Returns
    -------
    Tuple
        (chunk_id, result_chunk, y_start, y_end, x_start, x_end)
        result_chunk is None if processing fails
        
    Notes
    -----
    Uses pre-allocated numpy arrays and loads images one at a time to reduce
    memory footprint. Missing/failed images are filled with zeros.
    """

    chunk_id, images, y_start, y_end, x_start, x_end = args
    
    try:
        chunk_height = y_end - y_start
        chunk_width = x_end - x_start
        num_images = len(images)
        
        # Pre-allocate array for chunk data from all images
        chunk_data = np.zeros((num_images, chunk_height, chunk_width), dtype=np.float32)
        
        # Load chunk data from each image
        for img_idx, img_path in enumerate(images):
            img = get_image(img_path)
            if img is not None:
                # Check dimensions to avoid indexing errors
                if img.shape[0] >= y_end and img.shape[1] >= x_end:
                    chunk_data[img_idx] = img[y_start:y_end, x_start:x_end]
                else:
                    msg(f"Warning: Image {img_path} has incompatible dimensions {img.shape}")
                del img
            # If img is None, chunk_data[img_idx] remains zeros by default initialization
        
        # Compute median across images for this chunk
        result_chunk = bn.median(chunk_data, axis=0)
        
        # Clean up
        del chunk_data
        
        return (chunk_id, result_chunk, y_start, y_end, x_start, x_end)
        
    except Exception as e:
        msg(f"Error processing chunk {chunk_id}: {e}")
        return (chunk_id, None, y_start, y_end, x_start, x_end)


def compute_median_chunked_threaded(images: List[str], 
                                  resources: Dict[str, Any],
                                  chunk_size: Optional[int] = None,
                                  max_images: int = 1024,
                                  n_threads: Optional[int] = None) -> np.ndarray:
    """
    Compute median using threaded I/O for images that cannot fit in memory.
    
    Parameters
    ----------
    images : List[str]
        List of image file paths
    chunk_size : Optional[int]
        Manual chunk size override
    max_images : int
        Maximum number of images to use
    n_threads : Optional[int]
        Number of threads for parallel I/O
        
    Returns
    -------
    np.ndarray
        Median image
    """

    first_img = get_image(images[0])
    if first_img is None:
        raise ValueError("Could not load first image")
    
    ny, nx = first_img.shape
    del first_img
    
    # Set optimal threading parameters
    if n_threads is None:
        n_threads = resources['optimal_threads_io']
    
    if chunk_size is None:
        chunk_y, chunk_x, _ = calculate_optimal_chunk_size(
            ny, nx, len(images), resources, n_threads, worker_type='auto'
        )
    else:
        chunk_y = chunk_x = chunk_size
    
    msg(f"Processing {ny}x{nx} image using threaded I/O with {n_threads} threads")
    msg(f"Using chunk size: {chunk_y}x{chunk_x} pixels") 
    
    # Generate chunk tasks
    chunk_tasks = []
    chunk_id = 0
    for i in range(0, ny, chunk_y):
        for j in range(0, nx, chunk_x):
            y_start, y_end = i, min(i + chunk_y, ny)
            x_start, x_end = j, min(j + chunk_x, nx)
            chunk_tasks.append((chunk_id, images, y_start, y_end, x_start, x_end))
            chunk_id += 1
    
    result = np.zeros((ny, nx), dtype=np.float32)
    
    start_time = time.time()
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = [executor.submit(process_chunk_threaded, task) for task in chunk_tasks]
        
        for i, future in enumerate(as_completed(futures)):
            chunk_id, result_chunk, y_start, y_end, x_start, x_end = future.result()
            if result_chunk is not None:
                result[y_start:y_end, x_start:x_end] = result_chunk
            
            progress = 100 * (i + 1) / len(futures)
            elapsed = time.time() - start_time
            rate = (i + 1) / elapsed
            eta = (len(futures) - i - 1) / rate if rate > 0 else 0
            msg(f"Progress: {i+1}/{len(futures)} chunks ({progress:.1f}%) | "
                f"Rate: {rate:.1f} chunks/s | ETA: {eta:.0f}s")
    
    elapsed = time.time() - start_time
    msg(f"Threaded median computation completed in {elapsed:.1f}s")
    
    return result


def compute_median_chunked(images: List[str], 
                         resources: Dict[str, Any],
                         chunk_size: Optional[int] = None,
                         max_images: int = 1024) -> np.ndarray:
    """
    Intelligently compute median using optimal strategy based on available resources.
    
    Parameters
    ----------
    images : List[str]
        List of image file paths
    chunk_size : Optional[int]
        Manual chunk size override
    max_images : int
        Maximum number of images to use
        
    Returns
    -------
    np.ndarray
        Median image
    """    
    # Limit images if necessary
    if len(images) > max_images:
        msg(f"Limiting to {max_images} images from {len(images)} available")
        indices = np.linspace(0, len(images) - 1, max_images, dtype=int)
        images = [images[i] for i in indices]
    
    # Determine strategy based on available memory
    can_preload = can_load_all_images_in_memory(images, resources)
    
    if can_preload:
        msg("Strategy: Pre-loading all images for shared memory processing")
        
        # Get dimensions from first image
        first_img = get_image(images[0])
        if first_img is None:
            raise ValueError("Could not load first image to determine dimensions")
        
        ny, nx = first_img.shape
        num_images = len(images)
        
        # Try shared memory approach first
        try:
            # Create unique shared memory name to avoid conflicts
            shm_name = f"homogenize_beams_{uuid.uuid4().hex[:8]}_{os.getpid()}"
            
            # Create shared memory for the image stack
            shape = (num_images, ny, nx)
            nbytes = num_images * ny * nx * 4  # 4 bytes per float32
            msg(f"Creating shared memory stack: {num_images} images of {ny}x{nx} pixels ({nbytes / (1024**3):.2f} GB)")
            
            # Create shared memory block with unique name
            shm = shared_memory.SharedMemory(create=True, size=nbytes, name=shm_name)
            
            # Create numpy array view of shared memory
            shared_image_stack = np.ndarray(shape, dtype=np.float32, buffer=shm.buf)
            
            # Load first image directly into shared memory
            shared_image_stack[0] = first_img.astype(np.float32)
            loaded_count = 1
            
            # Load remaining images directly into shared memory
            for i in range(1, len(images)):
                img = get_image(images[i])
                if img is not None:
                    if img.shape != (ny, nx):
                        msg(f"Warning: Image {images[i]} has different dimensions {img.shape}, expected {(ny, nx)}")
                        continue
                    shared_image_stack[loaded_count] = img.astype(np.float32)
                    loaded_count += 1
                
                if (i + 1) % 10 == 0:
                    msg(f"Loaded {loaded_count}/{len(images)} images into shared memory")
            
            if loaded_count == 0:
                raise ValueError("Could not load any images")
            
            # Trim the shared memory view if some images failed to load
            if loaded_count < num_images:
                msg(f"Trimming shared memory view from {num_images} to {loaded_count} images due to load failures")
                shared_image_stack = shared_image_stack[:loaded_count]
            
            # Compute median from shared memory
            try:
                result = compute_median_preloaded_shared(shm_name, shared_image_stack.shape, resources, chunk_size)
                return result
            except Exception as median_error:
                msg(f"Shared memory median computation failed: {median_error}")
                raise  # Re-raise to trigger the outer fallback
            finally:
                # Clean up shared memory after processing
                try:
                    shm.close()
                    shm.unlink()
                    msg("Shared memory cleaned up successfully")
                except Exception as cleanup_error:
                    msg(f"Warning: Shared memory cleanup failed: {cleanup_error}")
            
        except Exception as e:
            msg(f"Shared memory approach failed: {e}")
            msg("Falling back to simple array approach...")
            
            # Clean up shared memory if it was created
            try:
                if 'shm' in locals():
                    shm.close()
                    shm.unlink()
            except:
                pass
            
            # Simple fallback: Load all images into regular array and compute median
            msg(f"Loading {num_images} images into regular array: {ny}x{nx} pixels")
            image_stack = np.zeros((num_images, ny, nx), dtype=np.float32)
            
            # Load first image into the stack
            image_stack[0] = first_img.astype(np.float32)
            loaded_count = 1
            
            # Load remaining images
            for i in range(1, len(images)):
                img = get_image(images[i])
                if img is not None:
                    if img.shape != (ny, nx):
                        msg(f"Warning: Image {images[i]} has different dimensions {img.shape}, expected {(ny, nx)}")
                        continue
                    image_stack[loaded_count] = img.astype(np.float32)
                    loaded_count += 1
                
                if (i + 1) % 10 == 0:
                    msg(f"Loaded {loaded_count}/{i+1} images into array")
            
            if loaded_count == 0:
                raise ValueError("Could not load any images")
            
            # Trim the array if some images failed to load
            if loaded_count < num_images:
                msg(f"Trimming array from {num_images} to {loaded_count} images due to load failures")
                image_stack = image_stack[:loaded_count]
            
            # Simple median computation
            msg("Computing median using simple array approach...")
            start_time = time.time()
            result = bn.median(image_stack, axis=0)
            elapsed = time.time() - start_time
            msg(f"Simple median computation completed in {elapsed:.1f}s")
            
            return result
    else:
        msg("Strategy: Using threaded I/O for memory-efficient processing")
        return compute_median_chunked_threaded(images, resources, chunk_size, max_images)

# ============================================================================
# IMAGE PROCESSING AND HOMOGENIZATION
# ============================================================================

def get_identifiers(prefix: str) -> Tuple[List[str], List[str]]:
    """
    Extract image identifiers with enhanced logic for time series handling.
    
    Parameters
    ----------
    prefix : str
        Input prefix for image search
        
    Returns
    -------
    Tuple[List[str], List[str]]
        (beam_identifiers, image_identifiers)
    """
    msg('Analyzing image naming patterns and grouping strategies')
    
    image_arr = sorted(glob.glob(f'{prefix}-t*'))
    
    if image_arr and not cfg.WSC_HOMOGENIZETIME:
        # Time-resolved processing
        suffix = np.unique([im.split(f'{prefix}-t')[-1].split('-')[0] for im in image_arr])
        beam_ids = [f'{prefix}-t{s}' for s in suffix]
        image_ids = [f'{prefix}-t{s}' for s in suffix]
        msg(f'Time-resolved mode: Found {len(suffix)} time intervals')
        return beam_ids, image_ids
    elif image_arr and cfg.WSC_HOMOGENIZETIME:
        # Time-averaged processing
        suffix = np.unique([im.split(f'{prefix}-t')[-1].split('-')[0] for im in image_arr])
        beam_ids = [prefix]
        image_ids = [f'{prefix}-t{s}' for s in suffix]
        msg(f'Time-averaged mode: {len(suffix)} intervals, single beam solution')
        return beam_ids, image_ids
    else:
        # Standard processing
        msg('Standard mode: Single identifier processing')
        return [prefix], [prefix]

def process_psf(psf: str, beam: Tuple[float, float, float], identifier: str) -> Dict[str, Any]:
    """
    Process a single PSF file: generate the homogenization kernel via pypher and
    apply it to all associated Stokes image/residual pairs.

    This is the worker function for parallel execution in homogenize_images().
    It must remain a top-level function so that ProcessPoolExecutor can pickle it.

    Parameters
    ----------
    psf : str
        Path to the PSF FITS file for this channel.
    beam : Tuple[float, float, float]
        Target beam parameters (semi-major, semi-minor, position_angle) as
        returned by get_homogenized_beam().  Position angle is in radians.
    identifier : str
        Base identifier for the image set, used only for log context.

    Returns
    -------
    Dict[str, Any]
        Status dictionary with keys:
          'psf'       - the input psf path
          'success'   - bool, True if the PSF was processed without fatal error
          'processed' - int, 1 if completed (including already-skipped), else 0
          'failed'    - int, 1 if a fatal error occurred, else 0
          'message'   - human-readable summary string
    """
    a, b, p = beam[0], beam[1], beam[2]

    # Derive sibling filenames
    psf_header = fits.getheader(psf)
    psf_zoom = psf.replace('psf.fits', 'psf.zoom.fits')
    psf_new  = psf.replace('psf.fits', 'psf.new.fits')
    kernel   = psf.replace('psf.fits', 'psf.kernel.fits')

    sky_to_pix  = (abs(psf_header.get('CDELT1'))) ** (-1)
    fwhm_to_sig = (8 * np.log(2)) ** (-0.5)
    psf_size    = 101

    # Safety check: skip flagged channels (WSCLEAN sets BMAJ=0 for flagged channels)
    if psf_header['BMAJ'] <= 1e-14:
        return {
            'psf': psf, 'success': False, 'processed': 0, 'failed': 1,
            'message': f'Skipping {psf} - invalid beam (BMAJ <= 1e-14)'
        }

    # If kernel already exists, skip generation and go straight to convolution
    if os.path.exists(kernel):
        msg(f'Kernel exists, skipping pypher for {os.path.basename(psf)}')
    else:
        # --- Build original (zoomed) PSF ---
        psf_zoom_header = psf_header.copy()
        a_pix = psf_zoom_header['BMAJ'] * sky_to_pix * fwhm_to_sig
        b_pix = psf_zoom_header['BMIN'] * sky_to_pix * fwhm_to_sig
        p_pix = np.radians(psf_zoom_header['BPA']) + 0.5 * np.pi

        psf_zoom_header['NAXIS1'] = psf_size
        psf_zoom_header['NAXIS2'] = psf_size
        psf_zoom_image = Gaussian2DKernel(
            x_stddev=a_pix, y_stddev=b_pix, theta=p_pix,
            x_size=psf_size, y_size=psf_size, mode='center'
        ).array
        psf_zoom_image /= np.amax(psf_zoom_image)
        fits.PrimaryHDU(data=psf_zoom_image, header=psf_zoom_header).writeto(psf_zoom, overwrite=True)

        # --- Build target (new) PSF ---
        a_pix = a * sky_to_pix
        b_pix = b * sky_to_pix

        psf_new_header = psf_header.copy()
        psf_new_header['NAXIS1'] = psf_size
        psf_new_header['NAXIS2'] = psf_size
        psf_new_header['BMAJ']   = a / fwhm_to_sig
        psf_new_header['BMIN']   = b / fwhm_to_sig
        psf_new_header['BPA']    = np.degrees(p - 0.5 * np.pi)
        psf_new_image = Gaussian2DKernel(
            x_stddev=a_pix, y_stddev=b_pix, theta=p,
            x_size=psf_size, y_size=psf_size, mode='center'
        ).array
        psf_new_image /= np.amax(psf_new_image)
        fits.PrimaryHDU(data=psf_new_image, header=psf_new_header).writeto(psf_new, overwrite=True)

        # --- Run pypher to generate the homogenization kernel ---
        try:
            msg(f'Generating convolution kernel using pypher for {os.path.basename(psf)}')
            subprocess.run(
                [f'pypher {psf_zoom} {psf_new} {kernel}'],
                shell=True, check=True
            )
            msg(f'Successfully created kernel: {os.path.basename(kernel)}')
        except subprocess.CalledProcessError as e:
            return {
                'psf': psf, 'success': False, 'processed': 0, 'failed': 1,
                'message': f'pypher failed for {os.path.basename(psf)}: {e}'
            }

    # Re-derive psf_new_image for the convolution step (needed whether kernel
    # was just built or already existed on disk)
    a_pix = a * sky_to_pix
    b_pix = b * sky_to_pix
    psf_new_image = Gaussian2DKernel(
        x_stddev=a_pix, y_stddev=b_pix, theta=p,
        x_size=psf_size, y_size=psf_size, mode='center'
    ).array
    psf_new_image /= np.amax(psf_new_image)

    # --- Apply kernel to all associated Stokes images ---
    images = glob.glob(psf.split('-psf.fits')[0] + '*-image.fits')
    images = [im for im in images if not ('-Plin-' in im or '-Ptot-' in im)]

    for im in images:
        try:
            image_name    = im.replace('image.fits', 'image.homogenized.fits')
            residual      = im.replace('image.fits', 'residual.fits')
            residual_name = im.replace('image.fits', 'residual.homogenized.fits')
            model         = im.replace('image.fits', 'model.fits')

            model_image    = get_image(model)
            residual_image = get_image(residual)
            kernel_image   = get_image(kernel)

            if model_image is None or residual_image is None or kernel_image is None:
                msg(f'Warning: Could not load required images for {im}, skipping')
                continue

            # Convolve model with new PSF, residual with homogenization kernel
            image_new  = scipy.signal.fftconvolve(model_image,    psf_new_image, mode='same')
            image_rms  = scipy.signal.fftconvolve(residual_image, kernel_image,  mode='same')
            image_new += image_rms

            header = fits.getheader(im)
            header['BMAJ']    = a / fwhm_to_sig
            header['BMIN']    = b / fwhm_to_sig
            header['BPA']     = np.degrees(p - 0.5 * np.pi)
            header['HISTORY'] = 'Beam homogenized using pypher kernel'

            create_fits(np.asarray(image_new, dtype=np.float32), header, image_name)
            create_fits(np.asarray(image_rms, dtype=np.float32), header, residual_name)

            msg(f'Successfully processed: {os.path.basename(image_name)}')

        except Exception as e:
            msg(f'Error processing image {im}: {e}')
            continue

    return {
        'psf': psf, 'success': True, 'processed': 1, 'failed': 0,
        'message': f'Completed: {os.path.basename(psf)}'
    }


def homogenize_images(identifier: str, beam: Tuple[float, float, float], 
                     resources: Dict[str, Any], sigma_threshold: float = 3.0) -> None:
    """
    Homogenize the synthesized beam across a set of radio interferometry images.
    
    This function takes a collection of images with varying beam sizes and convolves
    them to a common, larger beam size to enable direct comparison and analysis.
    The process involves:
    1. Filtering out PSF files with abnormal beam sizes (outliers)
    2. Creating convolution kernels to transform each original beam to the target beam
    3. Applying these kernels to both model and residual components of each image
    4. Generating homogenized versions of all images with consistent resolution
    5. Creating pseudo-MFS (Multi-Frequency Synthesis) images if not already present
    
    Parameters
    ----------
    identifier : str
        Base identifier/prefix for the image set (e.g., 'source_name-t001').
        Used to locate all related PSF and image files with this prefix.
    beam : tuple or array-like
        Target beam parameters as (major_axis, minor_axis, position_angle).
        - major_axis : float
            Semi-major axis of target beam in degrees
        - minor_axis : float  
            Semi-minor axis of target beam in degrees
        - position_angle : float
            Position angle of target beam in radians (measured from +x axis)
    resources : Dict[str, Any]
        System resource information
    sigma_threshold : float
        Sigma threshold for PSF outlier rejection (default: 3.0)
    
    Returns
    -------
    None
        Function operates by side-effect, creating new homogenized FITS files:
        - *-image.homogenized.fits : Homogenized total intensity images
        - *-residual.homogenized.fits : Homogenized residual images  
        - *-MFS*-image.homogenized.fits : Pseudo-MFS images (if MFS PSF absent)
    
    Notes
    -----
    - Requires pypher package for PSF convolution kernel generation
    - Only processes channels with valid beam information (BMAJ > 1e-14)
    - Automatically detects and processes all Stokes parameters present
    - Memory-efficient chunked processing for large images via compute_median_chunked()
    - Original images and headers are preserved; only homogenized versions created
    - Filters out PSF files with abnormal beam sizes based on frequency normalization
    
    File Dependencies
    -----------------
    Input files (must exist):
        - {identifier}*-psf.fits : Point Spread Function files
        - {identifier}*-image.fits : Total intensity image files
        - {identifier}*-model.fits : Model component files  
        - {identifier}*-residual.fits : Residual image files
    
    Output files (created):
        - {identifier}*-image.homogenized.fits : Beam-homogenized images
        - {identifier}*-residual.homogenized.fits : Beam-homogenized residuals
        - {identifier}-MFS*-image.homogenized.fits : Pseudo-MFS images
    """

    # Split out the beam components
    a, b, p = beam[0], beam[1], beam[2]
    
    msg(f"Homogenizing images for identifier: {identifier}")
    msg(f"Target beam: semi-major={a:.2e}, semi-minor={b:.2e}, angle={np.degrees(p):.1f}°")
            
    # Get list of PSF files
    psfs = sorted(glob.glob(f'{identifier}-*-psf.fits'))
    if not psfs:
        msg(f"Warning: No PSF files found for {identifier}")
        return
    
    # Filter out outlier PSFs based on frequency-normalized beam sizes
    good_psfs, rejected_psfs = filter_beam_outliers(psfs, sigma_threshold)
    
    if not good_psfs:
        msg(f"Error: All PSF files rejected as outliers for {identifier}")
        return
    
    msg(f"Processing {len(good_psfs)} good PSF files (rejected {len(rejected_psfs)} outliers)")
    
    processed_count = 0
    failed_count = 0

    # Determine worker count from available memory and CPU count.
    # Each worker runs two fftconvolve calls; budget ~14 array-equivalents
    # (inputs + FFT working buffers + outputs) plus 100 MB process overhead.
    available_gb = get_available_memory_gb()
    sample_image = good_psfs[0].replace('-psf.fits', '-I-image.fits')
    if not os.path.exists(sample_image):
        # Fall back to any image associated with the first PSF
        candidates = glob.glob(good_psfs[0].split('-psf.fits')[0] + '*-image.fits')
        sample_image = candidates[0] if candidates else None

    if sample_image:
        with fits.open(sample_image) as hdul:
            array_bytes = hdul[0].data.nbytes
        per_worker_gb = (14 * array_bytes + 100 * 1024**2) / 1024**3
    else:
        per_worker_gb = 2.0
        msg('Warning: Could not find sample image for memory estimation, assuming 2.00 GB/worker')

    mem_cap  = max(1, int(available_gb // per_worker_gb))
    # --- SLURM-aware CPU cap ---
    # SLURM_CPUS_ON_NODE is intentionally excluded: it reports the total node
    # CPU count, not the allocation, and would defeat the purpose of this check.
    import re as _re
    _cpus_per_task = os.environ.get('SLURM_CPUS_PER_TASK')
    _ntasks        = os.environ.get('SLURM_NTASKS')
    _job_cpus      = os.environ.get('SLURM_JOB_CPUS_PER_NODE', '')
    _job_cpus_m    = _re.match(r'(\d+)', _job_cpus)

    if _cpus_per_task and _ntasks:
        cpu_cap = int(_cpus_per_task) * int(_ntasks)
        _cpu_src = f'SLURM_CPUS_PER_TASK={_cpus_per_task} x SLURM_NTASKS={_ntasks}'
    elif _cpus_per_task:
        cpu_cap = int(_cpus_per_task)
        _cpu_src = f'SLURM_CPUS_PER_TASK={_cpus_per_task}'
    elif _job_cpus_m:
        cpu_cap = int(_job_cpus_m.group(1))
        _cpu_src = f'SLURM_JOB_CPUS_PER_NODE={_job_cpus} -> {cpu_cap}'
    elif _ntasks:
        cpu_cap = int(_ntasks)
        _cpu_src = f'SLURM_NTASKS={_ntasks}'
    else:
        cpu_cap = os.cpu_count() or 1
        _cpu_src = f'os.cpu_count() (no SLURM allocation detected)'

    n_workers = min(mem_cap, cpu_cap, len(good_psfs))

    msg(f'Memory available      : {available_gb:.2f} GB  [{"SLURM" if any(os.environ.get(v) for v in ("SLURM_MEM_PER_NODE","SLURM_MEM_PER_CPU")) else "psutil/proc"}]')
    msg(f'Est. memory/worker    : {per_worker_gb:.2f} GB')
    msg(f'Workers (mem cap)     : {mem_cap}')
    msg(f'Workers (CPU cap)     : {cpu_cap}  [{_cpu_src}]')
    msg(f'Workers (selected)    : {n_workers}  [min(mem={mem_cap}, cpu={cpu_cap}, images={len(good_psfs)})]')

    # Dispatch all PSFs in parallel; each worker is fully independent
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(process_psf, psf, beam, identifier): psf
            for psf in good_psfs
        }
        for future in as_completed(futures):
            try:
                result = future.result()
                msg(result['message'])
                processed_count += result['processed']
                failed_count    += result['failed']
            except Exception as e:
                msg(f'ERROR processing {futures[future]}: {e}')
                failed_count += 1

    msg(f"Homogenization summary for {identifier}:")
    msg(f"  Total PSF files found: {len(psfs)}")
    msg(f"  Good PSF files (passed filter): {len(good_psfs)}")
    msg(f"  Rejected PSF files (outliers): {len(rejected_psfs)}")

    # Check if MFS image exists, if not make a pseudo MFS image
    if glob.glob(f'{identifier}-MFS-psf.fits') == []:
          
        # Figure out what the Stokes parameters are included in the images
        stokes = []
        for stoke in ['-I', '-Q', '-U', '-V']:
            if glob.glob(f'{identifier}*{stoke}-*') != []:
                stokes.append(stoke)
        if stokes == []:        
            stokes = ['']

        # Make a MFS image for each stokes parameter
        for stoke in stokes:
            msg(f'Making a (pseudo-) Stokes{stoke} MFS image; take these images with a grain of salt')
            suffix = f'{stoke}-image.homogenized.fits'
            # suffix = f'{stoke}-image.fits' # REMOVE
            images = [im.replace('-psf.fits', suffix)  for im in good_psfs if os.path.exists(im.replace('-psf.fits', suffix))]

            if not images:
                msg(f'No homogenized images found for Stokes {stoke}')
                continue

            freq = []
            for k, im in enumerate(images):
                try:
                    freq.append(fits.getheader(im)['CRVAL3'])
                    if k == 0:
                        header = fits.getheader(im)
                except Exception as e:
                    msg(f'Warning: Could not read frequency from {im}: {e}')

            # Adopt median values for each pixel and output MFS image
            msg(f'Computing median for {len(images)} images using chunked processing')
            data = compute_median_chunked(images, resources=resources)
 
            if freq:
                header['CRVAL3'] = np.nanmean(freq)
            header['HISTORY'] = f'Pseudo-MFS created from {len(images)} homogenized images'
            mfs_name = f'{identifier}-MFS{stoke}-image.homogenized.fits'
            create_fits(data, header, mfs_name)
            msg(f'Created pseudo-MFS: {os.path.basename(mfs_name)}')

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main() -> None:
    """Enhanced main function with comprehensive resource management."""
    
    # Display system information
    system_info = get_system_resources()
    
    # Parse command line arguments
    if len(sys.argv) not in [2, 3, 4, 5]:
        msg('ERROR: Usage: python homogenize_beams.py <prefix> [max_images] [chunk_size] [sigma_threshold]')
        msg('  prefix: Image prefix (required)')
        msg('  max_images: Maximum number of images to use (optional, default: 1024)')
        msg('  chunk_size: Manual chunk size override (optional, auto-detect if not provided)')
        msg('  sigma_threshold: Outlier rejection threshold in sigma (optional, default: 3.0)')
        sys.exit(1)
    
    prefix = sys.argv[1]
    max_images = 1024
    chunk_size = None
    sigma_threshold = 3.0  
    
    if len(sys.argv) >= 3:
        try:
            max_images = int(sys.argv[2])
            msg(f"Using maximum {max_images} images (user specified)")
        except ValueError:
            msg("Warning: Invalid max_images, using default 1024")
            max_images = 1024
    
    if len(sys.argv) >= 4:
        try:
            chunk_size = int(sys.argv[3])
            msg(f"Using chunk size {chunk_size}x{chunk_size} (user specified)")
        except ValueError:
            msg("Warning: Invalid chunk_size, using auto-detection")
            chunk_size = None
    
    if len(sys.argv) >= 5:
        try:
            sigma_threshold = float(sys.argv[4])
            msg(f"Using sigma threshold {sigma_threshold} for outlier rejection")
        except ValueError:
            msg("Warning: Invalid sigma_threshold, using default 3.0")
            sigma_threshold = 3.0
            
    try:
        # Get image identifiers
        beam_identifiers, image_identifiers = get_identifiers(prefix)
        msg(f"Found {len(beam_identifiers)} beam groups and {len(image_identifiers)} image groups")
        
        # Compute homogenized beams
        beams = []
        for identifier in beam_identifiers:
            msg(f'Computing minimum enclosing ellipse for: {identifier}')
            try:
                beam = get_homogenized_beam(identifier)
                beams.append(beam)
                msg(f'Successfully computed beam parameters for {identifier}')
            except Exception as e:
                msg(f'Error computing beam for {identifier}: {e}')
                raise
        
        msg(f'Computed {len(beams)} homogenized beam(s)')
        
        # Homogenize images
        for k, identifier in enumerate(image_identifiers):
            msg(f'Homogenizing images for: {identifier} ({k+1}/{len(image_identifiers)})')
            try:
                beam = beams[0] if len(beams) == 1 else beams[k]
                homogenize_images(identifier, beam, system_info)
                msg(f'Successfully homogenized images for {identifier}')
            except Exception as e:
                msg(f'Error homogenizing images for {identifier}: {e}')
                continue
        
        msg('Beam homogenization completed successfully')
        
    except Exception as e:
        msg(f'Fatal error in main execution: {e}')
        sys.exit(1)

if __name__ == '__main__':
    # Guard required for ProcessPoolExecutor on macOS/Windows (spawn start method).
    # Safe no-op on Linux (fork), but good practice regardless.
    main()
