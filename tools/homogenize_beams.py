#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk
# fraser.cowie@physics.ox.ac.uk

import os
import sys
import glob
import subprocess
import gc
import time
import scipy
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from scipy.spatial import ConvexHull
import psutil

import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

# Helper functions
def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt, flush=True)

# Fits manipulation functions
def flush_fits(newimage,newheader, fitsfile):
    f = fits.open(fitsfile,mode='update')
    input_hdu = f[0]
    input_hdu.header = newheader
    if len(input_hdu.data.shape) == 2:
            input_hdu.data[:,:] = newimage
    elif len(input_hdu.data.shape) == 3:
            input_hdu.data[0,:,:] = newimage
    else:
            input_hdu.data[0,0,:,:] = newimage
    f.flush()

def create_fits(newimage, newheader, fitsfile):
    """
    Create a new FITS file with the provided image data and header.
    
    Parameters
    ----------
    newimage : numpy.ndarray
        Image data to write to the FITS file
    newheader : astropy.io.fits.Header
        FITS header to use for the new file
    fitsfile : str
        Path for the new FITS file to create
        
    Notes
    -----
    This function creates a brand new FITS file, unlike flush_fits which 
    updates an existing file. The output is always 4D to maintain proper
    radio astronomy FITS coordinate structure.
    """
    # Always create 4D data structure for radio astronomy FITS
    if len(newimage.shape) == 2:
        # 2D image - add Stokes and frequency axes
        data_to_write = newimage.reshape(1, 1, newimage.shape[0], newimage.shape[1])
    elif len(newimage.shape) == 3:
        # 3D image - add Stokes axis
        data_to_write = newimage.reshape(1, newimage.shape[0], newimage.shape[1], newimage.shape[2])
    else:
        # 4D image - use as is
        data_to_write = newimage
    
    # Create the HDU with data and header
    hdu = fits.PrimaryHDU(data=data_to_write, header=newheader)
    
    # Write to file
    hdu.writeto(fitsfile, overwrite=True)
    
    msg(f'Created FITS file: {fitsfile}')

def get_image(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
            image = np.array(input_hdu.data[:,:])
    elif len(input_hdu.data.shape) == 3:
            image = np.array(input_hdu.data[0,:,:])
    else:
            image = np.array(input_hdu.data[0,0,:,:])
    return image

# Ellipse functions are based on the Welzl recursion algorithmn
# Taken from https://github.com/dorshaviv/lowner-john-ellipse 

def ellipse_from_boundary5(S):
    """
    Compute the unique ellipse that passes through 5 boundary points.

    Arguments:
        S: an array of shape (5,2) containing points in R2 as row
            vectors, which are on the boundary of the desired ellipse.

    Returns:
        an ellipse given by a tuple (c, a, b, t), where c = (x, y) is the
            center, a and b are the major and minor radii, and t is the
            rotation angle.
    """
    assert S.shape == (5, 2)

    # find parameters of ellipse given in the form
    # s0 * x ** 2 + s1 * y ** 2 + 2 * s2 * x * y + s3 * x + s4 * y + 1 = 0.

    # build linear system of equations:
    x = S[:, 0]
    y = S[:, 1]
    A = np.column_stack((x**2, y**2, 2 * x * y, x, y))

    # if A is close to singular, then at least 3 points are colinear, in which
    # case an ellipse is not unique, then we give up on this ellipse
    if is_singular(A):
        return None

    # solve system of equations
    sol = np.linalg.solve(A, -np.ones(S.shape[0]))

    # find ellipse center
    c = np.linalg.solve(-2 * np.array([[sol[0], sol[2]], [sol[2], sol[1]]]), sol[3:5])

    # solve for the matrix F (ellipse representation in center form)
    A = np.vstack(
        [
            np.hstack([np.eye(3), -np.array([[sol[0], sol[2], sol[1]]]).T]),
            np.array([c[0] ** 2, 2 * c[0] * c[1], c[1] ** 2, -1]),
        ]
    )
    s = np.linalg.solve(A, np.array([0, 0, 0, 1]))
    F = np.array([[s[0], s[1]], [s[1], s[2]]])

    return center_form_to_geometric(F, c)


def ellipse_from_boundary4(S):
    """
    Compute the smallest ellipse that passes through 4 boundary points,
    based on the algorithm by:
    B. W. Silverman and D. M. Titterington, "Minimum covering ellipses,"
    SIAM Journal on Scientific and Statistical Computing 1, no. 4 (1980):
    401-409.

    Arguments:
        S: an array of shape (4,2) containing points in R2 as row
            vectors, which are on the boundary of the desired ellipse.

    Returns:
        an ellipse given by a tuple (c, a, b, t), where c = (x, y) is the
            center, a and b are the major and minor radii, and t is the
            rotation angle. This ellipse is the ellipse with the smallest
            area that passes through the 4 points.
    """
    assert S.shape == (4, 2)

    # sort coordinates in clockwise order
    Sc = S - np.mean(S, axis=0)
    angles = np.arctan2(Sc[:, 1], Sc[:, 0])
    S = S[np.argsort(-angles), :]

    # find intersection point of diagonals
    A = np.column_stack([S[2, :] - S[0, :], S[1, :] - S[3, :]])

    # if A is singular, give up on this ellipse
    if is_singular(A):
        return None

    b = S[1, :] - S[0, :]
    s = np.linalg.solve(A, b)
    diag_intersect = S[0, :] + s[0] * (S[2, :] - S[0, :])

    # shift to origin
    S = S - diag_intersect

    # rotate so one diagonal is parallel to x-axis
    AC = S[2, :] - S[0, :]
    theta = np.arctan2(AC[1], AC[0])
    rot_mat = np.array(
        [[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]]
    )
    S = rot_mat.dot(S.T).T

    # shear parallel to x-axis to make diagonals perpendicular
    m = (S[1, 0] - S[3, 0]) / (S[3, 1] - S[1, 1])
    shear_mat = np.array([[1, m], [0, 1]], dtype=float)
    S = shear_mat.dot(S.T).T

    # make the quadrilateral cyclic (i.e. all vertices lie on a circle)
    b = np.linalg.norm(S, axis=1)
    d = b[1] * b[3] / (b[2] * b[0])
    stretch_mat = np.diag(np.array([d**0.25, d**-0.25], dtype=float))
    S = stretch_mat.dot(S.T).T

    # compute optimal swing angle by solving cubic equation
    a = np.linalg.norm(S, axis=1)
    coeff = np.zeros(4)
    coeff[0] = -4 * a[1] ** 2 * a[2] * a[0]
    coeff[1] = -4 * a[1] * (a[2] - a[0]) * (a[1] ** 2 - a[2] * a[0])
    coeff[2] = (
        3 * a[1] ** 2 * (a[1] ** 2 + a[2] ** 2)
        - 8 * a[1] ** 2 * a[2] * a[0]
        + 3 * (a[1] ** 2 + a[2] ** 2) * a[0] ** 2
    )
    coeff[3] = coeff[1] / 2.0
    rts = np.roots(coeff)
    # take the unique root in the interval (-1, 1)
    rts = rts[(-1 < rts) & (rts < 1)]
    theta = np.arcsin(np.real(rts[0]))

    # apply transformation D_theta
    D_mat = np.array(
        [
            [np.cos(theta) ** -0.5, np.sin(theta) * np.cos(theta) ** -0.5],
            [0, np.cos(theta) ** 0.5],
        ]
    )
    S = D_mat.dot(S.T).T

    # find enclosing circle
    boundary = S[:-1, :]  # only 3 points are needed
    A = np.vstack([-2 * boundary.T, np.ones(boundary.shape[0])]).T
    if is_singular(A):
        return None
    b = -np.sum(boundary**2, axis=1)
    s = np.linalg.solve(A, b)

    # circle parameters (center and radius)
    circle_c = s[:2]
    circle_r = np.sqrt(np.sum(circle_c**2) - s[2])

    # total affine transform that was applied
    T_mat = D_mat.dot(stretch_mat).dot(shear_mat).dot(rot_mat)

    # find original ellipse parameters (in center form)
    ellipse_c = np.linalg.solve(T_mat, circle_c) + diag_intersect
    ellipse_F = T_mat.T.dot(T_mat) / circle_r**2

    return center_form_to_geometric(ellipse_F, ellipse_c)


def ellipse_from_boundary3(S):
    """
    Compute the smallest ellipse that passes through 3 boundary points.

    Arguments:
        S: an array of shape (3,2) containing points in R2 as row
            vectors, which are on the boundary of the desired ellipse.

    Returns:
        an ellipse given by a tuple (c, a, b, t), where c = (x, y) is the
            center, a and b are the major and minor radii, and t is the
            rotation angle. This ellipse is the ellipse with the smallest
            area that passes through the 3 points.
    """
    assert S.shape == (3, 2)

    # centroid
    c = np.mean(S, axis=0)

    # shift points
    Sc = S - c

    # if Sc is close to singular, then the 3 points are colinear, in which
    # case an ellipse is not unique, then we give up on this ellipse
    if is_singular(Sc):
        return None

    # ellipse matrix (center form)
    F = 1.5 * np.linalg.inv(Sc.T.dot(Sc))

    return center_form_to_geometric(F, c)

def center_form_to_geometric(F, c):
    """
    Convert ellipse represented in centre form:
        (x - c)^T * F * (x - c) = 1
    to geometrical representation, i.e. center, major-axis, minor-axis, and
    rotation angle.

    Arguments:
        F: array of shape (2,2), the matrix in the ellipse representation.
        c: array of length 2, the ellipse center.

    Returns:
        a tuple (c, a, b, t), where c = (x, y) is the center, a and
            b are the major and minor radii, and t is the rotation angle.
    """

    # extract a, b, and t from F by finding eigenvalues and eigenvectors
    w, V = np.linalg.eigh(F)

    # the eigenvalues are 1/a**2 and 1/b**2
    # The eigenvectors form the rotation matrix with angle t

    # if one the eigenvalues is not positive, the ellipse is degenerate
    if w[0] <= 0 or w[1] <= 0:
        return None

    # we assume without loss of generality 0 < t < pi.
    # V[1, 0] = sin(t), therefore it must be non-negative:
    if V[1, 0] < 0:
        V[:, 0] = -V[:, 0]

    # find t
    t = np.arccos(V[0, 0])  # V[0, 0] = cos(t)

    return c, 1 / np.sqrt(w[0]), 1 / np.sqrt(w[1]), t

def is_singular(A):
    """Checks if matrix is close to singular.

    Args:
        A: matrix

    Returns:
        bool: True if A is close to singular.
    """

    return np.linalg.cond(A) >= 1 / np.finfo(float).eps


def sample_ellipse(ellipse, num_pts = 100, endpoint=True):
    """
    Uniformly sample points on an ellipse.

    Arguments:
        ellipse: a tuple (a, b, p), assuming center is (0,0); a and
            b are the major and minor radii, and p is the rotation angle from +x direction.
        num_pts: number of points to sample.
        endpoint: boolean. If True, repeat first point at the end (used for
            plotting).

    Returns:
        x: an array of shape (num_pts, 2) containing the sampled points as row
            vectors.
    """

    # extract ellipse parameters
    a, b, p = ellipse

    # rotation matrix
    rot_mat = np.array([[np.cos(p), -np.sin(p)], [np.sin(p), np.cos(p)]])

    # array of angles uniformly chosen between 0 and 2 * pi
    theta = np.linspace(0, 2 * np.pi, num_pts, endpoint=endpoint)

    # points on an ellipse with axis a, b before rotation and shift
    z = np.column_stack((a * np.cos(theta), b * np.sin(theta)))

    # rotate points by angle t and shift to center c
    x = rot_mat.dot(z.T).T

    return x


def plot_ellipse(ellipse, num_pts=100, str="-"):
    """
    Plot ellipse.

    Arguments:
        ellipse: a tuple (a, b, p), assuming center is (0,0); a and
            b are the major and minor radii, and t is the rotation angle.
        num_pts: number of points to sample the ellipse and plot.
        str: plot string to be passed to plot function.
    """

    # if ellipse is empty, do nothing
    if ellipse is None:
        return

    # sample points on ellipse
    x = sample_ellipse(ellipse, num_pts)

    # plot ellipse
    plt.plot(x[:, 0], x[:, 1], str, label='Minimum enclosing ellipse', zorder = 100000)


def is_in_ellipse(point, ellipse):
    """
    Check if a point is contained in an ellipse.

    Arguments:
        point: array of length 2 representing a point in R2.
        ellipse: a tuple (c, a, b, t), where c = (x, y) is the center, a and
            b are the major and minor radii, and t is the rotation angle.

    Returns:
        bool: True if point is in ellipse, False otherwise.
    """

    # if the ellipse is empty, return False
    if ellipse is None:
        return False

    # extract ellipse parameters
    c, a, b, t = ellipse

    # shift point by center of ellipse
    v = point - c

    # rotation matrix, by angle t
    rot_mat = np.array([[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]])

    # matrix F parametrizing ellipse in center form:
    # (x - c)^T * F * (x - c) = 1
    F = rot_mat.T.dot(np.diag(1 / np.array([a, b], dtype=float) ** 2)).dot(rot_mat)

    return v.T.dot(F.dot(v)) <= 1

def welzl(interior, boundary=np.zeros((0, 2))):
    """
    Find the smallest ellipse containing a set of points in the interior and
    another on its boundary. To find the smallest ellipse
    containing a set of points without given boundary points, the function can
    be called with the second argument empty (default usage).

    Arguments:
        interior: an array containing points in R2 as row vectors, representing
            points to be contained in the desired ellipse.
        boundary: an array containing points in R2 as row vectors, representing
            points to be on the boundary of the desired ellipse.

    Returns:
        an ellipse given by a tuple (c, a, b, t), where c = (x, y) is the
            centre, a and b are the major and minor radii, and t is the
            rotation angle.
    """    
    # Stopping condition: stop either if the interior is empty or if there are 5
    # points on the boundary, in which case a unique ellipse is determined.
    if interior.shape[0] == 0 or boundary.shape[0] == 5:
        
        # if the boundary has only 2 points, an ellipse is degenerate
        if boundary.shape[0] <= 2:
            return None
        
        # Call primitive functions to compute the smallest ellipse going through
        # 3, 4, or 5 points.
        elif boundary.shape[0] == 3:
            return ellipse_from_boundary3(boundary)
        elif boundary.shape[0] == 4:
            return ellipse_from_boundary4(boundary)
        else:
            return ellipse_from_boundary5(boundary)

    # choose an arbitrary point in the interior set
    i = np.random.randint(interior.shape[0])
    p = interior[i, :]

    # Remove point from interior set
    interior_wo_p = np.delete(interior, i, 0)
    
    # Recursively call the function to find the smallest ellipse containing
    # the interior points without p
    ellipse = welzl(interior_wo_p, boundary) # <-- Ends here first. This call only returns something upon a fail of the above if statement
    
    # If p is in this ellipse, then this ellipse is also the smallest ellipse containing the interior
    if is_in_ellipse(p, ellipse):
        return ellipse #<--- Loop ends here when all of the points are in the ellipse

    # If not, then p must be on the boundary of the smallest ellipse
    else:
        ellipse = welzl(interior_wo_p, np.vstack([boundary, p]))
        return ellipse

def convexhull(images):
    """
    Read in channelized (in frequency) images outputted by wsclean and
    calculate the convex hull vertices, these will then be used to find the
    smallest enclosing ellipse for beam homogenization
    
    Arguments:
        images: a string or list of strings containing image names to make the convex hull
    Returns:
        hull_points: a convexhull of points for minimum ellipse fitting
    """

    c = []
    a = []
    b = []
    p = []

    # Append beam shapes, noting that the beam is given an angle between 0 and 2pi measure CCW from north
    for im in sorted(images):
        header = fits.getheader(im)
        a.append(header['BMAJ'] / (8 * np.log(2)) ** 0.5)
        b.append(header['BMIN'] / (8 * np.log(2)) ** 0.5)
        p.append(np.radians(header['BPA']) + np.pi * 0.5)

    # Generate a 100 point ellipse for each image
    all_points = []
    for ellipse in np.array([a, b, p]).T:
        all_points.append(sample_ellipse(ellipse, endpoint=False))
    all_points = np.vstack(all_points)

    # Get and return the ConvexHull
    hull = ConvexHull(all_points)
    hull_points = all_points[hull.vertices]
        
    return hull_points, all_points


def get_homogenized_beam(identifier):
    """
    Function to run the Welzl algorithmn to get the smallest enclosing
    ellipse. Steps include: (i) get convex hull; (ii) run welzl on convex hull; 
    (iii) plot enclosing ellipse.
    Arguments:
        identifier: a string pointing to a set of standard WSCLEAN images
    Returns:
        ellipse: semi-major, semi-minor, position angle of the restoring ellipse 
    """

    images = glob.glob(f'{identifier}*[!MFS]-psf.fits')

    # Get hull points
    hull_points, all_points = convexhull(images)

    # Get Smallest confining ellipse
    ellipse = welzl(hull_points)

    # Plot the homogenized beam to make sure it succesfully enclosed the points
    fig, ax = plt.subplots()
    fig.tight_layout()
    ax.set_ylabel('Dec')
    ax.set_xlabel('RA')
    ax.scatter(hull_points[:,0], hull_points[:,1], color='r', label='Convex hull points', zorder=10000, marker='o', s = 8)
    ax.scatter(all_points[:,0], all_points[:,1], color='g', label ='All beam-sampled points', marker='s', s = 8)
    ax.set_xlim(np.amin(hull_points) * 1.3, np.amax(hull_points) * 1.3)
    ax.set_ylim(np.amin(hull_points) * 1.3, np.amax(hull_points) * 1.3)
    plot_ellipse(ellipse[1:], str="k--")
    ax.legend()
    plt.savefig(f'{identifier.replace("IMAGES/","/").replace("INTERVALS/","/")}_ellipse.png')
    plt.clf()
    plt.close()

    # return ellipse paramters
    msg(f'Centre of ellipse is {ellipse[0]}; this should be very close to zero!')
    return ellipse[1:]
    

# Functions to fix names and actually run spatial homogenization
def get_identifiers(prefix):
    """
    Read in the input prefix and return the image name identifiers
    for each image (or set of images).

    Arguments:
        prefix: a string of the input prefix for the wsclean imaging call
    Returns:
        identifiers: an array containing the identifiers, where each identifier
        corresponds to a set of channelized images where we will calculate
        a single beam and homogenize the resolution
    """

    
    # Get identifiers should
    msg('Getting naming identifier "groups"')
    image_arr = sorted(glob.glob(f'{prefix}-t*'))
    if  image_arr != [] and not cfg.WSC_HOMOGENIZETIME:
        suffix = np.unique([im.split(f'{prefix}-t')[-1].split('-')[0] for im in image_arr])
        return [f'{prefix}-t{s}' for s in suffix] , [f'{prefix}-t{s}' for s in suffix]
    elif image_arr != [] and cfg.WSC_HOMOGENIZETIME:
        suffix = np.unique([im.split(f'{prefix}-t')[-1].split('-')[0] for im in image_arr])
        return [prefix] , [f'{prefix}-t{s}' for s in suffix]
    else:
        return [prefix], [prefix]

def homogenize_images(identifier, beam):
    """
    Homogenize the synthesized beam across a set of radio interferometry images.
    
    This function takes a collection of images with varying beam sizes and convolves
    them to a common, larger beam size to enable direct comparison and analysis.
    The process involves:
    1. Creating convolution kernels to transform each original beam to the target beam
    2. Applying these kernels to both model and residual components of each image
    3. Generating homogenized versions of all images with consistent resolution
    4. Creating pseudo-MFS (Multi-Frequency Synthesis) images if not already present
    
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
    
    Examples
    --------
    >>> # Homogenize images to a common 5"x3" beam at 45° position angle
    >>> target_beam = (5.0/3600, 3.0/3600, np.radians(45))
    >>> homogenize_images('1934-638_L_IMAGES/1934-638-t001', target_beam)
    
    Raises
    ------
    FileNotFoundError
        If required PSF or image files are not found
    subprocess.CalledProcessError
        If pypher convolution kernel generation fails
    """

    # Split out the beam components
    a, b, p = beam[0], beam[1], beam[2]
            
    # Get list of image names
    psfs = sorted(glob.glob(f'{identifier}*-psf.fits'))
    good_images = []
    
    # Iterate through psfs (these are Stokes independant)
    for psf in psfs[:]:

        # Define various names for psf
        psf_header = fits.getheader(psf)
        psf_zoom = psf.replace('psf.fits', 'psf.zoom.fits') 
        psf_new   = psf.replace('psf.fits', 'psf.new.fits') 
        kernel     = psf.replace('psf.fits', 'psf.kernel.fits') 

        # If kernel exists skip this psf
        if os.path.exists(kernel):
            msg(f'Skipping {psf} as homogenization kernel already exists; delete {kernel} to re-run')
            continue
        
        # Get sky to pixel and FWHM to SIGMA conversions
        sky_to_pix    = (abs(psf_header.get('CDELT1'))) ** (-1)
        fwhm_to_sig = (8 * np.log(2)) ** (-0.5)

        # Make PSF images
        psf_size = 101      

        # This will check if channel is flagged (WSCLEAN assigns the BMAJ as 0 for these channels)
        if psf_header['BMAJ'] > 1e-14:

            good_images.append(True)

            # Original PSF
            psf_zoom_header = psf_header.copy()
            a_pix = psf_zoom_header['BMAJ'] * sky_to_pix * fwhm_to_sig
            b_pix = psf_zoom_header['BMIN'] * sky_to_pix * fwhm_to_sig
            p_pix = np.radians(psf_zoom_header['BPA']) + 0.5 * np.pi
        
            psf_zoom_header['NAXIS1'] = psf_size
            psf_zoom_header['NAXIS2'] = psf_size
            psf_zoom_image = Gaussian2DKernel(x_stddev = a_pix, y_stddev = b_pix, theta = p_pix, x_size = psf_size , y_size = psf_size, mode='center').array
            psf_zoom_image /= np.amax(psf_zoom_image)

            psf_zoom_hdul = fits.PrimaryHDU(data=psf_zoom_image, header=psf_zoom_header)
            psf_zoom_hdul.writeto(psf_zoom, overwrite=True)

            # New PSF
            a_pix = a * sky_to_pix
            b_pix = b * sky_to_pix

            psf_new_header = psf_header.copy()
            psf_new_header['NAXIS1'] = psf_size
            psf_new_header['NAXIS2'] = psf_size
            psf_new_header['BMAJ'] = a / fwhm_to_sig
            psf_new_header['BMIN'] = b / fwhm_to_sig
            psf_new_header['BPA'] = np.degrees(p - 0.5 * np.pi)
            psf_new_image = Gaussian2DKernel(x_stddev=a_pix, y_stddev = b_pix, theta = p, x_size = psf_size , y_size = psf_size, mode='center').array
            psf_new_image /= np.amax(psf_new_image)

            psf_new_hdul = fits.PrimaryHDU(data=psf_new_image, header=psf_new_header)
            psf_new_hdul.writeto(psf_new, overwrite=True)

            # Run pypher to generate a homogenization kernel
            subprocess.run(["pypher {} {} {}".format(psf_zoom, psf_new, kernel)], shell=True)        

            # Iterate through the images
            images = glob.glob(psf.split('-psf.fits')[0] + '*-image.fits')
            for im in images:

                # Get various names
                image_name =  im.replace('image.fits', 'image.homogenized.fits') 
                residual = im.replace('image.fits', 'residual.fits')
                residual_name = im.replace('image.fits', 'residual.homogenized.fits') 
                model    = im.replace('image.fits', 'model.fits')

                # Convolve model image
                image_new = scipy.signal.fftconvolve(get_image(model), psf_new_image, mode='same')

                # Apply kernel to residual 
                image_rms = scipy.signal.fftconvolve(get_image(residual), get_image(kernel), mode='same')
                image_new += image_rms

                # Save outputs
                header = fits.getheader(im)
                header['BMAJ'] = a / fwhm_to_sig
                header['BMIN'] = b / fwhm_to_sig
                header['BPA'] = np.degrees(p - 0.5 * np.pi)

                image_fits = fits.PrimaryHDU(data=image_new, header=header)
                image_fits.writeto(image_name, overwrite=True)

                residual_fits = fits.PrimaryHDU(data=image_rms, header=header)
                residual_fits.writeto(residual_name, overwrite=True)

        else:
            good_images.append(False)

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
            z = 0
            data = []
            freq  = []
            images = sorted(glob.glob(f'{identifier}-[!MFS]*{stoke}-image.homogenized.fits'))
            for k, im in enumerate(images[:]):
                freq.append(fits.getheader(im)['CRVAL3'])
                if z == 0:
                    header = fits.getheader(im)
                    z += 1

            # Adopt median values for each pixel and output MFS image
            msg(f'Computing median for {len(images)} images using chunked processing')
            data = compute_median_chunked(images)  # Will auto-calculate optimal chunk size
 
            header['CRVAL3'] = np.nanmean(freq)
            mfs_name = f'{identifier}-MFS{stoke}-image.homogenized.fits'
            create_fits(data, header, mfs_name)

def calculate_optimal_chunk_size(ny, nx, num_images, max_memory_pixels=10240*10240*32):
    """
    Calculate optimal chunk size for memory-efficient median computation.
    
    This function determines chunk sizes that are clean divisors of the image 
    dimensions while staying within memory constraints. The chunk size is 
    calculated to ensure total memory usage doesn't exceed the equivalent of
    processing 10240x10240x32 images simultaneously.
    
    Parameters
    ----------948
    ny : int
        Image height in pixels
    nx : int
        Image width in pixels  
    num_images : int
        Number of images to process simultaneously
    max_memory_pixels : int, optional
        Maximum memory usage in pixel equivalents. Default is 10240*10240*32
        
    Returns
    -------
    chunk_y : int
        Optimal chunk height (clean divisor of ny)
    chunk_x : int
        Optimal chunk width (clean divisor of nx)
    estimated_memory_pixels : int
        Estimated memory usage in pixel equivalents
    """
    
    # Include overhead factor for median computation and temporary arrays
    overhead_factor = 1.5
    
    # Calculate maximum chunk area we can afford
    max_chunk_pixels = max_memory_pixels / (overhead_factor * (num_images + 1))
    
    msg(f'Target max chunk area: {max_chunk_pixels:.0f} pixels')
    msg(f'Image dimensions: {ny}x{nx}, Number of images: {num_images}')
    
    # Generate possible divisors for each dimension
    def get_divisors(n):
        """Get all divisors of n in descending order"""
        divisors = []
        for i in range(1, int(np.sqrt(n)) + 1):
            if n % i == 0:
                divisors.append(n // i)  # Larger divisor first
                if i != n // i:
                    divisors.append(i)
        return sorted(divisors, reverse=True)
    
    y_divisors = get_divisors(ny)
    x_divisors = get_divisors(nx)
    
    msg(f'Y divisors: {y_divisors[:10]}...' if len(y_divisors) > 10 else f'Y divisors: {y_divisors}')
    msg(f'X divisors: {x_divisors[:10]}...' if len(x_divisors) > 10 else f'X divisors: {x_divisors}')
    
    # Find the best chunk size combination, favoring square chunks
    best_chunk_y, best_chunk_x = 64, 64  # Minimum chunk size
    best_chunk_area = best_chunk_y * best_chunk_x
    best_square_score = float('inf')  # Lower is better for square-ness
    
    for chunk_y in y_divisors:
        for chunk_x in x_divisors:
            chunk_area = chunk_y * chunk_x
            
            # Skip if chunk is too small (inefficient) or too large (memory)
            if chunk_area < 64*64 or chunk_area > max_chunk_pixels:
                continue
            
            # Calculate "square score" - how close to square the chunk is
            # Perfect square has score = 1.0
            aspect_ratio = max(chunk_y, chunk_x) / min(chunk_y, chunk_x)
            square_score = aspect_ratio
            
            # Prioritize chunks that are:
            # 1. Large enough to be efficient
            # 2. More square-like (lower aspect ratio)
            # 3. Fit within memory constraints
            
            is_better = False
            
            # If this chunk is significantly more square, prefer it
            if square_score < best_square_score * 0.8:  # 20% better square-ness
                is_better = True
            # If square-ness is similar, prefer larger area
            elif abs(square_score - best_square_score) / best_square_score < 0.2:  # Within 20%
                if chunk_area > best_chunk_area:
                    is_better = True
            
            if is_better:
                best_chunk_y = chunk_y
                best_chunk_x = chunk_x
                best_chunk_area = chunk_area
                best_square_score = square_score
    
    # Ensure we don't exceed image dimensions
    best_chunk_y = min(best_chunk_y, ny)
    best_chunk_x = min(best_chunk_x, nx)
    
    # Calculate actual memory usage
    estimated_memory_pixels = int(overhead_factor * (num_images + 1) * best_chunk_y * best_chunk_x)
    
    # Calculate number of chunks needed
    n_chunks_y = ny // best_chunk_y
    n_chunks_x = nx // best_chunk_x
    total_chunks = n_chunks_y * n_chunks_x
    
    msg(f'Optimal chunk size: {best_chunk_y}x{best_chunk_x} pixels')
    msg(f'Chunk area: {best_chunk_area} pixels')
    msg(f'Aspect ratio: {max(best_chunk_y, best_chunk_x)/min(best_chunk_y, best_chunk_x):.2f} (1.0 = perfect square)')
    msg(f'Number of chunks: {n_chunks_y}x{n_chunks_x} = {total_chunks} total')
    msg(f'Estimated memory usage: {estimated_memory_pixels} pixel equivalents')
    msg(f'Memory efficiency: {estimated_memory_pixels/max_memory_pixels*100:.1f}% of maximum')
    
    # Verify the chunks divide evenly
    assert ny % best_chunk_y == 0, f"Chunk height {best_chunk_y} doesn't divide image height {ny}"
    assert nx % best_chunk_x == 0, f"Chunk width {best_chunk_x} doesn't divide image width {nx}"
    
    return best_chunk_y, best_chunk_x, estimated_memory_pixels

def compute_median_chunked(images, chunk_size=None, max_images=512):
    """Compute median pixel-wise using spatial chunks to reduce memory usage"""
    
    # Subsample images if we have too many - include first/last + evenly spaced
    num_total_images = len(images)
    if num_total_images > max_images:
        if num_total_images == 2:
            # Special case: only 2 images, use both
            selected_images = images
        else:
            # Always include first and last
            selected_indices = [0]  # First image
            
            # Calculate how many middle images we can take
            remaining_slots = max_images - 2  # Reserve slots for first and last
            
            if remaining_slots > 0 and num_total_images > 2:
                # Space middle images evenly between first and last
                middle_start = 1
                middle_end = num_total_images - 1
                middle_span = middle_end - middle_start
                
                if remaining_slots >= middle_span:
                    # Can take all middle images
                    middle_indices = list(range(middle_start, middle_end))
                else:
                    # Need to subsample middle images
                    step = middle_span / remaining_slots
                    middle_indices = [int(middle_start + i * step) for i in range(remaining_slots)]
                
                selected_indices.extend(middle_indices)
            
            # Always include last image
            selected_indices.append(num_total_images - 1)
            
            # Remove duplicates and sort
            selected_indices = sorted(list(set(selected_indices)))
            
            selected_images = [images[i] for i in selected_indices]
        
        msg(f'Subsampled from {num_total_images} to {len(selected_images)} images')
        msg(f'Using images at indices: {selected_indices[:5]}...{selected_indices[-5:] if len(selected_indices) > 10 else selected_indices[5:]}')
        images = selected_images
    else:
        msg(f'Using all {num_total_images} images')
        
    # Get image dimensions from first image
    first_img = get_image(images[0])
    ny, nx = first_img.shape
    num_images = len(images)
    result = np.zeros_like(first_img)
    
    # Calculate optimal chunk size if not provided
    if chunk_size is None:
        chunk_y, chunk_x, estimated_memory = calculate_optimal_chunk_size(ny, nx, num_images)
        msg(f'Using optimal chunk size: {chunk_y}x{chunk_x}')
    else:
        # Use provided square chunk size with safety checks
        msg(f'Using provided chunk size: {chunk_size}x{chunk_size}')
        
        # Adjust chunk sizes to fit image dimensions evenly
        # This ensures we don't leave out any pixels
        n_chunks_y = max(1, int(np.ceil(ny / chunk_size)))
        n_chunks_x = max(1, int(np.ceil(nx / chunk_size)))
        
        # Recalculate actual chunk sizes to fit image dimensions
        chunk_y = int(np.ceil(ny / n_chunks_y))
        chunk_x = int(np.ceil(nx / n_chunks_x))
        
        msg(f'Adjusted chunk size to fit image: {chunk_y}x{chunk_x}')
    
    # Calculate number of chunks (should divide evenly with optimal sizing)
    n_chunks_y = ny // chunk_y if chunk_size is None else max(1, int(np.ceil(ny / chunk_y)))
    n_chunks_x = nx // chunk_x if chunk_size is None else max(1, int(np.ceil(nx / chunk_x)))
    
    msg(f'Processing {ny}x{nx} image in {chunk_y}x{chunk_x} chunks')
    msg(f'Total chunks: {n_chunks_y} x {n_chunks_x} = {n_chunks_y * n_chunks_x}')
    
    total_chunks = n_chunks_y * n_chunks_x
    processed_chunks = 0
    
    # Process in chunks with safety checks for pixel coverage
    for i in range(n_chunks_y):
        if chunk_size is None:
            # Optimal sizing - exact division
            y_start = i * chunk_y
            y_end = y_start + chunk_y
        else:
            # Provided sizing - adaptive with bounds checking
            y_start = i * chunk_y
            y_end = min(y_start + chunk_y, ny)  # Ensure we don't exceed image bounds
        
        for j in range(n_chunks_x):
            if chunk_size is None:
                # Optimal sizing - exact division
                x_start = j * chunk_x
                x_end = x_start + chunk_x
            else:
                # Provided sizing - adaptive with bounds checking
                x_start = j * chunk_x
                x_end = min(x_start + chunk_x, nx)  # Ensure we don't exceed image bounds
            
            # Pre-allocate numpy array for chunk data - FIXED FOR MEMORY EFFICIENCY
            chunk_height = y_end - y_start
            chunk_width = x_end - x_start
            chunk_data = np.zeros((num_images, chunk_height, chunk_width), dtype=np.float32)
            
            # Load chunk from all images directly into pre-allocated array
            for img_idx, im in enumerate(images[:]):
                img = get_image(im)
                chunk_data[img_idx] = img[y_start:y_end, x_start:x_end]
            
            # Compute median for this chunk along axis 0 (image dimension)
            result[y_start:y_end, x_start:x_end] = np.median(chunk_data, axis=0)
            
            # Progress reporting
            processed_chunks += 1
            if processed_chunks % 1 == 0 or processed_chunks == total_chunks:
                progress = 100 * processed_chunks / total_chunks
                memory_info = psutil.virtual_memory()
                msg(f'Progress: {processed_chunks}/{total_chunks} chunks ({progress:.1f}%) | '
                    f'Memory: {memory_info.percent:.1f}% used')
            
            # Free memory immediately after processing each chunk
            del chunk_data
            
            # Force garbage collection every few chunks
            if processed_chunks % 5 == 0:
                gc.collect()    
    # Verify we processed all pixels
    msg(f'Processed all {ny * nx} pixels successfully')
    return result

def main():

    # Read in prefix return error if missing it
    if len(sys.argv) != 2:
        msg('ERROR: Missing image prefix')
        sys.exit()
    prefix = sys.argv[-1]

    # Correct naming conventions
    beam_identifiers, image_identifiers = get_identifiers(prefix)
    beams = []

    # Get the beam for each identifier
    for identifier in beam_identifiers:
        msg(f'Solving for smallest enclosing ellipse for image set given by: {identifier}')
        beams.append(get_homogenized_beam(identifier))

    msg('Extracted homogenized beam')
    for k, identifier in enumerate(image_identifiers[:]):    
        msg(f'Homogenizing beam of image set given by: {identifier}')
        if len(beams) == 1:
            beam = beams[0]
        else:
            beam = beams[k]
        homogenize_images(identifier, beam)

            
if __name__ in '__main__':
    main()
