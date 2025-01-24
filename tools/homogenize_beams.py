#!/usr/bin/env python
# andrew.hughes@physics.ox.ac.uk

import os
import sys
import glob
import subprocess
import time
import scipy
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from scipy.spatial import ConvexHull

import os.path as o
sys.path.append(o.abspath(o.join(o.dirname(sys.modules[__name__].__file__), "..")))

from oxkat import config as cfg

# Helper functions
def msg(txt):
    stamp = time.strftime(' %Y-%m-%d %H:%M:%S | ')
    print(stamp+txt)

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

        # If kernel exists delete is
        if os.path.exists(kernel):
            os.remove(kernel)

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
            images = sorted(glob.glob(f'{identifier}-[!MFS]*{stoke}-image.fits'))
            for k, im in enumerate(images[:]):
                freq.append(fits.getheader(im)['CRVAL3'])
                if good_images[k]:
                    if z == 0:
                        header = fits.getheader(im)
                        z += 1
                    data.append(get_image(im))

            # Adopt median values for each pixel and output MFS image
            data = np.median(data, axis = 0)
            header['CRVAL3'] = np.nanmean(freq)
            mfs_name = f'{identifier}-MFS{stoke}-image.homogenized.fits'
            mfs_fits = fits.PrimaryHDU(data=data, header=header)
            mfs_fits.writeto(mfs_name, overwrite=True)


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

    for k, identifier in enumerate(image_identifiers[:]):    
        msg(f'Homogenizing beam of image set given by: {identifier}')
        if len(beams) == 1:
            beam = beams[0]
        else:
            beam = beams[k]
        homogenize_images(identifier, beam)


            
if __name__ in '__main__':
    main()
