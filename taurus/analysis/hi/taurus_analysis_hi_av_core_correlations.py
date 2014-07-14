#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the taurus molecular cloud.
'''

import pyfits as pf
import numpy as np


''' Plotting Functions
'''

''' Calculations
'''

def calculate_nhi(cube=None, velocity_axis=None, velocity_range=[],
        return_nhi_error=True, noise_cube=None,
        velocity_noise_range=[90,100], Tsys=30., header=None,
        fits_filename=None, fits_error_filename=None, verbose=True):

    ''' Calculates an N(HI) image given a velocity range within which to
    include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to
    fits_filename : str
        If specified, and a header is provided, the nhi image will be written.
    header : pyfits.Header
        Header from cube.

    '''

    import numpy as np

    # Calculate NHI from cube if set
    if cube is not None and velocity_axis is not None:
        image = np.empty((cube.shape[1],
                          cube.shape[2]))
        image[:,:] = np.NaN
        indices = np.where((velocity_axis > velocity_range[0]) & \
                (velocity_axis < velocity_range[1]))[0]
        image[:,:] = cube[indices,:,:].sum(axis=0)
        # Calculate image error
        if return_nhi_error:
            image_error = np.empty((cube.shape[1],
                              cube.shape[2]))
            image_error[:,:] = np.NaN
            image_error[:,:] = (noise_cube[indices,:,:]**2).sum(axis=0)**0.5

    # NHI in units of 1e20 cm^-2
    nhi_image = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2
    nhi_image = image * 1.823e-2

    if fits_filename is not None and header is not None:
        if verbose:
            print('Writing N(HI) image to FITS file %s' % fits_filename)
        header['BUNIT'] = '1e20 cm^-2'
        header.remove('CDELT3')
        header.remove('CRVAL3')
        header.remove('CRPIX3')
        header.remove('CTYPE3')
        header.remove('NAXIS3')
        header['NAXIS'] = 2

        pf.writeto(fits_filename, image*1.823e-2, header = header, clobber =
                True, output_verify = 'fix')

    if fits_error_filename is not None and header is not None:
        if verbose:
            print('Writing N(HI) error image to FITS file %s' % fits_filename)

        pf.writeto(fits_error_filename, image_error * 1.823e-2, header =
                header, clobber = True, output_verify = 'fix')

    if return_nhi_error:
        nhi_image_error = np.ma.array(image_error,
                mask=np.isnan(image_error)) * 1.823e-2
        nhi_image_error = image_error * 1.823e-2
        return nhi_image, nhi_image_error
    else:
        return nhi_image

def calculate_noise_cube(cube=None, velocity_axis=None,
            velocity_noise_range=[-110,-90,90,110], header=None, Tsys=30.,
            filename=None):

    """ Calcuates noise envelopes for each pixel in a cube
    """

    import numpy as np
    import pyfits as pf

    noise_cube = np.zeros(cube.shape)
    for i in range(cube.shape[1]):
        for j in range(cube.shape[2]):
            profile = cube[:,i,j]
            noise = calculate_noise(profile, velocity_axis,
                    velocity_noise_range)
            #noise = 0.1 # measured in line free region
            noise_cube[:,i,j] = calculate_noise_scale(Tsys,
                    profile, noise=noise)

    if filename is not None:
        pf.writeto(filename, noise_cube, header=header)

    return noise_cube

def calculate_noise(profile, velocity_axis, velocity_range):
    """ Calculates rms noise of Tile profile given velocity ranges.
    """
    import numpy as np

    std = 0

    # calculate noises for each individual region
    for i in range(len(velocity_range) / 2):
        velMin = velocity_range[2*i + 0]
        velMax = velocity_range[2*i + 1]

        noise_region = np.where((velocity_axis >= velMin) & \
                        (velocity_axis <= velMax))

        std += np.std(profile[noise_region])

    std /= len(velocity_range) / 2
    return std

def calculate_noise_scale(Tsys, profile, noise=None):
    """ Creates an array for scaling the noise by (Tsys + Tb) / Tsys
    """
    import numpy as np
    n = np.zeros(len(profile))
    n = (Tsys + profile) / Tsys * noise

    return n

def correlate_hi_av(hi_cube=None, hi_velocity_axis=None,
        hi_noise_cube=None, av_image=None, av_image_error=None,
        velocity_centers=None, velocity_widths=None, return_correlations=True):

    '''
    Parameters
    ----------

    Returns
    -------
    hi_vel_range : tuple
        Lower and upper bound of HI velocity range in km/s which provides the
        best correlated N(HI) distribution with Av.
    correlations : array-like, optional
        Array of Pearson correlation coefficients corresponding to each
        permutation through the velocity centers and velocity widths.

    '''

    import numpy as np
    from scipy.stats import pearsonr

    # calculate the velocity ranges given a set of centers and widths
    velocity_ranges = np.zeros(shape=[len(velocity_centers) * \
            len(velocity_widths),2])
    count = 0
    for i, center in enumerate(velocity_centers):
        for j, width in enumerate(velocity_widths):
            velocity_ranges[count,0] = center - width/2.
            velocity_ranges[count,1] = center + width/2.
            count += 1

    # calculate the correlation coefficient for each velocity range
    correlations = np.zeros(velocity_ranges.shape[0])
    pvalues = np.zeros(velocity_ranges.shape[0])

    for i, velocity_range in enumerate(velocity_ranges):
        nhi_image_temp, nhi_image_error = calculate_nhi(cube=hi_cube,
                velocity_axis=hi_velocity_axis,
                velocity_range=velocity_range,
                noise_cube=hi_noise_cube)

        nhi_image = np.ma.array(nhi_image_temp,
                                mask=np.isnan(nhi_image_temp))

        # Select pixels with Av > 1.0 mag and Av_SNR > 5.0.
        # Av > 1.0 mag is used to avoid too low Av.
        # 1.0 mag corresponds to SNR = 1 / 0.2 ~ 5
        # (see Table 2 of Ridge et al. 2006).
        indices = np.where((nhi_image_temp == nhi_image_temp) & \
                (av_image == av_image))

        nhi_image_corr = nhi_image_temp[indices]
        av_image_corr = av_image[indices] #- 0.8 # subtract background of 0.8
        # Use Pearson's correlation test to compare images
        correlations[i] = pearsonr(nhi_image_corr.ravel(),
                av_image_corr.ravel())[0]

        # Shows progress each 10%
        #total = float(velocity_ranges.shape[0])
        #abs_step = int((total * 1)/100) or 1
        #if i and not i % abs_step:
        #    print "{0:.2%} processed".format(i/total)

    correlations = np.ma.array(correlations,
            mask=(correlations != correlations))

    # Find best-correlating velocity range
    best_corr = correlations.max()
    best_corr_index = np.where(correlations == best_corr)
    best_corr_vel_range = velocity_ranges[best_corr_index][0]
    best_corr_vel_range = best_corr_vel_range.tolist()

    if not return_correlations:
    	return best_corr_vel_range, best_corr
    else:
    	return best_corr_vel_range, best_corr, correlations

''' DS9 Region and Coordinate Functions
'''

def convert_core_coordinates(cores, header):

    for core in cores:
        cores[core].update({'box_pixel': 0})
        cores[core].update({'center_pixel': 0})

        #box_wcs = cores[core]['box_wcs']
        #box_pixel = len(box_wcs) * [0,]
        center_wcs = cores[core]['center_wcs']

        # convert centers to pixel coords
        cores[core]['center_pixel'] = get_pix_coords(ra=center_wcs[0],
                                                     dec=center_wcs[1],
                                                     header=header)[:2].tolist()
        # convert box corners to pixel coords
        #for i in range(len(box_wcs)/2):
        #    pixels = get_pix_coords(ra=box_wcs[2*i], dec=box_wcs[2*i + 1],
        #            header=header)
        #    box_pixel[2*i], box_pixel[2*i + 1] = int(pixels[0]), int(pixels[1])
        #cores[core]['box_pixel'] = box_pixel

    return cores

def load_fits(filename,return_header=False):
    ''' Loads a fits file.
    '''

    import pyfits as pf

    f = pf.open(filename)
    if return_header:
        return f[0].data,f[0].header
    else:
        return f[0].data

def get_sub_image(image, indices):

    return image[indices[1]:indices[3],
            indices[0]:indices[2]]

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec), or Ra in degrees
    and dec in degrees.
    '''

    import pywcsgrid2 as wcs
    import pywcs

    # convert to degrees if ra and dec are array-like
    try:
        if len(ra) == 3 and len(dec) == 3:
            ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)
        else:
            raise ValueError('RA and Dec must be in (hrs,min,sec) and' + \
                    ' (deg,arcmin,arcsec) or in degrees.')
    except TypeError:
        ra_deg, dec_deg = ra, dec

    wcs_header = pywcs.WCS(header)
    pix_coords = wcs_header.wcs_sky2pix([[ra_deg, dec_deg, 0]], 0)[0]

    return pix_coords

def hrs2degs(ra=None, dec=None):
    ''' Ra and dec tuples in hrs min sec and deg arcmin arcsec.
    '''

    ra_deg = 15*(ra[0] + ra[1]/60. + ra[2]/3600.)
    dec_deg = dec[0] + dec[1]/60. + dec[2]/3600.

    return (ra_deg, dec_deg)

def read_ds9_region(filename):

    ''' Converts DS9 region file into format for plotting region.

    Need the following format:
        angle : degrees
        xy : pixels
        width : pixels
        height : pixels

    Region file provides following format:
        # Region file format: DS9 version 4.1
        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        fk5
        box(4:17:04.740,+29:20:31.32,5854.33",11972.7",130) # text={test}

    pyregion module reads DS9 regions:
    http://leejjoon.github.io/pyregion/users/overview.html


    '''

    # Import external modules
    import pyregion as pyr

    # Read region file
    try:
        region = pyr.open(filename)
    except IOError:
        return None

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    return region[0].coord_list

def load_ds9_region(cores, filename_base = 'taurus_av_boxes_', header=None):

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    for core in cores:
        region = read_ds9_region(filename_base + core + '.reg')
        if region is not None:
            box_center_pixel = get_pix_coords(ra = region[0],
                                              dec = region[1],
                                              header = header)
            box_center_pixel = (int(box_center_pixel[1]),
                    int(box_center_pixel[0]))
            box_height = region[2] / header['CDELT1']
            box_width = region[3] / header['CDELT2']
            cores[core].update({'box_center_pix': box_center_pixel})
            cores[core].update({'box_width': box_width})
            cores[core].update({'box_height': box_height})
            cores[core].update({'box_angle': region[4]})

    return cores

'''
The main script
'''

def main():

    import grid
    import numpy as np
    import numpy
    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)
    import mygeometry as myg
    reload(myg)
    from mycoords import make_velocity_axis
    import json

    # parameters used in script
    # -------------------------
    # HI velocity integration range
    # Determine HI integration velocity by CO or correlation with Av?
    hi_av_correlation = True
    velocity_centers = np.linspace(-5, 10, 15)
    velocity_widths = np.linspace(5, 20, 15)

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/cores/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/hi/'
    co_dir = '/d/bip3/ezbc/taurus/data/co/'
    core_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'

    # load Planck Av and GALFA HI images, on same grid
    av_data_planck, av_header = load_fits(av_dir + \
                'taurus_av_planck_5arcmin.fits',
            return_header=True)

    av_error_data_planck, av_error_header = load_fits(av_dir + \
                'taurus_av_error_planck_5arcmin.fits',
            return_header=True)

    hi_data, h = load_fits(hi_dir + \
                'taurus_hi_galfa_cube_regrid_planckres.fits',
            return_header=True)

    # make the velocity axis
    velocity_axis = make_velocity_axis(h)

    # Plot NHI vs. Av for a given velocity range
    noise_cube_filename = 'taurus_hi_galfa_cube_regrid_planckres_noise.fits'
    if not path.isfile(hi_dir + noise_cube_filename):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=velocity_axis,
                velocity_noise_range=[90,110], header=h, Tsys=30.,
                filename=hi_dir + noise_cube_filename)
    else:
        noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    # define core properties
    with open(core_dir + 'taurus_core_properties.txt', 'r') as f:
        cores = json.load(f)

    cores = convert_core_coordinates(cores, h)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'taurus_av_boxes_',
            header = h)

    for core in cores:
        print('\nCalculating for core %s' % core)

        # Grab the mask
        mask = myg.get_polygon_mask(av_data_planck,
                cores[core]['box_vertices_rotated'])

        indices = mask == 1

        hi_data_sub = np.copy(hi_data)
        hi_data_sub[:, indices] = np.NaN
        av_data_sub = np.copy(av_data_planck)
        av_error_data_sub = np.copy(av_error_data_planck)

        # Correlate each core region Av and N(HI) for velocity ranges
        hi_vel_range, correlation_coeff, correlations = \
                correlate_hi_av(hi_cube=hi_data_sub,
                    hi_velocity_axis=velocity_axis,
                    hi_noise_cube=noise_cube,
                    av_image=av_data_sub,
                    av_image_error=av_error_data_sub,
                    velocity_centers=velocity_centers,
                    velocity_widths=velocity_widths,
                    return_correlations=True)

        print('HI velocity integration range:')
        print('%.1f to %.1f km/s' % (hi_vel_range[0], hi_vel_range[1]))
        print('With correlation of:')
        print('%.2f' % correlation_coeff)

        cores[core]['hi_velocity_range'] = hi_vel_range
        cores[core]['correlation_coeff'] = correlation_coeff

    with open(core_dir + 'taurus_core_properties.txt', 'w') as f:
        json.dump(cores, f)

if __name__ == '__main__':
    main()


