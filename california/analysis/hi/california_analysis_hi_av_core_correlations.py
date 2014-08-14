#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the california molecular cloud.
'''

import pyfits as pf
import numpy as np


''' Plotting Functions
'''
def plot_correlations(correlations,velocity_centers,velocity_widths,
        savedir='./', filename=None,show=True, returnimage=False):
    ''' Plots a heat map of correlation values as a function of velocity width
    and velocity center.
    '''

    # Import external modules
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 8
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale * 3 / 4.0,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              #'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)


    fig = plt.figure(figsize=(3,2))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="3%",
                 cbar_size='6%',
                 axes_pad=0,
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    correlations_image = np.empty((velocity_centers.shape[0],
                                   velocity_widths.shape[0]))
    correlations_image[:,:] = np.NaN
    count = 0
    try:
        for i, center in enumerate(velocity_centers):
            for j, width in enumerate(velocity_widths):
                correlations_image[i,j] = correlations[count]
                count += 1
    except IndexError:
        print(' plot_correlations: O-d array input, cannot proceed')

    image = np.ma.array(correlations_image, mask=np.isnan(correlations_image))

    ax = imagegrid[0]

    ax.set_xlabel('Velocity Width (km/s)')
    ax.set_ylabel('Velocity Center (km/s)')

    #ax.set_xticks(np.arange(0,velocity_widths.shape[0],1)[::5],
    #        velocity_centers[::5])

    plt.rc('text', usetex=False)
    im = ax.imshow(image, interpolation='nearest', origin='lower',
            extent=[velocity_widths[0],velocity_widths[-1],
                    velocity_centers[0],velocity_centers[-1]],
            cmap=plt.cm.gray)
    cb = ax.cax.colorbar(im)

    #cb.set_clim(vmin=0.)
    # Write label to colorbar
    cb.set_label_text(r'Correlation coefficient')

    fractions = np.array([0.95, 0.85, 0.75, 0.65,])
    levels = fractions * image.max()

    cs = ax.contour(image, levels=levels, origin='lower',
            extent=[velocity_widths[0],velocity_widths[-1],
                    velocity_centers[0],velocity_centers[-1]],
            colors='k'
            )

    # Define a class that forces representation of float to look a certain way
    # This remove trailing zero so '1.0' becomes '1'
    class nf(float):
         def __repr__(self):
             str = '%.1f' % (self.__float__(),)
             if str[-1]=='0':
                 return '%.0f' % self.__float__()
             else:
                 return '%.1f' % self.__float__()

    # Recast levels to new class
    cs.levels = [nf(val) for val in fractions*100.0]

    #fmt = {}
    #for level, fraction in zip(cs.levels, fractions):
    #    fmt[level] = fraction
    fmt = '%r %%'

    ax.clabel(cs, cs.levels, fmt=fmt, fontsize=9, inline=1)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
    	plt.draw()
        plt.show()
    if returnimage:
        return correlations_image

''' Calculations
'''

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
    from scipy.stats import kendalltau
    from myimage_analysis import calculate_nhi
    from scipy import signal

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
        nhi_image_error_corr = nhi_image_error[indices]
        av_image_corr = av_image[indices]
        av_image_error_corr = av_image_error[indices]

        # Normalize images by their mean
        #nhi_image_corr = (nhi_image_corr - nhi_image_corr.mean()) / \
        #                  nhi_image_error_corr
        #av_image_corr = (av_image_corr - av_image_corr.mean()) / \
        #                 av_image_error_corr

        #correlations[i] = np.sum(np.abs(nhi_image_corr - av_image_corr))

        # Use Pearson's correlation test to compare images
        #correlations[i] = pearsonr(nhi_image_corr.ravel(),
        #        av_image_corr.ravel())[0]
        correlations[i], pvalues[i] = kendalltau(nhi_image_corr.ravel(),
                                                 av_image_corr.ravel())

        # Shows progress each 10%
        total = float(correlations.shape[0])
        abs_step = int((total * 1)/100) or 1
        if i and not i % abs_step:
            print "\t{0:.2%} processed".format(i/total)

    #correlations /= correlations.max()

    #correlations = 1.0 / correlations

    correlations = np.ma.array(correlations,
            mask=(correlations != correlations))
    pvalues = np.ma.array(pvalues,
            mask=np.isnan(pvalues))

    plot_correlations(correlations,
                      velocity_centers,
                      velocity_widths,
                      show=0,
                      returnimage=False,
                      savedir='/usr/users/ezbc/Desktop/',
                      filename='correlation.png')

    # Find best-correlating velocity range
    #best_corr = correlations.max()
    #best_corr_index = np.where(correlations == best_corr)
    #best_corr_vel_range = velocity_ranges[best_corr_index][0]
    #best_corr_vel_range = best_corr_vel_range.tolist()

    correlations_image = np.empty((velocity_centers.shape[0],
                                   velocity_widths.shape[0]))
    correlations_image[:,:] = np.NaN
    count = 0
    for i, center in enumerate(velocity_centers):
        for j, width in enumerate(velocity_widths):
            correlations_image[i,j] = correlations[count]
            count += 1

    max_index = np.where(correlations_image == correlations_image.max())

    # Define parameter resolutions
    delta_center = velocity_centers[1] - velocity_centers[0]
    delta_width = velocity_widths[1] - velocity_widths[0]

    center_corr = np.sum(correlations_image, axis=1) * delta_center
    center_corr = correlations_image[:, max_index[1][0]]
    center_confint = threshold_area(velocity_centers,
                                         center_corr,
                                         area_fraction=0.68)

    width_corr = np.sum(correlations_image, axis=0) * delta_width
    width_corr = correlations_image[max_index[0][0], :]
    width_confint = threshold_area(velocity_widths,
                                        width_corr,
                                        area_fraction=0.68)

    print('Velocity widths = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(width_confint[0],
                                                    width_confint[2],
                                                    np.abs(width_confint[1])))
    print('Velocity centers = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(center_confint[0],
                                                    center_confint[2],
                                                    np.abs(center_confint[1])))

    # Write PDF

    center = center_confint[0]
    upper_lim = (center_confint[0] + width_confint[0]/2.)
    lower_lim = (center_confint[0] - width_confint[0]/2.)
    upper_lim_error = (center_confint[2]**2 + width_confint[2]**2)**0.5
    lower_lim_error = (center_confint[1]**2 + width_confint[1]**2)**0.5

    vel_range_confint = (lower_lim, upper_lim, lower_lim_error, upper_lim_error)

    '''
    if not return_correlations:
    	return best_corr_vel_range, best_corr
    else:
    	return best_corr_vel_range, best_corr, correlations
    '''

    if not return_correlations:
        return vel_range_confint
    else:
        return vel_range_confint, correlations, center_corr, width_corr

def threshold_area(x, y, area_fraction=0.68):

    '''
    Finds the limits of a 1D array which includes a given fraction of the
    integrated data.

    Parameters
    ----------
    data : array-like
        1D array.
    area_fraction : float
        Fraction of area.

    Returns
    -------
    limits : tuple
        Lower and upper bound including fraction of area.

    '''

    import numpy as np
    from scipy.integrate import simps as integrate

    # Step for lowering threshold
    step = (np.max(y) - np.median(y)) / 100.0

    # initial threshold
    threshold = np.max(y) - step
    threshold_area = 0.0

    # area under whole function
    area = integrate(y, x)

    # Stop when the area below the threshold is greater than the max area
    while threshold_area < area * area_fraction:

        threshold_indices = np.where(y > threshold)[0]

        try:
            bounds_indices = (threshold_indices[0], threshold_indices[-1])
        except IndexError:
            bounds_indices = ()

        try:
            threshold_area = integrate(y[bounds_indices[0]:bounds_indices[1]],
                                       x[bounds_indices[0]:bounds_indices[1]])
            threshold_area += threshold * (x[bounds_indices[1]] - \
                                           x[bounds_indices[0]])
        except IndexError:
            threshold_area = 0

        threshold -= step

    x_peak = x[y == y.max()][0]
    low_error, up_error = x_peak - x[bounds_indices[0]], \
                          x[bounds_indices[1]] - x_peak

    return (x_peak, low_error, up_error)

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

def load_ds9_region(cores, filename_base = 'california_av_boxes_', header=None):

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
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    # parameters used in script
    # -------------------------
    # HI velocity integration range
    # Determine HI integration velocity by CO or correlation with Av?
    hi_av_correlation = True
    velocity_centers = np.arange(-10, 15, 1)
    velocity_widths = np.arange(1, 60, 2)
    #velocity_centers = np.linspace(-10, 10, 7)
    #velocity_widths = np.linspace(2, 30, 5)
    #velocity_centers = np.linspace(-10, 10, 14)
    #velocity_widths = np.linspace(2, 30, 15)
    #velocity_centers = np.linspace(-20, 20, 15)
    #velocity_widths = np.linspace(2, 40, 15)

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/california/figures/cores/'
    av_dir = '/d/bip3/ezbc/california/data/av/'
    hi_dir = '/d/bip3/ezbc/california/data/hi/'
    co_dir = '/d/bip3/ezbc/california/data/co/'
    core_dir = '/d/bip3/ezbc/california/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/california/data/python_output/ds9_regions/'

    # load Planck Av and GALFA HI images, on same grid
    av_data_planck, av_header = load_fits(av_dir + \
                'california_av_planck_5arcmin.fits',
            return_header=True)

    av_error_data_planck, av_error_header = load_fits(av_dir + \
                'california_av_error_planck_5arcmin.fits',
            return_header=True)

    hi_data, h = load_fits(hi_dir + \
                'california_hi_galfa_cube_regrid_planckres.fits',
            return_header=True)

    # make the velocity axis
    velocity_axis = make_velocity_axis(h)

    # Plot NHI vs. Av for a given velocity range
    noise_cube_filename = 'california_hi_galfa_cube_regrid_planckres_noise.fits'
    if not path.isfile(hi_dir + noise_cube_filename):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=velocity_axis,
                velocity_noise_range=[90,110], header=h, Tsys=30.,
                filename=hi_dir + noise_cube_filename)
    else:
        noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    # define core properties
    with open(core_dir + 'california_core_properties.txt', 'r') as f:
        cores = json.load(f)

    cores = convert_core_coordinates(cores, h)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'california_av_boxes_',
            header = h)

    for core in cores:
        print('\nCalculating for core %s' % core)

        # Grab the mask
        mask = myg.get_polygon_mask(av_data_planck,
                cores[core]['box_vertices_rotated'])

        #indices = mask == 1
        indices = ((mask == 0) &\
                   (av_data_planck / av_error_data_planck > 5))

        #hi_data_sub = np.copy(hi_data)
        #hi_data_sub[:, indices] = np.NaN
        hi_data_sub = np.copy(hi_data[:, indices])
        noise_cube_sub = np.copy(noise_cube[:, indices])
        av_data_sub = np.copy(av_data_planck[indices])
        av_error_data_sub = np.copy(av_error_data_planck[indices])

        # Correlate each core region Av and N(HI) for velocity ranges
        vel_range_confint, correlations, center_corr, width_corr = \
                correlate_hi_av(hi_cube=hi_data_sub,
                                hi_velocity_axis=velocity_axis,
                                hi_noise_cube=noise_cube_sub,
                                av_image=av_data_sub,
                                av_image_error=av_error_data_sub,
                                velocity_centers=velocity_centers,
                                velocity_widths=velocity_widths,
                                return_correlations=True)

        print('HI velocity integration range:')
        print('%.1f to %.1f km/s' % (vel_range_confint[0],
                                     vel_range_confint[1]))

        #print('Median center + width correlations')
        #print(np.median(center_corr))
        #print(np.median(width_corr))

        cores[core]['hi_velocity_range'] = vel_range_confint[0:2]
        cores[core]['hi_velocity_range_error'] = vel_range_confint[2:]
        cores[core]['center_corr'] = center_corr.tolist()
        cores[core]['width_corr'] = width_corr.tolist()
        cores[core]['vel_centers'] = velocity_centers.tolist()
        cores[core]['vel_widths'] = velocity_widths.tolist()
        #cores[core]['correlation_coeff'] = correlation_coeff

    with open(core_dir + 'california_core_properties.txt', 'w') as f:
        json.dump(cores, f)

if __name__ == '__main__':
    main()








