#!/usr/bin/python

''' Calculates the N(HI) map for california

'''

import numpy as np
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')

''' Plotting Functions
'''

def plot_av_vs_nhi(nhi, av, av_error=None, limits=None,
        fit=True, savedir='', filename=None, show=True, fit_params=None,
        contour_plot=True,
        scale=('linear','linear'), title = '', gridsize=(100,100), std=None):

    # import external modules
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from matplotlib import cm
    from astroML.plotting import scatter_contour

    # set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    #plt.rcdefaults()

    # color map
    cmap = plt.cm.gnuplot

    # color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    params = {'axes.color_cycle': color_cycle, # colors of different plots
             }
    #plt.rcparams.update(params)

    # Create figure instance
    fig = plt.figure(figsize=(3.6, 3.6))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1, 1),
                 ngrids=1,
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 #cbar_mode='single',
                 cbar_pad=0.1,
                 cbar_size=0.2,
                 )

    # Drop the NaNs from the images
    if type(av_error) is float:
        indices = np.where((av == av) &\
                           (nhi == nhi)
                           )

    if type(av_error) is np.ndarray or \
            type(av_error) is np.ma.core.MaskedArray:
        indices = np.where((av == av) &\
                           (nhi == nhi) &\
                           (av_error == av_error)
                           )

    av_nonans = av[indices]
    nhi_nonans = nhi[indices]

    # Fix error data types
    if type(av_error) is np.ndarray:
        av_error_nonans = av_error[indices]
    else:
        av_error_nonans = np.array(av_error[indices])

    # Create plot
    ax = imagegrid[0]

    if limits is not None:
        contour_range = ((limits[0], limits[1]),
                         (limits[2], limits[3]))
    else:
        contour_range = limits

    print contour_plot
    if contour_plot:
        l1 = scatter_contour(nhi_nonans.ravel(),
                             av_nonans.ravel(),
                             threshold=2,
                             log_counts=0,
                             levels=5,
                             ax=ax,
                             histogram2d_args=dict(bins=30,
                                    range=contour_range),
                             plot_args=dict(marker='o',
                                            linestyle='none',
                                            color='black',
                                            alpha=0.7,
                                            markersize=2),
                             contour_args=dict(
                                               cmap=plt.cm.gray_r,
                                               #cmap=cmap,
                                               ),
                             )
    else:
        image = ax.errorbar(nhi_nonans.ravel(),
                av_nonans.ravel(),
                xerr=(av_error_nonans.ravel()),
                alpha=0.2,
                color='k',
                marker='^',
                ecolor='k',
                linestyle='None',
                markersize=3
                )

    # Plot sensitivies
    #av_limit = np.median(av_errors[0])
    #ax.axvline(av_limit, color='k', linestyle='--')

    # Plot 1 to 1 pline
    if 1:
        p = np.polyfit(nhi_nonans.ravel(), av_nonans.ravel(), deg=1)
        x_fit = np.linspace(-10, 100, 100)
        y_fit = fit_params['dgr'] * x_fit + fit_params['intercept']
        y_poly_fit = p[0] * x_fit + p[1]
        ax.plot(x_fit,
                y_fit,
                #color='0.5',
                linestyle='--',
                linewidth=2,
                alpha=1)
        ax.plot(x_fit,
                y_poly_fit,
                color='r',
                linestyle='--',
                linewidth=2,
                alpha=1)
        print('fitted intercept = {0:.1f} mag'.format(p[1]))
        std = None
        if std is not None:
            ax.fill_between((-std, 10),
                            (0, 10 + std),
                            (-2*std, 10 - std),
                            color='0.2',
                            alpha=0.2,
                            edgecolor='none'
                            )

    # Annotations
    anno_xpos = 0.95

    ax.set_xscale(scale[0], nonposx = 'clip')
    ax.set_yscale(scale[1], nonposy = 'clip')

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'$N($H$\textsc{i}) \times\,10^{20}$ cm$^{-3}$')
    ax.set_ylabel(r'$A_V$ [mag]')
    #ax.set_title(core_names[i])

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        fig.show()

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
        center_pixel = get_pix_coords(ra=center_wcs[0],
                                      dec=center_wcs[1],
                                      header=header)[:2]
        cores[core]['center_pixel'] = center_pixel.tolist()

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

def main(dgr=None, vel_range=None, vel_range_type='single', region=None,
        av_data_type='planck', use_binned_images=False):
    ''' Executes script.

    Parameters
    ----------
    dgr : float
        If None, pulls best-fit value from properties.
    vel_range : tuple
        If None, pulls best-fit value from properties.
    '''

    # import external modules
    import pyfits as pf
    import numpy as np
    from mycoords import make_velocity_axis
    import mygeometry as myg
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    import json

    # Script parameters
    # -----------------
    if use_binned_images:
        bin_string = '_bin'
    else:
        bin_string = ''

    # Name of noise cube
    noise_cube_filename = \
            'california_hi_galfa_cube_regrid_planckres_noise' + bin_string + \
            '.fits'

    # Name of property files results are written to
    prop_file = 'california_global_properties_' + av_data_type + '_scaled'

    # Regions, regions to edit the global properties with
    if region == 1:
        region_limit = {'wcs' : (((5, 10, 0), (19, 0, 0)),
                                 ((4, 30, 0), (27, 0, 0))),
                          'pixel' : ()
                         }
    elif region == 2:
        region_limit = {'wcs' : (((4, 30, 0), (19, 0, 0)),
                                 ((3, 50, 0), (29, 0, 0))),
                          'pixel' : ()
                        }
    elif region == 3:
        region_limit = {'wcs' : (((4, 30, 0), (29, 0, 0)),
                                 ((3, 50, 0), (33, 0, 0))),
                          'pixel' : ()
                        }
    else:
    	region_limit = None

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/california/figures/av/'
    av_dir = '/d/bip3/ezbc/california/data/av/'
    hi_dir = '/d/bip3/ezbc/california/data/hi/'
    co_dir = '/d/bip3/ezbc/california/data/co/'
    core_dir = '/d/bip3/ezbc/california/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/california/data/python_output/'
    region_dir = '/d/bip3/ezbc/california/data/python_output/ds9_regions/'

    # load Planck Av and GALFA HI images, on same grid
    if av_data_type == 'lee12_2mass':
    	print('\nLoading Lee+12 data...')
        av_image, av_header = load_fits(av_dir + \
                    'california_av_lee12_2mass_regrid_planckres' + bin_string + \
                    '.fits',
                return_header=True)
        av_image_error = 0.1 * np.ones(av_image.shape)
    elif av_data_type == 'lee12_iris':
    	print('\nLoading Lee+12 data...')
        av_image, av_header = load_fits(av_dir + \
                    'california_av_lee12_iris_regrid_planckres' + bin_string + \
                    '.fits',
                return_header=True)
        av_image_error = 0.1 * np.ones(av_image.shape)
    else:
    	print('\nLoading Planck data...')
        av_image, av_header = load_fits(av_dir + \
                    'california_av_planck_5arcmin' + bin_string + \
                    '.fits',
                return_header=True)

        av_image_error, av_error_header = load_fits(av_dir + \
                    'california_av_error_planck_5arcmin' + bin_string + \
                    '.fits',
                return_header=True)

    hi_cube, hi_header = load_fits(hi_dir + \
                'california_hi_galfa_cube_regrid_planckres' + bin_string + \
                '.fits',
            return_header=True)

    hi_noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    if not use_binned_images:
        co_data, co_header = load_fits(co_dir + \
                    'california_co_cfa_cube_regrid_planckres' + bin_string + \
                    '.fits',
                return_header=True)

    # Load global properties of cloud
    # global properties written from script
    # 'av/california_analysis_global_properties.txt'
    if region is not None:
        likelihood_filename += '_region{0:.0f}'.format(region)
        results_filename += '_region{0:.0f}'.format(region)

    print('\nReading global parameter file\n' + prop_file + '.txt')
    with open(property_dir + prop_file + '.txt', 'r') as f:
        props = json.load(f)

    if vel_range is not None:
        props['hi_velocity_range'] = vel_range
    else:
        vel_width = props['hi_velocity_width_max']['value']
        vel_center = np.array(props['hi_velocity_center']['value'])
        vel_center = -4.0
        vel_range = (vel_center - vel_width / 2.0,
                     vel_center + vel_width / 2.0)
    if dgr is not None:
        props['dust2gas_ratio_max']['value'] = dgr
    else:
        dgr = props['dust2gas_ratio_max']['value']
        intercept = props['intercept_max']['value']

    fit_params = {}
    fit_params['dgr'] = dgr
    fit_params['intercept'] = intercept

    # define core properties
    with open(core_dir + 'california_core_properties.txt', 'r') as f:
        cores = json.load(f)

    # make velocity axis for hi cube
    velocity_axis = make_velocity_axis(hi_header)

    if not use_binned_images:
        # make velocity axis for co cube
        co_velocity_axis = make_velocity_axis(co_header)

    # Write core coordinates in pixels
    cores = convert_core_coordinates(cores, hi_header)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'california_av_boxes_',
            header = hi_header)

    # create nhi image
    nhi_image = calculate_nhi(cube=hi_cube,
            velocity_axis=velocity_axis,
            velocity_range=vel_range,
            header=hi_header,
            noise_cube=hi_noise_cube)

    # create model av map
    av_model = nhi_image * dgr

    if vel_range_type == 'single':
        print('\nHI velocity integration range:')
        print('%.1f to %.1f km/s' % (vel_range[0],
                                     vel_range[1]))
    elif vel_range_type == 'multiple':
        print('\nHI velocity integration ranges:')
        for i in xrange(0, vel_range.shape[0]):
            print('%.1f to %.1f km/s' % (vel_range[i, 0],
                                         vel_range[i, 1]))

    print('\nDGR:')
    print('%.2f x 10^-20 cm^2 mag' % (dgr))

    print('\nIntercept:')
    print('%.2f mag' % (intercept))

    # Get mask and mask images
    mask = np.asarray(props['mask' + bin_string])

    mask_images = 1

    if mask_images:
        av_image[mask] = np.nan
        nhi_image[mask] = np.nan
        av_image_error[mask] = np.nan
        av_model[mask] = np.nan

    indices = ((np.isnan(av_model)) & \
               (np.isnan(av_image)) & \
               (np.isnan(av_image_error)))

    if 1:
        import matplotlib.pyplot as plt
        plt.imshow(av_image)
        plt.show()

    print('\nTotal number of pixels after masking = ' + str(props['npix']))

    # Plot
    figure_types = ['png', 'pdf']
    for figure_type in figure_types:
        if region is None:
            filename = 'california_av_vs_nhi_' + av_data_type + bin_string

        filename = figure_dir + filename + '.' + figure_type

        print('\nSaving Av model image to \n' + filename)

        plot_av_vs_nhi(nhi_image,
                av_image,
                av_error=av_image_error,
                #limits=[10**-1, 10**1.9, 10**0, 10**1.7],
                fit_params=fit_params,
                limits=[5,40,-0.2,2],
                #limits=[0,30,0,10],
                gridsize=(10,10),
                #scale=('log', 'log'),
                #scale=('linear', 'linear'),
                filename=filename,
                contour_plot=not use_binned_images,
                std=0.22,
                )

if __name__ == '__main__':
    # Use Planck dust Av map or Kainulainen 2009 optical extinction Av map?

    #main(av_data_type='planck', use_binned_images=False)
    main(av_data_type='planck', use_binned_images=True)
    #main(av_data_type='planck', vel_range=[-18, 25])
    #main(av_data_type='planck', use_binned_images=True)


