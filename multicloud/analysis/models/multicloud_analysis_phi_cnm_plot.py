#!/usr/bin/python

'''

Plots Av image for all three clouds.

'''

import pyfits as pf
import numpy as np

''' Plotting Functions
'''

def plot_phi_cnm_vs_glat(glats, phi_cnms, glat_errors=None,
        phi_cnm_errors=None, limits=None, savedir='./', filename=None,
        show=True, scale=['linear', 'linear'], title='', core_names=''):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale * 3 / 4.0,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': True,
              'figure.figsize': (8, 8),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # Create plot
    ax = fig.add_subplot(111)

    phi_cnm_errors = np.copy(phi_cnm_errors).T

    image = ax.errorbar(
                    phi_cnms,
                    glats,
                    xerr=(phi_cnm_errors),
                    alpha=0.5,
                    color='k',
                    marker='^',ecolor='k',linestyle='none',
                    markersize=4
                    )

    ax.set_xscale(scale[0], nonposx = 'clip')
    ax.set_yscale(scale[1], nonposy = 'clip')

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'$\phi_{\rm CNM}$',)
    ax.set_ylabel(r'Galactic Latitude (deg)',)

    if title is not None:
        fig.suptitle(title.capitalize(), fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()

def plot_T_cnm_vs_glat(cloud_dict, limits=None, savedir='./', filename=None,
        show=True, scale=['linear', 'linear'], title='', core_names=''):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 4)]
    font_scale = 9
    line_weight = 600
    font_weight = 600
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'axes.weight': line_weight,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale,
              'xtick.labelsize': font_scale,
              'xtick.weight': line_weight,
              'ytick.labelsize': font_scale,
              'ytick.weight': line_weight,
              'font.weight': font_weight,
              'axes.labelweight': font_weight,
              'text.usetex': True,
              #'font.family': 'sans-serif',
              'figure.figsize': (3.6, 3.6),
              'figure.dpi': 600,
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # Create plot
    ax = fig.add_subplot(111)

    markers = ['^', 's', 'o']

    for i, cloud in enumerate(cloud_dict):

        glats = []
        phi_cnms = []
        phi_cnm_errors = []
        T_cnms = []
        T_cnm_errors = []
        core_plot_count = 0

        for j, core in enumerate(cloud_dict[cloud]['cores']):
            core_dict = cloud_dict[cloud]['cores'][core]
            glat = core_dict['center_wcs']['gal'][1]
            try:
                phi_cnm = core_dict['krumholz_results']['phi_cnm']
                phi_cnm_error = core_dict['krumholz_results']['phi_cnm_error']
                T_cnm = core_dict['krumholz_results']['T_cnm']
                T_cnm_error = np.array(core_dict['krumholz_results']['T_cnm_error'])
                T_cnm_error = np.copy((T_cnm_error,)).T

                if core_plot_count == 0:
                    label = cloud.capitalize()
                else:
                    label = None

                image = ax.errorbar(
                                T_cnm,
                                glat,
                                xerr=T_cnm_error,
                                alpha=0.8,
                                color=color_cycle[i],
                                marker=markers[i],
                                ecolor=color_cycle[i],
                                linestyle='none',
                                markersize=4,
                                label=label,
                                )
                core_plot_count += 1

            except KeyError:
                pass

    ax.set_xscale(scale[0], nonposx = 'clip')
    ax.set_yscale(scale[1], nonposy = 'clip')

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'$T_{\rm CNM}$ [K]',)
    ax.set_ylabel(r'Galactic Latitude [deg]',)

    ax.legend(loc='best', numpoints=1)

    if title is not None:
        fig.suptitle(title.capitalize(), fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        fig.show()

''' Calculations
'''

def create_box(core_pos, width, height, core_rel_pos=0.1):

    '''
    Parameters
    ----------
    core_pos : array-like
        x and y pixel coordinates of core
    width : int, float
        Width of box along x axis.
    height : int, float
        Height of box along y axis.
    core_rel_pos : float, optional
        Core position in box along y-axis as a fraction of the height with the
        origin in the south.

    Returns
    -------
    box_vertices : numpy.array
        4 x 2 array with box pixel vertices, starting from lower-left going
        clockwise.

    '''

    box_vertices = np.empty((4, 2))

    # x-coords
    box_vertices[:2, 0] = core_pos[0] - width / 2.
    box_vertices[2:4, 0] = core_pos[0] + width / 2.

    # y-coords
    offset = height * core_rel_pos
    box_vertices[1:3, 1] = core_pos[1] + height - offset
    box_vertices[(0, 3), 1] = core_pos[1] - offset

    return box_vertices

def rotate_box(box_vertices, anchor, angle):

    '''
    Parameters
    ----------
    box_vertices : numpy.array
        4 x 2 array with box pixel vertices, starting from lower-left going
        clockwise.
    anchor : tuple
        x and y coordinates of pivot point.
    angle : float
        Angle to rotate polygon clockwise from North.

    Returns
    -------
    box_vertices_rotated : numpy.array
        4 x 2 array with rotated box pixel vertices.
    '''

    from mygeometry import rotate_polygon

    box_vertices_rotated = rotate_polygon(box_vertices, anchor, angle)

    return box_vertices_rotated

def get_radial_profile(image, center=None, stddev=False, binsize=1,
        mask=None, weights=None):

    ''' Calculates radial profiles of an image at the center.

    '''

    # import external modules
    import numpy as np
    from scipy.optimize import curve_fit
    from agpy import azimuthalAverage as radial_average

    if stddev and weights is not None:
        weights=None

    result = radial_average(image, binsize=binsize, center=center,
            stddev=stddev, mask=mask, interpnan=False, returnradii=True,
            weights=weights)

    return result

def fit_profile(radii, profile, function, sigma=None):

    ''' Fits a radial profile with a power function A * radius**alpha where A
    and alpha are the fitted constants.
    '''

    # import external modules
    from scipy.optimize import curve_fit
    import numpy as np

    profile_fit = curve_fit(function, radii, profile, sigma=sigma,
            maxfev=1000000,)

    return profile_fit

def write_fitting_results(cloud_dict, filename):

    results_list = []

    for cloud in cloud_dict:
        for core in cloud_dict[cloud]['cores']:
            results_row = []
            core_dict = cloud_dict[cloud]['cores'][core]
            results_row.append(cloud)
            results_row.append(core)
            results_row.append(core_dict['phi_cnm'])
            results_row.append(core_dict['phi_cnm_error'][0])
            results_row.append(core_dict['phi_cnm_error'][1])
            results_row.append(core_dict['phi_mol'])
            results_row.append(core_dict['phi_mol_error'][0])
            results_row.append(core_dict['phi_mol_error'][1])
            results_row.append(core_dict['T_cnm'])
            results_row.append(core_dict['T_cnm_error'][0])
            results_row.append(core_dict['T_cnm_error'][1])
            results_row.append(core_dict['Z'])
            results_row.append(core_dict['Z_error'][0])
            results_row.append(core_dict['Z_error'][1])
            results_row.append(core_dict['center_wcs']['eq'][0])
            results_row.append(core_dict['center_wcs']['eq'][1])
            results_row.append(core_dict['center_wcs']['gal'][0])
            results_row.append(core_dict['center_wcs']['gal'][1])

            results_list.append(results_row)

    results_header = ('{0:s}\t{0:s}\t'.format(cloud, core) + \
                      'phi_cnm\tphi_cnm_herr\tphi_cnm_lerr\tphi_mol\t' + \
                      'phi_mol_herr\tphi_mol_lerr\tT_cnm\tT_cnm_herr\t' + \
                      'T_cnm_lerr\tZ\tZ_herr\tZ_lerr\tRA\tDec\tGlat\tGlong')

    results_array = np.asarray(results_list, dtype='string')

    print('Saving fitting results as {0:s}'.format(filename))

    np.savetxt(filename, results_array, header=results_header, delimiter='\t',
               fmt='%s')

def calc_core_coords(cloud_dict):

    ''' Calculates galactic coordinates for cores.

    '''

    for cloud in cloud_dict:
        #cloud_dict[cloud] = convert_core_coordinates(cloud_dict[cloud]['cores'],
        #                                cloud_dict[cloud]['av_header'])

        for core in cloud_dict[cloud]['cores']:
            core_dict = dict(cloud_dict[cloud]['cores'][core])
            center_wcs = np.copy(core_dict['center_wcs'])
            core_dict['center_wcs'] = {}
            eq_coords = hrs2degs(ra=center_wcs[0],
                                 dec=center_wcs[1])

            core_dict['center_wcs']['eq'] = eq_coords
            core_dict['center_wcs']['gal'] = \
                    switch_coords(eq_coords[0], eq_coords[1])

            cloud_dict[cloud]['cores'][core] = core_dict

    return cloud_dict

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

def convert_limit_coordinates(clouds):

    for cloud in clouds:
        clouds[cloud].update({'limit_pixels': []})

        header = clouds[cloud]['av_header']

        limit_wcs = clouds[cloud]['limit_wcs']

        for limits in limit_wcs:
            # convert centers to pixel coords
            limit_pixels = get_pix_coords(ra=limits[0],
                                         dec=limits[1],
                                         header=header)[:2].tolist()

            clouds[cloud]['limit_pixels'].append(limit_pixels)

    return clouds

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

def switch_coords(x_coords, y_coords, coord_type='equatorial'):

    ''' Switches coordinates between equatorial and galactic.

    Parameters
    ----------
    x_coords, y_coords : array-like
        N-dimensional x and y coordinates
    coord_type : str
        Coordinate system of x_coords and y_coords. Options are 'equatorial'
        and 'galactic'. Default is to switch between 'equatorial' to 'galactic'.

    Returns
    -------
    x_coords_sw, y_coords_sw : array-like
        N-dimensional x and y coordinates in the switched coordinate system.

    '''

    from astropy.coordinates import ICRS as eq
    from astropy.coordinates import Galactic as gal
    from astropy import units

    # Convert coordinates to arrays
    x_coords, y_coords = np.copy(x_coords), np.copy(y_coords)

    if coord_type.lower() == 'galactic':
        coords = gal(l=x_coords,
                  b=y_coords,
                  unit=(units.degree, units.degree)
                  )
        x_coords_sw = coords.icrs.ra.deg
        y_coords_sw = coords.icrs.dec.deg
    elif coord_type.lower() == 'equatorial':
        coords = eq(ra=x_coords,
                  dec=y_coords,
                  unit=(units.degree, units.degree)
                  )
        x_coords_sw = coords.galactic.l.deg
        y_coords_sw = coords.galactic.b.deg

    return x_coords_sw, y_coords_sw

'''
The main script
'''

def main():

    import numpy as np
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json
    #from astropy.io import fits
    import pyfits as fits

    # define directory locations
    figure_dir = '/d/bip3/ezbc/multicloud/figures/models/'
    output_dir = '/d/bip3/ezbc/multicloud/data/python_output/'

    cloud_dict = {'taurus' : {},
                  'perseus' : {},
                  'california' : {},
                  }

    # load Planck Av and GALFA HI images, on same grid
    for cloud in cloud_dict:
        file_dir = '/d/bip3/ezbc/{0:s}/data/av/'.format(cloud)
        av_data, av_header = fits.getdata(file_dir + \
                    '{0:s}_av_planck_5arcmin.fits'.format(cloud),
                    header=True)
        av_error_data, av_error_header = fits.getdata(file_dir + \
                '{0:s}_av_error_planck_5arcmin.fits'.format(cloud),
                header=True)

        hi_dir = '/d/bip3/ezbc/{0:s}/data/hi/'.format(cloud)
        hi_data, hi_header = fits.getdata(hi_dir + \
                '{0:s}_hi_galfa_cube_regrid_planckres.fits'.format(cloud),
                header=True)

        cloud_dict[cloud]['av_data'] = av_data
        cloud_dict[cloud]['av_header'] = av_header
        cloud_dict[cloud]['av_error_data'] = av_error_data
        cloud_dict[cloud]['av_error_header'] = av_error_header
        cloud_dict[cloud]['hi_data'] = hi_data
        cloud_dict[cloud]['hi_header'] = hi_header

        # define core properties
        with open('/d/bip3/ezbc/{0:s}/data/python_output/'.format(cloud) + \
                  'core_properties/{0:s}_core_properties.txt'.format(cloud),
                  'r') as f:
             cores = json.load(f)

        cores = convert_core_coordinates(cores, av_header)

        cores = load_ds9_region(cores,
                filename_base = '/d/bip3/ezbc/{0:s}/data/'.format(cloud) + \
                                'python_output/ds9_regions/taurus_av_boxes_',
                header = av_header)

        cloud_dict[cloud]['cores'] = cores

    cloud_dict['taurus']['limit_wcs'] = (((4, 51, 0), (21, 54, 0)),
                                         ((4, 5, 0), (30, 25, 0)))
    cloud_dict['perseus']['limit_wcs'] = (((3, 58, 0), (27, 6, 0)),
                                          ((3, 20, 0), (35, 0, 0)))
    cloud_dict['california']['limit_wcs'] = (((4, 49, 0), (31, 54, 0)),
                                             ((3, 40, 0), (42, 24, 0)))
    cloud_dict['taurus']['plot_cores'] = ('L1495A', 'L1495')
    cloud_dict['perseus']['plot_cores'] = ('NGC1333', 'B5')
    cloud_dict['california']['plot_cores'] = ('L1482', 'L1456')

    cloud_dict = convert_limit_coordinates(cloud_dict)

    # Calculate galactic coordinates
    cloud_dict = calc_core_coords(cloud_dict)

    # Write out cloud_dict
    #write_fitting_results(cloud_dict, output_dir + 'k09_fitting_results.txt')

    # Plot
    figure_types = ['png',]# 'pdf']
    for figure_type in figure_types:
        if 0:
            plot_phi_cnm_vs_glat(cloud_dict,
                                 savedir=figure_dir,
                                 filename='multicloud_phi_cnm_vs_glat' + \
                                    '.{0:s}'.format(figure_type),
                                 show=0,
                                 #limits=[5, 45, -6, -22]
                                 )

        plot_T_cnm_vs_glat(cloud_dict,
                           savedir=figure_dir,
                           filename='multicloud_T_cnm_vs_glat' + \
                              '.{0:s}'.format(figure_type),
                           show=0,
                           #limits=[5, 45, -6, -22]
                           )

if __name__ == '__main__':
    main()



