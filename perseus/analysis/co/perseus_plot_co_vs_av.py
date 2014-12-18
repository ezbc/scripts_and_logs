#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the perseus molecular cloud.
'''

import pyfits as pf
import numpy as np

''' Plotting Functions
'''

def plot_av_image(av_image=None, header=None, title=None, limits=None,
        savedir='./', filename=None, show=True, av_mask=None, co_mask=None,
        av_threshold=None, av_thresholds=None):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    import matplotlib.cm as cm
    import matplotlib.lines as mlines

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    fontScale = 15
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (8, 7),
              'figure.titlesize': fontScale
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    nrows_ncols=(1,1)
    ngrids=1

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=1,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # create axes
    ax = imagegrid[0]
    cmap = cm.pink # colormap
    cmap = cm.gray # colormap
    # show the image
    im = ax.imshow(av_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',)
    ax.set_ylabel('Declination (J2000)',)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Write label to colorbar
    cb.set_label_text(r'A$_V$ (Mag)',)

    # Show contour masks
    cs_co = ax.contour(co_mask,
                       levels=(bad_pix,),
                       origin='lower',
                       colors='r',
                       linestyles='-')
    if av_mask is not None:
        cs_av = ax.contour(av_mask,
                           levels=(bad_pix,),
                           origin='lower',
                           colors='c',
                           linestyles='solid')

    # Legend
    co_line = mlines.Line2D([], [],
                color='r',
                linestyle='--',
                label=r'CO threshold = 2$\times \sigma_{\rm CO}$')
    if av_thresholds is not None:
    	colors = cm.rainbow(np.linspace(0, 1, len(av_thresholds)))
        for i, av_threshold in enumerate(av_thresholds):
            label = r'$A_V$ threshold = {0:1f} mag'.format(av_threshold)
            av_line = mlines.Line2D([], [],
                        color=colors[i],
                        linestyle='solid',
                        label=label)
    elif av_threshold is not None and av_mask is not None:
        label = r'$A_V$ threshold = {0:1f} mag'.format(av_threshold)
        av_line = mlines.Line2D([], [],
                    color='c',
                    linestyle='solid',
                    label=label)

    #ax.clabel(cs_co, inline=1, fontsize=10)
    #ax.clabel(cs_av, inline=1, fontsize=10)

    #cs_co.collections.set_label(r'CO threshold = 2$\times \sigma_{\rm CO}$')
    #cs_av.collections.set_label(r'$A_V$ threshold = ' + \
    #                             '{0:1f} mag'.format(av_threshold))
    lines = [cs_co.collections[0], cs_av.collections[0]]
    labels = [r'CO threshold = 2$\times\ \sigma_{\rm CO}$',
              r'$A_V$ threshold = {0:.1f} mag'.format(av_threshold)]
    ax.legend(lines, labels,)
              #bbox_to_anchor=(1,1),
              #loc=3,
              #ncol=2,
              #mode="expand",
              #borderaxespad=0.)

    if title is not None:
        fig.suptitle(title, fontsize=fontScale)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        fig.show()

def plot_co_vs_av(co_images, av_images, co_error_images=None,
        av_error_images=None, limits=None, fit=True, savedir='./',
        filename=None, show=True, scale='linear', title = ''):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from myscience.krumholz09 import calc_T_cnm
    from matplotlib import cm

    n = int(np.ceil(len(av_images)**0.5))
    if n**2 - n > len(av_images):
        nrows = n - 1
        ncols = n
        y_scaling = 1.0 - 1.0 / n
    else:
        nrows, ncols = n, n
        y_scaling = 1.0

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
              'figure.figsize': (5, 5 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=len(av_images),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 cbar_mode='single',
                 cbar_pad=0.1,
                 cbar_size=0.2,
                 )

    # Cycle through lists
    for i in xrange(len(av_images)):
        av = av_images[i]
        co = co_images[i]
        av_error = av_error_images[i]
        co_error = co_error_images[i]

        # Drop the NaNs from the images
        if type(av_error) is float:
            indices = np.where((av == av) &\
                               (co == co)&\
                               (co > 0) &\
                               (av > 0))

        if type(av_error) is np.ndarray or \
                type(av_error) is np.ma.core.MaskedArray or \
                type(co_error) is np.ndarray or \
                type(co_error) is np.ma.core.MaskedArray:
            indices = np.where((av == av) &\
                               (co == co) &\
                               (co_error == co_error) &\
                               (av_error == av_error) &\
                               (co > 0) &\
                               (av > 0))

        av_nonans = av[indices]
        co_nonans = co[indices]

        # Fix error data types
        if type(av_error) is np.ndarray:
            av_error_nonans = av_error[indices]
        else:
            av_error_nonans = np.array(av_error[indices])

        if type(co_error) is np.ndarray or \
                type(co_error) is np.ma.core.MaskedArray:
            co_error_nonans = co_error[indices]
        else:
            co_error_nonans = co_error * \
                    np.ones(co[indices].shape)

        # Create plot
        ax = imagegrid[i]

        image = ax.hexbin(av_nonans.ravel(),
            co_nonans.ravel(),
            norm=matplotlib.colors.LogNorm(),
            mincnt=1,
            yscale=scale[0],
            xscale=scale[0],
            gridsize=(20,100),
            cmap = cm.Greys)
        cb = ax.cax.colorbar(image,)
        # Write label to colorbar
        cb.set_label_text('Bin Counts',)

        # Plot sensitivies
        co_limit = co_error_images[0]
        av_limit = np.median(av_error_images[0])
        print 'av_limit', av_limit
        ax.axhline(co_limit, color='k', linestyle='--')
        ax.axvline(av_limit, color='k', linestyle='--')

        '''
        # Plot with error bars
        image = ax.errorbar(co_nonans.ravel(),
                av_nonans.ravel(),
                xerr=(co_error_nonans.ravel()),
                yerr=(av_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='None',
                markersize=3
                )
        '''

        # Annotations
        anno_xpos = 0.95

        '''
        if phi_cnm_list is not None and Z_list is not None:
            if phi_cnm_error_list is None and Z_error_list is not None:
                ax.annotate(r'$\phi_{\rm CNM}$ = {0:.2f}\n'.format(phi_cnm) + \
                            r'Z = {0:.2f} Z$_\odot$'.format(Z),
                        xytext=(anno_xpos, 0.05),
                        xy=(anno_xpos, 0.05),
                        textcoords='axes fraction',
                        xycoords='axes fraction',
                        color='k',
                        bbox=dict(boxstyle='round',
                                  facecolor='w',
                                  alpha=0.3),
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        )
            else:
            	T_cnm = calc_T_cnm(phi_cnm, Z=Z)
            	T_cnm_error = []
            	T_cnm_error.append(\
            	        T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[0], Z=Z))
            	T_cnm_error.append(\
            	        T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[1], Z=Z))

                phi_cnm_text = r'\noindent$\phi_{\rm CNM}$ =' + \
                               r' %.2f' % (phi_cnm) + \
                               r'$^{+%.2f}_{-%.2f}$ \\' % (phi_cnm_error[0],
                                                         phi_cnm_error[1])
                T_cnm_text = r'\noindent T$_{\rm CNM}$ =' + \
                             r' %.2f' % (T_cnm) + \
                             r'$^{+%.2f}_{-%.2f}$ \\' % (T_cnm_error[0],
                                                         T_cnm_error[1])
                if Z_error == (0.0, 0.0):
                	Z_text = r'Z = %.1f Z$_\odot$ \\' % (Z)
                else:
                	Z_text = r'Z = %.2f' % (Z) + \
                    r'$^{+%.2f}_{-%.2f}$ Z$_\odot$ \\' % (Z_error[0],
                                                          Z_error[1])
                if phi_mol_error == (0.0, 0.0):
                    phi_mol_text = r'\noindent$\phi_{\rm mol}$ = ' + \
                                     '%.1f' % (phi_mol)
                else:
                    phi_mol_text = r'\noindent$\phi_{\rm mol}$ =' + \
                                r' %.2f' % (phi_mol) + \
                                r'$^{+%.2f}_{-%.2f}$ \\' % (phi_mol_error[0],
                                                         phi_mol_error[1])

                ax.annotate(phi_cnm_text + T_cnm_text + Z_text + phi_mol_text,
                        xytext=(anno_xpos, 0.05),
                        xy=(anno_xpos, 0.05),
                        textcoords='axes fraction',
                        xycoords='axes fraction',
                        size=font_scale*3/4.0,
                        color='k',
                        bbox=dict(boxstyle='round',
                                  facecolor='w',
                                  alpha=1),
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        )
        '''

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'$A_V$ (mag)')
        ax.set_ylabel(r'I$_{\rm CO}$ K km/s')
        #ax.set_title(core_names[i])

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
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

def derive_ideal_box(av_image, cores_dict, box_width, box_height,
        av_image_error=None, core_rel_pos=0.1, angle_res=1.0):

    import mygeometry as myg

    """
    Parameters
    ----------
    angle_res : float
        Resolution with which to rotate each new box in degrees. 1.0 degree
        gives 360 different box orientations.


    """

    angle_grid = np.arange(0, 360, angle_res)
    box_dict = {}

    for core in cores_dict:
        print('Calculating optimal angle for core {:s}'.format(core))

        # axes are reversed
        core_pos = cores_dict[core]['center_pixel'][::-1]

        box_vertices = create_box(core_pos, box_width, box_height,
                core_rel_pos=core_rel_pos)

        gradient_sums = np.zeros((len(angle_grid)))

        for i, angle in enumerate(angle_grid):
            box_vertices_rotated = rotate_box(box_vertices, core_pos, angle)

            mask = myg.get_polygon_mask(av_image, box_vertices_rotated)

            av_image_masked = np.copy(av_image)

            # extract radial profile weighted by SNR
            radii, profile = get_radial_profile(av_image, binsize=3,
                    center=core_pos[::-1],
                    weights=av_image_error,
                    mask=mask
                    )

            if angle == 90:
                av_image_masked = np.copy(av_image)
                mask = myg.get_polygon_mask(av_image_masked, box_vertices)
                av_image_masked[mask==0]=np.NaN

            indices = np.where((radii == radii) & \
                               (profile == profile))
            profile, radii = profile[indices], radii[indices]

            # steeper gradients will have smaller sums
            gradient_sum = np.sum(np.gradient(profile, radii))
            gradient_sums[i] = gradient_sum

        # find steepest profile and recreate the box mask
        angle_ideal = angle_grid[gradient_sums == np.min(gradient_sums)][0]

        box_vertices_rotated = rotate_box(box_vertices, core_pos, angle_ideal)

        box_dict[core] = {}
        box_dict[core]['box_vertices_rotated'] = box_vertices_rotated

    return box_dict

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

def calc_co_noise(co_mom0, prop_dict):

    co_noise_region = []

    # Append pixels from each region to CO region map
    for region in prop_dict['co_noise_limits']['pixel']:
        co_noise_region.append(co_mom0[region[0][1]:region[1][1],
                                       region[0][0]:region[1][0]])

    # Calc noise
    noise = 0.0
    for region in co_noise_region:
    	std = np.std(np.array(region)[~np.isnan(region)])
    	noise += std

    # Take average of stds
    noise = noise / len(co_noise_region)

    return noise

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

def convert_limit_coordinates(prop_dict,
        coords=('region_limit', 'co_noise_limits'), header=None):

    # Initialize pixel keys
    for coord in coords:
        prop_dict[coord].update({'pixel': []})

        if coord == 'region_limit':
            limit_wcs = prop_dict[coord]['wcs']

            for limits in limit_wcs:
                # convert centers to pixel coords
                limit_pixels = get_pix_coords(ra=limits[0],
                                             dec=limits[1],
                                             header=header)[:2].tolist()

                prop_dict[coord]['pixel'].append(limit_pixels[0])
                prop_dict[coord]['pixel'].append(limit_pixels[1])
        elif coord == 'co_noise_limits':
            region_limits = prop_dict[coord]['wcs']

            # Cycle through each region, convert WCS limits to pixels
            for region in region_limits:
                region_pixels = []
                for limits in region:
                    # convert centers to pixel coords
                    limit_pixels = get_pix_coords(ra=limits[0],
                                                  dec=limits[1],
                                                  header=header)[:2].tolist()
                    region_pixels.append(limit_pixels)

                # Append individual regions back to CO noise
                prop_dict[coord]['pixel'].append(region_pixels)

    return prop_dict

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

def load_ds9_region(cores, filename_base = 'perseus_av_boxes_', header=None):

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
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json

    # parameters used in script
    # -------------------------
    # Pixel value of masks
    global bad_pix
    bad_pix = -1e8

    global_property_file = 'perseus_global_properties.txt'

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/perseus/figures/co/'
    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    hi_dir = '/d/bip3/ezbc/perseus/data/hi/'
    co_dir = '/d/bip3/ezbc/perseus/data/co/'
    core_dir = '/d/bip3/ezbc/perseus/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/perseus/data/python_output/'
    region_dir = '/d/bip3/ezbc/perseus/data/python_output/ds9_regions/'
    likelihood_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'

    # load Planck Av and CfA CO images, on same grid
    av_data, av_header = load_fits(av_dir + \
                'perseus_av_planck_5arcmin.fits',
            return_header=True)

    av_error_data, av_error_header = load_fits(av_dir + \
                'perseus_av_error_planck_5arcmin.fits',
            return_header=True)

    co_data, co_header = load_fits(co_dir + \
                'perseus_co_cfa_cube_regrid_planckres.fits',
            return_header=True)

    # Load global properties
    with open(property_dir + global_property_file, 'r') as f:
        global_props = json.load(f)
    global_props = convert_limit_coordinates(global_props, header=av_header)

    # Create moment 0 map of CO
    co_mom0 = np.sum(co_data, axis=0)

    # Extract pixels where only CO noise is present
    co_noise = calc_co_noise(co_mom0, global_props)

    # Plot
    figure_types = ['png',]
    for figure_type in figure_types:
        plot_co_vs_av((co_mom0,),
                (av_data,),
                co_error_images=(co_noise,),
                av_error_images=(av_error_data,),
                limits=[10**-1, 10**1.9, 10**0, 10**2],
                #limits=[0,5,0,20],
                savedir=figure_dir,
                scale=('log', 'log'),
                #scale=('linear', 'linear'),
                filename='perseus_co_vs_av.%s' % \
                        figure_type,
                show = False,
                )

if __name__ == '__main__':
    main()



