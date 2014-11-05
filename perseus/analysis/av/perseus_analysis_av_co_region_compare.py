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
                       linestyles='--')
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

def convert_limit_coordinates(noise_dict):

    noise_dict.update({'limit_pixels': []})

    header = noise_dict['av_header']

    limit_wcs = noise_dict['limit_wcs']

    for limits in limit_wcs:
        # convert centers to pixel coords
        limit_pixels = get_pix_coords(ra=limits[0],
                                     dec=limits[1],
                                     header=header)[:2].tolist()

        noise_dict['limit_pixels'].append(limit_pixels[0])
        noise_dict['limit_pixels'].append(limit_pixels[1])

    return noise_dict

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
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json

    # parameters used in script
    # -------------------------
    # WCS limits of boxed region without emission
    noise_dict = {}
    noise_region = noise_dict['limit_wcs'] = (((3, 34, 40), (33, 3, 0)),
                                             ((3, 21, 0), (36, 38, 0)))
    prop_dict = {}
    prop_dict['limit_wcs'] = (((3, 58, 0), (27, 6, 0)),
                              ((3, 20, 0), (35, 0, 0)))
    prop_dict['limit_wcs'] = (((3, 58, 0), (26, 6, 0)),
                              ((3, 0, 0), (35, 0, 0)))

    # Pixel value of masks
    global bad_pix
    bad_pix = -1e8

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/perseus/figures/maps/'
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
    noise_dict['av_header'] = av_header
    prop_dict['av_header'] = av_header
    prop_dict = convert_limit_coordinates(prop_dict)

    av_error_data, av_error_header = load_fits(av_dir + \
                'perseus_av_error_planck_5arcmin.fits',
            return_header=True)

    co_data, co_header = load_fits(co_dir + \
                'perseus_co_cfa_cube_regrid_planckres.fits',
            return_header=True)

    # Create moment 0 map of CO
    co_mom0 = np.sum(co_data, axis=0)

    # Extract pixels where only CO noise is present
    noise_dict = convert_limit_coordinates(noise_dict)
    limits = noise_dict['limit_pixels']
    co_mom0_noise_region = co_mom0[limits[0]:limits[2], limits[1]:limits[3]]
    std = np.std(co_mom0_noise_region[~np.isnan(co_mom0_noise_region)])

    # Extract locations of diffuse gas traced by lack of CO emission
    co_indices = ((co_mom0 < std * 2.0) & \
                  (av_data == av_data))
    av_image_co_masked = np.copy(av_data)
    av_image_co_masked[co_indices] = bad_pix

    # Extract locations of diffuse gas traced by low Av. Raise Av threshold
    # until the same number of solely atomic gas column pixels are recovered
    # for both the Av and CO priors
    num_co_pix = co_indices[co_indices].size
    num_av_pix = 0
    av_threshold = 0.0 # mag
    while num_av_pix <= num_co_pix:
    	av_threshold += 0.05
    	av_indices = ((av_data < av_threshold) & \
    	              (co_mom0 == co_mom0))
    	num_av_pix = av_indices[av_indices].size

    av_image_av_masked = np.copy(av_data)
    av_image_av_masked[av_data < av_threshold] = bad_pix

    print('\nAv threshold = {0:.1f} mag'.format(av_threshold))

    # Plot
    figure_types = ['png', 'pdf']
    for figure_type in figure_types:
        plot_av_image(av_image=av_data,
                header=av_header,
                av_mask=av_image_av_masked,
                co_mask=av_image_co_masked,
                av_threshold=av_threshold,
                savedir=figure_dir,
                limits=prop_dict['limit_pixels'],
                filename='perseus_av_co_masks_map.%s' % \
                        figure_type,
                show=0)

        plot_av_image(av_image=av_data,
                header=av_header,
                co_mask=av_image_co_masked,
                av_thresholds=av_thresholds,
                savedir=figure_dir,
                limits=prop_dict['limit_pixels'],
                filename='perseus_av_co_masks_map.%s' % \
                        figure_type,
                show=0)

    '''
    # define core properties
    with open(core_dir + 'perseus_core_properties.txt', 'r') as f:
        cores = json.load(f)

    cores = convert_core_coordinates(cores, av_header)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'perseus_av_boxes_',
            header = av_header)

    av_image_list = []
    av_image_error_list = []
    core_name_list = []

    box_dict = derive_ideal_box(av_data, cores, box_width, box_height,
            core_rel_pos=0.1, angle_res=10., av_image_error=av_error_data)

    for core in cores:
        cores[core]['box_vertices_rotated'] = \
            box_dict[core]['box_vertices_rotated'].tolist()
        try:
            cores[core]['center_pixel'] = cores[core]['center_pixel'].tolist()
        except AttributeError:
            cores[core]['center_pixel'] = cores[core]['center_pixel']

    with open(core_dir + 'perseus_core_properties.txt', 'w') as f:
        json.dump(cores, f)

    for core in cores:
        mask = myg.get_polygon_mask(av_data,
                cores[core]['box_vertices_rotated'])

        av_data_mask = np.copy(av_data)
        av_data_mask[mask == 0] = np.NaN

    # Define limits for plotting the map
    noise_dict = {}
    noise_dict['limit_wcs'] = (((3, 58, 0), (27, 6, 0)),
                              ((3, 20, 0), (35, 0, 0)))
    noise_dict['limit_wcs'] = (((3, 58, 0), (26, 6, 0)),
                              ((3, 0, 0), (35, 0, 0)))
    noise_dict['av_header'] = av_header

    noise_dict = convert_limit_coordinates(noise_dict)

    # Plot
    figure_types = ['pdf', 'png']
    for figure_type in figure_types:
        plot_av_image(av_image=av_data, header=av_header,
                boxes=True, cores=cores, #limits=[50,37,200,160],
                #title=r'perseus: A$_V$ map with core boxed-regions.',
                savedir=figure_dir,
                limits=noise_dict['limit_pixels'],
                filename='perseus_av_cores_map.%s' % \
                        figure_type,
                show=0)
    '''



if __name__ == '__main__':
    main()


