#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the taurus molecular cloud.
'''

import pyfits as pf
import numpy as np

''' Plotting Functions
'''

def plot_av_image(av_image=None, header=None, cores=None, title=None,
        limits=None, boxes=False, savedir='./', filename=None, show=True):

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

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    font_scale = 9
    line_weight = 600
    font_weight = 600
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'axes.weight': line_weight,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'xtick.weight': line_weight,
              'ytick.labelsize': font_scale,
              'ytick.weight': line_weight,
              'font.weight': font_weight,
              'axes.labelweight': font_weight,
              'text.usetex': True,
              #'font.family': 'sans-serif',
              'figure.figsize': (3.3, 3.3),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
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
                 cbar_size='5%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # create axes
    ax = imagegrid[0]

    # show the image
    im = ax.imshow(av_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=0,
            vmax=15,
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
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

    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    for core in cores:
        pix_coords = cores[core]['center_pixel']

        anno_color = (0.3, 0.5, 1)

        ax.scatter(pix_coords[0],pix_coords[1],
                color='w',
                s=100,
                marker='+',
                linewidths=1.5)

        ax.annotate(core,
                xy=[pix_coords[0], pix_coords[1]],
                xytext=(5,5),
                textcoords='offset points',
                fontsize=font_scale,
                color='w')

        if boxes:
            vertices = np.copy(cores[core]['poly_verts']['pixel'])
            #[:, ::-1]
            rect = ax.add_patch(Polygon(
                    vertices[:, ::-1],
                    facecolor='none',
                    edgecolor='w'))

    if title is not None:
        fig.suptitle(title, fontsize=fontScale)
    if filename is not None:
        #plt.tight_layout()
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        fig.show()

''' Calculations
'''

def create_wedge(core_pos, radius, angle, core_rel_pos=0.1):

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
    wedge_vertices : numpy.array
        4 x 2 array with box pixel vertices, starting from lower-left going
        clockwise.

    '''

    from matplotlib.patches import Circle, Wedge, Polygon

    center_pos = (core_pos[0] - core_rel_pos * radius,
                  core_pos[1])

    wedge_vertices = Wedge(center_pos, radius, -angle/2., angle/2.).get_verts()

    return wedge_vertices

def rotate_wedge(wedge_vertices, anchor, angle):

    '''
    Parameters
    ----------
    wedge_vertices : numpy.array
        4 x 2 array with box pixel vertices, starting from lower-left going
        clockwise.
    anchor : tuple
        x and y coordinates of pivot point.
    angle : float
        Angle to rotate polygon clockwise from North.

    Returns
    -------
    wedge_vertices_rotated : numpy.array
        4 x 2 array with rotated box pixel vertices.
    '''

    from mygeometry import rotate_polygon

    wedge_vertices_rotated = rotate_polygon(wedge_vertices, anchor, angle)

    return wedge_vertices_rotated

def derive_ideal_wedge(av_image, cores_dict, wedge_angle, wedge_radius,
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
    wedge_dict = {}

    for core in cores_dict:
        print('Calculating optimal angle for core {:s}'.format(core))

        # axes are reversed
        core_pos = cores_dict[core]['center_pixel'][::-1]

        wedge_vertices = create_wedge(core_pos, wedge_radius, wedge_angle,
                core_rel_pos=core_rel_pos)

        gradient_sums = np.zeros((len(angle_grid)))

        for i, angle in enumerate(angle_grid):
            wedge_vertices_rotated = rotate_wedge(wedge_vertices, core_pos, angle)

            mask = myg.get_polygon_mask(av_image, wedge_vertices_rotated)

            av_image_masked = np.copy(av_image)

            # extract radial profile weighted by SNR
            radii, profile = get_radial_profile(av_image, binsize=3,
                    center=core_pos[::-1],
                    weights=av_image_error,
                    mask=mask
                    )

            if angle == 90:
                av_image_masked = np.copy(av_image)
                mask = myg.get_polygon_mask(av_image_masked, wedge_vertices)
                av_image_masked[mask==0]=np.NaN

            indices = np.where((radii == radii) & \
                               (profile == profile))
            profile, radii = profile[indices], radii[indices]

            # steeper gradients will have smaller sums
            gradient_sum = np.sum(np.gradient(profile, radii))
            gradient_sums[i] = gradient_sum

        # find steepest profile and recreate the box mask
        angle_ideal = angle_grid[gradient_sums == np.min(gradient_sums)][0]

        wedge_vertices_rotated = rotate_wedge(wedge_vertices, core_pos, angle_ideal)

        wedge_dict[core] = {}
        wedge_dict[core]['wedge_vertices_rotated'] = wedge_vertices_rotated

    return wedge_dict

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

    return region

def load_ds9_region(cores, filename_base = 'taurus_av_boxes_', header=None):

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    regions = read_ds9_region(filename_base + '.reg')

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        core = tag[tag.find('{')+1:tag.find('}')]

        # Format vertices to be 2 x N array
        poly_verts = []
        for i in xrange(0, len(region.coord_list)/2):
            poly_verts.append((region.coord_list[2*i],
                               region.coord_list[2*i+1]))

        poly_verts_pix = []
        for i in xrange(0, len(poly_verts)):
            poly_verts_pix.append(get_pix_coords(ra=poly_verts[i][0],
                                            dec=poly_verts[i][1],
                                            header=header)[:-1][::-1].tolist())

        cores[core]['poly_verts'] = {}
        cores[core]['poly_verts']['wcs'] = poly_verts
        cores[core]['poly_verts']['pixel'] = poly_verts_pix

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
    # wedge should be a few tens of pc.
    # D = 300 pc
    # res = 5'
    # d/pix = 0.43 pc/pix
    wedge_angle = 40.0 # degrees
    wedge_radius = 10.0 / 0.43 # pixels,
    core_rel_pos = 0.15 # fraction of radius core is within wedge

    # Name of property files
    global_property_file = 'taurus_global_properties.txt'

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/maps/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/hi/'
    co_dir = '/d/bip3/ezbc/taurus/data/co/'
    core_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'
    property_dir = '/d/bip3/ezbc/taurus/data/python_output/'

    # load Planck Av and GALFA HI images, on same grid
    av_data, av_header = load_fits(av_dir + \
                'taurus_av_planck_5arcmin.fits',
            return_header=True)

    av_error_data, av_error_header = load_fits(av_dir + \
                'taurus_av_error_planck_5arcmin.fits',
            return_header=True)

    # av_data[dec, ra], axes are switched

    # define core properties
    with open(core_dir + 'taurus_core_properties.txt', 'r') as f:
        cores = json.load(f)

    cores = convert_core_coordinates(cores, av_header)

    cores = load_ds9_region(cores,
                            filename_base = region_dir + 'taurus_av_poly_cores',
                            header = av_header)

    # Open core properties
    with open(core_dir + 'taurus_core_properties.txt', 'w') as f:
        json.dump(cores, f)

    # Open file with WCS region limits
    with open(property_dir + global_property_file, 'r') as f:
        global_props = json.load(f)

    global_props = convert_limit_coordinates(global_props, header=av_header)

    # Plot
    figure_types = ['pdf', 'png']
    for figure_type in figure_types:
        plot_av_image(av_image=av_data,
                      header=av_header,
                      boxes=True,
                      cores=cores,
                      savedir=figure_dir,
                      limits=global_props['region_limit']['pixel'],
                      filename='taurus_av_cores_map.' + figure_type,
                      show=0)

if __name__ == '__main__':
    main()



