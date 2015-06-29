#!/usr/bin/python

''' Calculates the N(HI) map for multicloud

'''

import numpy as np
import matplotlib
matplotlib.use('Agg')


import pyfits as pf
import numpy as np
import warnings
warnings.filterwarnings('ignore')


''' Plotting Functions
'''

def plot_av_images(av_image_2mass=None, header=None, contour_image=None,
        av_image_planck=None, cores=None, props=None, regions=None, title=None,
        limits=None, contours=None, boxes=False, filename=None,
        show=True, hi_vlimits=None, av_vlimits=None, av_diff_vlimits=None,):

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
    plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    params = {
              'figure.figsize': (7.3, 9),
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(3,1),
                 ngrids=3,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0.3,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # ------------------
    # 2MASS image
    # ------------------
    # create axes
    ax = imagegrid[0]

    # show the image
    im = ax.imshow(av_image_2mass,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=av_vlimits[0],
            vmax=av_vlimits[1],
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)
    ax.set_title('2MASS $A_V$', fontsize=9)

    ax.locator_params(nbins=6)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')

    # Write label to colorbar
    cb.set_label_text(r'$A_V$ [mag]',)

    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    if type(cores) is dict:
        for core in cores:
            pix_coords = cores[core]['center_pixel']

            anno_color = (0.3, 0.5, 1)

            ax.scatter(pix_coords[0],pix_coords[1],
                    color=anno_color,
                    s=200,
                    marker='+',
                    linewidths=2)

            ax.annotate(core,
                    xy=[pix_coords[0], pix_coords[1]],
                    xytext=(5,5),
                    textcoords='offset points',
                    color=anno_color)

            if boxes:
                rect = ax.add_patch(Polygon(
                    cores[core]['box_vertices'][:, ::-1],
                        facecolor='none',
                        edgecolor=anno_color))

    if regions is not None:
        for region in regions:
            vertices = np.copy(regions[region]['poly_verts']['pixel'])
            rect = ax.add_patch(Polygon(
                    vertices[:, ::-1],
                    facecolor='none',
                    edgecolor='w'))

    # ------------------
    # Planck image
    # ------------------
    # create axes
    ax = imagegrid[1]
    # show the image
    im = ax.imshow(av_image_planck,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=av_vlimits[0],
            vmax=av_vlimits[1],
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)
    ax.set_title('Planck $A_V$', fontsize=9)

    ax.locator_params(nbins=6)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])


    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Write label to colorbar
    cb.set_label_text(r'$A_V$ [mag]',)

    # ------------------
    # Division image
    # ------------------
    # create axes
    ax = imagegrid[2]

    # show the image
    im = ax.imshow(av_image_planck / av_image_2mass,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=av_diff_vlimits[0],
            vmax=av_diff_vlimits[1],
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)
    ax.set_title('Planck $A_V$ / 2MASS $A_V$', fontsize=9)

    ax.locator_params(nbins=6)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Write label to colorbar
    #cb.set_label_text(r'$A_V$ [mag]',)


    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    if type(cores) is dict:
        for core in cores:
            pix_coords = cores[core]['center_pixel']

            anno_color = (0.3, 0.5, 1)

            if boxes:
                rect = ax.add_patch(Polygon(
                    cores[core]['box_vertices'][:, ::-1],
                        facecolor='none',
                        edgecolor='w'))

    if regions is not None:
        for region in regions:
            vertices = np.copy(regions[region]['poly_verts']['pixel'])
            rect = ax.add_patch(Polygon(
                    vertices[:, ::-1],
                    facecolor='none',
                    edgecolor='w'))

        if props is not None:
            for region in props['region_name_pos']:
                if region == 'taurus1':
                    region = 'taurus 1'
                if region == 'taurus2':
                    region = 'taurus 2'
                ax.annotate(region.capitalize(),
                            xy=props['region_name_pos'][region]['pixel'],
                            xytext=(0,0),
                            textcoords='offset points',
                            color='w',
                            fontsize=font_scale*7.0/9.0,
                            zorder=10)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()

def plot_av(av_x, av_y,
        av_x_error=None, av_y_error=None, limits=None,
        savedir='./', filename=None, show=True, scale='linear',
        axis_labels=(None, None),
        returnimage=False, hess_binsize=None, title='', plot_type='scatter',
        color_scale = 'linear', errorbar = None, errorbar_pos = None):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib import cm

    # Drop the NaNs from the images
    indices = np.where((av_x == av_x) &\
                       (av_y == av_y) &\
                       (av_x > 0) &\
                       (av_y > 0))

    try:
        av_x_nonans = av_x[indices]
        av_y_nonans = av_y[indices]

        if type(av_y_error) is float:
            av_y_error_nonans = sd_image_error * \
                    np.ones(av_y[indices].shape)
        else:
            av_y_error_nonans = sd_image_error[indices]

        if type(av_x_error) is np.ndarray:
            av_x_error_nonans = av_x_error[indices]
        else:
            av_x_error_nonans = av_x_error * \
                    np.ones(av_x[indices].shape)
    except NameError:
        no_errors = True

    # Create figure
    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    fontScale = 12
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
              'figure.figsize': (8, 8),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)


    # Create figure
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    if av_y_error is None:
        if plot_type is 'hexbin':
            if color_scale == 'linear':
                image = ax.hexbin(av_x_nonans.ravel(),
                        av_y_nonans.ravel(),
                        mincnt=1,
                        yscale=scale,
                        xscale=scale)
            elif color_scale == 'log':
                image = ax.hexbin(av_x_nonans.ravel(),
                    av_y_nonans.ravel(),
                    norm=matplotlib.colors.LogNorm(),
                    mincnt=1,
                    yscale=scale,
                    xscale=scale,
                    cmap = cm.gist_stern)
                for i, pos in enumerate(errorbar_pos):
                	error_bar = ax.errorbar(pos[0], pos[1],
                            xerr = (errorbar[0]),
                            yerr = (errorbar[1]),
                            color = 'c')
                ax.set_xscale(scale, nonposx = 'clip')
                ax.set_yscale(scale, nonposy = 'clip')
            # Adjust color bar of density plot
            cb = plt.colorbar(image)
            cb.set_label('Bin Counts')
        elif plot_type is 'scatter':
            image = ax.scatter(av_x_nonans.ravel(),
                    av_y_nonans.ravel(),
                    alpha=0.3,
                    color='k'
                    )
            ax.set_xscale(scale)
            ax.set_yscale(scale)
    else:
        image = ax.errorbar(av_x_nonans.ravel(),
                av_y_nonans.ravel(),
                xerr=(av_x_error_nonans.ravel()),
                yerr=(av_y_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=2
                )

        ax.set_xscale(scale)
        ax.set_yscale(scale)

    if limits is not None:
    	ax.set_xlim(limits[0],limits[1])
    	ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    #ax.set_xlabel(r'Planck A$_{\rm V}$ (mag)')
    #ax.set_ylabel(r'Lee+12 A$_{\rm V}$ (mag)')
    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    ax.set_title(title)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

''' DS9 Region and Coordinate Functions
'''

def convert_limit_coordinates(prop_dict,
        coords=('region_limit', 'co_noise_limits', 'plot_limit'), header=None):

    # Initialize pixel keys
    for coord in coords:

        if coord == 'region_limit' or coord == 'plot_limit':
            prop_dict[coord].update({'pixel': []})
            limit_wcs = prop_dict[coord]['wcs']

            for limits in limit_wcs:
                # convert centers to pixel coords
                limit_pixels = get_pix_coords(ra=limits[0],
                                             dec=limits[1],
                                             header=header)[:2].tolist()

                prop_dict[coord]['pixel'].append(limit_pixels[0])
                prop_dict[coord]['pixel'].append(limit_pixels[1])
        elif coord == 'co_noise_limits':
            prop_dict[coord].update({'pixel': []})
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
        elif coord == 'region_name_pos':

            # convert centers to pixel coords
            for region in prop_dict[coord]:
                prop_dict[coord][region].update({'pixel': []})

                coord_wcs = prop_dict[coord][region]['wcs']

                coord_pixel = get_pix_coords(ra=coord_wcs[0],
                                             dec=coord_wcs[1],
                                             header=header)[:2].tolist()

                prop_dict[coord][region]['pixel'].append(coord_pixel[0])
                prop_dict[coord][region]['pixel'].append(coord_pixel[1])

    return prop_dict

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
    from astropy.io import fits

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

def load_ds9_region(props, filename=None, header=None):

    import pyregion as pyr

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    props['regions'] = {}


    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        region_name = tag[tag.find('text={')+6:tag.find('}')].lower()

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

        props['regions'][region_name] = {}
        props['regions'][region_name]['poly_verts'] = {}
        props['regions'][region_name]['poly_verts']['wcs'] = poly_verts
        props['regions'][region_name]['poly_verts']['pixel'] = poly_verts_pix

    return props

'''
The main script
'''

def main(dgr=None, vel_range=None, vel_range_type='single', region=None,
        av_data_type='planck'):
    ''' Executes script.

    Parameters
    ----------
    dgr : float
        If None, pulls best-fit value from properties.
    vel_range : tuple
        If None, pulls best-fit value from properties.
    '''

    # import external modules
    #import pyfits as fits
    from astropy.io import fits
    import numpy as np
    from mycoords import make_velocity_axis
    import mygeometry as myg
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    import json
    from os import system,path

    # Script parameters
    # -----------------
    # Name of noise cube
    noise_cube_filename = 'multicloud_hi_galfa_cube_regrid_planckres_noise.fits'

    # Use Planck dust Av map or Kainulainen 2009 optical extinction Av map?
    # options are 'planck' or 'lee12'
    #av_data_type = 'lee12'
    #av_data_type = 'planck'

    # Global parameter file
    prop_file = 'multicloud_global_properties'

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
    output_dir = '/d/bip3/ezbc/multicloud/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/multicloud/figures/maps/'
    av_dir = '/d/bip3/ezbc/multicloud/data/av/'
    hi_dir = '/d/bip3/ezbc/multicloud/data/hi/'
    co_dir = '/d/bip3/ezbc/multicloud/data/co/'
    core_dir = '/d/bip3/ezbc/multicloud/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'

    # load Planck Av and GALFA HI images, on same grid
    print('\nLoading Planck data...')
    av_image_planck, av_header_planck = load_fits(av_dir + \
                'multicloud_av_planck_5arcmin.fits',
            return_header=True)

    av_image_2mass, av_header_2mass = load_fits(av_dir + \
                'multicloud_av_k09_nan_regrid_planckres.fits',
            return_header=True)

    # Prepare data products
    # ---------------------
    # Load global properties of cloud
    # global properties written from script
    # 'av/multicloud_analysis_global_properties.txt'
    if region is not None:
        likelihood_filename += '_region{0:.0f}'.format(region)
        results_filename += '_region{0:.0f}'.format(region)
    with open(property_dir + prop_file + '.txt', 'r') as f:
        props = json.load(f)

    # Change WCS coords to pixel coords of images
    props = convert_limit_coordinates(props,
                                      header=av_header_planck,
                                      coords=('region_limit',
                                              'co_noise_limits',
                                              'plot_limit',
                                              'region_name_pos'))

    # Plot
    figure_types = ['png', 'pdf']
    for figure_type in figure_types:

        filename = figure_dir + 'multicloud_av_2mass_planck_compare_map.' + \
            figure_type

        print('\nSaving Av model image to \n' + filename)

        plot_av_images(av_image_2mass=av_image_2mass,
                       av_image_planck=av_image_planck,
                       header=av_header_planck,
                       props=props,
                       av_vlimits=(-0.1,16),
                       av_diff_vlimits=(0.3,2.1),
                       limits=[50,50, 500, 370],
                       filename=filename,
                       show=False)

        filename = figure_dir + 'multicloud_av_2mass_planck_plot.' + \
            figure_type

        print('\nSaving Av vs Av to \n' + filename)

        plot_av(av_image_planck, av_image_2mass,
                plot_type = 'hexbin',
                color_scale = 'log',
                scale = 'log',
                limits = (0.2, 30, 0.2, 30),
                axis_labels=(r'Planck $A_V$ [mag]', '2MASS $A_V$ [mag]'),
                errorbar = (0.1, 0.3),
                errorbar_pos = ((2., 2.),
                                (0.5, 0.5),
                                (9., 9.)),
                savedir = '',
                filename = filename,
                show = False,
                )

if __name__ == '__main__':

    main()


