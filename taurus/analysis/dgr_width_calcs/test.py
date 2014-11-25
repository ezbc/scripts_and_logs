#!/usr/bin/python
''' Plotting Functions
'''
def plot_likelihoods(likelihoods,velocity_centers,velocity_widths,
        filename=None,show=True, returnimage=False):
    ''' Plots a heat map of likelihoodelation values as a function of velocity width
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

    # Unravel the likelihoods if raveled
    if len(likelihoods.shape) == 1:
        likelihoods = np.empty((velocity_centers.shape[0],
                                       velocity_widths.shape[0]))
        likelihoods[:,:] = np.NaN
        count = 0
        try:
            for i, center in enumerate(velocity_centers):
                for j, width in enumerate(velocity_widths):
                    likelihoods[i,j] = likelihoods[count]
                    count += 1
        except IndexError:
            print(' plot_likelihoods: O-d array input, cannot proceed')
    else:
           likelihoods = likelihoods

    image = np.ma.array(likelihoods, mask=np.isnan(likelihoods))

    ax = imagegrid[0]

    ax.set_xlabel('Velocity Width (km/s)')
    ax.set_ylabel('Velocity Center (km/s)')

    #ax.set_xticks(np.arange(0,velocity_widths.shape[0],1)[::5],
    #        velocity_centers[::5])

    plt.rc('text', usetex=False)
    im = ax.imshow(image, interpolation='nearest', origin='lower',
            extent=[velocity_widths[0],velocity_widths[-1],
                    velocity_centers[0],velocity_centers[-1]],
            cmap=plt.cm.gist_stern,
            #cmap=plt.cm.gray,
            #norm=matplotlib.colors.LogNorm(),
            )
    cb = ax.cax.colorbar(im)

    #cb.set_clim(vmin=0.)
    # Write label to colorbar
    cb.set_label_text(r'log L')

    fractions = np.array([0.95, 0.68])
    levels = (1 + fractions * image.min())

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
        plt.savefig(filename,bbox_inches='tight')
    if show:
        plt.draw()
        plt.show()
    if returnimage:
        return likelihoods

def plot_likelihoods_hist(global_props, filename=None, show=True,
        returnimage=False, plot_axes=('centers', 'widths'),
        contour_confs=None):

    ''' Plots a heat map of likelihoodelation values as a function of velocity width
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
              'text.usetex': False,
              #'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    fig, ax_image = plt.subplots(figsize=(6,6))

    if plot_axes[0] == 'centers':
        x_extent = x_grid[0], x_grid[-1]
        ax_image.set_xlabel(r'Velocity Center (km/s)')
        x_sum_axes = (1, 2)
        y_pdf_label = r'Centers PDF'
    if plot_axes[1] == 'centers':
        y_extent = y_grid[0], y_grid[-1]
        ax_image.set_ylabel(r'Velocity Center (km/s)')
        y_sum_axes = (1, 2)
        x_pdf_label = r'Centers PDF'
    if plot_axes[0] == 'widths':
    	x_grid = global_props['vel_widths']
    	x_confint = global_props['width_confint']
        x_extent = x_grid[0], x_grid[-1]
        ax_image.set_xlabel(r'Velocity Width (km/s)')
        x_sum_axes = (0, 2)
        y_pdf_label = r'Width PDF'
        x_limits = (x_grid[0], x_grid[-1])
    if plot_axes[1] == 'widths':
        y_extent = y_grid[0], y_grid[-1]
        ax_image.set_ylabel(r'Velocity Width (km/s)')
        y_sum_axes = (0, 2)
        x_pdf_label = r'Width PDF'
    if plot_axes[0] == 'dgrs':
        x_extent = x_grid[0], x_grid[-1]
        ax_image.set_xlabel(r'DGR (10$^{-20}$ cm$^2$ mag$^1$)')
        x_sum_axes = (0, 1)
        y_pdf_label = r'DGR PDF'
    if plot_axes[1] == 'dgrs':
    	y_grid = global_props['dgrs']
    	y_confint = global_props['dgr_confint']
        y_extent = y_grid[0], y_grid[-1]
        ax_image.set_ylabel(r'DGR (10$^{-20}$ cm$^2$ mag$^1$)')
        y_sum_axes = (0, 1)
        x_pdf_label = r'DGR PDF'
        y_limits = (y_grid[0], y_grid[-1])

    # Create axes
    sum_axes = np.array((x_sum_axes, y_sum_axes))
    sum_axis = np.argmax(np.bincount(np.ravel(sum_axes)))

    # Mask NaNs
    likelihoods = global_props['likelihoods']
    image = np.ma.array(likelihoods, mask=np.isnan(likelihoods))

    # Create likelihood image
    image = np.sum(likelihoods, axis=sum_axis) / np.sum(likelihoods)

    # Derive marginal distributions of both centers and widths
    x_sum = np.sum(likelihoods, axis=x_sum_axes)
    x_pdf = x_sum / np.sum(x_sum)
    y_sum = np.sum(likelihoods, axis=y_sum_axes)
    y_pdf = y_sum / np.sum(y_sum)

    extent = np.ravel(np.array((x_extent, y_extent)))

    #plt.rc('text', usetex=False)
    im = ax_image.imshow(image.T, interpolation='nearest', origin='lower',
            extent=extent,
            #cmap=plt.cm.gist_stern,
            #cmap=plt.cm.gray,
            cmap=plt.cm.binary,
            #norm=matplotlib.colors.LogNorm(),
            aspect='auto',
            )

    show_pdfs = 1

    if show_pdfs:
        divider = make_axes_locatable(ax_image)
        ax_pdf_x = divider.append_axes("top", 1, pad=0.1, sharex=ax_image)
        ax_pdf_y  = divider.append_axes("right", 1, pad=0.1,
                sharey=ax_image)

        # make some labels invisible
        plt.setp(ax_pdf_x.get_xticklabels() + \
                 ax_pdf_y.get_yticklabels(),
                 visible=False)

        ax_pdf_x.plot(x_grid,
                      x_pdf,
                      color='k',
                      drawstyle='steps-post',
                      linewidth=2,
                      )

        ax_pdf_y.plot(y_pdf,
                      y_grid,
                      color='k',
                      drawstyle='steps-post',
                      linewidth=2,
                      )

        #axHistx.axis["bottom"].major_ticklabels.set_visible(False)

        # Tick marks on the pdf?
        pdf_ticks = False

        for tl in ax_pdf_x.get_xticklabels():
            tl.set_visible(False)

        if pdf_ticks:
            wmax = x_pdf.max()
            ticks = [0, 0.5*wmax, 1.0*wmax]
            tick_labels = ['{0:.1f}'.format(ticks[0]),
                           '{0:.1f}'.format(ticks[1]),
                           '{0:.1f}'.format(ticks[2]),
                            ]
            ax_pdf_x.set_yticks(ticks)
            ax_pdf_x.set_yticklabels(tick_labels)
        else:
            for tl in ax_pdf_x.get_yticklabels():
                tl.set_visible(False)

        ax_pdf_x.set_ylabel(y_pdf_label)

        for tl in ax_pdf_y.get_yticklabels():
            tl.set_visible(False)
        if pdf_ticks:
            cmax = y_pdf.max()
            ticks = [0, 0.5*cmax, 1.0*cmax]
            tick_labels = ['{0:.1f}'.format(ticks[0]),
                           '{0:.1f}'.format(ticks[1]),
                           '{0:.1f}'.format(ticks[2]),
                            ]
            ax_pdf_y.set_xticks(ticks)
            ax_pdf_y.set_xticklabels(tick_labels)
        else:
            for tl in ax_pdf_y.get_xticklabels():
                tl.set_visible(False)

        ax_pdf_y.set_xlabel(x_pdf_label)

        # Show confidence limits
        if y_confint is not None:
            ax_pdf_y.axhspan(y_confint[0] - y_confint[1],
                             y_confint[0] + y_confint[2],
                             color='k',
                             linewidth=1,
                             alpha=0.2)
            ax_pdf_y.axhline(y_confint[0],
                             color='k',
                             linestyle='--',
                             linewidth=3,
                             alpha=1)
        if x_confint is not None:
            ax_pdf_x.axvspan(x_confint[0] - x_confint[1],
                                 x_confint[0] + x_confint[2],
                                  color='k',
                                 linewidth=1,
                                  alpha=0.2)
            ax_pdf_x.axvline(x_confint[0],
                                 color='k',
                                 linestyle='--',
                                 linewidth=3,
                                 alpha=1)

    #cb.set_clim(vmin=0.)
    # Write label to colorbar
    #cb.set_label_text(r'log L')

    # Plot contours
    if contour_confs is not None:

        fractions = (1.0 - np.asarray(contour_confs))
        levels = (fractions * image.max())

        cs = ax_image.contour(image.T, levels=levels, origin='lower',
                extent=extent,
                colors='k'
                )

        # Define a class that forces representation of float to look a certain
        # way This remove trailing zero so '1.0' becomes '1'
        class nf(float):
             def __repr__(self):
                 str = '%.1f' % (self.__float__(),)
                 if str[-1]=='0':
                     return '%.0f' % self.__float__()
                 else:
                     return '%.1f' % self.__float__()

        # Recast levels to new class
        cs.levels = [nf(val) for val in np.asarray(contour_confs)*100.0]

        #fmt = {}
        #for level, fraction in zip(cs.levels, fractions):
        #    fmt[level] = fraction
        fmt = '%r %%'

        ax_image.clabel(cs, cs.levels, fmt=fmt, fontsize=9, inline=1)

    try:
        ax_image.set_xlim(x_limits)
        ax_image.set_ylim(y_limits)
    except UnboundLocalError:
        pass

    if 0:
    #if npix is not None or av_threshold is not None:
    	text = ''
        if npix is not None:
            text += r'N$_{\rm pix}$ = ' + \
                     '{0:.0f}'.format(npix)
            if av_threshold is not None:
            	text += '\n'
        if av_threshold is not None:
            text += r'$A_V$ threshold = {0:.1f} mag'.format(av_threshold)
            text += '\n'
        text += r'DGR = {0:.2f} '.format(y_confint[0]) + \
                r'$\times$ 10$^{-20}$ (cm$^2$ mag$^1$)'
        text += '\n'
        text += r'Velocity width = {0:.2f} '.format(x_confint[0]) + \
                r'km/s'
        ax_image.annotate(text,
                xytext=(0.95, 0.95),
                xy=(0.95, 0.95),
                textcoords='axes fraction',
                xycoords='axes fraction',
                color='k',
                fontsize=font_scale*0.75,
                bbox=dict(boxstyle='round',
                          facecolor='w',
                          alpha=0.3),
                horizontalalignment='right',
                verticalalignment='top',
                )

    if filename is not None:
        plt.draw()
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.draw()
        plt.show()
    if returnimage:
        return likelihoods

def plot_av_image(av_image=None, header=None, title=None,
        limits=None, savedir='./', filename=None, show=True):

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
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 15
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (8, 7),
              'figure.titlesize': font_scale
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
    cmap = cm.jet # colormap
    # show the image
    im = ax.imshow(av_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            #norm=matplotlib.colors.LogNorm()
            vmin=0,
            vmax=1.4
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

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
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

        if coord == 'region_limit' or coord == 'plot_limit':
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

# ==============================================================================
# Script begins
# ==============================================================================

region = 1

import grid
import numpy as np
import numpy
from os import system,path
import mygeometry as myg
from mycoords import make_velocity_axis
import json
from myimage_analysis import calculate_nhi, calculate_noise_cube, \
    calculate_sd, calculate_nh2, calculate_nh2_error
from multiprocessing import Pool

global hi_cube
global hi_velocity_axis
global hi_noise_cube
global av_image
global av_image_error

# parameters used in script
# -------------------------
# HI velocity integration range
# Determine HI integration velocity by CO or likelihoodelation with Av?
hi_av_likelihoodelation = True

center_vary = False
width_vary = True
dgr_vary = True

# Check if likelihood file already written, rewrite?
clobber = 1

# Confidence of parameter errors
conf = 0.68
# Confidence of contour levels
contour_confs = (0.68, 0.95)

# Course, large grid or fine, small grid?
grid_res = 'fine'
grid_res = 'course'

# Use Av+CO mask or only CO?
av_and_co_mask = True

# Derive CO mask? If co_thres = None, co_thres will be 2 * std(co)
co_thres = 6.00 # K km/s

# Threshold of Av below which we expect only atomic gas, in mag
av_thres = 1.4

# Center method, HI center or CO center?
vel_center_method = 'co'

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

# Results and fits filenames
if av_and_co_mask:
    likelihood_filename = 'taurus_nhi_av_likelihoods_co_' + \
                          'av{0:.1f}mag'.format(av_thres) + \
                          '_region{0:.0f}'.format(region)
    results_filename = 'taurus_likelihood_co_' + \
                       'av{0:.1f}mag'.format(av_thres) + \
                       '_region{0:.0f}'.format(region)
else:
    likelihood_filename = 'taurus_nhi_av_likelihoods_co_only'
    results_filename = 'taurus_likelihood_co_only'

# Name of property files results are written to
global_property_file = 'taurus_global_properties'
core_property_file = 'taurus_core_properties.txt'

# Name of noise cube
noise_cube_filename = 'taurus_hi_galfa_cube_regrid_planckres_noise.fits'

# Define ranges of parameters
if center_vary and width_vary and dgr_vary:
    likelihood_filename += '_width_dgr_center'
    results_filename += '_width_dgr_center'

    velocity_centers = np.arange(-15, 30, 1)
    velocity_widths = np.arange(1, 80, 1)
    dgrs = np.arange(1e-2, 1, 2e-2)
elif not center_vary and width_vary and dgr_vary:
    if grid_res == 'course':
        likelihood_filename += '_dgr_width_lowres'
        results_filename += '_dgr_width_lowres'
        velocity_centers = np.arange(5, 6, 1)
        velocity_widths = np.arange(0.1, 40, 5*0.16667)
        dgrs = np.arange(0.07, 0.4, 1e-2)
    elif grid_res == 'fine':
        likelihood_filename += '_dgr_width_highres'
        results_filename += '_dgr_width_highres'
        velocity_centers = np.arange(5, 6, 1)
        velocity_widths = np.arange(1, 100, 0.16667)
        dgrs = np.arange(0.15, 0.4, 1e-3)
        velocity_widths = np.arange(1, 15, 0.16667)
        dgrs = np.arange(0.1, 0.9, 3e-3)
        velocity_widths = np.arange(1, 60, 0.16667)
        dgrs = np.arange(0.01, 0.5, 2e-3)
        #velocity_widths = np.arange(1, 40, 1)
        #dgrs = np.arange(0.15, 0.4, 1e-1)
elif center_vary and width_vary and not dgr_vary:
    likelihood_filename += '_width_center'
    results_filename += '_width_center'

    velocity_centers = np.arange(-15, 30, 1)
    velocity_widths = np.arange(1, 80, 1)
    dgrs = np.arange(1.1e-1, 1.2e-1, 0.1e-1)
elif not center_vary and width_vary and not dgr_vary:
    likelihood_filename += '_width'
    results_filename += '_width'

    velocity_centers = np.arange(5, 6, 1)
    velocity_widths = np.arange(1, 80, 1)
    dgrs = np.arange(1.1e-1, 1.2e-1, 0.1e-1)

# define directory locations
# --------------------------
output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
figure_dir = '/d/bip3/ezbc/taurus/figures/hi_velocity_range/'
av_dir = '/d/bip3/ezbc/taurus/data/av/'
hi_dir = '/d/bip3/ezbc/taurus/data/hi/'
co_dir = '/d/bip3/ezbc/taurus/data/co/'
core_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
property_dir = '/d/bip3/ezbc/taurus/data/python_output/'
region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'
likelihood_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'

# load Planck Av and GALFA HI images, on same grid
av_data, av_header = load_fits(av_dir + \
            'taurus_av_planck_5arcmin.fits',
        return_header=True)

av_data_error, av_error_header = load_fits(av_dir + \
            'taurus_av_error_planck_5arcmin.fits',
        return_header=True)

hi_data, h = load_fits(hi_dir + \
            'taurus_hi_galfa_cube_regrid_planckres.fits',
        return_header=True)

co_data, co_header = load_fits(co_dir + \
            'taurus_co_cfa_cube_regrid_planckres.fits',
        return_header=True)

# make the velocity axis
velocity_axis = make_velocity_axis(h)
co_velocity_axis = make_velocity_axis(co_header)

# Plot NHI vs. Av for a given velocity range
if not path.isfile(hi_dir + noise_cube_filename):
    noise_cube = calculate_noise_cube(cube=hi_data,
            velocity_axis=velocity_axis,
            velocity_noise_range=[90,110], header=h, Tsys=30.,
            filename=hi_dir + noise_cube_filename)
else:
    noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
        return_header=True)

# define core properties
with open(core_dir + core_property_file, 'r') as f:
    cores = json.load(f)
with open(property_dir + global_property_file + '.txt', 'r') as f:
    global_props = json.load(f)

# Change region limits if region is defined
if region is not None:
    global_props['region_limit'] = region_limit

# Change WCS coords to pixel coords of images
global_props = convert_limit_coordinates(global_props,
        header=av_header,
        coords=('region_limit', 'co_noise_limits', 'plot_limit'))

# Derive relevant region
pix = global_props['region_limit']['pixel']
region_vertices = ((pix[1], pix[0]),
                   (pix[1], pix[2]),
                   (pix[3], pix[2]),
                   (pix[3], pix[0])
                   )

# block offregion
region_mask = myg.get_polygon_mask(av_data, region_vertices)

print('\nRegion size = ' + \
      '{0:.0f} pix'.format(region_mask[region_mask == 1].size))

# Set velocity center as CO peak
co_data_nonans = np.copy(co_data)
co_data_nonans[np.isnan(co_data_nonans)] = 0.0
co_data_nonans[:,region_mask == 1] = 0
if vel_center_method.lower() == 'co':
    co_spectrum = np.sum(co_data_nonans, axis=(1,2))
    co_avg_vel = np.average(co_velocity_axis, weights=co_spectrum)
    co_peak_vel = co_velocity_axis[co_spectrum == np.max(co_spectrum)]
    #velocity_centers = np.arange(co_peak_vel, co_peak_vel + 1, 1)
    velocity_centers = np.arange(co_avg_vel, co_avg_vel + 1, 1)
    vel_center_image = None
elif vel_center_method.lower() == 'hi':
    vel_center_image = np.zeros(hi_cube.shape[1:])
    for i in xrange(0, hi_cube.shape[0]):
        for j in xrange(0, hi_cube.shape[1]):
            hi_spectrum = hi_cube[:, i, j]
            hi_spectrum[np.isnan(hi_spectrum)] = 0.0
            vel_center_image[i,j] = \
                    velocity_axis[hi_spectrum == np.max(hi_spectrum)]


print('\nVelocity center from CO = ' +\
        '{0:.2f} km/s'.format(velocity_centers[0]))

# Create mask where CO is present
core_mask = np.zeros(av_data.shape)
#for core in cores:
#    # Grab the mask
#    core_mask += myg.get_polygon_mask(av_data,
#            cores[core]['box_vertices_rotated'])

# Calc moment 0 map of CO
co_mom0 = np.sum(co_data_nonans, axis=0)

# calc noise without any emission if CO threshold not already set
if co_thres is None:
    co_noise = calc_co_noise(co_mom0, global_props)
    co_thres = 2.0 * co_noise

# Get indices which trace only atomic gas, i.e., no CO emission
if av_and_co_mask:
    indices = (((co_mom0 < co_thres) & \
                (av_data < av_thres)) & \
                (region_mask == 1))
elif not av_and_co_mask:
    indices = ((co_mom0 < co_thres) & \
               (region_mask == 1))
    av_thres = None

# Write mask of pixels not used
mask = ~indices

# Mask global data with CO indices
hi_data_sub = np.copy(hi_data[:, indices])
noise_cube_sub = np.copy(noise_cube[:, indices])
av_data_sub = np.copy(av_data[indices])
av_error_data_sub = np.copy(av_data_error[indices])

# import matplotlib.pyplot as plt
# av_plot_data = np.copy(av_data)
# av_plot_data[~indices] = np.nan
# plt.imshow(av_plot_data, origin='lower')
# plt.contour(co_mom0, levels=(6, 12, 24), origin='lower')
# plt.show()
# plt.clf()
# plt.close()

# Plot the masked image
av_data_masked = np.copy(av_data)
av_data_masked[~indices] = np.nan
figure_types = ['png',]
for figure_type in figure_types:
    plot_av_image(av_image=av_data_masked, header=av_header,
            savedir=figure_dir + '../maps/',
            limits=global_props['plot_limit']['pixel'],
            filename='taurus_co_av_masked_map_' + \
                     'region{0:.0f}'.format(region) + \
                     figure_type,
            show=0)

# Set global variables
hi_cube = hi_data_sub
hi_velocity_axis = velocity_axis
hi_noise_cube = noise_cube_sub
av_image = av_data_sub
av_image_error = av_error_data_sub

# Define filename for plotting results
results_filename = figure_dir + results_filename

# ==============================================================================
# Populate likelihood space
# ==============================================================================

from myimage_analysis import calculate_nhi

velocity_center = 6.22
velocity_widths = np.linspace(14, 17, 10)
dgrs = np.linspace(0.08, 0.09, 20)
f = 0

#logLs = -np.inf * np.ones(len(velocity_widths) * len(dgrs) * len(fs))
logLs = -np.inf * np.ones(len(velocity_widths) * len(dgrs))

count = 0
for velocity_width in velocity_widths:
    for dgr in dgrs:
        velocity_range = (velocity_center - velocity_width / 2.,
                          velocity_center + velocity_width / 2.)

        nhi_image_temp = \
                calculate_nhi(cube=hi_cube,
                    velocity_axis=hi_velocity_axis,
                    velocity_range=velocity_range,
                    return_nhi_error=False)

        # Avoid NaNs
        indices = np.where((nhi_image_temp == nhi_image_temp) & \
                           (av_image == av_image) & \
                           (nhi_image_temp > 0))

        # Create 1D images without nans
        nhi_image_likelihood = nhi_image_temp[indices]
        av_image_likelihood = av_image[indices]

        # Create model of Av with N(HI) and DGR
        av_image_model = nhi_image_likelihood * dgr

        model = av_image_model
        data = av_image_likelihood
        error = av_image_error[indices]

        sn = np.sqrt(error**2 + f**2 * model**2)

        logL = -0.5 * np.sum((data - model)**2 / (sn**2) + \
                             np.log(2*np.pi*sn**2))

        logLs[count] = logL
        count += 1

# Normalize likelihoods
logLs -= np.max(logLs)
Ls = np.exp(logLs)
Ls /= np.sum(Ls)

# Get best-fit parameters
count = 0
max_L = 0
for velocity_width in velocity_widths:
    for dgr in dgrs:
        if Ls[count] > max_L:
            max_L = Ls[count]
            max_f = f
            max_width = velocity_width
            max_dgr = dgr
        count += 1

print ''
print 'max likelihood = ', max_L
print 'max width = ', max_width
print 'max dgr = ', max_dgr
print 'max f = ', max_f

# Sort likelihoods and print top ten best likelihoods
Ls_sorted = np.sort(Ls)
print ''
print 'Top ten likelihoods = ', Ls_sorted[-10:]

# ==============================================================================
# Calculate the error in the best-fit model
# ==============================================================================

velocity_range = (velocity_center - max_width / 2.,
                  velocity_center + max_width / 2.)

nhi_image_temp = \
        calculate_nhi(cube=hi_cube,
            velocity_axis=hi_velocity_axis,
            velocity_range=velocity_range,
            return_nhi_error=False)

# Avoid NaNs
indices = np.where((nhi_image_temp == nhi_image_temp) & \
                   (av_image == av_image) & \
                   (nhi_image_temp > 0))

# Create 1D images without nans
nhi_image_likelihood = nhi_image_temp[indices]
av_image_likelihood = av_image[indices]

# Create model of Av with N(HI) and DGR
av_image_model = nhi_image_likelihood * max_dgr

model = av_image_model
data = av_image_likelihood

print 'STD =', np.std(model - data)

error = np.std(model - data)

count = 0
for velocity_width in velocity_widths:
    for dgr in dgrs:
        velocity_range = (velocity_center - velocity_width / 2.,
                          velocity_center + velocity_width / 2.)

        nhi_image_temp = \
                calculate_nhi(cube=hi_cube,
                    velocity_axis=hi_velocity_axis,
                    velocity_range=velocity_range,
                    return_nhi_error=False)

        # Avoid NaNs
        indices = np.where((nhi_image_temp == nhi_image_temp) & \
                           (av_image == av_image) & \
                           (nhi_image_temp > 0))

        # Create 1D images without nans
        nhi_image_likelihood = nhi_image_temp[indices]
        av_image_likelihood = av_image[indices]

        # Create model of Av with N(HI) and DGR
        av_image_model = nhi_image_likelihood * dgr

        model = av_image_model
        data = av_image_likelihood

        sn = np.sqrt(error**2 + f**2 * model**2)

        logL = -0.5 * np.sum((data - model)**2 / (sn**2) + \
                       np.log(2*np.pi*sn**2))

        logLs[count] = logL
        count += 1

# Normalize likelihoods
logLs -= np.max(logLs)
Ls = np.exp(logLs)
Ls /= np.sum(Ls)

# Get best-fit parameters
count = 0
max_L = 0
for velocity_width in velocity_widths:
    for dgr in dgrs:
        if Ls[count] > max_L:
            max_L = Ls[count]
            max_width = velocity_width
            max_dgr = dgr
        count += 1

print ''
print 'max likelihood = ', max_L
print 'max width = ', max_width
print 'max dgr = ', max_dgr
print 'max f = ', max_f

# Sort likelihoods and print top ten best likelihoods
Ls_sorted = np.sort(Ls)
print ''
print 'Top ten likelihoods = ', Ls_sorted[-10:]


