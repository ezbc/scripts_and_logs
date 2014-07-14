#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the Taurus molecular cloud.
'''

import pyfits as pf
import numpy as np

def plot_pdf(av_image, limits=None, savedir='./', filename=None, show=True,
        scale=(0,0), n_bins=50, fit_gaussian=True, returnimage=False, title=''):

    ''' Plots a probability distribution function of an image.

    '''


    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from scipy.optimize import curve_fit

    # Drop the NaNs from the images
    indices = np.where(av_image == av_image)
    av_image_nonans = av_image[indices]

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
              'figure.figsize': (6, 6),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Derive the histograms
    bin_edges = np.logspace(-3, 3, num=n_bins, base=np.e)
    n = np.zeros(n_bins - 1)

    for i in xrange(n_bins - 1):
        bin_count = len(av_image_nonans[(av_image_nonans > bin_edges[i]) & \
                                    (av_image_nonans < bin_edges[i + 1])])
        n[i] = bin_count

    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

    n = np.append(n, 0)
    bin_centers = np.append(bin_centers,
            bin_centers[-1] + (bin_centers[-1] - bin_centers[-2]))

    # Normalize the bins
    mean_loc = np.argmin(np.abs(bin_centers - av_image_nonans.mean()))
    bin_centers = np.log(bin_centers / bin_centers[mean_loc])

    ax.errorbar(
        bin_centers,
        n,
        yerr = n**0.5,
        marker = '.',
        color = 'k',
        drawstyle = 'steps-mid'
    )

    # Fit a gausssian to the distribution
    if fit_gaussian:
        def gauss(x, a, x0, sigma):
            return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

        indices = np.where((bin_centers > -1.5) & \
                                       (bin_centers < 1))

        bin_centers_crop, n_crop = bin_centers[indices], n[indices]

        popt, pcov = curve_fit(gauss,
                               bin_centers_crop,
                               n_crop,
                               p0=[200, 0, 1],
                               maxfev=1000000)
        ax.plot(bin_centers,
                gauss(bin_centers, *popt),
                color = 'r')

    try:
        if scale[0] == 0:
            x_scale = 'linear'
        elif scale[0] == 1:
            x_scale = 'log'
        if scale[1] == 0:
            y_scale = 'linear'
        elif scale[1] == 1:
            y_scale = 'log'
    except IndexError('Scale must be tuple with 2 integer elements.'):
        pass

    ax.set_xscale(x_scale, nonposx = 'clip')
    ax.set_yscale(y_scale, nonposy = 'clip')


    if limits is not None:
    	ax.set_xlim(limits[0],limits[1])
    	ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'ln(A$_{\rm V}$ / $\bar{\rm A}_{\rm V}$ (mag))',)
    ax.set_ylabel(r'N',)
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()

def plot_pdfs(av_images, limits=None, savedir='./', filename=None, show=True,
        scale=(0,0), n_bins=50, fit_gaussian=True, returnimage=False,
        title='', core_names=''):

    ''' Plots a probability distribution function of an image.

    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from scipy.optimize import curve_fit
    from mpl_toolkits.axes_grid1 import ImageGrid

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
              'figure.figsize': (10, 10),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    n = int(np.ceil(len(av_images)**0.5))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(n, n),
                 ngrids=len(av_images),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    for j, av_image in enumerate(av_images):
        ax = imagegrid[j]

        # Drop the NaNs from the images
        indices = np.where(av_image == av_image)
        av_image_nonans = av_image[indices]

        # Derive the histograms
        bin_edges = np.logspace(-3, 3, num=n_bins, base=np.e)
        n = np.zeros(n_bins - 1)

        for i in xrange(n_bins - 1):
            bin_count = len(av_image_nonans[(av_image_nonans > bin_edges[i]) & \
                                        (av_image_nonans < bin_edges[i + 1])])
            n[i] = bin_count

        bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

        n = np.append(n, 0)
        bin_centers = np.append(bin_centers,
                bin_centers[-1] + (bin_centers[-1] - bin_centers[-2]))

        # Normalize the bins
        mean_loc = np.argmin(np.abs(bin_centers - av_image_nonans.mean()))
        bin_centers = np.log(bin_centers / bin_centers[mean_loc])

        ax.errorbar(
            bin_centers,
            n,
            yerr = n**0.5,
            marker = '.',
            color = 'k',
            drawstyle = 'steps-mid'
        )

        # Fit a gausssian to the distribution
        if fit_gaussian:
            def gauss(x, a, x0, sigma):
                return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

            indices = np.where((bin_centers > -1.5) & \
                                           (bin_centers < 1))

            indices = np.where(bin_centers == bin_centers)

            bin_centers_crop, n_crop = bin_centers[indices], n[indices]

            popt, pcov = curve_fit(gauss,
                                   bin_centers_crop,
                                   n_crop,
                                   p0=[200, 0, 1],
                                   maxfev=1000000)
            ax.plot(bin_centers,
                    gauss(bin_centers, *popt),
                    color = 'r')

        try:
            if scale[0] == 0:
                x_scale = 'linear'
            elif scale[0] == 1:
                x_scale = 'log'
            if scale[1] == 0:
                y_scale = 'linear'
            elif scale[1] == 1:
                y_scale = 'log'
        except IndexError('Scale must be tuple with 2 integer elements.'):
            pass

        ax.set_xscale(x_scale, nonposx = 'clip')
        ax.set_yscale(y_scale, nonposy = 'clip')


        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'ln(A$_{\rm V}$ / $\bar{\rm A}_{\rm V}$)',)
        ax.set_ylabel(r'N',)
        ax.set_title(core_names[j])
        ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()

def convert_core_coordinates(cores, header):

    for core in cores:
        cores[core].update({'box_pixel': 0})
        cores[core].update({'center_pixel': 0})

    	box_wcs = cores[core]['box_wcs']
    	box_pixel = len(box_wcs) * [0,]
    	center_wcs = cores[core]['center_wcs']

        # convert centers to pixel coords
        cores[core]['center_pixel'] = get_pix_coords(ra=center_wcs[0],
                                                     dec=center_wcs[1],
                                                     header=header)
        # convert box corners to pixel coords
        for i in range(len(box_wcs)/2):
        	pixels = get_pix_coords(ra=box_wcs[2*i], dec=box_wcs[2*i + 1],
        	        header=header)
        	box_pixel[2*i], box_pixel[2*i + 1] = int(pixels[0]), int(pixels[1])
        cores[core]['box_pixel'] = box_pixel

    return cores

''' DS9 Region and Coordinate Functions
'''

def convert_core_coordinates(cores, header):

    for core in cores:
        cores[core].update({'box_pixel': 0})
        cores[core].update({'center_pixel': 0})

        box_wcs = cores[core]['box_wcs']
        box_pixel = len(box_wcs) * [0,]
        center_wcs = cores[core]['center_wcs']

        # convert centers to pixel coords
        cores[core]['center_pixel'] = get_pix_coords(ra=center_wcs[0],
                                                     dec=center_wcs[1],
                                                     header=header)[:2]
        # convert box corners to pixel coords
        for i in range(len(box_wcs)/2):
            pixels = get_pix_coords(ra=box_wcs[2*i], dec=box_wcs[2*i + 1],
                    header=header)
            box_pixel[2*i], box_pixel[2*i + 1] = int(pixels[0]), int(pixels[1])
        cores[core]['box_pixel'] = box_pixel

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
    region = pyr.open(filename)

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
        box_center_pixel = get_pix_coords(ra = region[0],
                                          dec = region[1],
                                          header = header)
        box_center_pixel = (int(box_center_pixel[1]), int(box_center_pixel[0]))
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
    import myclumpfinder as clump_finder
    reload(clump_finder)
    import mygeometry as myg
    reload(myg)
    import json

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/cores/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    core_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'

    # load 2mass Av and GALFA HI images, on same grid
    av_data_k09, av_header = load_fits(av_dir + 'taurus_av_k09_regrid_planckres.fits',
            return_header=True)

    av_data_k09_orig, av_header = load_fits(av_dir + \
                                          'taurus_av_kainulainen2009.fits',
                                       return_header=True)

    av_data_k09_orig[av_data_k09_orig == -1] = np.NaN

    # load Planck Av and GALFA HI images, on same grid
    av_data_planck, h = load_fits(av_dir + \
                'taurus_av_planck_5arcmin.fits',
            return_header=True)

    # define core properties
    with open(core_dir + 'taurus_core_properties.txt', 'r') as f:
        cores = json.load(f)

    cores = convert_core_coordinates(cores, h)

    if True:
        print('Calculating PDF for global map')

        figure_types = ['pdf', 'png']
        for figure_type in figure_types:
            if 1:
                plot_pdf(av_data_planck,
                        limits = [-2, 3.1, 1, 10**5],
                        savedir=figure_dir + 'panel_cores/',
                        scale=(0,1),
                        filename='taurus_av_pdf_global_planck.%s' % \
                                figure_type,
                        title=r'A$_{\rm V}$ PDF of ' + \
                                'Taurus',
                        show=False)

                plot_pdf(av_data_k09,
                        limits = [-2, 3.1, 1, 10**5],
                        savedir=figure_dir + 'panel_cores/',
                        scale=(0,1),
                        filename='taurus_av_pdf_global_k09.%s' % \
                                figure_type,
                        title=r'A$_{\rm V}$ PDF of ' + \
                                'Taurus',
                        show=False)
            if 1:
                # Create PDF with original Av image from Jouni
                plot_pdfs((av_data_k09_orig, av_data_k09),
                        limits = [-2, 3.1, 1, 10**5],
                        savedir=figure_dir + 'panel_cores/',
                        scale=(0,1),
                        core_names=('Original', 'Planck Regrid'),
                        filename='taurus_av_pdf_global_k09_orig_vs_regrid.%s'%\
                                figure_type,
                        title=r'A$_{\rm V}$ PDF of Taurus',
                        show=False)

        cores = load_ds9_region(cores,
                filename_base = region_dir + 'taurus_av_boxes_',
                header = h)

        av_data_list_k09 = []
        av_data_list_planck = []
        core_name_list = []

        for core in cores:
            print('Calculating PDF for core %s' % core)

            # Grab the mask from the DS9 regions
            xy = cores[core]['box_center_pix']
            box_width = cores[core]['box_width']
            box_height = cores[core]['box_height']
            box_angle = cores[core]['box_angle']
            mask = myg.get_rectangular_mask(av_data_planck,
        	        xy[0], xy[1],
                    width = box_width,
                    height = box_height,
                    angle = box_angle)

            # Get indices where there is no mask, and extract those pixels
            indices = np.where(mask == 1)
            av_data_k09_sub = av_data_k09[indices]
            av_data_planck_sub = av_data_planck[indices]

            if 1:
                figure_types = ['pdf', 'png']
                for figure_type in figure_types:
                    plot_pdf(av_data_planck_sub,
                            limits = [-1,4,1,1000],
                            savedir=figure_dir + 'individual_cores/',
                            scale=(0,1),
                            filename='taurus_av_pdf_' + core + '_planck.%s' % \
                                    figure_type,
                            title=r'A$_{\rm V}$ PDF of ' + \
                                    'Taurus Core ' + core,
                            show=False)

            av_data_list_k09.append(av_data_k09_sub)
            av_data_list_planck.append(av_data_planck_sub)
            core_name_list.append(core)

        if 1:
            for figure_type in figure_types:
                plot_pdfs(av_data_list_k09,
                          limits = [-2,2.9,1,100],
                          savedir=figure_dir + 'panel_cores/',
                          scale=(0,1),
                          filename='taurus_av_pdfs_k09.%s' % figure_type,
                          title=r'A$_{\rm V}$ PDF of Taurus Core ',
                          core_names=core_name_list,
                          show=False)

                plot_pdfs(av_data_list_planck,
                          limits = [-2,2.9,1,100],
                          savedir=figure_dir + 'panel_cores/',
                          scale=(0,1),
                          filename='taurus_av_pdfs_planck.%s' % figure_type,
                          title=r'A$_{\rm V}$ PDF of Taurus Core ',
                          core_names=core_name_list,
                          show=False)

if __name__ == '__main__':
    main()



