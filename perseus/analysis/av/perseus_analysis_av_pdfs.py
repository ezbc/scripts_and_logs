#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the perseus molecular cloud.
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

    bin_centers = np.log(bin_centers)

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
        popt, pcov = curve_fit(gauss, bin_centers, n, p0=[200, 1, 2])
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
    ax.set_xlabel(r'ln(A$_{\rm V}$ (mag))',)
    ax.set_ylabel(r'N',)
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()

def convert_core_coordinates(cores, header):

    for core in cores:
        cores[core].update({'box_pixel': 0})
        cores[core].update({'center_pixel': 0})

        try:
            box_wcs = cores[core]['box_wcs']
            box_pixel = len(box_wcs) * [0,]
            # convert box corners to pixel coords
            for i in range(len(box_wcs)/2):
                pixels = get_pix_coords(ra=box_wcs[2*i], dec=box_wcs[2*i + 1],
                        header=header)
                box_pixel[2*i], box_pixel[2*i + 1] = int(pixels[0]),\
                        int(pixels[1])
            cores[core]['box_pixel'] = box_pixel
        except TypeError:
            do_nothing = True


    	center_wcs = cores[core]['center_wcs']

        # convert centers to pixel coords
        cores[core]['center_pixel'] = get_pix_coords(ra=center_wcs[0],
                                                     dec=center_wcs[1],
                                                     header=header)

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

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec).
    '''

    import pywcsgrid2 as wcs
    import pywcs

    # convert to degrees
    if type(ra) is tuple and type(dec) is tuple:
        ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)
    else:
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

def load_ds9_region(cores, filename_base = 'perseus_av_boxes_', header=None):

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

def main():

    import grid
    import numpy as np
    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)
    import mygeometry as myg
    reload(myg)


    # define directory locations
    output_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/perseus/figures/cores/'
    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    hi_dir = '/d/bip3/ezbc/perseus/data/galfa/'
    core_dir = output_dir + 'core_arrays/'
    region_dir = '/d/bip3/ezbc/perseus/data/python_output/ds9_regions/'

    # load Planck Av and GALFA HI images, on same grid
    av_data_planck, h = load_fits(av_dir + \
                'perseus_planck_av_regrid.fits',
            return_header=True)

    cores = {'IC348':
                {'center_wcs': [(3, 44, 0), (32, 8, 0)],
                 'map': None,
                 'threshold': None,
                 'box_wcs': [(3,46,13), (26,3,24), (3,43,4), (32,25,41)],
                 },
             'NGC1333':
                {'center_wcs': [(3, 29, 11), (31, 16, 53)],
                 'map': None,
                 'threshold': None,
                 'box_wcs': None,
                 },
             'B4':
                {'center_wcs': [(3, 44, 18), (32, 05, 20)],
                 'map': None,
                 'threshold': None,
                 'box_wcs': None,
                 },
             'B5':
                {'center_wcs': [(3, 47, 34), (32, 48, 17)],
                 'map': None,
                 'threshold': None,
                 'box_wcs': None,
                 },
             #'':
             #   {'center_wcs': [],
             #    'map': None,
             #    'threshold': None,
             #    'box_wcs': None,
             #    },
            }


    # define core properties
    cores = convert_core_coordinates(cores, h)


    if True:
        print('Calculating PDF for global map')

        figure_types = ['pdf', 'png']
        for figure_type in figure_types:
            plot_pdf(av_data_planck,
                    limits = [-2, 3.1, 1, 10**5],
                    savedir=figure_dir,
                    scale=(0,1),
                    filename='perseus_av_pdf_global_planck.%s' % \
                            figure_type,
                    title=r'A$_{\rm V}$ PDF of ' + \
                            'perseus',
                    show=False)

        cores = load_ds9_region(cores,
                filename_base = region_dir + 'perseus_av_boxes_',
                header = h)

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
            av_data_planck_sub = av_data_planck[indices]

            figure_types = ['pdf', 'png']
            for figure_type in figure_types:
                plot_pdf(av_data_planck_sub,
                        limits = [-1,4,1,1000],
                        savedir=figure_dir,
                        scale=(0,1),
                        filename='perseus_av_pdf' + core + '_planck.%s' % \
                                figure_type,
                        title=r'A$_{\rm V}$ PDF of ' + \
                                'perseus Core ' + core,
                        show=False)

if __name__ == '__main__':
    main()



