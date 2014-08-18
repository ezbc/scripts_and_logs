#!/usr/bin/python

''' Calculates the N(H) / Av correlation for the taurus molecular cloud. Uses
Xco factor from Paradis et al. (2012) A&A, 543, 103 to calculate N(H2).
Integrates GALFA HI image to determine N(HI).

'''


def plot_av_vs_nhi(nhi_image, av_image, limits=None,
        savedir='./', filename=None, show=False, scale=['linear', 'linear'],
        returnimage=False, hess_binsize=None, title='', plot_type='hexbin',
        color_scale='linear'):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib import cm
    from mpl_toolkits.axes_grid1 import ImageGrid

    # Drop the NaNs from the images
    indices = np.where((nhi_image == nhi_image) &\
                       (av_image == av_image) &\
                       (nhi_image > 0) &\
                       (av_image > 0))

    try:
        nhi_image_nonans = nhi_image[indices]
        av_image_nonans = av_image[indices]

        if type(av_image_error) is float:
            av_image_error_nonans = sd_image_error * \
                    np.ones(av_image[indices].shape)
        else:
            av_image_error_nonans = sd_image_error[indices]

        if type(nhi_image_error) is np.ndarray:
            nhi_image_error_nonans = nhi_image_error[indices]
        else:
            nhi_image_error_nonans = nhi_image_error * \
                    np.ones(nhi_image[indices].shape)
    except NameError:
        no_errors = True

    # Create figure
    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    fig_size = (4,4)
    font_scale = 10
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
              'figure.figsize': fig_size,
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    if plot_type == 'scatter':
    	cbar_mode = 'None'
    else:
    	cbar_mode = 'single'

    # Create figure
    plt.clf()
    fig = plt.figure()
    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1, 1),
                 ngrids=1,
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 cbar_mode=cbar_mode,
                 cbar_pad=0.1,
                 cbar_size=0.2)

    ax = imagegrid[0]

    if plot_type is 'hexbin':
        if color_scale == 'linear':
            image = ax.hexbin(nhi_image_nonans.ravel(),
                    av_image_nonans.ravel(),
                    mincnt=1,
                    xscale=scale[0],
                    yscale=scale[1],
                    cmap = cm.gist_stern)
            cb = ax.cax.colorbar(image,)
            # Write label to colorbar
            cb.set_label_text('Bin Counts',)
        elif color_scale == 'log':
            image = ax.hexbin(nhi_image_nonans.ravel(),
                av_image_nonans.ravel(),
                norm=matplotlib.colors.LogNorm(),
                mincnt=1,
                xscale=scale[0],
                yscale=scale[1],
                gridsize=(100,200),
                cmap = cm.gist_stern)
            cb = ax.cax.colorbar(image,)
            # Write label to colorbar
            cb.set_label_text('Bin Counts',)
        # Adjust color bar of density plot
        #cb = image.colorbar(image)
        #cb.set_label('Bin Counts')
    elif plot_type is 'scatter':
        image = ax.scatter(nhi_image_nonans.ravel(),
                av_image_nonans.ravel(),
                alpha=0.3,
                color='k'
                )
        ax.set_xscale(scale[0])
        ax.set_yscale(scale[1])

    if limits is not None:
    	ax.set_xlim(limits[0],limits[1])
    	ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'$N(HI)$ (10$^{20}$ cm$^{-2}$)')
    ax.set_ylabel(r'A$_{\rm V}$ (mag)')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_av_vs_nhi_grid(nhi_images, av_images, nhi_error_images=None,
        av_error_images=None, limits=None, savedir='./', filename=None,
        show=False, scale=['linear', 'linear'], returnimage=False,
        hess_binsize=None, title='', plot_type='hexbin', color_scale='linear'):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid

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
              'figure.figsize': (8, 8 * y_scaling),
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
                 share_all=True)

    # Cycle through lists
    for i in xrange(len(av_images)):
        av = av_images[i]
        nhi = nhi_images[i]
        av_error = av_error_images[i]
        nhi_error = nhi_error_images[i]
        #av_fit = av_fits[i]
        #nhi_fit = nhi_fits[i]

        # Drop the NaNs from the images
        if type(av_error) is float:
            indices = np.where((av == av) &\
                               (nhi == nhi)&\
                               (nhi > 0) &\
                               (av > 0))

        if type(av_error) is np.ndarray or \
                type(av_error) is np.ma.core.MaskedArray or \
                type(nhi_error) is np.ndarray or \
                type(nhi_error) is np.ma.core.MaskedArray:
            indices = np.where((av == av) &\
                               (nhi == nhi) &\
                               (nhi_error == nhi_error) &\
                               (av_error == av_error) &\
                               (nhi > 0) &\
                               (av > 0))

        av_nonans = av[indices]
        nhi_nonans = nhi[indices]

        if type(av_error) is np.ndarray:
            av_error_nonans = av_error[indices]
        else:
            av_error_nonans = np.array(av_error[indices])

        if type(nhi_error) is np.ndarray or \
                type(nhi_error) is np.ma.core.MaskedArray:
            nhi_error_nonans = nhi_error[indices]
        else:
            nhi_error_nonans = nhi_error * \
                    np.ones(nhi[indices].shape)

                # Create plot
        ax = imagegrid[i]

        image = ax.errorbar(nhi_nonans.ravel(),
                        av_nonans.ravel(),
                        xerr=(nhi_error_nonans.ravel()),
                        yerr=(av_error_nonans.ravel()),
                        alpha=0.3,
                        color='k',
                        marker='^',ecolor='k',linestyle='none',
                        markersize=4
                        )

        #if av_fit is not None:
        #    ax.plot(nhi_fit, av_fit,
        #            color = 'r')

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
                                  alpha=0.5),
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        )
            else:
                ax.annotate(r'\noindent$\phi_{\rm CNM}$ =' + \
                            r' %.2f' % (phi_cnm) + \
                            r'$^{+%.2f}_{-%.2f}$ \\' % (phi_cnm_error[0],
                                                     phi_cnm_error[1]) + \
                            r'Z = %.2f' % (Z) + \
                            r'$^{+%.2f}_{-%.2f}$ Z$_\odot$' % (Z_error[0],
                                                               Z_error[1]) + \
                            r'',
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
        ax.set_xlabel(r'$N(HI)$ (10$^{20}$ cm$^{-2}$)')
        ax.set_ylabel(r'A$_{\rm V}$ (mag)')
        ax.set_title(title)
        ax.grid(True)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
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
                                                     header=header)[:2]
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

def main():

    import grid
    import numpy as np
    from myimage_analysis import calculate_nhi
    from mycoords import make_velocity_axis
    import pyfits as pf
    import mygeometry as myg
    import json

    # parameters used in script
    # -------------------------
    # Regions
    # Options are 'ds9' or 'av_gradient'
    box_method = 'av_gradient'

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/dgr/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/hi/'
    co_dir = '/d/bip3/ezbc/taurus/data/cfa/'
    core_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'
    property_dir = '/d/bip3/ezbc/taurus/data/python_output/'

    av_data_planck, planck_header = pf.getdata(av_dir + \
                'taurus_av_planck_5arcmin.fits',
            header=True)
    av_data_error_planck, planck_header = pf.getdata(av_dir + \
                'taurus_av_error_planck_5arcmin.fits',
            header=True)

    # load GALFA HI
    hi_data, hi_header = pf.getdata(hi_dir + \
            'taurus_hi_galfa_cube_regrid_planckres.fits',
            header=True)
    velocity_axis = make_velocity_axis(hi_header)

    noise_cube, noise_header = pf.getdata(hi_dir + \
            'taurus_hi_galfa_cube_regrid_planckres_noise.fits', header=True)

    # define core properties
    with open(core_dir + 'taurus_core_properties.txt', 'r') as f:
        cores = json.load(f)

    cores = convert_core_coordinates(cores, planck_header)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'taurus_av_boxes_',
            header = planck_header)

    # Initialize lists
    av_images = []
    av_error_images = []
    nhi_images = []
    nhi_error_images = []

    for core in cores:
        print('\nCalculating for core %s' % core)
        if box_method == 'ds9':
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
        elif box_method == 'av_gradient':
            mask = myg.get_polygon_mask(av_data_planck,
                    cores[core]['box_vertices_rotated'])
        else:
        	raise ValueError('Method for boxes is either ds9 or av_gradient')

        indices = mask == 1

        # Get only the relevant pixels to decrease computation time
        hi_data_sub = np.copy(hi_data[:, indices])
        noise_cube_sub = np.copy(noise_cube[:, indices])
        av_data_planck_sub = np.copy(av_data_planck[indices])
        av_data_error_planck_sub = np.copy(av_data_error_planck[indices])

        # Derive N(HI) image
        nhi_image, nhi_image_error = calculate_nhi(cube=hi_data_sub,
                velocity_axis=velocity_axis,
                noise_cube=noise_cube_sub,
                velocity_range=cores[core]['hi_velocity_range'])

        nhi_images.append(nhi_image)
        nhi_error_images.append(nhi_image_error)
        av_images.append(av_data_planck_sub)
        av_error_images.append(av_data_error_planck_sub)

    plot_av_vs_nhi_grid(nhi_images,
                        av_images,
                        av_error_images=av_error_images,
                        nhi_error_images=nhi_error_images,
                        #limits=[0,14, 0,10],
                        scale=['linear', 'log'],
                        savedir=figure_dir,
                        plot_type='scatter',
                        filename='taurus_av_vs_nhi_panels.png',
                        color_scale='linear')

    # Derive N(HI) image
    nhi_image, nhi_image_error = calculate_nhi(cube=hi_data,
            velocity_axis=velocity_axis,
            noise_cube=noise_cube,
            velocity_range=cores[core]['hi_velocity_range'])

    # Plot correlation, similar to Figure 3 of Paradis et al. (2012)
    plot_av_vs_nhi(nhi_image,
            av_data_planck,
            savedir=figure_dir,
            limits=[0, 15, 0, 3],
            scale=['linear', 'linear'],
            filename='taurus_av_vs_nhi_global.png',
            color_scale='linear')

if __name__ == '__main__':
    main()


