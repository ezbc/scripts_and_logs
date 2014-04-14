#!/usr/bin/python

''' Calculates the CfA CO 2nd moment (velocity dispersion) in Perseus MC
'''

import pyfits as pf
import numpy as np

def plot_sd_vs_av(sd_image, av_image,
        sd_image_error=None, av_image_error=None, limits=None,
        savedir='./', filename=None, show=True, scale='linear',
        returnimage=False, hess_binsize=None, title='', plot_type='hexbin'):
    ''' Plots N(HI) as a function of Av for individual pixels in an N(HI) image
    and an Av image.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    # Drop the NaNs from the images
    indices = np.where((sd_image == sd_image) &\
                       (av_image == av_image)&\
                       (av_image > 0) &\
                       (sd_image > -5))

    sd_image_nonans = sd_image[indices]
    av_image_nonans = av_image[indices]

    if type(sd_image_error) is np.ndarray:
        sd_image_error_nonans = sd_image_error[indices]
    else:
        sd_image_error_nonans = np.array(sd_image_error[indices])

    if type(av_image_error) is np.ndarray:
        av_image_error_nonans = av_image_error[indices]
    else:
        av_image_error_nonans = av_image_error * \
                np.ones(av_image[indices].shape)

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
    if sd_image_error is None:
        if plot_type is 'hexbin':
            image = ax.hexbin(av_image_nonans.ravel(),
                    sd_image_nonans.ravel(),
                    norm=matplotlib.colors.LogNorm(),
                    mincnt=1,
                    yscale='log',
                    xscale='log')
            # Adjust color bar of density plot
            cb = plt.colorbar(image)
            cb.set_label('Counts')
        elif plot_type is 'scatter':
            image = ax.scatter(av_image_nonans.ravel(),
                    sd_image_nonans.ravel(),
                    alpha=0.3,
                    color='k'
                    )
            ax.set_xscale('log')
            ax.set_yscale('log')
    else:
        image = ax.errorbar(av_image_nonans.ravel(),
                sd_image_nonans.ravel(),
                xerr=(av_image_error_nonans.ravel()),
                yerr=(sd_image_error_nonans.ravel()),
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
    ax.set_xlabel('A$_v$ (mag)',)
    ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def calculate_nhi(cube = None, velocity_axis = None, velocity_range = [],
        return_nhi_error = True, noise_cube = None,
        velocity_noise_range=[90,100], Tsys = 30., header = None,
        fits_filename = None, fits_error_filename = None, verbose = True):

    ''' Calculates an N(HI) image given a velocity range within which to
    include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to
    fits_filename : str
        If specified, and a header is provided, the nhi image will be written.
    header : pyfits.Header
        Header from cube.

    '''

    import numpy as np

    # Calculate NHI from cube if set
    if cube is not None and velocity_axis is not None:
        image = np.empty((cube.shape[1],
                          cube.shape[2]))
        image[:,:] = np.NaN
        indices = np.where((velocity_axis > velocity_range[0]) & \
                (velocity_axis < velocity_range[1]))[0]
        image[:,:] = cube[indices,:,:].sum(axis=0)
        # Calculate image error
        if return_nhi_error:
            image_error = np.empty((cube.shape[1],
                              cube.shape[2]))
            image_error[:,:] = np.NaN
            image_error[:,:] = (noise_cube[indices,:,:]**2).sum(axis=0)**0.5

    # NHI in units of 1e20 cm^-2
    nhi_image = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2


    if fits_filename is not None and header is not None:
        if verbose:
        	print('Writing N(HI) image to FITS file %s' % fits_filename)
        header['BUNIT'] = '1e20 cm^-2'
        header.remove('CDELT3')
        header.remove('CRVAL3')
        header.remove('CRPIX3')
        header.remove('CTYPE3')
        header.remove('NAXIS3')
        header['NAXIS'] = 2

        pf.writeto(fits_filename, image*1.823e-2, header = header, clobber =
                True, output_verify = 'fix')

    if fits_error_filename is not None and header is not None:
        if verbose:
        	print('Writing N(HI) error image to FITS file %s' % fits_filename)

        pf.writeto(fits_error_filename, image_error * 1.823e-2, header =
                header, clobber = True, output_verify = 'fix')

    if return_nhi_error:
        nhi_image_error = np.ma.array(image_error,
                mask=np.isnan(image_error)) * 1.823e-2
        return nhi_image, nhi_image_error
    else:
        return nhi_image

def calculate_sd(image, sd_factor=1/1.25):

    ''' Calculates a surface density image given a velocity range within which
    to include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to '''

    import numpy as np

    # NHI in units of 1e20 cm^-2
    sd_image = image * sd_factor

    return sd_image

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

def load_fits(filename,return_header=False):
    ''' Loads a fits file.
    '''

    import pyfits as pf

    f = pf.open(filename)
    if return_header:
        return f[0].data,f[0].header
    else:
        return f[0].data

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

def plot_mom2_image(mom2_image=None, header=None, contour_image=None,
        cores=None, title=None, limits=None,
        contours=None, boxes=False, savedir='./', filename=None, show=True):

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

    # Create figure instance
    fig = plt.figure()

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    # create axes
    ax = imagegrid[0]
    cmap = cm.winter # colormap
    cmap = cm.gist_stern # colormap
    # show the image
    im = ax.imshow(mom2_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,)

    # Asthetics
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',)
    ax.set_ylabel('Declination (J2000)',)
    if title is not None:
    	ax.set_title(title)

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
    cb.set_label_text(r'km/s',)

    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    core_anno_color = 'c'
    if cores is not None:
        for core in cores:
            pix_coords = cores[core]['center_pixel']

            ax.scatter(pix_coords[0],pix_coords[1],
                    color=core_anno_color,
                    s=200, marker='+',
                    linewidths=2)

            ax.annotate(core,
                    xy=[pix_coords[0], pix_coords[1]],
                    xytext=(5,5),
                    textcoords='offset points',
                    color=core_anno_color,
                    )

            if boxes:
                rect = ax.add_patch(Polygon(
                    cores[core]['box_vertices'][:, ::-1],
                        facecolor='none',
                        edgecolor='k' ))

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        fig.show()

def main():

    import grid
    import numpy as np
    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)
    import mygeometry as myg
    reload(myg)
    from mycoords import make_velocity_axis
    import mymath
    reload(mymath)

    # define directory locations
    output_dir = '/d/bip3/ezbc/perseus/data/python_output/co_dispersion/'
    figure_dir = '/d/bip3/ezbc/perseus/figures/'
    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    hi_dir = '/d/bip3/ezbc/perseus/data/galfa/'
    cfa_dir = '/d/bip3/ezbc/perseus/data/cfa/'
    core_dir = output_dir + 'core_arrays/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'

    # load 2mass Av and GALFA HI images, on same grid
    cfa_data, cfa_header = load_fits(cfa_dir + \
                'perseus_cfa_cube_galfa_regrid.fits',
            return_header=True)

    #av_data += - 0.4 # subtracts background of 0.4 mags
    hi_data, hi_header = load_fits(hi_dir + \
                'perseus_galfa_cube_bin_3.7arcmin.fits',
            return_header = True)

    # make the velocity axis
    hi_velocity_axis = make_velocity_axis(hi_header)
    cfa_velocity_axis = make_velocity_axis(cfa_header)

    cfa_mom2 = mymath.calc_moment(cfa_data, moment = 2, spectral_axis = 0)

    # define core properties
    cores = {'L1495':
                {'center_wcs': [(4,14,0), (28, 11, 0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(4,16,30), (27,44,30), (4,5,20), (28,28,33)]
                 },
             'L1495A':
                {'center_wcs': [(4,18,0), (28,23., 0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(4,28,23),(28,12,50),(4,16,23),(29,46,5)],
                 },
             'B213':
                {'center_wcs': [(4, 19, 0), (27, 15,0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(4,22,27), (26,45,47),(4,5,25),(27,18,48)],
                },
             'B220':
                {'center_wcs': [(4, 41, 0.), (26,7,0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,47,49),(25,31,13),(4,40,37),(27,31,17)],
                 },
             'L1527':
                {'center_wcs': [(4, 39, 0.), (25,47, 0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,40,13), (24,46,38), (4,34,35), (25,56,7)],
                 },
             'B215':
                {'center_wcs': [(4, 23, 0), (25, 3, 0)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,24,51), (22,36,7), (4,20,54), (25,26,31)],
                 },
             'L1524':
                {'center_wcs': [(4,29,0.), (24,31.,0)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,31,0), (22,4,6), (4,25,33), (25,0,55)],
                 }
                }

    cores = convert_core_coordinates(cores, hi_header)

    if True:
        hsd_limits =[0.1,300]
        hisd_limits = [2,20]
        av_limits =[0.01,100]
        nhi_limits = [2,20]

        cores = load_ds9_region(cores,
                filename_base = region_dir + 'taurus_av_boxes_',
                header = hi_header)

        # save cores for later
        if 0:
            mask = np.zeros((cfa_mom2.shape))
            for core in cores:
                print('Calculating for core %s' % core)

                av_limits = [0.01,100]

                # Grab the mask
                xy = cores[core]['box_center_pix']
                box_width = cores[core]['box_width']
                box_height = cores[core]['box_height']
                box_angle = cores[core]['box_angle']
                mask += myg.get_rectangular_mask(cfa_mom2,
                        xy[0], xy[1],
                        width = box_width,
                        height = box_height,
                        angle = box_angle)

                cores[core]['box_vertices'] = myg.get_rect(
                            xy[0], xy[1],
                            width = box_width,
                            height = box_height,
                            angle = box_angle,)

                indices = np.where(mask == 1)

            mask[mask > 1] = 1
            cfa_mom2[mask == 0] = np.nan

        # currently have variance, need dispersion
        cfa_mom2 = cfa_mom2**0.5


        # Plot
        plot_mom2_image(mom2_image = cfa_mom2,
                header = cfa_header,
                boxes = False,
                #limits=[128,37,308,206],
                title = r'Perseus: $\sigma_{\rm CO}$ map with core ' + \
                        'boxed-regions.',
                savedir = figure_dir,
                filename='perseus_co_veldisp_map.png',
                show=True)

if __name__ == '__main__':
    main()



