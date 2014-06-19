#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the Taurus molecular cloud.
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

        ax.set_xscale('linear')
        ax.set_yscale(scale)

    if limits is not None:
    	ax.set_xlim(limits[0],limits[1])
    	ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'A$_{\rm V}$ (mag)',)
    ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_sd_vs_av_grid(sd_images, av_images,
        sd_image_errors=None, av_image_errors=None, limits=None,
        savedir='./', filename=None, show=True, scale=('linear', 'linear'),
        returnimage=False, title=None, core_names='', plot_type='scatter',
        hexbin_range=None):
    ''' Plots N(HI) as a function of Av for individual pixels in an N(HI) image
    and an Av image.
    '''

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

    if plot_type == 'scatter':
        imagegrid = ImageGrid(fig, (1,1,1),
                     nrows_ncols=(n, n),
                     ngrids=len(av_images),
                     axes_pad=0,
                     aspect=False,
                     label_mode='L',
                     share_all=True)
    elif plot_type == 'hexbin':
        imagegrid = ImageGrid(fig, (1,1,1),
                              nrows_ncols=(n, n),
                              ngrids=len(av_images),
                              cbar_mode='single',
                              cbar_location='right',
                              cbar_pad=0.1,
                              cbar_size=0.25,
                              #cbar_set_cax=True,
                              axes_pad=0,
                              aspect=False,
                              label_mode='L',
                              share_all=True)

    for i in xrange(len(av_images)):
    	sd_image = sd_images[i]
    	av_image = av_images[i]
    	sd_image_error = sd_image_errors[i]
    	av_image_error = av_image_errors[i]

        # Drop the NaNs from the images
        indices = np.where((sd_image == sd_image) &\
                           (av_image == av_image)&\
                           (av_image > 0) &\
                           (sd_image > 0))

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

        # Create plot
        ax = imagegrid[i]

        if plot_type == 'scatter':
            image = ax.errorbar(av_image_nonans.ravel(),
                    sd_image_nonans.ravel(),
                    xerr=(av_image_error_nonans.ravel()),
                    yerr=(sd_image_error_nonans.ravel()),
                    alpha=0.1,
                    color='k',
                    marker='^',ecolor='k',linestyle='none',
                    markersize=2
                    )
            ax.set_xscale(scale[0], nonposx = 'clip')
            ax.set_yscale(scale[1], nonposy = 'clip')
        elif plot_type == 'hexbin':
            image = ax.hexbin(av_image_nonans.ravel(),
                sd_image_nonans.ravel(),
                #norm=matplotlib.colors.LogNorm(), # logscale?
                mincnt=1,
                yscale=scale[0],
                xscale=scale[1],
                cmap = matplotlib.cm.gist_stern,
                vmin=hexbin_range[0],
                vmax=hexbin_range[1])

            cbar = imagegrid[0].cax.colorbar(image)
            cbar.set_label_text('Bin Counts')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'A$_{\rm V}$ (mag)',)
        ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
        ax.annotate(core_names[i].capitalize(),
                xytext=(0.7, 0.9),
                xy=(0.7, 0.9),
                textcoords='axes fraction',
                xycoords='axes fraction',
                color='k'
                )
        ax.grid(True)

    if title is not None:
    	fig.suptitle(title, fontsize=fontScale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
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
    from mycoords import make_velocity_axis

    if velocity_axis is None:
        velocity_axis = make_velocity_axis(header)

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

def calculate_noise_cube(cube=None, velocity_axis=None,
            velocity_noise_range=[-110,-90,90,110], header=None, Tsys=30.,
            filename=None):

    """ Calcuates noise envelopes for each pixel in a cube
    """

    import numpy as np
    import pyfits as pf
    from mycoords import make_velocity_axis

    if velocity_axis is None:
        velocity_axis = make_velocity_axis(header)


    noise_cube = np.zeros(cube.shape)
    for i in range(cube.shape[1]):
        for j in range(cube.shape[2]):
            profile = cube[:,i,j]
            noise = calculate_noise(profile, velocity_axis,
                    velocity_noise_range)
            #noise = 0.1 # measured in line free region
            noise_cube[:,i,j] = calculate_noise_scale(Tsys,
                    profile, noise=noise)

    if filename is not None:
        pf.writeto(filename, noise_cube, header=header)

    return noise_cube

def calculate_noise(profile, velocity_axis, velocity_range):
    """ Calculates rms noise of Tile profile given velocity ranges.
    """
    import numpy as np

    std = 0

    # calculate noises for each individual region
    for i in range(len(velocity_range) / 2):
        velMin = velocity_range[2*i + 0]
        velMax = velocity_range[2*i + 1]

        noise_region = np.where((velocity_axis >= velMin) & \
                        (velocity_axis <= velMax))

        std += np.std(profile[noise_region])

    std /= len(velocity_range) / 2
    return std

def calculate_noise_scale(Tsys, profile, noise=None):
    """ Creates an array for scaling the noise by (Tsys + Tb) / Tsys
    """
    import numpy as np
    n = np.zeros(len(profile))
    n = (Tsys + profile) / Tsys * noise

    return n

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

def calculate_nh2(nhi_image = None, av_image = None, dgr = 1.1e-1):

    ''' Calculates the total gas column density given N(HI), A_v and a
    dust-to-gas ratio.

    Parameters
    ----------
    '''

    import numpy as np

    nh_image = 0.5 * (av_image / dgr - nhi_image)

    return nh_image

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

def main():

    import grid
    import numpy as np
    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)
    import mygeometry as myg
    reload(myg)



    cloud_list = ('taurus', 'perseus', 'california')
    data_dict = {}

    for i, cloud in enumerate(cloud_list):
        # define directory locations
        output_dir = '/d/bip3/ezbc/%s/data/python_output/nhi_av/' % cloud
        figure_dir = '/d/bip3/ezbc/multicloud/figures/'
        av_dir = '/d/bip3/ezbc/%s/data/av/' % cloud
        hi_dir = '/d/bip3/ezbc/%s/data/hi/' % cloud
        core_dir = output_dir + 'core_arrays/'
        region_dir = '/d/bip3/ezbc/%s/data/python_output/ds9_regions/' % cloud

        cloud_dict = {}

        # Load av maps from Taurus, California, and Perseus
        cloud_dict['av_data'], cloud_dict['av_header'] = load_fits(
                '/d/bip3/ezbc/%s/data/av/%s_av_planck_5arcmin.fits' % \
                        (cloud, cloud),
                return_header=True)

        cloud_dict['av_data_error'] = load_fits(
                '/d/bip3/ezbc/%s/data/av/%s_av_error_planck_5arcmin.fits' % \
                        (cloud, cloud),
                return_header=False)

        # Load HI maps from Taurus, California, and Perseus
        cloud_dict['hi_data'], cloud_dict['hi_header'] = load_fits(
                '/d/bip3/ezbc/{0}/data/hi/{0}_hi_galfa_'.format(cloud) + \
                        'cube_regrid_planckres.fits',
                return_header=True)


        # Plot NHI vs. Av for a given velocity range
        filename = '%s_hi_noise_galfa_cube_regrid_planckres.fits' % cloud
        if not path.isfile(hi_dir + filename):
            cloud_dict['noise_cube'] = \
                calculate_noise_cube(cube=cloud_dict['hi_data'],
                        velocity_noise_range=[90,110],
                        header=cloud_dict['hi_header'],
                        Tsys=30.,
                        filename=hi_dir + filename)
        else:
            cloud_dict['noise_cube'], cloud_dict['noise_header'] = \
                    load_fits(hi_dir + filename,
                        return_header=True)

        # calculate nhi and error maps, write nhi map to fits file
        cloud_dict['nhi_image'], cloud_dict['nhi_image_error'] = \
                calculate_nhi(cube=cloud_dict['hi_data'],
                    noise_cube=cloud_dict['noise_cube'],
                    velocity_range=[-100,100],
                    return_nhi_error=True,
                    fits_filename=hi_dir + '%s_nhi_galfa_5arcmin.fits' % cloud,
                    fits_error_filename=hi_dir + \
                            '%s_nhi_error_galfa_5arcmin.fits' % cloud,
                    header=cloud_dict['hi_header'])

        # convert to column density to surface density
        cloud_dict['hi_sd_image'] = calculate_sd(cloud_dict['nhi_image'],
                sd_factor=1/1.25)
        cloud_dict['hi_sd_image_error'] =\
                calculate_sd(cloud_dict['nhi_image_error'],\
                    sd_factor=1/1.25)

        # Write the cloud properties to the main dictionary
        data_dict[cloud] = cloud_dict


    sd_image_list = []
    av_image_list = []
    core_name_list = []
    sd_image_error_list = []
    av_image_error_list = []

    for cloud in data_dict:
        sd_image_list.append(data_dict[cloud]['hi_sd_image'])
        av_image_list.append(data_dict[cloud]['av_data'])
        core_name_list.append(cloud)
        sd_image_error_list.append(data_dict[cloud]['hi_sd_image_error'])
        av_image_error_list.append(data_dict[cloud]['av_data_error'])

    #figure_types = ['png', 'pdf']
    figure_types = ('png',)
    for figure_type in figure_types:
        # scatter plot
        plot_sd_vs_av_grid(sd_image_list,
                        av_image_list,
                        sd_image_errors = sd_image_error_list,
                        av_image_errors = av_image_error_list,
                        limits = [-3,40,4,32],
                        savedir=figure_dir,
                        scale=('linear', 'linear'),
                        filename='multicloud_sd_vs_av_scatter_planck.%s'%\
                                figure_type,
                        show=False,
                        core_names=core_name_list,
                        title=r'Global $\Sigma_{HI}$ vs. Planck A$_{\rm V}$' + \
                                ' of GMCs')

        # density plot
        plot_sd_vs_av_grid(sd_image_list,
                        av_image_list,
                        sd_image_errors = sd_image_error_list,
                        av_image_errors = av_image_error_list,
                        limits = [-3,40,4,32],
                        savedir=figure_dir,
                        scale=('linear', 'linear'),
                        plot_type='hexbin',
                        hexbin_range=(0,240),
                        filename='multicloud_sd_vs_av_hexbin_planck.%s'%\
                                figure_type,
                        show=False,
                        core_names=core_name_list,
                        title=r'Global $\Sigma_{HI}$ vs. Planck A$_{\rm V}$' + \
                                ' of GMCs')


if __name__ == '__main__':
    main()



