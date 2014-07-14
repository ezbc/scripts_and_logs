#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the california molecular cloud.
'''

import pyfits as pf
import numpy as np

def plot_sd_vs_av(sd_image, av_image,
        sd_image_error=None, av_image_error=None, limits=None,
        fit = True,
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

    # Perform fit
    if fit:
        av = av_image_nonans.ravel()
        nhi = sd_image_nonans.ravel()

        nhi_fit, av_extended = fit_krumholz(av, nhi, (0, 200))

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
        if fit:
            ax.plot(av_extended, nhi_fit)

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

def plot_hi_vs_h(hi_sd, h_sd, hi_sd_error=None, h_sd_error=None,
        hi_sd_fit = None, h_sd_fit = None,
        limits = None, fit = True, savedir='./', filename=None, show=True,
        scale = 'linear', title=''):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    # Drop the NaNs from the images
    if type(hi_sd_error) is float:
        indices = np.where((hi_sd == hi_sd) &\
                           (h_sd == h_sd)&\
                           (h_sd > 0) &\
                           (hi_sd > 0))

    if type(hi_sd_error) is np.ndarray or \
            type(hi_sd_error) is np.ma.core.MaskedArray or \
            type(h_sd_error) is np.ndarray or \
            type(h_sd_error) is np.ma.core.MaskedArray:
        indices = np.where((hi_sd == hi_sd) &\
                           (h_sd == h_sd) &\
                           (h_sd_error == h_sd_error) &\
                           (hi_sd_error == hi_sd_error) &\
                           (h_sd > 0) &\
                           (hi_sd > 0))

    hi_sd_nonans = hi_sd[indices]
    h_sd_nonans = h_sd[indices]

    if type(hi_sd_error) is np.ndarray:
        hi_sd_error_nonans = hi_sd_error[indices]
    else:
        hi_sd_error_nonans = np.array(hi_sd_error[indices])

    if type(h_sd_error) is np.ndarray or \
            type(h_sd_error) is np.ma.core.MaskedArray:
        h_sd_error_nonans = h_sd_error[indices]
    else:
        h_sd_error_nonans = h_sd_error * \
                np.ones(h_sd[indices].shape)

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
    image = ax.errorbar(h_sd_nonans.ravel(),
            hi_sd_nonans.ravel(),
            xerr=(h_sd_error_nonans.ravel()),
            yerr=(hi_sd_error_nonans.ravel()),
            alpha=0.3,
            color='k',
            marker='^',ecolor='k',linestyle='none',
            markersize=2
            )
    if hi_sd_fit is not None:
        ax.plot(h_sd_fit, hi_sd_fit,
                color = 'r')

    ax.set_xscale(scale)
    ax.set_yscale(scale)

    if limits is not None:
    	ax.set_xlim(limits[0],limits[1])
    	ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
    ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()

def plot_rh2_vs_h(rh2, h_sd, rh2_error=None, h_sd_error=None,
        rh2_fit = None, h_sd_fit = None,
        limits = None, fit = True, savedir = './', filename = None, show = True,
        scale = 'linear', title = ''):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib


    # Drop the NaNs from the images
    if type(rh2_error) is float:
        indices = np.where((rh2 == rh2) &\
                           (h_sd == h_sd)&\
                           (h_sd > 0) &\
                           (rh2 > -5))

    if type(rh2_error) is np.ndarray or \
            type(rh2_error) is np.ma.core.MaskedArray or \
            type(h_sd_error) is np.ndarray or \
            type(h_sd_error) is np.ma.core.MaskedArray:
        indices = np.where((rh2 == rh2) &\
                           (h_sd == h_sd) &\
                           (h_sd_error == h_sd_error) &\
                           (rh2_error == rh2_error) &\
                           (h_sd > 0) &\
                           (rh2 > 0))

    rh2_nonans = rh2[indices]
    h_sd_nonans = h_sd[indices]

    if type(rh2_error) is np.ndarray:
        rh2_error_nonans = rh2_error[indices]
    else:
        rh2_error_nonans = np.array(rh2_error[indices])

    if type(h_sd_error) is np.ndarray or \
            type(h_sd_error) is np.ma.core.MaskedArray:
        h_sd_error_nonans = h_sd_error[indices]
    else:
        h_sd_error_nonans = h_sd_error * \
                np.ones(h_sd[indices].shape)

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
    if 1:
        image = ax.errorbar(h_sd_nonans.ravel(),
                rh2_nonans.ravel(),
                xerr=(h_sd_error_nonans.ravel()),
                yerr=(rh2_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=3
                )
    if rh2_fit is not None:
        ax.plot(h_sd_fit, rh2_fit,
                color = 'r')

    ax.set_xscale(scale)
    ax.set_yscale(scale)

    if limits is not None:
    	ax.set_xlim(limits[0],limits[1])
    	ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
    ax.set_ylabel(r'R$_{H2}$',)
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()

def krumholz_eq(h_sd, phi_CNM = None,
        Z = 1.0, # metallicity
        a = 0.2, # ?
        f_diss = 0.1, # fraction of absorbing H2 which disociates
        phi_mol = 10.0, # molecular gas fraction
        mu_H = 2.3e-24, # molecular weight of H, g
        ):
    '''
    Krumholz et al. 2008
    '''

    # Constants
    c = 3.0e10 # speed of light, cm/s

    # solar values
    sigma_d_solar = 1e-21 # solar dust grain cross section, cm^2
    R_d_solar = 10**-16.5 # solar cloud radius, cm
    E_0_solar = 7.5e-4 # Radiation field, erg/s

    # cloud values
    sigma_d = sigma_d_solar * Z # dust grain cross section, cm^2
    R_d = R_d_solar * Z # cloud radius, cm

    # normalized radiation field strength, EQ 7
    chi = ((f_diss * sigma_d_solar * c * E_0_solar) \
            * (1.0 + (3.1 * Z**0.365))) \
            / (31.0 * phi_CNM * R_d_solar)

    # dust-adjusted radiation field, EQ 10
    psi = chi * (2.5 + chi) / (2.5 + (chi * np.e))

    # cloud optical depth, EQ 21
    tau_c = (3.0 * h_sd * sigma_d) / (4.0 * (3.1 * Z**0.365) * mu_H)

    # cloud optical depth, EQ 21
    tau_c = 0.067 * Z * h_sd

    # ratio of molecular to atomic fraction, EQ 17 Lee et al. 2012
    R_H2 = 4 * tau_c / (3 * psi) \
            * (1+ 0.8 * psi * phi_mol \
                / (4 * tau_c + 3 * (phi_mol - 1) * psi)) -1

    return R_H2

def krumholz_eq_simple(h_sd, phi_CNM):

   return krumholz_eq(h_sd, phi_CNM = phi_CNM,
        Z = 1.0, # metallicity
        a = 0.2, # ?
        f_diss = 0.1, # fraction of absorbing H2 which disociates
        phi_mol = 10.0, # molecular gas fraction
        mu_H = 2.3e-24, # molecular weight of H, g
        )

def fit_krumholz(h_sd, rh2, h_sd_extent, p0 = 10, return_params = False):

    from scipy.optimize import curve_fit

    # fit the krumholz model, choose first element of tuple = parameters
    rh2_fit_params = curve_fit(krumholz_eq_simple, h_sd, rh2,
            maxfev=10000000, p0 = p0)[0]

    # Create large array of h_sd
    h_sd_extended = np.linspace(h_sd_extent[0], h_sd_extent[1], 1e5)
    rh2_fit = krumholz_eq_simple(h_sd_extended, rh2_fit_params)

    if return_params:
        return rh2_fit, h_sd_extended, rh2_fit_params
    else:
        return rh2_fit, h_sd_extended

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

def calculate_noise_cube(cube=None, velocity_axis=None,
            velocity_noise_range=[-110,-90,90,110], header=None, Tsys=30.,
            filename=None):

    """ Calcuates noise envelopes for each pixel in a cube
    """

    import numpy as np
    import pyfits as pf

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

def calculate_sd(image, sd_factor=1/1.25, units = 'galactic'):

    ''' Calculates a surface density image given a velocity range within which
    to include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to
    units : str
        'galactic' = M_sun / pc^-2
        'cgs' = g / cm^-2

    '''

    import numpy as np

    # NHI in units of 1e20 cm^-2
    hi_image = image * sd_factor

    if units == 'cgs':
        hi_image *= 1.981e33 / 3e19

    return hi_image

def calculate_nh2(nhi_image = None, av_image = None, dgr = 5.3e-2):

    ''' Calculates the total gas column density given N(HI), A_v and a
    dust-to-gas ratio.

    Parameters
    ----------
    nhi_image in units of 10^20 cm^-2
    '''

    import numpy as np

    nh_image = 0.5 * (av_image / dgr - nhi_image)

    return nh_image

def calculate_correlation(SpectralGrid=None,cube=None,velocity_axis=None,
        av_image=None, velocity_centers=[], velocity_widths=[],av_noise=0.5,
        av_SNR=None):
    ''' Calculates the correlation coefficient between either a SpectralGrid or
    a cube and an Av image.

    Parameters
    ----------
    cube : array-like, optional

    '''

    # Import external modules
    import grid
    import numpy as np
    from scipy.stats import pearsonr

    # calculate the velocity ranges given a set of centers and widths
    velocity_ranges = np.zeros(shape=[len(velocity_centers) * \
            len(velocity_widths),2])
    count = 0
    for i, center in enumerate(velocity_centers):
        for j, width in enumerate(velocity_widths):
            velocity_ranges[count,0] = center - width/2.
            velocity_ranges[count,1] = center + width/2.
            count += 1

    # calculate the correlation coefficient for each velocity range
    correlations = np.zeros(velocity_ranges.shape[0])
    pvalues = np.zeros(velocity_ranges.shape[0])

    # calc
    for i in range(velocity_ranges.shape[0]):
        nhi_image_temp = calculate_nhi(SpectralGrid=SpectralGrid,
                cube=cube,
                velocity_axis=velocity_axis,
                velocity_range=velocity_ranges[i])
        nhi_image = np.ma.array(nhi_image_temp,
                                mask=np.isnan(nhi_image_temp))

        # Select pixels with Av > 1.0 mag and Av_SNR > 5.0.
        # Av > 1.0 mag is used to avoid too low Av.
        # 1.0 mag corresponds to SNR = 1 / 0.2 ~ 5
        # (see Table 2 of Ridge et al. 2006).
        indices = np.where((nhi_image == nhi_image) & \
                (av_image == av_image) & \
                (av_image < 2))# & \
                #(av_image > 5*av_noise))

        nhi_image_corr = nhi_image[indices]
        av_image_corr = av_image[indices] #- 0.8 # subtract background of 0.8
        # Use Pearson's correlation test to compare images
        correlations[i] = pearsonr(nhi_image_corr.ravel(),
                av_image_corr.ravel())[0]

        # Shows progress each 10%
        total = float(velocity_ranges.shape[0])
        abs_step = int((total * 1)/100) or 1
        if i and not i % abs_step:
                 print "{0:.2%} processed".format(i/total)
    return correlations

def convert_core_coordinates(cores, header):

    for core in cores:
        cores[core].update({'center_pixel': 0})
    	center_wcs = cores[core]['center_wcs']
        # convert centers to pixel coords
        cores[core]['center_pixel'] = get_pix_coords(ra=center_wcs[0],
                                                     dec=center_wcs[1],
                                                     header=header)

        try:
            cores[core].update({'box_pixel': 0})

            box_wcs = cores[core]['box_wcs']
            box_pixel = len(box_wcs) * [0,]

            # convert box corners to pixel coords
            for i in range(len(box_wcs)/2):
                pixels = get_pix_coords(ra=box_wcs[2*i], dec=box_wcs[2*i + 1],
                        header=header)
                box_pixel[2*i], box_pixel[2*i + 1] = int(pixels[0]), \
                    int(pixels[1])
            cores[core]['box_pixel'] = box_pixel
        except KeyError:
            print('WARNING: One or more cores do not have boxes.')

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

def load_ds9_region(cores, filename_base = 'california_av_boxes_', header=None):

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

    import numpy as np
    from os import system,path
    import mygeometry as myg
    reload(myg)


    # define directory locations
    if 1:
        output_dir = '/home/ezbc/Desktop/'
        figure_dir = '/home/ezbc/research/figures/'
        av_dir = '/home/ezbc/research/data/california/'
        hi_dir = '/home/ezbc/research/data/california/'
        core_dir = output_dir + 'core_arrays/'
        region_dir = '/home/ezbc/research/data/california/'
    if 0:
        output_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'
        figure_dir = '/d/bip3/ezbc/california/figures/'
        av_dir = '/d/bip3/ezbc/california/data/av/'
        hi_dir = '/d/bip3/ezbc/california/data/galfa/'
        core_dir = output_dir + 'core_arrays/'
        region_dir = '/d/bip3/ezbc/california/data/python_output/ds9_regions/'

    # load 2mass Av and GALFA HI images, on same grid
    av_data, av_header = load_fits(av_dir + 'california_planck_av_regrid.fits',
            return_header=True)

    #av_data += - 0.4 # subtracts background of 0.4 mags
    hi_data,h = load_fits(hi_dir + 'california_galfa_cube_bin_3.7arcmin.fits',
            return_header=True)

    # make the velocity axis
    velocity_axis = (np.arange(h['NAXIS3']) - h['CRPIX3'] + 1) * h['CDELT3'] + \
            h['CRVAL3']
    velocity_axis /= 1000.

    # Plot NHI vs. Av for a given velocity range
    noise_cube_filename = 'california_galfa_cube_bin_3.7arcmin_noise.fits'
    if not path.isfile(hi_dir + noise_cube_filename):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=velocity_axis,
                velocity_noise_range=[90,110], header=h, Tsys=30.,
                filename=hi_dir + noise_cube_filename)
    else:
        noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    # calculate nhi and error maps, write nhi map to fits file
    nhi_image, nhi_image_error = calculate_nhi(cube=hi_data,
        velocity_axis=velocity_axis,
        noise_cube = noise_cube,
        velocity_range=[-100,100],
        return_nhi_error=True,
        fits_filename = hi_dir + 'california_galfa_nhi_3.7arcmin.fits',
        fits_error_filename = hi_dir + 'california_galfa_nhi_error_3.7arcmin.fits',
        header = h)

    # calculate N(H2) maps
    nh2_image = calculate_nh2(nhi_image = nhi_image,
            av_image = av_data, dgr = 5.3e-2)
    nh2_image_error = calculate_nh2(nhi_image = nhi_image_error,
            av_image = 0.1, dgr = 5.3e-2)

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi_image, sd_factor=1./1.25,)
    hi_sd_image_error = calculate_sd(nhi_image_error, sd_factor=1./1.25)
    h2_sd_image = calculate_sd(nh2_image, sd_factor=10./6.25,)
    h2_sd_image_error = calculate_sd(nh2_image_error, sd_factor=10./6.25)
    h_sd_image = hi_sd_image + h2_sd_image
    h_sd_image_error = (hi_sd_image_error**2 + h2_sd_image_error**2)**0.5

    # Write ratio between H2 and HI
    rh2_image = h2_sd_image / hi_sd_image
    rh2_image_error = rh2_image * (hi_sd_image_error**2 / hi_sd_image**2 \
                 + h2_sd_image_error**2 / h2_sd_image**2)**0.5

    # define core properties
    cores = {'L1482':
                {'center_wcs': [(4, 29, 41), (35, 48, 41)],
                 'map': None,
                 'threshold': None,
                 },
              'L1483':
                {'center_wcs': [(4, 34, 57), (36, 18, 12)],
                 'map': None,
                 'threshold': None,
                 },
              'L1478':
                {'center_wcs': [(4, 25, 7), (37, 13, 0)],
                 'map': None,
                 'threshold': None,
                 },
              'L1434':
                {'center_wcs': [(3, 50, 51), (35, 15, 10)],
                 'map': None,
                 'threshold': None,
                 },
              'L1503':
                {'center_wcs': [(4, 40, 27), (29, 57, 12)],
                 'map': None,
                 'threshold': None,
                 },
              'L1507':
                {'center_wcs': [(4, 42, 51), (29, 44, 47)],
                 'map': None,
                 'threshold': None,
                 },
                }

    cores = convert_core_coordinates(cores, h)

    if True:
        hsd_limits =[0.1,300]
        hisd_limits = [4,50]
        av_limits =[0.1,50]
        nhi_limits = [2,20]

        cores = load_ds9_region(cores,
                filename_base = region_dir + 'california_av_boxes_',
                header = h)

        print('core\t phi_CNM')

        for core in cores:

            av_limits =[0.01,100]

            # Grab the mask
            xy = cores[core]['box_center_pix']
            box_width = cores[core]['box_width']
            box_height = cores[core]['box_height']
            box_angle = cores[core]['box_angle']
            mask = myg.get_rectangular_mask(nhi_image,
        	        xy[0], xy[1],
                    width = box_width,
                    height = box_height,
                    angle = box_angle)
            # apply mask and avoid NaNs
            indices = np.where((mask == 1) &\
                    (hi_sd_image == hi_sd_image) &\
                    (h2_sd_image == h2_sd_image) &\
                    (h_sd_image == h_sd_image))

            # Apply mask to images
            av_data_sub = av_data[indices]
            hi_sd_image_sub = hi_sd_image[indices]
            hi_sd_image_error_sub = hi_sd_image_error[indices]
            h2_sd_image_sub = h2_sd_image[indices]
            h2_sd_image_error_sub = h2_sd_image_error[indices]
            h_sd_image_sub = h_sd_image[indices]
            h_sd_image_error_sub = h_sd_image_error[indices]
            rh2_image_sub = rh2_image[indices]
            rh2_image_error_sub = rh2_image_error[indices]

            # Unravel images to single dimension
            hi_sd_ravel = hi_sd_image_sub.ravel()
            h2_sd_ravel = h2_sd_image_sub.ravel()
            hi_sd_error_ravel = hi_sd_image_error_sub.ravel()
            h2_sd_error_ravel = h2_sd_image_error_sub.ravel()
            h_sd_ravel = h_sd_image_sub.ravel()
            h_sd_error_ravel = h_sd_image_error_sub.ravel()
            rh2_ravel = rh2_image_sub.ravel()
            rh2_error_ravel = rh2_image_error_sub.ravel()

            # write indices for only ratios > 0
            indices = np.where(rh2_ravel > 0)

            # Fit to krumholz model, init guess of phi_CNM = 10
            rh2_fit, h_sd_fit, phi_cnm = fit_krumholz(h_sd_ravel[indices],
                    rh2_ravel[indices],
                    (0.001, 1000, 1e6),
                    p0 = 10,
                    return_params = True)

            print('%s\t %.2f' % (core, phi_cnm))
            # see eq 6 of krumholz+09
            # phi_cnm is the number density of the CNM over the minimum number
            # density required for pressure balance
            # the lower phi_cnm values than for perseus mean that california
            # has a more diffuse CNM

            '''
            rh2 = fh2 / fhi
            rh2 = (2*(1-hi_sd) / (hi_sd+2*h2_sd)) / (hi_sd / (hi_sd+2*h2_sd))
            rh2 = (2*(1-hi_sd) / hi_sd
            rh2 = 2/hi_sd - 2
            rh2 / 2 + 1 = 1 / hi_sd
            hi_sd = 1 / (rh2/2 + 1)

            '''

            hi_sd_fit = h_sd_fit / (rh2_fit/2. + 1)

            if True:
                plot_rh2_vs_h(rh2_ravel, h_sd_ravel,
                        rh2_error = rh2_error_ravel,
                        h_sd_error = h_sd_error_ravel,
                        rh2_fit = rh2_fit,
                        h_sd_fit = h_sd_fit,
                        limits = [0.5, 1000, 10**-3, 10**2],
                        savedir = figure_dir,
                        scale = 'log',
                        filename = 'california_rh2_vs_hsd_' + core + '_box.png',
                        title = r'$R_{\rm H2}$ vs. $\Sigma_{\rm HI}$'\
                                + ' of California Core ' + core,
                        show = False)

                plot_hi_vs_h(hi_sd_ravel, h_sd_ravel,
                        hi_sd_error = hi_sd_error_ravel,
                        h_sd_error = h_sd_error_ravel,
                        hi_sd_fit = hi_sd_fit,
                        h_sd_fit = h_sd_fit,
                        limits = [0.5, 1000, 10**0, 10**2],
                        savedir = figure_dir,
                        scale = 'log',
                        filename = 'california_hisd_vs_hsd_' + core + \
                                '_box.png',
                        title = r'$\Sigma_{\rm HI}$ vs. $\Sigma_{\rm H}$'\
                                + ' of California Core ' + core,
                        show = False)


if __name__ == '__main__':
    main()



