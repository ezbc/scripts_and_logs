#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the Taurus molecular cloud.
'''

def plot_correlations(correlations,velocity_centers,velocity_widths,
        savedir='./', filename=None,show=True, returnimage=False):
    ''' Plots a heat map of correlation values as a function of velocity width
    and velocity center.
    '''

    # Import external modules
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8,2.5))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="0.3%",
                 cbar_size='6%',
                 axes_pad=0,
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    correlations_image = np.empty((velocity_centers.shape[0],
                                   velocity_widths.shape[0]))
    correlations_image[:,:] = np.NaN
    count = 0
    try:
        for i, center in enumerate(velocity_centers):
            for j, width in enumerate(velocity_widths):
                correlations_image[i,j] = correlations[count]
                count += 1
    except IndexError:
        print(' plot_correlations: O-d array input, cannot proceed')

    image = np.ma.array(correlations_image,mask=np.isnan(correlations_image))

    ax = imagegrid[0]
    ax.set_xlabel('Velocity Width (km/s)',
              size = 'small',
              family='serif')
    ax.set_ylabel('Velocity Center (km/s)',
              size = 'small',
              family='serif')

    #ax.set_xticks(np.arange(0,velocity_widths.shape[0],1)[::5],
    #        velocity_centers[::5])

    plt.rc('text', usetex=True)
    im = ax.imshow(image, interpolation='nearest',origin='lower',
            extent=[velocity_widths[0],velocity_widths[-1],
                velocity_centers[0],velocity_centers[-1]],
            cmap=plt.jet())
    cb = ax.cax.colorbar(im)

    #cb.set_clim(vmin=0.)
    # Write label to colorbar
    cb.set_label_text(r'Correlation coefficient',
                   size='small',
                   family='serif')

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_center_velocity(correlations, velocity_centers, velocity_widths,
        velocity_center=None, savedir='./', filename=None, show=True,
        returnimage=False):
    ''' Plots correlation vs. velocity width for a single center velocity.
    '''

    # Import external modules
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8,8))

    correlations_image = np.empty((velocity_centers.shape[0],
                                   velocity_widths.shape[0]))
    correlations_image[:,:] = np.NaN
    count = 0
    for i, center in enumerate(velocity_centers):
        for j, width in enumerate(velocity_widths):
            correlations_image[i,j] = correlations[count]
            count += 1

    ax = fig.add_subplot(111)
    ax.set_xlabel('Velocity Width (km/s)',
              size = 'small',
              family='serif')
    ax.set_ylabel('Correlation Coefficient',
              size = 'small',
              family='serif')
    ax.set_title('Center velocity ' + str(velocity_center) + ' km/s')

    center_index = np.where(velocity_centers == velocity_center)[0][0]

    ax.plot(velocity_widths,
            correlations_image[center_index,:])
    # fig.show()

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_nhi_vs_av(nhi_image, av_image,
        nhi_image_error=None, av_image_error=None, limits=None,
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
    indices = np.where((nhi_image == nhi_image) &\
                       (av_image == av_image)&\
                       (av_image > 0) &\
                       (nhi_image > -5))

    nhi_image_nonans = nhi_image[indices]
    av_image_nonans = av_image[indices]

    if type(nhi_image_error) is np.ndarray:
        nhi_image_error_nonans = nhi_image_error[indices]
    else:
        nhi_image_error_nonans = np.array(nhi_image_error[indices])

    if type(av_image_error) is np.ndarray:
        av_image_error_nonans = av_image_error[indices]
    else:
        av_image_error_nonans = av_image_error * \
                np.ones(av_image[indices].shape)

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    if nhi_image_error is None:
        if plot_type is 'hexbin':
            image = ax.hexbin(av_image_nonans.ravel(),
                    nhi_image_nonans.ravel(),
                    norm=matplotlib.colors.LogNorm(),
                    mincnt=1,
                    yscale='log',
                    xscale='log')
            # Adjust color bar of density plot
            cb = plt.colorbar(image)
            cb.set_label('Counts')
        elif plot_type is 'scatter':
            image = ax.scatter(av_image_nonans.ravel(),
                    nhi_image_nonans.ravel(),
                    alpha=0.3,
                    color='k'
                    )
            ax.set_xscale('log')
            ax.set_yscale('log')
    else:
        image = ax.errorbar(av_image_nonans.ravel(),
                nhi_image_nonans.ravel(),
                xerr=(av_image_error_nonans.ravel()),
                yerr=(nhi_image_error_nonans.ravel()),
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
    ax.set_xlabel('A$_v$ (mag)',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'N(HI) (1 $\times 10^{20}$ cm$^{-2}$)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_hisd_vs_hsd(hi_sd_image, h_sd_image,
        hi_sd_image_error=None, h_sd_image_error=None, limits=None,
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
    indices = np.where((hi_sd_image == hi_sd_image) &\
                       (h_sd_image == h_sd_image)&\
                       (h_sd_image > 0) &\
                       (hi_sd_image > -5))

    hi_sd_image_nonans = hi_sd_image[indices]
    h_sd_image_nonans = h_sd_image[indices]

    if type(hi_sd_image_error) is np.ndarray:
        hi_sd_image_error_nonans = hi_sd_image_error[indices]
    else:
        hi_sd_image_error_nonans = np.array(hi_sd_image_error[indices])

    if type(h_sd_image_error) is np.ndarray:
        h_sd_image_error_nonans = h_sd_image_error[indices]
    elif type(h_sd_image_error) is np.ma.core.MaskedArray:
        #h_sd_image_error_nonans = np.copy(h_sd_image_error[indices])
        h_sd_image_error_nonans = h_sd_image_error[indices]
        print 'nope'
    else:
        h_sd_image_error_nonans = h_sd_image_error * \
                np.ones(h_sd_image[indices].shape)

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    if hi_sd_image_error is None:
        if plot_type is 'hexbin':
            image = ax.hexbin(h_sd_image_nonans.ravel(),
                    hi_sd_image_nonans.ravel(),
                    norm=matplotlib.colors.LogNorm(),
                    mincnt=1,
                    yscale='log',
                    xscale='log')
            # Adjust color bar of density plot
            cb = plt.colorbar(image)
            cb.set_label('Counts')
        elif plot_type is 'scatter':
            image = ax.scatter(h_sd_image_nonans.ravel(),
                    hi_sd_image_nonans.ravel(),
                    alpha=0.3,
                    color='k'
                    )
            ax.set_xscale('log')
            ax.set_yscale('log')
    else:
        image = ax.errorbar(h_sd_image_nonans.ravel(),
                hi_sd_image_nonans.ravel(),
                xerr=(h_sd_image_error_nonans.ravel()),
                yerr=(hi_sd_image_error_nonans.ravel()),
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
    ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{H2}$ (M$_\odot$ / pc$^2$)',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

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

    # Create figure
    fig = plt.figure(figsize=(8,8))
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
    ax.set_xlabel('A$_v$ (mag)',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def calculate_NHI(cube=None, velocity_axis=None, SpectralGrid=None,
        velocity_range=[], return_nhi_error=True, noise_cube=None,
        velocity_noise_range=[-110,-90,90,100],
        Tsys=30.):
    ''' Calculates an N(HI) image given a velocity range within which to
    include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to '''

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
        nhi_image_temp = calculate_NHI(SpectralGrid=SpectralGrid,
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

def hrs2degs(ra=None, dec=None):
    ''' Ra and dec tuples in hrs min sec and deg arcmin arcsec.
    '''

    ra_deg = 15*(ra[0] + ra[1]/60. + ra[2]/3600.)
    dec_deg = dec[0] + dec[1]/60. + dec[2]/3600.

    return (ra_deg, dec_deg)

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in degrees.
    '''

    import pywcsgrid2 as wcs
    import pywcs

    ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)

    wcs_header = pywcs.WCS(header)
    pix_coords = wcs_header.wcs_sky2pix([[ra_deg, dec_deg, 0]], 0)[0]

    return pix_coords

def main():

    import grid
    import numpy as np
    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    core_dir = output_dir + 'core_arrays/'

    # load 2mass Av and GALFA HI images, on same grid
    av_data, av_header = load_fits(av_dir + 'taurus_av_k09_regrid.fits',
            return_header=True)
    # load Av image from goldsmith: Pineda et al. 2010, ApJ, 721, 686
    av_data_goldsmith = load_fits(av_dir + \
            'taurus_av_p10_regrid.fits')

    #av_data += - 0.4 # subtracts background of 0.4 mags
    hi_data,h = load_fits(hi_dir + 'taurus_galfa_cube_bin_3.7arcmin.fits',
            return_header=True)

    # make the velocity axis
    velocity_axis = (np.arange(h['NAXIS3']) - h['CRPIX3'] + 1) * h['CDELT3'] + \
            h['CRVAL3']
    velocity_axis /= 1000.

    # Plot NHI vs. Av for a given velocity range
    noise_cube_filename = 'taurus_galfa_cube_bin_3.7arcmin_noise.fits'
    if not path.isfile(hi_dir + noise_cube_filename):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=velocity_axis,
                velocity_noise_range=[90,110], header=h, Tsys=30.,
                filename=hi_dir + noise_cube_filename)
    else:
        noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    # calculate maps
    nhi_image, nhi_image_error = calculate_NHI(cube=hi_data,
        velocity_axis=velocity_axis, noise_cube = noise_cube,
        velocity_range=[-100,100], return_nhi_error=True)

    nh2_image = calculate_nh2(nhi_image = nhi_image,
            av_image = av_data, dgr = 1.1e-1)
    nh2_image_error = calculate_nh2(nhi_image = nhi_image_error,
            av_image = 0.1, dgr = 1.1e-1)

    hi_sd_image = calculate_sd(nhi_image, sd_factor=1/1.25)
    hi_sd_image_error = calculate_sd(nhi_image_error, sd_factor=1/1.25)

    h2_sd_image = calculate_sd(nh2_image, sd_factor=1/6.25)
    h2_sd_image_error = calculate_sd(nh2_image_error, sd_factor=1/6.25)

    h_sd_image = av_data / (1.25 * 0.11) # DGR = 1.1e-12 mag / cm^-2
    h_sd_image_error = 0.1 / (1.25 * 0.11)


    cores = {'L1495':
                {'wcs_position': [15*(4+14/60.), 28+11/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [get_pix_coords(ra=(4,16,30.031),
                                        dec=(27,44,30),
                                        header=h),
                         get_pix_coords(ra=(4,5,20),
                                        dec=(28,28,33),
                                        header=h)]
                 },
             'L1495A':
                {'wcs_position': [15*(4+18/60.), 28+23/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [get_pix_coords(ra=(4,28,23),
                                        dec=(28,12,50),
                                        header=h),
                         get_pix_coords(ra=(4,16,23),
                                        dec=(29,46,5),
                                        header=h)],
                 },
             'B213':
                {'wcs_position': [15*(4+19/60.), 27+15/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [get_pix_coords(ra=(4,22,27),
                                        dec=(26,45,47),
                                        header=h),
                         get_pix_coords(ra=(4,5,25),
                                        dec=(27,18,48),
                                        header=h)],
                },
             'B220':
                {'wcs_position': [15*(4+41/60.), 26+7/60., 0],
                 'map': None,
                 'threshold': 7,
                 'box': [get_pix_coords(ra=(4,47,49),
                                        dec=(25,31,13),
                                        header=h),
                         get_pix_coords(ra=(4,40,37),
                                        dec=(27,31,17),
                                        header=h)],
                 },
             'L1527':
                {'wcs_position': [15*(4+39/60.), 25+47/60., 0],
                 'map': None,
                 'threshold': 7,
                 'box': [get_pix_coords(ra=(4,40,13),
                                        dec=(24,46,38),
                                        header=h),
                         get_pix_coords(ra=(4,34,35),
                                        dec=(25,56,7),
                                        header=h)],
                 },
             'B215':
                {'wcs_position': [15*(4+23/60.), 25+3/60., 0],
                 'map': None,
                 'threshold': 3,
                 'box': [get_pix_coords(ra=(4,24,51),
                                        dec=(22,36,7),
                                        header=h),
                         get_pix_coords(ra=(4,20,54),
                                        dec=(25,26,31),
                                        header=h)],
                 },
             'L1524':
                {'wcs_position': [15*(4+29/60.), 24+31/60., 0],
                 'map': None,
                 'threshold': 3,
                 'box': [get_pix_coords(ra=(4,31,0),
                                        dec=(22,4,6),
                                        header=h),
                         get_pix_coords(ra=(4,25,33),
                                        dec=(25,0,55),
                                        header=h)],
                 }
                }

    for core in cores:
    	box = cores[core]['box']
    	cores[core]['box'] = (int(box[0][0]),int(box[0][1]),
    	        int(box[1][0]),int(box[1][1]))

    # calculate correlation for cores
    if False:
        limits =[1,22,6,16]

        for core in cores:

            print('Calculating N(HI) vs. Av for core %s' % core)

            core_map = np.load(core_dir + core + '.npy')
            plot_nhi_vs_av(nhi_image,core_map,
                    nhi_image_error = nhi_image_error,
                    av_image_error = 0.1,
                    limits = limits,
                    savedir=figure_dir,
                    plot_type='scatter',
                    filename='taurus_nhi_vs_av_' + core + '_small.png',
                    title=r'N(HI) vs. A$_v$ of Taurus Core ' + core,
                    show=False)

            plot_sd_vs_av(sd_image, core_map,
                    sd_image_error = sd_image_error,
                    av_image_error = 0.1,
                    limits = limits,
                    savedir=figure_dir,
                    plot_type='scatter',
                    filename='taurus_sd_vs_av_' + core + '_small.png',
                    title=r'$\Sigma_{HI}$ vs. A$_v$ of Taurus Core ' + core,
                    show = False)

    if True:
        hsd_limits =[0.1,300]
        hisd_limits = [2,20]
        av_limits =[0.01,100]
        nhi_limits = [2,20]

        for core in cores:
            print('Calculating for core %s' % core)

            indices = cores[core]['box']
            nhi_image_sub = get_sub_image(nhi_image, indices)
            nhi_image_error_sub = get_sub_image(nhi_image_error, indices)
            av_data_sub = get_sub_image(av_data, indices)
            hi_sd_image_sub = get_sub_image(hi_sd_image, indices)
            hi_sd_image_error_sub = get_sub_image(hi_sd_image_error, indices)
            av_limits =[0.01,100]

            if core is 'B213':
            	import mygeometry as myg
            	reload(myg)
            	angle = 40.
                xy = (115, 223)
                width = 20
                height = 40

                # Grab the mask
                indices = myg.get_rectangular_mask(nhi_image, xy[0], xy[1],
                        width=width, height=height, angle=angle,
                        return_indices=True)

                av_data_sub = av_data[indices]
                hi_sd_image_sub = hi_sd_image[indices]
                hi_sd_image_error_sub = hi_sd_image_error[indices]

            if False:
                plot_nhi_vs_av(nhi_image_sub,av_data_sub,
                        nhi_image_error = nhi_image_error_sub,
                        av_image_error = 0.1,
                        limits = [0.01,100,2,20],
                        savedir=figure_dir,
                        plot_type='scatter',
                        scale='log',
                        filename='taurus_nhi_vs_av_' + core + '_box.png',
                        title=r'N(HI) vs. A$_v$ of Taurus Core ' + core,
                        show=False)


            if True:
                plot_sd_vs_av(hi_sd_image_sub, av_data_sub,
                        sd_image_error = hi_sd_image_error_sub,
                        av_image_error = 0.1,
                        limits = [0.01,100,2,20],
                        savedir=figure_dir,
                        plot_type='scatter',
                        scale='log',
                        filename='taurus_sd_vs_av_' + core + '_box.png',
                        title=r'$\Sigma_{HI}$ vs. A$_v$ of Taurus Core ' + core,
                        show=False)

            if False:
                h_sd_image_sub = get_sub_image(h_sd_image, indices)
                #h_sd_image_error_sub = get_sub_image(h_sd_image_error, indices)
                plot_hisd_vs_hsd(hi_sd_image_sub, h_sd_image_sub,
                        h_sd_image_error = h_sd_image_error,
                        hi_sd_image_error = hi_sd_image_error_sub,
                        limits = [0.1,300,2,20],
                        savedir=figure_dir,
                        plot_type='scatter',
                        scale='log',
                        filename='taurus_hisd_vs_hsd_' + core + '_box.png',
                        title=r'$\Sigma_{HI}$ vs. $\Sigma_{HI}$ + ' + \
                                '$\Sigma_{H2}$ of Taurus Core ' + core,
                        show=False)

if __name__ == '__main__':
    main()



