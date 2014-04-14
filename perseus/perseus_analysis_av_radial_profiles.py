#!/usr/bin/python

''' Calculates the Av profiles Taurus molecular cloud.
'''

def plot_profile(radii, profile, limits=None, savedir='./', filename=None,
        show=True, scale='linear', title='', profile_errors=None,
        profile_fit_params=None, profile_fit_function=None):

    ''' Plots N(HI) as a function of Av for individual pixels in an N(HI) image
    and an Av image.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    fontScale = 20
    params = {#'backend': 'png',
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
             }
    plt.rcParams.update(params)

    # Create figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(radii, profile,
            yerr=profile_errors,
            color='k',
            markersize=5,
            marker='s',
            linestyle='None',)
    if profile_fit_params is not None:
    	profile_fit_radii = np.linspace(limits[0], limits[1], 1000)
    	profile_fit_data = profile_fit_function(profile_fit_radii,
    	        *profile_fit_params)
    	ax.plot(profile_fit_radii, profile_fit_data, color='k')

    if limits is not None:
    	ax.set_xlim(limits[0],limits[1])
    	ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_ylabel('A$_V$ (mag)',)
    ax.set_xlabel(r'Radius (pc)',)
    ax.set_title(title)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid(True)



    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()

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

def get_sub_image(image,indices):

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
            stddev=stddev, mask=mask, interpnan=True, returnradii=True,
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
            maxfev=10000)

    return profile_fit

def print_fit_params(cores, A_p, pho_c, R_flat, p, filename=None):

    ''' Prints statistics about each core.
    '''

    import os

    #print('A_p\t pho_c\t R_flat\t p ')
    #print(core_name + ':')
    #print('%.2f\t %.2f \t %.2f \t %.2f \n' % \
    #        (A_p, pho_c, R_flat, p))

    if filename is not None:
        os.system('rm -rf ' + filename)
        f = open(filename, 'w')
        f.write('Core\t A_p\t pho_c\t R_flat\t p \n')
        for i, core in enumerate(cores):
        	f.write('%s\t %.2f\t %.2f \t %.2f \t %.2f \n' % \
                    (core, A_p[i], pho_c[i], R_flat[i], p[i]))
        f.close()

def main():

    import grid
    import numpy as np
    from os import system,path
    import myclumpfinder as clump_finder

    # define directory locations
    output_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/perseus/figures/'
    av_dir = '/d/bip3/ezbc/perseus/data/2mass/'
    hi_dir = '/d/bip3/ezbc/perseus/data/galfa/'
    core_dir = output_dir + 'core_arrays/'

    # load 2mass Av and GALFA HI images, on same grid
    av_image, h = load_fits(av_dir + '2mass_av_lee12_nocal.fits',
            return_header=True)

    cores = {'IC348':
                {'center_wcs': [(3,44,0), (32, 8, 0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(3,46,13), (26,03,24),(3,43,4),(32,25,41)]
                 },
              'B5':
                {'center_wcs': [(03, 47, 53), (32, 50, 8)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(3,58,46),(32,13,24), (3,45,4), (34,50,41),]
                 },
               'NGC 1333':
                {'center_wcs': [(03, 28, 57), (31, 22, 47)],
                 'map': None,
                 'threshold': 4,
                 'box_wcs': [(3,30,54), (30,58,1), (3,10,50), (33,49,45)]
                 },
                'B1':
                {'center_wcs': [(3,36,21), (31,13,4)],
                 'map': None,
                 'threshold': 4,
                 'box_wcs': [(3,38,52), (27,3,2), (3,34,48), (31,38,8)]}
                }

    cores = convert_core_coordinates(cores, h)

    if True:
        limits = [0.1, 10, 1, 10]

        A_p = []
        pho_c = []
        R_flat = []
        p = []

        for core in cores:
            print('Calculating for core %s' % core)

            # plotting radial profiles only within boxes
            #av_image_sub = get_sub_image(av_image, cores[core]['box_pixel'])
            # change center pixel to correspond to sub image
            #pix = cores[core]['center_pixel']
            #pix = (pix[0] - cores[core]['box_pixel'][0],
            #    pix[1] - cores[core]['box_pixel'][1],)
            #cores[core]['center_pixel'] = pix

            av_image_sub = av_image

            # extract radial profile weighted by SNR
            radii, profile = get_radial_profile(av_image_sub, binsize=3,
                    center=cores[core]['center_pixel'],
                    weights=av_image_sub / 0.3)
            # extract std
            radii, profile_std = get_radial_profile(av_image_sub, binsize=3,
                    center=cores[core]['center_pixel'], stddev=True,
                    weights=av_image_sub / 0.3)

            # convert radii from degrees to parsecs
            radii_arcmin = radii * h['CD2_2'] * 60 * 60. # radii in arcminutes
            radii_pc = radii_arcmin * 300 / 206265. # radii in parsecs

            # extract radii from within the limits
            indices = np.where(radii_pc < limits[1])
            radii_pc = radii_pc[indices]
            profile = profile[indices]
            profile_std = profile_std[indices]

            # fit profile with power function
            def function(radius, A_p, pho_c, R_flat, p):
                return A_p * pho_c * R_flat / \
                        (1 + (radius / R_flat)**2)**(p/2. - 0.5)
                #return A_p * radius**p

            profile_fit_params = fit_profile(radii_pc, profile, function,
                    sigma=profile / profile_std)[0]

            # plot the radial profile
            plot_profile(radii_pc, profile,
                    profile_errors = profile_std,
                    limits = limits,
                    profile_fit_params = profile_fit_params,
                    profile_fit_function = function,
                    savedir=figure_dir,
                    filename = 'perseus_profile_av_' + core + '.png',
                    title=r'Radial A$_V$ Profile of Taurus Core ' + core,
                    show = False)

            A_p.append(profile_fit_params[0])
            pho_c.append(profile_fit_params[1])
            R_flat.append(profile_fit_params[2])
            p.append(profile_fit_params[3])

        print_fit_params(cores, A_p, pho_c, R_flat, p,
                filename=output_dir + 'core_profile_fit_data.txt')

if __name__ == '__main__':
    main()



