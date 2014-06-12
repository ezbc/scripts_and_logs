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
    ax.set_yscale('linear', nonposy = 'clip')
    ax.set_xscale('linear', nonposx = 'clip')
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
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
            stddev=stddev, mask=mask, interpnan=False, returnradii=True,
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
            maxfev=1000000,)

    return profile_fit

def print_fit_params(cores, A_p, pho_c, R_flat, p, filename=None):

    ''' Prints statistics about each core.
    '''

    import os

    if filename is None:
        for i, core in enumerate(cores):
            print('A_p\t pho_c\t R_flat\t p ')
            print(core + ':')
            print('%.2f\t %.2f \t %.2f \t %.2f \n' % \
                    (A_p[i], pho_c[i], R_flat[i], p[i]))

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
    import mygeometry as myg

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/cores/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'
    core_dir = output_dir + 'core_arrays/'

    # load 2mass Av and GALFA HI images, on same grid
    av_image, h = load_fits(av_dir + 'taurus_planck_av_regrid.fits',
            return_header=True)


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

    cores = convert_core_coordinates(cores, h)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'taurus_av_boxes_',
            header = h)

    if True:
        #limits = [0.1, 100, 0.1, 100] # x-log limits
        limits = [0, 20, -1, 30] # x-linear limits

        A_p = []
        pho_c = []
        R_flat = []
        p = []

        for core in cores:
            print('Calculating for core %s' % core)

            # Grab the mask from the DS9 regions
            xy = cores[core]['box_center_pix']
            box_width = cores[core]['box_width']
            box_height = cores[core]['box_height']
            box_angle = cores[core]['box_angle']
            mask = myg.get_rectangular_mask(av_image,
        	        xy[0], xy[1],
                    width = box_width,
                    height = box_height,
                    angle = box_angle)

            # Get indices where there is no mask, and extract those pixels
            indices = np.where(mask == 1)

            # plotting radial profiles only within boxes
            #av_image_sub = get_sub_image(av_image, cores[core]['box_pixel'])
            # change center pixel to correspond to sub image
            #pix = cores[core]['center_pixel']
            #pix = (pix[0] - cores[core]['box_pixel'][0],
            #    pix[1] - cores[core]['box_pixel'][1],)
            #cores[core]['center_pixel'] = pix

            av_image_sub = av_image
            #av_image_sub[mask == 0] = np.NaN
            temp_mask = np.copy(mask)
            temp_mask[mask == 0] = 1
            temp_mask[mask == 1] = 0

            import matplotlib.pyplot as plt

            plt.clf()
            plt.imshow(np.ma.array(av_image_sub, mask=temp_mask))
            plt.savefig('/usr/users/ezbc/Desktop/map%s.png' % core)
            plt.clf()

            pix = cores[core]['center_pixel']

            # extract radial profile weighted by SNR
            radii, profile = get_radial_profile(av_image_sub, binsize=3,
                    center=pix,
                    weights=av_image_sub / 0.3,
                    mask=mask
                    )

            # extract std
            radii, profile_std = get_radial_profile(av_image_sub, binsize=3,
                    center=pix,
                    stddev=True,
                    weights=av_image_sub / 0.3,
                    #mask=mask
                    )

            # convert radii from degrees to parsecs
            radii_arcmin = radii * h['CDELT2'] * 60 * 60. # radii in arcminutes
            radii_pc = radii_arcmin * 300 / 206265. # radii in parsecs

            # extract radii from within the limits
            indices = np.where((radii_pc < limits[1]) & \
                               (profile == profile) & \
                               (profile_std == profile_std))
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
            figure_types = ['.pdf', '.png']
            for figure_type in figure_types:
                plot_profile(radii_pc, profile,
                        profile_errors = profile_std,
                        limits = limits,
                        profile_fit_params = profile_fit_params,
                        profile_fit_function = function,
                        savedir=figure_dir,
                        filename = 'taurus_profile_av_' + core + figure_type,
                        title=r'Radial A$_V$ Profile of Taurus Core ' + core,
                        show = False)

            A_p.append(profile_fit_params[0])
            pho_c.append(profile_fit_params[1])
            R_flat.append(profile_fit_params[2])
            p.append(profile_fit_params[3])

        print_fit_params(cores, A_p, pho_c, R_flat, p,
                filename=output_dir + 'core_profile_fit_data.txt')

        print_fit_params(cores, A_p, pho_c, R_flat, p)

if __name__ == '__main__':
    main()



