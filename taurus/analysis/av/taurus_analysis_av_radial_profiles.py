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

def plot_profile_grid(radii_list, profile_list, limits=None, savedir='./',
        filename=None, show=True, scale=('linear', 'linear'), title='',
        profile_errors_list=None, profile_fit_params_list=None,
        profile_fit_function=None, core_names=None):

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
    fontScale = 12
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
    n = int(np.ceil(len(radii_list)**0.5))
    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(n, n),
                 ngrids=len(radii_list),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    for i in xrange(len(radii_list)):
        radii = radii_list[i]
        profile = profile_list[i]
        profile_errors = profile_errors_list[i]
        profile_fit_params = profile_fit_params_list[i]

    	ax = imagegrid[i]
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
        ax.set_title(core_names[i])
        ax.set_yscale(scale[0], nonposy = 'clip')
        ax.set_xscale(scale[1], nonposx = 'clip')
        ax.grid(True)

    if title:
        fig.suptitle(title, fontsize=1.5*fontScale)
    if filename is not None:
        plt.savefig(savedir + filename)#, bbox_inches='tight')
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
    av_image, h = load_fits(av_dir + 'taurus_av_planck_5arcmin.fits',
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
        limits = [0, 20, -1, 25] # x-linear limits

        # Initialize fit params
        A_p = []
        pho_c = []
        R_flat = []
        p = []

        # Initialize data lists
        radii_pc_list = []
        profile_list = []
        profile_std_list = []
        profile_fit_params_list = []
        core_names_list = []

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

            av_image_sub = np.copy(av_image)
            #av_image_sub[mask == 0] = np.NaN
            av_image_sub = np.ma.array(av_image, mask=(mask == 0))

            # to check the positions of the boxes, uncomment the following
            #import matplotlib.pyplot as plt
            #plt.clf()
            #plt.imshow(np.ma.array(av_image_sub, mask=temp_mask))
            #plt.savefig('/usr/users/ezbc/Desktop/map%s.png' % core)
            #plt.clf()

            pix = cores[core]['center_pixel']

            # extract radial profile weighted by SNR
            radii, profile = get_radial_profile(av_image, binsize=3,
                    center=pix,
                    weights=av_image / 0.3,
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
                        savedir=figure_dir + 'individual_cores/',
                        filename = 'taurus_profile_av_' + core + figure_type,
                        title=r'Radial A$_V$ Profile of Taurus Core ' + core,
                        show = False)

            A_p.append(profile_fit_params[0])
            pho_c.append(profile_fit_params[1])
            R_flat.append(profile_fit_params[2])
            p.append(profile_fit_params[3])

            radii_pc_list.append(radii_pc)
            profile_list.append(profile)
            profile_std_list.append(profile_std)
            profile_fit_params_list.append(profile_fit_params)
            core_names_list.append(core)

        for figure_type in figure_types:
            plot_profile_grid(radii_pc_list, profile_list,
                    profile_errors_list = profile_std_list,
                    limits = limits,
                    profile_fit_params_list = profile_fit_params_list,
                    profile_fit_function = function,
                    savedir=figure_dir + 'panel_cores/',
                    filename = 'taurus_profile_av_cores_planck' + figure_type,
                    title=r'Radial A$_V$ Profiles of Taurus Cores',
                    core_names=core_names_list,
                    show = False)


        print_fit_params(cores, A_p, pho_c, R_flat, p,
                filename=output_dir + 'core_profile_fit_data.txt')

        print_fit_params(cores, A_p, pho_c, R_flat, p)

if __name__ == '__main__':
    main()



