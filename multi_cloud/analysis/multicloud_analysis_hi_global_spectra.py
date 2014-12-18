#!/usr/bin/python

''' Calculates the N(HI) map for Taurus

'''

from astropy.io import fits
import numpy as np

def calc_global_spectrum(hi_cube=None, statistic='average'):

    hi_cube_nonans = np.ma.array(hi_cube, mask=(hi_cube != hi_cube))

    n_pix = hi_cube.shape[1] * hi_cube.shape[2]

    if statistic == 'average':
        spectrum = np.sum(hi_cube_nonans, axis=(1,2)) / n_pix
    elif statistic == 'median':
        spectrum = np.zeros(hi_cube.shape[0])
        for i in xrange(len(spectrum)):
            cube = hi_cube_nonans[i, :, :].ravel()
            spectrum[i] = np.median(cube[~np.isnan(cube)])
    elif statistic == 'std':
        spectrum = np.zeros(hi_cube.shape[0])
        for i in xrange(len(spectrum)):
            spectrum[i] = np.std(hi_cube_nonans[i, :, :].ravel())

    return spectrum

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

def plot_spectra_grid(hi_velocity_axis_list=None, hi_spectrum_list=None,
        co_spectrum_list=None, co_velocity_axis_list=None, hi_std_list=None,
        title=None, limits=None, savedir='./', filename=None, show=True,
        spectra_names='',):

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
    #colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    fontScale = 15
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
              'figure.figsize': (15, 7),
              'figure.titlesize': fontScale
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # Determine number of plots on each axis
    n = int(np.ceil(len(hi_velocity_axis_list)**0.5))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(n, n),
                 ngrids=len(hi_velocity_axis_list),
                 axes_pad=0,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    for i in xrange(len(hi_velocity_axis_list)):

        # create axes
        ax = imagegrid[i]

        ax.plot(hi_velocity_axis_list[i],
                hi_spectrum_list[i],
                color='k',
                linestyle='--',
                label='median HI',
                drawstyle = 'steps-mid'
                )

        if hi_std_list is not None:
            ax.plot(hi_velocity_axis_list[i],
                    hi_std_list[i],
                    color='k',
                    linestyle='-',
                    label=r'$\sigma_{\rm HI}$',
                    drawstyle='steps-mid'
                    )

        if co_velocity_axis_list is not None:
            ax.plot(co_velocity_axis_list[i],
                    co_spectrum_list[i] * 10.,
                    color='r',
                    label=r'median $^{12}$CO X 10',
                    drawstyle = 'steps-mid'
                    )

        ax.annotate(spectra_names[i].capitalize(),
                xytext=(0.1, 0.9),
                xy=(0.1, 0.9),
                textcoords='axes fraction',
                xycoords='axes fraction',
                color='k'
                )

        ax.set_xlabel('Velocity (km/s)')
        ax.set_ylabel(r'T$_b$ (K)')

        # plot limits
        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

    if co_velocity_axis_list is not None:
        # Single legend
        ax.legend(bbox_to_anchor=(1.5, 0.2),
                loc='lower right',
                borderaxespad=0.)

    if title is not None:
        fig.suptitle(title, fontsize=fontScale)
    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        fig.show()

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
    ''' Executes script.
    '''
    import numpy as np
    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)
    import mygeometry as myg
    reload(myg)
    from mycoords import make_velocity_axis


    cloud_list = ('taurus', 'perseus', 'california', 'perseus')
    data_dict = {}

    for i, cloud in enumerate(cloud_list):
        # define directory locations
        output_dir = '/d/bip3/ezbc/%s/data/python_output/nhi_av/' % cloud
        figure_dir = '/d/bip3/ezbc/multicloud/figures/'
        av_dir = '/d/bip3/ezbc/%s/data/av/' % cloud
        hi_dir = '/d/bip3/ezbc/%s/data/hi/' % cloud
        co_dir = '/d/bip3/ezbc/%s/data/co/' % cloud
        core_dir = output_dir + 'core_arrays/'
        region_dir = '/d/bip3/ezbc/%s/data/python_output/ds9_regions/' % cloud

        cloud_dict = {}

        # Load HI maps from Taurus, California, and Perseus
        cloud_dict['hi_data'], cloud_dict['hi_header'] = fits.getdata(
                '/d/bip3/ezbc/{0}/data/hi/{0}_hi_galfa_'.format(cloud) + \
                        'cube_regrid_planckres.fits',
                header=True)

        # Load CO maps from Taurus, California, and Perseus
        cloud_dict['co_data'], cloud_dict['co_header'] = fits.getdata(
                    co_dir + '{0}_co_cfa_'.format(cloud) + \
                        'cube_regrid_planckres.fits',
                header=True)

        # sum along spatial axes
        cloud_dict['hi_spectrum'] = calc_global_spectrum(
                                        hi_cube=cloud_dict['hi_data'],
                                        statistic='median'
                                        )
        cloud_dict['hi_std'] = calc_global_spectrum(
                                        hi_cube=cloud_dict['hi_data'],
                                        statistic='std'
                                        )
        cloud_dict['co_spectrum'] = calc_global_spectrum(
                                        hi_cube=cloud_dict['co_data'],
                                        statistic='median'
                                        )

        # Calculate velocity
        cloud_dict['hi_velocity_axis'] = make_velocity_axis(
                                        cloud_dict['hi_header'])
        cloud_dict['co_velocity_axis'] = make_velocity_axis(
                                        cloud_dict['co_header'])

        data_dict[cloud] = cloud_dict

    hi_velocity_axis_list = []
    co_velocity_axis_list = []
    hi_spectrum_list = []
    hi_std_list = []
    co_spectrum_list = []
    spectra_names = []

    for cloud in data_dict:
        hi_velocity_axis_list.append(data_dict[cloud]['hi_velocity_axis'])
        hi_spectrum_list.append(data_dict[cloud]['hi_spectrum'])
        hi_std_list.append(data_dict[cloud]['hi_std'])
        co_velocity_axis_list.append(data_dict[cloud]['co_velocity_axis'])
        co_spectrum_list.append(data_dict[cloud]['co_spectrum'])
        spectra_names.append(cloud)

    figure_types = ('png',)
    for figure_type in figure_types:
        limits = [-75, 40, -1, 55]

        # scatter plot
        plot_spectra_grid(hi_velocity_axis_list=hi_velocity_axis_list,
                        hi_spectrum_list=hi_spectrum_list,
                        co_velocity_axis_list=co_velocity_axis_list,
                        co_spectrum_list=co_spectrum_list,
                        hi_std_list=hi_std_list,
                        savedir=figure_dir,
                        limits=limits,
                        filename='multicloud_hi_spectrum_global.%s' % \
                                figure_type,
                        show=False,
                        spectra_names=spectra_names,
                        title=r'Global HI Profile' + \
                                ' of GMCs')


if __name__ == '__main__':
    main()


