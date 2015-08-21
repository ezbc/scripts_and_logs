#!/usr/bin/python

''' Calculates the N(HI) map for Taurus

'''

from astropy.io import fits
import pyfits as fits
import numpy as np

''' Plots
'''

def plot_spectra_grid(hi_velocity_axis_list=None, hi_spectrum_list=None,
        hi_vel_range_list=None, co_spectrum_list=None,
        co_velocity_axis_list=None, hi_std_list=None, title=None, limits=None,
        savedir='./', filename=None, show=True, spectra_names='',):

    # Import external modules
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Create figure instance
    fig = plt.figure(figsize=(3.6, 3.6))
    ax = fig.add_subplot(111)

    myplt.set_color_cycle(num_colors=3)

    for i in xrange(len(hi_velocity_axis_list)):

        # create axes
        spectra_names[i] = spectra_names[i].replace('1', ' 1')
        spectra_names[i] = spectra_names[i].replace('2', ' 2')

        ax.plot(co_velocity_axis_list[i],
                co_spectrum_list[i],
                #color='k',
                drawstyle = 'steps-mid',
                label=spectra_names[i].capitalize(),
                )

        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel(r'T$_b$ [K]')

    ax.legend(loc='upper left')

    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight', dpi=100)


''' Calculations
'''

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


''' DS9 Region and Coordinate Functions
'''

def convert_limit_coordinates(prop_dict,
        coords=('region_limit', 'co_noise_limits', 'plot_limit'), header=None):

    # Initialize pixel keys
    for coord in coords:
        prop_dict[coord].update({'pixel': []})

        if coord in ('region_limit',
                     'plot_limit',
                     'region_limit_bin',
                     'plot_limit_bin'):
            limit_wcs = prop_dict[coord]['wcs']

            for limits in limit_wcs:
                # convert centers to pixel coords
                limit_pixels = get_pix_coords(ra=limits[0],
                                             dec=limits[1],
                                             header=header)[:2].tolist()

                prop_dict[coord]['pixel'].append(limit_pixels[0])
                prop_dict[coord]['pixel'].append(limit_pixels[1])
        elif coord == 'co_noise_limits':
            region_limits = prop_dict[coord]['wcs']

            # Cycle through each region, convert WCS limits to pixels
            for region in region_limits:
                region_pixels = []
                for limits in region:
                    # convert centers to pixel coords
                    limit_pixels = get_pix_coords(ra=limits[0],
                                                  dec=limits[1],
                                                  header=header)[:2].tolist()
                    region_pixels.append(limit_pixels)

                # Append individual regions back to CO noise
                prop_dict[coord]['pixel'].append(region_pixels)

    return prop_dict

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

def load_ds9_region(props, filename=None, header=None):

    import pyregion as pyr

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    props['regions'] = {}


    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        region_name = tag[tag.find('text={')+6:tag.find('}')].lower()

        # Format vertices to be 2 x N array
        poly_verts = []
        for i in xrange(0, len(region.coord_list)/2):
            poly_verts.append((region.coord_list[2*i],
                               region.coord_list[2*i+1]))

        poly_verts_pix = []
        for i in xrange(0, len(poly_verts)):
            poly_verts_pix.append(get_pix_coords(ra=poly_verts[i][0],
                                            dec=poly_verts[i][1],
                                            header=header)[:-1][::-1].tolist())

        props['regions'][region_name] = {}
        props['regions'][region_name]['poly_verts'] = {}
        props['regions'][region_name]['poly_verts']['wcs'] = poly_verts
        props['regions'][region_name]['poly_verts']['pixel'] = poly_verts_pix

    return props

def calc_data():
    ''' Executes script.
    '''
    import numpy as np
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json


    cloud_list = ('taurus',
                  #'taurus1',
                  #'taurus2',
                  'california',
                  'perseus')
    data_dict = {}

    for i, cloud in enumerate(cloud_list):

        if cloud == 'taurus1' or cloud == 'taurus2':
            cloud_name = 'taurus'
        else:
            cloud_name = cloud

        # define directory locations
        output_dir = '/d/bip3/ezbc/%s/data/python_output/nhi_av/' % cloud_name
        figure_dir = '/d/bip3/ezbc/multicloud/figures/spectra/'
        av_dir = '/d/bip3/ezbc/%s/data/av/' % cloud_name
        hi_dir = '/d/bip3/ezbc/%s/data/hi/' % cloud_name
        co_dir = '/d/bip3/ezbc/%s/data/co/' % cloud_name
        core_dir = output_dir + 'core_arrays/'
        region_dir = '/d/bip3/ezbc/%s/data/' % cloud_name + \
                     'python_output/ds9_regions/'

        # global property filename
        global_property_filename = '{0}_global_properties'.format(cloud)
        property_dir = '/d/bip3/ezbc/{0}/data/python_output/'.format(cloud_name)

        av_data_type = 'planck'

        cloud_dict = {}

        # Load HI maps from Taurus, California, and Perseus
        hi_data, hi_header = fits.getdata(
                '/d/bip3/ezbc/{0}/data/hi/{0}_hi_galfa_'.format(cloud_name) + \
                        'cube_regrid_planckres_bin.fits',
                header=True)

        # Load CO maps from Taurus, California, and Perseus
        co_data, co_header = fits.getdata(
                    co_dir + '{0}_co_cfa_'.format(cloud_name) + \
                        'cube_regrid_planckres_bin.fits',
                header=True)
        av_data, av_header = \
                fits.getdata(av_dir + \
                             '{0}_av_planck_tau353_5arcmin_bin.fits'.format(cloud_name),
                             header=True)

        # Mask out region
        with open(property_dir + \
                  global_property_filename + '_' + av_data_type + \
                  '_scaled.txt', 'r') as f:
            global_props = json.load(f)

        # Load cloud division regions from ds9
        global_props = load_ds9_region(global_props,
                        filename='/d/bip3/ezbc/multicloud/data/' + \
                        'python_output/multicloud_divisions.reg',
                                header=av_header)

        # Derive relevant region
        region_vertices = \
            np.array(global_props['regions'][cloud]['poly_verts']['pixel'])
        region_mask = np.logical_not(myg.get_polygon_mask(av_data,
                                                          region_vertices))
        if 0:
            import matplotlib.pyplot as plt
            plt.imshow(np.ma.array(av_data,
                                   mask=region_mask), origin='lower')
            plt.colorbar()
            plt.show()
        hi_data[:, region_mask] = np.nan
        co_data[:, region_mask] = np.nan

        import matplotlib.pyplot as plt
        plt.close(); plt.clf();
        co = np.copy(co_data[30,:,:])
        co[region_mask] = np.nan
        plt.imshow(co, origin='lower')
        plt.savefig('/usr/users/ezbc/Desktop/comap_' + cloud_name + '.png')


        # sum along spatial axes
        cloud_dict['hi_spectrum'] = calc_global_spectrum(
                                        hi_cube=hi_data,
                                        statistic='median'
                                        )
        cloud_dict['hi_std'] = calc_global_spectrum(
                                        hi_cube=hi_data,
                                        statistic='std'
                                        )
        cloud_dict['co_spectrum'] = calc_global_spectrum(
                                        hi_cube=co_data,
                                        statistic='median'
                                        )
        vel_center = global_props['hi_velocity_width']['value']
        vel_width = global_props['hi_velocity_center']['value']
        #cloud_dict['hi_vel_range'] = (vel_center + vel_width / 2.0,
        #                              vel_center - vel_width / 2.0)
        #cloud_dict['hi_vel_range'] = global_props['hi_velocity_range_conf']
        cloud_dict['hi_vel_range'] = global_props['hi_velocity_range']

        # Calculate velocity
        cloud_dict['hi_velocity_axis'] = make_velocity_axis(
                                        hi_header)
        cloud_dict['co_velocity_axis'] = make_velocity_axis(
                                        co_header)

        data_dict[cloud] = cloud_dict

    return data_dict

def main():

    import numpy as np
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json
    import pickle

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    figure_dir = '/d/bip3/ezbc/multicloud/figures/spectra/'

    clobber_results = 1

    if clobber_results:
        data_dict = calc_data()
        with open(output_dir + 'multicloud_hi_spectra.pickle', 'w') as f:
            pickle.dump(data_dict, f)
    else:
        with open(output_dir + 'multicloud_hi_spectra.pickle', 'r') as f:
            data_dict = pickle.load(f)

    hi_velocity_axis_list = []
    co_velocity_axis_list = []
    hi_spectrum_list = []
    hi_vel_range_list = []
    hi_std_list = []
    co_spectrum_list = []
    spectra_names = []

    cloud_list = ('taurus',
                  #'taurus1',
                  #'taurus2',
                  'california',
                  'perseus')

    for cloud in cloud_list:
        hi_velocity_axis_list.append(data_dict[cloud]['hi_velocity_axis'])
        hi_spectrum_list.append(data_dict[cloud]['hi_spectrum'])
        hi_std_list.append(data_dict[cloud]['hi_std'])
        hi_vel_range_list.append(data_dict[cloud]['hi_vel_range'])
        co_velocity_axis_list.append(data_dict[cloud]['co_velocity_axis'])
        co_spectrum_list.append(data_dict[cloud]['co_spectrum'])
        spectra_names.append(cloud)

    figure_types = ('png', 'pdf')
    for figure_type in figure_types:
        #limits = [-40, 40, -9, 61]
        limits = None

        # scatter plot
        plot_spectra_grid(hi_velocity_axis_list=hi_velocity_axis_list,
                        hi_spectrum_list=hi_spectrum_list,
                        co_velocity_axis_list=co_velocity_axis_list,
                        co_spectrum_list=co_spectrum_list,
                        hi_std_list=hi_std_list,
                        hi_vel_range_list=hi_vel_range_list,
                        savedir=figure_dir,
                        limits=limits,
                        filename='multicloud_co_spectrum_global.%s' % \
                                figure_type,
                        show=False,
                        spectra_names=spectra_names,
                        #title=r'Global HI Profile' + \
                        #        ' of GMCs'
                        )


if __name__ == '__main__':
    main()


