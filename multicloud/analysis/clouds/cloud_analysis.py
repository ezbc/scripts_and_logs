#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import cloudpy

class KeyboardInterruptError(Exception): pass

def fit_background(av_data, background_mask=None, background_dim=1):

    from scipy.interpolate import interp2d
    from scipy.interpolate import SmoothBivariateSpline as spline

    if background_mask is None:
        background_mask = np.zeros(av_data.shape)

    if background_dim == 1:
        background = np.nanmean(av_data[~background_mask])

    if background_dim == 2:
        #av_data = np.ma.array(av_data, mask=background_mask)

        loc = np.where(~background_mask)

        x = loc[0]
        y = loc[1]

        z = av_data[~background_mask]

        assert z.shape == x.shape

        bbox = [0, av_data.shape[0], 0, av_data.shape[1]]

        result_interp = spline(x, y, z, bbox=bbox, kx=1, ky=1)

        x_grid, y_grid = np.where(av_data)

        background_flat = result_interp.ev(x_grid, y_grid)

        background = np.reshape(background_flat, av_data.shape)

    return background

def perform_background_subtraction(av_filename, background_mask=None,
        background_dim=1, background_filename=None, background_init=None,
        background_region_filename=None):

    # Import external modules
    # -----------------------
    from myio import check_file
    from mycoords import convert_limit_coordinates, get_pix_coords, \
                         hrs2degs, load_ds9_region
    from myimage_analysis import fit_background
    import pyfits as fits
    import mygeometry as myg

    av_data, av_header = fits.getdata(av_filename,
                                      header=True)

    if background_init is not None:
        av_data = av_data - background_init

    file_exists = check_file(background_filename, clobber=False)

    if not file_exists:
        props = {}

        print('writing new background')

        # Load background regions from ds9
        props = load_ds9_region(props,
                                filename=background_region_filename,
                                header=av_header,
                                key='background_regions')

        # Derive relevant region
        background_mask = np.ones(av_data.shape)
        for background_region in props['background_regions']:
            background_vertices = \
              props['background_regions']\
                   [background_region]['poly_verts']['pixel']

            # block off region
            background_mask_temp = ~np.logical_not(myg.get_polygon_mask(av_data,
                                                background_vertices))

            background_mask[background_mask_temp] = 0
        background_mask = ~np.logical_not(background_mask)

        # Fit the background
        background = fit_background(av_data, background_mask,
                background_dim=background_dim)

        fits.writeto(background_filename, background, av_header)
    else:
        background = fits.getdata(background_filename)

    return background

def plot_mask_map(cloud_results, load_original_av=True):

    filename = \
            cloud_results['figure_dir'] + 'diagnostics/maps/' + \
            cloud_results['filename_extension'] + '_mask_map.png'

    print('\nWriting mask map to \n' + filename)

    cloudpy.plot_mask_map(cloud_results['cloud'],
                      filename=filename,
                      load_original_av=load_original_av,
                      #limits=[10, 20, -1, 1],
                      )

def plot_likelihoods(cloud_results):

    likelihood_filename_base = \
            cloud_results['figure_dir'] + 'likelihood/' + \
            cloud_results['filename_extension'] + '_likelihood_'

    if 0:
        cloudpy.plot_likelihoods_hist(cloud=cloud_results['cloud'],
                          filename=likelihood_filename_base + 'wd.png',
                          )

    cloudpy.plot_likelihoods_hist(cloud=cloud_results['cloud'],
                      props=cloud_results['props'],
                      plot_axes=('widths', 'dgrs'),
                      show=0,
                      returnimage=False,
                      filename=likelihood_filename_base + 'wd.png',
                      #limits=[0, 100, 0.0, 0.3],
                      logscale=False,
                      )
    cloudpy.plot_likelihoods_hist(cloud=cloud_results['cloud'],
                      props=cloud_results['props'],
                      plot_axes=('widths', 'intercepts'),
                      show=0,
                      returnimage=False,
                      filename=likelihood_filename_base + 'wi.png',
                      #limits=[0, 100, -2, 2],
                      logscale=False,
                      )

def plot_dgr_intercept_progression(cloud_results):

    filename = \
            cloud_results['figure_dir'] + 'diagnostics/' + \
            cloud_results['filename_extension'] + '_dgr_intercept_progress.png'

    cloudpy.plot_dgr_intercept_progression(cloud_results['cloud'],
                      filename=filename,
                      #limits=[10, 20, -1, 1],
                      )

def plot_av_vs_nhi(cloud_results):

    filename_base = \
            cloud_results['figure_dir'] + 'diagnostics/' + \
            cloud_results['filename_extension'] + '_av_vs_nhi'

    cloud = cloud_results['cloud']
    props = cloud.props
    fit_params = {
                  'dgr': props['dust2gas_ratio']['value'],
                  'intercept': props['intercept']['value']}

    from astropy.io import fits
    from myimage_analysis import calculate_nhi

    if 0:
        av_data = fits.getdata(cloud.av_filename)
        if cloud.av_error_filename is not None:
            av_error_data = fits.getdata(cloud.av_error_filename)
        else:
            av_error_data = np.ones(av_data.shape) * cloud.av_error
        hi_data = fits.getdata(cloud.hi_filename)

    av_data = fits.getdata(cloud.av_filename_bin)
    av_error_data = fits.getdata(cloud.av_error_filename_bin)
    if cloud.av_error_filename_bin is not None:
        av_error_data = fits.getdata(cloud.av_error_filename_bin)
    else:
        av_error_data = np.ones(av_data.shape) * cloud.av_error
    hi_data = fits.getdata(cloud.hi_filename_bin)

    nhi_image = calculate_nhi(cube=hi_data,
                        velocity_axis=cloud.hi_vel_axis,
                        velocity_range=props['hi_velocity_range_max']['value'],
                        )

    nhi_image[nhi_image < 0] = np.nan

    if cloud.av_background is not None:
        av_data = av_data - cloud.av_background

    cloudpy.plot_av_vs_nhi(nhi_image[~cloud.mask],
                      av_data[~cloud.mask],
                      av_error=av_error_data[~cloud.mask],
                      filename=filename_base + '_masked.png',
                      fit_params=fit_params,
                      gridsize=(10,10),
                      #limits=[3,20, 0, 3],
                      std=0.22,
                      contour_plot=0,
                      #limits=[10, 20, -1, 1],
                      )
    if 1:
        av_data = fits.getdata(cloud.av_filename)
        if cloud.av_error_filename is not None:
            av_error_data = fits.getdata(cloud.av_error_filename)
        else:
            av_error_data = np.ones(av_data.shape) * cloud.av_error
        hi_data = fits.getdata(cloud.hi_filename)

    nhi_image = calculate_nhi(cube=hi_data,
                        velocity_axis=cloud.hi_vel_axis,
                        velocity_range=props['hi_velocity_range_max']['value'],
                        )

    nhi_image[nhi_image < 0] = np.nan

    cloudpy.plot_av_vs_nhi(nhi_image,
                      av_data,
                      av_error=av_error_data,
                      filename=filename_base + '.png',
                      fit_params=fit_params,
                      gridsize=(10,10),
                      #limits=[1,20, 0, 4],
                      std=0.22,
                      contour_plot=True,
                      #limits=[10, 20, -1, 1],
                      )

def plot_nh2_vs_nhi(cloud_results):

    filename_base = \
            cloud_results['figure_dir'] + 'diagnostics/' + \
            cloud_results['filename_extension'] + '_nh2_vs_nhi'

    cloud = cloud_results['cloud']
    props = cloud.props
    fit_params = {
                  'dgr': props['dust2gas_ratio']['value'],
                  'intercept': props['intercept']['value']}

    from astropy.io import fits
    from myimage_analysis import calculate_nhi

    av_data = fits.getdata(cloud.av_filename)
    if cloud.av_error_filename is not None:
        av_error_data = fits.getdata(cloud.av_error_filename)
    else:
        av_error_data = np.ones(av_data.shape) * cloud.av_error
    hi_data = fits.getdata(cloud.hi_filename)

    nhi_image = calculate_nhi(cube=hi_data,
                        velocity_axis=cloud.hi_vel_axis,
                        velocity_range=props['hi_velocity_range_max']['value'],
                        )

    nhi_image[nhi_image < 0] = np.nan

    nh2_image = 0.5 * ((av_data - fit_params['intercept']) / \
                        fit_params['dgr'] - nhi_image)

    cloudpy.plot_nh2_vs_nhi(nhi_image,
                      nh2_image,
                      filename=filename_base + '.png',
                      fit_params=fit_params,
                      gridsize=(10,10),
                      #limits=[3,20, 0, 3],
                      #std=0.22,
                      contour_plot=1,
                      #limits=[10, 20, -1, 1],
                      )

def make_map_movie(cloud_results):

    filename = \
            cloud_results['figure_dir'] + 'diagnostics/' + \
            cloud_results['filename_extension'] + '_residual_maps.gif'

    cloudpy.plot_map_movie(cloud_results['cloud'], filename=filename)

def make_residual_hist_movie(cloud_results):

    filename = \
            cloud_results['figure_dir'] + 'diagnostics/' + \
            cloud_results['filename_extension'] + '_residual_hists.gif'

    cloudpy.plot_residual_hist_movie(cloud_results['cloud'], filename=filename)

def run_cloud_analysis(args,
    cloud_name = 'perseus',
    region = None,
    load = False,
    data_type = 'planck'):

    if 1:
        cloud_name = args['cloud_name']
        region = args['region']
        load = args['load']
        data_type = args['data_type']
        background_subtract = args['background_subtract']

    # define directory locations
    # --------------------------
    figure_dir = \
        '/d/bip3/ezbc/' + cloud_name + '/figures/'
    av_dir = '/d/bip3/ezbc/' + cloud_name + '/data/av/'
    hi_dir = '/d/bip3/ezbc/' + cloud_name + '/data/hi/'
    co_dir = '/d/bip3/ezbc/' + cloud_name + '/data/co/'
    core_dir = \
       '/d/bip3/ezbc/' + cloud_name + '/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/' + cloud_name + '/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    background_region_dir = '/d/bip3/ezbc/' + cloud_name + \
                            '/data/python_output/ds9_regions/'
    likelihood_dir = \
            '/d/bip3/ezbc/' + cloud_name + '/data/python_output/nhi_av/'

    # define filenames
    prop_filename = property_dir + \
       cloud_name + '_global_properties.txt'
    hi_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres.fits'
    hi_error_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres_noise.fits'
    if cloud_name == 'perseus' and data_type == 'lee12':
        av_filename = av_dir + \
           cloud_name + '_av_lee12_iris_regrid_planckres.fits'
        av_error_filename = None
        av_error = 0.1
        if background_subtract:
            av_background = 0.5
        else:
            av_background = None
    if data_type == 'planck':
        av_filename = av_dir + \
           cloud_name + '_av_planck_tau353_5arcmin.fits'
        av_error_filename = av_dir + \
           cloud_name + '_av_error_planck_tau353_5arcmin.fits'
        av_error = None

        if 0:
            av_error_filename = None
            av_error = 0.1

        av_background = 0.0
    if data_type == 'k09':
        av_filename = av_dir + \
           cloud_name + '_av_regrid_planckres.fits'

        av_error_filename = None
        av_error = 0.1

        av_background = 0.0

    # Name of diagnostic files
    if background_subtract:
        background_name = '_backsub'
    else:
        background_name = ''

    if args['bin_image']:
        bin_name = '_binned'
    else:
        bin_name = ''

    if args['use_weights']:
        weights_name = '_weights'
        weights_filename = av_dir + \
           cloud_name + '_bin_weights.fits'
    else:
        weights_name = ''
        weights_filename = None

    filename_extension = cloud_name + '_' + data_type + background_name + \
            bin_name + weights_name + '_' + args['likelihood_resolution'] + \
            'res'

    # Plot args
    residual_hist_filename_base = figure_dir + 'diagnostics/residuals/' + \
                                  filename_extension + '_residual_hist'
    residual_map_filename_base = figure_dir + 'diagnostics/residuals/' + \
                                  filename_extension + '_residual_map'
    likelihood_filename_base = figure_dir + 'diagnostics/likelihoods/' + \
                                  filename_extension + '_likelihood'
    av_bin_map_filename_base = figure_dir + 'diagnostics/maps/' + \
                                  filename_extension + '_bin_map'

    plot_args = {
            'residual_hist_filename_base': residual_hist_filename_base,
            'residual_map_filename_base': residual_map_filename_base,
            'likelihood_filename_base': likelihood_filename_base,
            'av_bin_map_filename_base' : av_bin_map_filename_base,
            }

    if 1:
        import os
        os.system('rm -rf ' + hi_error_filename)

    region_filename = region_dir + 'multicloud_divisions.reg'


    cloud_filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/' + \
            filename_extension + \
            '.pickle'
    props_filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/' + \
            filename_extension + \
            '_props.pickle'

    # background filenames
    background_filename = av_dir + filename_extension + '_background.fits'
    background_region_filename = background_region_dir + cloud_name + \
                                 '_background.reg'

    diagnostic_filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/diagnostics/' + \
            filename_extension + '_diagnostic.txt'

    if args['likelihood_resolution'] == 'fine':
        width_grid = np.arange(1, 100, 2*0.16667)
        dgr_grid = np.arange(0.000, 0.3, 2e-4)
        intercept_grid = np.arange(-2, 2, 0.005)
    elif args['likelihood_resolution'] == 'coarse':
        width_grid = np.arange(1, 100, 2*0.16667)
        #width_grid = np.arange(100, 101, 1)
        dgr_grid = np.arange(0.000, 0.3, 3e-3)
        intercept_grid = np.arange(-2, 2, 0.1)
        #intercept_grid = np.arange(-0.5, 0.5, 0.01)
        #intercept_grid = np.array((0.9,))
    else:
        raise ValueError('likelihood_resolution should be either ' + \
                         '"coarse" or "fine"')

    # Define number of pixels in each bin
    # size of binned pixel in degrees * number of arcmin / degree * number of
    # arcmin / pixel
    binsize = 0.5 * 60.0 / 5.0

    if cloud_name == 'taurus':
        region = 'taurus'
    else:
        region = cloud_name

    if not load:
        print('\n\nPerforming analysis on ' + cloud_name)

        if cloud_name == 'california':
            if args['background_subtract']:
                print('\nPerforming a background subtraction...')
                av_filename_back = av_filename.replace('.fits', '_bin.fits')
                av_background = perform_background_subtraction(av_filename_back,
                       background_dim=2,
                       #background_init=0.9,
                       background_filename=\
                               background_filename,
                       background_region_filename=\
                               background_region_filename)
                #intercept_grid = np.arange(0, 1, 1)
                av_background = 0.9
                #results['av_background'] = av_background

        cloud = cloudpy.Cloud(av_filename,
                              hi_filename,
                              av_error_filename=av_error_filename,
                              av_error=av_error,
                              av_background=av_background,
                              hi_error_filename=hi_error_filename,
                              cloud_prop_filename=prop_filename,
                              dgr_grid=dgr_grid,
                              intercept_grid=intercept_grid,
                              width_grid=width_grid,
                              residual_width_scale=3,
                              threshold_delta_dgr=0.001,
                              hi_noise_vel_range=[90,110],
                              vel_range_diff_thres=5,
                              init_vel_range=[-50, 50],
                              verbose=True,
                              clobber_likelihoods=True,
                              binsize=binsize,
                              use_bin_weights=args['use_weights'],
                              use_only_binned_data=args['bin_image'],
                              pixel_mask_increase_fraction=0.05,
                              binned_data_filename_ext=\
                                args['binned_data_filename_ext'],
                              weights_filename=weights_filename,
                              #diagnostic_filename=diagnostic_filename,
                              av_mask_threshold=args['av_mask_threshold'],
                              plot_args=plot_args,
                              )

        cloud.run_analysis(region_filename=region_filename,
                           region=region)

        print('\nSaving cloud to file \n' + cloud_filename)
        cloudpy.save(cloud, cloud_filename)
        cloudpy.save(cloud.props, props_filename)
        props = cloud.props
    else:
        if not args['load_props']:
            print('\nLoading cloud from file \n' + cloud_filename)
            cloud = cloudpy.load(cloud_filename)
            cloudpy.save(cloud.props, props_filename)
            props = cloud.props
        else:
            print('\nLoading cloud props from file \n' + props_filename)
            cloud = None
            props = cloudpy.load(props_filename)

    results = {}
    results['cloud'] = cloud
    results['props'] = props
    results['filename_extension'] = filename_extension
    results['figure_dir'] = figure_dir

    return results

def main():

    from Queue import Queue
    from threading import Thread
    import multiprocessing

    #load_clouds = 0

    results = {}
    #clouds['perseus'] = run_cloud_analysis('perseus')
    #clouds['california'] = run_cloud_analysis('california')
    #clouds['taurus'] = run_cloud_analysis('taurus')

    clouds = (
              #'perseus',
              #'taurus',
              'california',
              )

    for cloud in clouds:
        args = {'cloud_name':cloud,
                'region':None,
                'load': 0,
                'load_props': 0,
                'data_type': 'planck',
                #'data_type': 'lee12',
                'background_subtract': 0,
                'bin_image': True,
                'use_weights': False,
                'av_mask_threshold': None,
                'binned_data_filename_ext': '_bin',
                'likelihood_resolution': 'fine',
                #'likelihood_resolution': 'coarse',
                }

        results[cloud] = run_cloud_analysis(args)

        print('\nDGR max = {0:e}'.format(results[cloud]['cloud']\
                                       .props['dust2gas_ratio_max']['value']))
        print('Intercept max = {0:e}'.format(results[cloud]['cloud']\
                                       .props['intercept_max']['value']))
        print('HI width max = {0:e}'.format(results[cloud]['cloud']\
                                .props['hi_velocity_width_max']['value']))

        plot_likelihoods(results[cloud])
        plot_mask_map(results[cloud])
        plot_dgr_intercept_progression(results[cloud])
        plot_av_vs_nhi(results[cloud])
        plot_nh2_vs_nhi(results[cloud])
        make_map_movie(results[cloud])
        make_residual_hist_movie(results[cloud])



if __name__ == '__main__':
    main()



