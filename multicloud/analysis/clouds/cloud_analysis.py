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

        #print av_data.size, z.shape
        assert z.shape == x.shape

        bbox = [0, av_data.shape[0], 0, av_data.shape[1]]

        result_interp = spline(x, y, z, bbox=bbox, kx=1, ky=1)

        x_grid, y_grid = np.where(av_data)

        background_flat = result_interp.ev(x_grid, y_grid)

        background = np.reshape(background_flat, av_data.shape)

    return background

def setup_background_subtraction():

    pass

def plot_likelihoods(cloud):

    likelihood_filename_base = \
            '/d/bip3/ezbc/perseus/figures/likelihood/' + \
            'perseus_likelihood_lee12_'

    cloudpy.plot_likelihoods_hist(cloud=cloud,
                      plot_axes=('widths', 'dgrs'),
                      show=0,
                      returnimage=False,
                      filename=likelihood_filename_base + 'wd.png',
                      limits=[10, 20, 0.05, 0.15],
                      )
    cloudpy.plot_likelihoods_hist(cloud=cloud,
                      plot_axes=('widths', 'intercepts'),
                      show=0,
                      returnimage=False,
                      filename=likelihood_filename_base + 'wi.png',
                      limits=[10, 20, -1, 1],
                      )

def plot_dgr_intercept_progression(cloud):

    filename = \
            '/d/bip3/ezbc/perseus/figures/diagnostics/' + \
            'perseus_dgr_intercept_progress_lee12.png'

    cloudpy.plot_dgr_intercept_progression(cloud,
                      filename=filename,
                      #limits=[10, 20, -1, 1],
                      )

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
    output_dir = '/d/bip3/ezbc/' + cloud_name + '/data/python_output/nhi_av/'
    figure_dir = \
        '/d/bip3/ezbc/' + cloud_name + '/figures/'
    av_dir = '/d/bip3/ezbc/' + cloud_name + '/data/av/'
    hi_dir = '/d/bip3/ezbc/' + cloud_name + '/data/hi/'
    co_dir = '/d/bip3/ezbc/' + cloud_name + '/data/co/'
    core_dir = \
       '/d/bip3/ezbc/' + cloud_name + '/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/' + cloud_name + '/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
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
        av_error = None,
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

    diagnostic_filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/diagnostics/' + \
            filename_extension + '_diagnostic.txt'

    if args['likelihood_resolution'] == 'fine':
        width_grid = np.arange(1, 75, 2*0.16667)
        dgr_grid = np.arange(0.000, 0.2, 2e-4)
        intercept_grid = np.arange(-2, 2, 0.01)
    elif args['likelihood_resolution'] == 'coarse':
        width_grid = np.arange(1, 50, 2*0.16667)
        dgr_grid = np.arange(0.000, 0.2, 3e-3)
        intercept_grid = np.arange(-1, 2, 0.1)
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
        print('Performing analysis on ' + cloud_name)

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
                              residual_width_scale=1.5,
                              threshold_delta_dgr=0.001,
                              hi_noise_vel_range=[90,110],
                              vel_range_diff_thres=2,
                              init_vel_range=[-40,30],
                              verbose=True,
                              clobber_likelihoods=True,
                              binsize=binsize,
                              use_bin_weights=args['use_weights'],
                              use_only_binned_data=args['bin_image'],
                              pixel_mask_increase_fraction=0.03,
                              binned_data_filename_ext=\
                                args['binned_data_filename_ext'],
                              weights_filename=weights_filename,
                              #diagnostic_filename=diagnostic_filename,
                              plot_args=plot_args,
                              )

        cloud.run_analysis(region_filename=region_filename,
                           region=region)

        cloudpy.save(cloud, cloud_filename)
    else:
        print('Loading cloud from file \n' + cloud_filename)
        cloud = cloudpy.load(cloud_filename)

    results = {}
    results['cloud'] = cloud
    results['filename_extension'] = filename_extension
    results['figure_dir'] = figure_dir

    return results

def main():

    from Queue import Queue
    from threading import Thread
    import multiprocessing

    load_clouds = 0

    q = Queue(maxsize=0)
    num_threads = 10

    results = {}
    #clouds['perseus'] = run_cloud_analysis('perseus')
    #clouds['california'] = run_cloud_analysis('california')
    #clouds['taurus'] = run_cloud_analysis('taurus')

    clouds = (
              #'perseus',
              'california',
              'taurus'
              )

    if 1:
        for cloud in clouds:
            args = {'cloud_name':cloud,
                    'region':None,
                    'load': 0,
                    'data_type': 'planck',
                    'background_subtract': False,
                    'bin_image': True,
                    'use_weights': False,
                    'binned_data_filename_ext': '_bin',
                    'likelihood_resolution': 'fine',
                    }

            results[cloud] = run_cloud_analysis(args)

    # Define a worker function for multiprocessing
    if 0:
        def worker(args):
            """thread worker function"""
            try:
                results[cloud] = run_cloud_analysis(args)
            except KeyboardInterrupt:
                raise KeyboardInterruptError()

            return results

        jobs = []
        clouds = (
                  'perseus',
                  #'california',
                  #'taurus'
                  )
        data_type = 'lee12'

        for cloud in clouds:
            args = {'cloud_name':cloud,
                    'region':None,
                    'load': True,
                    'data_type': 'lee12'}

            try:
                p = multiprocessing.Process(target=worker, args=(args,))
                jobs.append(p)
                p.start()
                p.join()
            except KeyboardInterrupt():
                p.terminate()

        p.terminate()

    plot_likelihoods(results['perseus'], args)
    plot_dgr_intercept_progression(results['perseus'], args)



if __name__ == '__main__':
    main()



