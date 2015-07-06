#!/usr/bin/python

import numpy as np
import cloudpy

class KeyboardInterruptError(Exception): pass

def run_cloud_analysis(cloud_name='', region=None):


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
    av_filename = av_dir + \
       cloud_name + '_av_planck_tau353_5arcmin.fits'
    av_error_filename = av_dir + \
       cloud_name + '_av_error_planck_tau353_5arcmin.fits'
    hi_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres.fits'
    hi_error_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres_noise.fits'

    # Plot args
    residual_hist_filename_base = figure_dir + 'diagnostics/residuals/' + \
                                  cloud_name + '_residual_hist'
    residual_map_filename_base = figure_dir + 'diagnostics/residuals/' + \
                                  cloud_name + '_residual_map'
    likelihood_filename_base = figure_dir + 'diagnostics/likelihoods/' + \
                                  cloud_name + '_likelihood'

    plot_args = {
            'residual_hist_filename_base': residual_hist_filename_base,
            'residual_map_filename_base': residual_map_filename_base,
            'likelihood_filename_base': likelihood_filename_base,
            }

    if 1:
        import os
        os.system('rm -rf ' + hi_error_filename)

    region_filename = region_dir + 'multicloud_divisions.reg'
    cloud_filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/' + \
            cloud_name + '.pickle'
    diagnostic_filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/diagnostics/' + \
            cloud_name + '_diagnostic.txt'

    width_grid = np.arange(1, 75, 2*0.16667)
    dgr_grid = np.arange(0.001, 0.4, 3e-3)
    intercept_grid = np.arange(-2, 2, 0.1)
    intercept_grid = np.arange(-1, 1, 0.1)

    width_grid = width_grid
    dgr_grid = dgr_grid
    intercept_grid = intercept_grid

    # Define number of pixels in each bin
    binsize = 1.0 * 60.0 / 5.0

    if cloud_name == 'taurus':
        region = 'taurus'
    else:
        region = cloud_name

    print('Performing analysis on ' + cloud_name)

    cloud = cloudpy.Cloud(av_filename,
                          hi_filename,
                          av_error_filename=av_error_filename,
                          hi_error_filename=hi_error_filename,
                          cloud_prop_filename=prop_filename,
                          dgr_grid=dgr_grid,
                          intercept_grid=intercept_grid,
                          width_grid=width_grid,
                          residual_width_scale=3.0,
                          threshold_delta_dgr=0.01,
                          hi_noise_vel_range=[90,110],
                          vel_range_diff_thres=2,
                          init_vel_range=[-40,30],
                          verbose=True,
                          clobber_likelihoods=True,
                          binsize=binsize,
                          diagnostic_filename=diagnostic_filename,
                          plot_args=plot_args,
                          )

    cloud.run_analysis(region_filename=region_filename,
                       region=region)

    cloudpy.save(cloud, cloud_filename)

def main():

    from Queue import Queue
    from threading import Thread
    import multiprocessing


    q = Queue(maxsize=0)
    num_threads = 10

    clouds = {}
    #clouds['perseus'] = run_cloud_analysis('perseus')
    #clouds['california'] = run_cloud_analysis('california')
    #clouds['taurus'] = run_cloud_analysis('taurus')

    # Define a worker function for multiprocessing
    if 1:
        def worker(cloud):
            """thread worker function"""
            try:
                clouds[cloud] = run_cloud_analysis(cloud)
            except KeyboardInterrupt:
                raise KeyboardInterruptError()

            return clouds

        jobs = []
        for cloud in ('perseus', 'california', 'taurus'):
            try:
                p = multiprocessing.Process(target=worker, args=(cloud,))
                jobs.append(p)
                p.start()
            except KeyboardInterrupt():
                p.terminate()

    #for cloud in clouds:
    #    print cloud.props

    #clouds['taurus'] = run_cloud_analysis('taurus')



if __name__ == '__main__':
    main()



