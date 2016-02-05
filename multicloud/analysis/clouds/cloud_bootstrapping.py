#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import myimage_analysis as myia
from multiprocessing.queues import Queue
import mygeometry as myg
import scipy

'''
Plotting
'''

from local_module_plotting import *

'''
Data Prep functions
'''

from local_module_dataprep import *

'''
Multiprocessing functions
'''

from local_module_multiprocessing import *

'''
Region functions
'''

from local_module_regions import *

'''
Modeling Functions
'''

from local_module_fitting import *

'''
Bootstrapping functions
'''

from local_module_bootstrapping import *

'''
Main function
'''
def get_results(global_args):

    import myio

    print('\nPerforming analysis on ' + global_args['cloud_name'])
    print('=======================' + '=' * len(global_args['cloud_name']))

    # Get the results filename
    filename_base, global_args = create_filename_base(global_args)
    print('\n\tFilename base = \n\t' + filename_base)
    results_dir =  '/d/bip3/ezbc/multicloud/data/python_output/'
    results_filename = results_dir + \
               'bootstrap_results/' + filename_base + \
               '_bootstrap_results.pickle'
    global_args['results_filename'] = results_filename

    exists = myio.check_file(results_filename)

    # either load or perform analysis
    if global_args['load'] and exists:
        print('\n\tLoading results...')
        results_dict = load_results(global_args['results_filename'])

    else:
        results_dict = run_cloud_analysis(global_args)

    # derive col dens images and statistics on MC sim
    add_results_analysis(results_dict)

    # calculate errors on dgrs and intercept
    results_dict['params_summary'] = calc_param_errors(results_dict)

    return results_dict

def run_cloud_analysis(global_args,):

    from astropy.io import fits
    from myimage_analysis import calculate_nhi, calc_region_mask
    import myimage_analysis as myia
    from mycoords import make_velocity_axis
    from mystats import calc_symmetric_error, calc_logL
    import os
    import myio
    import pickle
    import mystats

    cloud_name = global_args['cloud_name']
    region = global_args['region']
    load = global_args['load']
    data_type = global_args['data_type']
    background_subtract = global_args['background_subtract']


    # define directory locations
    # --------------------------
    figure_dir = \
        '/d/bip3/ezbc/multicloud/figures/'
    av_dir = '/d/bip3/ezbc/' + cloud_name + '/data/av/'
    dust_temp_dir = '/d/bip3/ezbc/' + cloud_name + '/data/dust_temp/'
    hi_dir = '/d/bip3/ezbc/' + cloud_name + '/data/hi/'
    co_dir = '/d/bip3/ezbc/' + cloud_name + '/data/co/'
    core_dir = \
       '/d/bip3/ezbc/' + cloud_name + '/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/' + cloud_name + '/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'
    background_region_dir = '/d/bip3/ezbc/' + cloud_name + \
                            '/data/python_output/ds9_regions/'
    results_dir =  '/d/bip3/ezbc/multicloud/data/python_output/'

    # define filenames
    prop_filename = property_dir + \
       cloud_name + '_global_properties.txt'
    hi_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres.fits'
    hi_error_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres_noise.fits'
    co_filename = co_dir + \
       cloud_name + '_co_cfa_cube_regrid_planckres.fits'

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
            av_error = 1
        av_background = None
        av_ref_filename = av_dir + \
           cloud_name + '_av_k09_regrid_planckres.fits'
    if cloud_name == 'perseus' and data_type == 'planck_lee12mask':
        av_filename = av_dir + \
           cloud_name + '_av_planck_tau353_5arcmin_lee12mask.fits'
        av_error_filename = av_dir + \
           cloud_name + '_av_error_planck_tau353_5arcmin.fits'
        av_error = None
        av_background = None
    if data_type == 'k09':
        av_filename = av_dir + \
           cloud_name + '_av_k09_regrid_planckres.fits'

        av_error_filename = None
        av_error = 0.4

        av_background = 0.0

    # Get the filename base to differentiate between different parameters
    filename_base, global_args = create_filename_base(global_args)

    # set up plotting variables
    plot_kwargs = {
                   'figure_dir': figure_dir,
                   'cloud_name': cloud_name,
                   'filename_base': filename_base,
                   'plot_diagnostics': global_args['plot_diagnostics'],
                   #'av_nhi_contour': av_nhi_contour,
                   'av_nhi_contour': True,
                   'av_nhi_limits': [0, 20, -1, 9],
                   #'av_nhi_limits': None,
                    }

    # Load data
    if global_args['bin_image']:
        av_filename = av_filename.replace('.fits', '_bin.fits')
        if av_error_filename is not None:
            av_error_filename = av_error_filename.replace('.fits', '_bin.fits')
        hi_filename = hi_filename.replace('.fits', '_bin.fits')
        av_nhi_contour = False
    else:
        av_nhi_contour = True

    av_data, av_header = fits.getdata(av_filename, header=True)
    av_data_ref, av_header = fits.getdata(av_ref_filename, header=True)
    if av_error_filename is not None:
        av_error_data, av_error_header = fits.getdata(av_error_filename,
                                                      header=True)
    else:
        av_error_data = av_error * np.ones(av_data.shape)

    # mask data
    region_filename = region_dir + 'multicloud_divisions.reg'
    region_mask = calc_region_mask(region_filename,
                                   av_data,
                                   av_header,
                                   region_name=global_args['region_name'])

    av_data[region_mask] = np.nan
    av_data_ref[region_mask] = np.nan

    # Scale the data to the 2MASS K+09 data
    scale_kwargs = scale_av_with_refav(av_data, av_data_ref, av_error_data)
    av_data_backsub = av_data - scale_kwargs['intercept']
    avg_scalar = (scale_kwargs['av_scalar'] + 1) / 2.0
    av_data_backsub_scaled = av_data_backsub / avg_scalar

    if 1:
        print('\n\tAv scaling stats:')
        print('\t\tSlope: ' + \
              '{0:.2f} +/- {1:.2f}'.format(scale_kwargs['av_scalar'],
                                           scale_kwargs['av_scalar_error']))
        print('\t\tIntercept: ' + \
              '{0:.2f} +/- {1:.2f}'.format(scale_kwargs['intercept'],
                                           scale_kwargs['intercept_error']))


    # Load HI and CO cubes
    hi_data, hi_header = fits.getdata(hi_filename, header=True)
    co_data, co_header = fits.getdata(co_filename, header=True)


    hi_data[:, region_mask] = np.nan
    co_data[:, region_mask] = np.nan

    hi_vel_axis = make_velocity_axis(hi_header)
    co_vel_axis = make_velocity_axis(co_header)

    # Load HI error
    if global_args['clobber_hi_error']:
        print('\n\tCalculating HI noise cube...')
        os.system('rm -rf ' + hi_error_filename)
        hi_data_error = \
            myia.calculate_noise_cube(cube=hi_data,
                                      velocity_axis=hi_vel_axis,
                                      velocity_noise_range=[-110,-90, 90,110],
                                      Tsys=30.0,
                                      filename=hi_error_filename)
    else:
        hi_data_error = fits.getdata(hi_error_filename)


    # Derive N(HI)
    # -------------------------------------------------------------------------
    # get fit kwargs
    gauss_fit_kwargs, ncomps_in_cloud = get_gauss_fit_kwargs(global_args)

    # derive spectra or load
    spectra_filename = results_dir + 'spectra/' + global_args['cloud_name'] + \
            '_spectra.pickle'
    load_spectra = myio.check_file(spectra_filename,
                                   clobber=global_args['clobber_spectra'])
    if load_spectra:
        hi_spectrum, hi_std_spectrum, co_spectrum = \
                myio.load_pickle(spectra_filename)
    else:
        print('\n\tCalculating spectra...')
        if global_args['smooth_hi_to_co_res']:
            from astropy.convolution import Gaussian2DKernel, convolve
            # Create kernel
            # one pix = 5 arcmin, need 8.4 arcmin for CO res
            # The beamsize is the FWHM. The convolution kernel needs the
            # standard deviation
            hi_res = 1.0
            co_res = 8.4 / 5.0
            width = (co_res**2 - hi_res**2)**0.5
            std = width / 2.355
            g = Gaussian2DKernel(width)

            # Convolve data
            hi_data_co_res = np.zeros(hi_data.shape)
            for i in xrange(hi_data.shape[0]):
                hi_data_co_res[i, :, :] = \
                    convolve(hi_data[i, :, :], g, boundary='extend')

        hi_spectrum = myia.calc_spectrum(hi_data_co_res)
        hi_std_spectrum = myia.calc_spectrum(hi_data_co_res,
                                             statistic=np.nanstd)
        co_spectrum = myia.calc_spectrum(co_data)
        myio.save_pickle(spectra_filename,
                         (hi_spectrum, hi_std_spectrum, co_spectrum))

    if global_args['hi_range_calc'] == 'gaussian':
        velocity_range, gauss_fits, comp_num, hi_range_error = \
                calc_hi_vel_range(hi_spectrum,
                                  hi_vel_axis,
                                  gauss_fit_kwargs,
                                  co_spectrum=co_spectrum,
                                  co_vel_axis=co_vel_axis,
                                  ncomps=ncomps_in_cloud,
                                  )
        global_args['vel_range_error'] = hi_range_error
    else:
        velocity_range = [-5, 15]
        gauss_fits = None
        comp_num = None

    hi_range_kwargs = {
                       'velocity_range': velocity_range,
                       'gauss_fits': gauss_fits,
                       'comp_num': comp_num,
                       'hi_range_error': hi_range_error,
                       'vel_range': velocity_range,
                       'gauss_fit_kwargs': gauss_fit_kwargs,
                       }

    # plot the results
    filename = plot_kwargs['figure_dir'] + \
               'spectra/' + plot_kwargs['filename_base'] + \
               '_spectra.png'
    plot_spectra(hi_spectrum,
                 hi_vel_axis,
                 hi_std_spectrum=hi_std_spectrum,
                 gauss_fits=gauss_fits,
                 comp_num=comp_num,
                 co_spectrum=co_spectrum,
                 co_vel_axis=co_vel_axis,
                 vel_range=velocity_range,
                 filename=filename,
                 limits=[-50, 30, -10, 70],
                 )

    if 0:
        print('\n\tVelocity range = ' + \
              '{0:.1f} to {1:.1f}'.format(*velocity_range))

    # use the vel range to derive N(HI)
    nhi_image, nhi_image_error = \
        calculate_nhi(cube=hi_data,
                      velocity_axis=hi_vel_axis,
                      velocity_range=velocity_range,
                      noise_cube=hi_data_error,
                      return_nhi_error=True,
                      )
    nhi_error_median = scipy.stats.nanmedian(nhi_image_error, axis=None)
    print('\n\tMedian N(HI) error = ' + \
          '{0:.2f} x 10^20 cm^-2'.format(nhi_error_median))

    nhi_image_background = calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=(-100,velocity_range[0]),
                              )
    nhi_image_background += calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=(velocity_range[1],100),
                              )

    # mask for erroneous pixels
    nhi_image[nhi_image < 0] = np.nan
    nhi_image_error[nhi_image_error < 0] = np.nan
    nhi_image_background[nhi_image_background < 0] = np.nan

    if not global_args['use_background']:
        nhi_image_background = None


    # Write filenames
    filenames = {
                 'region_filename': region_filename,
                 'av_filename': av_filename,
                 'av_error_filename': av_error_filename,
                 'av_ref_filename': av_ref_filename,
                 'hi_filename': hi_filename,
                 'hi_error_filename': hi_error_filename,
                 'co_filename': co_filename,
                 'results_dir': results_dir,
                 'figure_dir': figure_dir,
                 }

    # Collect data
    data = {
            'av_data': av_data,
            'av_error_data': av_error_data,
            'av_data_ref': av_data_ref,
            'hi_data': hi_data,
            'hi_data_error': hi_data_error,
            'hi_vel_axis': hi_vel_axis,
            'co_data': co_data,
            'co_vel_axis': co_vel_axis,
            'av_header': av_header,
            'av_error_header': av_error_header,
            'hi_header': hi_header,
            'co_header': co_header,
            }

    # Collect data products
    data_products = {
                     'av_data_backsub': av_data_backsub,
                     'av': av_data_backsub_scaled,
                     'av_error': av_error_data,
                     'nhi': nhi_image,
                     'nhi_error': nhi_image_error,
                     'nhi_image_background': nhi_image_background,
                     'region_mask': region_mask,
                     'scale_kwargs': scale_kwargs,
                     'hi_spectrum': hi_spectrum,
                     'hi_std_spectrum': hi_std_spectrum,
                     'co_spectrum': co_spectrum,
                     'hi_range_kwargs': hi_range_kwargs,
                     }

    # Get model fitting params
    model_fitting = get_model_fit_kwargs(cloud_name,
                                         vary_phi_g=global_args['vary_phi_g'])
    model_fitting['sternberg_params']['radiation_type'] = \
            global_args['radiation_type']

    # Get cores params
    cores = get_core_properties(data, cloud_name)

    # Calculate average dust temperature of each core
    #cores = add_core_temps(data, )

    # get the cores in the cloud
    cores_to_plot = get_cores_to_plot()
    cores_to_plot = trim_cores_to_plot(cores, cores_to_plot)

    if 0:
        import sys
        sys.exit()

    global_args['ss_model_kwargs'] = {}
    global_args['ss_model_kwargs']['cores'] = cores
    global_args['ss_model_kwargs']['cores_to_plot'] = cores_to_plot
    global_args['ss_model_kwargs']['model_kwargs'] = model_fitting

    # Bootstrap data
    # -------------------------------------------------------------------------

    # crop hi_data to be a reasonable size
    hi_data_crop, hi_vel_axis_crop = myia.crop_cube(hi_data,
                                                    hi_vel_axis,
                                                    [-20, 30])
    hi_data_error_crop, hi_vel_axis_crop = myia.crop_cube(hi_data_error,
                                                    hi_vel_axis,
                                                    [-20, 30])

    bootstrap_filename = results_dir + filename_base + '_bootresults.npy'
    results_filename = results_dir + \
               'bootstrap_results/' + filename_base + \
               '_bootstrap_results.pickle'

    # Bootstrap residuals of best fitting models
    if global_args['bootstrap_fit_residuals']:
        print('\n\tBeginning residual bootstrapping...')
        resid_mc_results = \
            bootstrap_residuals(av_data_backsub,
                           nhi_image=nhi_image,
                           nhi_image_error=nhi_image_error,
                           av_error_data=av_error_data,
                           nhi_image_background=nhi_image_background,
                           plot_kwargs=plot_kwargs,
                           hi_data=hi_data_crop,
                           hi_data_error=hi_data_error_crop,
                           vel_axis=hi_vel_axis_crop,
                           vel_range=velocity_range,
                           vel_range_error=2,
                           av_reference=av_data_ref,
                           use_intercept=global_args['use_intercept'],
                           num_bootstraps=global_args['num_resid_bootstraps'],
                           scale_kwargs=scale_kwargs,
                           sim_hi_error=global_args['sim_hi_error'],
                           ss_model_kwargs=global_args['ss_model_kwargs'],
                           multiprocess=global_args['multiprocess'],
                           rotate_cores=global_args['rotate_cores'],
                           )
    else:
        resid_mc_results = None

    print('\n\tBeginning bootstrap monte carlo...')
    # Perform bootsrapping
    boot_result, mc_results = \
        bootstrap_fits(av_data_backsub,
                       nhi_image=nhi_image,
                       nhi_image_error=nhi_image_error,
                       av_error_data=av_error_data,
                       nhi_image_background=nhi_image_background,
                       plot_kwargs=plot_kwargs,
                       hi_data=hi_data_crop,
                       hi_data_error=hi_data_error_crop,
                       vel_axis=hi_vel_axis_crop,
                       vel_range=velocity_range,
                       vel_range_error=hi_range_error,
                       av_reference=av_data_ref,
                       use_intercept=global_args['use_intercept'],
                       num_bootstraps=global_args['num_bootstraps'],
                       scale_kwargs=scale_kwargs,
                       sim_hi_error=global_args['sim_hi_error'],
                       ss_model_kwargs=global_args['ss_model_kwargs'],
                       multiprocess=global_args['multiprocess'],
                       rotate_cores=global_args['rotate_cores'],
                       calc_median_error=global_args['calculate_median_error']
                       )
    np.save(bootstrap_filename, boot_result)

    results_dict = {'boot_result': boot_result,
                    'data': data,
                    'data_products': data_products,
                    'global_args': global_args,
                    'plot_kwargs': plot_kwargs,
                    'filenames': filenames,
                    'mc_results': mc_results,
                    'resid_mc_results': resid_mc_results,
                    }

    print('\n\tSaving results...')
    save_results(results_dict, global_args['results_filename'])
    #results_dict = load_results(global_args['results_filename'])

    return results_dict

def main():

    import itertools

    results = {}

    clouds = (
              'taurus',
              'california',
              'perseus',
              )

    data_types = (
                  'planck',
                  #'lee12',
                  #'planck_lee12mask',
                  #'k09',
                  )
    recalculate_likelihoods = (
                               #True,
                               False,
                               )
    bin_image = (
                 #True,
                 False,
                 )

    init_vel_width = (#20,
                      20,
                      #200,
                      )
    fixed_width = (
                   #20,
                   #'gaussfit',
                   None,
                   #50,
                   )

    hi_range_calc = ('gaussian',
                     #'std',
                     )

    use_intercept = (
                     False,
                     #True,
                     )

    use_background = (
                      False,
                      #True,
                      )
    av_mask_threshold = (
                         None,
                         #1.2,
                         #2.5,
                         )

    regions = (None,
               #'1',
               #'2'
               )

    subtract_comps = (#True,
                      False,
                      )

    radiation_type = (#'beamed',
                      'isotropic',
                      )

    rotate_cores = (
                    False,
                    #True,
                    )

    vary_phi_g = (
                    #True,
                    False,
                    )

    elements = (clouds, data_types, recalculate_likelihoods, bin_image,
            init_vel_width, fixed_width, use_intercept, av_mask_threshold,
            regions, subtract_comps, use_background, hi_range_calc,
            radiation_type, rotate_cores, vary_phi_g)

    permutations = list(itertools.product(*elements))

    print('Number of permutations to run: ' + str(len(permutations)))

    #for cloud in clouds:
    for permutation in permutations:
        global_args = {
                'cloud_name':permutation[0],
                'load': 0,
                'load_props': 0,
                'data_type' : permutation[1],
                'background_subtract': 0,
                'recalculate_likelihoods': permutation[2],
                'bin_image': permutation[3],
                'use_weights': 0,
                'init_vel_width': permutation[4],
                'fixed_width': permutation[5],
                'use_intercept': permutation[6],
                'av_mask_threshold': permutation[7],
                'binned_data_filename_ext': '_bin',
                'likelihood_resolution': 'coarse',
                'region': permutation[8],
                'subtract_comps': permutation[9],
                'plot_diagnostics': 0,
                'use_background': permutation[10],
                'clobber_spectra': 0,
                'smooth_hi_to_co_res': 1,
                'clobber_hi_error': 0,
                'sim_hi_error': True,
                'hi_range_calc': permutation[11],
                #'num_bootstraps': 10000,
                'num_bootstraps': 10,
                'num_resid_bootstraps': 100,
                'bootstrap_fit_residuals': False,
                'calculate_median_error': False,
                'multiprocess': 1,
                'radiation_type': permutation[12],
                'rotate_cores': permutation[13],
                'vary_phi_g': permutation[14],
                }
        run_analysis = False
        if global_args['data_type'] in ('planck_lee12mask', 'lee12'):
            if global_args['cloud_name'] == 'perseus':
                run_analysis = True
        else:
            if global_args['cloud_name'] == 'california':
                if global_args['region'] is None:
                    run_analysis = True
            else:
                run_analysis = True

        if run_analysis:
            results[global_args['cloud_name']] = \
                    get_results(global_args)

            print('\n\tPlotting')
            plot_results(results[global_args['cloud_name']])

    plot_multicloud_results(results)

if __name__ == '__main__':
    main()



