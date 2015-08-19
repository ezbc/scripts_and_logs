#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import cloudpy
import matplotlib.pyplot as plt

global debugging
debugging = True

def mask_nans(arrays, return_mask=False):

    """ Masks any positions where any array in the list has a NaN.

    Parameters
    ----------
    arrays : tuple
        Tuple of arrays with same dimension.

    """

    mask = np.zeros(arrays[0].shape, dtype=bool)

    for array in arrays:
        mask[array != array] = 1

    masked_arrays = []
    for array in arrays:
        if isinstance(array, np.ndarray):
            masked_arrays.append(array[~mask])
        else:
            masked_arrays.append(array)

    if return_mask:
        return masked_arrays, mask
    else:
        return masked_arrays

def fit_model(av, nhi, av_error=None, algebraic=False, nhi_background=None):

    from lmfit import minimize, Parameters

    if nhi_background is None:
        use_background = False
    else:
        use_background = True

    # Use linear algebra?
    if algebraic:
        b = av
        A = np.array([nhi, np.ones(nhi.shape)]).T

        # weights
        if av_error is not None:
            W = 1.0 / av_error**2
        else:
            W = np.ones(av.shape)

        A = np.array([nhi, np.ones(nhi.shape)]).T

        params = np.dot(np.linalg.pinv(A), b)

    else:
        # Set parameter limits and initial guesses
        print('Fitting...')
        params = Parameters()
        params.add('dgr_cloud',
                   value=0.1,
                   min=0.0,
                   max=0.5,
                   )
        params.add('dgr_background',
                   value=0.0,
                   min=0.0,
                   max=0.5,
                   vary=use_background,
                   )
        params.add('intercept',
                   value=0.0,
                   min=-3,
                   max=3,
                   )

        #bin_edges = residuals_crop
        #counts = np.ones(residuals_crop.size - 1)

        def norm(params, av, nhi, av_error=None, nhi_background=None,):
            if nhi_background is None:
                nhi_background = np.zeros(nhi.shape)
            if av_error is None:
                av_error = np.ones(av.shape)


            model = params['dgr_cloud'] * nhi + \
                    params['dgr_background'] * nhi_background + \
                    params['intercept']

            norm = np.sum((av - model)**2 * (1.0/av_error**2)) / \
                   np.sum(1.0/av_error**2)
            #norm = np.sum((av - model)**2)
            return norm

        #print('fitting')
        # Perform the fit!
        result = minimize(norm,
                          params,
                          args=(av, nhi, av_error, nhi_background),
                          #method='lbfgsb',
                          method='nelder',
                          )

        dgr_cloud = params['dgr_cloud'].value
        dgr_background = params['dgr_background'].value
        intercept = params['intercept'].value

        if debugging:
            print('dgr = ', dgr_cloud)
            print('dgr background = ', dgr_background)
            print('intercept = ', intercept)
            plt.close(); plt.clf()
            background = dgr_background * nhi_background
            plt.plot(nhi, av - background,
                     linestyle='',
                     marker='o',
                     alpha=0.1,
                     markersize=1)
            xfit = np.linspace(0,50)
            plt.plot(xfit, dgr_cloud * xfit + intercept)
            plt.xlim(0, 22)
            plt.ylim(0, 15)
            plt.savefig('/usr/users/ezbc/Desktop/avfit.png')

        return (dgr_cloud, dgr_background, intercept)


def bootstrap_fits(av_data, nhi_image, av_error_data=None,
        nhi_image_background=None, num_bootstraps=1000):

    from astropy.stats import bootstrap

    if av_error_data is None:
        av_error_data = np.ones(av_data.size)

    # mask for nans, arrays will be 1D
    av, av_error, nhi, nhi_back = \
            mask_nans((av_data, av_error_data, nhi_image, nhi_image_background))

    # bootstrap
    for i in xrange(num_bootstraps):
        # get bootstrap indices and apply them to each dataset
        boot_indices = np.random.choice(av.size, size=av.size)
        av_boot = av[boot_indices]
        av_error_boot = av_error[boot_indices]
        nhi_boot = nhi[boot_indices]
        nhi_back_boot = nhi_back[boot_indices]

        # fit the bootstrapped data
        fit_results = fit_model(av_boot, nhi_boot, av_error=av_error_boot,
                                nhi_background=nhi_back_boot)

def run_cloud_analysis(args,
    cloud_name = 'perseus',
    region = None,
    load = False,
    data_type = 'planck'):


    from astropy.io import fits
    from myimage_analysis import calculate_nhi
    from mycoords import make_velocity_axis

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

    # define filename
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

    # Load data
    av_filename = av_filename.replace('.fits', '_bin.fits')
    av_error_filename = av_error_filename.replace('.fits', '_bin.fits')
    hi_filename = hi_filename.replace('.fits', '_bin.fits')

    av_data, av_header = fits.getdata(av_filename, header=True)
    if av_error_filename is not None:
        av_error_data, av_error_header = fits.getdata(av_error_filename,
                                                      header=True)
    else:
        av_error_data = av_error * np.ones(av_data.shape)
    hi_data, hi_header = fits.getdata(hi_filename, header=True)
    co_data, co_header = fits.getdata(co_filename, header=True)

    hi_vel_axis = make_velocity_axis(hi_header)
    nhi_image = calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=(-5,15),
                              )
    nhi_image_background = calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=(-100,-5),
                              )
    nhi_image_background += calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=(15,100),
                              )
    nhi_image[nhi_image < 0] = np.nan
    nhi_image_background[nhi_image_background < 0] = np.nan

    # Perform bootsrapping
    boot_result = bootstrap_fits(av_data,
                                 nhi_image,
                                 av_error_data=av_error_data,
                                 nhi_image_background=nhi_image_background,
                                 )

    return results

def main():

    import itertools

    results = {}

    clouds = (
              'california',
              'taurus',
              'perseus',
              )

    data_types = (
                  #'planck_lee12mask',
                  'planck',
                  #'lee12',
                  #'k09',
                  )
    recalculate_likelihoods = (
                               #True,
                               False,
                               )
    bin_image = (
                 True,
                 #False,
                 )
    init_vel_width = (#20,
                      20,
                      #200,
                      )
    fixed_width = (
                   #20,
                   'gaussfit',
                   #None,
                   #50,
                   )
    use_intercept = (
                     True,
                     #False,
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

    subtract_comps = (True,
                      #False,
                      )

    elements = (clouds, data_types, recalculate_likelihoods, bin_image,
            init_vel_width, fixed_width, use_intercept, av_mask_threshold,
            regions, subtract_comps)

    permutations = list(itertools.product(*elements))

    print('Number of permutations to run: ' + str(len(permutations)))

    #for cloud in clouds:
    for permutation in permutations:
        args = {'cloud_name':permutation[0],
                'load': 1,
                'load_props': 0,
                #'data_type': 'planck',
                #'data_type': 'k09',
                #'data_type': 'planck_lee12mask',
                #'data_type': 'lee12',
                'data_type' : permutation[1],
                'background_subtract': 0,
                'recalculate_likelihoods': permutation[2],
                'bin_image': permutation[3],
                'use_weights': 0,
                'init_vel_width': permutation[4],
                #'fixed_width': 20,
                'fixed_width': permutation[5],
                'use_intercept': permutation[6],
                'av_mask_threshold': permutation[7],
                #'av_mask_threshold': 1.2,
                'binned_data_filename_ext': '_bin',
                #'likelihood_resolution': 'fine',
                'likelihood_resolution': 'coarse',
                'region': permutation[8],
                'subtract_comps': permutation[9],
                }
        run_analysis = False
        if args['data_type'] in ('planck_lee12mask', 'lee12'):
            if args['cloud_name'] == 'perseus':
                run_analysis = True
        else:
            if args['cloud_name'] == 'california':
                if args['region'] is None:
                    run_analysis = True
            else:
                run_analysis = True

        if run_analysis:
            results[args['cloud_name']] = run_cloud_analysis(args)

if __name__ == '__main__':
    main()



