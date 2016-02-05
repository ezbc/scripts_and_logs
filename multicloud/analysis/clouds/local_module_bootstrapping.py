#!/usr/bin/python


import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import myimage_analysis as myia
from multiprocessing.queues import Queue
import mygeometry as myg
import scipy



def create_cloud_model(av, nhi_background, dgr_background,):

    if nhi_background is None:
        return av

    return av - dgr_background * nhi_background

def create_background_model(av, nhi_cloud, dgr_cloud):

    return av - dgr_cloud * nhi_cloud

def create_filename_base(global_args):

    # Name of diagnostic files
    if global_args['background_subtract']:
        background_name = '_backsub'
    else:
        background_name = ''

    if global_args['bin_image']:
        bin_name = '_binned'
        global_args['bin_procedure'] = 'all'
    else:
        bin_name = ''
        global_args['bin_procedure'] = 'none'
    if global_args['fixed_width'] is None:
        width_name = ''
        init_vel_width = global_args['init_vel_width']
        vel_center_gauss_fit_kwargs = None
    else:
        if global_args['fixed_width'] == 'gaussfit':
            if global_args['cloud_name'] == 'perseus':
                guesses = (28, 3, 5,
                           2, -20, 20)
                ncomps = 2
            elif global_args['cloud_name'] == 'taurus':
                guesses = (28, 3, 5,
                           5, -30, 20,
                           3, -15, 5,
                           )
                ncomps = 3
            elif global_args['cloud_name'] == 'california':
                guesses = (50, 3, 5,
                           20, -10, 10,
                           3, -45, 10,
                           #2, -20, 20,
                           )
                ncomps = 3
            vel_center_gauss_fit_kwargs = {'guesses': guesses,
                                           'ncomps': ncomps,
                                           #'width_scale': 2,
                                           }
        else:
            vel_center_gauss_fit_kwargs = None
        width_name = '_fixedwidth'
        init_vel_width = global_args['fixed_width']
    if global_args['use_weights']:
        weights_name = '_weights'
        weights_filename = av_dir + \
           global_args['cloud_name'] + '_binweights.fits'
    else:
        weights_name = ''
        weights_filename = None
    if global_args['region'] is None:
        region_name = ''
        global_args['region_name'] = global_args['cloud_name']
    else:
        region_name = '_region' + global_args['region']
        global_args['region_name'] = global_args['cloud_name'] + global_args['region']
    if global_args['av_mask_threshold'] is not None:
        avthres_name = '_avthres'
    else:
        avthres_name = ''
    if not global_args['use_intercept']:
        intercept_name = '_noint'
    else:
        intercept_name = ''
    if global_args['recalculate_likelihoods']:
        error_name = '_errorrecalc'
    else:
        error_name = ''
    if global_args['subtract_comps']:
        compsub_name = '_compsub'
    else:
        compsub_name = ''
    if global_args['use_background']:
        backdgr_name = '_backdgr'
    else:
        backdgr_name = ''
    if global_args['hi_range_calc'] == 'gaussian':
        hi_range_name = 'gaussrange'
    else:
        hi_range_name = 'stdrange'
    if global_args['rotate_cores']:
        rotate_cores_name = '_rotatedcores'
    else:
        rotate_cores_name = ''
    if global_args['vary_phi_g']:
        vary_phi_g_name = '_varyphig'
    else:
        vary_phi_g_name = ''

    filename_extension = global_args['cloud_name'] + '_' + global_args['data_type'] + \
            background_name + \
            bin_name + weights_name + \
            region_name + width_name + avthres_name + \
            intercept_name + error_name + compsub_name + backdgr_name + \
            '_' + hi_range_name + '_' + global_args['radiation_type'] + \
            rotate_cores_name + vary_phi_g_name

    return filename_extension, global_args

def mask_nans(arrays, return_mask=False):

    """ Masks any positions where any array in the list has a NaN.

    Parameters
    ----------
    arrays : tuple
        Tuple of arrays. The mask will be the shape of the first array. The
        last axes of the rest of the arrays will be masked.

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

def simulate_noise(av, av_error):

    ''' Simulates noise of Av data

    Possible noise contributions:
        + uncertainty in dust params, tau_353, Beta, T

        + CIB background - section 4.2, Planck 2013

        + background subtraction

        + variation in dust temperature, e.g. difference between radiance and
        tau_353

    '''

    #av_error = (av_error**2 + (0.07 * av)**2)**0.5

    # get bootstrap indices and apply them to each dataset
    np.random.seed()

    # empirical uncertainty from comparison with Schlegel 98
    #sigma_ = 0.003*3.1

    av_noise_sim = np.random.normal(0, scale=av_error)

    # empirical uncertainty from comparison with Schlegel 98
    #av_noise_sim += np.random.normal(0, scale=0.003*3.1)

    return av_noise_sim

def simulate_rescaling(av, scalar=1.0):

    denominator = np.random.uniform(low=1, high=scalar)
    av_rescale = av / denominator

    return av_rescale, denominator

def simulate_background_error(av, scale=1.0):

    # bias from 2MASS image, see lombardi et al. (2009), end of section 2
    av_bias = 0.2
    scale = (scale**2 + (av_bias)**2)**0.5

    background = np.random.normal(0, scale=scale)

    av_background_sim = av + background

    return av_background_sim, background

def simulate_nhi(hi_data, vel_axis, vel_range, vel_range_error,
        hi_data_error=None):

    vel_range_sim = [vel_range[0] + np.random.normal(0, scale=vel_range_error),
                     vel_range[1] + np.random.normal(0, scale=vel_range_error)]

    if hi_data_error is None:
        nhi_sim = myia.calculate_nhi(cube=hi_data,
                                  velocity_axis=vel_axis,
                                  velocity_range=vel_range_sim,
                                  )
        return nhi_sim.ravel(), vel_range_sim
    else:
        nhi_sim, nhi_sim_error = \
            myia.calculate_nhi(cube=hi_data,
                               velocity_axis=vel_axis,
                               velocity_range=vel_range_sim,
                               noise_cube=hi_data_error,
                               return_nhi_error=True,
                               )
        return nhi_sim.ravel(), nhi_sim_error.ravel(), vel_range_sim

def get_rotated_core_indices(core, mask, corename=None, iteration=None):

    import mygeometry as myg

    # Get random angle
    angle = np.random.uniform(low=-1, high=1) * 90.0 # deg

    # rotate vertices
    # ---------------
    vertices = np.array(core['poly_verts']['pixel'])
    center = core['center_pixel']
    vertices_rotated = myg.rotate_wedge(vertices,
                                        center,
                                        angle)

    # Get mask
    # --------
    core_mask = np.logical_not(myg.get_polygon_mask(mask,
                                                    vertices_rotated))

    if 0:
        import matplotlib.pyplot as plt
        plt.close(); plt.clf()
        plt.imshow(core_mask, origin='lower')
        plt.savefig('/d/bip3/ezbc/scratch/' + corename + '_mask' + \
                    str(iteration) + '.png')

    # Map the mask to the unraveled mask and get the indices
    # ------------------------------------------------------
    core_indices = np.where(core_mask == 0)[0]

    return core_indices

def bootstrap_worker(global_args, i):

    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    av = global_args['av']
    av_error = global_args['av_error']
    nhi = global_args['nhi']
    nhi_back = global_args['nhi_back']
    hi_data = global_args['hi_data']
    hi_data_error = global_args['hi_data_error']
    mask = global_args['mask']
    vel_axis = global_args['vel_axis']
    vel_range = global_args['vel_range']
    vel_range_error = global_args['vel_range_error']
    init_guesses = global_args['init_guesses']
    plot_kwargs = global_args['plot_kwargs']
    use_intercept = global_args['use_intercept']
    probabilities = global_args['probabilities']
    av_scalar = global_args['scale_kwargs']['av_scalar']
    intercept_error = global_args['scale_kwargs']['intercept_error']
    model_kwargs = global_args['ss_model_kwargs']
    rotate_cores = global_args['rotate_cores']

    #i = global_args['i']
    #np.random.seed()

    #queue = global_args['queue']


    # Create simulated data
    # -------------------------------------------------------------------------
    # add random noise
    av_sim = av + simulate_noise(av, av_error)

    # rescale the data somewhere between Planck and 2MASS:
    # rescaled Planck = Planck / beta where beta is between 1.0 and 1.4
    av_sim, av_scalar_sim = simulate_rescaling(av_sim, scalar=av_scalar)

    # remove background
    if 0:
        av_sim, av_background_sim = \
                simulate_background_error(av_sim,
                                          scale=intercept_error)
    else:
        av_background_sim = 0.0

    # calculate N(HI)
    if global_args['sim_hi_error']:
        nhi_sim, nhi_sim_error, vel_range_sim = \
            simulate_nhi(hi_data,
                         vel_axis,
                         vel_range,
                         vel_range_error,
                         hi_data_error=hi_data_error,
                         )
    else:
        nhi_sim = nhi

    if global_args['calculate_median_error']:
        error_dict = {}
        error_dict['av_sim'] = av_sim
        error_dict['nhi_sim'] = nhi_sim
    else:
        error_dict = None

    # Bootstrap data
    # -------------------------------------------------------------------------
    boot_indices = np.random.choice(av.size, size=av.size,)# p=probabilities)

    av_boot = av_sim[boot_indices]
    av_error_boot = av_error[boot_indices]
    nhi_boot = nhi_sim[boot_indices]
    nhi_error_boot = nhi_sim_error[boot_indices]

    if nhi_back is not None:
        nhi_back_boot = nhi_back[boot_indices]
    else:
        nhi_back_boot = None

    # for plotting
    plot_kwargs['bootstrap_num'] = i

    # Fit the bootstrapped data
    # -------------------------------------------------------------------------
    av_model_results = fit_av_model(av_boot,
                            nhi_boot,
                            av_error=av_error_boot,
                            nhi_error=nhi_error_boot,
                            nhi_background=nhi_back_boot,
                            init_guesses=init_guesses,
                            plot_kwargs=plot_kwargs,
                            use_intercept=use_intercept)

    # Calculate N(H2), then HI + H2 surf dens, fit relationship
    # -------------------------------------------------------------------------
    # calculate N(H2) maps
    nh2_image = calculate_nh2(nhi_image=nhi_sim,
                              av_image=av,
                              dgr=av_model_results['dgr_cloud'])
    nh2_image_error = calculate_nh2(nhi_image=nhi_sim_error,
                              av_image=av_error,
                              dgr=av_model_results['dgr_cloud'])

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi_sim,
                               sd_factor=1/1.25)
    hi_sd_image_error = calculate_sd(nhi_sim_error,
                                     sd_factor=1/1.25)

    h2_sd_image = calculate_sd(nh2_image,
                               sd_factor=1/0.625)
    h2_sd_image_error = calculate_sd(nh2_image_error,
                               sd_factor=1/0.625)

    h_sd_image = hi_sd_image + h2_sd_image
    h_sd_image_error = (hi_sd_image_error**2 + h2_sd_image_error**2)**0.5

    # Write ratio between H2 and HI
    rh2_image = h2_sd_image / hi_sd_image
    rh2_image_error = rh2_image * (h2_sd_image_error**2 / h2_sd_image**2 + \
                      h_sd_image_error**2 / h_sd_image**2)**0.5

    model_kwargs = global_args['ss_model_kwargs']['model_kwargs']
    cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']
    cores = global_args['ss_model_kwargs']['cores']

    ss_model_results = {}

    # scale phi_g by relative to the galactic DGR
    # Galactic DGR = 5.3 x 10^-22 cm^2 mag
    # our DGR in units of 10^-20 cm^2 mag
    # Galactic DGR in our units = 5.3 x 10^-22 / 10^-20 = 5.3 x 10^-2 = 0.053
    if 1:
        DGR = av_model_results['dgr_cloud']
        #print 'DGR = ', DGR
        phi_g = DGR / 0.053
        Z = DGR / 0.053
        #print 'phi_g', phi_g
        new_model_kwargs = dict(model_kwargs)
        new_model_kwargs['sternberg_params']['guesses'][2] = phi_g
        #new_model_kwargs['sternberg_params']['guesses'][1] = Z

    # cycle through each core, bootstrapping the pixels
    for core in cores_to_plot:
        if rotate_cores:
            core_indices = get_rotated_core_indices(cores[core],
                                                    mask,
                                                    corename=core,
                                                    iteration=i,
                                                    )
            #if core == 'G168.54-6.22':
            #    print np.sum(core_indices)
        else:
            # grab the indices of the core in the unraveled array
            core_indices = cores[core]['indices']

        if 0:
            assert av[core_indices] == cores[core]['test_pix']


        if 0:
            print('\n\tRegion size = ' + \
                  '{0:.0f} pix'.format(core_indices.size))

        # Bootstrap core indices
        #core_boot_indices = core_indices[index_ints]
        np.random.seed()
        core_boot_indices = np.random.choice(core_indices.size,
                                             size=core_indices.size,)
        # get bootstrapped pixels
        #h_sd_core = h_sd_image[core_boot_indices]
        #rh2_core = rh2_image[core_boot_indices]
        h_sd_core = h_sd_image[core_indices]
        h_sd_core_error = h_sd_image_error[core_indices]
        rh2_core = rh2_image[core_indices]
        rh2_core_error = rh2_image_error[core_indices]

        # mask negative ratios
        if 1:
            mask_rh2 = (rh2_core < 0) | (np.isnan(rh2_core))
            rh2_core = rh2_core[~mask_rh2]
            rh2_core_error = rh2_core_error[~mask_rh2]
            h_sd_core = h_sd_core[~mask_rh2]
            h_sd_core_error = h_sd_core_error[~mask_rh2]

        G0 = model_kwargs['krumholz_params']['G0'] + \
             np.random.normal(model_kwargs['krumholz_params']['G0_error'])
        ss_model_result = \
            fit_steady_state_models(h_sd_core,
                                    rh2_core,
                                    rh2_error=rh2_core_error,
                                    h_sd_error=h_sd_core_error,
                                    model_kwargs=model_kwargs,
                                    G0=G0,
                                    )

        if plot_kwargs['plot_diagnostics']:
            filename = plot_kwargs['figure_dir'] + \
                       'diagnostics/models/' + plot_kwargs['filename_base'] + \
                       '_rh2_vs_h_bootstrap' + \
                       '{0:03d}.png'.format(plot_kwargs['bootstrap_num'])
            plot_rh2_vs_h_diagnostic(h_sd_core,
                                     rh2_core,
                                     h_sd_error=h_sd_core_error,
                                     rh2_error=rh2_core_error,
                                     model_results=ss_model_result,
                                     filename=filename)

        #print '\npost fit phi_g:'
        #print ss_model_result['sternberg_results']['phi_g']

        # get HI transition result
        add_hi_transition_calc(ss_model_result)
        ss_model_results[core] = ss_model_result

        phi_cnm = ss_model_result['krumholz_results']['phi_cnm']
        #if phi_cnm < 0:
            #print 'core = ', core, ' phi_cnm =', phi_cnm

    # Write results
    # -------------------------------------------------------------------------
    #global_args['init_guesses'] = av_model_result

    # Write results
    mc_results = {}
    mc_results['data_params'] = {'av_background_sim': av_background_sim,
                                 'vel_range_sim': vel_range_sim,
                                 'av_scalar_sim': av_scalar_sim}
    mc_results['av_model_results'] = av_model_results
    mc_results['ss_model_results'] = ss_model_results
    mc_results['sim_images'] = error_dict

    # Plot distribution and fit
    #if plot_kwargs['plot_diagnostics']:
    if 0:
        dgr_cloud = av_model_results['dgr_cloud']
        dgr_background = av_model_results['dgr_background']
        intercept = av_model_results['intercept']

        filename = plot_kwargs['figure_dir'] + \
                   'diagnostics/av_nhi/' + plot_kwargs['filename_base'] + \
                   '_av_vs_nhi_bootstrap' + \
                   '{0:03d}.png'.format(plot_kwargs['bootstrap_num'])
        av_cloud = create_cloud_model(av_boot,
                                     nhi_back_boot,
                                     dgr_background,)
        av_background = create_background_model(av_boot,
                                     nhi_boot,
                                     dgr_cloud,)
        if nhi_back_boot is not None:
            #nhi_total = nhi_boot + nhi_back_boot
            #nhi_total = np.hstack((nhi_boot, nhi_back_boot))
            #av_boot = np.hstack((av_cloud, av_background))
            #av_images = (av_boot, av_cloud, av_background)
            av_images = (av_cloud, av_background)
            #nhi_images = (nhi_total, nhi_boot, nhi_back_boot)
            nhi_images = (nhi_boot, nhi_back_boot)
        else:
            nhi_total = nhi_boot
            av_images = (av_boot,)
            nhi_images = (nhi_total,)

        fit_params = {
                      'dgr_cloud': dgr_cloud,
                      'dgr_cloud_error': (0, 0),
                      'dgr_background': dgr_background,
                      'dgr_background_error': (0, 0),
                      'intercept': intercept,
                      'intercept_error': (0,0),
                      }

        #print('plotting')
        plot_av_vs_nhi(av_images,
                       nhi_images,
                       av_error=av_error_boot,
                       fit_params=fit_params,
                       contour_plot=plot_kwargs['av_nhi_contour'],
                       limits=plot_kwargs['av_nhi_limits'],
                       filename=filename,
                       )
    #queue.put(result)
    return mc_results

def bootstrap_worker_wrapper(args, i):

    import sys
    import traceback

    try:
        output = bootstrap_worker(args, i)
        return output
    except Exception as error:
        # capture the exception and bundle the traceback
        # in a string, then raise new exception with the string traceback
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

def bootstrap_fits(av_data, nhi_image=None, hi_data=None,
        nhi_image_error=None, hi_data_error=None, vel_axis=None,
        vel_range=None, vel_range_error=1, av_error_data=None,
        av_reference=None, nhi_image_background=None, num_bootstraps=100,
        plot_kwargs=None, scale_kwargs=None, use_intercept=True,
        sim_hi_error=False, ss_model_kwargs=None, multiprocess=True,
        rotate_cores=False, calc_median_error=False):

    import multiprocessing as mp
    import sys
    import traceback

    if av_error_data is None:
        av_error_data = np.ones(av_data.size)

    # mask for nans, arrays will be 1D
    (av, av_error, nhi, nhi_error, nhi_back, ), mask = \
        mask_nans((av_data, av_error_data, nhi_image, nhi_image_error, nhi_image_background),
                   return_mask=True)
    hi_data = hi_data[:, ~mask]
    hi_data_error = hi_data_error[:, ~mask]
    cores = ss_model_kwargs['cores']
    for core in cores:
        if cores[core]['mask'] is not None:
            cores[core]['mask_raveled'] = cores[core]['mask'][~mask]
            cores[core]['indices'] = \
                np.where(cores[core]['mask_raveled'] == 0)[0]
        else:
            cores[core]['indices'] = None

    probabilities = 1.0 / av_error**2
    probabilities /= np.nansum(probabilities)

    # for plotting
    plot_kwargs['num_bootstraps'] = num_bootstraps

    # initialize array for storing output
    boot_results = np.empty((3, num_bootstraps))
    init_guesses = [0.05, 0.05, 0.0] # dgr_cloud, dgr_background, intercept

    # Prep arguments
    global_args = {}
    global_args['av'] = av
    global_args['av_unmasked'] = av_data
    global_args['av_error'] = av_error
    global_args['nhi'] = nhi
    global_args['rotate_cores'] = rotate_cores
    global_args['mask'] = mask
    global_args['nhi_back'] = nhi_back
    global_args['init_guesses'] = init_guesses
    global_args['plot_kwargs'] = plot_kwargs
    global_args['use_intercept'] = use_intercept
    global_args['probabilities'] = probabilities
    global_args['scale_kwargs'] = scale_kwargs
    global_args['sim_hi_error'] = sim_hi_error
    global_args['ss_model_kwargs'] = ss_model_kwargs
    global_args['calculate_median_error'] = calc_median_error

    # initialize errors on av and nhi
    error_dict = {}
    error_dict['nhi_error'] =  0.0
    error_dict['av_error'] =  0.0
    error_dict['av_sim_images'] = []
    error_dict['nhi_sim_images'] = []
    global_args['error_dict'] = error_dict

    if sim_hi_error:
        global_args['hi_data'] = hi_data
        global_args['hi_data_error'] = hi_data_error
        global_args['vel_axis'] = vel_axis
        global_args['vel_range'] = vel_range
        global_args['vel_range_error'] = vel_range_error
    else:
        global_args['hi_data'] = None
        global_args['hi_data_error'] = None
        global_args['vel_axis'] = None
        global_args['vel_range'] = None
        global_args['vel_range_error'] = None

    # Prep multiprocessing
    queue = mp.Queue(10)
    pool = mp.Pool()
    processes = []

    # bootstrap
    if multiprocess:
        for i in xrange(num_bootstraps):
            processes.append(pool.apply_async(bootstrap_worker_wrapper,
                                              args=(global_args,i,)))
        pool.close()
        pool.join()

        # Get the results
        mc_results = collect_bootstrap_results(processes, ss_model_kwargs)

        for i in xrange(len(processes)):
            result = processes[i].get()
            boot_results[:, i] = result['av_model_results'].values()
    else:
        for i in xrange(num_bootstraps):
            processes.append(bootstrap_worker(global_args, i))

        mc_results = collect_bootstrap_results(processes, ss_model_kwargs,
                                               multiprocess=False)

        for i in xrange(len(processes)):
            result = processes[i]
            boot_results[:, i] = result['av_model_results'].values()

    return boot_results, mc_results

def residual_worker(global_args, core):

    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    rh2 = global_args['rh2']
    rh2_error = global_args['rh2_error']
    h_sd = global_args['h_sd']
    nboot = global_args['nboot']
    av = global_args['av']
    av_error = global_args['av_error']
    nhi = global_args['nhi']
    nhi_back = global_args['nhi_back']
    hi_data = global_args['hi_data']
    mask = global_args['mask']
    vel_axis = global_args['vel_axis']
    vel_range = global_args['vel_range']
    vel_range_error = global_args['vel_range_error']
    init_guesses = global_args['init_guesses']
    plot_kwargs = global_args['plot_kwargs']
    use_intercept = global_args['use_intercept']
    probabilities = global_args['probabilities']
    av_scalar = global_args['scale_kwargs']['av_scalar']
    intercept_error = global_args['scale_kwargs']['intercept_error']
    model_kwargs = global_args['ss_model_kwargs']
    rotate_cores = global_args['rotate_cores']

    model_kwargs = global_args['ss_model_kwargs']['model_kwargs']
    cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']
    cores = global_args['ss_model_kwargs']['cores']

    # cycle through each core, bootstrapping the pixels
    if rotate_cores:
        core_indices = get_rotated_core_indices(cores[core],
                                                mask,
                                                corename=core,
                                                iteration=i,
                                                )
        #if core == 'G168.54-6.22':
        #    print np.sum(core_indices)
    else:
        # grab the indices of the core in the unraveled array
        core_indices = cores[core]['indices']

    h_sd_core = h_sd[core_indices]
    rh2_core = rh2[core_indices]
    rh2_core_error = rh2_error[core_indices]

    # mask negative ratios
    mask_rh2 = (rh2_core < 0) | (np.isnan(rh2_core))
    rh2_core = rh2_core[~mask_rh2]
    rh2_core_error = rh2_core_error[~mask_rh2]
    h_sd_core = h_sd_core[~mask_rh2]

    ss_model_result = \
        fit_steady_state_models(h_sd_core,
                                rh2_core,
                                rh2_error=rh2_core_error,
                                model_kwargs=model_kwargs,
                                bootstrap_residuals=True,
                                nboot=nboot,
                                )



    # get HI transition result
    add_hi_transition_calc(ss_model_result)

    return ss_model_result, core

def residual_worker_wrapper(args, core):

    import sys
    import traceback

    try:
        output = residual_worker(args, core)
        return output
    except Exception as error:
        # capture the exception and bundle the traceback
        # in a string, then raise new exception with the string traceback
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

def bootstrap_residuals(av_data, nhi_image=None, hi_data=None, vel_axis=None,
        vel_range=None, vel_range_error=1, av_error_data=None,
        av_reference=None, nhi_image_background=None, num_bootstraps=100,
        plot_kwargs=None, scale_kwargs=None, use_intercept=True,
        sim_hi_error=False, ss_model_kwargs=None, multiprocess=True,
        rotate_cores=False,):

    import multiprocessing as mp
    import sys
    import traceback
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    if av_error_data is None:
        av_error_data = np.ones(av_data.size)

    # mask for nans, arrays will be 1D
    (av, av_error, nhi, nhi_back, ), mask = \
        mask_nans((av_data, av_error_data, nhi_image, nhi_image_background),
                   return_mask=True)
    hi_data = hi_data[:, ~mask]
    cores = ss_model_kwargs['cores']
    for core in cores:
        if cores[core]['mask'] is not None:
            cores[core]['mask_raveled'] = cores[core]['mask'][~mask]
            cores[core]['indices'] = \
                np.where(cores[core]['mask_raveled'] == 0)[0]
        else:
            cores[core]['indices'] = None

    probabilities = 1.0 / av_error**2
    probabilities /= np.nansum(probabilities)

    # for plotting
    plot_kwargs['num_bootstraps'] = num_bootstraps

    # initialize array for storing output
    boot_results = np.empty((3, num_bootstraps))
    init_guesses = [0.05, 0.05, 0.0] # dgr_cloud, dgr_background, intercept

    # calculate N(HI)
    nhi = myia.calculate_nhi(cube=hi_data,
                             velocity_axis=vel_axis,
                             velocity_range=vel_range,
                             )

    # Fit the data
    # -------------------------------------------------------------------------
    av_model_results = fit_av_model(av,
                            nhi,
                            av_error=av_error,
                            nhi_error=nhi_error,
                            nhi_background=nhi_back,
                            init_guesses=init_guesses,
                            plot_kwargs=plot_kwargs,
                            use_intercept=use_intercept)

    # Calculate N(H2), then HI + H2 surf dens, fit relationship
    # -------------------------------------------------------------------------
    # calculate N(H2) maps
    nh2_image = calculate_nh2(nhi_image=nhi,
                              av_image=av,
                              dgr=av_model_results['dgr_cloud'])
    nh2_image_error = calculate_nh2(nhi_image=nhi,
                              av_image=av_error,
                              dgr=av_model_results['dgr_cloud'])

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi,
                               sd_factor=1/1.25)
    hi_sd_image_error = 0.02 * hi_sd_image

    h2_sd_image = calculate_sd(nh2_image,
                               sd_factor=1/0.625)
    h2_sd_image_error = calculate_sd(nh2_image_error,
                               sd_factor=1/0.625)

    h_sd_image = hi_sd_image + h2_sd_image
    h_sd_image_error = (hi_sd_image_error**2 + h2_sd_image_error**2)**0.5

    # Write ratio between H2 and HI
    rh2_image = h2_sd_image / hi_sd_image
    rh2_image_error = rh2_image * (h2_sd_image_error**2 / h2_sd_image**2 + \
                      h_sd_image_error**2 / h_sd_image**2)**0.5

    # Prep arguments
    global_args = {}
    global_args['av'] = av
    global_args['av_unmasked'] = av_data
    global_args['av_error'] = av_error
    global_args['rh2'] = rh2_image
    global_args['rh2_error'] = rh2_image_error
    global_args['h_sd'] = h_sd_image
    global_args['nhi'] = nhi
    global_args['rotate_cores'] = rotate_cores
    global_args['mask'] = mask
    global_args['nboot'] = num_bootstraps
    global_args['nhi_back'] = nhi_back
    global_args['init_guesses'] = init_guesses
    global_args['plot_kwargs'] = plot_kwargs
    global_args['use_intercept'] = use_intercept
    global_args['probabilities'] = probabilities
    global_args['scale_kwargs'] = scale_kwargs
    global_args['sim_hi_error'] = sim_hi_error
    global_args['ss_model_kwargs'] = ss_model_kwargs
    if sim_hi_error:
        global_args['hi_data'] = hi_data
        global_args['vel_axis'] = vel_axis
        global_args['vel_range'] = vel_range
        global_args['vel_range_error'] = vel_range_error
        print 'vel_range_error', vel_range_error
    else:
        global_args['hi_data'] = None
        global_args['vel_axis'] = None
        global_args['vel_range'] = None
        global_args['vel_range_error'] = None
    cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']

    # Prep multiprocessing
    queue = mp.Queue(10)
    pool = mp.Pool()
    processes = []

    # bootstrap
    if multiprocess:
        for core in cores_to_plot:
            processes.append(pool.apply_async(residual_worker_wrapper,
                                              args=(global_args,core)))
        pool.close()
        pool.join()

        # Get the results
        resid_results = collect_residual_results(processes, ss_model_kwargs)

        for i in xrange(len(processes)):
            result = processes[i].get()
    else:
        for i in xrange(num_bootstraps):
            processes.append(residual_worker(global_args, i))

        resid_results = collect_residual_results(processes, ss_model_kwargs,
                                               multiprocess=False)

    return resid_results

def collect_residual_results(processes, ss_model_kwargs, multiprocess=True):

    empty = lambda: np.empty(len(processes))
    mc_analysis = {}
    mc_analysis['cores'] = {}

    for i in xrange(len(processes)):
        #result = queue.get())

        if multiprocess:
            result = processes[i].get()
        else:
            result = processes[i]

        model_result, core = result
        mc_analysis[core] = {}

        for model in model_result:
            mc_analysis[core][model] = {}
            model_dict = \
                mc_analysis[core][model]
            for param in model_result[model]:
                param_result = \
                    model_result[model][param]
                model_dict[param] = param_result

    return mc_analysis

def collect_bootstrap_results(processes, ss_model_kwargs, multiprocess=True):

    empty = lambda: np.empty(len(processes))
    mc_results = {}
    mc_results['av_model_results'] = {}
    mc_results['ss_model_results'] = {}
    mc_results['av_model_results']['dgr'] = empty()
    mc_results['ss_model_results']['cores'] = {}
    mc_results['sim_images'] = {}
    mc_results['sim_images']['av_sim'] = []
    mc_results['sim_images']['nhi_sim'] = []
    for core in ss_model_kwargs['cores_to_plot']:
        mc_results['ss_model_results']['cores'][core] = {}
        core_dict = mc_results['ss_model_results']['cores'][core]
        core_dict['krumholz_results'] = \
                {'phi_cnm': empty(),
                 'Z': empty(),
                 'phi_mol': empty(),
                 'hi_transition': empty(),
                 }
        core_dict['sternberg_results'] = \
                {'alphaG': empty(),
                 'Z': empty(),
                 'phi_g': empty(),
                 'hi_transition': empty(),
                 }
    mc_results['data_params'] = \
            {'av_background_sim': empty(),
             'vel_range_sim': np.empty((len(processes), 2)),
             'av_scalar_sim': empty()}

    mc_results['sim_images'] = {'av_sim': [],
                                'nhi_sim': []}

    for i in xrange(len(processes)):
        #result = queue.get())

        if multiprocess:
            result = processes[i].get()
        else:
            result = processes[i]

        for data_param in mc_results['data_params']:
            mc_results['data_params'][data_param][i] = \
                result['data_params'][data_param]
        mc_results['av_model_results']['dgr'][i] = \
            result['av_model_results']['dgr_cloud']

        if result['sim_images'] is not None:
            mc_results['sim_images']['av_sim']\
                    .append(result['sim_images']['av_sim'])
            mc_results['sim_images']['nhi_sim']\
                    .append(result['sim_images']['nhi_sim'])
        else:
            mc_results['sim_images'] = None

        for core in result['ss_model_results']:
            for model in result['ss_model_results'][core]:
                for param in result['ss_model_results'][core][model]:
                    param_result = \
                        result['ss_model_results'][core][model][param]
                    model_dict = \
                        mc_results['ss_model_results']['cores'][core][model]
                    model_dict[param][i] = param_result



    return mc_results

def fit_av_with_refav(av_data, av_reference, av_error_data):

    import scipy as sp

    # Mask AV for Av > 10 mag. The 2MASS Av will saturate at high Av thus is
    # unreliable.

    nan_mask = (np.isnan(av_reference) | \
                np.isnan(av_data) | \
                np.isnan(av_error_data) | \
                (av_reference > 10.0))
    p, V = np.polyfit(av_reference[~nan_mask], av_data[~nan_mask], deg=1,
                   #w=1.0/av_error_data[~nan_mask]**2,
                   cov=True,
                   )
    av_scalar, intercept = p
    intercept_error = V[1, 1]

    # Perform residual bootstrapping to get errors on intercept
    bootindex = sp.random.random_integers
    nboot = 1000
    x_fit = av_reference[~nan_mask]
    y_fit = p[0] * x_fit  + p[1]
    residuals = av_data[~nan_mask] - y_fit
    fits = np.empty((2, nboot))

    weights = np.abs(1.0 / av_error_data[~nan_mask])
    weights /= np.nansum(weights)

    for i in xrange(nboot): # loop over n bootstrap samples from the resids
        if 0:
            boot_indices = bootindex(0,
                                     len(residuals)-1,
                                     len(residuals))
            residuals_bootstrapped = residuals[boot_indices]
            fits[:, i] = sp.polyfit(av_reference[~nan_mask],
                                    y_fit + residuals_bootstrapped,
                                    deg=1)
        else:
            boot_indices = bootindex(0,
                                     len(x_fit)-1,
                                     len(x_fit))
            x = x_fit[boot_indices] + np.random.normal(0, 0.4,
                                                       size=x_fit.size)
            y = av_data[~nan_mask][boot_indices] + \
                    np.random.normal(0,
                            av_error_data[~nan_mask][boot_indices])

            fits[:, i] = np.polyfit(x,
                                    y,
                                    deg=1,)

    intercept_error = np.std(fits[1])
    av_scalar_error = np.std(fits[0])

    return av_scalar, av_scalar_error, intercept, intercept_error

def scale_av_with_refav(av_data, av_reference, av_error_data, perform_mc=1):

    import scipy as sp
    import mystats


    if not perform_mc:
        av_scalar, av_scalar_error, intercept, intercept_error = \
                fit_av_with_refav(av_data, av_reference, av_error_data)
    else:
        if 1:
            av_scalar, av_scalar_error, intercept, intercept_error = \
                    fit_av_with_refav(av_data, av_reference, av_error_data)

            print 'av_scaling without mc:'
            print av_scalar, av_scalar_error, intercept, intercept_error

        N_mc = 10
        mc_array = np.empty((N_mc, 4))
        for i in xrange(N_mc):
            av_data_sim = av_data + np.random.normal(0.2,
                                                     size=av_data.size)
            mc_array[i] = \
                fit_av_with_refav(av_data_sim, av_reference, av_error_data)

            #if i < 4:
            if 1:
                c_cycle = ['c', 'b', 'm', 'y', 'r']
                import matplotlib.pyplot as plt
                c_cycle = plt.cm.copper(np.linspace(0,1,N_mc))
                plt.scatter(av_reference,
                            av_data_sim,
                            alpha=0.01,
                            marker='o',
                            color=c_cycle[i]
                            )
                x_fit = np.linspace(-10, 100, 1000)
                y_fit = mc_array[i, 0] * x_fit + mc_array[i, 2]
                plt.plot(x_fit,
                         y_fit,
                         linestyle='-',
                         color='k',
                         )
                plt.xlabel('2MASS Av')
                plt.ylabel('Planck Av')
                plt.xlim(-1,20)
                plt.ylim(-1,20)

                plt.savefig('/d/bip3/ezbc/multicloud/figures/av_scaling/' + \
                            'av_scaling_' + str(i) + '.png')

        av_scalar, av_scalar_error = mystats.calc_cdf_error(mc_array[:,0])
        intercept, intercept_error = mystats.calc_cdf_error(mc_array[:,2])

        if 1:
            import matplotlib.pyplot as plt
            plt.close(); plt.clf();
            for i in xrange(N_mc):
                c_cycle = ['c', 'b', 'm', 'y', 'r']
                x_fit = np.linspace(-10, 100, 1000)
                y_fit = mc_array[i, 0] * x_fit + mc_array[i, 2]
                plt.plot(x_fit,
                         y_fit,
                         linestyle='-',
                         color='k',
                         alpha=0.3,
                         )
                plt.xlabel('2MASS Av')
                plt.ylabel('Planck Av')
                plt.xlim(-1,20)
                plt.ylim(-1,20)

            plt.savefig('/d/bip3/ezbc/multicloud/figures/av_scaling/' + \
                        'av_slopes.png')

        print 'av_scaling with mc:'
        print av_scalar, av_scalar_error, intercept, intercept_error

    kwargs = {}
    kwargs['av_scalar'] = av_scalar
    kwargs['av_scalar_error'] = av_scalar_error
    kwargs['intercept'] = intercept
    kwargs['intercept_error'] = intercept_error

    #print 'Reference Av characteristcs:'
    #print '\t', av_scalar, av_scalar_error
    #print '\t', intercept, intercept_error
    #print ''

    return kwargs

def calc_hi_vel_range(hi_spectrum, hi_vel_axis, gauss_fit_kwargs,
        width_scale=2, co_spectrum=None, co_vel_axis=None, ncomps=1,):

    '''
    Discussion of HI width error from Min:
    (3) For the uncertainty in the HI velocity range, you take the absolute
    difference between the Gaussian HI component center and the median CO peak
    velocity.  In fact, this is sort of the minimum uncertainty we expect.  The
    uncertainty mostly comes from our adoption of +/- 2 sigma(HI) to determine
    the HI velocity range.  In this sense, the maximum uncertainty we expect
    would be (HI velocity range we determine by using +/- 2 sigma(HI) as
    thresholds) - (HI velocity range over which CO emission appears).
    Therefore, a more realistic estimate for the uncertainty in the HI velocity
    range would be, e.g., ("minimum" error + "maximum" error) / 2.  This
    estimate would be a better measure in particular for Taurus, where CO
    emission appears over a relatively small velocity range.

    '''

    from scipy.stats import nanmedian
    from myfitting import fit_gaussians

    co_width = np.copy(gauss_fit_kwargs['co_width'])
    del(gauss_fit_kwargs['co_width'])

    hi_fits = fit_gaussians(hi_vel_axis,
            hi_spectrum, **gauss_fit_kwargs)

    # use either the gaussian closest to the CO peak, or the tallest gaussian
    if co_spectrum is not None:
        co_peak_vel = co_vel_axis[co_spectrum == np.nanmax(co_spectrum)][0]
        vel_diffs_center = []
        vel_diffs_edge = []
        for i, param in enumerate(hi_fits[2]):
            vel_center = param[1]
            width = param[2]
            vel_center_diff = np.abs(vel_center - co_peak_vel)
            vel_diffs_center.append(vel_center_diff)
            #vel_diffs_edge.append(vel_center_diff + 2.0 * width)
            vel_diffs_edge.append(4.0 * width)

        # get the velocity range
        vel_diffs_center = np.asarray(vel_diffs_center)
        vel_diffs_edge = np.asarray(vel_diffs_edge)

        sort_indices = np.argsort(vel_diffs_center)
        cloud_comp_num = np.asarray(sort_indices[:ncomps])
        velocity_range = [np.inf, -np.inf]
        for i in cloud_comp_num:
            # get the component
            param = hi_fits[2][i]
            vel_center = param[1]
            width = param[2]

            # set absolute bounds if the component bounds extend beyond
            upper_vel = vel_center + width * width_scale
            lower_vel = vel_center - width * width_scale
            if upper_vel > velocity_range[1]:
                velocity_range[1] = upper_vel
            if lower_vel < velocity_range[0]:
                velocity_range[0] = lower_vel

        # the offset between the co and fitted gaussians will be the HI error
        hi_width_error_min = np.max(vel_diffs_center[cloud_comp_num])

        # max error
        hi_width_error_max = np.max(vel_diffs_edge[cloud_comp_num]) - \
                             co_width
                             #gauss_fit_kwargs['co_width']
        hi_width_error_max = np.abs(velocity_range[1] - velocity_range[0]) - \
                             co_width

        hi_width_error = (hi_width_error_min + hi_width_error_max) / 2.0

        #print 'min error', hi_width_error_min, 'km/s'
        #print 'max error', hi_width_error_max, 'km/s'
        #print 'avg error', hi_width_error, 'km/s'

    else:
        amp_max = -np.Inf
        for i, param in enumerate(hi_fits[2]):
            if param[0] > amp_max:
                amp_max = param[0]
                vel_center = param[1]
                width = param[2] * 4
                cloud_comp_num = i

        velocity_range = [vel_center - width * width_scale,
                          vel_center + width * width_scale]

    return velocity_range, hi_fits, cloud_comp_num, hi_width_error

def get_gauss_fit_kwargs(global_args):
    if global_args['cloud_name'] == 'perseus':
        guesses = (28, 3, 5,
                   2, -20, 20)
        ncomps = 2
        ncomps_in_cloud = 1
        co_width = np.abs(10.0 - (-1.0))
    elif global_args['cloud_name'] == 'taurus':
        guesses = (35, 3, 5,
                   5, -15, 20,
                   #3, -2, 2,
                   )
        ncomps = 2
        ncomps_in_cloud = 1
        co_width = np.abs(8.0 - (4.0))
    elif global_args['cloud_name'] == 'california':
        guesses = (50, 3, 5,
                   10, -3, 3,
                   12, -10, 10,
                   3, -30, 10,
                   #2, -20, 20,
                   )
        ncomps = 4
        ncomps_in_cloud = 2
        co_width = np.abs(8.0 - (-4.0))
    gauss_fit_kwargs = {'guesses': guesses,
                        'ncomps': ncomps,
                        'co_width': co_width,
                        #'width_scale': 2,
                        }

    return gauss_fit_kwargs, ncomps_in_cloud

def calc_param_errors(results_dict):

    import mystats

    boot_result = results_dict['boot_result']
    global_args = results_dict['global_args']

    dgr_cloud, dgr_background, intercept = np.mean(boot_result, axis=1)
    dgr_cloud_error, dgr_background_error, intercept_error = \
                np.std(boot_result, axis=1)

    # Calculate conf interval
    dgrs = boot_result[0]
    dgr_cloud, dgr_cloud_error = mystats.calc_cdf_error(dgrs,
                                                        alpha=0.32)

    if global_args['use_background']:
        dgrs = boot_result[1]
        dgr_background_error = (mean - conf_int_a[0], conf_int_a[1] - mean)
    else:
        dgr_background_error = (0.0, 0.0)
        dgr_background = 0.0

    if global_args['use_intercept']:
        intercepts = boot_result[2]
        intercept, intercept_error = mystats.calc_cdf_error(intercepts,
                                                            alpha=0.32)
    else:
        intercept_error = (0.0, 0.0)
        intercept = 0.0

    fit_params = {
                  'dgr_cloud': dgr_cloud,
                  'dgr_cloud_error': dgr_cloud_error,#, dgr_cloud_error),
                  'dgr_background': dgr_background,
                  'dgr_background_error': dgr_background_error,
                  'intercept': intercept,
                  'intercept_error': intercept_error,
                  }

    return fit_params

def load_results(filename, load_fits=True):

    import pickle
    from astropy.io import fits

    with open(filename, 'rb') as input:
        results = pickle.load(input)
    input.close()

    if load_fits:
        results['data']['av_data'], results['data']['av_header'] = \
                fits.getdata(results['filenames']['av_filename'],
                             header=True)
        results['data']['av_error_data'], results['data']['av_error_header'] = \
                fits.getdata(results['filenames']['av_error_filename'],
                             header=True)
        results['data']['av_data_ref'] = \
                fits.getdata(results['filenames']['av_ref_filename'])
        results['data']['hi_data'], results['data']['hi_header'] = \
                fits.getdata(results['filenames']['hi_filename'], header=True)
        results['data']['hi_data_error'] = \
                fits.getdata(results['filenames']['hi_error_filename'])
        results['data']['co_data'], results['data']['co_header'] = \
                fits.getdata(results['filenames']['co_filename'], header=True)

    return results

def save_results(results_dict, filename, write_fits=False):

    import pickle
    from astropy.io import fits

    if not write_fits:
        results_dict['data']['av_data'] = None
        results_dict['data']['av_error_data'] = None
        results_dict['data']['av_data_ref'] = None
        results_dict['data']['hi_data'] = None
        results_dict['data']['hi_data_error'] = None
        results_dict['data']['co_data'] = None

    with open(filename, 'wb') as output:
        pickle.dump(results_dict, output)
    output.close()

def get_model_fit_kwargs(cloud_name, vary_phi_g=False):

    '''

    '''
    vary_alphaG = True # Vary alphaG in S+14 fit?
    vary_Z = False # Vary metallicity in S+14 fit?
    vary_phi_g = vary_phi_g # Vary phi_g in S+14 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    guesses=[10.0, 1.0, 1] # Guesses for (alphaG, Z, phi_g)
    h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

    # Monte carlo results file bases
    results_filename = '/d/bip3/ezbc/multicloud/' + \
                       '/data/python_output/' + \
                       'monte_carlo_results/' + \
                       cloud_name + '_mc_results_' + \
                       'planck' + '_'

    sternberg_params = {}
    sternberg_params['param_vary'] = [vary_alphaG, vary_Z, vary_phi_g]
    sternberg_params['error_method'] = error_method
    sternberg_params['alpha'] = alpha
    sternberg_params['guesses'] = guesses
    sternberg_params['h_sd_fit_range'] = h_sd_fit_range
    sternberg_params['results_filename'] = results_filename
    sternberg_params['parameters'] = ['alphaG', 'Z', 'phi_g']

    # Krumholz Parameters
    # --------------------
    vary_phi_cnm = True # Vary phi_cnm in K+09 fit?
    vary_Z = False # Vary metallicity in K+09 fit?
    vary_phi_mol = False # Vary phi_mol in K+09 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    guesses=[8.0, 1.0, 10.0] # Guesses for (phi_cnm, Z, phi_mol)
    h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

    krumholz_params = {}
    krumholz_params['param_vary'] = [vary_phi_cnm, vary_Z, vary_phi_mol]
    krumholz_params['error_method'] = error_method
    krumholz_params['alpha'] = alpha
    krumholz_params['guesses'] = guesses
    krumholz_params['h_sd_fit_range'] = h_sd_fit_range
    krumholz_params['results_filename'] = results_filename
    krumholz_params['parameters'] = ['phi_cnm', 'Z', 'phi_mol']
    if cloud_name == 'taurus':
        G0 = 0.6
        G0_error = 0.1
    elif cloud_name == 'california':
        G0 = 1.0
        G0_error = 0.2
    elif cloud_name == 'perseus':
        G0 = 0.7
        G0_error = 0.2
    krumholz_params['G0'] = G0
    krumholz_params['G0_error'] = G0_error

    results = {}
    results['results_filename'] = results_filename
    results['krumholz_params'] = krumholz_params
    results['sternberg_params'] = sternberg_params

    return results

def add_coldens_images(data_products, mc_analysis, mc_results):

    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    nhi = data_products['nhi']
    nhi_error = data_products['nhi_error']
    av = data_products['av']
    av_error = data_products['av_error']
    dgr = mc_analysis['dgr']
    dgr_error = np.mean(np.abs(mc_analysis['dgr_error']))

    nh2_image = calculate_nh2(nhi_image=nhi,
                              av_image=av,
                              dgr=dgr)

    # nh2 = (av * dgr - nhi) / 2
    # nh2_error = ((nh2 * ((av_error / av)**2 + (dgr_error / dgr)**2)**0.5)**2 \
    #              - nhi_error**2.0)**0.5 / 2
    av_comp_error = nh2_image * ((av_error / av)**2 + (dgr_error / dgr)**2)**0.5
    nh2_image_error = (av_comp_error**2 + nhi_error**2)**0.5 / 2.0

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi,
                               sd_factor=1/1.25)
    hi_sd_image_error = calculate_sd(nhi_error,
                                     sd_factor=1/1.25)

    h2_sd_image = calculate_sd(nh2_image,
                               sd_factor=1/0.625)
    h2_sd_image_error = calculate_sd(nh2_image_error,
                               sd_factor=1/0.625)

    h_sd_image = hi_sd_image + h2_sd_image
    h_sd_image_error = (hi_sd_image_error**2 + h2_sd_image_error**2)**0.5

    # Write ratio between H2 and HI
    rh2_image = h2_sd_image / hi_sd_image
    rh2_image_error = rh2_image * (h2_sd_image_error**2 / h2_sd_image**2 + \
                      h_sd_image_error**2 / h_sd_image**2)**0.5

    data_products['nh2'] = nh2_image
    data_products['h2_sd'] = h2_sd_image
    data_products['hi_sd'] = hi_sd_image
    data_products['h_sd'] = h_sd_image
    data_products['rh2'] = rh2_image
    data_products['nh2_error'] = nh2_image_error
    data_products['h2_sd_error'] = h2_sd_image_error
    data_products['hi_sd_error'] = hi_sd_image_error
    data_products['h_sd_error'] = h_sd_image_error
    data_products['rh2_error'] = rh2_image_error

    # Median sim errors
    # --------------------------------------------------------------------------
    if mc_results['sim_images'] is not None:
        av_error = np.median(np.nanstd(mc_results['sim_images']['av_sim']))
        nhi_error = np.median(np.nanstd(mc_results['sim_images']['nhi_sim']))

        # nh2 = (av * dgr - nhi) / 2
        # nh2_error = ((nh2 * ((av_error / av)**2 + (dgr_error / dgr)**2)**0.5)**2 \
        #              - nhi_error**2.0)**0.5 / 2
        av_comp_error = nh2_image * ((av_error / av)**2 + (dgr_error / dgr)**2)**0.
        nh2_image_error = (av_comp_error**2 + nhi_error**2)**0.5 / 2.0

        # convert to column density to surface density
        hi_sd_image = calculate_sd(nhi,
                                   sd_factor=1/1.25)
        hi_sd_image_error = calculate_sd(nhi_error,
                                         sd_factor=1/1.25)

        h2_sd_image = calculate_sd(nh2_image,
                                   sd_factor=1/0.625)
        h2_sd_image_error = calculate_sd(nh2_image_error,
                                   sd_factor=1/0.625)

        h_sd_image = hi_sd_image + h2_sd_image
        h_sd_image_error = (hi_sd_image_error**2 + h2_sd_image_error**2)**0.5

        # Write ratio between H2 and HI
        rh2_image = h2_sd_image / hi_sd_image
        rh2_image_error = rh2_image * (h2_sd_image_error**2 / h2_sd_image**2 + \
                          h_sd_image_error**2 / h_sd_image**2)**0.5

        data_products['h_sd_median_error'] = \
            scipy.stats.nanmedian(h_sd_image_error.ravel())
        data_products['hi_sd_median_error'] = \
                scipy.stats.nanmedian(hi_sd_image_error.ravel())
        data_products['rh2_median_error'] = \
                scipy.stats.nanmedian(rh2_image_error.ravel())
    else:
        data_products['h_sd_median_error'] = None
        data_products['hi_sd_median_error'] = None
        data_products['rh2_median_error'] = None

def calc_mc_analysis(mc_results, resid_results, data_products):

    import mystats

    mc_analysis = {}
    core_analysis = {}

    # Calculate conf intervals on parameters
    # -------------------------------------------------------------------------
    # DGR
    dgrs = mc_results['av_model_results']['dgr']
    dgr, dgr_error = mystats.calc_cdf_error(dgrs,
                                            alpha=0.32)

    if 1:
        import matplotlib.pyplot as plt
        import myplotting as myplt
        plt.clf(); plt.close()
        myplt.plot_cdf(dgrs)
        plt.axvline(dgr, linewidth=2)
        plt.axvline(dgr - dgr_error[0], linewidth=1, linestyle='--')
        plt.axvline(dgr_error[1] + dgr, linewidth=1, linestyle='--')
        plt.savefig('/d/bip3/ezbc/scratch/dgr_cdf.png')

    # Model params fit for each core:
    for core in mc_results['ss_model_results']['cores']:
        core_analysis[core] = {}
        for model in mc_results['ss_model_results']['cores'][core]:
            params = {}
            core_analysis[core][model] = {}
            for param in mc_results['ss_model_results']['cores'][core][model]:
                param_results = \
                    mc_results['ss_model_results']['cores'][core][model][param]
                param_value, param_error = \
                    mystats.calc_cdf_error(param_results,
                                           alpha=0.32)

                # add fit error
                if param != 'hi_transition' and resid_results is not None:
                    fit_error = \
                        np.array(resid_results[core][model][param + '_error'])

                    param_error = np.array(param_error)
                    #print ''
                    #print 'before', param, param_error / param_value
                    #print 'fit_error', fit_error / param_value
                    #print 'mc error', param_error / param_value
                    if 0:
                        param_error = np.sqrt(fit_error**2 + \
                                              np.array(param_error)**2)
                    #print 'after', param, param_error / param_value


                core_analysis[core][model][param] = param_value
                core_analysis[core][model][param + '_error'] = param_error

                params[param] = param_value

                #if param == 'phi_g':
                    #print '\nanalysis'
                    #print param_value
                    #print param_error


                if 0:
                    if param == 'phi_g' or param == 'alphaG':
                        import matplotlib.pyplot as plt
                        import mystats
                        plt.close(); plt.clf()
                        cdf = mystats.calc_cdf(param_results)
                        plt.plot(np.sort(param_results), cdf)
                        plt.xlabel(param)
                        plt.savefig('/d/bip3/ezbc/scratch/'+core+\
                                    '_' + param + '_cdf.png')

                    print('core ' + core + ': ' + param  + \
                          ' = {0:.2f}'.format(param_value) + \
                          '+{0:.2f} / -{1:.2f}'.format(*param_error))

            if 'sternberg' in model:
                model_fits = calc_sternberg(params,
                                          h_sd_extent=(0.001, 1000),
                                          return_fractions=False,
                                          return_hisd=True)
            elif 'krumholz' in model:
                model_fits = calc_krumholz(params,
                                          h_sd_extent=(0.001, 1000),
                                          return_fractions=False,
                                          return_hisd=True)

            core_analysis[core][model]['rh2_fit'] = model_fits[0]
            core_analysis[core][model]['hsd_fit'] = model_fits[1]
            core_analysis[core][model]['hisd_fit'] = model_fits[2]

    mc_analysis = {'dgr': dgr,
                   'dgr_error': dgr_error,
                   'cores': core_analysis,
                   }

    # comment

    return mc_analysis

def refit_data(h_sd, rh2, h_sd_error=None, rh2_error=None, model_kwargs=None):

    data_array = h_sd, rh2, h_sd_error, rh2_error
    h_sd, rh2, h_sd_error, rh2_error = mask_nans(data_array)

    G0 = model_kwargs['krumholz_params']['G0'] + \
         np.random.normal(model_kwargs['krumholz_params']['G0_error'])
    ss_model_result = \
        fit_steady_state_models(h_sd.ravel(),
                                rh2.ravel(),
                                rh2_error=rh2_error.ravel(),
                                h_sd_error=h_sd_error.ravel(),
                                model_kwargs=model_kwargs,
                                G0=G0,
                                )
    fitted_models = {}
    for model in ss_model_result:
        params = {}
        for param in ss_model_result[model]:
            params[param] = ss_model_result[model][param]

        if 'sternberg' in model:
            model_fits = calc_sternberg(params,
                                      h_sd_extent=(0.001, 200),
                                      return_fractions=False,
                                      return_hisd=True)
        elif 'krumholz' in model:
            model_fits = calc_krumholz(params,
                                      h_sd_extent=(0.001, 200),
                                      return_fractions=False,
                                      return_hisd=True)

        fitted_models[model] = {}
        fits = fitted_models[model]
        fits['rh2'] = model_fits[0]
        fits['h_sd'] = model_fits[1]
        fits['hi_sd'] = model_fits[2]

    return fitted_models

def add_results_analysis(results_dict):

    if 0:
        for core in results_dict['mc_analysis']['cores']:
            for key in results_dict['mc_analysis']['cores'][core]:
                print results_dict['mc_analysis']['cores'][core][key]

    # calculate statistics of bootstrapped model values
    results_dict['mc_analysis'] = calc_mc_analysis(results_dict['mc_results'],
                                            results_dict['resid_mc_results'],
                                            results_dict['data_products'])

    # derive N(H2), N(H) etc...
    add_coldens_images(results_dict['data_products'],
                       results_dict['mc_analysis'],
                       results_dict['mc_results'])


