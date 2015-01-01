
''' Calculations
'''

def search_likelihoods(mesh):

    from myimage_analysis import calculate_nhi

    # calculate the likelihoodelation coefficient for each velocity range

    # Progress bar parameters
    #total = float(likelihoods.size)
    #count = 0

    try:
        velocity_center = mesh[0]
        velocity_width = mesh[1]
        dgr = mesh[2]

        velocity_range = (velocity_center - velocity_width / 2.,
                          velocity_center + velocity_width / 2.)

        nhi_image_temp, nhi_image_error = \
                calculate_nhi(cube=hi_cube,
                    velocity_axis=hi_velocity_axis,
                    velocity_range=velocity_range,
                    noise_cube=hi_noise_cube)

        # Avoid NaNs
        indices = np.where((nhi_image_temp == nhi_image_temp) & \
                           (av_image == av_image) & \
                           (nhi_image_temp > 0))

        nhi_image_likelihood = nhi_image_temp[indices]
        nhi_image_error_likelihood = nhi_image_error[indices]
        av_image_likelihood = av_image[indices]
        if type(av_image_error) != float:
            av_image_error_likelihood = np.median(av_image_error[indices])
        else:
            av_image_error_likelihood = np.median(av_image_error)

        # Create model of Av with N(HI) and DGR
        av_image_model = nhi_image_likelihood * dgr

        logL = calc_logL(av_image_model,
                         av_image_likelihood,
                         data_error=av_image_error_likelihood)

        likelihood = logL

        return likelihood
    except KeyboardInterrupt:
        raise KeyboardInterruptError

def setup_likelihood_mesh(velocity_centers, velocity_widths, dgrs):
    mesh = np.meshgrid(velocity_centers, velocity_widths, dgrs)
    mesh = np.vstack(mesh).reshape(3, -1).T

    return mesh

def reshape_likelihoods(likelihoods, velocity_centers=None,
        velocity_widths=None, dgrs=None):

    likelihoods = np.asarray(likelihoods).T.reshape(len(velocity_centers),
                                                    len(velocity_widths),
                                                    len(dgrs))

    return likelihoods

def calc_likelihood_hi_av(#hi_cube=None, hi_velocity_axis=None,
        #hi_noise_cube=None, av_image=None, av_image_error=None,
        velocity_centers=None, velocity_widths=None, return_likelihoods=True,
        dgrs=None, plot_results=True, results_filename='',
        vel_center_image=None, likelihood_filename=None, clobber=False,
        conf=0.68, contour_confs=None, multithread=False):

    '''
    Parameters
    ----------

    Returns
    -------
    hi_vel_range : tuple
        Lower and upper bound of HI velocity range in km/s which provides the
        best likelihoodelated N(HI) distribution with Av.
    likelihoods : array-like, optional
        Array of Pearson likelihoodelation coefficients likelihoodesponding to each
        permutation through the velocity centers and velocity widths.

    '''

    import numpy as np
    from scipy.stats import pearsonr
    from scipy.stats import kendalltau
    from myimage_analysis import calculate_nhi
    from scipy import signal
    from os import path
    from astropy.io import fits
    import multiprocessing
    import itertools

    # Check if likelihood grid should be derived
    if likelihood_filename is not None:
        if not path.isfile(likelihood_filename):
            perform_mle = True
            write_mle = True
        elif clobber:
            perform_mle = True
            write_mle = True
        else:
            perform_mle = False
            write_mle = False
    # If no filename provided, do not read file and do not write file
    else:
        write_mle = False
        perform_mle = True

    if perform_mle:
        if multithread:
            print('\nUsing multithreading in likelihood claculation...')
            # calculate the velocity ranges given a set of centers and widths
            velocity_ranges = np.zeros(shape=[len(velocity_centers) * \
                    len(velocity_widths),2])
            count = 0
            for i, center in enumerate(velocity_centers):
                for j, width in enumerate(velocity_widths):
                    velocity_ranges[count, 0] = center - width/2.
                    velocity_ranges[count, 1] = center + width/2.
                    count += 1

            # Set up iterable whereby each row contains the parameter values
            mesh = setup_likelihood_mesh(velocity_centers,
                                         velocity_widths,
                                         dgrs)

            # Use multiple processors to iterate through parameter
            # combinations

            p = multiprocessing.Pool()
            likelihoods = p.map(search_likelihoods, mesh)
            p.close()

            # reshape likelihoods
            likelihoods = reshape_likelihoods(likelihoods,
                                velocity_centers=velocity_centers,
                                velocity_widths=velocity_widths,
                                dgrs=dgrs)


        else:
            # calculate the likelihoodelation coefficient for each velocity
            # range
            likelihoods = np.zeros((len(velocity_centers),
                                     len(velocity_widths),
                                     len(dgrs)))

            # Progress bar parameters
            total = float(likelihoods.size)
            count = 0

            for i, velocity_center in enumerate(velocity_centers):
                for j, velocity_width in enumerate(velocity_widths):
                    # Construct N(HI) image outside of DGR loop, then apply
                    # DGRs in loop

                    if vel_center_image is None:
                        velocity_range = (velocity_center-velocity_width / 2.,
                                          velocity_center+velocity_width / 2.)
                    else:
                        velocity_range = (vel_center_image - velocity_width/2.,
                                          vel_center_image + velocity_width/2.)

                    if 0:
                        velocity_range = np.array(velocity_range)
                        import matplotlib.pyplot as plt
                        plt.close(); plt.clf()
                        plt.imshow(velocity_range[0, :, :], origin='lower')
                        plt.show()

                    nhi_image_temp = \
                            calculate_nhi(cube=hi_cube,
                                velocity_axis=hi_velocity_axis,
                                velocity_range=velocity_range,
                                )

                    # Avoid NaNs
                    indices = np.where((nhi_image_temp == nhi_image_temp)&\
                                       (av_image == av_image))

                    nhi_image_likelihood = nhi_image_temp[indices]
                    av_image_likelihood = av_image[indices]
                    if type(av_image_error) != float:
                        av_image_error_likelihood = av_image_error[indices]
                    else:
                        av_image_error_likelihood = av_image_error

                    for k, dgr in enumerate(dgrs):
                        # Create model of Av with N(HI) and DGR
                        av_image_model = nhi_image_likelihood * dgr

                        logL = calc_logL(av_image_model,
                                         av_image_likelihood,
                                         data_error=av_image_error_likelihood)

                        likelihoods[i, j, k] = logL

                        # Shows progress each 10%
                        count += 1
                        abs_step = int((total * 1)/10) or 10
                        if count and not count % abs_step:
                            print "\t{0:.0%} processed".format(count/total)

        #likelihoods /= (av_image_model.size)

        # Normalize the log likelihoods
        likelihoods -= likelihoods.max()

        # Convert to likelihoods
        likelihoods = np.exp(likelihoods)

        # Normalize the likelihoods
        likelihoods = likelihoods / np.nansum(likelihoods)

        # Write out fits file of likelihoods
        if write_mle:
            write_mle_tofits(filename=likelihood_filename,
                             velocity_centers=velocity_centers,
                             velocity_widths=velocity_widths,
                             dgrs=dgrs,
                             likelihoods=likelihoods,
                             clobber=clobber)

    # Load file of likelihoods
    elif not perform_mle:
        print('Reading likelihood grid file:')
        print(likelihood_filename)

        hdu = fits.open(likelihood_filename)
        likelihoods = hdu[0].data

        if len(velocity_centers) != likelihoods.shape[0] or \
            len(velocity_widths) != likelihoods.shape[1]:
            raise ValueError('Specified parameter grid not the same as in' + \
                    'loaded data likelihoods.')

        likelihoods = np.ma.array(likelihoods,
                mask=(likelihoods != likelihoods))

    # Define parameter resolutions
    #delta_center = velocity_centers[1] - velocity_centers[0]
    #delta_width = velocity_widths[1] - velocity_widths[0]

    # Derive marginal distributions of both centers and widths
    center_likelihood = np.sum(likelihoods, axis=(1,2)) / \
            np.sum(likelihoods)
    width_likelihood = np.sum(likelihoods, axis=(0,2)) / \
            np.sum(likelihoods)
    dgr_likelihood = np.sum(likelihoods, axis=(0,1)) / \
            np.sum(likelihoods)

    # Derive confidence intervals of parameters
    center_confint = threshold_area(velocity_centers,
                                    center_likelihood,
                                    area_fraction=conf)
    width_confint = threshold_area(velocity_widths,
                                   width_likelihood,
                                   area_fraction=conf)
    dgr_confint = threshold_area(dgrs,
                                 dgr_likelihood,
                                 area_fraction=conf)

    # Get values of best-fit model parameters
    max_loc = np.where(likelihoods == np.max(likelihoods))
    center_max = velocity_centers[max_loc[0]][0]
    width_max = velocity_widths[max_loc[1]][0]
    dgr_max = dgrs[max_loc[2]][0]

    print('\nVelocity widths = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(width_confint[0],
                                                    width_confint[2],
                                                    np.abs(width_confint[1])))
    print('\nVelocity centers = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(center_confint[0],
                                                    center_confint[2],
                                                    np.abs(center_confint[1])))
    print('\nDGRs = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} 10^20 cm^2 mag'.format(dgr_confint[0],
                                                    dgr_confint[2],
                                                    np.abs(dgr_confint[1])))

    # Write PDF
    center = center_confint[0]
    upper_lim = (center_confint[0] + width_confint[0]/2.)
    lower_lim = (center_confint[0] - width_confint[0]/2.)
    upper_lim_error = (center_confint[2]**2 + width_confint[2]**2)**0.5
    lower_lim_error = (center_confint[1]**2 + width_confint[1]**2)**0.5

    vel_range_confint = (lower_lim, upper_lim, lower_lim_error,
            upper_lim_error)

    '''
    if plot_results:
        plot_likelihoods_hist(likelihoods,
                              velocity_widths,
                              velocity_centers,
                              x_confint=width_confint,
                              y_confint=center_confint,
                              plot_axes=('widths', 'centers'),
                              show=0,
                              returnimage=False,
                              filename=results_filename + '_wc.png',
                              contour_confs=contour_confs)
        plot_likelihoods_hist(likelihoods,
                              velocity_centers,
                              dgrs,
                              x_confint=center_confint,
                              y_confint=dgr_confint,
                              plot_axes=('centers', 'dgrs'),
                              show=0,
                              returnimage=False,
                              filename=results_filename + '_cd.png',
                              contour_confs=contour_confs)
        plot_likelihoods_hist(likelihoods,
                              velocity_widths,
                              dgrs,
                              x_confint=width_confint,
                              y_confint=dgr_confint,
                              plot_axes=('widths', 'dgrs'),
                              show=0,
                              returnimage=False,
                              filename=results_filename + '_wd.png',
                              contour_confs=contour_confs)
    '''

    if not return_likelihoods:
        return vel_range_confint, dgr_confint
    else:
        return (vel_range_confint, width_confint, dgr_confint, likelihoods,
            center_likelihood, width_likelihood, dgr_likelihood, center_max,
            width_max, dgr_max)

def calc_logL(model, data, data_error=None):

    '''
    Calculates log likelihood

    http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node2.html

    '''

    if data_error is None:
        data_error = np.std(data)

    logL = -np.sum((data - model)**2 / (2 * (data_error)**2))

    return logL

def calc_model_chisq(params, av_image=None, av_image_error=None, hi_cube=None,
        hi_velocity_axis=None, hi_noise_cube=None, dgr=None):

    from myimage_analysis import calculate_nhi

    velocity_range = params['low_vel'].value, params['high_vel'].value

    nhi_image_temp, nhi_image_error = calculate_nhi(cube=hi_cube,
            velocity_axis=hi_velocity_axis,
            velocity_range=velocity_range,
            noise_cube=hi_noise_cube)

    # Select pixels with Av > 1.0 mag and Av_SNR > 5.0.
    # Av > 1.0 mag is used to avoid too low Av.
    # 1.0 mag likelihoodesponds to SNR = 1 / 0.2 ~ 5
    # (see Table 2 of Ridge et al. 2006).
    indices = np.where((nhi_image_temp == nhi_image_temp) & \
                       (av_image == av_image) & \
                       (av_image_error == av_image_error))

    nhi_image_likelihood = nhi_image_temp[indices]
    nhi_image_error_likelihood = nhi_image_error[indices]
    av_image_data = av_image[indices]
    av_image_error_data = av_image_error[indices]

    # Create model image
    av_image_model = nhi_image_likelihood * dgr
    av_image_error_model = nhi_image_error_likelihood * dgr

    chisq = np.sum((av_image_data - av_image_model)**2 / \
                   (av_image_error_data**2))

    return chisq

def fit_hi_vel_range(guesses=None, av_image=None, av_image_error=None,
        hi_cube=None, hi_velocity_axis=None, hi_noise_cube=None, dgr=None):

    from scipy.optimize import curve_fit
    from scipy import stats
    from lmfit import minimize, Parameters, report_fit, report_errors
    from pprint import pprint

    params = Parameters()
    params.add('low_vel',
               value=guesses[0],
               min=-100,
               max=100,
               vary=True)
    params.add('high_vel',
               value=guesses[1],
               min=-100,
               max=100,
               vary=True)

    result = minimize(calc_model_chisq,
                      params,
                      kws={'av_image': av_image,
                           'av_image_error': av_image_error,
                           'hi_cube': hi_cube,
                           'hi_velocity_axis': hi_velocity_axis,
                           'hi_noise_cube': hi_noise_cube,
                           'dgr': dgr},
                      #method='BFGS',
                      method='anneal',
                      #method='powell',
                      #method='SLSQP',
                      options={'gtol': 1e-6,
                               'disp': True,
                               'maxiter' : 1e9}
                      )

    report_fit(params)
    report_errors(params)
    print result.values
    print result.errorbars
    #print(result.__dict__)
    #print(dir(result))
    #print result.vars
    #print help(result)
    print result.view_items()
    print result.Bopt


    #report_fit(result)

def threshold_area(x, y, area_fraction=0.68):

    '''
    Finds the limits of a 1D array which includes a given fraction of the
    integrated data.

    Parameters
    ----------
    data : array-like
        1D array.
    area_fraction : float
        Fraction of area.

    Returns
    -------
    limits : tuple
        Lower and upper bound including fraction of area.

    '''

    import numpy as np
    from scipy.integrate import simps as integrate

    # Check if size of data
    if x.size == 1:
        return x[0], 0, 0

    # Step for lowering threshold
    step = (np.max(y) - np.median(y)) / 10000.0

    # initial threshold
    threshold = np.max(y) - step
    threshold_area = 0.0

    # area under whole function
    area = integrate(y, x)

    # Stop when the area below the threshold is greater than the max area
    while threshold_area < area * area_fraction and threshold > 0:

        threshold_indices = np.where(y > threshold)[0]

        try:
            bounds_indices = (threshold_indices[0], threshold_indices[-1])
        except IndexError:
            bounds_indices = ()

        try:
            threshold_area = integrate(y[bounds_indices[0]:bounds_indices[1]],
                                       x[bounds_indices[0]:bounds_indices[1]])
            threshold_area += threshold * (x[bounds_indices[1]] - \
                                           x[bounds_indices[0]])
        except IndexError:
            threshold_area = 0

        threshold -= step

    if threshold < 0:
        bounds_indices = (0, len(x) - 1)

    x_peak = x[y == y.max()][0]
    low_error, up_error = x_peak - x[bounds_indices[0]], \
                          x[bounds_indices[1]] - x_peak

    return (x_peak, low_error, up_error)

def write_mle_tofits(filename='', velocity_centers=None,
        velocity_widths=None, dgrs=None, likelihoods=None, clobber=False):

    from astropy.io import fits

    print('Writing likelihood grid to file:')
    print(filename)

    header = fits.Header()
    header['NAXIS'] = 3
    header['CTYPE1'] = 'CENTERS'
    header['CTYPE2'] = 'WIDTHS'
    header['CTYPE3'] = 'DGR'
    header['CRPIX1'] = 0
    header['CRPIX2'] = 0
    header['CRPIX3'] = 0
    header['CRVAL1'] = velocity_centers[0]
    header['CRVAL2'] = velocity_widths[0]
    header['CRVAL3'] = dgrs[0]
    try:
        header['CDELT1'] = velocity_centers[1] - velocity_centers[0]
    except IndexError:
        header['CDELT1'] = 1
    try:
        header['CDELT2'] = velocity_widths[1] - velocity_widths[0]
    except IndexError:
        header['CDELT2'] = 1
    try:
        header['CDELT3'] = dgrs[1] - dgrs[0]
    except IndexError:
        header['CDELT3'] = 1

    fits.writeto(filename,
                 likelihoods,
                 header,
                 clobber=clobber)

def gauss(x, width, amp, x0):
    import numpy as np

    return amp * np.exp(-(x - x0)**2 / (2 * width**2))

def get_residual_mask(residuals, resid_width_scale=3.0, plot_progress=False):

    '''

    '''

    import numpy as np
    from scipy.optimize import curve_fit

    # Fit the rising portion of the residuals
    residuals_crop = residuals[(residuals < 0) & ~np.isnan(residuals)]

    counts, bin_edges = np.histogram(np.ravel(residuals_crop),
                                     bins=100,
                                     )
    fit_params = curve_fit(gauss,
                           bin_edges[:-1],
                           counts,
                           p0=(2, np.nanmax(counts), 0),
                           maxfev=10000,
                           )[0]

    # Include only residuals within 3 sigma
    mask = residuals > resid_width_scale * np.abs(fit_params[0])

    if plot_progress:
        import matplotlib.pyplot as plt
        x_fit = np.linspace(-10,
                            10,
                            1000)

        y_fit = gauss(x_fit, *fit_params)

        counts, bin_edges = \
            np.histogram(np.ravel(residuals[~np.isnan(residuals)]),
                                     bins=1000,
                                     )

        bin_edges_ext = np.zeros(len(counts) + 1)
        counts_ext = np.zeros(len(counts) + 1)

        bin_edges_ext[0] = bin_edges[0] - (bin_edges[1] - bin_edges[0])
        bin_edges_ext[1:] = bin_edges[:-1]
        counts_ext[0] = 0
        counts_ext[1:] = counts

        counts_ext /= np.nanmax(counts_ext)
        y_fit /= np.max(y_fit)

        plt.close;plt.clf()
        plt.plot(bin_edges_ext, counts_ext, drawstyle='steps-mid')
        plt.plot(x_fit, y_fit, color='r')
        plt.xlim([np.nanmin(bin_edges_ext),4])
        plt.ylim([-0.1, 1.1])
        plt.axvline(resid_width_scale * np.abs(fit_params[0]),
                    color='k',
                    linestyle='--',
                    linewidth=3)
        plt.xlabel(r'Residual $A_V$ [mag]')
        plt.ylabel('Normalized PDF')
        plt.show()

    return mask

def main(av_data_type='planck'):

    import numpy as np
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    from astropy.io import fits
    import matplotlib.pyplot as plt

    # Check if likelihood file already written, rewrite?
    clobber = 1

    # Confidence of parameter errors
    conf = 0.68
    # Confidence of contour levels
    contour_confs = (0.68, 0.95)

    likelihood_filename = 'taurus_likelihood_{0:s}'.format(av_data_type)
    results_filename = 'taurus_likelihood_{0:s}'.format(av_data_type)

    # Threshold for converging DGR
    threshold_delta_dgr = 0.00005

    # define directory locations
    # --------------------------
    output_dir = '/home/ezbc/research/data/taurus/python_output/nhi_av/'
    figure_dir = \
        '/d/bip3/ezbc/taurus/figures/hi_velocity_range/'
    av_dir = '/home/ezbc/research/data/taurus/av/'
    hi_dir = '/home/ezbc/research/data/taurus/hi/'
    co_dir = '/home/ezbc/research/data/taurus/co/'
    core_dir = '/home/ezbc/research/data/taurus/python_output/core_properties/'
    property_dir = '/home/ezbc/research/data/taurus/python_output/'
    region_dir = '/home/ezbc/research/data/taurus/python_output/ds9_regions/'
    likelihood_dir = '/home/ezbc/research/data/taurus/python_output/nhi_av/'


    # Load data
    # ---------
    av_data, av_header = fits.getdata(av_dir + \
                                       'taurus_av_planck_5arcmin.fits',
                                       header=True)

    av_data_error, av_error_header = fits.getdata(av_dir + \
                'taurus_av_error_planck_5arcmin.fits',
            header=True)

    hi_data, hi_header = fits.getdata(hi_dir + \
                'taurus_hi_galfa_cube_regrid_planckres.fits',
            header=True)

    # make the velocity axes
    hi_vel_axis = make_velocity_axis(hi_header)

    # Velocity range over which to integrate HI
    vel_range = (0, 10)

    # Make Av model
    # -------------
    nhi_image = calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=vel_range,
                              return_nhi_error=False,
                              )

    #plt.clf(); plt.close()
    #plt.imshow(nhi_image, origin='lower')
    #plt.show()

    # Mask out nans and high-valued pixels
    mask = ((av_data > 30.0) | \
            np.isnan(av_data) | \
            np.isnan(av_data_error) | \
            (av_data_error == 0) | \
            np.isnan(nhi_image))

    # solve for DGR using linear least squares
    print('\nSolving for DGR...')

    delta_dgr = 1e10
    dgr = 1e10
    while delta_dgr > threshold_delta_dgr:
        A = np.array((np.ravel(nhi_image[~mask] / av_data_error[~mask]),))
        b = np.array((np.ravel(av_data[~mask] / av_data_error[~mask]),))
        A = np.matrix(A).T
        b = np.matrix(b).T
        #dgr = np.dot(np.linalg.pinv(A), b)
        dgr_new = (np.linalg.pinv(A) * b)[0, 0]

        # Create model with the DGR
        print('\nDGR = {0:.2} 10^20 cm^2 mag'.format(dgr))
        av_image_model = nhi_image * dgr_new

        residuals = av_data - av_image_model

        # Include only residuals which are white noise
        mask_new = get_residual_mask(residuals,
                resid_width_scale=2.0, plot_progress=0)

        # Mask non-white noise, i.e. correlated residuals.
        mask += mask_new

        npix = mask.size - np.sum(mask)
        print('Number of non-masked pixels = {0:.0f}'.format(npix))

        # Reset while loop conditions
        delta_dgr = np.abs(dgr - dgr_new)
        dgr = dgr_new

    plt.clf(); plt.close()
    nhi_image_copy = np.copy(nhi_image)
    nhi_image_copy[mask] = np.nan
    av_image_copy = np.copy(av_data)
    resid_image = av_image_copy - nhi_image_copy * dgr
    plt.imshow(resid_image, origin='lower')
    plt.title(r'$A_V$ Data - Model')
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    main()









