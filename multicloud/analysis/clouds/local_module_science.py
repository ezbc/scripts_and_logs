#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import myimage_analysis as myia
from multiprocessing.queues import Queue
import mygeometry as myg
import scipy
import pickle

def scale_dust_areas(DGR, model_kwargs):

    ''' Scales the K+09 and S+14 dust cross sections.
    '''

    phi_g = DGR / 0.053
    sigma_d = DGR / 0.053 * 1.9
    #phi_g = DGR / 0.053 / 1.9
    #sigma_d = DGR / 0.053
    new_model_kwargs = dict(model_kwargs)
    new_model_kwargs['sternberg_params']['guesses'][2] = phi_g
    new_model_kwargs['krumholz_params']['guesses'][2] = sigma_d

    return new_model_kwargs

def fit_steady_state_models(h_sd, rh2, model_kwargs, rh2_error=None,
        h_sd_error=None, bootstrap_residuals=False, nboot=100, G0=1.0,
        odr_fit=False,):

    from myscience import bialy16

    # Fit R_H2
    #---------
    sternberg_params = model_kwargs['sternberg_params']
    sternberg_results = {}
    krumholz_params = model_kwargs['krumholz_params']
    krumholz_results = {}

    # Fit to sternberg model
    if rh2.size > 3:
        result = \
            fit_sternberg(h_sd,
                          rh2,
                          guesses=sternberg_params['guesses'],
                          vary=sternberg_params['param_vary'],
                          radiation_type=sternberg_params['radiation_type'],
                          bootstrap_residuals=bootstrap_residuals,
                          nboot=nboot,
                          rh2_error=np.abs(rh2_error),
                          h_sd_error=np.abs(h_sd_error),
                          odr_fit=odr_fit,
                          )



        if bootstrap_residuals:
            alphaG, alphaG_error, Z_s14, Z_s14_error, phi_g, phi_g_error = \
                    result
            sternberg_results['alphaG_error'] = alphaG_error
            sternberg_results['Z_error'] = Z_s14_error
            sternberg_results['phi_g_error'] = phi_g_error
        else:
            alphaG, Z_s14, phi_g = result
            alphaG_error, Z_s14_error, phi_g_error = 3*[np.nan]


        # Fit to krumholz model
        result = \
            fit_krumholz(h_sd,
                         rh2,
                         guesses=krumholz_params['guesses'],
                         vary=krumholz_params['param_vary'],
                         bootstrap_residuals=bootstrap_residuals,
                         nboot=nboot,
                         rh2_error=np.abs(rh2_error),
                         h_sd_error=np.abs(h_sd_error),
                         odr_fit=odr_fit,
                         )


        if bootstrap_residuals:
            phi_cnm, phi_cnm_error,Z_k09, Z_k09_error,sigma_d, sigma_d_error= \
                    result
            krumholz_results['phi_cnm_error'] = phi_cnm_error
            krumholz_results['Z_error'] = Z_k09_error
            krumholz_results['sigma_d_error'] = sigma_d_error
        else:
            phi_cnm, Z_k09, sigma_d = result
            phi_cnm_error, Z_k09_error, sigma_d_error = 3*[np.nan]
    else:
        alphaG, Z_s14, phi_g, phi_cnm, Z_k09, sigma_d = 6 * [np.nan]
        alphaG_error, Z_s14_error, phi_g_error, \
        phi_cnm_error, Z_k09_error, sigma_d_error = 6*[np.nan]

    # keep results
    sternberg_results['alphaG'] = alphaG
    sternberg_results['Z'] = Z_s14
    sternberg_results['phi_g'] = phi_g

    # sigma_g = 1.9 * 10^-21 phi_g * Z cm^2
    # sigma_g,gal is relative to galactic, and so is phi_g in this case, so
    # sigma_g,gal propto phi_g.
    sigma_g = phi_g
    sternberg_results['sf_threshold'] = bialy16.calc_sf_threshold(alphaG,
                                                                  sigma_g)

    # keep results
    krumholz_results['phi_cnm'] = phi_cnm
    krumholz_results['Z'] = Z_k09
    krumholz_results['sigma_d'] = sigma_d


    #print('sigma_d, phi_g', sigma_d, phi_g)
    #print('phi_cnm, alphaG', phi_cnm, alphaG)

    # see eq 6 of sternberg+09
    # alphaG is the number density of the CNM over the minimum number
    # density required for pressure balance
    # the lower alphaG values than for taurus mean that taurus
    # has a more diffuse CNM

    # By fitting the model to the observed R_H2 vs total H, you
    # basically constrained psi in Equation (35) of sternberg+09.  This
    # means that you can calculate f_H2 for a given total hydrogen
    # surface density.  In this case, H2 surface density = f_H2 *
    # total hydrogen surface density HI surface density = (1 - f_HI) *
    # total hydrogen surface density

    results = {}
    results['sternberg_results'] = sternberg_results
    results['krumholz_results'] = krumholz_results

    return results

def fit_krumholz(h_sd, rh2, guesses=[10.0, 1.0, 1.0], h_sd_error=None,
        rh2_error=None, verbose=False, vary=[True, True, True],
        bootstrap_residuals=False, nboot=100, G0=1.0, odr_fit=False):

    '''
    Parameters
    ----------
    h_sd : array-like
        Hydrogen surface density in units of solar mass per parsec**2
    rh2 : array-like
        Ratio between molecular and atomic hydrogen masses.
    guesses : None, scalar, or M-length sequence.
        Initial guess for the parameters. See scipy.optimize.curve_fit.
    rh2_error : bool
        Error in rh2 parameter. Calculates a more accurate chi^2 statistic
    G0 : float
        radiation field

    Returns
    -------
    rh2_fit_params : array-like, optional
        Model parameter fits.

    '''

    from scipy.optimize import curve_fit
    from scipy import stats
    from lmfit import minimize, Parameters, report_fit
    from myscience import krumholz09 as k09
    import mystats
    import scipy.odr as odr

    if not odr_fit:
        h_sd_error = None
        rh2_error = None

    phi_cnm, Z, sigma_d = guesses
    def odr_func(phi_cnm, h_sd):

        return k09.calc_rh2(h_sd, phi_cnm, Z=Z, sigma_d=sigma_d,
                            return_fractions=False)

    model = odr.Model(odr_func)
    data = odr.RealData(h_sd, rh2, sx=h_sd_error, sy=rh2_error)
    odr_instance = odr.ODR(data, model, beta0=[phi_cnm,])
    output = odr_instance.run()
    phi_cnm = output.beta[0]

    rh2_fit_params = (phi_cnm, Z, sigma_d)

    return rh2_fit_params

def fit_sternberg(h_sd, rh2, guesses=[10.0, 1.0, 10.0], rh2_error=None,
        h_sd_error=None, verbose=False, vary=[True, True, True],
        radiation_type='isotropic', bootstrap_residuals=False, nboot=100,
        odr_fit=True):

    '''
    Parameters
    ----------
    h_sd : array-like
        Hydrogen surface density in units of solar mass per parsec**2
    rh2 : array-like
        Ratio between molecular and atomic hydrogen masses.
    guesses : None, scalar, or M-length sequence.
        Initial guess for the parameters. See scipy.optimize.curve_fit.
    rh2_error : bool
        Error in rh2 parameter. Calculates a more accurate chi^2 statistic

    Returns
    -------
    rh2_fit_params : array-like, optional
        Model parameter fits.


    Residual bootstrapping:
    http://stats.stackexchange.com/questions/67519/bootstrapping-residuals-am-i-doing-it-right


    '''

    from scipy.optimize import curve_fit
    from scipy import stats
    from lmfit import minimize, Parameters, report_fit, Minimizer
    import lmfit
    from myscience import sternberg14 as s14
    import mystats
    import scipy.odr as odr

    #if not odr_fit:
    if 0:
        def calc_residual(params, h_sd, rh2):
            alphaG = params['alphaG'].value
            phi_g = params['phi_g'].value
            Z = params['Z'].value

            rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False)

            residual = rh2 - rh2_model

            return residual

        # Set parameter limits and initial guesses
        params = Parameters()
        params.add('alphaG',
                   value=guesses[0],
                   min=0.1,
                   max=500,
                   vary=vary[0])
        params.add('phi_g',
                   value=guesses[2],
                   min=0.01,
                   max=10,
                   vary=vary[2])
        params.add('Z',
                   value=guesses[1],
                   min=0.1,
                   max=4,
                   vary=vary[1])

        # Perform the fit!
        result = minimize(calc_residual,
                          params,
                          args=(h_sd, rh2),
                          #method='leastsq',
                          method='anneal',
                          )
        rh2_fit_params = (params['alphaG'].value, params['Z'].value,
                    params['phi_g'].value)
    else:

        if not odr_fit:
            h_sd_error = None
            rh2_error = None

        alphaG, Z, phi_g = guesses
        def odr_func(alphaG, h_sd):
            return s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False)

        #h_sd_error, rh2_error = None, None
        model = odr.Model(odr_func)
        data = odr.RealData(h_sd, rh2, sx=h_sd_error, sy=rh2_error)
        odr_instance = odr.ODR(data, model, beta0=[alphaG,])
        output = odr_instance.run()
        alphaG = output.beta[0]

        rh2_fit_params = (alphaG, Z, phi_g)

    return rh2_fit_params

def calc_coldens_products(nhi, av, dgr, nhi_error=0.0, av_error=0.0,
        dgr_error=0.0, ):

    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    # Calculate N(H2) and error
    # ---------------------------------------------------------------------------
    nh2 = calculate_nh2(nhi_image=nhi,
                              av_image=av,
                              dgr=dgr)

    # nh2 = (av / dgr - nhi) / 2
    # nh2_error = ((nh2 * ((av_error / av)**2 + (dgr_error / dgr)**2)**0.5)**2 \
    #              - nhi_error**2.0)**0.5 / 2
    comp_av_error = av_error / dgr
    comp_dgr_error = - av / dgr**2 * dgr_error
    comp_nhi_error = nhi_error
    nh2_error = 0.5 * (comp_av_error**2 + comp_dgr_error**2 + \
                      comp_nhi_error**2)**0.5

    # Convert to column density to surface density
    # ---------------------------------------------------------------------------
    hi_sd = calculate_sd(nhi,
                               sd_factor=1/1.25)
    hi_sd_error = calculate_sd(nhi_error,
                                     sd_factor=1/1.25)

    h2_sd = calculate_sd(nh2,
                               sd_factor=1/0.625)
    h2_sd_error = calculate_sd(nh2_error,
                               sd_factor=1/0.625)

    h_sd = hi_sd + h2_sd
    h_sd_error = (hi_sd_error**2 + h2_sd_error**2)**0.5

    # Write ratio between H2 and HI
    # ---------------------------------------------------------------------------
    rh2 = h2_sd / hi_sd
    comp_hi_error = - h2_sd / hi_sd**2 * hi_sd_error
    comp_h2_error = h2_sd_error / hi_sd
    rh2_error = (comp_hi_error**2 + comp_h2_error**2)**0.5

    return ((nhi, nh2, hi_sd, h2_sd, h_sd, rh2),
            (nhi_error, nh2_error, hi_sd_error, h2_sd_error, h_sd_error,
                rh2_error))



