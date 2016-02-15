#!/usr/bin/python


import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import myimage_analysis as myia
from multiprocessing.queues import Queue
import mygeometry as myg
import scipy
from local_module_plotting import *
from local_module_dataprep import *
from local_module_multiprocessing import *
from local_module_regions import *
from local_module_fitting import *
from local_module_bootstrapping import *



def fit_av_model(av, nhi, av_error=None, nhi_error=None, algebraic=False,
        nhi_background=None, plot_kwargs=None, init_guesses=[0.05, 0.05, 0],
        use_intercept=True, return_fit=False, fit_method='odr',
        odr_fit=True,):

    from lmfit import minimize, Parameters
    import lmfit
    import scipy.odr as odr

    if nhi_background is None:
        use_background = False
        init_guesses[1] = 0.0
        nhi_background = 0.0
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
    elif fit_method is not 'odr':
        # Set parameter limits and initial guesses
        params = Parameters()
        params.add('dgr_cloud',
                   value=init_guesses[0],
                   min=-0.5,
                   max=1,
                   )
        params.add('dgr_background',
                   value=init_guesses[1],
                   min=0.0,
                   max=1,
                   vary=use_background,
                   )
        params.add('intercept',
                   value=init_guesses[2],
                   min=-5,
                   max=5,
                   vary=use_intercept,
                   )

        #bin_edges = residuals_crop
        #counts = np.ones(residuals_crop.size - 1)

        def norm(params, av, nhi, av_error=None, nhi_background=None,):
            if nhi_background is None:
                nhi_background = 0.0
            if av_error is None:
                av_error = np.ones(av.shape)


            model = params['dgr_cloud'] * nhi + \
                    params['dgr_background'] * nhi_background + \
                    params['intercept']

            #if fit_method == 'leastsq':
                #norm = np.sum((av - model)**2 * (1.0/av_error**2)) / \
                #       np.sum(1.0/av_error**2)
            norm = np.sum((av - model)**2 * (1.0/av_error**2)) / \
                    np.sum(1.0/av_error**2)
            #norm = np.sum((av - model)**2)
            #else:
            #    norm = (av - model)

            return norm

        #print('fitting')
        # Perform the fit!
        result = minimize(norm,
                          params,
                          args=(av, nhi, av_error, nhi_background),
                          #method='leastsq',
                          #method=fit_method,
                          method='nelder',
                          )

        #print lmfit.report_fit(params)
        #print lmfit.printfuncs.report_ci(lmfit.conf_interval(result))
        #print lmfit.conf_interval(result)

        dgr_cloud = params['dgr_cloud'].value
        dgr_background = params['dgr_background'].value
        intercept = params['intercept'].value

        #if debugging:
        if 0:
            print('dgr = ', dgr_cloud)
            print('dgr background = ', dgr_background)
            print('intercept = ', intercept)
        if 0:
            plt.close(); plt.clf()
            background = dgr_background * nhi_background
            plt.errorbar(nhi, av - background,
                     yerr=(av_error),
                     linestyle='',
                     marker='o',
                     alpha=0.1,
                     markersize=2)
            xfit = np.linspace(0,50)
            plt.plot(xfit, dgr_cloud * xfit + intercept)
            plt.xlim(0, 22)
            plt.ylim(0, 15)
            plt.xlabel(r'N(HI)')
            plt.ylabel(r'$A_V$')
            plt.savefig(plot_kwargs['figure_dir'] + \
                        'diagnostics/av_nhi/' + plot_kwargs['filename_base']+ \
                        '_avnhi_bootsrap' + \
                        '{0:03d}.png'.format(plot_kwargs['bootstrap_num']))
    else:
        def odr_func(dgr, nhi):
            model = dgr * nhi
            return model

        model = odr.Model(odr_func)
        data = odr.RealData(nhi, av, sx=nhi_error, sy=av_error)
        odr_instance = odr.ODR(data, model, beta0=[0.1,])
        output = odr_instance.run()
        dgr_cloud = output.beta[0]

        dgr_background, intercept = 0.0, 0.0

        #print 'dgr =', dgr_cloud

    results = {'dgr_cloud': dgr_cloud,
              'dgr_background': dgr_background,
              'intercept': intercept}

    if return_fit:
        return (result, params), (dgr_cloud, dgr_background, intercept)
    else:
        return results

def fit_steady_state_models(h_sd, rh2, model_kwargs, rh2_error=None,
        h_sd_error=None, bootstrap_residuals=False, nboot=100, G0=1.0,
        odr_fit=False,):

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
                          rh2_error=rh2_error,
                          h_sd_error=h_sd_error,
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
                         G0=G0,
                         rh2_error=rh2_error,
                         h_sd_error=h_sd_error,
                         odr_fit=odr_fit,
                         )


        if bootstrap_residuals:
            phi_cnm, phi_cnm_error,Z_k09, Z_k09_error,phi_mol, phi_mol_error= \
                    result
            krumholz_results['phi_cnm_error'] = phi_cnm_error
            krumholz_results['Z_error'] = Z_k09_error
            krumholz_results['phi_mol_error'] = phi_mol_error
        else:
            phi_cnm, Z_k09, phi_mol = result
            phi_cnm_error, Z_k09_error, phi_mol_error = 3*[np.nan]
    else:
        alphaG, Z_s14, phi_g, phi_cnm, Z_k09, phi_mol = 6 * [np.nan]
        alphaG_error, Z_s14_error, phi_g_error, \
        phi_cnm_error, Z_k09_error, phi_mol_error = 6*[np.nan]

    # keep results
    sternberg_results['alphaG'] = alphaG
    sternberg_results['Z'] = Z_s14
    sternberg_results['phi_g'] = phi_g

    # keep results
    krumholz_results['phi_cnm'] = phi_cnm
    krumholz_results['Z'] = Z_k09
    krumholz_results['phi_mol'] = phi_mol


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

def add_hi_transition_calc(ss_model_result):

    h_sd_fit = np.linspace(0, 100, 1000)

    # To get HI transition, calculate model fits, then find where RH2 = 1
    for model_name in ss_model_result:
        model = ss_model_result[model_name]

        params = {}

        if 'sternberg' in model_name:
            params['phi_g'] = model['phi_g']
            params['Z'] = model['Z']
            params['alphaG'] = model['alphaG']
            model_fits = calc_sternberg(params,
                                      h_sd=h_sd_fit,
                                      return_fractions=False,
                                      return_hisd=False,
                                      )
        elif 'krumholz' in model_name:
            params['phi_cnm'] = model['phi_cnm']
            params['Z'] = model['Z']
            params['phi_mol'] = model['phi_mol']
            model_fits = calc_krumholz(params,
                                      h_sd=h_sd_fit,
                                      return_fractions=False,
                                      return_hisd=False,
                                      )

        rh2_fit = model_fits[0]

        try:
            if np.isnan(np.sum(rh2_fit)):
                hi_transition = np.nan
            else:
                # when R_H2 = 1, HI-to-H2 transition
                hi_transition = np.interp(1, rh2_fit, h_sd_fit) / 2.0

        except ValueError:
            hi_transition = np.nan

        model['hi_transition'] = hi_transition

def calc_krumholz(params, h_sd_extent=(0.001, 500), return_fractions=True,
        return_hisd=False, h_sd=None):

    '''
    Parameters
    ----------
    phi_cnm, Z : float
        Phi_cnm and Z parameters for Krumholz model.
    h_sd_extent : tuple
        Lower and upper bound of hydrogen surface densities with which to
        build the output model array.
    return_fractions : bool
        Return f_H2 and f_HI?

    Returns
    -------
    rh2_fit : array-like
        Model ratio between molecular and atomic hydrogen masses.
    h_sd_extended : list
        Model hydrogen surface density in units of solar mass per parsec**2.
    f_H2, f_HI : array-like, optional
        f_H2 = mass fraction of molecular hydrogen
        f_HI = mass fraction of atomic hydrogen

    '''

    from scipy import stats
    from myscience import krumholz09 as k09

    # Get large array of h_sd
    if h_sd is None:
        h_sd = np.linspace(h_sd_extent[0], h_sd_extent[1], 1e2)

    params = [params['phi_cnm'], params['Z'], params['phi_mol']]
    if params[0] <= 0 or np.isnan(params[0]):
        rh2_fit, f_H2, f_HI = np.empty(1), np.empty(1), np.empty(1)
    else:
        rh2_fit, f_H2, f_HI = k09.calc_rh2(h_sd,
                                           phi_cnm=params[0],
                                           Z=params[1],
                                           phi_mol=params[2],
                                           return_fractions=True)

    output = [rh2_fit, h_sd]

    if return_fractions:
        output.append(f_H2)
        output.append(f_HI)
    if return_hisd:
        #hi_sd = f_HI * h_sd
        hi_sd = (1 - f_H2) * h_sd
        output.append(hi_sd)

    return output

def fit_krumholz(h_sd, rh2, guesses=[10.0, 1.0, 10.0], h_sd_error=None,
        rh2_error=None, verbose=False, vary=[True, True, True],
        bootstrap_residuals=False, nboot=100, G0=1.0, odr_fit=True):

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
        def chisq(params, h_sd, rh2):
            phi_cnm = params['phi_cnm'].value
            phi_mol = params['phi_mol'].value
            Z = params['Z'].value

            rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol)

            chisq = np.sum(np.abs(rh2 - rh2_model))

            return chisq

        def calc_residual(params, h_sd, rh2, G0):
            phi_cnm = params['phi_cnm'].value
            phi_mol = params['phi_mol'].value
            Z = params['Z'].value

            rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol, G_0=G0)

            residual = rh2 - rh2_model

            return residual

        # Set parameter limits and initial guesses
        params = Parameters()
        params.add('phi_cnm',
                   value=guesses[0],
                   min=0.001,
                   max=1000,
                   vary=vary[0])
        params.add('phi_mol',
                   value=guesses[2],
                   min=1,
                   max=20,
                   vary=vary[2])
        params.add('Z',
                   value=guesses[1],
                   min=0.1,
                   max=4,
                   vary=vary[1])

        # Perform the fit!
        result = minimize(calc_residual,
                          params,
                          args=(h_sd, rh2, G0),
                          method='leastsq')

        if bootstrap_residuals:
            def resample_residuals(residuals):
                return np.random.choice(residuals,
                                        size=residuals.size,
                                        replace=True)
            phi_cnm = params['phi_cnm'].value
            phi_mol = params['phi_mol'].value
            Z = params['Z'].value
            rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol)
            residual = rh2 - rh2_model

            empty = np.empty(nboot)
            param_dict = {}
            param_names = ('phi_cnm', 'Z', 'phi_mol')
            for param_name in param_names:
                param_dict[param_name] = empty.copy()

            for i in xrange(nboot):
                rh2_resampled = rh2_model + resample_residuals(residual)
                result = minimize(calc_residual,
                                  params,
                                  args=(h_sd, rh2_resampled),
                                  #method='anneal',
                                  method='leastsq',
                                  )

                for param_name in param_names:
                    param_dict[param_name][i] = params[param_name].value

            rh2_fit_params = []
            for param_name in param_names:
                conf = mystats.calc_cdf_error(param_dict[param_name])
                rh2_fit_params.append(conf[0])
                rh2_fit_params.append(conf[1])
                #print param_name, conf
        else:
            rh2_fit_params = (params['phi_cnm'].value, params['Z'].value,
                    params['phi_mol'].value)
    elif odr_fit:

        phi_cnm, Z, phi_mol = guesses
        def odr_func(phi_cnm, h_sd):

            return k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol,
                                return_fractions=False)

        model = odr.Model(odr_func)
        data = odr.RealData(h_sd, rh2, sx=h_sd_error, sy=rh2_error)
        odr_instance = odr.ODR(data, model, beta0=[phi_cnm,])
        output = odr_instance.run()
        phi_cnm = output.beta[0]

        rh2_fit_params = (phi_cnm, Z, phi_mol)

    return rh2_fit_params

def analyze_krumholz_model(krumholz_results):

    ''' Calculates various properties of Krumholz model from fitted phi_cnm
    values.

    '''

    from myscience.krumholz09 import calc_T_cnm

    # Calculate T_cnm from Krumholz et al. (2009) Eq 19
    phi_cnm, phi_cnm_error, Z = \
        [krumholz_results[key] for key in ('phi_cnm', 'phi_cnm_error', 'Z')]
    T_cnm = calc_T_cnm(phi_cnm, Z=Z)
    T_cnm_error = []
    T_cnm_error.append(\
            T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[0], Z=Z))
    T_cnm_error.append(\
            T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[1], Z=Z))

    krumholz_results.update({'T_cnm': T_cnm,
                             'T_cnm_error': T_cnm_error})

    # Get fitted surface density ratios...
    params = [krumholz_results[param] for param in \
              krumholz_results['parameters']]

    rh2_fit, h_sd_fit, f_H2, f_HI = \
           calc_krumholz(params=params,
                          h_sd_extent=krumholz_results['h_sd_fit_range'],
                          return_fractions=True)

    krumholz_results.update({'rh2_fit': rh2_fit,
                             'h_sd_fit' : h_sd_fit,
                             'hi_sd_fit' : f_HI * h_sd_fit,
                             'f_H2' : f_H2,
                             'f_HI' : f_HI,})

    return krumholz_results

def calc_sternberg(params, h_sd_extent=(0.001, 500), return_fractions=True,
        return_hisd=False, h_sd=None):

    '''
    Parameters
    ----------
    alphaG, Z : float
        alphaG and Z parameters for sternberg model.
    h_sd_extent : tuple
        Lower and upper bound of hydrogen surface densities with which to
        build the output model array.
    return_fractions : bool
        Return f_H2 and f_HI?

    Returns
    -------
    rh2_fit : array-like
        Model ratio between molecular and atomic hydrogen masses.
    h_sd_extended : list
        Model hydrogen surface density in units of solar mass per parsec**2.
    f_H2, f_HI : array-like, optional
        f_H2 = mass fraction of molecular hydrogen
        f_HI = mass fraction of atomic hydrogen

    '''

    from scipy import stats
    from myscience import sternberg14 as s14

    # Create large array of h_sd
    if h_sd is None:
        h_sd = np.linspace(h_sd_extent[0], h_sd_extent[1], 1e2)

    params = [params['alphaG'], params['Z'], params['phi_g']]
    if params[0] <= 0 or np.isnan(params[0]):
        rh2_fit, f_H2, f_HI = np.empty(1), np.empty(1), np.empty(1)
    else:
        rh2_fit, f_H2, f_HI = s14.calc_rh2(h_sd,
                                           alphaG=params[0],
                                           Z=params[1],
                                           phi_g=params[2],
                                           return_fractions=True)

    output = [rh2_fit, h_sd]

    if return_fractions:
        output.append(f_H2)
        output.append(f_HI)
    if return_hisd:
        hi_sd = f_HI * h_sd
        output.append(hi_sd)

    return output

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

    if not odr_fit:
        def chisq(params, h_sd, rh2):
            alphaG = params['alphaG'].value
            phi_g = params['phi_g'].value
            Z = params['Z'].value

            rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False,
                                     radiation_type=radiation_type)

            chisq = np.sum(np.abs(rh2 - rh2_model))

            return chisq

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

        if bootstrap_residuals:
            def resample_residuals(residuals):
                return np.random.choice(residuals,
                                        size=residuals.size,
                                        replace=True)
            alphaG = params['alphaG'].value
            phi_g = params['phi_g'].value
            Z = params['Z'].value
            rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False)
            residual = rh2 - rh2_model

            empty = np.empty(nboot)
            param_dict = {}
            param_names = ('alphaG', 'Z', 'phi_g',)
            for param_name in param_names:
                param_dict[param_name] = empty.copy()

            for i in xrange(nboot):
                rh2_resampled = rh2_model + resample_residuals(residual)
                result = minimize(calc_residual,
                                  params,
                                  args=(h_sd, rh2_resampled),
                                  #method='leastsq',
                                  method='anneal',
                                  )

                for param_name in param_names:
                    param_dict[param_name][i] = params[param_name].value

            rh2_fit_params = []
            for param_name in param_names:
                conf = mystats.calc_cdf_error(param_dict[param_name])
                rh2_fit_params.append(conf[0])
                rh2_fit_params.append(conf[1])
                #print param_name, conf
        else:
            rh2_fit_params = (params['alphaG'].value, params['Z'].value,
                    params['phi_g'].value)
    elif odr_fit:

        alphaG, Z, phi_g = guesses
        def odr_func(alphaG, h_sd):

            return s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False)

        model = odr.Model(odr_func)
        data = odr.RealData(h_sd, rh2, sx=h_sd_error, sy=rh2_error)
        odr_instance = odr.ODR(data, model, beta0=[alphaG,])
        output = odr_instance.run()
        alphaG = output.beta[0]

        rh2_fit_params = (alphaG, Z, phi_g)

    return rh2_fit_params

def analyze_sternberg_model(sternberg_results):

    ''' Calculates various properties of Krumholz model from fitted phi_cnm
    values.

    '''

    # Get fitted surface density ratios...
    params = [sternberg_results[param] for param in \
              sternberg_results['parameters']]

    rh2_fit, h_sd_fit, f_H2, f_HI = \
           calc_sternberg(params=params,
                          h_sd_extent=sternberg_results['h_sd_fit_range'],
                          return_fractions=True)

    sternberg_results.update({'rh2_fit': rh2_fit,
                              'h_sd_fit' : h_sd_fit,
                              'hi_sd_fit' : (1 - f_H2) * h_sd_fit,
                              'f_H2' : f_H2,
                              'f_HI' : f_HI,})

    return sternberg_results


