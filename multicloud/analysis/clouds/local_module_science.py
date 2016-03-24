#!/usr/bin/python

def scale_phi_g(DGR, model_kwargs):

    phi_g = DGR / 0.053
    sigma_d = DGR / 0.053
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

