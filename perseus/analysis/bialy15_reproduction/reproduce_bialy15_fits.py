#!/usr/bin/python

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def get_filename(region_name):

    DIR_DATA = '/d/bip3/ezbc/perseus/data/bialy15/'

    return DIR_DATA + 'Lee15_' + region_name + '_data.dat'

def plot_rh2_vs_h(data_dict, fits_dict):

    DIR_FIGURE = '/d/bip3/ezbc/perseus/figures/bialy15_reproduction/'

    for region_name in data_dict:
        data = data_dict[region_name]
        fits = fits_dict[region_name]

        filename = region_name + '_rh2_vs_h.png'

        plt.close(); plt.clf();
        plt.scatter(data['H'].values,
                    data['RH2'].values,
                    color='k',
                    alpha=0.5,
                    )
        plt.plot(fits['h_resampled'],
                 fits['rh2_resampled'],
                 color='r',
                 alpha=0.5,
                 linewidth=2,
                 )
        plt.yscale('log')
        plt.ylim([0.001, 100])
        plt.xlim([0, 80])
        text = region_name + '\n' + \
               r'$\alpha$G = {0:.1f}'.format(fits['alphaG']) + '\n' + \
               r'n = {0:.1f} [cm$^{{-3}}$]'.format(fits['n'])
        plt.annotate(text,
                     xy=(0.7,0.1),
                     xycoords='axes fraction',
                     )
        plt.savefig(DIR_FIGURE + filename, dpi=100)
        #plt.show()

def print_values(fits_dict):

    DIR_DATA = '/d/bip3/ezbc/perseus/data/bialy15/'
    f = open(DIR_DATA + 'fitting_results.txt', 'w')
    print('Region\tMy aG\tB+15 aG\tMy n\tB+15 n')
    f.write('Region\tMy aG\tB+15 aG\tMy n\tB+15 n\r')
    for region_name in fits_dict:
        fits = fits_dict[region_name]

        print(region_name + \
              '\t{0:.1f}\t{1:.1f}'.format(fits['alphaG'], fits['alphaG_B15']) + \
              '\t{0:.1f}\t{1:.1f}'.format(fits['n'], fits['n_B15']))

        f.write(region_name + \
                '\t{0:.1f}\t{1:.1f}'.format(fits['alphaG'],
                                              fits['alphaG_B15']) + \
                '\t{0:.1f}\t{1:.1f}\r'.format(fits['n'], fits['n_B15']))

    f.close()

def load_data():

    '''
    Column descriptions:
        column (1): HI surface density
        column (2): 1-sigma uncertainty in the HI surface density
        column (3): RH2 = H2 surface density / HI surface density
        column (4): 1-sigma uncertainty in the R_H2
        column (5): HI + H2 surface density
        column (6): 1-sigma uncertainty in the HI + H2 surface density

    '''

    # Define constants
    REGION_NAMES = ['B1', 'B1E', 'B5', 'IC348', 'NGC1333',]
    COLUMNS = ['HI', 'HI_ERROR', 'RH2', 'RH2_ERROR', 'H', 'H_ERROR']
    row_verify = np.array([6.506969e+00,
                           6.927944e-01,
                           -1.646940e-01,
                           -1.316955e-01,
                           5.435310e+00,
                           1.591130e+00,])
    #
    data_dict = {}
    for region_name in REGION_NAMES:

        # get the name of the file
        filename_data = get_filename(region_name)

        # read the file
        data_region = pd.read_csv(filename_data,
                                  delimiter=',',
                                  header=None,
                                  names=COLUMNS,
                                  engine='python',
                                  )
        if region_name == 'B1':
            assert np.sum(data_region.iloc[0].values - row_verify) < 1e6

        # assign
        data_dict[region_name] = data_region

    return data_dict

def fit_data(data_dict):

    # initialize new dict for fits
    fits_dict = {}
    guesses = {'B1': 26.1,
               'B1E': 23.8,
               'B5': 17.7,
               'IC348': 23.2,
               'NGC1333': 47.0}

    for region_name in data_dict:

        # Extract the data
        data = data_dict[region_name]
        h_sd = data['H'].values
        rh2 = data['RH2'].values
        h_sd_error = data['H_ERROR'].values
        rh2_error = data['RH2_ERROR'].values

        # fit a model to the data with weighted least squares
        alphaG, rh2_model, rh2_resampled, h_resampled = \
            fit_sternberg(h_sd,
                          rh2,
                          guesses=[guesses[region_name],],
                          rh2_error=rh2_error,
                          h_sd_error=h_sd_error,
                          phi_g=2.0,
                          Z=1.0,
                          radiation_type='isotropic',
                          )

        # write the fits to the fits dict
        fits_dict[region_name] = {}
        fits_dict[region_name]['alphaG'] = alphaG
        fits_dict[region_name]['alphaG_B15'] = guesses[region_name]
        fits_dict[region_name]['rh2_model'] = rh2_model
        fits_dict[region_name]['rh2_resampled'] = rh2_resampled
        fits_dict[region_name]['h_resampled'] = h_resampled

    return fits_dict

def fit_sternberg(h_sd, rh2, guesses=[10,], rh2_error=None,
        h_sd_error=None, radiation_type='beamed', phi_g=2.0, Z=1.0,):

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

    import sys
    sys.path.insert(0, '/usr/users/ezbc/.local/lib64/python2.7/site-packages/')
    from scipy.optimize import curve_fit, minimize
    from scipy import stats
    import scipy
    #from lmfit import minimize, Parameters, report_fit, Minimizer
    #import lmfit
    from myscience import sternberg14 as s14
    import mystats
    import scipy.odr as odr

    if any(np.isnan(h_sd)):
        raise ValueError('HSD should not have nans')

    mask = (rh2 < 0.0) | (np.isnan(h_sd))
    mask = np.isnan(h_sd)
    h_sd = h_sd[~mask]
    rh2 = rh2[~mask]
    rh2_error = rh2_error[~mask]

    if 0:
        def calc_residual(params, h_sd, rh2, rh2_error=None):
            alphaG = params['alphaG'].value

            rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False,
                                     radiation_type=radiation_type,
                                     )
            if rh2_error is not None:
                residual = (rh2 - rh2_model) / rh2_error
            else:
                residual = rh2 - rh2_model

            return residual

        # Set parameter limits and initial guesses
        params = Parameters()
        params.add('alphaG',
                   value=guesses[0],
                   min=0.01,
                   max=1000,
                   vary=True)

        # Perform the fit!
        result = minimize(calc_residual,
                          params,
                          args=(h_sd, rh2, rh2_error),
                          method='leastsq',
                          #method='anneal',
                          is_weighted=True,
                          )

        # best-fit parameter
        alphaG = params['alphaG'].value
    elif 0:
        def calc_model(h_sd, alphaG):
            rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False,
                                     radiation_type=radiation_type,
                                     )

            return rh2_model

        popt, pcov = curve_fit(calc_model,
                               h_sd,
                               rh2,
                               p0=guesses,
                               sigma=rh2_error,
                               method='lm',
                               )

        alphaG = popt[0]
    elif 1:
        def calc_model(alphaG):
            if 0:
                plt.scatter(h_sd, rh2)
                #plt.yscale('log')
                #plt.ylim([1e-3, 1e2])
                plt.annotate(alphaG, xy=(0.8, 0.1), xycoords='axes fraction')
                plt.show()

            rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False,
                                     radiation_type=radiation_type,
                                     )

            #mask = rh2_model < 0
            #print np.sum(rh2_model[~mask] < 0)
            mask = np.zeros(np.size(rh2), dtype=bool)

            chisq = np.nansum((rh2[~mask] - rh2_model[~mask])**2 / \
                    rh2_error[~mask]**2)

            return chisq

        res = minimize(calc_model,
                       guesses[0],
                       #method='powell',
                       #method='nelder-mead',
                       method='L-BFGS-B',
                       bounds=[[0, None],],
                       options={'disp':False},
                       )

        alphaG = res.x[0]
    else:
        def odr_func(alphaG, h_sd,):
            if 0:
                plt.scatter(h_sd, rh2)
                plt.yscale('log')
                plt.ylim([1e-3, 1e2])
                plt.annotate(alphaG, xy=(0.8, 0.1), xycoords='axes fraction')
                plt.show()
            rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False,
                                     radiation_type=radiation_type,
                                     )
            return rh2_model

        model = odr.Model(odr_func)
        data = odr.RealData(h_sd, rh2, sy=rh2_error)
        odr_instance = odr.ODR(data, model, beta0=guesses)
        output = odr_instance.run()
        alphaG = output.beta[0]

    # calculate model for observed H
    rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                             return_fractions=False,
                             radiation_type=radiation_type,
                             )

    # calculate model for resampled H
    h_sd_resampled = np.logspace(-3, 3, 1000)
    rh2_resampled = s14.calc_rh2(h_sd_resampled, alphaG, Z, phi_g=phi_g,
                             return_fractions=False,
                             radiation_type=radiation_type,
                             )

    return alphaG, rh2_model, rh2_resampled, h_sd_resampled

def calc_n(fits_dict):

    from myscience.sternberg14 import calc_n_H

    for region_name in fits_dict:
        fits = fits_dict[region_name]

        fits['n'] = calc_n_H(I_UV=1.0,
                             alphaG=fits['alphaG'],
                             phi_g=2.0,
                             Z_g=1.0,
                             calc_errors=False,
                             )
        fits['n_B15'] = calc_n_H(I_UV=1.0,
                             alphaG=fits['alphaG_B15'],
                             phi_g=2.0,
                             Z_g=1.0,
                             calc_errors=False,
                             )

    return fits_dict

def main():

    # load the data as a dict of pandas dataframes
    data = load_data()

    # fit the data with the sternberg model
    fits = fit_data(data)

    # calculate the volume densities
    fits = calc_n(fits)

    print_values(fits)

    # plot the data
    plot_rh2_vs_h(data, fits)

if __name__ == '__main__':
    main()

