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

'''
Plotting
'''
import local_module_plotting as lm_plt

def print_av_error_stats(av, av_error):

    from scipy.stats import nanmedian

    error_above = nanmedian(av_error[av > 5])
    error_below = nanmedian(av_error[av <= 5])
    print('\n\tMedian Av error below 5 mag = {0:.2f} mag'.format(error_below))
    print('\n\tMedian Av error above 5 mag = {0:.2f} mag'.format(error_above))

def calc_model_plot_fit(analysis, model='krumholz', hsd=None,):

    if 'sternberg' in model:
        alphaG = analysis['alphaG']
        alphaG_low = alphaG - analysis['alphaG_error'][0]
        alphaG_high = alphaG + analysis['alphaG_error'][1]

        params = {'alphaG': alphaG,
                  'phi_g': analysis['phi_g'],
                  'Z': analysis['Z'],
                  }

        h_sd, hi_sd = calc_sternberg(params,
                                  h_sd_extent=(0, 100),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[1:3]

        params = {'alphaG': alphaG_low,
                  'phi_g': analysis['phi_g'],
                  'Z': analysis['Z'],
                  }
        hi_sd_low = calc_sternberg(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[2]

        params = {'alphaG': alphaG_high,
                  'phi_g': analysis['phi_g'],
                  'Z': analysis['Z'],
                  }
        hi_sd_high = calc_sternberg(params,
                                  h_sd_extent=(0, 100),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[2]
    elif 'krumholz' in model:
        phi_cnm = analysis['phi_cnm']
        phi_cnm_low = phi_cnm - analysis['phi_cnm_error'][0]
        phi_cnm_high = phi_cnm + analysis['phi_cnm_error'][1]

        params = {'phi_cnm': phi_cnm,
                  'phi_mol': analysis['phi_mol'],
                  'Z': analysis['Z'],
                  }
        h_sd, hi_sd = calc_krumholz(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[1:3]

        params = {'phi_cnm': phi_cnm_high,
                  'phi_mol': analysis['phi_mol'],
                  'Z': analysis['Z'],
                  }
        hi_sd_low = calc_krumholz(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[2]

        params = {'phi_cnm': phi_cnm_low,
                  'phi_mol': analysis['phi_mol'],
                  'Z': analysis['Z'],
                  }
        hi_sd_high = calc_krumholz(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[2]

    return h_sd, hi_sd, hi_sd_low, hi_sd_high

def plot_multicloud_results(results):

    print('\nPlotting multicloud results...')

    # Collect Data
    # =========================================================================

    spectra_list = []
    hi_range_kwargs_list = []
    av_list = []
    av_error_list = []
    nhi_list = []
    nh2_list = []
    rh2_list = []
    hsd_list = []
    hisd_list = []
    nhi_error_list = []
    nh2_error_list = []
    rh2_error_list = []
    hsd_error_list = []
    hisd_error_list = []
    hsd_median_error_list = []
    hisd_median_error_list = []
    hisd_cores_list = []
    hsd_cores_list = []
    rh2_cores_list = []
    hisd_error_cores_list = []
    hsd_error_cores_list = []
    hisd_median_error_cores_list = []
    hsd_median_error_cores_list = []
    rh2_error_cores_list = []
    model_results_list = []
    model_analysis_list = []
    model_analysis_dict = {}
    dgr_list = []
    dgr_error_list = []
    fit_params_list = []
    cloud_name_list = []
    core_names_list = []
    core_list = []
    cloud_model_fits_list = []
    stats_list = {'krumholz_results': {'sum_of_resid': [], 'BIC': []},
            'sternberg_results': {'sum_of_resid': [], 'BIC': []}}
    for i, cloud_name in enumerate(results):
        results_dict = results[cloud_name]
        figure_dir = results_dict['filenames']['figure_dir']
        results_dir = results_dict['filenames']['results_dir']
        plot_kwargs = results_dict['plot_kwargs']
        data_products = results_dict['data_products']
        spectra_list.append((data_products['hi_spectrum'],
                             data_products['hi_std_spectrum'],
                             data_products['co_spectrum'],)
                             )
        cloud_name_list.append(cloud_name)
        hi_range_kwargs_list.append(data_products['hi_range_kwargs'])
        #av_list.append(data_products['av_data_backsub'])
        av_list.append(data_products['av'])
        av_error_list.append(results_dict['data']['av_error_data'])
        nhi_list.append(data_products['nhi'])
        nh2_list.append(data_products['nh2'])
        rh2_list.append(data_products['rh2'])
        hsd_list.append(data_products['h_sd'])
        hisd_list.append(data_products['hi_sd'])
        nhi_error_list.append(data_products['nhi_error'])
        nh2_error_list.append(data_products['nh2_error'])
        rh2_error_list.append(data_products['rh2_error'])
        hsd_error_list.append(data_products['h_sd_error'])
        hisd_error_list.append(data_products['hi_sd_error'])
        hsd_median_error_list.append(data_products['h_sd_median_error'])
        hisd_median_error_list.append(data_products['hi_sd_median_error'])
        dgr_list.append(results_dict['mc_analysis']['dgr'])
        dgr_error_list.append(results_dict['mc_analysis']['dgr_error'])
        #fit_params_list.append(results_dict['params_summary'])
        model_results_list.append(results_dict['mc_results']['ss_model_results'])
        model_analysis_list.append(results_dict['mc_analysis'])
        model_analysis_dict[cloud_name] = results_dict['mc_analysis']

        # Modeling and cores
        global_args = results_dict['global_args']
        cores = global_args['ss_model_kwargs']['cores']
        cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']
        model_kwargs = global_args['ss_model_kwargs']['model_kwargs']
        rh2_core_list = []
        hsd_core_list = []
        hisd_core_list = []
        rh2_error_core_list = []
        hsd_error_core_list = []
        hsd_median_error_core_list = []
        hisd_error_core_list = []
        hisd_median_error_core_list = []
        core_names = []
        model_fits_list = []
        core_list.append(cores)
        for j, core in enumerate(cores_to_plot):
            core_indices = cores[core]['indices_orig']
            core_names.append(core)
            hisd_core_list.append(hisd_list[i][core_indices])
            hsd_core_list.append(hsd_list[i][core_indices])
            rh2_core_list.append(rh2_list[i][core_indices])
            hisd_error_core_list.append(hisd_error_list[i][core_indices])
            hsd_error_core_list.append(hsd_error_list[i][core_indices])
            hisd_median_error_core_list.append(hisd_median_error_list[i])
            hsd_median_error_core_list.append(hsd_median_error_list[i])
            rh2_error_core_list.append(rh2_error_list[i][core_indices])

            model_fits_list.append(refit_data(hsd_core_list[j],
                                              rh2_core_list[j],
                                              h_sd_error=hsd_error_core_list[j],
                                              rh2_error=rh2_error_core_list[j],
                                              model_kwargs=model_kwargs,
                                              )
                                   )

            # get the residual sum of squares
            for model in model_fits_list[j]:
                fits = refit_data(hsd_core_list[j],
                                  rh2_core_list[j],
                                  h_sd_error=hsd_error_core_list[j],
                                  rh2_error=rh2_error_core_list[j],
                                  h_sd_fit=hsd_core_list[j],
                                  model_kwargs=model_kwargs,
                                  )
                stats_list[model]['sum_of_resid'].append(\
                        fits[model]['sum_of_resid'])
                stats_list[model]['BIC'].append(\
                        fits[model]['BIC'])


            if 0:
                rh2_copy = rh2_list[i].copy()
                rh2_copy[core_indices] = 1000
                plt.imshow(rh2_copy, origin='lower')
                plt.savefig('/d/bip3/ezbc/scratch/core_' + core + '.png')

        cloud_model_fits_list.append(model_fits_list)
        core_names_list.append(core_names)
        hisd_cores_list.append(hisd_core_list)
        hsd_cores_list.append(hsd_core_list)
        rh2_cores_list.append(rh2_core_list)
        hisd_error_cores_list.append(hisd_error_core_list)
        hsd_error_cores_list.append(hsd_error_core_list)
        hisd_median_error_cores_list.append(hisd_median_error_core_list)
        hsd_median_error_cores_list.append(hsd_median_error_core_list)

    # calculate statistics on hi
    hi_dict = calc_hi_statistics(cloud_name_list, core_names_list,
                                 hisd_cores_list, hsd_cores_list,
                                 rh2_cores_list, model_analysis_dict,
                                 filename=global_args['filename_hi_props'],
                                 )

    # plot CDFs of diffuse fraction LOS
    lm_plt.plot_diffusefraction_cdfs(hi_dict)

    # Print results
    # =========================================================================
    # Write results to a
    #print_av_error_stats(av_list[0], av_error_list[0])

    print_BIC_results(stats_list, core_names_list)

    filename = results_dir + 'tables/multicloud_model_params.tex'
    write_model_params_table(model_analysis_dict,
                             filename,
                             models=('krumholz','sternberg'))

    filename = results_dir + 'tables/multicloud_hi_core_properties.tex'
    write_core_HI_table(hi_dict,
                             filename,
                             )

    # Write param summary to dataframe for ease of use
    filename = results_dir + 'tables/multicloud_model_params.csv'
    write_param_csv(model_analysis_dict,
                    core_list,
                    cloud_name_list,
                    filename,
                    hi_dict=hi_dict,
                    )

    # Write nhi properties
    filename = results_dir + 'tables/nhi_properties.csv'
    write_nhi_properties_csv(nhi_list,
                             cloud_name_list,
                             filename,
                             )

    # Write param summary to dataframe for ease of use
    filename = results_dir + 'tables/multicloud_model_summary.pickle'
    write_fit_summary_dict(model_analysis_dict,
                           core_list,
                           cloud_name_list,
                           filename,
                           hi_dict=hi_dict,
                           )

    # Write table for
    filename = results_dir + 'tables/multicloud_hi_transitions.csv'
    hi_trans_dict = collect_hi_transition_results(model_analysis_list,
                                                  cloud_name_list,
                                                  filename=filename,)

    # Plot the results
    # =========================================================================

    # Plot HI spectra
    # -------------------------------------------------------------------------
    hi_vel_axis = results_dict['data']['hi_vel_axis']
    co_vel_axis = results_dict['data']['co_vel_axis']


    filetypes = ['png', 'pdf']
    for filetype in filetypes:
        filename = plot_kwargs['figure_dir'] + \
                   'spectra/multicloud_spectra.' + filetype

        lm_plt.plot_spectra_grid(spectra_list,
                     hi_range_kwargs_list=hi_range_kwargs_list,
                     names_list=cloud_name_list,
                     hi_vel_axis=hi_vel_axis,
                     co_vel_axis=co_vel_axis,
                     filename=filename,
                     limits=[-30, 30, -10, 59],
                     )

        filename = results_dir + 'tables/multicloud_vel_ranges.tex'
        write_hi_vel_range_table(cloud_name_list,
                                 hi_range_kwargs_list,
                                 filename,
                                 dgr_list=dgr_list,
                                 dgr_error_list=dgr_error_list)

        # Plot N(HI) vs. Av
        # ---------------------------------------------------------------------
        filename = plot_kwargs['figure_dir'] + \
                   'av_nhi/multicloud_av_vs_nhi.' + filetype
        lm_plt.plot_av_vs_nhi_grid(av_list,
                            nhi_list,
                            av_error_list=av_error_list,
                            fit_params_list=model_analysis_list,
                            names_list=cloud_name_list,
                            filename=filename,
                            levels=(0.99, 0.98, 0.95, 0.86, 0.59),
                            poly_fit=False,
                            limits=[-2, 19, 0, 20]
                            )

        # Plot Av PDF
        # ---------------------------------------------------------------------
        filename = plot_kwargs['figure_dir'] + \
                   'pdfs/multicloud_pdfs.' + filetype
        lm_plt.plot_pdf_grid(av_list,
                  nhi_list=nhi_list,
                  nh2_list=nh2_list,
                  dgr_list=dgr_list,
                  hi_trans_dict=hi_trans_dict,
                  #limits = [-4,3,1,10000],
                  names_list=cloud_name_list,
                  #limits = [0.07,14,7,6000],
                  limits=[10**-2, 3 * 10**1, 5 * 10**-4, 3 * 10**-1],
                  scale=(1,1),
                  filename=filename,
                  #core_names=core_name_list,
                  show=False)

        # Plot HI vs H
        # ---------------------------------------------------------------------
        levels = (0.9, 0.8, 0.6, 0.3)
        ncols = 2
        for i, cloud in enumerate(cloud_name_list):
            core_names = core_names_list[i]
            print('\n\tPlotting Models')
            # RH2 vs. H SD for L1478
            # -----------------------------------------------------------------
            if cloud == 'taurus':
                single_core = 'G169.32-16.17'
                filename = plot_kwargs['figure_dir'] + \
                           'models/' + cloud + '_rh2_vs_hsd_' + \
                           single_core + '.' + filetype
                index = core_names.index(single_core)
                lm_plt.plot_rh2_vs_h((hsd_cores_list[i][index],),
                              (hisd_cores_list[i][index],),
                              core_names=(core_names[index],),
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              ylimits=[10**-3, 10**2],
                              levels=levels,
                              #scale=('log', 'linear'),
                              scale=('linear', 'log'),
                              filename=filename,
                              ncols=ncols
                              )

                filename = plot_kwargs['figure_dir'] + \
                           'models/' + cloud + '_rh2_vs_hsd_' + \
                           single_core + '.ps'
                lm_plt.plot_rh2_vs_h((hsd_cores_list[i][index],),
                              (hisd_cores_list[i][index],),
                              core_names=(core_names[index],),
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              ylimits=[10**-3, 10**2],
                              levels=levels,
                              #scale=('log', 'linear'),
                              scale=('linear', 'log'),
                              filename=filename,
                              ncols=ncols
                              )

            # HI SD vs. H SD
            # -----------------------------------------------------------------
            filename = plot_kwargs['figure_dir'] + \
                       'models/' + cloud + '_hisd_vs_hsd.' + filetype
            lm_plt.plot_hi_vs_h_grid(hsd_cores_list[i],
                              hisd_cores_list[i],
                              core_names=core_names,
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              model_fits=cloud_model_fits_list[i],
                              hsd_error_list=hsd_median_error_cores_list[i],
                              hisd_error_list=hisd_median_error_cores_list[i],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              levels=levels,
                              #scale=('log', 'linear'),
                              #scale=('linear', 'log'),
                              scale=('linear', 'linear'),
                              filename=filename,
                              ncols=ncols
                              )
            # HI SD vs. H2 SD
            # -----------------------------------------------------------------
            filename = plot_kwargs['figure_dir'] + \
                       'models/' + cloud + '_hisd_vs_h2sd.' + filetype
            lm_plt.plot_hi_vs_h2_grid(hsd_cores_list[i],
                              hisd_cores_list[i],
                              core_names=core_names,
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              model_fits=cloud_model_fits_list[i],
                              hsd_error_list=hsd_median_error_cores_list[i],
                              hisd_error_list=hisd_median_error_cores_list[i],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              levels=levels,
                              #scale=('log', 'linear'),
                              #scale=('linear', 'log'),
                              scale=('linear', 'linear'),
                              filename=filename,
                              ncols=ncols
                              )
            # H2 SD vs. H SD
            # -----------------------------------------------------------------
            filename = plot_kwargs['figure_dir'] + \
                       'models/' + cloud + '_h2sd_vs_hsd.' + filetype
            lm_plt.plot_h2_vs_h_grid(hsd_cores_list[i],
                              hisd_cores_list[i],
                              core_names=core_names,
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              model_fits=cloud_model_fits_list[i],
                              hsd_error_list=hsd_median_error_cores_list[i],
                              hisd_error_list=hisd_median_error_cores_list[i],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              levels=levels,
                              #scale=('log', 'linear'),
                              #scale=('linear', 'log'),
                              scale=('linear', 'linear'),
                              filename=filename,
                              ncols=ncols
                              )
            # HI CDF
            # -----------------------------------------------------------------
            filename = plot_kwargs['figure_dir'] + \
                       'models/' + cloud + '_hisd_cdf.' + filetype
            lm_plt.plot_hi_cdf_grid(hsd_cores_list[i],
                              hisd_cores_list[i],
                              core_names=core_names,
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              #xlimits=[-1, 20],
                              levels=levels,
                              #scale=('log', 'linear'),
                              #scale=('linear', 'log'),
                              scale=('linear', 'linear'),
                              filename=filename,
                              ncols=ncols
                              )

            # RH2 vs. H SD
            # -----------------------------------------------------------------
            filename = plot_kwargs['figure_dir'] + \
                       'models/' + cloud + '_rh2_vs_hsd.' + filetype
            lm_plt.plot_rh2_vs_h_grid(hsd_cores_list[i],
                              hisd_cores_list[i],
                              core_names=core_names,
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              model_fits=cloud_model_fits_list[i],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              ylimits=[10**-3, 10**2],
                              levels=levels,
                              #scale=('log', 'linear'),
                              scale=('linear', 'log'),
                              filename=filename,
                              ncols=ncols
                              )

        if 0:
            filename = plot_kwargs['figure_dir'] + \
                       'av_nhi/multicloud_av_vs_nhi_log.' + filetype
            lm_plt.plot_av_vs_nhi_grid(av_list,
                                nhi_list,
                                av_error_list=av_error_list,
                                fit_params_list=fit_params_list,
                                names_list=cloud_name_list,
                                filename=filename,
                                levels=(0.99, 0.98, 0.95, 0.86, 0.59),
                                poly_fit=True,
                                scale=['log','log'],
                                #limits=[2, 20, -13, 20]
                                limits=[0.1, 100, 0.01, 100]
                                )

            filename = plot_kwargs['figure_dir'] + \
                       'av_nhi/multicloud_av_vs_nhi_scatter.' + filetype
            lm_plt.plot_av_vs_nhi_grid(av_list,
                                nhi_list,
                                av_error_list=av_error_list,
                                fit_params_list=fit_params_list,
                                names_list=cloud_name_list,
                                filename=filename,
                                levels=(0.99, 0.98, 0.95, 0.86, 0.59),
                                poly_fit=True,
                                contour=False,
                                limits=[2, 20, -6, 10]
                                )

'''
Data Prep functions
'''

def print_BIC_results(stats_list, core_names):

    bayes_factor = np.array(stats_list['krumholz_results']['BIC']) - \
                   np.array(stats_list['sternberg_results']['BIC'])

    # cores gathers by cloud, unravel
    core_names = np.ravel(core_names)

    print('Median difference in BIC between k09 and s14 models:')
    print(np.median(bayes_factor))

    print('Bayes Factor for each core:')
    print(bayes_factor)

    import matplotlib.pyplot as plt
    import myplotting as myplt
    import mystats

    print('Confidence intervals on Bayes Factor')
    print(mystats.calc_cdf_error(bayes_factor, alpha=0.5))
    print(np.sort(bayes_factor))

    print('Core with max BF of {0:.0f}'.format(np.max(bayes_factor)))
    print(core_names[np.argmax(bayes_factor)])
    print('Core with min BF of {0:.0f}'.format(np.min(bayes_factor)))
    print(core_names[np.argmin(bayes_factor)])

    plt.close(); plt.clf()
    myplt.plot_cdf(bayes_factor)
    plt.xlabel('K+09 - S+14')
    plt.title('Bayes Factor CDF')
    plt.savefig('/d/bip3/ezbc/multicloud/figures/models/bayes_factor_cdf.png')


def calc_hi_statistics(cloud_name_list, core_names_list,
                                 hisd_cores_list, h_sd_cores_list,
                                 rh2_cores_list, model_analysis_dict,
                                 filename=None):

    hi_dict = {}
    for i, cloud in enumerate(cloud_name_list):

        # initialize hi properties
        hi_dict[cloud] = {}
        hi_dict[cloud]['cores'] = []
        hi_dict[cloud]['hi_sd_mean'] = []
        hi_dict[cloud]['hi_sd_median'] = []
        hi_dict[cloud]['hi_sd_std'] = []
        hi_dict[cloud]['fraction_LOS_diffuse'] = []
        hi_dict[cloud]['fraction_LOS_sf'] = []
        hi_dict[cloud]['sf_threshold'] = []

        # add a row for each core
        for j, core in enumerate(core_names_list[i]):
            hi = hisd_cores_list[i][j]
            h = h_sd_cores_list[i][j]
            rh2 = rh2_cores_list[i][j]
            hi_dict[cloud]['cores'].append(core_names_list[i][j])
            hi_dict[cloud]['hi_sd_mean'].append(np.nanmean(hi))
            hi_dict[cloud]['hi_sd_median'].append(scipy.stats.nanmedian(hi))
            hi_dict[cloud]['hi_sd_std'].append(np.nanstd(hi))

            # frac of diffuse LOS
            indices_diffuse = np.where(rh2 < 1.0)[0]
            frac_diffuse = np.size(indices_diffuse) / float(np.size(rh2))
            hi_dict[cloud]['fraction_LOS_diffuse'].append(frac_diffuse)

            #frac above sf threshold
            model = model_analysis_dict[cloud]['cores'][core]
            sf_threshold = model['sternberg_results']['sf_threshold']
            hi_dict[cloud]['sf_threshold'].append(sf_threshold)
            indices_sf = np.where(h > sf_threshold)[0]
            frac_sf = np.size(indices_sf) / float(np.size(h))
            hi_dict[cloud]['fraction_LOS_sf'].append(frac_sf)

    # save the dict?
    if filename is not None:
        pickle.dump(hi_dict, open(filename, 'w'),)

    return hi_dict

def collect_hi_transition_results(model_analysis_list, cloud_list,
        filename=None):

    import pandas as pd

    hi_trans_dict = {}

    for i, cloud in enumerate(cloud_list):
        for core_name in model_analysis_list[i]['cores']:
            core = model_analysis_list[i]['cores'][core_name]
            hi_trans_k09 = \
                core['krumholz_results']['hi_transition']
            #print hi_trans_k09
            hi_trans_k09_error = \
                core['krumholz_results']['hi_transition_error']
            hi_trans_s14 = \
                core['sternberg_results']['hi_transition']
            hi_trans_s14_error = \
                core['sternberg_results']['hi_transition_error']

            hi_trans_dict[core_name] = \
                {'cloud': cloud,
                 'k09_transition': hi_trans_k09,
                 'k09_transition_error': hi_trans_k09_error,
                 's14_transition': hi_trans_s14,
                 's14_transition_error': hi_trans_s14_error,
                 }

    return hi_trans_dict

def print_dict_keys(d):

    for key in d:
        print(key)
        if type(d[key]) is dict:
            print('--')
            print_dict_keys(d[key])

def write_fit_summary_dict(mc_analysis_dict, core_list, cloud_name_list,
        filename, hi_dict=None,):

    import pandas as pd
    import pickle

    d = {}

    # Collect parameter names for each model for each core
    for i, cloud in enumerate(cloud_name_list):
        mc_analysis = mc_analysis_dict[cloud]
        core_names = np.sort(mc_analysis['cores'].keys())

        cores = core_list[i]

        # each core correspond to a row
        for cloud_row, core_name in enumerate(core_names):
            core = mc_analysis['cores'][core_name]

            core_props = cores[core_name]

            d[core_name] = {}
            core_new = d[core_name]
            core_new['cloud'] = cloud
            ra_deg = core_props['center_wcs'][0]
            dec_deg = core_props['center_wcs'][1]
            #ra_deg = 15*(ra[0] + ra[1] / 60. + ra[2] / 3600.)
            #dec_deg = dec[0] + dec[1] / 60. + dec[2] / 3600.
            core_new['ra'] = ra_deg
            core_new['dec'] = dec_deg
            try:
                core_new['temp'] = core_props['temp']
                core_new['temp_error'] = core_props['temp_error']
            except KeyError:
                core_new['temp'], core_new['temp_error'] = 17, 1
            core_new['region_vertices'] = core_props['poly_verts']['wcs']

            if hi_dict is not None:
                hi_cloud_dict = hi_dict[cloud]
                index = np.where(hi_cloud_dict['cores'] == core)[0]
                keys = hi_cloud_dict.keys()

                # remove cores from dict
                keys.pop(keys.index('cores'))

                hi_core_dict = {}
                for param in hi_cloud_dict:
                    hi_core_dict[param] = hi_cloud_dict[param]

                core_new['hi_props'] = hi_core_dict

            # append model params and errors to row
            for model in ('krumholz', 'sternberg'):
                if model == 'krumholz':
                    params_to_write = ['phi_cnm', 'Z', 'phi_mol',
                                       'hi_transition']
                else:
                    params_to_write = ['alphaG', 'Z', 'phi_g',
                                       'hi_transition']
                core_new[model] = {}
                for j, param_name in enumerate(params_to_write):
                    param = \
                        core[model + '_results'][param_name]
                    param_error = \
                        core[model + '_results'][param_name + '_error']

                    core_new[model][param_name] = param
                    #core_new[model][param_name + '_error_low'] = \
                    #        param_error[0]
                    #core_new[model][param_name + '_error_high'] = \
                    #        param_error[1]
                    core_new[model][param_name + '_error'] = param_error

    with open(filename, 'wb') as f:
        pickle.dump(d, f)

def write_param_csv(mc_analysis_dict, core_list, cloud_name_list, filename,
        hi_dict=None,):

    import pandas as pd

    d = {}

    #print_dict_keys(mc_analysis_dict)
    params_to_write = ['phi_cnm', 'alphaG', 'hi_transition']

    d['cloud'] = []
    d['core'] = []
    d['ra'] = []
    d['dec'] = []
    d['region_vertices'] = []
    if hi_dict is not None:
        d['fraction_diffuse'] = []
    for param in params_to_write:
        d[param] = []
        d[param + '_error_low'] = []
        d[param + '_error_high'] = []

    # Collect parameter names for each model for each core
    for i, cloud in enumerate(cloud_name_list):
        mc_analysis = mc_analysis_dict[cloud]
        core_names = np.sort(mc_analysis['cores'].keys())

        cores = core_list[i]

        # each core correspond to a row
        for cloud_row, core_name in enumerate(core_names):
            core = mc_analysis['cores'][core_name]

            core_props = cores[core_name]

            d['cloud'].append(cloud)
            d['core'].append(core_name)
            ra_deg = core_props['center_wcs'][0]
            dec_deg = core_props['center_wcs'][1]
            #ra_deg = 15*(ra[0] + ra[1] / 60. + ra[2] / 3600.)
            #dec_deg = dec[0] + dec[1] / 60. + dec[2] / 3600.
            d['ra'].append(ra_deg)
            d['dec'].append(dec_deg)
            d['region_vertices'].append(core_props['poly_verts']['wcs'])

            if hi_dict is not None:
                hi_cloud_dict = hi_dict[cloud]
                index = hi_cloud_dict['cores'].index(core_name)
                d['fraction_diffuse'].append(\
                    hi_cloud_dict['fraction_LOS_diffuse'][index])

            # append model params and errors to row
            for model in ('krumholz', 'sternberg'):
                if model == 'krumholz':
                    params_to_write = ['phi_cnm', 'hi_transition']
                else:
                    params_to_write = ['alphaG',]
                for i, param_name in enumerate(params_to_write):
                    param = \
                        core[model + '_results'][param_name]
                    param_error = \
                        core[model + '_results'][param_name + '_error']

                    d[param_name].append(param)
                    d[param_name + '_error_low'].append(param_error[0])
                    d[param_name + '_error_high'].append(param_error[1])

            #print d[param_name]

    # Create dataframe and write it!
    df = pd.DataFrame(data=d,)
    df.to_csv(filename,
              sep=',',
              columns=('cloud',
                       'core',
                       'ra',
                       'dec',
                       'phi_cnm',
                       'phi_cnm_error_low',
                       'phi_cnm_error_high',
                       'alphaG',
                       'alphaG_error_low',
                       'alphaG_error_high',
                       'hi_trans',
                       'hi_trans_error_low',
                       'hi_trans_error_high',
                       ),
              index=False,
              )

    df.save(filename.replace('csv', 'pickle'))

def write_hi_vel_range_table(names_list, hi_range_kwargs_list, filename,
        dgr_list=None, dgr_error_list=None,):

    # Open file to be appended
    f = open(filename, 'wb')

    if 0:
        for i in xrange(0, len(names_list)):
            cloud_name = names_list[i]
            hi_range_kwargs = hi_range_kwargs_list[i]
            vel_range = hi_range_kwargs['vel_range']
            vel_range_error = hi_range_kwargs['hi_range_error']

            row_text = cloud_name.capitalize()

            row_text = add_row_element(row_text,
                            vel_range,
                            text_format='[{0:.0f}, {1:.0f}]')
            row_text = add_row_element(row_text,
                            vel_range_error,
                            text_format='{0:.0f}')
            row_text += ' \\\\[0.1cm] \n'

            f.write(row_text)
    else:
        row_text = ''

        for cloud in names_list:
            row_text = add_row_element(row_text, cloud.capitalize())
        row_text += ' \\\\[0.1cm] \n'
        f.write(row_text)

        if dgr_list is not None:
            row_text = r'DGR'
            for i in xrange(len(names_list)):
                dgr = dgr_list[i]
                dgr_error = dgr_error_list[i]

                row_text = \
                    add_row_element(row_text,
                              (dgr*100, dgr_error[1]*100, dgr_error[0]*100),
                              text_format=r'{0:.1f}$^{{+{1:.1f}}}_{{-{2:.1f}}}$')
            row_text += ' \\\\[0.1cm] \n'
            f.write(row_text)

            # write next row with dgr units
            row_text = r'[10$^{-22}$ cm$^{2}$ mag] & & & '
            row_text += ' \\\\[0.1cm] \n'
            f.write(row_text)

        # write hi vel range row
        row_text = r'\hi\ Range'
        for i in xrange(len(names_list)):
            hi_range_kwargs = hi_range_kwargs_list[i]
            vel_range = hi_range_kwargs['vel_range']
            vel_range_error = hi_range_kwargs['hi_range_error']

            row_text = add_row_element(row_text,
                            (vel_range[0], vel_range[1], vel_range_error),
                            text_format=r'[{0:.0f}, {1:.0f}]\,$\pm$\,{2:.0f}')
        row_text += ' \\\\[0.1cm] \n'
        f.write(row_text)

        # write second row for units
        row_text = r'[\kms, \kms] & & & '
        row_text += ' \\\\[0.1cm] \n'
        f.write(row_text)

    f.close()

def write_nhi_properties_csv(nhi_list, cloud_name_list, filename):

    import pandas as pd

    d = {}
    d['cloud'] = []
    d['nhi_std'] = []
    d['nhi_median'] = []
    d['nhi_mean'] = []
    d['nhi_min'] = []
    d['nhi_max'] = []

    # Collect parameter names for each model for each core
    for i, cloud in enumerate(cloud_name_list):
        d['cloud'].append(cloud)

        nhi = nhi_list[i]
        nhi = nhi[~np.isnan(nhi)]
        d['nhi_std'].append(np.std(nhi))
        d['nhi_median'].append(np.median(nhi))
        d['nhi_mean'].append(np.mean(nhi))
        d['nhi_min'].append(np.min(nhi))
        d['nhi_max'].append(np.max(nhi))

    # Create dataframe and write it!
    df = pd.DataFrame(data=d,)
    df.to_csv(filename,
              sep='\t',
              columns=('cloud',
                       'nhi_std',
                       'nhi_median',
                       'nhi_mean',
                       'nhi_min',
                       'nhi_max',
                       ),
              float_format='%.1f',
              index=False,
              )

    df.save(filename.replace('csv', 'tsv'))
    df.save(filename.replace('csv', 'pickle'))

def write_model_params_table(mc_analysis_dict, filename, models=('krumholz',)):

    # Open file to be appended
    f = open(filename, 'wb')

    text_param_format ='{0:.1f}$^{{+{1:.1f}}}_{{-{2:.1f}}}$'

    #print_dict_keys(mc_analysis_dict)
    params_to_write = ['phi_cnm', 'alphaG']

    # Collect parameter names for each model for each core
    for cloud in ('california', 'perseus', 'taurus'):
        mc_analysis = mc_analysis_dict[cloud]
        core_names = np.sort(mc_analysis['cores'].keys())
        # each core correspond to a row
        for cloud_row, core_name in enumerate(core_names):
            core = mc_analysis['cores'][core_name]
            if cloud_row == 0:
                row_text = cloud.capitalize()
            else:
                row_text = ''
            row_text = add_row_element(row_text,
                                       core_name)

            # append model params and errors to row
            for model in ('sternberg', 'krumholz',):
                if model == 'krumholz':
                    params_to_write = ['phi_cnm', 'hi_transition']
                else:
                    params_to_write = ['alphaG',]
                for i, param_name in enumerate(params_to_write):
                    param = \
                        core[model + '_results'][param_name]
                    param_error = \
                        core[model + '_results'][param_name + '_error']

                    param_info = (param, param_error[1], param_error[0])

                    #if param_name == 'alphaG':
                        #print core_name, param_info

                    row_text = \
                        add_row_element(row_text,
                                        param_info,
                                        text_format=text_param_format)


            row_text += ' \\\\[0.1cm] \n'
            if cloud_row == len(mc_analysis['cores']) - 1 \
                and cloud != 'taurus':
                row_text += '\hline  \\\\[-0.2cm] \n'
            elif cloud_row == len(mc_analysis['cores']) - 1 and \
                    cloud == 'taurus':
                row_text.replace(r'\\[0.1cm] \n', '')

            f.write(row_text)

    f.close()

def write_core_HI_table(hi_dict, filename,):

    # Open file to be appended
    f = open(filename, 'wb')

    text_param_format_sd ='{0:.1f}'
    text_param_format_frac ='{0:.0f}'

    params_to_write = ['hi_sd_mean', 'hi_sd_median', 'hi_sd_std',
    'fraction_LOS_diffuse', 'fraction_LOS_sf', 'sf_threshold']

    # Collect parameter names for each model for each core
    for cloud in ('california', 'perseus', 'taurus'):
        core_dict = hi_dict[cloud]
        core_names = np.sort(core_dict['cores'])
        # each core correspond to a row
        for cloud_row, core_name in enumerate(core_names):

            if cloud_row == 0:
                row_text = cloud.capitalize()
            else:
                row_text = ''

            row_text = add_row_element(row_text,
                                       core_name)

            for i, param_name in enumerate(params_to_write):
                param = \
                    core_dict[param_name][cloud_row]

                param_info = param

                if 'fraction' in param_name:
                    text_param_format = text_param_format_frac
                    #print param_info
                    param_info = param_info * 100.0
                else:
                    text_param_format = text_param_format_sd


                #if param_name == 'alphaG':
                    #print core_name, param_info

                row_text = \
                    add_row_element(row_text,
                                    param_info,
                                    text_format=text_param_format)


            row_text += ' \\\\[0.1cm] \n'
            if cloud_row == len(core_names) - 1 \
                and cloud != 'taurus':
                row_text += '\hline  \\\\[-0.2cm] \n'
            elif cloud_row == len(core_names) - 1 and \
                    cloud == 'taurus':
                row_text.replace(r'\\[0.1cm] \n', '')

            f.write(row_text)

    f.close()

def add_row_element(row_text, element, text_format='{0:s}'):

    if type(element) is list or type(element) is tuple:
        return row_text + ' & ' + text_format.format(*element)
    else:
        return row_text + ' & ' + text_format.format(element)

def write_coldens_maps(results_dict, global_args):

    from astropy.io import fits

    FILENAME_HI = '/d/bip3/ezbc/multicloud/data/nhi/' + \
                  global_args['cloud_name'] + '_hisurfdens.fits'
    FILENAME_H2 = '/d/bip3/ezbc/multicloud/data/nh2/' + \
                  global_args['cloud_name'] + '_h2surfdens.fits'
    FILENAME_HI_ERROR = '/d/bip3/ezbc/multicloud/data/nhi/' + \
                  global_args['cloud_name'] + '_hisurfdens_error.fits'
    FILENAME_H2_ERROR = '/d/bip3/ezbc/multicloud/data/nh2/' + \
                  global_args['cloud_name'] + '_h2surfdens_error.fits'

    # get av_header
    header = results_dict['data']['av_header'].copy()
    header['BUNIT'] = 'Msun.pc2'

    hi_sd_image = results_dict['data_products']['hi_sd']
    h2_sd_image = results_dict['data_products']['h2_sd']
    hi_sd_error_image = results_dict['data_products']['hi_sd_error']
    h2_sd_error_image = results_dict['data_products']['h2_sd_error']

    print('median h2 error:', scipy.stats.nanmedian(h2_sd_error_image.ravel()))
    print('median hi error:', scipy.stats.nanmedian(hi_sd_error_image.ravel()))

    fits.writeto(FILENAME_HI, hi_sd_image, header=header, clobber=True)
    fits.writeto(FILENAME_H2, h2_sd_image, header=header, clobber=True)
    fits.writeto(FILENAME_HI_ERROR, hi_sd_error_image, header=header, clobber=True)
    fits.writeto(FILENAME_H2_ERROR, h2_sd_error_image, header=header, clobber=True)

'''
Multiprocessing functions
'''
class QueueGet(Queue):
    """Queue which will retry if interrupted with EINTR."""
    def get( block=True, timeout=None):
        return retry_on_eintr(Queue.get,  block, timeout)

def retry_on_eintr(function, *args, **kw):
    from multiprocessing.queues import Queue
    import errno

    while True:
        try:
            return function(*args, **kw)
        except IOError, e:
            if e.errno == errno.EINTR:
                continue
            else:
                raise

def _my_queue_get(queue, block=True, timeout=None):
    import errno
    while True:
        try:
            return queue.get(block, timeout)
        except IOError, e:
            if e.errno != errno.EINTR:
                raise

class KeyboardInterruptError(Exception): pass

'''
Region functions
'''
def read_ds9_region(filename):

    ''' Converts DS9 region file into format for plotting region.

    Need the following format:
        angle : degrees
        xy : pixels
        width : pixels
        height : pixels

    Region file provides following format:
        # Region file format: DS9 version 4.1
        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        fk5
        box(4:17:04.740,+29:20:31.32,5854.33",11972.7",130) # text={test}

    pyregion module reads DS9 regions:
    http://leejjoon.github.io/pyregion/users/overview.html


    '''

    # Import external modules
    import pyregion as pyr

    # Read region file
    try:
        region = pyr.open(filename)
    except IOError:
        return None

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    return region

def load_ds9_core_region(cores, filename='',
        header=None):

    from myimage_analysis import get_pix_coords

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    regions = read_ds9_region(filename)

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        core = tag[tag.find('{')+1:tag.find('}')]

        if core in cores:
            # Format vertices to be 2 x N array
            poly_verts = []
            for i in xrange(0, len(region.coord_list)/2):
                poly_verts.append((region.coord_list[2*i],
                                   region.coord_list[2*i+1]))

            poly_verts_pix = []
            for i in xrange(0, len(poly_verts)):
                poly_verts_pix.append(get_pix_coords(ra=poly_verts[i][0],
                                          dec=poly_verts[i][1],
                                          header=header)[:-1][::-1].tolist())

            cores[core]['poly_verts'] = {}
            cores[core]['poly_verts']['wcs'] = poly_verts
            cores[core]['poly_verts']['pixel'] = poly_verts_pix

    return cores

def load_wedge_core_region(cores, filename='', header=None):

    from myimage_analysis import get_pix_coords
    import pickle
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from astropy.io import fits
    from astropy.wcs import WCS

    # Create WCS object
    wcs_header = WCS(header)

    with open(filename, 'rb') as f:
        region_dict = pickle.load(f)

    for core_name in region_dict:
        if core_name in cores:
            core = region_dict[core_name]
            # Format vertices to be 2 x N array
            poly_verts = np.array((core['ra'], core['dec']))

            # Make a galactic coords object and convert to Ra/dec
            coords_fk5 = SkyCoord(core['ra'] * u.deg,
                                  core['dec'] * u.deg,
                                  frame='fk5',
                                  )
            # convert to pixel
            coords_pixel = np.array(coords_fk5.to_pixel(wcs_header))

            #print coords_pixel
            #print coords_pixel.shape

            # write data to dataframe
            poly_verts_pix = np.array((coords_pixel[1], coords_pixel[0])).T

            #print poly_verts_pix.shape

            #poly_verts_pix = np.array((core['ypix'], core['xpix'])).T

            cores[core_name]['poly_verts'] = {}
            cores[core_name]['poly_verts']['wcs'] = poly_verts
            cores[core_name]['poly_verts']['pixel'] = poly_verts_pix

    return cores

def get_cores_to_plot():

    '''

    '''

    # Which cores to include in analysis?
    cores_to_keep = [ # N(H2) cores
                     'G168.54-6.22',
                     'G168.12-6.42',
                     'G166.91-7.76',
                     'G165.36-7.51',
                     'G164.70-7.63',
                     'G164.65-8.12',
                     'G165.71-9.15',
                     'G164.99-8.60',
                     'G164.26-8.39',
                     'G164.18-8.84',
                     'G174.40-13.45',
                     'G174.70-15.47',
                     'G174.05-15.82',
                     'G171.49-14.91',
                     'G172.93-16.73',
                     'G171.00-15.80',
                     'G172.12-16.94',
                     'G171.14-17.57',
                     'G169.32-16.17',
                     'G168.10-16.38',
                     'G160.49-16.81',
                     'G160.46-17.99',
                     'G159.80-18.49',
                     'G160.14-19.08',
                     'G160.53-19.73',
                     'G159.19-20.11',
                     'G159.17-21.09',
                     'G158.39-20.72',
                     'G158.89-21.60',
                     'G158.26-21.81',
                     ]

    if 0: # Random cores
        cores_to_keep = [
                 'G166.83-8.68',
                 'G168.82-6.37',
                 'G168.05-7.01',
                 'G164.16-8.46',
                 'G165.23-8.78',
                 'G167.06-7.77',
                 'G168.12-6.42',
                 'G167.58-6.64',
                 'G164.70-7.63',
                 'G166.35-8.77',
                 'G166.73-15.06',
                 'G173.08-16.50',
                 'G172.74-14.53',
                 'G169.44-16.18',
                 'G173.86-17.65',
                 'G173.71-13.91',
                 'G171.75-14.18',
                 'G173.70-15.21',
                 'G170.28-19.48',
                 'G171.00-15.80',
                 'G158.23-20.15',
                 'G159.01-22.19',
                 'G159.19-20.11',
                 'G157.12-23.49',
                 'G160.10-19.90',
                 'G160.34-18.42',
                 'G158.40-21.86',
                 'G159.79-21.32',
                 'G158.89-21.60',
                 'G159.51-18.41',
                ]

    if 0:
        cores_to_keep = [# taur
                         'L1495',
                         'L1495A',
                         'B213',
                         'L1498',
                         'B215',
                         'B18',
                         'B217',
                         'B220-1',
                         'B220-2',
                         'L1521',
                         'L1524',
                         'L1527-1',
                         'L1527-2',
                         # Calif
                         'L1536',
                         'L1483-1',
                         'L1483-2',
                         'L1482-1',
                         'L1482-2',
                         'L1478-1',
                         'L1478-2',
                         'L1456',
                         'NGC1579',
                         #'L1545',
                         #'L1517',
                         #'L1512',
                         #'L1523',
                         #'L1512',
                         # Pers
                         'B5',
                         'IC348',
                         'B1E',
                         'B1',
                         'NGC1333',
                         'B4',
                         'B3',
                         'L1455',
                         'L1448',
                         ]

    return cores_to_keep

def get_core_properties(data_dict, cloud_name):

    from myimage_analysis import load_ds9_region, get_pix_coords
    import json

    box_method = 'ds9'
    core_dir = '/d/bip3/ezbc/multicloud/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'

    header = data_dict['av_header']
    if 0:
        # define core properties
        with open(core_dir + cloud_name + '_core_properties.txt', 'r') as f:
            cores = json.load(f)

        #cores = convert_core_coordinates(cores, header)

        # convert center WCS coordinate to pixel
        for core in cores:
            cores[core].update({'box_pixel': 0})
            cores[core].update({'center_pixel': 0})

            center_wcs = cores[core]['center_wcs']

            # convert centers to pixel coords
            center_pixel = get_pix_coords(ra=center_wcs[0],
                                          dec=center_wcs[1],
                                          header=header)[:2]
            cores[core]['center_pixel'] = center_pixel

    else:
        filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
            'planck11_coldclumps.pickle'

        import pickle
        with open(filename, 'rb') as f:
            core_sample = pickle.load(f)
        #core_sample = pd.load(filename)
        cores = {}
        df_cores = core_sample[cloud_name]
        for name in df_cores['Name'].values:
            cores[name] = {}
            ra = df_cores[df_cores['Name'] == name]['ra'].values[0]
            dec = df_cores[df_cores['Name'] == name]['dec'].values[0]
            cores[name]['center_wcs'] = (ra, dec)
            cores[name]['temp'] = \
                df_cores[df_cores['Name'] == name]['temp'].values[0]
            cores[name]['temp_error'] = \
                df_cores[df_cores['Name'] == name]['temp_error'].values[0]

            #print cores[name]['temp_error']

            # convert centers to pixel coords
            center_pixel = get_pix_coords(ra=ra,
                                          dec=dec,
                                          header=header)[:2]
            cores[name]['center_pixel'] = center_pixel[::-1]

    # load the bounding regions
    if 0:
        region_filename = region_dir + 'multicloud_coldclump_divisions.reg'
        cores = load_ds9_core_region(cores,
                                filename=region_filename,
                                header=header)
    else:
        filename = region_dir + 'multicloud_divisions_coldcore_wedges.pickle'
        cores = load_wedge_core_region(cores,
                                       filename=filename,
                                       header=header)

    # add indices of core in data
    add_core_mask(cores, data_dict['av_data'])

    return cores

def trim_cores_to_plot(cores, cores_to_plot):

    # Trim down cores to keep list to include only cores for cloud
    cores_to_keep_old = list(cores_to_plot)
    for core in cores_to_keep_old:
        if core not in cores or 'poly_verts' not in cores[core]:
            cores_to_plot.remove(core)

    return cores_to_plot

def add_core_mask(cores, data):

    for core in cores:
        try:
            vertices = cores[core]['poly_verts']['pixel']

            mask = np.logical_not(myg.get_polygon_mask(data,
                                                       vertices))

            cores[core]['indices_orig'] = np.where(mask == 0)

            cores[core]['mask'] = mask

            #print 'npix in core ' + core + ':'
            #print np.where(mask == 0)[0].size

        except KeyError:
            cores[core]['mask'] = None
            cores[core]['indices_orig'] = None

'''
Modeling Functions
'''

def refit_data(h_sd, rh2, h_sd_error=None, rh2_error=None, model_kwargs=None,
        h_sd_fit=None):

    import mystats

    #if h_sd_fit is not None and np.size(h_sd_fit) == np.size(h_sd):
        #data_array = h_sd, rh2, h_sd_error, rh2_error, h_sd_fit
        #h_sd, rh2, h_sd_error, rh2_error, h_sd_fit = mask_nans(data_array)
    #else:
    if 1:
        data_array = h_sd, rh2, h_sd_error, rh2_error
        h_sd, rh2, h_sd_error, rh2_error = mask_nans(data_array)

    if 0:
        print 'h_sd', h_sd
        print 'h_sd_error', h_sd_error
        print 'rh2', rh2
        print 'rh2_error', rh2_error

    ss_model_result = \
        fit_steady_state_models(h_sd.ravel(),
                                rh2.ravel(),
                                rh2_error=rh2_error.ravel(),
                                h_sd_error=h_sd_error.ravel(),
                                model_kwargs=model_kwargs,
                                )
    fitted_models = {}
    for model in ss_model_result:
        params = {}
        for param in ss_model_result[model]:
            params[param] = ss_model_result[model][param]

        if 0:
            params['phi_cnm'] = 3
            params['alphaG'] = 10
            params['phi_g'] = 3
            params['phi_mol'] = 10
            params['Z'] = 1



        # Calculate sum of residuals
        if 'sternberg' in model:
            model_fits = calc_sternberg(params,
                                      h_sd_extent=(0.001, 200),
                                      h_sd=h_sd,
                                      return_fractions=False,
                                      return_hisd=True)
        elif 'krumholz' in model:
            model_fits = calc_krumholz(params,
                                      h_sd_extent=(0.001, 200),
                                      h_sd=h_sd,
                                      return_fractions=False,
                                      return_hisd=True)

        # calculate resid sum of squares and log likelihood
        fitted_models[model] = {}
        fits = fitted_models[model]
        fits['sum_of_resid'] = np.sum((rh2 - model_fits[0])**2)
        fits['logL'] = mystats.calc_logL(model_fits[0], rh2, rh2_error)


        # use fitted h_sd
        # ---------------
        if 'sternberg' in model:
            model_fits = calc_sternberg(params,
                                      h_sd_extent=(0.001, 200),
                                      h_sd=h_sd_fit,
                                      return_fractions=False,
                                      return_hisd=True)
        elif 'krumholz' in model:
            model_fits = calc_krumholz(params,
                                      h_sd_extent=(0.001, 200),
                                      h_sd=h_sd_fit,
                                      return_fractions=False,
                                      return_hisd=True)

        fits['rh2'] = model_fits[0]
        fits['h_sd'] = model_fits[1]
        fits['hi_sd'] = model_fits[2]

        # calculate bayesian information criterion
        k = 1 # number of parameters
        N = np.size(rh2)
        BIC = k * np.log(N) - 2 * fits['logL']
        fits['BIC'] = BIC

    return fitted_models

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

    # sigma_g = 1.9 * 10^-21 phi_g * Z cm^2
    # sigma_g,gal is relative to galactic, and so is phi_g in this case, so
    # sigma_g,gal propto phi_g.
    sigma_g = phi_g
    sternberg_results['sf_threshold'] = bialy16.calc_sf_threshold(alphaG,
                                                                  sigma_g)

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
        h_sd = np.linspace(h_sd_extent[0], h_sd_extent[1], 1e3)

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

    #if not odr_fit:
    if 0:
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

            rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol,)

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

    if not odr_fit:
        h_sd_error = None
        rh2_error = None

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
        h_sd = np.linspace(h_sd_extent[0], h_sd_extent[1], 1e3)

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
                          method='leastsq',
                          #method='anneal',
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

    if not odr_fit:
        h_sd_error = None
        rh2_error = None

    alphaG, Z, phi_g = guesses
    def odr_func(alphaG, h_sd):
        return s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                 return_fractions=False)

    h_sd_error, rh2_error = None, None
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




'''
Bootstrapping functions
'''
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

    bootstrap_name = '{0:.0f}mcsim'.format(global_args['num_bootstraps'])

    filename_extension = global_args['cloud_name'] + '_' + global_args['data_type'] + \
            background_name + \
            bin_name + weights_name + \
            region_name + width_name + avthres_name + \
            intercept_name + error_name + compsub_name + backdgr_name + \
            '_' + hi_range_name + '_' + global_args['radiation_type'] + \
            rotate_cores_name + vary_phi_g_name + '_' + \
            bootstrap_name

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
    global_args['av_error'] = np.abs(av_error)
    global_args['rh2'] = rh2_image
    global_args['rh2_error'] = np.abs(rh2_image_error)
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
    global_args['sim_hi_error'] = np.abs(sim_hi_error)
    global_args['ss_model_kwargs'] = ss_model_kwargs
    if sim_hi_error:
        global_args['hi_data'] = hi_data
        global_args['vel_axis'] = vel_axis
        global_args['vel_range'] = vel_range
        global_args['vel_range_error'] = vel_range_error
        #print 'vel_range_error', vel_range_error
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
                 'sf_threshold': empty(),
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

def scale_av_with_refav(av_data, av_reference, av_error_data, perform_mc=0):

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
                                                     size=av_data.shape)
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
    data_products['nh2_error'] = np.abs(nh2_image_error)
    data_products['h2_sd_error'] = np.abs(h2_sd_image_error)
    data_products['hi_sd_error'] = np.abs(hi_sd_image_error)
    data_products['h_sd_error'] = np.abs(h_sd_image_error)
    data_products['rh2_error'] = np.abs(rh2_image_error)

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

    # write N(HI) and N(H2) images for each cloud
    write_coldens_maps(results_dict, global_args)

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
    lm_plt.plot_spectra(hi_spectrum,
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

    FILENAME_HI_DICT = '/d/bip3/ezbc/multicloud/data/' + \
                        'python_output/' + \
                        'core_properties/hi_properties.pickle'

    print('Number of permutations to run: ' + str(len(permutations)))

    #for cloud in clouds:
    for permutation in permutations:
        global_args = {
                'cloud_name':permutation[0],
                'load': 1,
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
                'num_bootstraps': 10000,
                'num_resid_bootstraps': 100,
                'bootstrap_fit_residuals': False,
                'calculate_median_error': False,
                'multiprocess': 1,
                'radiation_type': permutation[12],
                'rotate_cores': permutation[13],
                'vary_phi_g': permutation[14],
                'filename_hi_props': FILENAME_HI_DICT,
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

            save_results(results[global_args['cloud_name']],
                         global_args['results_filename'])

            print('\n\tPlotting')
            lm_plt.plot_results(results[global_args['cloud_name']])

    plot_multicloud_results(results)

if __name__ == '__main__':
    main()



