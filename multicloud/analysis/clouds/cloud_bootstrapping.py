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
import local_module_science as lm_science

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
                  'sigma_d': analysis['sigma_d'],
                  'Z': analysis['Z'],
                  }
        h_sd, hi_sd = calc_krumholz(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[1:3]

        params = {'phi_cnm': phi_cnm_high,
                  'sigma_d': analysis['sigma_d'],
                  'Z': analysis['Z'],
                  }
        hi_sd_low = calc_krumholz(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[2]

        params = {'phi_cnm': phi_cnm_low,
                  'sigma_d': analysis['sigma_d'],
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
    h2sd_list = []
    hisd_list = []
    nhi_error_list = []
    nh2_error_list = []
    rh2_error_list = []
    hsd_error_list = []
    h2sd_error_list = []
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
    stats_list = {'krumholz_results':  {'sum_of_resid': [],
                                        'chisq_reduced': [],
                                        'BIC': [],
                                        },
                  'sternberg_results':  {'sum_of_resid': [],
                                        'chisq_reduced': [],
                                        'BIC': [],
                                        },
                  }
    for i, cloud_name in enumerate(('california', 'perseus', 'taurus')):
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
        av_error_list.append(data_products['av_error'])
        nhi_list.append(data_products['nhi'])
        nh2_list.append(data_products['nh2'])
        rh2_list.append(data_products['rh2'])
        hsd_list.append(data_products['h_sd'])
        h2sd_list.append(data_products['h2_sd'])
        hisd_list.append(data_products['hi_sd'])
        nhi_error_list.append(data_products['nhi_error'])
        nh2_error_list.append(data_products['nh2_error'])
        rh2_error_list.append(data_products['rh2_error'])
        hsd_error_list.append(data_products['h_sd_error'])
        h2sd_error_list.append(data_products['h2_sd_error'])
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
        model_kwargs = lm_science.scale_dust_areas(dgr_list[i], model_kwargs)
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
            core_names.append(core)

            # get the pixels for the core
            core_indices = cores[core]['indices_orig']

            # get the data for each core
            hisd = hisd_list[i][core_indices]
            hisd_error = hisd_error_list[i][core_indices]
            hsd = hsd_list[i][core_indices]
            hsd_error = hsd_error_list[i][core_indices]
            rh2 = rh2_list[i][core_indices]
            rh2_error = rh2_error_list[i][core_indices]


            # mask the data based on nans and negative RH2
            data_list = (hisd, hisd_error, hsd, hsd_error, rh2, rh2_error)
            [hisd, hisd_error, hsd, hsd_error, rh2, rh2_error] = \
                mask_nans(data_list)
            data_list = [hisd, hisd_error, hsd, hsd_error, rh2, rh2_error]


            if 0:
                print('number of neg rh2 points:', np.sum(rh2 < 0))

            hisd_core_list.append(hisd)
            hsd_core_list.append(hsd)
            rh2_core_list.append(rh2)
            hisd_error_core_list.append(hisd_error)
            hsd_error_core_list.append(hsd_error)
            rh2_error_core_list.append(rh2_error)
            hisd_median_error_core_list.append(hisd_median_error_list[i])
            hsd_median_error_core_list.append(hsd_median_error_list[i])

            model_fits_list.append(refit_data(hsd_core_list[j],
                                              rh2_core_list[j],
                                              h_sd_error=hsd_error_core_list[j],
                                              rh2_error=rh2_error_core_list[j],
                                              model_kwargs=model_kwargs,
                                  hi_sd_error=hisd_error_core_list[j],
                                  hi_sd=hisd_core_list[j],
                                  odr_fitting=global_args['odr_fitting'],
                                              )
                                   )

            if 0:
                print('hi sd error',
                        scipy.stats.nanmedian(hisd_error_core_list[j].ravel()))
            if 0:
                if core in ('G174.70-15.47',):
                    print('')
                    print(core)
                    print('Number of pixels: {0:.0f}'.format(len(hisd)))

            # get the residual sum of squares
            # -------------------------------------------------------------------
            fits = refit_data(hsd_core_list[j],
                              rh2_core_list[j],
                              rh2_error=rh2_error_core_list[j],
                              h_sd_error=hsd_error_core_list[j],
                              h_sd_fit=hsd_core_list[j],
                              hi_sd_error=hisd_error_core_list[j],
                              hi_sd=hisd_core_list[j],
                              model_kwargs=model_kwargs,
                              odr_fitting=global_args['odr_fitting'],
                              )
            for model in model_fits_list[j]:
                if core in ('G160.53-19.73','G172.93-16.73','G164.70-7.63'):
                    rh2 = rh2_core_list[j]
                    rh2_error = rh2_error_core_list[j]

                    if 0:
                        print('core', core)
                        print('rh2 median:',scipy.stats.nanmedian(rh2))
                        print('rh2 error median:',scipy.stats.nanmedian(rh2_error))
                        print('model fits')
                        print(np.sum((rh2 - \
                            model_fits_list[j]['krumholz_results']['rh2_ind'])**2))

                stats_list[model]['sum_of_resid'].append(\
                        fits[model]['sum_of_resid'])
                stats_list[model]['chisq_reduced'].append(\
                        fits[model]['chisq_reduced'])
                stats_list[model]['BIC'].append(\
                        fits[model]['BIC'])
                #print 'K+09 core params for ' + core
                #print fits['krumholz_results']['params']
                #print 'S+14 core params for ' + core
                #print fits['sternberg_results']['params']

            if 0:
                rh2_copy = rh2_list[i].copy()
                rh2_copy[core_indices] = 1000
                plt.imshow(rh2_copy, origin='lower')
                plt.savefig('/d/bip3/ezbc/scratch/core_' + core + '.png')


            write_core_values_to_csv(core,
                                     hisd=hisd_core_list[j],
                                     hisd_error=hisd_error_core_list[j],
                                     rh2=rh2_core_list[j],
                                     rh2_error=rh2_error_core_list[j],
                                     hsd=hsd_core_list[j],
                                     hsd_error=hsd_error_core_list[j],
                                     fits=model_fits_list[j],
                                     )

        #print('hi sd median error', hisd_median_error_cores_list)

        cloud_model_fits_list.append(model_fits_list)
        core_names_list.append(core_names)
        hisd_cores_list.append(hisd_core_list)
        hsd_cores_list.append(hsd_core_list)
        rh2_cores_list.append(rh2_core_list)
        hisd_error_cores_list.append(hisd_error_core_list)
        hsd_error_cores_list.append(hsd_error_core_list)
        hisd_median_error_cores_list.append(hisd_median_error_core_list)
        hsd_median_error_cores_list.append(hsd_median_error_core_list)

    # print difference between "observed" data and median simulated parameters
    print_modelparam_diff(core_names_list,
                          cloud_model_fits_list,
                          model_analysis_list)

    # calculate statistics on hi
    hi_dict = calc_hi_statistics(cloud_name_list, core_names_list,
                                 hisd_cores_list, hsd_cores_list,
                                 rh2_cores_list, model_analysis_dict,
                                 filename=global_args['filename_hi_props'],
                                 stats_list=stats_list,
                                 )

    # plot CDFs of diffuse fraction LOS
    lm_plt.plot_diffusefraction_cdfs(hi_dict)

    # Print results
    # =========================================================================
    # Write results to a
    #print_av_error_stats(av_list[0], av_error_list[0])

    #print_BIC_results(stats_list, core_names_list)
    #print_RSS_results(stats_list, core_names_list)
    filename = results_dir + 'tables/cloud_obs_stats.csv'
    write_cloud_stats(cloud_name_list,
                      filename,
                      av_list=av_list,
                      rh2_list=rh2_list,
                      hisd_list=hisd_list,
                      hisd_error_list=hisd_error_list,
                      hsd_list=hsd_list,
                      h2sd_error_list=h2sd_error_list,
                      h2sd_list=h2sd_list,
                      nhi_list=nhi_list,
                      nh2_list=nh2_list,
                      nhi_error_list=nhi_error_list,
                      nh2_error_list=nh2_error_list,
                      )

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

    filename = results_dir + 'tables/multicloud_vel_ranges.tex'
    write_hi_vel_range_table(cloud_name_list,
                             hi_range_kwargs_list,
                             filename,
                             dgr_list=dgr_list,
                             dgr_error_list=dgr_error_list)

    # write dust properties table
    filename = results_dir + 'tables/dust_properties.tex'
    write_dust_props(cloud_name_list,
                     filename,
                     dgr_list=dgr_list,
                     dgr_error_list=dgr_error_list)


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

        # Plot RH2 vs. HSD
        # ---------------------------------------------------------------------
        filename = plot_kwargs['figure_dir'] + \
                   'av_nhi/multicloud_rh2_vs_hsd.' + filetype
        lm_plt.plot_rh2_vs_hsd_cloud(hsd_list,
                            rh2_list,
                            names_list=cloud_name_list,
                            filename=filename,
                            levels=(0.999, 0.99, 0.95, 0.87, 0.61,),
                            #levels=7,
                            #limits=[-5, 100, -1.5, 20],
                            scale=['linear', 'linear'],
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
                                #limit=[2, 20, -13, 20]
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

def write_cloud_stats(cloud_name_list,
                      filename,
                      av_list=None,
                      rh2_list=None,
                      hisd_list=None,
                      hsd_list=None,
                      hisd_error_list=None,
                      h2sd_error_list=None,
                      h2sd_list=None,
                      nhi_list=None,
                      nh2_list=None,
                      nhi_error_list=None,
                      nh2_error_list=None,
                      ):

    from myscience import calc_physical_scale

    f = open(filename, 'w')
    distances = {'california': 450.0,
                 'perseus': 240.0,
                 'taurus': 150.0,
                 }

    for i, cloud in enumerate(cloud_name_list):

        # prep data
        h2sd = h2sd_list[i]
        nh2 = nh2_list[i]
        h2sd_error = h2sd_error_list[i]
        hsd = hsd_list[i]
        hisd = hisd_list[i]
        av = av_list[i]
        (h2sd, nh2, h2sd_error, hsd, hisd, av) =\
            mask_nans((h2sd, nh2, h2sd_error, hsd,
                       hisd, av))

        h2sd_neg_frac = \
            float(np.nansum(h2sd < 0.0)) / np.size(h2sd)

        if any(hsd < h2sd):
            raise ValueError('HSD cannot be less than H2SD')

        # goldsmith et al. (2008) found half of mass of H2 below N(H2) < 2.1 x
        # 10^21 cm^-2.
        threshold = 2.1*10.0 # cm^-2
        nh2_diffuse_frac = \
            float(np.nansum(nh2 < threshold)) / np.size(nh2)
        hsd_diffuse_sum = np.nansum(hsd[nh2 < threshold])
        hsd_diffuse_frac = np.nansum(h2sd[nh2 < threshold]) / np.nansum(hsd)
        h2sd_diffuse_frac = np.nansum(h2sd[nh2 < threshold]) / np.nansum(h2sd)

        # calculate mass of cloud
        D = distances[cloud]
        physical_scale = calc_physical_scale(5.0 * 60.0, D)
        threshold_aks = [0.1, 0.2, -1e10, 0.8]
        hi_masses = []
        h2_masses = []
        region_sizes = []
        for threshold_ak in threshold_aks:
            # convert from A_K to A_V:
            # scalar reference: http://adsabs.harvard.edu/abs/2005ApJ...619..931I
            mask = av < 8.8 * threshold_ak
            pixel_size = physical_scale**2 # pc^2
            region_sizes.append(physical_scale**2 * np.size(hisd[~mask]))# pc^2
            h2_masses.append(pixel_size * np.sum(h2sd[~mask]) / 10**3)
            hi_masses.append(pixel_size * np.sum(hisd[~mask]) / 10**3)

        # calculate median errors
        h2sd_error_median = scipy.stats.nanmedian(np.ravel(h2sd_error_list[i]))
        hisd_error_median = scipy.stats.nanmedian(np.ravel(hisd_error_list[i]))
        nhi_error_median = scipy.stats.nanmedian(np.ravel(nhi_error_list[i]))
        nh2_error_median = scipy.stats.nanmedian(np.ravel(nh2_error_list[i]))

        h2sd_negerror_frac = \
            float(np.nansum(h2sd < -h2sd_error)) / np.size(h2sd)
        h2sd_negerror_3sig_frac = \
            float(np.nansum(h2sd < -3*h2sd_error)) / np.size(h2sd)

        f.write('cloud: ' + cloud)
        f.write("\r")
        # Fractions
        # -----------------------------------------------------------------------
        f.write('fraction N(H2) below 2.1*10^21 cm^-2: ' + \
                '{0:.2f}'.format(nh2_diffuse_frac))
        f.write("\r")
        f.write('Total surf dense for N(H2) below 2.1*10^21 cm^-2 ' + \
                '[Msun /pc^2]: ' + \
                '{0:.2f}'.format(hsd_diffuse_sum))
        f.write("\r")
        f.write('Fraction of total mass where N(H2) below 2.1*10^21 cm^-2 ' + \
                ': ' + \
                '{0:.2f}'.format(hsd_diffuse_frac))
        f.write("\r")
        f.write('Fraction of H2 mass where N(H2) below 2.1*10^21 cm^-2 ' + \
                ': ' + \
                '{0:.2f}'.format(h2sd_diffuse_frac))
        f.write("\r")
        f.write('fraction H2 below 1 sigma: {0:.2f}'.format(h2sd_negerror_frac))
        f.write("\r")
        f.write('fraction H2 below 3 sigma: ' + \
                '{0:.2f}'.format(h2sd_negerror_3sig_frac))
        f.write("\r")
        f.write('fraction negative H2: {0:.2f}'.format(h2sd_neg_frac))
        f.write("\r")

        # Masses
        # -----------------------------------------------------------------------
        for i, threshold_ak in enumerate(threshold_aks):
            if threshold_ak > 0:
                f.write('H2 mass for Ak > {0:.1f}'.format(threshold_ak) + \
                        'mag [10^3 Msun]: ' + \
                        '{0:.2f}'.format(h2_masses[i]))
                f.write("\r")
                f.write('HI mass for Ak > {0:.1f}'.format(threshold_ak) + \
                        'mag [10^3 Msun]: ' + \
                        '{0:.2f}'.format(hi_masses[i]))
                f.write("\r")
                f.write('Total mass for Ak > {0:.1f}'.format(threshold_ak) + \
                        'mag [10^3 Msun]: ' + \
                        '{0:.2f}'.format(h2_masses[i] + hi_masses[i]))
                f.write("\r")
                f.write('Region size for Ak > {0:.1f} '.format(threshold_ak) + \
                        '[pc^2]: ' + \
                        '{0:.2f}'.format(region_sizes[i]))
                f.write("\r")
            else:
                f.write('H2 mass of whole region [10^3 Msun]: ' + \
                        '{0:.2f}'.format(h2_masses[i]))
                f.write("\r")
                f.write('HI mass of whole region [10^3 Msun]: ' + \
                        '{0:.2f}'.format(hi_masses[i]))
                f.write("\r")
                f.write('Total H mass of whole region [10^3 Msun]: ' + \
                        '{0:.2f}'.format(h2_masses[i] + hi_masses[i]))
                f.write("\r")
                f.write('Region size for whole region' + \
                        '[pc^2]: ' + \
                        '{0:.2f}'.format(region_sizes[i]))
                f.write("\r")

        # Uncertainties
        # -----------------------------------------------------------------------
        f.write('median random uncertainty of HI SD [msun / pc^2]: ' + \
                '{0:.2g} '.format(hisd_error_median))
        f.write("\r")
        f.write('median random uncertainty of H2 SD [msun / pc^2]: ' + \
                '{0:.2g} '.format(h2sd_error_median))
        f.write("\r")
        f.write('median random uncertainty of N(HI) [10^20 cm^-2]: ' + \
                '{0:.2g} '.format(nhi_error_median))
        f.write("\r")
        f.write('median random uncertainty of N(H2) [10^20 cm^-2]: ' + \
                '{0:.2g} '.format(nh2_error_median))
        f.write("\r")
        f.write("\r")

    f.close()


def collect_results(results_dict):

    clouds_dict = {}
    final_product = {}
    final_product['clouds'] = {}

    for i, cloud_name in enumerate(results_dict):
        results = results_dict[cloud_name]
        clouds_dict[cloud_name] = {}
        cloud = clouds_dict[cloud_name]

        cloud['filenames'] = results['filenames']
        cloud['plot_kwargs'] = results['plot_kwargs']
        cloud['data_products'] = results['data_products']
        cloud['mc_results'] = results_dict['mc_results']
        cloud['mc_analysis'] = results_dict['mc_analysis']
        cloud[''] = results['']

def write_core_values_to_csv(core_name, hisd=None, hsd=None, rh2=None,
        rh2_error=None, hisd_error=None, hsd_error=None, fits=None,):

    import pandas as pd

    # core values
    filename = '/d/bip3/ezbc/multicloud/data/python_output/tables/cores/' + \
               'properties_' + core_name + '.csv'

    data = [hisd, hisd_error, hsd, hsd_error, rh2, rh2_error,
            fits['sternberg_results']['hisd_ind'],
            fits['sternberg_results']['hsd_ind'],
            fits['sternberg_results']['rh2_ind'],
            fits['krumholz_results']['hisd_ind'],
            fits['krumholz_results']['hsd_ind'],
            fits['krumholz_results']['rh2_ind'],
            ]

    columns = ('HI_SD', 'HI_SD_ERR', 'H_SD', 'H_SD_ERR', 'RH2', 'RH2_ERR',
               'HI_SD_S14FIT', 'H_SD_S14FIT', 'RH2_S14FIT',
               'HI_SD_K09FIT', 'H_SD_K09FIT', 'RH2_K09FIT',
               )

    #if core_name in ('G160.53-19.73','G172.93-16.73','G164.70-7.63'):
        #print('rh2 median:',scipy.stats.nanmedian(rh2))
        #print('rh2 error median:',scipy.stats.nanmedian(rh2_error))
        #print('hsd median:',scipy.stats.nanmedian(hsd))
        #print('hsd error median:',scipy.stats.nanmedian(hsd_error))

    # mask nans from data
    data[:6] = mask_nans(data[:6])

    data_dict = {}
    for i, column in enumerate(columns):
        data_dict[column] = data[i]

    # WRite the data to a dataframe
    df = pd.DataFrame(data_dict, columns=columns)

    # csv
    df.to_csv(filename, index=False)


    # core fitted parameters
    filename = '/d/bip3/ezbc/multicloud/data/python_output/tables/cores/' + \
               'modelparams_' + core_name + '.csv'

    data = [
            fits['sternberg_results']['params']['alphaG'],
            fits['sternberg_results']['params']['phi_g'],
            fits['krumholz_results']['params']['phi_cnm'],
            fits['krumholz_results']['params']['sigma_d'],
            ]

    columns = ('ALHPHAG','PHI_G','PHI_CNM','SIGMA_D'
               )

    data_dict = {}
    for i, column in enumerate(columns):
        data_dict[column] = [data[i],]

    # WRite the data to a dataframe
    df = pd.DataFrame(data_dict, columns=columns)

    # csv
    df.to_csv(filename, index=False)

def print_modelparam_diff(cores_list, fits_list, model_analysis_list,):

    #print('\n\tDifferences between fitted parameters to observed data and MC' +\
            #'median data:')
    for i, cloud in enumerate(cores_list):
        for j, core_name in enumerate(cores_list[i]):
            core = model_analysis_list[i]['cores'][core_name]
            #print('\n\t' + core_name)

            for model in ('krumholz', 'sternberg'):
                #print fits_list[i][j][model + '_results']['params']
                for param_name in fits_list[i][j][model + '_results']['params']:
                    param_obs = \
                        fits_list[i][j][model + \
                        '_results']['params'][param_name]
                    param_med = core[model + '_results'][param_name]

                    diff = param_obs - param_med
                    #print('\t' + param_name + ': {0:.3f}'.format(diff))


def print_BIC_results(stats_list, core_names):

    bayes_factor = np.array(stats_list['krumholz_results']['BIC']) - \
                   np.array(stats_list['sternberg_results']['BIC'])

    chisq_s14 = np.array(stats_list['sternberg_results']['chisq_reduced'])
    chisq_k09 = np.array(stats_list['krumholz_results']['chisq_reduced'])

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

    print('Confidence intervals on reduced chi squared')
    print('krumholz')
    print(mystats.calc_cdf_error(chisq_k09, alpha=0.5))
    print(np.sort(chisq_k09))
    print('sternberg')
    print(mystats.calc_cdf_error(chisq_k09, alpha=0.5))
    print(np.sort(chisq_s14))
    print('diff')
    print(np.sort(chisq_k09 - chisq_s14))

    print('Core with max BF of {0:.0f}'.format(np.max(bayes_factor)))
    print(core_names[np.argmax(bayes_factor)])
    print('Core with min BF of {0:.0f}'.format(np.min(bayes_factor)))
    print(core_names[np.argmin(bayes_factor)])

    print('Core BF')
    index = np.where(core_names == 'G158.39-20.72')[0]
    print(core_names[index], bayes_factor[index])

    plt.close(); plt.clf()
    myplt.plot_cdf(bayes_factor)
    plt.xlabel('K+09 - S+14')
    plt.title('Bayes Factor CDF')
    plt.savefig('/d/bip3/ezbc/multicloud/figures/models/bayes_factor_cdf.png')

def print_RSS_results(stats_list, core_names):

    RSS_s14 = np.array(stats_list['sternberg_results']['sum_of_resid'])
    RSS_k09 = np.array(stats_list['krumholz_results']['sum_of_resid'])
    RSS_diff = RSS_k09 - RSS_s14

    # cores gathers by cloud, unravel
    core_names = np.ravel(core_names)

    print('Median difference in BIC between k09 and s14 models:')
    print(np.median(RSS_k09 - RSS_s14))

    import matplotlib.pyplot as plt
    import myplotting as myplt
    import mystats

    print('Confidence intervals on RSS K09 - S14 difference')
    print(mystats.calc_cdf_error(RSS_diff, alpha=0.5))
    print(np.sort(RSS_diff))

    print('Core with max BF of {0:.0f}'.format(np.max(RSS_diff)))
    print(core_names[np.argmax(RSS_diff)])
    print('Core with min BF of {0:.0f}'.format(np.min(RSS_diff)))
    print(core_names[np.argmin(RSS_diff)])

    print('Core BF')
    index = np.where(core_names == 'G158.39-20.72')[0]
    print(core_names[index], RSS_diff[index])

    plt.close(); plt.clf()
    myplt.plot_cdf(RSS_diff)
    plt.xlabel('K+09 - S+14')
    plt.title('Bayes Factor CDF')
    plt.savefig('/d/bip3/ezbc/multicloud/figures/models/RSS_diff_cdf.png')


def calc_hi_statistics(cloud_name_list, core_names_list,
                                 hisd_cores_list, h_sd_cores_list,
                                 rh2_cores_list, model_analysis_dict,
                                 filename=None, stats_list=None):

    hi_dict = {}
    for i, cloud in enumerate(cloud_name_list):

        # initialize hi properties
        hi_dict[cloud] = {}
        hi_dict[cloud]['cores'] = []
        hi_dict[cloud]['hi_sd_mean'] = []
        hi_dict[cloud]['hi_sd_median'] = []
        hi_dict[cloud]['hi_sd_std'] = []
        hi_dict[cloud]['fraction_LOS_diffuse'] = []
        hi_dict[cloud]['chisq_reduced_krumholz'] = []
        hi_dict[cloud]['chisq_reduced_sternberg'] = []
        hi_dict[cloud]['rh2_neg_fraction'] = []

        # add a row for each core
        for j, core in enumerate(core_names_list[i]):
            hi = hisd_cores_list[i][j]
            h = h_sd_cores_list[i][j]
            rh2 = rh2_cores_list[i][j]
            rh2_neg_fraction = np.sum(rh2 < 0) / float(np.size(rh2))
            hi_dict[cloud]['cores'].append(core)
            hi_dict[cloud]['hi_sd_mean'].append(np.nanmean(hi))
            hi_dict[cloud]['hi_sd_median'].append(scipy.stats.nanmedian(hi))
            hi_dict[cloud]['hi_sd_std'].append(np.nanstd(hi))
            hi_dict[cloud]['rh2_neg_fraction'].append(rh2_neg_fraction)

            # frac of diffuse LOS
            indices_diffuse = np.where(rh2 < 1.0)[0]
            frac_diffuse = np.size(indices_diffuse) / float(np.size(rh2))
            hi_dict[cloud]['fraction_LOS_diffuse'].append(frac_diffuse)

            #frac above sf threshold
            if 0:
                model = model_analysis_dict[cloud]['cores'][core]
                sf_threshold = model['sternberg_results']['sf_threshold']
                hi_dict[cloud]['sf_threshold'].append(sf_threshold)
                indices_sf = np.where(h > sf_threshold)[0]
                frac_sf = np.size(indices_sf) / float(np.size(h))
                hi_dict[cloud]['fraction_LOS_sf'].append(frac_sf)

            if stats_list is not None:
                hi_dict[cloud]['chisq_reduced_krumholz'].append( \
                    stats_list['krumholz_results']['chisq_reduced'][j])
                hi_dict[cloud]['chisq_reduced_sternberg'].append( \
                    stats_list['sternberg_results']['chisq_reduced'][j])

                if 0:
                    print 'chisq:'
                    print stats_list['sternberg_results']['chisq_reduced'][j]
                    print stats_list['krumholz_results']['chisq_reduced'][j]

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

            # Write HI properties
            if hi_dict is not None:
                hi_cloud_dict = hi_dict[cloud]
                index = hi_cloud_dict['cores'].index(core_name)
                keys = hi_cloud_dict.keys()

                # remove cores from dict
                #keys.pop(keys.index('cores'))

                hi_core_dict = {}
                for param in hi_cloud_dict:
                    #print hi_cloud_dict[param]
                    hi_core_dict[param] = hi_cloud_dict[param][index]

                core_new['hi_props'] = hi_core_dict

            # append model params and errors to row
            for model in ('krumholz', 'sternberg'):
                if model == 'krumholz':
                    params_to_write = ['phi_cnm', 'Z', 'sigma_d',
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

def write_dust_props(names_list, filename,
        dgr_list=None, dgr_error_list=None,):

    # Open file to be appended
    f = open(filename, 'wb')

    row_text = ''

    for cloud in names_list:
        row_text = add_row_element(row_text, cloud.capitalize())
    row_text += ' \\\\[0.1cm] \n'
    f.write(row_text)

    text_format = r'{0:.1f}$^{{+{1:.1f}}}_{{-{2:.1f}}}$'

    # DGR
    # -----------------------------------------------------------------------
    row_text = r'DGR'
    for i in xrange(len(names_list)):
        dgr = dgr_list[i]
        dgr_error = dgr_error_list[i]

        row_text = \
            add_row_element(row_text,
                      (dgr*100, dgr_error[1]*100, dgr_error[0]*100),
                      text_format=text_format)
    row_text += ' \\\\[0.1cm] \n'
    f.write(row_text)

    # phi_g
    # -----------------------------------------------------------------------
    row_text = r'$\phi_g$'
    for i in xrange(len(names_list)):
        dgr = dgr_list[i]
        dgr_error = dgr_error_list[i]

        row_text = \
            add_row_element(row_text,
                      (dgr*100/5.3,
                       dgr_error[1]*100/5.3,
                       dgr_error[0]*100/5.3),
                      text_format=r'{0:.1f}$^{{+{1:.1f}}}_{{-{2:.1f}}}$')
    row_text += ' \\\\[0.1cm] \n'
    f.write(row_text)

    # sigma_g
    # -----------------------------------------------------------------------
    row_text = r'$\sigma_g$'
    for i in xrange(len(names_list)):
        dgr = dgr_list[i]
        dgr_error = dgr_error_list[i]

        row_text = \
            add_row_element(row_text,
                      (dgr*100/5.3 * 1.9,
                       dgr_error[1]*100/5.3 * 1.9,
                       dgr_error[0]*100/5.3 * 1.9),
                      text_format=r'{0:.1f}$^{{+{1:.1f}}}_{{-{2:.1f}}}$')
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
    params_to_write = ['phi_cnm', 'Z', 'sigma_d', 'alphaG', 'Z', 'phi_g']

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
                    params_to_write = ['phi_cnm', 'Z', 'sigma_d',]
                else:
                    params_to_write = ['alphaG', 'Z', 'phi_g']
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
    text_param_format_chisq ='{0:.2g}'

    params_to_write = ['hi_sd_mean', 'hi_sd_median', 'hi_sd_std',
    'fraction_LOS_diffuse', 'rh2_neg_fraction'] #'chisq_reduced_krumholz', 'chisq_reduced_sternberg']

    # Collect parameter names for each model for each core
    row = 0
    for cloud in ('california', 'perseus', 'taurus'):
        core_dict = hi_dict[cloud]
        core_indices = np.argsort(core_dict['cores'])

        row_core = 0

        # each core correspond to a row
        for core_index in core_indices:

            core_name = core_dict['cores'][core_index]

            if row_core == 0:
                row_text = cloud.capitalize()
            else:
                row_text = ''

            row_text = add_row_element(row_text,
                                       core_name)

            for i, param_name in enumerate(params_to_write):
                param = \
                    core_dict[param_name][core_index]

                param_info = param

                if 'fraction' in param_name:
                    text_param_format = text_param_format_frac
                    #print param_info
                    param_info = param_info * 100.0
                elif 'chisq' in param_name:
                    text_param_format = text_param_format_chisq
                else:
                    text_param_format = text_param_format_sd

                #if param_name == 'alphaG':
                    #print core_name, param_info

                row_text = \
                    add_row_element(row_text,
                                    param_info,
                                    text_format=text_param_format)


            row_text += ' \\\\[0.1cm] \n'
            if row_core == len(core_indices) - 1 \
                and cloud != 'taurus':
                row_text += '\hline  \\\\[-0.2cm] \n'
            elif row_core == len(core_indices) - 1 and \
                    cloud == 'taurus':
                row_text.replace(r'\\[0.1cm] \n', '')

            f.write(row_text)

            row_core += 1
            row += 1

    f.close()

def add_row_element(row_text, element, text_format='{0:s}'):

    if type(element) is list or type(element) is tuple:
        return row_text + ' & ' + text_format.format(*element)
    else:
        return row_text + ' & ' + text_format.format(element)

def write_final_maps(results_dict, global_args):


    ''' Writes final data products.

    '''

    from astropy.io import fits

    DIR_FIGURES = '/d/bip3/ezbc/multicloud/data/final_products/'
    FILENAME_BASE = DIR_FIGURES + global_args['cloud_name']

    FILENAME_AV = FILENAME_BASE + '_av.fits'
    FILENAME_AV_ERROR = FILENAME_BASE + '_av_error.fits'
    FILENAME_HI = FILENAME_BASE + '_hisurfdens.fits'
    FILENAME_H2 = FILENAME_BASE + '_h2surfdens.fits'
    FILENAME_H = FILENAME_BASE + '_hsurfdens.fits'
    FILENAME_NHI = FILENAME_BASE + '_hicoldens.fits'
    FILENAME_NH2 = FILENAME_BASE + '_h2coldens.fits'
    FILENAME_NH = FILENAME_BASE + '_hcoldens.fits'
    FILENAME_RH2 = FILENAME_BASE + '_rh2.fits'
    FILENAME_HI_ERROR = FILENAME_BASE + '_hisurfdens_error.fits'
    FILENAME_H_ERROR = FILENAME_BASE + '_hsurfdens_error.fits'
    FILENAME_H2_ERROR = FILENAME_BASE + '_h2surfdens_error.fits'
    FILENAME_NHI_ERROR = FILENAME_BASE + '_hicoldens_error.fits'
    FILENAME_NH_ERROR = FILENAME_BASE + '_hcoldens_error.fits'
    FILENAME_NH2_ERROR = FILENAME_BASE + '_h2coldens_error.fits'
    FILENAME_RH2_ERROR = FILENAME_BASE + '_rh2_error.fits'

    # get av_header
    header = results_dict['data']['av_header'].copy()
    av_image = results_dict['data_products']['av']
    av_error_image = results_dict['data_products']['av_error']
    hi_sd_image = results_dict['data_products']['hi_sd']
    h_sd_image = results_dict['data_products']['h_sd']
    h2_sd_image = results_dict['data_products']['h2_sd']
    nhi_image = results_dict['data_products']['nhi']
    nh_image = results_dict['data_products']['nh']
    nh2_image = results_dict['data_products']['nh2']
    rh2_image = results_dict['data_products']['rh2']
    hi_sd_error_image = results_dict['data_products']['hi_sd_error']
    h_sd_error_image = results_dict['data_products']['h_sd_error']
    h2_sd_error_image = results_dict['data_products']['h2_sd_error']
    nhi_error_image = results_dict['data_products']['nhi_error']
    nh_error_image = results_dict['data_products']['nh_error']
    nh2_error_image = results_dict['data_products']['nh2_error']
    rh2_error_image = results_dict['data_products']['rh2_error']

    # Av map
    fits.writeto(FILENAME_AV, av_image, header=header, clobber=True)
    fits.writeto(FILENAME_AV_ERROR, av_error_image, header=header, clobber=True)

    # Surface density maps
    header['BUNIT'] = 'Msun/pc2'
    fits.writeto(FILENAME_HI, hi_sd_image, header=header, clobber=True)
    fits.writeto(FILENAME_HI_ERROR, hi_sd_error_image, header=header, clobber=True)
    fits.writeto(FILENAME_H2, h2_sd_image, header=header, clobber=True)
    fits.writeto(FILENAME_H2_ERROR, h2_sd_error_image, header=header, clobber=True)
    fits.writeto(FILENAME_H, h_sd_image, header=header, clobber=True)
    fits.writeto(FILENAME_H_ERROR, h_sd_error_image, header=header, clobber=True)

    # column density maps
    header['BUNIT'] = '10e20/cm2'
    fits.writeto(FILENAME_NHI, nhi_image, header=header, clobber=True)
    fits.writeto(FILENAME_NHI_ERROR, nhi_error_image, header=header, clobber=True)
    fits.writeto(FILENAME_NH2, nh2_image, header=header, clobber=True)
    fits.writeto(FILENAME_NH2_ERROR, nh2_error_image, header=header, clobber=True)
    fits.writeto(FILENAME_NH, nh_image, header=header, clobber=True)
    fits.writeto(FILENAME_NH_ERROR, nh_error_image, header=header, clobber=True)

    # ratio map
    header['BUNIT'] = ''
    fits.writeto(FILENAME_RH2, rh2_image, header=header, clobber=True)
    fits.writeto(FILENAME_RH2_ERROR, rh2_error_image, header=header, clobber=True)

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
        h_sd_fit=None, hi_sd=None, hi_sd_error=None, odr_fitting=False,):

    import mystats

    #if h_sd_fit is not None and np.size(h_sd_fit) == np.size(h_sd):
        #data_array = h_sd, rh2, h_sd_error, rh2_error, h_sd_fit
        #h_sd, rh2, h_sd_error, rh2_error, h_sd_fit = mask_nans(data_array)
    #else:
    if 1:
        data_array = h_sd, rh2, h_sd_error, rh2_error, hi_sd, hi_sd_error
        h_sd, rh2, h_sd_error, rh2_error, hi_sd, hi_sd_error = \
            mask_nans(data_array)

    ss_model_result = \
        lm_science.fit_steady_state_models(h_sd.ravel(),
                                rh2.ravel(),
                                rh2_error=rh2_error.ravel(),
                                h_sd_error=h_sd_error.ravel(),
                                model_kwargs=model_kwargs,
                                odr_fit=odr_fitting,
                                )
    fitted_models = {}
    radiation_type = model_kwargs['sternberg_params']['radiation_type']
    for model in ss_model_result:
        params = {}
        for param in ss_model_result[model]:
            params[param] = ss_model_result[model][param]

        # Calculate sum of residuals
        if 'sternberg' in model:
            model_fits = calc_sternberg(params,
                                      h_sd_extent=(0.001, 200),
                                      h_sd=h_sd,
                                      return_fractions=False,
                                      radiation_type=radiation_type,
                                      return_hisd=True)
        elif 'krumholz' in model:
            model_fits = calc_krumholz(params,
                                      h_sd_extent=(0.001, 200),
                                      h_sd=h_sd,
                                      return_fractions=False,
                                      return_hisd=True)

        # calculate resid sum of squares and log likelihood
        rh2_fit, hsd_fit, hisd_fit = model_fits

        fitted_models[model] = {}
        fits = fitted_models[model]
        fits['dof'] = np.size(rh2) - 1
        fits['sum_of_resid'] = np.sum((rh2 - model_fits[0])**2)
        fits['logL'] = mystats.calc_logL(model_fits[0], rh2, rh2_error)
        fits['logL'] = mystats.calc_logL(hisd_fit, hi_sd, hi_sd_error)
        fits['rh2_ind'] = rh2_fit
        fits['hsd_ind'] = hsd_fit
        fits['hisd_ind'] = hisd_fit
        try:
            fits['chisq_reduced'] = \
                mystats.calc_chisq(rh2_fit, rh2, rh2_error,
                                   dof=fits['dof'])


            #print('hi_sd_error', np.median(hi_sd_error))
            #fits['chisq_reduced'] = \
                #mystats.calc_chisq(hisd_fit, hi_sd, hi_sd_error,
                                   #dof=fits['dof'])
        except ValueError:
            fits['chisq_reduced'] = np.nan

        # use fitted h_sd
        # ---------------
        if 'sternberg' in model:
            model_fits = calc_sternberg(params,
                                      h_sd_extent=(0.001, 200),
                                      h_sd=h_sd_fit,
                                      return_fractions=False,
                                      radiation_type=radiation_type,
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
        fits['params'] = params

        # calculate bayesian information criterion
        k = 1 # number of parameters
        N = np.size(rh2)
        BIC = k * np.log(N) - 2 * fits['logL']
        #BIC = N * np.log(RSS / N) + k * np.log(N)
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

def add_hi_transition_calc(ss_model_result, model_kwargs):

    h_sd_fit = np.linspace(0, 100, 1000)

    # To get HI transition, calculate model fits, then find where RH2 = 1
    radiation_type = model_kwargs['sternberg_params']['radiation_type']
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
                                      radiation_type=radiation_type,
                                      return_hisd=False,
                                      )
        elif 'krumholz' in model_name:
            params['phi_cnm'] = model['phi_cnm']
            params['Z'] = model['Z']
            params['sigma_d'] = model['sigma_d']
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

    params = [params['phi_cnm'], params['Z'], params['sigma_d']]
    #print 'calc_krumholz params', params

    if params[0] <= 0 or np.isnan(params[0]):
        rh2_fit, f_H2, f_HI = np.empty(1), np.empty(1), np.empty(1)
    else:
        rh2_fit, f_H2, f_HI = k09.calc_rh2(h_sd,
                                           phi_cnm=params[0],
                                           Z=params[1],
                                           sigma_d=params[2],
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

def calc_sternberg(params, h_sd_extent=(0.001, 500), return_fractions=True,
        return_hisd=False, h_sd=None, radiation_type='beamed',):

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
    #print 'calc_sternberg params', params

    if params[0] <= 0 or np.isnan(params[0]):
        rh2_fit, f_H2, f_HI = np.empty(1), np.empty(1), np.empty(1)
    else:
        rh2_fit, f_H2, f_HI = s14.calc_rh2(h_sd,
                                           alphaG=params[0],
                                           Z=params[1],
                                           phi_g=params[2],
                                           radiation_type=radiation_type,
                                           return_fractions=True)

    output = [rh2_fit, h_sd]

    if return_fractions:
        output.append(f_H2)
        output.append(f_HI)
    if return_hisd:
        hi_sd = f_HI * h_sd
        output.append(hi_sd)

    return output


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
    if global_args['odr_fitting']:
        odr_name = '_odrfit'
    else:
        odr_name = ''
    if global_args['param_vary'] == (True, True, False):
        fit_name = '_fit-phicnm+alphaG+Z'
    elif global_args['param_vary'] == (False, True, False):
        fit_name = '_fitZ'
    else:
        fit_name = ''

    dust_cross_section_name = global_args['dust_cross_section_type']
    bootstrap_name = '{0:.0f}mcsim'.format(global_args['num_bootstraps'])

    filename_extension = global_args['cloud_name'] + '_' + global_args['data_type'] + \
            background_name + \
            bin_name + weights_name + \
            region_name + width_name + avthres_name + \
            intercept_name + error_name + compsub_name + backdgr_name + \
            '_' + hi_range_name + '_' + global_args['radiation_type'] + \
            rotate_cores_name + vary_phi_g_name + '_' + \
            bootstrap_name + odr_name + fit_name + '_' + dust_cross_section_name

    return filename_extension, global_args

def mask_nans(arrays, return_mask=False):

    """ Masks any positions where any array in the list has a NaN.

    Parameters
    ----------
    arrays : tuple
        Tuple of arrays. The mask will be the shape of the first array. The
        last axes of the rest of the arrays will be masked.

    """

    mask = np.zeros(np.shape(arrays[0]), dtype=bool)

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
    if 1:
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

    # Calculate N(H2), then HI + H2 surf dens
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
    # ---------------------------------------------------------------------------
    # Galactic DGR = 5.3 x 10^-22 cm^2 mag
    # our DGR in units of 10^-20 cm^2 mag
    # Galactic DGR in our units = 5.3 x 10^-22 / 10^-20 = 5.3 x 10^-2 = 0.053
    DGR = av_model_results['dgr_cloud']
    new_model_kwargs = lm_science.scale_dust_areas(DGR, model_kwargs,
            dust_cross_section_type=global_args['dust_cross_section_type'])

    # cycle through each core, bootstrapping the pixels
    # ---------------------------------------------------------------------------
    for core in cores_to_plot:
        if rotate_cores:
            core_indices = get_rotated_core_indices(cores[core],
                                                    mask,
                                                    corename=core,
                                                    iteration=i,
                                                    )
        else:
            # grab the indices of the core in the unraveled array
            core_indices = cores[core]['indices']

        if 0:
            assert av[core_indices] == cores[core]['test_pix']

        if 0:
            #if core in ('158.26-21.81','174.70-15.47',):
            print core
            print('\n\tRegion size = ' + \
                  '{0:.0f} pix'.format(core_indices.size))

        # Bootstrap core indices
        #core_boot_indices = core_indices[index_ints]
        np.random.seed()
        core_boot_indices = np.random.choice(core_indices.size,
                                             size=core_indices.size,)
        # get bootstrapped pixels
        h_sd_core = h_sd_image[core_indices]
        h_sd_core_error = h_sd_image_error[core_indices]
        rh2_core = rh2_image[core_indices]
        rh2_core_error = rh2_image_error[core_indices]

        #if any(rh2_core) < 0:
        if 0:
            print 'fraction of neg rh2', np.sum(rh2_core < 0) / \
                float(np.size(rh2_core))
            if np.sum(rh2_core < 0) > 0:
                print core
                print 'number of R(H2) less than 0', np.sum(rh2_core < 0)

        # mask negative ratios
        if 0:
            mask_rh2 = (rh2_core < 0) | (np.isnan(rh2_core))
            print('rh2 neg size', np.sum(mask_rh2))
            rh2_core = rh2_core[~mask_rh2]
            rh2_core_error = rh2_core_error[~mask_rh2]
            h_sd_core = h_sd_core[~mask_rh2]
            h_sd_core_error = h_sd_core_error[~mask_rh2]

        # Fit the models
        # -----------------------------------------------------------------------
        ss_model_result = \
            lm_science.fit_steady_state_models(h_sd_core,
                                               rh2_core,
                                               rh2_error=rh2_core_error,
                                               h_sd_error=h_sd_core_error,
                                               model_kwargs=new_model_kwargs,
                                               odr_fit=\
                                                   global_args['odr_fitting'],
                                               )

        if plot_kwargs['plot_diagnostics']:
            filename = plot_kwargs['figure_dir'] + \
                       'diagnostics/models/' + plot_kwargs['filename_base'] + \
                       '_rh2_vs_h_bootstrap' + \
                       '{0:03d}.png'.format(plot_kwargs['bootstrap_num'])
            lm_plt.plot_rh2_vs_h_diagnostic(h_sd_core,
                                     rh2_core,
                                     h_sd_error=h_sd_core_error,
                                     rh2_error=rh2_core_error,
                                     model_results=ss_model_result,
                                     filename=filename)

        # get HI transition result
        add_hi_transition_calc(ss_model_result, new_model_kwargs)
        ss_model_results[core] = ss_model_result

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
        rotate_cores=False, calc_median_error=False, odr_fitting=False,
        dust_cross_section_type='sternberg',):

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

    print av.size

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
    global_args['odr_fitting'] = odr_fitting
    global_args['dust_cross_section_type'] = dust_cross_section_type

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
        lm_science.fit_steady_state_models(h_sd_core,
                                rh2_core,
                                rh2_error=rh2_core_error,
                                model_kwargs=model_kwargs,
                                bootstrap_residuals=True,
                                nboot=nboot,
                                odr_fit=global_args['odr_fitting'],
                                )



    # get HI transition result
    add_hi_transition_calc(ss_model_result, model_kwargs)

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
                 'sigma_d': empty(),
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

def get_model_fit_kwargs(cloud_name, vary_phi_g=False,
        param_vary=(True,False,False),):

    '''

    '''
    vary_alphaG = param_vary[0] # Vary alphaG in S+14 fit?
    vary_Z = param_vary[1] # Vary metallicity in S+14 fit?
    vary_phi_g = param_vary[2] # Vary phi_g in S+14 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    guesses=[20, 1.0, 1] # Guesses for (alphaG, Z, phi_g)
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
    vary_phi_cnm = param_vary[0] # Vary phi_cnm in K+09 fit?
    vary_Z = param_vary[1] # Vary metallicity in K+09 fit?
    vary_sigma_d = param_vary[2] # Vary sigma_d in K+09 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    guesses=[4, 1.0, 1.0] # Guesses for (phi_cnm, Z, sigma_d)
    h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

    krumholz_params = {}
    krumholz_params['param_vary'] = [vary_phi_cnm, vary_Z, vary_sigma_d]
    krumholz_params['error_method'] = error_method
    krumholz_params['alpha'] = alpha
    krumholz_params['guesses'] = guesses
    krumholz_params['h_sd_fit_range'] = h_sd_fit_range
    krumholz_params['results_filename'] = results_filename
    krumholz_params['parameters'] = ['phi_cnm', 'Z', 'sigma_d']
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

    # Get data products
    # ---------------------------------------------------------------------------
    nhi = data_products['nhi']
    nhi_error = data_products['nhi_error']
    av = data_products['av']
    av_error = data_products['av_error']
    dgr = mc_analysis['dgr']
    #dgr_error = mc_analysis['dgr_std']
    dgr_error = np.mean(np.abs(mc_analysis['dgr_error']))

    # Calculate N(H2) and error
    # ---------------------------------------------------------------------------
    ((nhi,
      nh2,
      nh,
      hi_sd,
      h2_sd,
      h_sd,
      rh2),
     (nhi_error,
      nh2_error,
      nh_error,
      hi_sd_error,
      h2_sd_error,
      h_sd_error,
      rh2_error)) = \
        lm_science.calc_coldens_products(nhi, av, dgr, nhi_error=nhi_error,
                                         av_error=av_error, dgr_error=dgr_error)

    # Write products to dictionary
    # ---------------------------------------------------------------------------
    data_products['nh2'] = nh2
    data_products['nh'] = nh
    data_products['h2_sd'] = h2_sd
    data_products['hi_sd'] = hi_sd
    data_products['h_sd'] = h_sd
    data_products['rh2'] = rh2
    data_products['nh2_error'] = np.abs(nh2_error)
    data_products['nh_error'] = np.abs(nh_error)
    data_products['h2_sd_error'] = np.abs(h2_sd_error)
    data_products['hi_sd_error'] = np.abs(hi_sd_error)
    data_products['h_sd_error'] = np.abs(h_sd_error)
    data_products['rh2_error'] = np.abs(rh2_error)

    # Median sim errors
    # --------------------------------------------------------------------------
    if mc_results['sim_images'] is not None:
        av_error = np.median(np.nanstd(mc_results['sim_images']['av_sim']))
        nhi_error = nhi_error + \
            np.median(np.nanstd(mc_results['sim_images']['nhi_sim']))
        dgr_error = np.mean(np.abs(mc_analysis['dgr_error']))

        ((nhi,
          nh2,
          hi_sd,
          h2_sd,
          h_sd,
          rh2),
         (nhi_error,
          nh2_error,
          hi_sd_error,
          h2_sd_error,
          h_sd_error,
          rh2_error)) = \
            lm_science.calc_coldens_products(nhi, av, dgr, nhi_error=nhi_error,
                                             av_error=av_error, dgr_error=dgr_error)


        data_products['h_sd_median_error'] = \
            scipy.stats.nanmedian(h_sd_error.ravel())
        data_products['hi_sd_median_error'] = \
                scipy.stats.nanmedian(hi_sd_error.ravel())
        data_products['rh2_median_error'] = \
                scipy.stats.nanmedian(rh2_error.ravel())

        if 1:
            rh2 = h2_sd / hi_sd
            rh2_error_before = rh2 * (h2_sd_error**2 / \
                                     h2_sd**2 + \
                                     h_sd_error**2 / h_sd**2)**0.5
            mask = np.isnan(rh2_error_before)
            mask[np.isnan(rh2_error)] = 1
            print 'rh2 error before - after:', rh2_error_before[~mask] - \
                rh2_error[~mask]

            data_products['nh2'] = nh2
            data_products['h2_sd'] = h2_sd
            data_products['hi_sd'] = hi_sd
            data_products['h_sd'] = h_sd
            data_products['rh2'] = rh2
            data_products['nh2_error'] = np.abs(nh2_error)
            data_products['h2_sd_error'] = np.abs(h2_sd_error)
            data_products['hi_sd_error'] = np.abs(hi_sd_error)
            data_products['h_sd_error'] = np.abs(h_sd_error)
            data_products['rh2_error'] = np.abs(rh2_error)
            print('hi_sd_error', scipy.stats.nanmedian(hi_sd_error.ravel()))
            print('h2_sd_error', scipy.stats.nanmedian(h2_sd_error.ravel()))


    else:
        data_products['h_sd_median_error'] = None
        data_products['hi_sd_median_error'] = None
        data_products['rh2_median_error'] = None

def calc_mc_analysis(mc_results, resid_results, data_products, model_kwargs):

    import mystats

    mc_analysis = {}
    core_analysis = {}

    # Calculate conf intervals on parameters
    # -------------------------------------------------------------------------
    # DGR
    dgrs = mc_results['av_model_results']['dgr']
    dgr, dgr_error = mystats.calc_cdf_error(dgrs,
                                            alpha=0.32)
    dgr_std = np.std(dgrs)

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
    radiation_type = \
        model_kwargs['model_kwargs']['sternberg_params']['radiation_type']
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
                                          radiation_type=radiation_type,
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
                   'dgr_std': dgr_std,
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
                                            results_dict['data_products'],
                                            results_dict['ss_model_kwargs'],
                                            )

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
    write_final_maps(results_dict, global_args)

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

    # if region 1 or 2 specified remove numbers
    cloud_region_name = np.copy(cloud_name)
    for string in ['1', '2']:
        cloud_name = cloud_name.replace(string, '')

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
    # ---------------------------------------------------------------------------
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

        print('\n\tVelocity range and error:')
        print('\t\t[{0:.2f}, {1:.2f}] +/- {2:.2f} km/s'.format(velocity_range[0],
                                                           velocity_range[1],
                                                           hi_range_error))
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

    # mask for erroneous pixel at RA,dec of 4h37m0s, 29d40m0s.
    nhi_image[nhi_image < -10] = np.nan
    #nhi_image_error[nhi_image_error < -10] = np.nan
    #nhi_image_background[nhi_image_background < 10] = np.nan

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
                                         vary_phi_g=global_args['vary_phi_g'],
                                         param_vary=global_args['param_vary'],
                                         )
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
                       calc_median_error=global_args['calculate_median_error'],
                       odr_fitting=global_args['odr_fitting'],
                       dust_cross_section_type=\
                               global_args['dust_cross_section_type'],
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
                    'ss_model_kwargs': global_args['ss_model_kwargs'],
                    }

    print('\n\tSaving results...')
    save_results(results_dict, global_args['results_filename'])
    #results_dict = load_results(global_args['results_filename'])

    return results_dict

def main():

    import itertools

    results = {}

    clouds = (
              'california',
              'perseus',
              'taurus',
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
               #'2',
               #'1',
               )

    subtract_comps = (#True,
                      False,
                      )

    radiation_type = ('beamed',
                      #'isotropic',
                      )

    rotate_cores = (
                    False,
                    #True,
                    )

    vary_phi_g = (
                    #True,
                    False,
                    )

    # fit for (phi_cnm, Z, sigma_d) in K+09 model
    #         (alphaG,  Z, phi_g)   in S+14 model
    param_vary = (
                  (True, True, False),
                  #(True, False, False),
                 )

    num_bootstraps = (
                      20000,
                      #20001,
                      #20002,
                      #100,
                      #10,
                      )

    dust_cross_section = (
                          'sternberg',
                          #'krumholz',
                          )

    elements = (clouds, data_types, recalculate_likelihoods, bin_image,
            init_vel_width, fixed_width, use_intercept, av_mask_threshold,
            regions, subtract_comps, use_background, hi_range_calc,
            radiation_type, rotate_cores, vary_phi_g, param_vary,
            num_bootstraps,
            dust_cross_section,)

    permutations = list(itertools.product(*elements))

    FILENAME_HI_DICT = '/d/bip3/ezbc/multicloud/data/' + \
                        'python_output/' + \
                        'core_properties/hi_properties.pickle'

    print('Number of permutations to run: ' + str(len(permutations)))

    #for cloud in clouds:
    for permutation in permutations:
        global_args = {
                'load': 1,
                'cloud_name':permutation[0],
                'odr_fitting': False,
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
                'num_resid_bootstraps': 100,
                'bootstrap_fit_residuals': False,
                'calculate_median_error': 0,
                'multiprocess': 1,
                'radiation_type': permutation[12],
                'rotate_cores': permutation[13],
                'vary_phi_g': permutation[14],
                'filename_hi_props': FILENAME_HI_DICT,
                'param_vary': permutation[15],
                'num_bootstraps': permutation[16],
                'dust_cross_section_type': permutation[17],
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

    #collect_results(results)
    plot_multicloud_results(results)

if __name__ == '__main__':
    main()



