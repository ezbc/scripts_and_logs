#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import myimage_analysis as myia
from multiprocessing.queues import Queue
import mygeometry as myg
import scipy

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
        filename):

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

            # append model params and errors to row
            for model in ('krumholz', 'sternberg'):
                if model == 'krumholz':
                    params_to_write = ['phi_cnm', 'Z', 'phi_mol',
                                       'hi_transition']
                else:
                    params_to_write = ['alphaG', 'Z', 'phi_g',
                                       'hi_transition']
                core_new[model] = {}
                for i, param_name in enumerate(params_to_write):
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

def write_param_csv(mc_analysis_dict, core_list, cloud_name_list, filename):

    import pandas as pd

    d = {}

    #print_dict_keys(mc_analysis_dict)
    params_to_write = ['phi_cnm', 'alphaG', 'hi_transition']

    d['cloud'] = []
    d['core'] = []
    d['ra'] = []
    d['dec'] = []
    d['region_vertices'] = []
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

def write_hi_vel_range_table(names_list, hi_range_kwargs_list, filename):

    # Open file to be appended
    f = open(filename, 'wb')

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

def add_row_element(row_text, element, text_format='{0:s}'):

    if type(element) is list or type(element) is tuple:
        return row_text + ' & ' + text_format.format(*element)
    else:
        return row_text + ' & ' + text_format.format(element)


