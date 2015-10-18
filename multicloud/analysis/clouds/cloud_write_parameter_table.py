#!/usr/bin/python

import pickle
import numpy as np

def load_cores():

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = table_dir + 'tables/multicloud_model_summary.pickle'

    with open(filename, 'rb') as f:
        core_dict = pickle.load(f)

    return core_dict

def add_model_analysis(core_dict):

    from myscience import calc_radiation_field, calc_temperature
    from myscience.sternberg14 import calc_n_H

    ''' need the following parameters

    Cloud Name} &
        Core Name} &
        Alternate Name} &
        $T_{\rm dust}$} &
        $U$} &
        \aG} &
        \phicnm} &
        \hisd\ Threshold} &
        $n$} &
        $T_{\rm CNM}$} \\
    '''

    for core_name in core_dict:
        core = core_dict[core_name]

        core['rad_field'] = calc_radiation_field(core['temp'])
        core['rad_field_error'] = [0.1 * core['rad_field'],
                                   0.1 * core['rad_field']]

        core['n_H'] = calc_n_H(I_UV=core['rad_field'],
                               alphaG=core['sternberg']['alphaG'],
                               phi_g=core['sternberg']['phi_g'])
        core['n_H_error'] = [0.1 * core['n_H'], 0.1 * core['n_H']]

        core['T_H'] = calc_temperature(n_H=core['n_H'],
                                       pressure=3800.0)
        core['T_H_error'] = [0.1 * core['T_H'], 0.1 * core['T_H']]

        core['alt_name'] = get_alt_name(core_name)

def get_alt_name(core_name):

    if core_name == 'G164.99-8.60' or core_name == 'G165.71-9.15':
        alt_name = 'L1482'
    elif core_name == 'G166.91-7.76':
        alt_name = 'L1483'
    elif core_name == 'G158.39-20.72':
        alt_name = 'NGC1333'
    elif core_name == 'G158.89-21.60':
        alt_name = 'L1455'
    elif core_name == 'G159.19-20.11':
        alt_name = 'B1'
    elif core_name == 'G159.80-18.49':
        alt_name = 'B3'
    elif core_name == 'G160.46-17.99':
        alt_name = 'IC348'
    elif core_name == 'G160.49-16.81':
        alt_name = 'B5'
    elif core_name == 'G168.10-16.38':
        alt_name = 'L1495'
    elif core_name == 'G169.32-16.17':
        alt_name = 'B213'
    elif core_name == 'G171.00-15.80':
        alt_name = 'B217'
    elif core_name == 'G172.12-16.94':
        alt_name = 'B215'
    elif core_name == 'G172.93-16.73':
        alt_name = 'L1524'
    elif core_name == 'G174.40-13.45':
        alt_name = 'B220'
    elif core_name == 'G174.70-15.47':
        alt_name = 'B18'
    else:
        alt_name = ''

    return alt_name

def get_core_names():


    core_names = [
        'G164.18-8.84',
        'G164.26-8.39',
        'G164.65-8.12',
        'G164.70-7.63',
        'G164.99-8.60',
        'G165.36-7.51',
        'G165.71-9.15',
        'G166.91-7.76',
        'G168.12-6.42',
        'G168.54-6.22',
        'G158.26-21.81',
        'G158.39-20.72',
        'G158.89-21.60',
        'G159.17-21.09',
        'G159.19-20.11',
        'G159.80-18.49',
        'G160.14-19.08',
        'G160.46-17.99',
        'G160.49-16.81',
        'G160.53-19.73',
        'G168.10-16.38',
        'G169.32-16.17',
        'G171.00-15.80',
        'G171.14-17.57',
        'G171.49-14.91',
        'G172.12-16.94',
        'G172.93-16.73',
        'G174.05-15.82',
        'G174.40-13.45',
        'G174.70-15.47',
        ]

    return core_names

def write_model_params_table(core_dict):

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = table_dir + 'tables/modelparams_table.tex'

    # get ordered core names
    core_name_list = get_core_names()

    # Open file to be appended
    f = open(filename, 'wb')

    text_param_format ='{0:.1f}$^{{+{1:.1f}}}_{{-{2:.1f}}}$'

    #print_dict_keys(mc_analysis_dict)
    params_to_write = ['phi_cnm', 'alphaG']

    # Collect parameter names for each model for each core
    cloud_old = 'california'
    add_cloud = True
    for cloud_row, core_name in enumerate(core_name_list):
        core = core_dict[core_name]
        cloud = core['cloud']

        # add cloud name
        # -------------
        if add_cloud:
            row_text = cloud.capitalize()
        else:
            row_text = ''

        # add core name
        # -------------
        row_text = add_row_element(row_text,
                                   core_name)

        # add alternate core name
        # -------------
        row_text = add_row_element(row_text,
                                   core['alt_name'])

        # add dust temp
        # -------------
        param = core['temp']
        param_error = core['temp_error']
        param_info = (param, param_error[1], param_error[0])
        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format)

        # add rad field
        # -------------
        param = core['rad_field']
        param_error = core['rad_field_error']
        param_info = (param, param_error[1], param_error[0])
        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format)

        # append model params and errors
        # ------------------------------
        for model in ( 'krumholz', 'sternberg',):
            if model == 'krumholz':
                params_to_write = ['phi_cnm',]
            else:
                params_to_write = ['alphaG', 'hi_transition']

            for i, param_name in enumerate(params_to_write):
                param = \
                    core[model][param_name]
                param_error = \
                    core[model][param_name + '_error']

                param_info = (param, param_error[1], param_error[0])

                row_text = \
                    add_row_element(row_text,
                                    param_info,
                                    text_format=text_param_format)

        # add H vol dens
        # --------------
        param = core['n_H']
        param_error = core['n_H_error']
        param_info = (param, param_error[1], param_error[0])
        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format)

        # add HI temp
        # -----------
        param = core['T_H']
        param_error = core['T_H_error']
        param_info = (param, param_error[1], param_error[0])
        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format)

        # Finish row
        row_text += ' \\\\[0.1cm] \n'
        if cloud_old != cloud \
            and cloud != 'taurus':
            row_text += '\hline  \\\\[-0.2cm] \n'
            cloud_old = cloud
        elif cloud_old != cloud and \
                cloud == 'taurus':
            row_text.replace(r'\\[0.1cm] \n', '')
        else:
            add_cloud = False

        f.write(row_text)

    f.close()

def add_row_element(row_text, element, text_format='{0:s}'):

    if type(element) is list or type(element) is tuple:
        return row_text + ' & ' + text_format.format(*element)
    else:
        return row_text + ' & ' + text_format.format(element)

def main():

    # load core summary file
    core_dict = load_cores()

    add_model_analysis(core_dict)

    # write the latex table
    write_model_params_table(core_dict)

if __name__ == '__main__':
    main()
