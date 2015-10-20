#!/usr/bin/python

import pickle
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

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

        temp = core['dust_temp_avg']
        temp_error = core['dust_temp_error_avg']
        core['rad_field'], core['rad_field_error'] = \
                calc_radiation_field(temp,
                                     T_dust_error=temp_error)

        #rad_error = np.sort((calc_radiation_field(temp + temp_error),
        #                     calc_radiation_field(temp - temp_error)))
        #core['rad_field_error'] = rad_error
        for param in core['sternberg']:
            if type(core['sternberg'][param]) is tuple:
                core['sternberg'][param] = np.array(core['sternberg'][param])
        for param in core['krumholz']:
            if type(core['krumholz'][param]) is tuple:
                core['krumholz'][param] = np.array(core['krumholz'][param])

        # calculate n_H and error
        core['n_H'], core['n_H_error'] = \
                calc_n_H(I_UV=core['rad_field'],
                         alphaG=core['sternberg']['alphaG'],
                         phi_g=core['sternberg']['phi_g'],
                         #alphaG=10,
                         #phi_g=2,
                         I_UV_error=core['rad_field_error'],
                         alphaG_error=core['sternberg']['alphaG_error'],
                         phi_g_error=core['sternberg']['phi_g_error'],
                         )

        #if core['n_H'] <= 0:
        #    core['n_H'] = 10

        #core['n_H_error'] = [0.1 * core['n_H'], 0.1 * core['n_H']]

        core['T_H'] = calc_temperature(n_H=core['n_H'],
                                       pressure=3000.0) / 1000.0
        core['T_H_error'] = np.array([0.1 * core['T_H'], 0.1 * core['T_H']])
        T_H_error = np.sort((calc_temperature(core['n_H'] - \
                              core['n_H_error'][1]),
                             calc_temperature(core['n_H'] + \
                              core['n_H_error'][0])))
        core['T_H_error'] = T_H_error / 1000.0

        print core['n_H_error']

        core['alt_name'] = get_alt_name(core_name)

def add_dust_temps(core_dict):

    import myimage_analysis as myia
    import mygeometry as myg

    # Get the data
    # ------------
    dust_temp_dir = '/d/bip3/ezbc/multicloud/data/dust_temp/'
    temp_filename = dust_temp_dir + 'multicloud_dust_temp_5arcmin.fits'
    temp_error_filename = dust_temp_dir + \
            'multicloud_dust_temp_error_5arcmin.fits'
    temp_data, temp_header = fits.getdata(temp_filename, header=True)
    temp_error_data = fits.getdata(temp_error_filename)

    # Get the mask for each core
    # --------------------------
    # Create WCS object
    wcs_header = WCS(temp_header)
    for core_name in core_dict:
        vertices_wcs = core_dict[core_name]['region_vertices']

        # Format vertices to be 2 x N array
        #vertices_wcs = np.array((vertices_wcs[0], vertices_wcs[1]))

        # Make a galactic coords object and convert to Ra/dec
        coords_fk5 = SkyCoord(vertices_wcs[0] * u.deg,
                              vertices_wcs[1] * u.deg,
                              frame='fk5',
                              )
        # convert to pixel
        coords_pixel = np.array(coords_fk5.to_pixel(wcs_header))

        # write data to dataframe
        vertices_pix = np.array((coords_pixel[1], coords_pixel[0])).T

        core_dict[core_name]['region_vertices_pix'] = vertices_pix

        # Mask pixels outside of the region
        region_mask = np.logical_not(myg.get_polygon_mask(temp_data,
                                                          vertices_pix))
        core_dict[core_name]['region_mask'] = region_mask

        # Grab the temperatures
        core_dict[core_name]['dust_temps'] = temp_data[~region_mask]
        core_dict[core_name]['dust_temp_errors'] = temp_error_data[~region_mask]

        # Calculate average temp
        core_dict[core_name]['dust_temp_avg'] = \
            np.nanmean(temp_data[~region_mask])
        core_dict[core_name]['dust_temp_error_avg'] = \
            np.sqrt(np.nansum(temp_error_data[~region_mask]**2)) / \
                temp_error_data[~region_mask].size**0.5

    return core_dict

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
    cloud_old = ''
    for cloud_row, core_name in enumerate(core_name_list):
        core = core_dict[core_name]
        cloud = core['cloud']

        # add cloud name
        # -------------
        if cloud_old != cloud:
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
        param = core['dust_temp_avg']
        param_error = core['dust_temp_error_avg']
        param_error = [param_error, param_error]
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
                params_to_write = ['alphaG', 'phi_g', 'hi_transition']

            for i, param_name in enumerate(params_to_write):
                param = \
                    core[model][param_name]
                param_error = \
                    core[model][param_name + '_error']

                param_info = (param, param_error[1], param_error[0])

                print param_info

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
        if cloud_old != cloud and cloud != 'california':
            row_text = '\hline  \\\\[-0.2cm] \n' + row_text

        cloud_old = cloud

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

    # average dust temperatures over each core region
    add_dust_temps(core_dict)

    # Add model_analysis
    add_model_analysis(core_dict)

    # write the latex table
    write_model_params_table(core_dict)

if __name__ == '__main__':
    main()
