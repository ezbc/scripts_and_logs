import matplotlib
matplotlib.use('Agg')

import warnings
warnings.filterwarnings('ignore')

'''
Main Script
'''

def main():

    import numpy as np
    from os import path
    import json
    from pandas import DataFrame
    import io

    av_data_type = 'planck'

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/'

    figure_types = ['png', 'pdf']

    cloud_list = (  'perseus', 'california', 'taurus', 'taurus1', 'taurus2',)
    data_dict = {}

    for cloud in cloud_list:
        if cloud == 'taurus1' or cloud == 'taurus2':
            prop_dir = '/d/bip3/ezbc/taurus/data/python_output/'
            cloud_name = 'taurus'
        else:
            prop_dir = '/d/bip3/ezbc/' + cloud + '/data/python_output/'
            cloud_name = cloud

        property_filename = cloud + '_global_properties_planck_scaled'

        with open(prop_dir + property_filename + '.txt', 'r') as f:
            global_props = json.load(f)

        data_dict[cloud] = {}
        data_dict[cloud]['width'] = global_props['hi_velocity_width']['value']
        data_dict[cloud]['width_error'] = \
                global_props['hi_velocity_width_error']['value']
        data_dict[cloud]['dgr'] = global_props['dust2gas_ratio']['value']
        data_dict[cloud]['dgr_error'] = \
                global_props['dust2gas_ratio_error']['value']
        data_dict[cloud]['intercept'] = global_props['intercept']['value']
        data_dict[cloud]['intercept_error'] = \
                global_props['intercept_error']['value']


        # Get summary files
        property_dir = \
            '/d/bip3/ezbc/' + cloud_name + \
                '/data/python_output/residual_parameter_results/'

        summary_filename = property_dir + \
                        property_filename.replace('_scaled','') + \
                        '_residscale1.5' + \
                        '.txt'

        with open(summary_filename, 'r') as f:
            results_summary[cloud] = json.load(f)

    #with open(table_dir + 'multicloud_parameters.tex', 'wb') as f:

    #f = io.open(table_dir + 'multicloud_parameters.tex', 'w',
    #            newline='\n')   # newline='' means don't convert \n

    f = open(table_dir + 'multicloud_parameters.tex', 'wb')

    for cloud in cloud_list:
        row_data = (cloud.capitalize(),
                    data_dict[cloud]['width'],
                    data_dict[cloud]['intercept'],
                    data_dict[cloud]['dgr'])

        width_text = \
            '{0:.1f}'.format(data_dict[cloud]['width']) + \
            '$^{{+{0:.1f}}}'.format(data_dict[cloud]['width_error'][0]) + \
            '_{{-{0:.1f}}}$'.format(data_dict[cloud]['width_error'][1])
        dgr_text = \
            '{0:.2f}'.format(data_dict[cloud]['dgr']) + \
            '$^{{+{0:.2f}}}'.format(data_dict[cloud]['dgr_error'][0]) + \
            '_{{-{0:.2f}}}$'.format(data_dict[cloud]['dgr_error'][1])
        intercept_text = \
            '{0:.1f}'.format(data_dict[cloud]['intercept']) + \
            '$^{{+{0:.1f}}}'.format(data_dict[cloud]['intercept_error'][0]) + \
            '_{{-{0:.1f}}}$'.format(data_dict[cloud]['intercept_error'][1])

        cloud = cloud.replace('1', ' 1')
        cloud = cloud.replace('2', ' 2')

        f.write(cloud.capitalize() + ' & ' + width_text + ' & ' + \
                #intercept_text + ' & ' + \
                dgr_text + \
                ' \\\\ \n')

    f.close()


if __name__ == '__main__':
    main()









