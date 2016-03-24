#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import mygeometry as myg
import pickle
import warnings
warnings.filterwarnings('ignore')

def load_cores():

    import pickle

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = table_dir + 'tables/multicloud_model_summary.pickle'

    with open(filename, 'rb') as f:
        core_dict = pickle.load(f)

    if 'n_H' not in core_dict[core_dict.keys()[0]]:
        raise ValueError('Need interpretation, run' + \
                         'cloud_write_parameter_table.py first')

    return core_dict

def main():

    import matplotlib
    import local_module_plotting as lm_plt

    DIR_FIGURE = '/d/bip3/ezbc/multicloud/figures/'
    FILENAME_FIGURE = DIR_FIGURE + 'models/himean_vs_radfield'

    # get core dict
    core_dict = load_cores()

    # plot
    file_types = ['png', 'pdf']
    for file_type in file_types:
        lm_plt.plot_radfield_vs_himean(core_dict,
                                          #limits=[0.6, 1.8],#, 1, 25],
                                          filename=FILENAME_FIGURE + '.' + \
                                                   file_type)



if __name__ == '__main__':
    main()

