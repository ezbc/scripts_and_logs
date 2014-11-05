#!/usr/bin/python

def read_data(filename, data_start=3):

    from astropy.io import ascii

    with open(filename, 'r') as f:
        table = ascii.read(f,
                           data_start=data_start,
                           delimiter='\t',)

    return table

def plot_spin_temps(spin_temps_list, cloud_bins_list=10,
        global_bins=None, global_spin_temps=None, clouds=None,
        scales=('linear', 'linear'), filename=None, show=False):

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    for i in xrange(len(clouds)):
    	cloud = clouds[i]
    	cloud_bins = cloud_bins_list[i]
    	spin_temps = spin_temps_list[i]
    	color = ('c', 'b', 'g')[i]

        if type(cloud_bins) is int:
            ax.hist(spin_temps, color='k', bins=cloud_bins, alpha=0.2,
                    label=cloud, normed=True)
        else:
            ax.plot(cloud_bins, spin_temps, color=color, drawstyle='steps-mid',
                    label=cloud)

    if global_spin_temps is not None:
        if type(cloud_bins) is int:
            ax.hist(global_spin_temps, color='r', bins=global_bins, alpha=0.2,
                    label=r'Global', normed=True)
            ax.legend()
        else:
            ax.plot(global_bins, global_spin_temps, color='r',
                    drawstyle='steps-mid', linestyle='--', label=r'Global')
            ax.legend()

    ax.set_xlabel(r'T$_{s}$')
    ax.set_ylabel(r'Normalized Counts')

    ax.set_xscale(scales[0])
    ax.set_yscale(scales[1])

    if show:
        plt.show()
    if filename is not None:
        plt.savefig(filename)

    plt.close

def main():

    import numpy as np
    import numpy

    # Parameters in script
    # --------------------------------------------------------------------------
    clouds = ('taurus', 'perseus', 'california')
    dec_ranges = ((35, 20), (35, 24), (44, 35))
    ra_ranges = ((80, 60), (60, 40), (70, 40))
    tspin_threshold = 500.0

    # Data locations in script
    # --------------------------------------------------------------------------
    filedir = '/d/bip3/ezbc/multicloud/data/heiles03_data/'

    # Read the data
    source_table = read_data(filedir + 'heiles03_sources.tsv')
    param_table = read_data(filedir + 'heiles03_fit_params.tsv')

    # Extract headers and data
    source_cols = source_table.colnames
    source_data = np.asarray(source_table._data)
    param_cols = param_table.colnames
    param_data = np.asarray(param_table._data)

    # Extract source names and RA / Dec
    source_decs = np.empty(source_data.shape)
    source_ras = np.empty(source_data.shape)
    source_names = ['' for x in range(source_data.shape[0])]
    source_names = np.empty(source_data.shape, dtype=object)

    for i, source in enumerate(source_data):
        source_decs[i] = source[source_cols.index('_DE.icrs')]
        source_ras[i] = source[source_cols.index('_RA.icrs')]
        source_names[i] = source[source_cols.index('Name')]

    # Find T_cnm temperatures for each cloud
    cloud_counts_list = []
    cloud_bins_list = []
    cloud_t_cnm_lists = []

    for i, cloud in enumerate(clouds):
    	dec_range = dec_ranges[i]
    	ra_range = ra_ranges[i]
        # Choose only sources within RA and Dec range
        indices = np.where((source_decs > dec_range[1]) &\
                           (source_decs < dec_range[0]) &\
                           (source_ras > ra_range[1]) &\
                           (source_ras < ra_range[0])
                           )[0]

        cloud_decs = source_decs[indices]
        cloud_ras = source_ras[indices]
        cloud_sources = source_names[indices]

        # Get the spin temperatures of chosen sources
        t_cnm_list = []
        cloud_t_cnm_list = []
        source_list = []

        for i in xrange(len(param_data)):
            tspin = param_data[i][param_cols.index('Tspin')]
            tau_error = param_data[i][param_cols.index('e_tau')]
            if param_data[i][param_cols.index('Name')] in cloud_sources:
                if tspin < tspin_threshold and tau_error != 0.0:
                    cloud_t_cnm_list.append(tspin)
            if tspin < tspin_threshold:
                t_cnm_list.append(tspin)

        if cloud=='perseus':
        	bins=3
        else:
        	bins=8

        cloud_counts, cloud_bins = np.histogram(cloud_t_cnm_list, bins=bins)
        cloud_counts = cloud_counts / np.sum(cloud_counts, dtype=numpy.float)
        cloud_counts = np.append(cloud_counts, 0)

        cloud_counts_list.append(cloud_counts)
        cloud_bins_list.append(cloud_bins)
        cloud_t_cnm_lists.append(cloud_t_cnm_list)

        print cloud, len(cloud_t_cnm_list), ' ncomps'

    global_counts, global_bins = np.histogram(t_cnm_list, bins=20)
    global_counts = global_counts / np.sum(global_counts, dtype=numpy.float)
    global_counts = np.append(global_counts, 0)

    # Plot the spin temperatures
    plot_spin_temps(cloud_counts_list, cloud_bins_list=cloud_bins_list,
            global_spin_temps=global_counts, global_bins=global_bins,
            clouds=clouds,
            filename=\
                '/d/bip3/ezbc/multicloud/figures/heiles03_spin_temp_hist.png')


if __name__ == '__main__':
    main()

