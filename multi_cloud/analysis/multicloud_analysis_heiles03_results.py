#!/usr/bin/python

def read_data(filename, data_start=3):

    from astropy.io import ascii

    with open(filename, 'r') as f:
        table = ascii.read(f,
                           data_start=data_start,
                           delimiter='\t',)

    return table

def plot_spin_temps(spin_temps, taurus_bins=10,
        global_bins=None, global_spin_temps=None,
        scales=('linear', 'linear'), filename=None, show=False):

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    if type(taurus_bins) is int:
        ax.hist(spin_temps, color='k', bins=taurus_bins, alpha=0.2,
                label=r'Taurus', normed=True)
    else:
        ax.plot(taurus_bins, spin_temps, color='k', drawstyle='steps-mid',
                label=r'Taurus')


    if global_spin_temps is not None:
        if type(taurus_bins) is int:
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
    dec_range = (40, 20)
    ra_range = (70, 35)
    tspin_threshold = 150.0

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

    # Choose only sources within RA and Dec range
    indices = np.where((source_decs > dec_range[1]) &\
                       (source_decs < dec_range[0]) &\
                       (source_ras > ra_range[1]) &\
                       (source_ras < ra_range[0])
                       )[0]

    taurus_decs = source_decs[indices]
    taurus_ras = source_ras[indices]
    taurus_sources = source_names[indices]

    # Get the spin temperatures of chosen sources
    t_cnm_list = []
    taurus_t_cnm_list = []
    source_list = []

    for i in xrange(len(param_data)):
        tspin = param_data[i][param_cols.index('Tspin')]
        tau_error = param_data[i][param_cols.index('e_tau')]
        if param_data[i][param_cols.index('Name')] in taurus_sources:
            if tspin < tspin_threshold and tau_error != 0.0:
                taurus_t_cnm_list.append(tspin)
        if tspin < tspin_threshold:
            t_cnm_list.append(tspin)

    taurus_counts, taurus_bins = np.histogram(taurus_t_cnm_list, bins=8)
    global_counts, global_bins = np.histogram(t_cnm_list, bins=20)

    taurus_counts = taurus_counts / np.sum(taurus_counts, dtype=numpy.float)
    global_counts = global_counts / np.sum(global_counts, dtype=numpy.float)

    taurus_counts = np.append(taurus_counts, 0)
    global_counts = np.append(global_counts, 0)

    # Plot the spin temperatures
    plot_spin_temps(taurus_t_cnm_list, taurus_bins=8,
            global_spin_temps=t_cnm_list, global_bins=20,
            filename=\
                '/d/bip3/ezbc/multicloud/figures/heiles03_spin_temp_hist.png')
    plot_spin_temps(taurus_counts, taurus_bins=taurus_bins,
            global_spin_temps=global_counts, global_bins=global_bins,
            filename=\
                '/d/bip3/ezbc/multicloud/figures/heiles03_spin_temp_hist.png')


if __name__ == '__main__':
    main()

