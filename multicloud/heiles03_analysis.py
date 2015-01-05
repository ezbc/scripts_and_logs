#!/usr/bin/python

def read_data(filename, data_start=3):

    from astropy.io import ascii

    with open(filename, 'r') as f:
        table = ascii.read(f,
                           data_start=data_start,
                           delimiter='\t',)

    return table

def plot_spin_temps(spin_temps, bins=10, scales=('linear', 'linear'),
        filename=None, show=False):

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    ax.hist(spin_temps, color='k', bins=bins, alpha=0.2)

    ax.set_xlabel(r'T$_{s}$')
    ax.set_ylabel(r'Counts')

    ax.set_xscale(scales[0])
    ax.set_yscale(scales[1])

    if show:
        plt.show()
    if filename is not None:
    	plt.savefig(filename)

    plt.close

def main():

    import numpy as np

    filedir = '/d/bip3/ezbc/multicloud/data/'
    filedir = './'

    source_table = read_data(filedir + 'heiles03_sources.tsv')
    param_table = read_data(filedir + 'heiles03_fit_params.tsv')

    source_cols = source_table.colnames
    source_data = np.asarray(source_table._data)
    param_cols = param_table.colnames
    param_data = np.asarray(param_table._data)

    print param_cols
    print param_data[:10]

    dec_range = (45, 15)
    ra_range = (75, 30)

    source_decs = np.empty(source_data.shape)
    source_ras = np.empty(source_data.shape)
    source_names = ['' for x in range(source_data.shape[0])]
    source_names = np.empty(source_data.shape, dtype=object)

    for i, source in enumerate(source_data):
        source_decs[i] = source[source_cols.index('_DE.icrs')]
        source_ras[i] = source[source_cols.index('_RA.icrs')]
        source_names[i] = source[source_cols.index('Name')]

    indices = np.where((source_decs > dec_range[1]) &\
                       (source_decs < dec_range[0]) &\
                       (source_ras > ra_range[1]) &\
                       (source_ras < ra_range[0])
                       )[0]

    taurus_decs = source_decs[indices]
    taurus_ras = source_ras[indices]
    taurus_sources = source_names[indices]

    t_cnm_list = []
    source_list = []

    #print taurus_sources

    for i in xrange(len(param_data)):
        if param_data[i][param_cols.index('Name')] in taurus_sources:
            #print param_data[i][param_cols.index('Name')]
            tspin = param_data[i][param_cols.index('Tspin')]
            if tspin < 150.0:
                t_cnm_list.append(tspin)

    plot_spin_temps(t_cnm_list, bins=20,
            filename=\
                './heiles03_spin_temp_hist.png')


if __name__ == '__main__':
    main()

