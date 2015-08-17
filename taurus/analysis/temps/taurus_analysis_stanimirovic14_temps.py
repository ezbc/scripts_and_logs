
def gaussian(x, amp, sigma, x0):

    import numpy as np

    return amp * np.exp(-(x - x0)**2 / (2 * sigma)**2)

def read_data(filename, data_start=2):

    from astropy.io import ascii

    with open(filename, 'r') as f:
        table = ascii.read(f,
                           data_start=data_start,
                           delimiter='\t',
                           guess=True)

    return table

def main():

    import numpy as np
    import numpy
    from scipy.integrate import quad as integrate

    # Parameters in script
    # --------------------------------------------------------------------------
    clouds = ('taurus', 'perseus',)# 'california')
    dec_ranges = ((35, 20), (35, 24), (44, 35))
    ra_ranges = ((5.5, 3.9), (4, 3), (70, 40))
    tspin_threshold = 200.0

    # Data locations in script
    # --------------------------------------------------------------------------
    filedir = '/home/ezbc/research/data/taurus/python_output/'

    # Read the data
    source_table = read_data(filedir + 'taurus_stanimirovic14_sources.txt')
    param_table = read_data(filedir + 'taurus_stanimirovic14_temps.txt')

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
        source_dec = source[source_cols.index('Decl. (J2000)')]
        source_decs[i] = float(source_dec.split(':')[0])
        source_ra = source[source_cols.index('R.A. (J2000)')]
        source_ras[i] = float(source_ra.split(':')[0])
        source_names[i] = source[source_cols.index('Source')]

    # Find T_cnm temperatures for each cloud
    cloud_counts_list = []
    cloud_bins_list = []
    cloud_t_cnm_lists = []
    weights_list = []

    for i, cloud in enumerate(clouds):
    	dec_range = dec_ranges[i]
    	ra_range = ra_ranges[i]

        # Choose only sources within RA and Dec range
        indices = np.where((source_decs >= dec_range[1]) &\
                           (source_decs <= dec_range[0]) &\
                           (source_ras >= ra_range[1]) &\
                           (source_ras <= ra_range[0])
                           )

        cloud_decs = source_decs[indices]
        cloud_ras = source_ras[indices]
        cloud_sources = source_names[indices]

        # Get the spin temperatures of chosen sources
        t_cnm_list = []
        cloud_t_cnm_list = []
        source_list = []
        weights = []

        for j in xrange(len(param_data)):
            tspin = param_data[j][param_cols.index('Ts')]
            tau_error = param_data[j][param_cols.index('e_Ts')]
            if param_data[j][param_cols.index('Source')] in cloud_sources:
                if tspin < tspin_threshold and tau_error != 0.0:
                    cloud_t_cnm_list.append(tspin)

                    Tpeak = param_data[j][param_cols.index('TB')]
                    DelV = param_data[j][param_cols.index('DelV')]
                    VLSR = param_data[j][param_cols.index('VLSR')]

                    T_integ = integrate(gaussian,
                                        -1000,
                                        1000,
                                        args=(Tpeak,
                                              DelV/2.0,
                                              VLSR))

                    weights.append(T_integ[0])

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
        weights_list.append(weights)

        print cloud, np.median(cloud_t_cnm_list), 'Median CNM temp [K]'
        print cloud, np.average(cloud_t_cnm_list,
                                weights=weights), 'Intensity-weighted mean CNM temp [K]'



if __name__ == '__main__':
    main()
