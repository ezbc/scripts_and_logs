#!/usr/bin/python

import numpy as np

def gaussian(x, amp, sigma, x0):

    import numpy as np

    return amp * np.exp(-(x - x0)**2 / (2 * sigma)**2)

def read_data(filename, data_start=2, delimiter='\t'):

    from astropy.io import ascii

    with open(filename, 'r') as f:
        table = ascii.read(f,
                           data_start=data_start,
                           delimiter=delimiter,
                           guess=True)

    return table

def calc_avg_Ts(TB_dict, tau_dict, x,):

    from scipy.integrate import simps

    # initialize emtpy array
    avg_Ts_array = np.empty(len(TB_dict))

    # compute average temperature for each source
    for i, source in enumerate(TB_dict):
        TB = TB_dict[source]
        tau = tau_dict[source]

        mask_zero_val = (TB == 0) | (tau == 0)
        TB_nonzero = TB[~mask_zero_val]
        tau_nonzero = tau[~mask_zero_val]
        x_nonzero = x[~mask_zero_val]

        avg_Ts_array[i] = simps(TB_nonzero / (1.0 - np.exp(-tau_nonzero)),
                                x_nonzero)
        avg_Ts_array[i] = simps(TB_nonzero / tau_nonzero,
                                x_nonzero)

    return avg_Ts_array

def main():

    import numpy as np
    import numpy
    from scipy.integrate import quad as integrate
    import pickle

    # Parameters in script
    # --------------------------------------------------------------------------
    clouds = ('taurus', 'perseus', 'california')
    dec_ranges = ((37, 19), (37, 19), (37, 35))
    ra_ranges = ((5.5, 3.8), (3.8, 2.9), (5.5, 3.5))
    tspin_threshold = 500.0

    # Data locations in script
    # --------------------------------------------------------------------------
    filedir = '/d/bip3/ezbc/multicloud/data/cnm_data/stanimirovic14/'

    # Read the data
    source_table = read_data(filedir + 'stanimirovic14_sources.txt')
    param_table = read_data(filedir + 'stanimirovic14_temps.txt')
    avg_temp_table = read_data(filedir + 'stanimirovic14_avg_temps.txt',
                               data_start=0,
                               delimiter=' ')

    # Extract headers and data
    source_cols = source_table.colnames
    source_data = np.asarray(source_table._data)
    param_cols = param_table.colnames
    param_data = np.asarray(param_table._data)


    print('cnm temp columns', param_cols)
    print('avg temp columns', avg_temp_table.colnames)
    #avg_temp_cols = np.asarray(avg_temp_table.colnames)
    #avg_temp_data = np.asarray(avg_temp_table._data)

    avg_temp_sources = []
    avg_temp_data = []
    for row in avg_temp_table._data:
        avg_temp_sources.append(row[0])
        avg_temp_data.append(row[1])
    avg_temp_data = np.asarray(avg_temp_data)
    avg_temp_sources = np.asarray(avg_temp_sources)

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
    temp_dict = {}


    sources_taurus = ['4C+27.14',
                      '3C133',
                      '3C132',
                      '4C+25.14',
                      '3C108',
                      'B20400+25',]

    for i, cloud in enumerate(clouds):
        temp_dict[cloud] = {}

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
        t_int_list = []
        tk_list = []
        tk_error_list = []
        t_cnm_error_list = []
        avg_tspin_list = []
        avg_tspin_error_list = []
        cloud_t_cnm_list = []
        source_list = []
        weights = []
        TB_dict = {}
        tau_dict = {}
        velocities = np.linspace(-1000, 1000, 1000)

        for j in xrange(len(param_data)):
            tspin = param_data[j][param_cols.index('Ts')]
            TB = param_data[j][param_cols.index('TB')]
            tau = param_data[j][param_cols.index('tau')]
            tspin_error = param_data[j][param_cols.index('e_Ts')]
            source = param_data[j][param_cols.index('Source')]

            Tpeak = param_data[j][param_cols.index('TB')]
            Tk = param_data[j][param_cols.index('Tkmax')]
            Tk_error = Tk * 0.1
            DelV = param_data[j][param_cols.index('DelV')]
            VLSR = param_data[j][param_cols.index('VLSR')]

            if source in cloud_sources:
                if source not in source_list:
                    avg_temp_added = False
                else:
                    avg_temp_added = True
                source_list.append(source)

                #print source

                # add the source avg Tspin
                if 1:
                #if not avg_temp_added:
                    print ''
                    print source
                    print 'TB:', TB
                    print 'tau:', tau
                    print 'T_kinetic:', Tk
                    print 'sigma:', DelV/2.0
                    print 'vlsr:', VLSR
                    avg_temp = \
                        avg_temp_data[avg_temp_sources == source][0]

                    # for error see section 3.1 of Kim et al. (2014) who
                    # describe the harmonic mean temperature is within 10% of
                    # the observed optical depth weighted spin temperature. This
                    # may be due to multiple clouds with varying optical depths
                    # along the LOS.
                    avg_temp_error = 0.1 * avg_temp

                    avg_tspin_list.append(avg_temp)
                    avg_tspin_error_list.append(avg_temp_error)

                #avg_temp = np.mean(TB) / np.mean(tau)
                else:

                    # Calculate brightness temp and tau spectra
                    TB_spectrum = gaussian(velocities, TB, DelV/2.0, VLSR)
                    tau_spectrum = gaussian(velocities, tau, DelV/2.0, VLSR)

                    # add them to the source list
                    if source not in TB_dict:
                        TB_dict[source] = TB_spectrum
                        tau_dict[source] = tau_spectrum
                    else:
                        TB_dict[source] = \
                            TB_dict[source] + TB_spectrum
                        tau_dict[source] = \
                            tau_dict[source] + tau_spectrum

                #avg_tspin_list.append(avg_temp)

                if tspin < tspin_threshold and tspin_error != 0.0:
                    cloud_t_cnm_list.append(tspin)

                    T_integ = integrate(gaussian,
                                        -1000,
                                        1000,
                                        args=(Tpeak,
                                              DelV/2.0,
                                              VLSR))

                    #weights.append(T_integ[0])

                    if tspin < tspin_threshold and VLSR < 15 and VLSR > -5:
                        t_cnm_list.append(tspin)
                        t_cnm_error_list.append(tspin_error)
                        t_int_list.append(Tk)

                    if source in sources_taurus:
                        print source, tspin

                if VLSR < 15 and VLSR > -5:
                    tk_list.append(Tk)
                    tk_error_list.append(Tk_error)

        # convert TB and tau to average spin temperatue
        #avg_tspin_list = calc_avg_Ts(TB_dict, tau_dict, velocities)

        #print 'avg t spins:', avg_tspin_list

        temp_dict[cloud]['Ts_list'] = t_cnm_list
        temp_dict[cloud]['Ts_error_list'] = t_cnm_error_list
        #avg_tspin_list = []
        #for source in TB_dict:
        #    avg_tspin_list.append(TB_dict[source])
        temp_dict[cloud]['Ts_avg_list'] = avg_tspin_list
        temp_dict[cloud]['Ts_avg_error_list'] = avg_tspin_error_list
        temp_dict[cloud]['Tk_list'] = tk_list
        temp_dict[cloud]['Tk_error_list'] = tk_error_list

        if 0:
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

            avg = np.average(cloud_t_cnm_list,
                             weights=weights)
            med = np.median(cloud_t_cnm_list)

            print(cloud.capitalize())
            print('\tMedian CNM temp = ' + \
                  '{0:.2f} [K]'.format(med))
            print('\tIntensity-weighted mean CNM temp = ' + \
                  '{0:.2f} [K]\n'.format(avg))

    # Save data
    with open('/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
              'stanimirovic14_temps.pickle', 'wb') as f:
        pickle.dump(temp_dict, f)

    filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
            'stanimirovic14_temps.npy'
    np.savetxt(filename, np.array(t_cnm_list))
    filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
            'stanimirovic14_int_temps.npy'
    np.savetxt(filename, np.array(t_int_list))
    filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
            'stanimirovic14_temp_errors.npy'
    np.savetxt(filename, np.array(t_cnm_error_list))

if __name__ == '__main__':
    main()
