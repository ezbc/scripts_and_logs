#!/usr/bin/python

import pickle
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import myimage_analysis as myia
import mygeometry as myg
import mystats

def plot_dust_histogram(dust_temps, limits=None, filename=None):

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    font_scale = 9
    params = {
              'figure.figsize': (3.5, 6),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    myplt.set_color_cycle(num_colors=1, cmap_limits=[0, 0.6])

    # Create figure instance
    fig = plt.figure()

    axes = AxesGrid(fig, (1,1,1),
                    nrows_ncols=(3, 1),
                    ngrids=3,
                    axes_pad=0.1,
                    aspect=False,
                    label_mode='L',
                    share_all=True)


    def hist(data, normalize=True):
        # Derive the histograms
        n_bins = 300
        #bin_edges = np.logspace(-3, 3, num=n_bins, base=base)
        bin_edges = np.linspace(np.min(data), np.max(data),
                                num=n_bins)
        n = np.zeros(n_bins - 1)
        for i in xrange(n_bins - 1):
            bin_count = len(data[(data > bin_edges[i]) & \
                                 (data < bin_edges[i + 1])])
            n[i] = bin_count

        bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

        n = np.append(n, 0)
        bin_centers = np.append(bin_centers,
                bin_centers[-1] + (bin_centers[-1] - bin_centers[-2]))

        if limits is not None:
            bin_size = (bin_centers[-1] - bin_centers[-2])
            while bin_centers[-1] < limits[0][1]:
                print 'end bin:', bin_centers[-1]
                n = np.append(n, 0)
                bin_centers = np.append(bin_centers,
                                        bin_centers[-1] + bin_size)
            while bin_centers[0] > limits[0][0]:
                print ':', bin_centers[0]
                n = np.append(0, n)
                bin_centers = np.append(bin_centers[0] - bin_size,
                                        bin_centers,
                                        )

        # Normalize the bins
        if normalize:
            n /= np.max(n)

        return n, bin_centers

    for i, cloud_name in enumerate(['california', 'perseus', 'taurus']):

        cloud_temps = dust_temps[cloud_name]['dust_temps']

        n, bin_centers = hist(cloud_temps)

        ax = axes[i]
        ax.locator_params(nbins = 6)

        ax.errorbar(
            bin_centers,
            n,
            #yerr = n**0.5,
            marker = '',
            #label=r'',
            linewidth=1.5,
            color='k',
            drawstyle = 'steps-mid'
        )

        ax.annotate(cloud_name.capitalize(),
                    xytext=(0.96, 0.9),
                    xy=(0.96, 0.9),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=10,
                    color='k',
                    bbox=dict(boxstyle='square',
                              facecolor='w',
                              alpha=1),
                    horizontalalignment='right',
                    verticalalignment='top',
                    )
        # legend!
        if i == 2:
            ax.legend(loc='upper left')

        # plot limits
        if limits is not None:
            ax.set_xlim(limits[0][0],limits[0][1])
            ax.set_ylim(limits[1][0],limits[1][1])

        ax.set_xlabel(r'$T_{\rm dust}$ [K]')
        ax.set_ylabel('PDF')

    if filename is not None:
        plt.draw()
        plt.savefig(filename,
                    bbox_inches='tight', dpi=100)

def load_cores():

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = table_dir + 'tables/multicloud_model_summary.pickle'

    with open(filename, 'rb') as f:
        core_dict = pickle.load(f)

    return core_dict

def add_model_params(core_dict):

    from myscience import calc_radiation_field, calc_temperature
    from myscience.sternberg14 import calc_n_H
    import myscience.krumholz09 as myk09

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

    clouds_printed = []
    for core_name in core_dict:
        core = core_dict[core_name]

        temp = core['dust_temp_median']
        temp_error = core['dust_temp_median_error']
        core['rad_field'], core['rad_field_error'] = \
                calc_radiation_field(temp,
                                     T_dust_error=temp_error)

        cloud = core_dict[core_name]['cloud']
        if 0:
            if cloud not in clouds_printed:
                print ''
                print cloud + ' dust temp and error:'
                print core['dust_temp_median']
                print core['dust_temp_median_error']
                print 'rad field and error:'
                print core['rad_field']
                print core['rad_field_error']
                clouds_printed.append(cloud)

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
                calc_n_H(I_UV=core['rad_field_draine_median'],
                         alphaG=core['sternberg']['alphaG'],
                         phi_g=core['sternberg']['phi_g'],
                         Z_g=core['sternberg']['Z'],
                         I_UV_error=core['rad_field_draine_median_error'],
                         alphaG_error=core['sternberg']['alphaG_error'],
                         Z_g_error=core['sternberg']['Z_error'],
                         )

        #if core['n_H'] <= 0:
        #    core['n_H'] = 10

        #core['n_H_error'] = [0.1 * core['n_H'], 0.1 * core['n_H']]

        if 0:
            core['T_H'] = calc_temperature(n_H=core['n_H'],
                                           pressure=3000.0) / 1000.0
            core['T_H_error'] = np.array([0.1 * core['T_H'], 0.1 * core['T_H']])
            T_H_error = np.sort((calc_temperature(core['n_H'] - \
                                  core['n_H_error'][1]),
                                 calc_temperature(core['n_H'] + \
                                  core['n_H_error'][0])))
            core['T_H_error'] = T_H_error / 1000.0
            if abs(core['T_H_error'][0]) > core['T_H']:
                core['T_H_error'][0] = core['T_H']



        pressure = 3700.0
        pressure_error = (1200.0, 1200.0)
        core['T_H'], core['T_H_error'] = \
                calc_temperature(n_H=core['n_H'],
                                 pressure=pressure,
                                 pressure_error=pressure_error,
                                 n_H_error=core['n_H_error'])

        core['T_H'] /= 1000.0
        core['T_H_error'] = np.array(core['T_H_error']) / 1000.0

        # calculate krumholz params
        krumholz_pressure_calc = True
        if krumholz_pressure_calc:
            # get the minimum CNM density, calculate n_CNM with phi_CNM, then
            # calculate temperature given galactic pressure between WNM and CNM
            n_min, n_min_error = \
                myk09.calc_n_min(G_0=core['rad_field_habing_median'],
                               G_0_error=core['rad_field_habing_median_error'],
                               Z=1.0,
                               calc_error=True,
                               )
            phi_cnm = core['krumholz']['phi_cnm']
            phi_cnm_error = core['krumholz']['phi_cnm_error']
            core['n_cnm'] = n_min * phi_cnm
            core['n_cnm_error'] = ((phi_cnm * n_min_error)**2 + \
                                   (n_min * phi_cnm_error)**2)**0.5

            #print core['n_cnm_error']

            core['T_cnm'], core['T_cnm_error'] = \
                calc_temperature(n_H=core['n_cnm'],
                                 pressure=pressure,
                                 pressure_error=pressure_error,
                                 n_H_error=core['n_cnm_error'])

        else:
            core['T_cnm'], core['T_cnm_error'] = \
                    myk09.calc_T_cnm(core['krumholz']['phi_cnm'],
                        phi_cnm_error=core['krumholz']['phi_cnm_error'],
                        calc_error=True,
                        )

            core['n_cnm'], core['n_cnm_error'] = \
                    myk09.calc_n_cnm(G_0=(core['rad_field'] / 1.7),
                             G_0_error=(core['rad_field_error'] / 1.7),
                             T_cnm=core['T_cnm'],
                             T_cnm_error=core['T_cnm_error'],
                             calc_error=True,
                             )

        core['alt_name'] = get_alt_name(core_name)

def write_cloud_temps_table(cloud_temps):

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = table_dir + 'tables/cloud_params_table.tex'

    text_param_format ='{0:.2f}\,$\pm$\,{1:.2f}'

    params = [
        'dust_temp_median',
        'dust_beta_median',
        'rad_field_mathis_median',
        'rad_field_draine_median',
        'rad_field_habing_median',
        ]

    # latex names
    param_names_pretty = [r'T$_{\rm dust}$\,[K]',
                          r'$\beta$',
                          r'$U_{M83}\,[U_{M83,0}$]',
                          r'$I_{\rm UV}\,[I_{D,0}$]',
                          r'$G^\prime_{0}\,[G_{0}$]',
                          ]

    f = open(filename, 'wb')

    for i, param_name in enumerate(params):
        row_text = ''
        param_name_pretty = param_names_pretty[i]
        row_text += param_name_pretty

        for cloud in ('california', 'perseus', 'taurus'):

            # get the parameter
            param = cloud_temps[cloud][param_name]
            param_error = cloud_temps[cloud][param_name + '_error']

            row_text = add_row_element(row_text,
                                       (param, param_error),
                                       text_format=text_param_format,
                                       )

        # finish row
        row_text += ' \\\\[0.1cm] \n'

        f.write(row_text)

    f.close()

def run_mc_simulations(core_dict, wcs_header, temp_data, temp_error_data,
        beta_data, beta_error_data):

    from myscience import calc_radiation_field

    cloud_temps = {}
    for core_name in core_dict:
        # load cloud regions
        core_dict = add_cloud_region(core_dict)
        vertices_wcs = core_dict[core_name]['cloud_region_vertices'].T

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

        core_dict[core_name]['cloud_region_vertices_pix'] = vertices_pix

        # Mask pixels outside of the region
        region_mask = np.logical_not(myg.get_polygon_mask(temp_data,
                                                          vertices_pix))
        core_dict[core_name]['cloud_region_mask'] = region_mask

        # Grab the temperatures
        core_dict[core_name]['dust_temps'] = temp_data[~region_mask]
        core_dict[core_name]['dust_temp_errors'] = \
            temp_error_data[~region_mask]

        # adjust vertices to get errors on mean T_dust
        cloud = core_dict[core_name]['cloud']
        N_mc = 10
        temp_mc = np.empty(N_mc)
        temp_error_mc = np.empty(N_mc)
        beta_mc = np.empty(N_mc)
        beta_error_mc = np.empty(N_mc)
        rad_mc = np.empty(N_mc)
        rad_error_mc = np.empty(N_mc)
        if cloud not in cloud_temps:
            for j in xrange(N_mc):
                if j != 0:
                    new_vertices_wcs = vertices_wcs + \
                                       np.random.normal(scale=1.0,
                                                        size=vertices_wcs.shape)
                else:
                    new_vertices_wcs = vertices_wcs

                # Make a galactic coords object and convert to Ra/dec
                coords_fk5 = SkyCoord(new_vertices_wcs[0] * u.deg,
                                      new_vertices_wcs[1] * u.deg,
                                      frame='fk5',
                                      )

                # convert to pixel
                coords_pixel = np.array(coords_fk5.to_pixel(wcs_header))

                # write data to dataframe
                vertices_pix = np.array((coords_pixel[1],
                                         coords_pixel[0])).T

                # Mask pixels outside of the region
                region_mask = \
                        np.logical_not(myg.get_polygon_mask(temp_data,
                                                            vertices_pix))

                # Get the region's temperature
                if j == 0:
                    temps = temp_data[~region_mask]
                    betas = beta_data[~region_mask]
                    rads = calc_radiation_field(temps,
                                             beta=betas,
                                             )

                # simulate new observation of temperature and beta
                temp_sim = temp_data + np.random.normal(0,
                                                        scale=temp_error_data,)
                beta_sim = beta_data + np.random.normal(0,
                                                        scale=beta_error_data,)

                # Calculate the radiation field
                # -----------------------------
                rad_field = \
                    calc_radiation_field(temp_sim,
                                         beta=beta_sim,
                                         )

                # Grab the median values of temp, beta, and rad field
                temp_mc[j] = np.median(temp_sim[~region_mask])
                beta_mc[j] = np.median(beta_sim[~region_mask])
                rad_mc[j] = np.median(rad_field[~region_mask])

            # Calculate average temp
            #core_dict[core_name]['dust_temp_median'] = \
            #    np.nanmean(temp_data[~region_mask])
            #core_dict[core_name]['dust_temp_median_error'] = \
            #    np.sqrt(np.nansum(temp_error_data[~region_mask]**2)) / \
            #        temp_error_data[~region_mask].size
            dust_temp_median, mc_error = mystats.calc_cdf_error(temp_mc)
            dust_temp_median_error = np.mean(mc_error)
            dust_beta_median, mc_error = mystats.calc_cdf_error(beta_mc)
            dust_beta_median_error = np.mean(mc_error)
            rad_field_draine_median, mc_error = mystats.calc_cdf_error(rad_mc)
            rad_field_draine_median_error = np.mean(mc_error)

            # calculate habing field from draine:
            rad_field_habing_median = rad_field_draine_median / 1.71
            rad_field_habing_median_error = rad_field_draine_median_error / 1.71
            rad_field_mathis_median = rad_field_draine_median / 1.48
            rad_field_mathis_median_error = rad_field_draine_median_error / 1.48


            cloud_temps[cloud] = \
                    {
                     'dust_temp_median': dust_temp_median,
                     'dust_temp_median_error': dust_temp_median_error,
                     'dust_temps': temps,
                     'dust_beta_median': dust_beta_median,
                     'dust_beta_median_error': dust_beta_median_error,
                     'dust_betas': betas,
                     'rad_field_draine_median': rad_field_draine_median,
                     'rad_field_draine_median_error': \
                        rad_field_draine_median_error,
                     'rad_field_habing_median': rad_field_habing_median,
                     'rad_field_habing_median_error': \
                        rad_field_habing_median_error,
                     'rad_field_mathis_median': rad_field_mathis_median,
                     'rad_field_mathis_median_error': \
                        rad_field_mathis_median_error,
                     'rad_field_map': rads,
                     }

            for param_name in cloud_temps[cloud]:
                core_dict[core_name][param_name] = \
                    cloud_temps[cloud][param_name]
        else:
            core_dict[core_name]['dust_temp_median'] = \
                cloud_temps[cloud]['dust_temp_median']
            core_dict[core_name]['dust_temp_median_error'] = \
                cloud_temps[cloud]['dust_temp_median_error']

    return cloud_temps

def add_cloud_params(core_dict, cloud_average=True, load_results=0):

    # Get the data
    # ------------
    dust_temp_dir = '/d/bip3/ezbc/multicloud/data/dust_temp/'
    temp_filename = dust_temp_dir + 'multicloud_dust_temp_5arcmin.fits'
    temp_error_filename = dust_temp_dir + \
            'multicloud_dust_temp_error_5arcmin.fits'
    temp_data, temp_header = fits.getdata(temp_filename, header=True)
    temp_error_data = fits.getdata(temp_error_filename)

    table_temp_dir = '/d/bip3/ezbc/multicloud/data/python_output/dust_temps/'
    cloud_temp_filename = table_temp_dir + 'dust_temps.pickle'

    dust_hist_filename = '/d/bip3/ezbc/multicloud/figures/temps/' + \
                         'multicloud_dust_hist.png'

    # Define filenames
    # ------------
    DIR_DUST = '/d/bip3/ezbc/multicloud/data/dust_temp/'
    DIR_RAD = '/d/bip3/ezbc/multicloud/data/radiation_field/'
    FILENAME_TEMPERATURE = DIR_DUST + 'multicloud_dust_temp_5arcmin.fits'
    FILENAME_TEMPERATURE_ERROR = DIR_DUST + \
            'multicloud_dust_temp_error_5arcmin.fits'
    FILENAME_BETA = DIR_DUST + 'multicloud_dust_beta_5arcmin.fits'
    FILENAME_BETA_ERROR = DIR_DUST + 'multicloud_dust_beta_error_5arcmin.fits'
    FILENAME_RAD = DIR_RAD + 'multicloud_rad_field_5arcmin.fits'
    FILENAME_RAD_ERROR = DIR_RAD + 'multicloud_rad_field_error_5arcmin.fits'

    # Get the data
    # ------------
    temp_data, temp_header = fits.getdata(FILENAME_TEMPERATURE, header=True)
    temp_error_data = fits.getdata(FILENAME_TEMPERATURE_ERROR)
    beta_data, beta_header = fits.getdata(FILENAME_BETA, header=True)
    beta_error_data = fits.getdata(FILENAME_BETA_ERROR)

    # Get the mask for each core
    # --------------------------
    # Create WCS object
    wcs_header = WCS(temp_header)
    if cloud_average:
        if load_cloud_average:
            with open(cloud_temp_filename, 'rb') as f:
                cloud_temps = pickle.load(f)
                for core_name in core_dict:
                    cloud = core_dict[core_name]['cloud']
                    core_dict[core_name]['dust_temp_median'] = \
                        cloud_temps[cloud]['dust_temp_median']
                    core_dict[core_name]['dust_temp_median_error'] = \
                        cloud_temps[cloud]['dust_temp_median_error']
                    core_dict[core_name]['rad_field_draine_median'] = \
                        cloud_temps[cloud]['rad_field_draine_median']
                    core_dict[core_name]['rad_field_draine_median_error'] = \
                        cloud_temps[cloud]['rad_field_draine_median_error']
                    core_dict[core_name]['rad_field_habing_median'] = \
                        cloud_temps[cloud]['rad_field_habing_median']
                    core_dict[core_name]['rad_field_habing_median_error'] = \
                        cloud_temps[cloud]['rad_field_habing_median_error']
        else:
            cloud_temps = run_mc_simulations(core_dict,
                                               wcs_header,
                                               temp_data,
                                               temp_error_data,
                                               beta_data,
                                               beta_error_data,
                                               )


            with open(cloud_temp_filename, 'wb') as f:
                pickle.dump(cloud_temps, f)


            plot_dust_histogram(cloud_temps,
                                limits=[[14, 20], [-0.05, 1.05]],
                                filename=dust_hist_filename)
        write_cloud_temps_table(cloud_temps)

    else:
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
            core_dict[core_name]['dust_temp_errors'] = \
                temp_error_data[~region_mask]

            # Calculate average temp
            dust_temp, dust_temp_error = \
                mystats.calc_cdf_error(temp_data[~region_mask])
            core_dict[core_name]['dust_temp_median'] = \
                dust_temp
            core_dict[core_name]['dust_temp_median_error'] = \
                np.sqrt(np.nansum(temp_error_data[~region_mask]**2)) / \
                    temp_error_data[~region_mask].size

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
        alt_name = 'L1524'
    elif core_name == 'G171.14-17.57':
        alt_name = 'L1504'
    elif core_name == 'G172.93-16.73':
        alt_name = 'B18'
    elif core_name == 'G174.40-13.45':
        alt_name = 'B220'
    elif core_name == 'G174.70-15.47':
        alt_name = 'L1535'
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

    from mystats import sigfig

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = table_dir + 'tables/modelparams_table.tex'

    # get ordered core names
    core_name_list = get_core_names()

    # Open file to be appended
    f = open(filename, 'wb')

    text_param_format ='{0:.0f}$^{{+{1:.0f}}}_{{-{2:.0f}}}$'
    text_param_format_int ='{0:.0f}$^{{+{1:.0f}}}_{{-{2:.0f}}}$'

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

        if 0:
            # add dust temp
            # -------------
            param = core['dust_temp_median']
            param_error = core['dust_temp_median_error']
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
        if 0:
            for model in ( 'krumholz', 'sternberg',):
                if model == 'krumholz':
                    params_to_write = ['phi_cnm',]
                else:
                    params_to_write = ['alphaG',
                                       #'phi_g',
                                       'hi_transition']


                    #print '\nphi_g:'
                    #print core[model]['phi_g']
                    #print core[model]['phi_g_error']

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

        # Krumholz parameters
        # -------------------
        model = 'krumholz'
        param_name = 'phi_cnm'
        param = \
            core[model][param_name]
        param_error = \
            core[model][param_name + '_error']

        param_info = (param, param_error[1], param_error[0])

        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format)

        model = 'krumholz'
        param_name = 'hi_transition'
        param = \
            core[model][param_name]
        param_error = \
            core[model][param_name + '_error']

        param_info = (param, param_error[1], param_error[0])

        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format_int)

        # add H vol dens
        # --------------
        param = core['n_cnm']
        param_error = core['n_cnm_error']
        param_info = (param, param_error[1], param_error[0])
        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format_int)

        # add HI temp
        # -----------
        param = core['T_cnm']
        param_error = core['T_cnm_error']
        param_info = (param, param_error[1], param_error[0])
        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format_int)


        # Sternberg parameters
        # -------------------
        model = 'sternberg'
        param_name = 'alphaG'
        param = \
            core[model][param_name]
        param_error = \
            core[model][param_name + '_error']

        param_info = (param, param_error[1], param_error[0])

        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format)

        model = 'sternberg'
        param_name = 'hi_transition'
        param = \
            core[model][param_name]
        param_error = \
            core[model][param_name + '_error']

        param_info = (param, param_error[1], param_error[0])

        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format_int)

        # add H vol dens
        # --------------
        param = core['n_H']
        param_error = core['n_H_error']
        param_info = (param, param_error[1], param_error[0])
        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format_int)

        # add HI temp
        # -----------
        #param = sigfig(core['T_H'] * 10.0, 2) # convert from units of 1,000 K to 100 K
        #param_error = sigfig(core['T_H_error'] * 10.0, 2)
        param = core['T_H'] * 10.0 # convert from units of 1,000 K to 100 K
        param_error = core['T_H_error'] * 10.0
        param_info = (param, param_error[1], param_error[0])
        row_text = \
            add_row_element(row_text,
                            param_info,
                            text_format=text_param_format_int)


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

def add_cloud_region(core_dict):

    import pyregion as pyr

    data_dir = '/d/bip3/ezbc/multicloud/data/cold_clumps/'

    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'
    #filename = region_dir + 'multicloud_divisions_coldcore_selection.reg'
    filename = region_dir + 'multicloud_divisions.reg'

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        region_name = tag[tag.find('{')+1:tag.find('}')].lower()

        if region_name in ('taurus', 'california', 'perseus'):
            # Format vertices to be 2 x N array
            poly_verts = []
            for i in xrange(0, len(region.coord_list)/2):
                poly_verts.append((region.coord_list[2*i],
                                   region.coord_list[2*i+1]))

            for core in core_dict:
                if core_dict[core]['cloud'] == region_name:
                    core_dict[core]['cloud_region_vertices'] = \
                            np.array(poly_verts)

    return core_dict

def save_core_dict(core_dict):

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = table_dir + 'tables/multicloud_model_summary.pickle'

    with open(filename, 'wb') as f:
        pickle.dump(core_dict, f)

def main():

    LOAD_MC_RESULTS = 0

    # load core summary file
    core_dict = load_cores()

    # average dust temperatures over each core region
    add_cloud_params(core_dict, load_results=LOAD_MC_RESULTS)

    # Add model_analysis
    add_model_params(core_dict)

    # write the latex table
    write_model_params_table(core_dict)

    # save core_dict
    save_core_dict(core_dict)

if __name__ == '__main__':
    main()

