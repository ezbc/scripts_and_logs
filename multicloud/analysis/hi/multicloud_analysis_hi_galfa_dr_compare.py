#!/usr/bin/python

import numpy as np

def plot_spectra(hi_spectrum, hi_vel_axis, hi_std_spectrum=None,
        co_spectrum=None, co_vel_axis=None, gauss_fits=None, limits=None,
        filename='', vel_range=None, comp_num=None):

    ''' Plots a heat map of likelihoodelation values as a function of velocity
    width and velocity center.

    Parameters
    ----------
    cloud : cloudpy.Cloud
        If provided, properties taken from cloud.props.


    '''

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    font_scale = 9
    params = {
              'figure.figsize': (3.6, 3.6 / 1.618),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    axes = AxesGrid(fig, (1,1,1),
                    nrows_ncols=(1, 1),
                    ngrids=1,
                    axes_pad=0,
                    aspect=False,
                    label_mode='L',
                    share_all=True)

    ax = axes[0]

    ax.plot(hi_vel_axis,
            hi_spectrum,
            linewidth=1.5,
            label=r'Median H$\textsc{i}$',
            drawstyle = 'steps-mid'
            )

    if hi_std_spectrum is not None:
        ax.plot(hi_vel_axis,
                hi_std_spectrum,
                linewidth=1.5,
                linestyle='-.',
                label=r'$\sigma_{HI}$',
                )

    if gauss_fits is not None:
        ax.plot(hi_vel_axis, gauss_fits[0],
                alpha=0.4,
                linewidth=3,
                label='Fit',
                )

        label = 'Component'
        for i, comp in enumerate(gauss_fits[1]):
            if comp_num is not None and i in comp_num:
                linewidth = 2
            else:
                linewidth = 1
            ax.plot(hi_vel_axis, comp,
                    linewidth=linewidth,
                    linestyle='--',
                    color='k',
                    label=label,
                    )
            label = None

    if co_spectrum is not None:
        co_scalar = np.nanmax(hi_spectrum) / 2.0
        co_spectrum = co_spectrum / np.nanmax(co_spectrum) * co_scalar
        ax.plot(co_vel_axis,
                co_spectrum,
                #color='k',
                label=r'Median $^{12}$CO $\times$' + \
                       '{0:.0f}'.format(co_scalar),
                drawstyle = 'steps-mid'
                )

    ax.axvspan(vel_range[0],
               vel_range[1],
               alpha=0.3,
               color='k',
               )

    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    ax.legend(loc='upper left')
    ax.set_xlabel('Velocity [km/s]')
    ax.set_ylabel(r'T$_b$ [K]')

    if filename is not None:
        plt.draw()
        plt.savefig(filename, bbox_inches='tight', dpi=100)


def plot_nhi_image(nhi_image=None, header=None, contour_image=None,
        av_image=None, cores=None, props=None, regions=None, title=None,
        limits=None, contours=None, boxes=False, filename=None,
        show=True, hi_vlimits=[None,None], av_vlimits=[None,None],
        cb_text='N(HI)'):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    import pyfits as fits
    import matplotlib.pyplot as plt
    import myplotting as myplt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Set up plot aesthetics
    plt.clf(); plt.close()

    # Color map
    cmap = plt.cm.gnuplot
    #cmap = myplt.reverse_colormap(plt.cm.copper)
    cmap = plt.cm.copper

    # Create figure instance
    fig = plt.figure(figsize=(7.5, 3.6*2))

    if av_image is not None:
        nrows_ncols=(2,1)
        ngrids=2
    else:
        nrows_ncols=(1,1)
        ngrids=1

    axes = AxesGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0.1,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # ------------------
    # NHI image
    # ------------------
    # create axes
    ax = axes[0]

    # show the image
    im = ax.imshow(nhi_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=hi_vlimits[0],
            vmax=hi_vlimits[1],
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)

    ax.locator_params(nbins=6)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')

    # plot limits
    if limits is not None:
        limits_pix = myplt.convert_wcs_limits(limits, header, frame='fk5')
        ax.set_xlim(limits_pix[0],limits_pix[1])
        ax.set_ylim(limits_pix[2],limits_pix[3])

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Write label to colorbar
    #cb.set_label_text(r'$N$(H\textsc{i}) [10$^{20}$ cm$^{-2}$]',)
    if cb_text == 'N(HI)':
        cb.set_label_text(r'$N$(H\textsc{i}) [10$^{20}$ cm$^{-2}$]',)
    else:
        cb.set_label_text(cb_text,)

    # Convert sky to pix coordinates
    #wcs_header = pywcs.WCS(header)
    if type(cores) is dict:
        for core in cores:
            pix_coords = cores[core]['center_pixel']

            anno_color = (0.3, 0.5, 1)

            ax.scatter(pix_coords[0],pix_coords[1],
                    color=anno_color,
                    s=200,
                    marker='+',
                    linewidths=2)

            ax.annotate(core,
                    xy=[pix_coords[0], pix_coords[1]],
                    xytext=(5,5),
                    textcoords='offset points',
                    color=anno_color)

            if boxes:
                rect = ax.add_patch(Polygon(
                    cores[core]['box_vertices'][:, ::-1],
                        facecolor='none',
                        edgecolor=anno_color))

    if regions is not None:
        for region in props['region_name_pos']:
            vertices = np.copy(regions[region]['poly_verts']['pixel'])
            rect = ax.add_patch(Polygon(
                    vertices[:, ::-1],
                    facecolor='none',
                    edgecolor='w'))

    # ------------------
    # Av image
    # ------------------
    if av_image is not None:
        # create axes
        ax = axes[1]
        # show the image
        im = ax.imshow(av_image,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                vmin=av_vlimits[0],
                vmax=av_vlimits[1],
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension [J2000]',)
        ax.set_ylabel('Declination [J2000]',)

        ax.locator_params(nbins=6)

        # colorbar
        cb = ax.cax.colorbar(im)
        cmap.set_bad(color='w')

        # plot limits
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits, header, frame='fk5')
            ax.set_xlim(limits_pix[0],limits_pix[1])
            ax.set_ylim(limits_pix[2],limits_pix[3])

        ax.tick_params(axis='xy', which='major', colors='w')

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='r')

        # Write label to colorbar
        cb.set_label_text(r'$A_V$ [mag]',)

        # Convert sky to pix coordinates
        #wcs_header = pywcs.WCS(header)
        if type(cores) is dict:
            for core in cores:
                pix_coords = cores[core]['center_pixel']

                anno_color = (0.3, 0.5, 1)

                if boxes:
                    rect = ax.add_patch(Polygon(
                        cores[core]['box_vertices'][:, ::-1],
                            facecolor='none',
                            edgecolor='w'))

        if regions is not None:
            for region in props['region_name_pos']:
                vertices = np.copy(regions[region]['poly_verts']['pixel'])
                rect = ax.add_patch(Polygon(
                        vertices[:, ::-1],
                        facecolor='none',
                        edgecolor='w'))

            if props is not None:
                for region in props['region_name_pos']:
                    if region == 'taurus1':
                        region = 'taurus 1'
                    if region == 'taurus2':
                        region = 'taurus 2'
                    ax.annotate(region.capitalize(),
                                xy=props['region_name_pos'][region]['pixel'],
                                xytext=(0,0),
                                textcoords='offset points',
                                color='w',
                                fontsize=10,
                                zorder=10)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        fig.show()

def get_gauss_fit_kwargs(global_args):
    if global_args['cloud_name'] == 'perseus':
        guesses = (28, 3, 5,
                   2, -20, 20)
        ncomps = 2
        ncomps_in_cloud = 1
        co_width = np.abs(10.0 - (-1.0))
    elif global_args['cloud_name'] == 'taurus':
        guesses = (35, 3, 5,
                   5, -15, 20,
                   #3, -2, 2,
                   )
        ncomps = 2
        ncomps_in_cloud = 1
        co_width = np.abs(8.0 - (4.0))
    elif global_args['cloud_name'] == 'california':
        guesses = (50, 3, 5,
                   10, -3, 3,
                   12, -10, 10,
                   3, -30, 10,
                   #2, -20, 20,
                   )
        ncomps = 4
        ncomps_in_cloud = 2
        co_width = np.abs(8.0 - (-4.0))
    gauss_fit_kwargs = {'guesses': guesses,
                        'ncomps': ncomps,
                        'co_width': co_width,
                        #'width_scale': 2,
                        }

    return gauss_fit_kwargs, ncomps_in_cloud

def create_filename_base(global_args):

    # Name of diagnostic files
    if global_args['background_subtract']:
        background_name = '_backsub'
    else:
        background_name = ''

    if global_args['bin_image']:
        bin_name = '_binned'
        global_args['bin_procedure'] = 'all'
    else:
        bin_name = ''
        global_args['bin_procedure'] = 'none'
    if global_args['fixed_width'] is None:
        width_name = ''
        init_vel_width = global_args['init_vel_width']
        vel_center_gauss_fit_kwargs = None
    else:
        if global_args['fixed_width'] == 'gaussfit':
            if global_args['cloud_name'] == 'perseus':
                guesses = (28, 3, 5,
                           2, -20, 20)
                ncomps = 2
            elif global_args['cloud_name'] == 'taurus':
                guesses = (28, 3, 5,
                           5, -30, 20,
                           3, -15, 5,
                           )
                ncomps = 3
            elif global_args['cloud_name'] == 'california':
                guesses = (50, 3, 5,
                           20, -10, 10,
                           3, -45, 10,
                           #2, -20, 20,
                           )
                ncomps = 3
            vel_center_gauss_fit_kwargs = {'guesses': guesses,
                                           'ncomps': ncomps,
                                           #'width_scale': 2,
                                           }
        else:
            vel_center_gauss_fit_kwargs = None
        width_name = '_fixedwidth'
        init_vel_width = global_args['fixed_width']
    if global_args['use_weights']:
        weights_name = '_weights'
        weights_filename = av_dir + \
           global_args['cloud_name'] + '_binweights.fits'
    else:
        weights_name = ''
        weights_filename = None
    if global_args['region'] is None:
        region_name = ''
        global_args['region_name'] = global_args['cloud_name']
    else:
        region_name = '_region' + global_args['region']
        global_args['region_name'] = global_args['cloud_name'] + global_args['region']
    if global_args['av_mask_threshold'] is not None:
        avthres_name = '_avthres'
    else:
        avthres_name = ''
    if not global_args['use_intercept']:
        intercept_name = '_noint'
    else:
        intercept_name = ''
    if global_args['recalculate_likelihoods']:
        error_name = '_errorrecalc'
    else:
        error_name = ''
    if global_args['subtract_comps']:
        compsub_name = '_compsub'
    else:
        compsub_name = ''
    if global_args['use_background']:
        backdgr_name = '_backdgr'
    else:
        backdgr_name = ''
    if global_args['hi_range_calc'] == 'gaussian':
        hi_range_name = 'gaussrange'
    else:
        hi_range_name = 'stdrange'
    if global_args['rotate_cores']:
        rotate_cores_name = '_rotatedcores'
    else:
        rotate_cores_name = ''
    if global_args['vary_phi_g']:
        vary_phi_g_name = '_varyphig'
    else:
        vary_phi_g_name = ''

    filename_extension = global_args['cloud_name']

    return filename_extension, global_args

def get_results(global_args):

    import myio

    print('\nPerforming analysis on ' + global_args['cloud_name'])
    print('=======================' + '=' * len(global_args['cloud_name']))

    # Get the results filename
    filename_base, global_args = create_filename_base(global_args)
    print('\n\tFilename base = \n\t' + filename_base)
    results_dir =  '/d/bip3/ezbc/multicloud/data/python_output/'
    results_filename = results_dir + \
               'bootstrap_results/' + filename_base + \
               '_bootstrap_results.pickle'
    global_args['results_filename'] = results_filename

    exists = myio.check_file(results_filename)

    # either load or perform analysis
    if global_args['load'] and exists:
        print('\n\tLoading results...')
        results_dict = load_results(global_args['results_filename'])

    else:
        results_dict = run_cloud_analysis(global_args)

    # derive col dens images and statistics on MC sim
    #add_results_analysis(results_dict)

    # calculate errors on dgrs and intercept
    #results_dict['params_summary'] = calc_param_errors(results_dict)

    return results_dict

def calc_hi_vel_range(hi_spectrum, hi_vel_axis, gauss_fit_kwargs,
        width_scale=2, co_spectrum=None, co_vel_axis=None, ncomps=1,):

    '''
    Discussion of HI width error from Min:
    (3) For the uncertainty in the HI velocity range, you take the absolute
    difference between the Gaussian HI component center and the median CO peak
    velocity.  In fact, this is sort of the minimum uncertainty we expect.  The
    uncertainty mostly comes from our adoption of +/- 2 sigma(HI) to determine
    the HI velocity range.  In this sense, the maximum uncertainty we expect
    would be (HI velocity range we determine by using +/- 2 sigma(HI) as
    thresholds) - (HI velocity range over which CO emission appears).
    Therefore, a more realistic estimate for the uncertainty in the HI velocity
    range would be, e.g., ("minimum" error + "maximum" error) / 2.  This
    estimate would be a better measure in particular for Taurus, where CO
    emission appears over a relatively small velocity range.

    '''

    from scipy.stats import nanmedian
    from myfitting import fit_gaussians

    co_width = np.copy(gauss_fit_kwargs['co_width'])
    del(gauss_fit_kwargs['co_width'])

    hi_fits = fit_gaussians(hi_vel_axis,
            hi_spectrum, **gauss_fit_kwargs)

    gauss_fit_kwargs['co_width'] = co_width

    # use either the gaussian closest to the CO peak, or the tallest gaussian
    if co_spectrum is not None:
        co_peak_vel = co_vel_axis[co_spectrum == np.nanmax(co_spectrum)][0]
        vel_diffs_center = []
        vel_diffs_edge = []
        for i, param in enumerate(hi_fits[2]):
            vel_center = param[1]
            width = param[2]
            vel_center_diff = np.abs(vel_center - co_peak_vel)
            vel_diffs_center.append(vel_center_diff)
            #vel_diffs_edge.append(vel_center_diff + 2.0 * width)
            vel_diffs_edge.append(4.0 * width)

        # get the velocity range
        vel_diffs_center = np.asarray(vel_diffs_center)
        vel_diffs_edge = np.asarray(vel_diffs_edge)

        sort_indices = np.argsort(vel_diffs_center)
        cloud_comp_num = np.asarray(sort_indices[:ncomps])
        velocity_range = [np.inf, -np.inf]
        for i in cloud_comp_num:
            # get the component
            param = hi_fits[2][i]
            vel_center = param[1]
            width = param[2]

            # set absolute bounds if the component bounds extend beyond
            upper_vel = vel_center + width * width_scale
            lower_vel = vel_center - width * width_scale
            if upper_vel > velocity_range[1]:
                velocity_range[1] = upper_vel
            if lower_vel < velocity_range[0]:
                velocity_range[0] = lower_vel

        # the offset between the co and fitted gaussians will be the HI error
        hi_width_error_min = np.max(vel_diffs_center[cloud_comp_num])

        # max error
        hi_width_error_max = np.max(vel_diffs_edge[cloud_comp_num]) - \
                             co_width
                             #gauss_fit_kwargs['co_width']
        hi_width_error_max = np.abs(velocity_range[1] - velocity_range[0]) - \
                             co_width

        hi_width_error = (hi_width_error_min + hi_width_error_max) / 2.0

        #print 'min error', hi_width_error_min, 'km/s'
        #print 'max error', hi_width_error_max, 'km/s'
        #print 'avg error', hi_width_error, 'km/s'

    else:
        amp_max = -np.Inf
        for i, param in enumerate(hi_fits[2]):
            if param[0] > amp_max:
                amp_max = param[0]
                vel_center = param[1]
                width = param[2] * 4
                cloud_comp_num = i

        velocity_range = [vel_center - width * width_scale,
                          vel_center + width * width_scale]

    return velocity_range, hi_fits, cloud_comp_num, hi_width_error

def run_cloud_analysis(global_args,):

    from astropy.io import fits
    from myimage_analysis import calculate_nhi, calc_region_mask
    import myimage_analysis as myia
    from mycoords import make_velocity_axis
    from mystats import calc_symmetric_error, calc_logL
    import os
    import myio
    import pickle
    import mystats

    cloud_name = global_args['cloud_name']
    region = global_args['region']
    load = global_args['load']
    data_type = global_args['data_type']
    background_subtract = global_args['background_subtract']

    # define directory locations
    # --------------------------
    figure_dir = \
        '/d/bip3/ezbc/multicloud/figures/'
    av_dir = '/d/bip3/ezbc/' + cloud_name + '/data/av/'
    dust_temp_dir = '/d/bip3/ezbc/' + cloud_name + '/data/dust_temp/'
    hi_dir = '/d/bip3/ezbc/' + cloud_name + '/data/hi/'
    co_dir = '/d/bip3/ezbc/' + cloud_name + '/data/co/'
    core_dir = \
       '/d/bip3/ezbc/' + cloud_name + '/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/' + cloud_name + '/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'
    background_region_dir = '/d/bip3/ezbc/' + cloud_name + \
                            '/data/python_output/ds9_regions/'
    results_dir =  '/d/bip3/ezbc/multicloud/data/python_output/'

    av_filename = av_dir + \
       cloud_name + '_av_planck_tau353_5arcmin.fits'
    av_data, av_header = fits.getdata(av_filename, header=True)

    # define filenames
    prop_filename = property_dir + \
       cloud_name + '_global_properties.txt'
    hi_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres.fits'
    hi_dr1_filename = hi_dir + \
       cloud_name + '_hi_galfa_dr1_cube_regrid_planckres.fits'
    hi_error_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres_noise.fits'
    co_filename = co_dir + \
       cloud_name + '_co_cfa_cube_regrid_planckres.fits'

    # Get the filename base to differentiate between different parameters
    filename_base, global_args = create_filename_base(global_args)

    # set up plotting variables
    plot_kwargs = {
                   'figure_dir': figure_dir,
                   'cloud_name': cloud_name,
                   'filename_base': filename_base,
                   'plot_diagnostics': global_args['plot_diagnostics'],
                   #'av_nhi_contour': av_nhi_contour,
                   'av_nhi_contour': True,
                   'av_nhi_limits': [0, 20, -1, 9],
                   #'av_nhi_limits': None,
                    }


    # mask data
    region_filename = region_dir + 'multicloud_divisions.reg'
    region_mask = calc_region_mask(region_filename,
                                   av_data,
                                   av_header,
                                   region_name=global_args['region_name'])


    # Load HI and CO cubes
    hi_data, hi_header = fits.getdata(hi_filename, header=True)
    hi_dr1_data, hi_dr1_header = fits.getdata(hi_dr1_filename, header=True)
    co_data, co_header = fits.getdata(co_filename, header=True)


    hi_data[:, region_mask] = np.nan
    hi_dr1_data[:, region_mask] = np.nan
    co_data[:, region_mask] = np.nan

    hi_vel_axis = make_velocity_axis(hi_header)
    co_vel_axis = make_velocity_axis(co_header)

    # Load HI error
    if global_args['clobber_hi_error']:
        print('\n\tCalculating HI noise cube...')
        os.system('rm -rf ' + hi_error_filename)
        hi_data_error = \
            myia.calculate_noise_cube(cube=hi_data,
                                      velocity_axis=hi_vel_axis,
                                      velocity_noise_range=[-110,-90, 90,110],
                                      Tsys=30.0,
                                      filename=hi_error_filename)
    else:
        hi_data_error = fits.getdata(hi_error_filename)


    # Derive N(HI)
    # -------------------------------------------------------------------------
    # get fit kwargs
    gauss_fit_kwargs, ncomps_in_cloud = get_gauss_fit_kwargs(global_args)

    # derive spectra or load
    spectra_filename = results_dir + 'spectra/' + global_args['cloud_name'] + \
            '_spectra.pickle'
    spectra_dr1_filename = results_dir + 'spectra/' + \
                           global_args['cloud_name'] + \
                           '_spectra_dr1.pickle'
    load_spectra = myio.check_file(spectra_filename,
                                   clobber=global_args['clobber_spectra'])
    if load_spectra:
        hi_spectrum, hi_std_spectrum, co_spectrum = \
                myio.load_pickle(spectra_filename)
        hi_dr1_spectrum, hi_std_dr1_spectrum, co_spectrum = \
                myio.load_pickle(spectra_dr1_filename)
    else:
        print('\n\tCalculating spectra...')
        if global_args['smooth_hi_to_co_res']:
            from astropy.convolution import Gaussian2DKernel, convolve
            # Create kernel
            # one pix = 5 arcmin, need 8.4 arcmin for CO res
            # The beamsize is the FWHM. The convolution kernel needs the
            # standard deviation
            hi_res = 1.0
            co_res = 8.4 / 5.0
            width = (co_res**2 - hi_res**2)**0.5
            std = width / 2.355
            g = Gaussian2DKernel(width)

            # Convolve data
            hi_data_co_res = np.zeros(hi_data.shape)
            for i in xrange(hi_data.shape[0]):
                hi_data_co_res[i, :, :] = \
                    convolve(hi_data[i, :, :], g, boundary='extend')

            hi_dr1_data_co_res = np.zeros(hi_dr1_data.shape)
            for i in xrange(hi_dr1_data.shape[0]):
                hi_dr1_data_co_res[i, :, :] = \
                    convolve(hi_dr1_data[i, :, :], g, boundary='extend')

        hi_spectrum = myia.calc_spectrum(hi_data_co_res)
        hi_std_spectrum = myia.calc_spectrum(hi_data_co_res,
                                             statistic=np.nanstd)
        hi_dr1_spectrum = myia.calc_spectrum(hi_dr1_data_co_res)
        hi_std_dr1_spectrum = myia.calc_spectrum(hi_dr1_data_co_res,
                                             statistic=np.nanstd)
        co_spectrum = myia.calc_spectrum(co_data)
        myio.save_pickle(spectra_filename,
                         (hi_spectrum, hi_std_spectrum, co_spectrum))
        myio.save_pickle(spectra_dr1_filename,
                         (hi_dr1_spectrum, hi_std_dr1_spectrum, co_spectrum))

    if global_args['hi_range_calc'] == 'gaussian':
        velocity_range, gauss_fits, comp_num, hi_range_error = \
                calc_hi_vel_range(hi_spectrum,
                                  hi_vel_axis,
                                  gauss_fit_kwargs,
                                  co_spectrum=co_spectrum,
                                  co_vel_axis=co_vel_axis,
                                  ncomps=ncomps_in_cloud,
                                  )
        global_args['vel_range_error'] = hi_range_error
        velocity_range_dr1, gauss_fits_dr1, comp_num_dr1, hi_range_error_dr1 = \
                calc_hi_vel_range(hi_dr1_spectrum,
                                  hi_vel_axis,
                                  gauss_fit_kwargs,
                                  co_spectrum=co_spectrum,
                                  co_vel_axis=co_vel_axis,
                                  ncomps=ncomps_in_cloud,
                                  )
    else:
        velocity_range = [-5, 15]
        gauss_fits = None
        comp_num = None

    hi_range_kwargs = {
                       'velocity_range': velocity_range,
                       'gauss_fits': gauss_fits,
                       'comp_num': comp_num,
                       'hi_range_error': hi_range_error,
                       'vel_range': velocity_range,
                       'gauss_fit_kwargs': gauss_fit_kwargs,
                       }

    # plot the results
    # --------------------------------------------------------------------------
    filename = plot_kwargs['figure_dir'] + \
               'spectra/' + plot_kwargs['filename_base'] + \
               '_spectra_dr2.png'
    print('Saving\neog ' + filename + ' &')
    plot_spectra(hi_spectrum,
                 hi_vel_axis,
                 hi_std_spectrum=hi_std_spectrum,
                 gauss_fits=gauss_fits,
                 comp_num=comp_num,
                 co_spectrum=co_spectrum,
                 co_vel_axis=co_vel_axis,
                 vel_range=velocity_range,
                 filename=filename,
                 limits=[-50, 30, -10, 70],
                 )

    # DR1 data
    filename = plot_kwargs['figure_dir'] + \
               'spectra/' + plot_kwargs['filename_base'] + \
               '_spectra_dr1.png'
    print('Saving\neog ' + filename + ' &')
    plot_spectra(hi_dr1_spectrum,
                 hi_vel_axis,
                 hi_std_spectrum=hi_std_dr1_spectrum,
                 gauss_fits=gauss_fits_dr1,
                 comp_num=comp_num_dr1,
                 co_spectrum=co_spectrum,
                 co_vel_axis=co_vel_axis,
                 vel_range=velocity_range_dr1,
                 filename=filename,
                 limits=[-50, 30, -10, 70],
                 )

    # use the vel range to derive N(HI)
    nhi_image, nhi_image_error = \
        calculate_nhi(cube=hi_data,
                      velocity_axis=hi_vel_axis,
                      velocity_range=velocity_range,
                      noise_cube=hi_data_error,
                      return_nhi_error=True,
                      )
    # use the vel range to derive N(HI)
    nhi_image_dr1 = \
        calculate_nhi(cube=hi_dr1_data,
                      velocity_axis=hi_vel_axis,
                      velocity_range=velocity_range_dr1,
                      )

    # mask for erroneous pixels
    mask_nhi = (nhi_image < 0) & (nhi_image_dr1 < 0)
    nhi_image[mask_nhi] = np.nan
    nhi_image_dr1[mask_nhi] = np.nan

    # Plot residuals between nhi maps
    filename = plot_kwargs['figure_dir'] + \
               'maps/' + plot_kwargs['filename_base'] + \
               '_nhi_dr2_dr1_residuals.png'
    print('Saving\neog ' + filename + ' &')
    plot_nhi_image(nhi_image=nhi_image / nhi_image_dr1,
                   header=hi_header,
                   limits=None,
                   filename=filename,
                   show=0,
                   cb_text='DR2 / DR1'
                   #hi_vlimits=None,
                   )

def main():

    import itertools

    results = {}

    clouds = (
              'taurus',
              'california',
              'perseus',
              )

    data_types = (
                  'planck',
                  #'lee12',
                  #'planck_lee12mask',
                  #'k09',
                  )
    recalculate_likelihoods = (
                               #True,
                               False,
                               )
    bin_image = (
                 #True,
                 False,
                 )

    init_vel_width = (#20,
                      20,
                      #200,
                      )
    fixed_width = (
                   #20,
                   #'gaussfit',
                   None,
                   #50,
                   )

    hi_range_calc = ('gaussian',
                     #'std',
                     )

    use_intercept = (
                     False,
                     #True,
                     )

    use_background = (
                      False,
                      #True,
                      )
    av_mask_threshold = (
                         None,
                         #1.2,
                         #2.5,
                         )

    regions = (None,
               #'1',
               #'2'
               )

    subtract_comps = (#True,
                      False,
                      )

    radiation_type = (#'beamed',
                      'isotropic',
                      )

    rotate_cores = (
                    False,
                    #True,
                    )

    vary_phi_g = (
                    #True,
                    False,
                    )

    elements = (clouds, data_types, recalculate_likelihoods, bin_image,
            init_vel_width, fixed_width, use_intercept, av_mask_threshold,
            regions, subtract_comps, use_background, hi_range_calc,
            radiation_type, rotate_cores, vary_phi_g)

    permutations = list(itertools.product(*elements))

    print('Number of permutations to run: ' + str(len(permutations)))

    #for cloud in clouds:
    for permutation in permutations:
        global_args = {
                'cloud_name':permutation[0],
                'load': 0,
                'load_props': 0,
                'data_type' : permutation[1],
                'background_subtract': 0,
                'recalculate_likelihoods': permutation[2],
                'bin_image': permutation[3],
                'use_weights': 0,
                'init_vel_width': permutation[4],
                'fixed_width': permutation[5],
                'use_intercept': permutation[6],
                'av_mask_threshold': permutation[7],
                'binned_data_filename_ext': '_bin',
                'likelihood_resolution': 'coarse',
                'region': permutation[8],
                'subtract_comps': permutation[9],
                'plot_diagnostics': 0,
                'use_background': permutation[10],
                'clobber_spectra': 0,
                'smooth_hi_to_co_res': 1,
                'clobber_hi_error': 0,
                'sim_hi_error': True,
                'hi_range_calc': permutation[11],
                #'num_bootstraps': 10000,
                'num_bootstraps': 10,
                'num_resid_bootstraps': 100,
                'bootstrap_fit_residuals': False,
                'calculate_median_error': False,
                'multiprocess': 1,
                'radiation_type': permutation[12],
                'rotate_cores': permutation[13],
                'vary_phi_g': permutation[14],
                }
        run_analysis = False
        if global_args['data_type'] in ('planck_lee12mask', 'lee12'):
            if global_args['cloud_name'] == 'perseus':
                run_analysis = True
        else:
            if global_args['cloud_name'] == 'california':
                if global_args['region'] is None:
                    run_analysis = True
            else:
                run_analysis = True

        if run_analysis:
            results[global_args['cloud_name']] = \
                    get_results(global_args)

    plot_multicloud_results(results)

if __name__ == '__main__':
    main()


