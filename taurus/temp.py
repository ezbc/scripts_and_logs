

if True:
    import grid
    import numpy as np
    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    core_dir = output_dir + 'core_arrays/'

    # load 2mass Av and GALFA HI images, on same grid
    av_data, av_header = load_fits(av_dir + 'taurus_av_k09_regrid.fits',
            return_header=True)
    # load Av image from goldsmith: Pineda et al. 2010, ApJ, 721, 686
    av_data_goldsmith = load_fits(av_dir + \
            'taurus_av_p10_regrid.fits')

    #av_data += - 0.4 # subtracts background of 0.4 mags
    hi_data,h = load_fits(hi_dir + 'taurus_galfa_cube_bin_3.7arcmin.fits',
            return_header=True)

    # make the velocity axis
    velocity_axis = (np.arange(h['NAXIS3']) - h['CRPIX3'] + 1) * h['CDELT3'] + \
            h['CRVAL3']
    velocity_axis /= 1000.

    # Plot NHI vs. Av for a given velocity range
    noise_cube_filename = 'taurus_galfa_cube_bin_3.7arcmin_noise.fits'
    if not path.isfile(hi_dir + noise_cube_filename):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=velocity_axis,
                velocity_noise_range=[90,110], header=h, Tsys=30.,
                filename=hi_dir + noise_cube_filename)
    else:
        noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    nhi_image, nhi_image_error = calculate_NHI(cube=hi_data,
        velocity_axis=velocity_axis, noise_cube = noise_cube,
        velocity_range=[-5,15], return_nhi_error=True)

    nh2_image = calculate_nh2(nhi_image = nhi_image,
            av_image = av_data, dgr = 1.1e-1)
    nh2_image_error = calculate_nh2(nhi_image = nhi_image_error,
            av_image = 0.1, dgr = 1.1e-1)

    hi_sd_image = calculate_sd(nhi_image, sd_factor=1/1.25)
    hi_sd_image_error = calculate_sd(nhi_image_error, sd_factor=1/1.25)

    h2_sd_image = calculate_sd(nh2_image, sd_factor=1/6.25)
    h2_sd_image_error = calculate_sd(nh2_image_error, sd_factor=1/6.25)

    h_sd_image = hi_sd_image + h2_sd_image
    h_sd_image_error = (hi_sd_image_error**0.5 + h2_sd_image_error**0.5)**0.5

    cores = {'L1495':
                {'wcs_position': [15*(4+14/60.), 28+11/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [get_pix_coords(ra=(4,16,30.031),
                                        dec=(27,44,30),
                                        header=h),
                         get_pix_coords(ra=(4,5,20),
                                        dec=(28,28,33),
                                        header=h)]
                 },
             'L1495A':
                {'wcs_position': [15*(4+18/60.), 28+23/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [get_pix_coords(ra=(4,28,23),
                                        dec=(27,52,50),
                                        header=h),
                         get_pix_coords(ra=(4,16,23),
                                        dec=(29,46,5),
                                        header=h)],
                 },
             'B213':
                {'wcs_position': [15*(4+19/60.), 27+15/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [get_pix_coords(ra=(4,23,27),
                                        dec=(26,45,47),
                                        header=h),
                         get_pix_coords(ra=(4,5,25),
                                        dec=(27,38,48),
                                        header=h)],
                 },
             'B220':
                {'wcs_position': [15*(4+41/60.), 26+7/60., 0],
                 'map': None,
                 'threshold': 7,
                 'box': [get_pix_coords(ra=(4,47,49),
                                        dec=(25,31,13),
                                        header=h),
                         get_pix_coords(ra=(4,40,37),
                                        dec=(27,31,17),
                                        header=h)],
                 },
             'L1527':
                {'wcs_position': [15*(4+39/60.), 25+47/60., 0],
                 'map': None,
                 'threshold': 7,
                 'box': [get_pix_coords(ra=(4,40,13),
                                        dec=(24,46,38),
                                        header=h),
                         get_pix_coords(ra=(4,34,35),
                                        dec=(25,56,7),
                                        header=h)],
                 },
             'B215':
                {'wcs_position': [15*(4+23/60.), 25+3/60., 0],
                 'map': None,
                 'threshold': 3,
                 'box': [get_pix_coords(ra=(4,24,51),
                                        dec=(22,36,7),
                                        header=h),
                         get_pix_coords(ra=(4,20,54),
                                        dec=(25,26,31),
                                        header=h)],
                 },
             'L1524':
                {'wcs_position': [15*(4+29/60.), 24+31/60., 0],
                 'map': None,
                 'threshold': 3,
                 'box': [get_pix_coords(ra=(4,30,32),
                                        dec=(22,4,6),
                                        header=h),
                         get_pix_coords(ra=(4,25,33),
                                        dec=(25,0,55),
                                        header=h)],
                 }
                }

    # calculate correlation for cores
    if False:
        limits =[1,22,6,16]

        for core in cores:
        	core_map = np.load(core_dir + core + '.npy')
        	plot_nhi_vs_av(nhi_image,core_map,
                    nhi_image_error = nhi_image_error,
                    av_image_error = 0.1,
                    limits = limits,
                    savedir=figure_dir,
                    plot_type='scatter',
                    filename='taurus_nhi_vs_av_' + core + '_small.png',
                    title=r'N(HI) vs. A$_v$ of Taurus Core ' + core,
                    show=False)

        	plot_sd_vs_av(sd_image, core_map,
                    sd_image_error = sd_image_error,
                    av_image_error = 0.1,
                    limits = limits,
                    savedir=figure_dir,
                    plot_type='scatter',
                    filename='taurus_sd_vs_av_' + core + '_small.png',
                    title=r'$\Sigma_{HI}$ vs. A$_v$ of Taurus Core ' + core,
                    show = False)

    if True:
        limits =[1,22,1,16]
        limits = None

        for core in cores:
            indices = cores[core]['box']
            nhi_image_sub = get_sub_image(nhi_image, indices)
            nhi_image_error_sub = get_sub_image(nhi_image_error, indices)
            av_data_sub = get_sub_image(av_data, indices)
            plot_nhi_vs_av(nhi_image_sub,av_data_sub,
                    nhi_image_error = nhi_image_error_sub,
                    av_image_error = 0.1,
                    limits = limits,
                    savedir=figure_dir,
                    plot_type='scatter',
                    scale='log',
                    filename='taurus_nhi_vs_av_' + core + '_box.png',
                    title=r'N(HI) vs. A$_v$ of Taurus Core ' + core,
                    show=False)

            hi_sd_image_sub = get_sub_image(hi_sd_image, indices)
            hi_sd_image_error_sub = get_sub_image(hi_sd_image_error, indices)


            x0,x1,y0,y1 = int(indices[0][0]), int(indices[0][1]), \
                    int(indices[1][0]), int(indices[1][1])
            hi_sd_image_sub = hi_sd_image[x0:x1,y0:y1]



            print indices[1][0],indices[1][1]
            print hi_sd_image_sub.shape
            print nhi_image_sub.shape
            print indices
            print hi_sd_image.shape
            print hi_sd_image_sub[~np.isnan(hi_sd_image_sub)].max()
            print av_image_sub[~np.isnan(av_image_sub)].max()
            plot_sd_vs_av(hi_sd_image_sub, av_data_sub,
                    sd_image_error = hi_sd_image_error_sub,
                    av_image_error = 0.1,
                    limits = limits,
                    savedir=figure_dir,
                    plot_type='scatter',
                    scale='log',
                    filename='taurus_sd_vs_av_' + core + '_box.png',
                    title=r'$\Sigma_{HI}$ vs. A$_v$ of Taurus Core ' + core,
                    show=False)

            if False:
                h_sd_image_sub = get_sub_image(h_sd_image, indices)
                h_sd_image_error_sub = get_sub_image(h_sd_image_error, indices)
                plot_hisd_vs_hsd(hi_sd_image_sub, h_sd_image_sub,
                        h_sd_image_error = h_sd_image_error_sub,
                        hi_sd_image_error = hi_sd_image_error_sub,
                        limits = limits,
                        savedir=figure_dir,
                        plot_type='scatter',
                        scale='log',
                        filename='taurus_hisd_vs_hsd_' + core + '_box.png',
                        title=r'$\Sigma_{HI}$ vs. $\Sigma_{HI}$ + ' + \
                                '$\Sigma_{H2}$ of Taurus Core ' + core,
                        show=False)

    # Plot heat map of correlations
    if False:
        if correlations is not None:
            correlations_array = plot_correlations(correlations,
                    velocity_centers, velocity_widths,
                    returnimage=True,show=False)
        if False:
            # Print best-fit characteristics
            indices = np.where(cube_correlations_array == \
                    cube_correlations_array.max())
            print 'Maximum correlation values: '
            print str(velocity_centers[indices[0]][0]) + ' km/s center'
            print str(velocity_widths[indices[1]][0]) + ' km/s width'


