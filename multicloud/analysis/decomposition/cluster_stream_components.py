#!/usr/bin/python

import numpy as np
import pickle

def crop_results(results_dict, filename_crop, cropped_size=100, limits=None):

    # Get list positions of small region
    # 278 < glong < 282
    # -37 < glat < -35
    #CROP_LIMITS = [278, 282, -37, -35]

    wcs_pos = results_dict['positions']['wcs']

    indices = np.where(
                        (wcs_pos[:,0] < limits[1]) & # Glong
                        (wcs_pos[:,0] > limits[0]) & # Glong
                        (wcs_pos[:,1] > limits[2]) & # Glat
                        (wcs_pos[:,1] < limits[3]) # Glat
                        )[0]

    # Initialize cropped dict
    cropped_size = len(indices)
    results_dict_crop = {}
    results_dict_crop['velocity_axis'] = results_dict['velocity_axis'].copy()
    results_dict_crop['positions'] = {}
    results_dict_crop['positions']['wcs'] = np.empty((cropped_size, 2))
    results_dict_crop['positions']['pix'] = np.empty((cropped_size, 2))
    results_dict_crop['results'] = []

    # Get portion of full results
    count = 0
    for i in indices:
        wcs = results_dict['positions']['wcs'][i]
        pix = results_dict['positions']['pix'][i]
        results_dict_crop['positions']['wcs'][count] = wcs
        results_dict_crop['positions']['pix'][count] = pix
        results_dict_crop['results'].append(results_dict['results'][i])
        count += 1

    # Print gathered limits of cube
    wcs_array = results_dict_crop['positions']['wcs']
    print('\n\tLimits of cropped cube:')
    print('\tGlon: {0:.0f} to {1:.0f}'.format(np.min(wcs_array[:,0]),
                                            np.max(wcs_array[:,0])))
    print('\tGlat: {0:.0f} to {1:.0f}'.format(np.min(wcs_array[:,1]),
                                            np.max(wcs_array[:,1])))

    # Write the cropped results?
    if filename_crop is not None:
        pickle.dump(results_dict_crop, open(filename_crop, 'w'))

    return results_dict_crop

def crop_outliers(results_dict):

    results_dict['positions']['wcs'] = results_dict['positions']['wcs'].tolist()
    results_dict['positions']['pix'] = results_dict['positions']['pix'].tolist()

    for i in xrange(len(results_dict['results'])):
        pass

def reload_wcs_positions(results_dict):

    from mycoords import make_velocity_axis
    from mydecomposition import plot_nhi_maps, create_synthetic_cube
    import myimage_analysis as myia
    import pyfits as fits
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy import wcs

    # Plot names
    DIR_FIG = '/d/bip3/ezbc/magellanic_stream/figures/'
    FILENAME_FIG = DIR_FIG + 'nhi_map_data_synth.png'

    # Load HI Cube
    DIR_HI = '/d/bip3/ezbc/magellanic_stream/data/hi/'
    FILENAME_CUBE = 'gass_280_-45_1450212515.fits'
    FILENAME_CUBE_SYNTH = DIR_HI + 'cube_synth.npy'

    print('\nLoading data cube...')
    cube, header = fits.getdata(DIR_HI + FILENAME_CUBE, header=True)
    velocity_axis = make_velocity_axis(header)

    header['CUNIT3'] = 'm/s'
    header['CUNIT2'] = 'deg'
    header['CUNIT1'] = 'deg'
    header['CTYPE3'] = 'VOPT'
    header['SPECSYS'] = 'LSRK'

    shape = (cube.shape[1] * cube.shape[2], 2)

    # Create header object
    w = wcs.WCS(header)

    # add position for each spectrum
    for i in xrange(len(results_dict['results'])):
        coords_pix = results_dict['positions']['pix'][i][::-1]
        coords_wcs = w.wcs_pix2world(([coords_pix[0], coords_pix[1], 0],), 0)[0]
        results_dict['positions']['wcs'][i] = np.array(coords_wcs[:2])

    return results_dict

def get_data(load=True, filename=None):

    if not load:
        print('\nFormatting data...')
        # Load the data
        data_dict = {}
        add_data(data_dict)

        # Crop the cube velocities
        crop_data(data_dict)

        # Format data for decomposition
        format_data(data_dict, filename=filename)
    else:
        print('\nLoading data...')
        data_dict = pickle.load(open(filename, 'r'))

    return data_dict

def reformat_results(results_dict, glat_lim=[-np.inf, np.inf],
        glong_lim=[-np.inf,np.inf]):

    ''' Creates array of glong, glat, amp, center, width. Size N x 5.
    '''

    # Initialize
    n_spectra = len(results_dict['results'])
    results_array = []
    results_pix_pos = []

    for i in xrange(n_spectra):
        # positions
        glong = results_dict['positions']['wcs'][i, 0]
        glat = results_dict['positions']['wcs'][i, 1]
        xpix = results_dict['positions']['pix'][i, 0]
        ypix = results_dict['positions']['pix'][i, 1]

        # Include another row if Gaussians are present
        if results_dict['results'][i]['N_components'] > 0:
            fit_params = results_dict['results'][i]['best_fit_parameters']
            ncomps = len(fit_params) / 3
            amps = fit_params[0:ncomps]
            fwhms = fit_params[ncomps:ncomps*2]
            means = fit_params[ncomps*2:ncomps*3]

            # Add each Gaussian component parameters
            for j in xrange(len(amps)):
                keep = ((amps[j] < 100) & \
                        (fwhms[j] < 500) & \
                        (fwhms[j] > 0) & \
                        (glat > glat_lim[0]) & \
                        (glat < glat_lim[1]) & \
                        (glong > glong_lim[0]) & \
                        (glong < glong_lim[1]) & \
                        (means[j] < 500) & \
                        (means[j] > 0)
                        )

                #print amps[j], fwhms[j], means[j]

                if keep:
                    row = [glong, glat, amps[j], fwhms[j], means[j]]
                    results_array.append(row)

                    results_pix_pos.append([xpix, ypix])

    # Gather results
    results_reformatted = {}
    results_reformatted['data'] = np.array(results_array)
    results_reformatted['pos_pix'] = np.array(results_pix_pos)
    results_reformatted['velocity_axis'] = results_dict['velocity_axis'].copy()
    results_reformatted['header'] = ['glong', 'glat', 'amp', 'fwhm', 'mean']

    return results_reformatted

def plot_spectra(results_dict, data_dict):

    from mydecomposition import construct_spectrum

    #vel_axis = np.arange(0, 400, 1)

    vel_axis = results_dict['velocity_axis']

    for i in xrange(len(results_dict['results'])):
        results = results_dict['results'][i]

        # Construct spectrum
        if results['N_components'] > 0:
            spectrum = construct_spectrum(results['best_fit_parameters'],
                                          vel_axis)
        else:
            spectrum = np.zeros(len(vel_axis))

        # Plot scratch plot of fits
        if 1:
            import matplotlib.pyplot as plt
            plt.close(); plt.clf()
            plt.plot(vel_axis,
                     data_dict['data_list'][i])
            plt.plot(vel_axis,
                     spectrum,
                     alpha=0.5,
                     linewidth=3)
            plt.savefig('/d/bip3/ezbc/scratch/spectrum_fit_' + \
                        str(i) + '.png')

    # test

def plot_cluster_nhi_panels(results_ref=None, colors=None, limits=None,
        cube=None, header=None, load_synthetic_cube=False, show=False):

    from mycoords import make_velocity_axis
    from mydecomposition import plot_nhi_map_panels, create_synthetic_cube
    import myimage_analysis as myia
    import pyfits as fits

    # Plot names
    DIR_FIG = '/d/bip3/ezbc/magellanic_stream/figures/'
    FILENAME_FIG = DIR_FIG + 'nhi_maps_components.png'

    # Load HI Cube
    DIR_HI = '/d/bip3/ezbc/magellanic_stream/data_products/hi/'
    FILENAME_CUBE = 'gass_280_-45_1450212515.fits'
    FILENAME_CUBE_SYNTH_BASE = DIR_HI + 'cube_synth_comp'

    velocity_axis = make_velocity_axis(header)

    # Create N(HI) data
    nhi_data = myia.calculate_nhi(cube=cube,
                                  velocity_axis=velocity_axis,
                                  velocity_range=[100, 400],
                                  )

    # Create synthetic cube from fitted spectra
    velocity_axis = results_ref['velocity_axis']

    # get number of unique components
    component_colors = np.unique(colors)
    n_components = len(component_colors)

    nhi_list = []
    nhi_max = 0.0
    for i in xrange(n_components):
        if not load_synthetic_cube:
            print('\n\tCreating synthetic cube ' + str(i+1) + ' of ' + \
                   str(n_components))

            # get the relevant parameters
            indices = np.where(colors == component_colors[i])[0]
            pix_positions = results_ref['pos_pix'][indices]
            fit_params_list = results_ref['data'][indices, 2:]

            print('\n\t\tNumber of components in cube: ' + \
                  '{0:.0f}'.format(len(fit_params_list)))

            cube_synthetic = \
                create_synthetic_cube(pix_positions=pix_positions,
                                      velocity_axis=velocity_axis,
                                      fit_params_list=fit_params_list,
                                      cube_data=cube,
                                      )

            np.save(FILENAME_CUBE_SYNTH_BASE + str(i) + '.npy', cube_synthetic)
        else:
            print('\n\tLoading synthetic cube ' + str(i+1) + ' of ' + \
                   str(n_components))
            cube_synthetic = np.load(FILENAME_CUBE_SYNTH_BASE + str(i) + '.npy')

        # Create N(HI) synthetic
        nhi_synthetic = myia.calculate_nhi(cube=cube_synthetic,
                                           velocity_axis=velocity_axis,
                                           velocity_range=[100, 400],
                                           )

        nhi_list.append(nhi_synthetic)

        nhi_max_temp = np.max(nhi_synthetic)
        if nhi_max_temp > nhi_max:
            nhi_max = nhi_max_temp

    v_limits = [0, nhi_max]

    # Plot the maps together

    plot_nhi_map_panels(nhi_list,
                        header=header,
                        #limits=[278, -37, 282, -35],
                        limits=limits,
                        filename=FILENAME_FIG,
                        #nhi_vlimits=v_limits,
                        show=show
                        )

def get_cube():

    from mycoords import make_velocity_axis
    from mydecomposition import plot_nhi_maps, create_synthetic_cube
    import myimage_analysis as myia
    from astropy.io import fits
    DIR_FIG = '../../figures/'
    FILENAME_FIG = DIR_FIG + 'nhi_map_data_synth.png'

    # Load HI Cube
    DIR_HI = '../../data/hi/'
    FILENAME_CUBE = '../../data/hi/gass_280_-50_1452892031.fits'
    FILENAME_CUBE_SYNTH = DIR_HI + 'cube_synth.npy'

    print('\nLoading data cube...')
    cube_data, header = fits.getdata(DIR_HI + FILENAME_CUBE, header=True)

    header['CUNIT3'] = 'm/s'
    header['CUNIT2'] = 'deg'
    header['CUNIT1'] = 'deg'
    header['CTYPE3'] = 'VOPT'
    header['SPECSYS'] = 'LSRK'

    return cube_data, header

def plot_nhi_maps(results_dict, limits=None, cube_data=None, header=None,
        load_synthetic_cube=False):

    from mycoords import make_velocity_axis
    from mydecomposition import plot_nhi_maps, create_synthetic_cube
    import myimage_analysis as myia
    from astropy.io import fits

    # Plot names
    DIR_FIG = '/d/bip3/ezbc/magellanic_stream/figures/'
    FILENAME_FIG = DIR_FIG + 'nhi_map_data_synth.png'

    # Load HI Cube
    DIR_HI = '/d/bip3/ezbc/magellanic_stream/data_products/hi/'
    FILENAME_CUBE = 'gass_280_-45_1450212515.fits'
    FILENAME_CUBE_SYNTH = DIR_HI + 'cube_synth.npy'

    velocity_axis = make_velocity_axis(header)

    # Create N(HI) data
    nhi_data = myia.calculate_nhi(cube=cube_data,
                                  velocity_axis=velocity_axis,
                                  velocity_range=[100, 400],
                                  )

    # Create synthetic cube from fitted spectra
    velocity_axis = results_dict['velocity_axis']
    if not load_synthetic_cube:
        print('\nCreating synthetic cube...')
        cube_synthetic = create_synthetic_cube(results_dict, cube_data)

        np.save(FILENAME_CUBE_SYNTH, cube_synthetic)
    else:
        print('\nLoading synthetic cube...')
        cube_synthetic = np.load(FILENAME_CUBE_SYNTH)

    # Create N(HI) synthetic
    nhi_synthetic = myia.calculate_nhi(cube=cube_synthetic,
                                       velocity_axis=velocity_axis,
                                       velocity_range=[100, 400],
                                       )

    v_limits = [0, np.max(nhi_synthetic)]

    if 0:
        import matplotlib.pyplot as plt
        plt.close(); plt.clf()
        fig, axes = plt.subplots(2,1)
        axes[0].imshow(nhi_data, origin='lower')
        axes[1].imshow(nhi_synthetic, origin='lower')
        plt.show()

    print('\nPlotting N(HI) maps...')
    # Plot the maps together
    plot_nhi_maps(nhi_data,
                  nhi_synthetic,
                  header=header,
                  #limits=[278, -37, 282, -35],
                  limits=limits,
                  filename=FILENAME_FIG,
                  nhi_1_vlimits=v_limits,
                  nhi_2_vlimits=v_limits,
                  show=True
                  )

def get_PCA(data, n_components=3):

    from mydecomposition import perform_PCA

    data_reduced = perform_PCA(data, n_components=n_components)

    return data_reduced

def get_clusters(data, n_clusters=3, method='kmean'):

    from mydecomposition import cluster_data

    labels = cluster_data(data,
                          n_clusters=n_clusters,
                          method=method,
                          )

    return labels

def plot_cluster_data(data, colors=None, filename=None, show=False,
        labels=None, show_tick_labels=False, zlim=None, ylim=None, xlim=None):

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from myplotting import scatter_contour
    import matplotlib as mpl

    n_components = data.shape[1]

    colors = mpl.cm.rainbow(colors)

    plt.close();

    if labels is None:
        labels = ["1st eigenvector", "2nd eigenvector", "3rd eigenvector"]

    alpha = 1.0 / np.log10(data.shape[0])

    if n_components <= 2:
        fig = plt.figure(1, figsize=(4,4))
        plt.clf()
        ax = fig.add_subplot(111)
        #data[:,0] += np.abs(np.min(data[:,0]) + 1.0)

        ax.scatter(data[:, 0],
                   data[:, 1],
                   color=colors,
                   alpha=alpha)

        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        if not show_tick_labels:
            ax.xaxis.set_ticklabels([])
            ax.yaxis.set_ticklabels([])
        #ax.set_xscale('log')
        if 0:
            scatter_contour(data[:, 0],
                            data[:,1],
                                 threshold=2,
                                 log_counts=0,
                                 levels=5,
                                 ax=ax,
                                 histogram2d_args=dict(bins=50,),
                                 plot_args=dict(marker='o',
                                                linestyle='none',
                                                markeredgewidth=0,
                                                color='black',
                                                alpha=0.4,
                                                markersize=2.5,
                                                ),
                                 contour_args=dict(cmap=plt.cm.binary,)
                                 )

    elif n_components >= 3:
        plt.figure(2, figsize=(6, 6))
        plt.clf()
        # To getter a better understanding of interaction of the dimensions
        # plot the first three PCA dimensions
        fig = plt.figure(1, figsize=(8, 6))
        ax = Axes3D(fig, elev=20, azim=135)
        ax.scatter(data[:, 0], data[:, 1], data[:, 2], c=colors, alpha=alpha)
        #ax.set_title("First three PCA directions")
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_zlabel(labels[2])
        if zlim is not None:
            ax.set_zlim(zlim)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        if not show_tick_labels:
            ax.w_xaxis.set_ticklabels([])
            ax.w_yaxis.set_ticklabels([])
            ax.w_zaxis.set_ticklabels([])

    if filename is not None:
        plt.savefig(filename,)# bbox_inches='tight')
    if show:
        plt.show()

def main():

    ''' Script to cluster gaussian components found from AGD.

    '''

    import os
    from mydecomposition import get_decomposed_data, perform_PCA

    # Set the constants
    DIR_DECOMP = '../../data/decomposition/'
    DIR_FIGURE = '../../figures/'
    FILENAME_DATA = 'agd_gass_data.pickle'
    FILENAME_TRAIN = DIR_DECOMP + 'agd_gass_train.pickle'
    FILENAME_TRAIN_DECOMPOSED = \
        DIR_DECOMP + 'agd_gass_train_decomp.pickle'
    FILENAME_DECOMPOSED = DIR_DECOMP + 'agd_gass_decomp.pickle'
    FILENAME_DECOMPOSED_CROP = \
        DIR_DECOMP + 'agd_gass_decomp_crop.pickle'
    FILENAME_DECOMP_REFORMAT = \
        DIR_DECOMP + 'agd_gass_decomp_reformat.pickle'
    FILENAME_CLUSTERS = DIR_DECOMP + 'agd_gass_clusters.pickle'
    FILENAME_PLOT = DIR_FIGURE + 'pca.png'

    CROP_LIMITS = [285, 290, -39, -35]
    CROP_LIMITS = [285, 297, -43, -35]
    PLOT_NHI = 1
    SHOW_PLOTS = 0

    # Clustering constants
    N_PC = 5
    N_CLUSTERS = 20
    #CLUSTER_METHOD = 'kmeans'
    CLUSTER_METHOD = 'spectral'

    # Switch to working directory

    #os.chdir('/d/bip3/ezbc/magellanic_stream/scripts/gausspy_decomp/')
    os.chdir('/home/ezbc/research/magellanic_stream/scripts/gausspy_decomp/')
    #os.chdir('/d/bip3/ezbc/magellanic_stream/scripts/gausspy_decomp/')

    # Load cropped decomp data?
    LOAD_DECOMP_CROP = 0
    CLOBBER_CROPPED = 0
    LOAD_REFORMATTED_DATA = 0
    LOAD_SYNTHETIC_CUBE = 0
    LOAD_CLUSTERS = 0

    # Remove cropped file?
    if not LOAD_REFORMATTED_DATA:
        if CLOBBER_CROPPED:
            os.system('rm -rf ' + FILENAME_DECOMPOSED_CROP)

        # Load the results
        if os.path.isfile(FILENAME_DECOMPOSED_CROP) and LOAD_DECOMP_CROP:
            print('\nLoading cropped file...')
            results_dict = \
                get_decomposed_data(FILENAME_DATA,
                                    filename_decomposed=\
                                        FILENAME_DECOMPOSED_CROP,
                                    load=True,
                                    )
        elif LOAD_DECOMP_CROP:
            print('\nLoading full results to be cropped...')

            # crop decomposed data
            results_dict = \
                get_decomposed_data(FILENAME_DATA,
                                    filename_decomposed=FILENAME_DECOMPOSED,
                                    load=True,
                                    )

            if 1:
                reload_wcs_positions(results_dict)
                pickle.dump(results_dict, open(FILENAME_DATA, 'w'))

            print('\nCreating cropped file...')
            results_dict = crop_results(results_dict,
                                        FILENAME_DECOMPOSED_CROP,
                                        limits=CROP_LIMITS)
        else:
            print('\nLoading full cube results...')
            # load entire dataset
            results_dict = \
                get_decomposed_data(FILENAME_DATA,
                                    filename_decomposed=FILENAME_DECOMPOSED,
                                    load=True,
                                    )

    if 0:
        FILENAME_DATA = DIR_DECOMP + 'agd_gass_data.pickle'
        data_dict = get_data(load=1, filename=FILENAME_DATA)
        plot_spectra(results_dict, data_dict)

    if LOAD_REFORMATTED_DATA and os.path.isfile(FILENAME_DECOMP_REFORMAT):
        results_ref = pickle.load(open(FILENAME_DECOMP_REFORMAT, 'rb'))
    else:
        results_ref = reformat_results(results_dict)

        pickle.dump(results_ref, open(FILENAME_DECOMP_REFORMAT, 'wb'))

    print('\nPerforming PCA...')
    results_ref['data_reduced'] = get_PCA(results_ref['data'],
                                          n_components=N_PC)

    if LOAD_CLUSTERS and os.path.isfile(FILENAME_CLUSTERS):
        results_ref['cluster_labels'] = \
            pickle.load(open(FILENAME_CLUSTERS, 'rb'))
    else:
        print('\nClustering reduced components...')
        results_ref['cluster_labels'] = \
                get_clusters(results_ref['data_reduced'],
                             n_clusters=N_CLUSTERS,
                             method=CLUSTER_METHOD)
        pickle.dump(results_ref['cluster_labels'],
                    open(FILENAME_CLUSTERS, 'wb'))

    print('\nNumber of components to be ' + \
          'clustered: {0:.0f}'.format(len(results_ref['data'])))


    # Crop the data
    if 0:
        results_ref['data'] = \
            results_ref['data'][np.random.randint(results_ref['data'].shape[0],
                                                  size=500),
                                :]


    print('\nPlotting cluster analysis...')
    plot_cluster_data(results_ref['data_reduced'],
                      colors=results_ref['cluster_labels'],
                      filename=FILENAME_PLOT,
                      show=SHOW_PLOTS)

    print('\nPlotting cluster analysis...')
    plot_cluster_data(results_ref['data'][:, (0,4,3)],
                      colors=results_ref['cluster_labels'],
                      filename=FILENAME_PLOT.replace('.png','_data.png'),
                      #labels=['Glon [deg]', 'Glat [deg]', 'FWHM [km/s]',],
                      labels=['Glon [deg]', 'Velocity [km/s]', 'FWHM [km/s]',],
                      show_tick_labels=True,
                      zlim=[-10,100],
                      show=SHOW_PLOTS)

    # Plot the decomposed HI map with the data
    if PLOT_NHI and not LOAD_REFORMATTED_DATA:
        print('\nPlotting N(HI) maps...')
        cube, header = get_cube()

        plot_nhi_maps(results_dict,
                      limits=CROP_LIMITS,
                      cube_data=cube,
                      header=header)

    if PLOT_NHI:
        print('\nPlotting N(HI) maps...')
        cube, header = get_cube()
        plot_cluster_nhi_panels(results_ref=results_ref,
                                colors=results_ref['cluster_labels'],
                                cube=cube,
                                header=header,
                                limits=CROP_LIMITS,
                                load_synthetic_cube=LOAD_SYNTHETIC_CUBE,
                                show=SHOW_PLOTS)


if __name__ == '__main__':
    main()

