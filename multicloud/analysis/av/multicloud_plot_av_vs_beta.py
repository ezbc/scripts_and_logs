#!/usr/bin/python

''' Calculates the N(HI) map for multicloud

'''

import numpy as np
import matplotlib
matplotlib.use('Agg')

from astropy.io import fits
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from myimage_analysis import calc_region_mask


''' Plotting Functions
'''

def plot_av_vs_beta_grid(plot_dict,
         filename=None, levels=7,
        limits=None, poly_fit=False, contour=True,
        scale=['linear','linear']):

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
    import myplotting as myplt
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    font_scale = 9
    params = {
              'figure.figsize': (3.6, 8),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    axes = AxesGrid(fig, (1,1,1),
                    nrows_ncols=(3, 1),
                    ngrids=3,
                    axes_pad=0.1,
                    aspect=False,
                    label_mode='L',
                    share_all=True)

    for i, cloud_name in enumerate(plot_dict):
        x = plot_dict[cloud_name]['av']
        y = plot_dict[cloud_name]['beta']

        ax = axes[i]

        if 'log' not in scale:
            ax.locator_params(nbins = 6)

        # Drop the NaNs from the images
        indices = np.where((x == x) &\
                           (y == y)
                           )
        x_nonans = x[indices]
        y_nonans = y[indices]

        ax = axes[i]

        if contour:
            if i == 0:
                if limits is None:
                    xmin = np.min(x_nonans)
                    ymin = np.min(y_nonans)
                    xmax = np.max(x_nonans)
                    ymax = np.max(y_nonans)
                    xscalar = 0.15 * xmax
                    yscalar = 0.15 * ymax
                    limits = [xmin - xscalar, xmax + xscalar,
                              ymin - yscalar, ymax + yscalar]

            contour_range = ((limits[0], limits[1]),
                             (limits[2], limits[3]))

            cmap = myplt.truncate_colormap(plt.cm.binary, 0.2, 1, 1000)

            l1 = myplt.scatter_contour(x_nonans.ravel(),
                                 y_nonans.ravel(),
                                 threshold=2,
                                 log_counts=1,
                                 levels=levels,
                                 ax=ax,
                                 #errors=x_error_nonans.ravel(),
                                 histogram2d_args=dict(bins=40,
                                                       range=contour_range),
                                 plot_args=dict(marker='o',
                                                linestyle='none',
                                                color='black',
                                                alpha=0.3,
                                                markersize=2),
                                 contour_args=dict(
                                                   #cmap=plt.cm.binary,
                                                   cmap=cmap,
                                                   #cmap=cmap,
                                                   ),
                                 )
        else:
            image = ax.errorbar(
                                x_nonans.ravel()[::100],
                                y_nonans.ravel()[::100],
                                yerr=(x_error_nonans.ravel()[::100]),
                                alpha=0.2,
                                color='k',
                                marker='^',
                                ecolor='k',
                                linestyle='None',
                                markersize=3
                                )

        c_cycle = myplt.set_color_cycle(num_colors=2, cmap_limits=[0.5, 0.7])

        if poly_fit:
            from scipy.optimize import curve_fit

            p = np.polyfit(x_nonans, y_nonans, 1)

            x_fit = np.linspace(-10, 100, 100)
            y_poly_fit = p[0] * x_fit + p[1]
            ax.plot(x_fit,
                    y_poly_fit,
                    color=c_cycle[1],
                    linestyle='-',
                    linewidth=2,
                    label=\
                        'Polynomial fit: \n' + \
                        r'$\beta$ = {0:.2f}'.format(p[0] * 100.0) + \
                        r'$\frac{{A_V}}{{100 \rm mag}}$ + {0:.2f}'.format(p[1]) + \
                        r'',
                    alpha=0.7,
                    )
        # Annotations
        #anno_xpos = 0.95

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_ylabel(r'$\beta$')
        ax.set_xlabel(r'$A_V$ [mag]')

        if 1:
            loc = 'lower right'
        elif i == 0:
            loc = 'upper left'
        else:
            loc = 'lower left'

        ax.legend(loc=loc)

        ax.annotate(cloud_name.capitalize(),
                    #xytext=(0.96, 0.9),
                    xy=(0.96, 0.96),
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

    if filename is not None:
        plt.draw()
        plt.savefig(filename, bbox_inches='tight', dpi=400)

'''
The main script
'''

def load_results(filename, load_fits=True):

    import pickle
    from astropy.io import fits

    with open(filename, 'rb') as input:
        results = pickle.load(input)
    input.close()

    if load_fits:
        results['data']['av_data'], results['data']['av_header'] = \
                fits.getdata(results['filenames']['av_filename'],
                             header=True)
        results['data']['av_error_data'], results['data']['av_error_header'] = \
                fits.getdata(results['filenames']['av_error_filename'],
                             header=True)
        results['data']['av_data_ref'] = \
                fits.getdata(results['filenames']['av_ref_filename'])
        results['data']['hi_data'], results['data']['hi_header'] = \
                fits.getdata(results['filenames']['hi_filename'], header=True)
        results['data']['hi_data_error'] = \
                fits.getdata(results['filenames']['hi_error_filename'])
        results['data']['co_data'], results['data']['co_header'] = \
                fits.getdata(results['filenames']['co_filename'], header=True)

    return results

def main():

    # define constants
    DIR_FIGURES = '/d/bip3/ezbc/multicloud/figures/'
    DIR_RESULTS = '/d/bip3/ezbc/multicloud/data/python_output/bootstrap_results/'
    DIR_AV = '/d/bip3/ezbc/multicloud/data/av/'
    DIR_BETA = '/d/bip3/ezbc/multicloud/data/dust_temp/'
    DIR_REGION = '/d/bip3/ezbc/multicloud/data/python_output/regions/'

    FILENAME_EXT = '_planck_noint_gaussrange_isotropic_bootstrap_results.pickle'
    FILENAME_PLOT_BASE = DIR_FIGURES + 'dust/av_vs_beta'
    PLOT_FILETYPES = ['png', 'pdf']
    CLOUD_NAMES = ['california', 'perseus', 'taurus']

    beta, beta_header = fits.getdata(DIR_BETA + \
                'multicloud_dust_beta_5arcmin.fits',
            header=True)

    temp, temp_header = fits.getdata(DIR_BETA + \
                'multicloud_dust_temp_5arcmin.fits',
            header=True)
    av, av_header = fits.getdata(DIR_AV + \
                #'multicloud_av_planck_5arcmin.fits',
                'multicloud_av_k09_nan_regrid_planckres.fits',
            header=True)

    region_filename = DIR_REGION + 'multicloud_divisions.reg'

    # get the data needed to plot
    plot_dict = {}
    for cloud_name in CLOUD_NAMES:
        plot_dict[cloud_name] = {}
        cloud_dict = plot_dict[cloud_name]

        # get mask for data
        region_mask = calc_region_mask(region_filename,
                                       av,
                                       av_header,
                                       region_name=cloud_name)

        # load the analysis
        results = load_results(DIR_RESULTS + cloud_name + FILENAME_EXT)

        hi_sd = results['data_products']['hi_sd']
        h2_sd = results['data_products']['h2_sd']

        cloud_dict['contour_image'] = None# hi_sd
        cloud_dict['contours'] = [4, 8]
        cloud_dict['h2_sd'] = h2_sd
        cloud_dict['header'] = results['data']['av_header']
        cloud_dict['beta'] = beta[~region_mask]
        cloud_dict['av'] = av[~region_mask]

        if cloud_name == 'california':
            plot_dict[cloud_name]['limits'] = [75, 60, 30, 38]
        if cloud_name == 'perseus':
            plot_dict[cloud_name]['limits'] = [61, 45, 24, 36]
        if cloud_name == 'taurus':
            plot_dict[cloud_name]['limits'] = [75, 57, 20, 33]

    # plot the 3-panel H2 surface density map
    for filetype in PLOT_FILETYPES:
        plot_av_vs_beta_grid(plot_dict,
                             filename=FILENAME_PLOT_BASE + '.' + filetype,
                             poly_fit=True,
                             #vlimits=[-0.1, 50],
                             )


if __name__ == '__main__':
    main()

