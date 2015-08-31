#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import cloudpy
import matplotlib.pyplot as plt

global debugging
debugging = True
#debugging = False

def plot_planck_vs_2mass(av_k09, av_pl, filename=None, av_error=None,
        contour_plot=True, levels=10, plot_median=True, limits=None,
        scale=('linear','linear'), title = '', gridsize=(100,100),):

    # import external modules
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib import cm
    from astroML.plotting import scatter_contour
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid
    import myplotting as myplt

    # set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # color map
    myplt.set_color_cycle(num_colors=3)

    # Create figure instance
    fig = plt.figure(figsize=(3.6, 3.6))

    axes = AxesGrid(fig, (1,1,1),
                 nrows_ncols=(1, 1),
                 ngrids=1,
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 #cbar_mode='single',
                 cbar_pad=0.1,
                 cbar_size=0.2,
                 )

    # Drop the NaNs from the images
    if type(av_error) is float or av_error is None:
        indices = np.where((av_pl == av_pl) &\
                           (av_k09 == av_k09)
                           )

    av_pl = av_pl[indices]
    av_k09 = av_k09[indices]

    # Create plot
    ax = axes[0]

    if limits is None:
        xmin = np.min(av_k09)
        ymin = np.min(av_k09)
        xmax = np.max(av_pl)
        ymax = np.max(av_pl)
        xscalar = 0.15 * xmax
        yscalar = 0.15 * ymax
        limits = [xmin - xscalar, xmax + xscalar,
                  ymin - yscalar, ymax + yscalar]

    if contour_plot:

        contour_range = ((limits[0], limits[1]),
                         (limits[2], limits[3]))

        cmap = myplt.truncate_colormap(plt.cm.binary, 0.2, 1, 1000)

        l1 = myplt.scatter_contour(av_k09,
                             av_pl,
                             threshold=3,
                             log_counts=1,
                             levels=levels,
                             ax=ax,
                             histogram2d_args=dict(bins=30,
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
        image = ax.errorbar(nhi_nonans.ravel(),
                av_nonans.ravel(),
                yerr=(av_error_nonans.ravel()),
                alpha=0.2,
                color='k',
                marker='^',
                ecolor='k',
                linestyle='None',
                markersize=3
                )

    # Plot sensitivies
    #av_limit = np.median(av_errors[0])
    #ax.axvline(av_limit, color='k', linestyle='--')

    # Plot 1 to 1 pline
    if 1:
        x_fit = np.linspace(-5, 200, 1000)
        p, V = \
                np.polyfit(av_k09, av_pl, deg=1,
                           #w=np.abs(1.0/av_error_nonans.ravel()),
                           cov=True
                           )

        y_poly_fit = p[0] * x_fit + p[1]
    if 1:
        ax.plot(x_fit,
                y_poly_fit,
                color='r',
                linestyle='dashed',
                linewidth=2,
                label=\
                    'Poly fit: \n' + \
                    r'$A_{V,{\rm Planck}}$ = ' + '{0:.1f}'.format(p[0]) + \
                    r'$\times A_{V,{\rm K+09}}$' + \
                    ' + {0:.1f} mag'.format(p[1]),
                alpha=0.7,
                )
    if 0:
        from scipy.interpolate import UnivariateSpline

        spl = UnivariateSpline(av_k09, av_pl)
        ax.plot(x_fit, spl(x_fit), linewidth=3, alpha=0.5)


    if 0:
        print('Bootstrapping...')
        # bootstrap conf intervals
        import scipy as sp
        bootindex = sp.random.random_integers
        nboot = 20
        x_fit = av_k09
        y_fit = p[0] * x_fit  + p[1]
        r = av_pl - y_fit
        for i in xrange(nboot): # loop over n bootstrap samples from the resids
            pc = sp.polyfit(av_k09,  y_fit + r[bootindex(0, len(r)-1, len(r))], 1)
            ax.plot(x_fit, sp.polyval(pc,x_fit), 'k-', linewidth=2, alpha=1.0/float(nboot))


    # Annotations
    anno_xpos = 0.95

    ax.set_xscale(scale[0], nonposx = 'clip')
    ax.set_yscale(scale[1], nonposy = 'clip')

    ax.set_xlim(limits[0],limits[1])
    ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'$A_{V,{\rm K+09}}$ [mag]')
    ax.set_ylabel(r'$A_{V,{\rm Planck}}$ [mag]')
    ax.set_title(title)
    ax.legend(loc='best')

    if filename is not None:
        plt.savefig(filename)

def run_analysis(cloud_name):

    from astropy.io import fits
    from myimage_analysis import calculate_nhi, calc_region_mask
    from mycoords import make_velocity_axis
    from mystats import calc_symmetric_error, calc_logL
    import myio

    # define directory locations
    # --------------------------
    figure_dir = \
        '/d/bip3/ezbc/multicloud/figures/'
    av_dir = '/d/bip3/ezbc/' + cloud_name + '/data/av/'
    hi_dir = '/d/bip3/ezbc/' + cloud_name + '/data/hi/'
    co_dir = '/d/bip3/ezbc/' + cloud_name + '/data/co/'
    core_dir = \
       '/d/bip3/ezbc/' + cloud_name + '/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/' + cloud_name + '/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    background_region_dir = '/d/bip3/ezbc/' + cloud_name + \
                            '/data/python_output/ds9_regions/'
    likelihood_dir = \
            '/d/bip3/ezbc/' + cloud_name + '/data/python_output/nhi_av/'
    results_dir =  '/d/bip3/ezbc/multicloud/data/python_output/'

    # define filenames
    prop_filename = property_dir + \
       cloud_name + '_global_properties.txt'
    hi_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres.fits'
    hi_error_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres_noise.fits'
    co_filename = co_dir + \
       cloud_name + '_co_cfa_cube_regrid_planckres.fits'

    av_planck_filename = av_dir + \
       cloud_name + '_av_planck_tau353_5arcmin.fits'
    av_error_filename = av_dir + \
       cloud_name + '_av_error_planck_tau353_5arcmin.fits'
    av_k09_filename = av_dir + \
       cloud_name + '_av_k09_regrid_planckres.fits'

    av_error = 0.4

    av_background = 0.0

    # Get data
    av_data_pl, av_header_pl = fits.getdata(av_planck_filename, header=True)
    av_data_k09, av_header_k09 = fits.getdata(av_k09_filename, header=True)

    # mask data
    region_filename = region_dir + 'multicloud_divisions.reg'
    region_mask = calc_region_mask(region_filename,
                                   av_data_pl,
                                   av_header_pl,
                                   region_name=cloud_name)

    av_data_pl[region_mask] = np.nan
    av_data_k09[region_mask] = np.nan

    #if debugging:
    if 0:
        import matplotlib.pyplot as plt
        plt.close(); plt.clf()
        plt.imshow(av_data, origin='lower')
        plt.savefig('/usr/users/ezbc/Desktop/avmap.png')

    #hi_data, hi_header = fits.getdata(hi_filename, header=True)
    #co_data, co_header = fits.getdata(co_filename, header=True)
    #hi_vel_axis = make_velocity_axis(hi_header)

    filename = figure_dir + 'av/' + cloud_name + '_planck_vs_k09.png'
    plot_planck_vs_2mass(av_data_k09, av_data_pl,
                         filename=filename,
                         title=cloud_name,
                         limits=[-1, 20, -1, 20])

def main():

    cloud_names = ('perseus',
                   'california',
                   'taurus')

    for cloud in cloud_names:
        run_analysis(cloud)

if __name__ == '__main__':
    main()

