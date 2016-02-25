#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import myimage_analysis as myia
from multiprocessing.queues import Queue
import mygeometry as myg
import scipy
from local_module_plotting import *
from local_module_dataprep import *
from local_module_multiprocessing import *
from local_module_regions import *
from local_module_fitting import *
from local_module_bootstrapping import *

def plot_av_vs_nhi(av, nhi, fit_params=None, filename=None, av_error=None,
        contour_plot=True, levels=10, plot_median=True, limits=None,
        scale=('linear','linear'), title = '', gridsize=(100,100), std=None):

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
    #plt.rcdefaults()

    # color map
    cmap = plt.cm.gnuplot

    # color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    params = {'axes.color_cycle': color_cycle, # colors of different plots
             }
    #plt.rcparams.update(params)

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
        indices = np.where((av == av) &\
                           (nhi == nhi)
                           )
    elif type(av_error) is np.ndarray or \
            type(av_error) is np.ma.core.MaskedArray:
        indices = np.where((av == av) &\
                           (nhi == nhi) &\
                           (av_error == av_error) &\
                           (av_error > 0)
                           )
        av_error_nonans = av_error[indices]

    av_nonans = av[indices]
    nhi_nonans = nhi[indices]

    # Create plot
    ax = axes[0]

    if limits is None:
        xmin = np.min(nhi_nonans)
        ymin = np.min(av_nonans)
        xmax = np.max(nhi_nonans)
        ymax = np.max(av_nonans)
        xscalar = 0.15 * xmax
        yscalar = 0.15 * ymax
        limits = [xmin - xscalar, xmax + xscalar,
                  ymin - yscalar, ymax + yscalar]

    if contour_plot:
        contour_range = ((limits[0], limits[1]),
                         (limits[2], limits[3]))

        cmap = myplt.truncate_colormap(plt.cm.binary, 0.2, 1, 1000)

        l1 = myplt.scatter_contour(nhi_nonans.ravel(),
                             av_nonans.ravel(),
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

    if plot_median:
        from scipy.stats import nanmedian, binned_statistic
        x_median = np.arange(np.min(nhi_nonans), np.max(nhi_nonans), 0.3)
        x_median = np.arange(np.min(nhi_nonans), np.max(nhi_nonans), 1)
        #x_median = np.arange(6.5, 9, 0.3)
        y_median, x_median = binned_statistic(nhi_nonans, av_nonans,
                                    statistic=nanmedian,
                                    bins=x_median)[:2]
        x_median = x_median[:-1]
        x_median = x_median[~np.isnan(y_median)]
        y_median = y_median[~np.isnan(y_median)]

        ax.plot(x_median,
                y_median,
                alpha=1,
                color='r',
                marker='s',
                linestyle='None',
                label='Median value',
                markersize=4.5
                )


    # Plot sensitivies
    #av_limit = np.median(av_errors[0])
    #ax.axvline(av_limit, color='k', linestyle='--')

    # Plot 1 to 1 pline
    nhi_fit = np.linspace(0, 100, 2)
    av_fit = fit_params['dgr'] * nhi_fit + fit_params['intercept']
    if 'dgr_error' in fit_params:
        dgr_error_text = r'$^{+%.2f}_{-%.2f}$ ' % fit_params['dgr_error']
    else:
        dgr_error_text = ''
    if 'intercept_error' in fit_params:
        intercept_error_text = \
                r'$^{+%.2f}_{-%.2f}$' % fit_params['intercept_error'],
    else:
        intercept_error_text = ''

    ax.plot(nhi_fit,
            av_fit,
            #color='0.5',
            linestyle='--',
            linewidth=2,
            alpha=0.7,
            label=\
                'Fit: \n' + \
                'DGR = {0:.2f}'.format(fit_params['dgr']) + \
                dgr_error_text + \
                '\nIntercept = {0:.2f}'.format(fit_params['intercept']) + \
                intercept_error_text
            )

    # Annotations
    anno_xpos = 0.95

    ax.set_xscale(scale[0], nonposx = 'clip')
    ax.set_yscale(scale[1], nonposy = 'clip')

    ax.set_xlim(limits[0],limits[1])
    ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'$N($H$\textsc{i}) \times\,10^{20}$ cm$^{-2}$')
    ax.set_ylabel(r'$A_V$ [mag]')
    ax.set_title(title)
    ax.legend(loc='best')

    if filename is not None:
        plt.savefig(filename)

def plot_av_vs_nhi(av_grid, nhi_grid, fit_params=None, filename=None,
        av_error=None, contour_plot=True, levels=7, plot_median=True,
        limits=None, scale=('linear','linear'), title = '', gridsize=(100,100),
        std=None):

    # import external modules
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib import cm
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid
    import myplotting as myplt

    # set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    #plt.rcdefaults()

    # color map
    myplt.set_color_cycle(num_colors=2, cmap_limits=[0.2, 0.8])

    ngrids = len(av_grid)

    # Create figure instance
    if ngrids == 1:
        ysize = 3.6
    else:
        ysize = 7.5
    fig = plt.figure(figsize=(3.6, ysize))

    axes = AxesGrid(fig, (1,1,1),
                 nrows_ncols=(ngrids, 1,),
                 ngrids=ngrids,
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 #cbar_mode='single',
                 cbar_pad=0.1,
                 cbar_size=0.2,
                 )

    # Create plot
    for i, ax in enumerate(axes):

        av = av_grid[i]
        nhi = nhi_grid[i]

        # Drop the NaNs from the images
        if type(av_error) is float or av_error is None:
            indices = np.where((av == av) &\
                               (nhi == nhi)
                               )
        elif type(av_error) is np.ndarray or \
                type(av_error) is np.ma.core.MaskedArray:
            indices = np.where((av == av) &\
                               (nhi == nhi) &\
                               (av_error == av_error) &\
                               (av_error > 0)
                               )
            av_error_nonans = av_error[indices]

        av_nonans = av[indices]
        nhi_nonans = nhi[indices]

        ax = axes[i]

        if i == 0:
            if limits is None:
                xmin = np.min(nhi_nonans)
                ymin = np.min(av_nonans)
                xmax = np.max(nhi_nonans)
                ymax = np.max(av_nonans)
                xscalar = 0.15 * xmax
                yscalar = 0.15 * ymax
                limits = [xmin - xscalar, xmax + xscalar,
                          ymin - yscalar, ymax + yscalar]

        if contour_plot:
            contour_range = ((limits[0], limits[1]),
                             (limits[2], limits[3]))

            cmap = myplt.truncate_colormap(plt.cm.binary, 0.2, 1, 1000)

            l1 = myplt.scatter_contour(nhi_nonans.ravel(),
                                 av_nonans.ravel(),
                                 threshold=2,
                                 log_counts=1,
                                 levels=levels,
                                 ax=ax,
                                 #errors=av_error_nonans.ravel(),
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

        if plot_median:
            from scipy.stats import nanmedian, binned_statistic
            x_median = np.linspace(np.min(nhi_nonans), np.max(nhi_nonans), 6)
            #x_median = np.arange(6.5, 9, 0.3)
            y_median, x_median = binned_statistic(nhi_nonans, av_nonans,
                                        statistic=nanmedian,
                                        bins=x_median)[:2]
            x_median = x_median[:-1]
            x_median = x_median[~np.isnan(y_median)]
            y_median = y_median[~np.isnan(y_median)]
            if i == 0:
                label = 'Median value'
            else:
                label = ''

            ax.plot(x_median,
                    y_median,
                    alpha=1,
                    #color='r',
                    marker='s',
                    linestyle='None',
                    label=label,
                    markersize=4.5
                    )

        # Plot sensitivies
        #av_limit = np.median(av_errors[0])
        #ax.axvline(av_limit, color='k', linestyle='--')
        if 'dgr_error' in fit_params:
            dgr_error_text = r'$^{+%.2f}_{-%.2f}$ ' % fit_params['dgr_error']
        else:
            dgr_error_text = ''
        if 'intercept_error' in fit_params:
            intercept_error_text = \
                    r'$^{+%.2f}_{-%.2f}$' % fit_params['intercept_error'],
        else:
            intercept_error_text = ''

        # Plot 1 to 1 pline
        nhi_fit = np.linspace(0, 100, 1000)
        dgr_error_text = \
                r'$^{+%.2f}_{-%.2f}$ ' % fit_params['dgr_cloud_error']
        intercept_error_text = \
            r'$^{+%.2f}_{-%.2f}$ ' % fit_params['intercept_error']
        cloud_text = 'Cloud DGR = {0:.2f}'.format(fit_params['dgr_cloud']) + \
                dgr_error_text + \
                ''#+'\nIntercept = {0:.2f}'.format(fit_params['intercept']) + \
                #intercept_error_text
        #if i == 0:
        #    ax.set_title('Full line-of-sight')
        #    av_fit = fit_params['dgr_cloud'] * nhi_fit + \
        #             fit_params['dgr_background'] * nhi_fit + \
        #             fit_params['intercept']
        #    label=''
        #    if ngrids == 1:
        #        label = cloud_text
        if i == 0:
            ax.set_title('Cloud')
            av_fit = fit_params['dgr_cloud'] * nhi_fit + \
                     fit_params['intercept']
            label = cloud_text
        elif i == 1:
            ax.set_title('Background')
            av_fit = fit_params['dgr_background'] * nhi_fit + \
                     fit_params['intercept']
            dgr_error_text = \
                    r'$^{+%.2f}_{-%.2f}$ ' % fit_params['dgr_background_error']
            intercept_error_text = \
                r'$^{+%.2f}_{-%.2f}$ ' % fit_params['intercept_error']
            label = 'Background DGR = ' + \
                    '{0:.2f}'.format(fit_params['dgr_background']) + \
                    dgr_error_text + \
                    '\nIntercept = {0:.2f}'.format(fit_params['intercept']) + \
                    intercept_error_text

        ax.plot(nhi_fit,
                av_fit,
                #color='r',
                linestyle='--',
                linewidth=2,
                alpha=0.7,
                label=label,
                )

        # Annotations
        anno_xpos = 0.95

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'$N($H$\textsc{i}) \times\,10^{20}$ cm$^{-2}$')
        ax.set_ylabel(r'$A_V$ [mag]')
        ax.legend(loc='best')

    if filename is not None:
        plt.savefig(filename)

def plot_bootstrap_dist(dgrs, intercepts, limits=None, filename=None,
        levels=4, axis_labels=['',''], contour_plot=True):

    ''' Plots a heat map of likelihoodelation values as a function of velocity
    width and velocity center.

    Parameters
    ----------
    cloud : cloudpy.Cloud
        If provided, properties taken from cloud.props.


    '''

    # Import external modules
    import numpy as np
    import math
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import myplotting as myplt

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    font_scale = 9
    params = {
              'figure.figsize': (3.6, 3.6),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    fig, ax = plt.subplots()

    if contour_plot:
        if limits is None:
            xmin = np.min(dgrs)
            ymin = np.min(intercepts)
            xmax = np.max(dgrs)
            ymax = np.max(intercepts)
            xscalar = 0.15 * xmax
            yscalar = 0.15 * ymax
            limits = [xmin - xscalar, xmax + xscalar,
                      ymin - yscalar, ymax + yscalar]

        contour_range = ((limits[0], limits[1]),
                         (limits[2], limits[3]))

        cmap = myplt.truncate_colormap(plt.cm.binary, 0.2, 1, 1000)

        l1 = myplt.scatter_contour(dgrs.ravel(),
                             intercepts.ravel(),
                             threshold=4,
                             log_counts=1,
                             levels=[0.99, 0.95, 0.68],
                             ax=ax,
                             #errors=av_error_nonans.ravel(),
                             histogram2d_args=dict(bins=20,
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
        image = ax.plot(dgrs.ravel(),
                intercepts.ravel(),
                alpha=0.2,
                color='k',
                marker='o',
                linestyle='None',
                markersize=3
                )
    #ax.set_xscale(scale[0], nonposx = 'clip')
    #ax.set_yscale(scale[1], nonposy = 'clip')

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')

def plot_bootstrap_hist(dgrs, limits=None, filename=None,
        axis_label='', statistics=None):

    ''' Plots a heat map of likelihoodelation values as a function of velocity
    width and velocity center.

    Parameters
    ----------
    cloud : cloudpy.Cloud
        If provided, properties taken from cloud.props.


    '''

    # Import external modules
    import numpy as np
    import scipy
    import math
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import myplotting as myplt

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    font_scale = 9
    params = {
              'figure.figsize': (3.6, 3.6),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    fig, ax = plt.subplots()

    counts, bins = np.histogram(dgrs,bins=40)
    counts = counts / float(np.max(counts))

    if 1:
        ax.plot(bins[:-1],
                counts,
                color='k',
                drawstyle='steps-mid',
                label='PDF',
                linewidth=1.5,
                )

    # Compute the CDF
    dgrs_sorted = np.sort(dgrs)
    cdf = np.cumsum(dgrs_sorted)
    cdf = cdf / np.max(cdf)
    cdf = 1. * np.arange(len(dgrs_sorted)) / (len(dgrs_sorted) - 1)
    ax.plot(dgrs_sorted,
            cdf,
            color='k',
            linestyle='--',
            drawstyle='steps-mid',
            label='CDF',
            linewidth=1.5,
            )

    if statistics is not None:
        ax.axvspan(statistics[0] - statistics[1][0],
                   statistics[0] + statistics[1][1],
                   color='k',
                   linewidth=1,
                   alpha=0.2)
        ax.axvline(statistics[0],
                   color='k',
                   linestyle='--',
                   linewidth=3,
                   label='Median',
                   alpha=1)

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    ax.legend(loc='best')

    # Adjust asthetics
    ax.set_xlabel(axis_label)
    #ax.set_ylabel('Counts')
    ax.set_ylabel('CDF')
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')

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

def plot_spectra_grid(spectra_list, hi_range_kwargs_list=None,
        names_list=None, hi_vel_axis=None, co_vel_axis=None, filename=None,
        limits=None,):

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
              'figure.figsize': (3.5, 6),
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    myplt.set_color_cycle(num_colors=4, cmap_limits=[0, 0.9])

    # Create figure instance
    fig = plt.figure()

    axes = AxesGrid(fig, (1,1,1),
                    nrows_ncols=(3, 1),
                    ngrids=3,
                    axes_pad=0.1,
                    aspect=False,
                    label_mode='L',
                    share_all=True)

    for i in xrange(0, len(names_list)):
        hi_spectrum, hi_std_spectrum, co_spectrum = spectra_list[i]
        cloud_name = names_list[i]
        hi_range_kwargs = hi_range_kwargs_list[i]
        gauss_fits = hi_range_kwargs['gauss_fits']
        comp_num = hi_range_kwargs['comp_num']
        vel_range = hi_range_kwargs['vel_range']
        ax = axes[i]

        ax.locator_params(nbins = 6)

        ax.plot(hi_vel_axis,
                hi_spectrum,
                linewidth=1.5,
                label=r'Median H$\textsc{i}$',
                drawstyle = 'steps-mid'
                )

        if 0:
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
            for j, comp in enumerate(gauss_fits[1]):
                if comp_num is not None and j in comp_num:
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
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel(r'T$_b$ [K]')

    if filename is not None:
        plt.draw()
        plt.savefig(filename, bbox_inches='tight', dpi=100)

def plot_av_vs_nhi_grid(av_list, nhi_list, names_list=None,
        av_error_list=None, fit_params_list=None, filename=None, levels=7,
        limits=None, poly_fit=False, plot_median=True, contour=True,
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

    for i in xrange(0, len(names_list)):
        cloud_name = names_list[i]
        x = av_list[i]
        x_error = av_error_list[i]
        y = nhi_list[i]
        fit_params = fit_params_list[i]
        ax = axes[i]

        if 'log' not in scale:
            ax.locator_params(nbins = 6)

        # Drop the NaNs from the images
        if type(x_error) is float or x_error is None:
            indices = np.where((x == x) &\
                               (y == y)
                               )
            x_error_nonans = x_error
        elif type(x_error) is np.ndarray or \
                type(x_error) is np.ma.core.MaskedArray:
            indices = np.where((x == x) &\
                               (y == y) &\
                               (x_error == x_error) &\
                               (x_error > 0)
                               )
            x_error_nonans = x_error[indices]

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
        if plot_median:
            from scipy.stats import nanmedian, binned_statistic
            y_median = np.linspace(np.min(y_nonans), np.max(y_nonans), 6)
            #x_median = np.arange(6.5, 9, 0.3)
            x_median, y_median = binned_statistic(y_nonans, x_nonans,
                                        statistic=nanmedian,
                                        bins=y_median)[:2]
            y_median = y_median[:-1]
            x_median = x_median[~np.isnan(y_median)]
            y_median = y_median[~np.isnan(y_median)]
            if i == 0:
                label = 'Median value'
            else:
                label = ''


            if 1:
                ax.plot(x_median,
                        y_median,
                        #color='r',
                        color=c_cycle[0],
                        marker='s',
                        linestyle='None',
                        label=label,
                        zorder=1000,
                        #alpha=0.5,
                        markersize=4.5
                        )
            else:
                # Calculate median absolute deviation
                a = y_median
                c = 0.6745
                def statistic(a):
                    center = nanmedian
                    if callable(center):
                        axis = 0
                        center = np.apply_over_axes(center, a, axis)
                    return np.median((np.fabs(a-center))/c, axis=axis)
                y_median_error, _, = binned_statistic(y_nonans, x_nonans,
                                            statistic=statistic,
                                            bins=x_median)[:2]
                ax.plot(x_median,
                        y_median,
                        color='r',
                        ecolor='k',
                        marker='s',
                        linestyle='None',
                        label=label,
                        alpha=0.5,
                        markersize=4.5,
                        yerr=(y_median_error),
                        )

        if poly_fit:
            from scipy.optimize import curve_fit

            weights = np.abs(1.0 / x_error_nonans)
            weights /= np.sum(weights)
            def f(x, A): # this is your 'straight line' y=f(x)
                return A*x / weights

            b = x_nonans * weights
            A = np.array([y_nonans * weights,]).T
            p = [np.dot(np.linalg.pinv(A), b)[0], 0]

            #p, V = curve_fit(f, x_nonans, y_nonans, p0=0.15,)
            p = [p[0], 0]

            x_fit = np.linspace(-10, 100, 100)
            y_poly_fit = p[0] * x_fit + p[1]
            ax.plot(x_fit,
                    y_poly_fit,
                    color=c_cycle[1],
                    linestyle='dotted',
                    linewidth=2,
                    label=\
                        'Polynomial: \n' + \
                        'DGR = {0:.2f}'.format(p[0] * 100.0) + \
                        r' $\times\,10^{-22}$ cm$^{2}$ mag',
                    alpha=0.7,
                    )

        # Plot sensitivies
        #x_limit = np.median(x_errors[0])
        #ax.axvline(x_limit, color='k', linestyle='--')
        if 'dgr_error' in fit_params:
            dgr_error_text = r'$_{-%.1f}^{+%.1f}$ ' % fit_params['dgr_error']
        else:
            dgr_error_text = ''

        # Plot 1 to 1 pline
        y_fit = np.linspace(-10, 100, 1000)
        dgr_error_text = \
            r'$^{+%.1f}_{-%.1f}$ ' % (fit_params['dgr_error'][1] * 100.,
                                      fit_params['dgr_error'][0] * 100.)
        cloud_text = 'DGR:\n' + \
                     '{0:.1f}'.format(fit_params['dgr'] * 100.) + \
                dgr_error_text + \
                r' $\times\,10^{-22}$ cm$^{2}$ mag'

        x_fit = fit_params['dgr'] * y_fit
        label = cloud_text

        ax.plot(
                x_fit,
                y_fit,
                #color='#6B47B2',
                color=c_cycle[1],
                linestyle='--',
                linewidth=2,
                #alpha=0.8,
                label=label,
                )

        # Annotations
        #anno_xpos = 0.95

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_ylabel(r'$N($H$\textsc{i}) \times\,10^{20}$ cm$^{-2}$')
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
        plt.savefig(filename, bbox_inches='tight', dpi=100)

def plot_pdf_grid(av_list=None, nhi_list=None, nh2_list=None, dgr_list=None,
        names_list=None, limits=None, savedir='./', filename=None, show=True,
        scale=(0,0), n_bins=200, fit_gaussian=False, returnimage=False,
        title='', base=10.0, normalize=False, hi_trans_dict=None):

    ''' Plots a probability distribution function of an image.

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


    # create color cycle
    if hi_trans_dict is None:
        num_colors = 3
    else:
        num_colors = 3

    c_cycle = myplt.set_color_cycle(num_colors=num_colors,
                                    cmap_limits=[0, 0.8])

    for i in xrange(0, len(names_list)):
        cloud_name = names_list[i]
        av = av_list[i]
        nhi = nhi_list[i]
        nh2 = nh2_list[i]
        dgr = dgr_list[i]
        av_nhi = dgr * nhi
        av_nh2 = dgr * 2.0 * nh2
        #av_error = av_error_list[i]
        #nhi = nhi_list[i]
        #fit_params = fit_params_list[i]
        ax = axes[i]

        # Drop the NaNs from the images
        indices = np.where(av == av)
        av_nonans = av[indices]
        av_nhi_nonans = av_nhi[indices]
        av_nh2_nonans = av_nh2[indices]

        def hist(data):
            # Derive the histograms
            bin_edges = np.logspace(-3, 3, num=n_bins, base=base)
            n = np.zeros(n_bins - 1)
            for i in xrange(n_bins - 1):
                bin_count = len(data[(data > bin_edges[i]) & \
                                     (data < bin_edges[i + 1])])
                n[i] = bin_count

            bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])

            n = np.append(n, 0)
            bin_centers = np.append(bin_centers,
                    bin_centers[-1] + (bin_centers[-1] - bin_centers[-2]))

            # Normalize the bins
            if normalize:
                mean_loc = np.argmin(np.abs(bin_centers - data.mean()))
                bin_centers = np.log(bin_centers / bin_centers[mean_loc])

            return n, bin_centers

        n_av, bin_centers = hist(av_nonans)
        n_av_nhi, bin_centers = hist(av_nhi_nonans)
        n_av_nh2, bin_centers = hist(av_nh2_nonans)

        norm = np.sum(n_av)
        n_av /= norm
        n_av_nhi /= norm
        n_av_nh2 /= norm

        ax.errorbar(
            bin_centers,
            n_av,
            #yerr = n**0.5,
            marker = '',
            label=r'A$_V$',
            linewidth=1.5,
            color=c_cycle[0],
            drawstyle = 'steps-mid'
        )
        ax.errorbar(
            bin_centers,
            n_av_nhi,
            #yerr = n**0.5,
            marker = '',
            label=r'$N($H$\textsc{i}) \times$ DGR',
            #color = 'k',
            linewidth=1.5,
            color=c_cycle[1],
            drawstyle = 'steps-mid'
        )
        ax.errorbar(
            bin_centers,
            n_av_nh2,
            #yerr = n**0.5,
            label=r'$2\,N($H$_2) \times$ DGR',
            marker = '',
            linewidth=1.4,
            #color = 'k',
            color=c_cycle[2],
            drawstyle = 'steps-mid'
        )

        # Fit a gausssian to the distribution
        if fit_gaussian:
            def gauss(x, a, x0, sigma):
                return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

            indices = np.where((bin_centers > -1.5) & \
                                           (bin_centers < 1))

            bin_centers_crop, n_crop = bin_centers[indices], n[indices]

            popt, pcov = curve_fit(gauss,
                                   bin_centers_crop,
                                   n_crop,
                                   p0=[200, 0, 1],
                                   maxfev=1000000)
            ax.plot(bin_centers,
                    gauss(bin_centers, *popt),
                    color = 'r')

        if hi_trans_dict is not None:
            count = 0
            hi_trans_list = []
            hi_trans_error_list = []
            for core in hi_trans_dict:
                if hi_trans_dict[core]['cloud'] == cloud_name:
                    hi_trans = hi_trans_dict[core]['k09_transition']
                    hi_trans_error = \
                        np.array(hi_trans_dict[core]['k09_transition_error'])
                    nhi_trans = hi_trans * 1.25 * dgr
                    nhi_trans_error = hi_trans_error * 1.25 * dgr
                    #if count == 0:
                    #    label = r'H\,$\textsc{i}$-to-H$_2$ transition'
                    #else:
                    #    label = None
                    #count += 1

                    hi_trans_list.append(nhi_trans)
                    hi_trans_error_list.append(nhi_trans_error)

            label = r'H\,$\textsc{i}$-to-H$_2$ transition'

            ax.axvspan(np.min(hi_trans_list) - np.mean(hi_trans_error_list),
                       np.max(hi_trans_list) + np.mean(hi_trans_error_list),
                       alpha=0.4,
                       linewidth=0,
                       #color=c_cycle[3],
                       color='k',
                       edgecolor='none',
                       label=label,
                       )


        try:
            if scale[0] == 0:
                x_scale = 'linear'
            elif scale[0] == 1:
                x_scale = 'log'
            if scale[1] == 0:
                y_scale = 'linear'
            elif scale[1] == 1:
                y_scale = 'log'
        except IndexError('Scale must be tuple with 2 integer elements.'):
            pass

        ax.set_xscale(x_scale, nonposx = 'clip')
        ax.set_yscale(y_scale, nonposy = 'clip')

        if i == 0:
            ax.legend(loc='upper left')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        if normalize:
            if base == 10:
                ax.set_xlabel(r'log$_{10}$(A$_{\rm V}$ / $\bar{\rm A}_{\rm V}$ (mag))')
            elif base == np.e:
                ax.set_xlabel(r'ln(A$_{\rm V}$ / $\bar{\rm A}_{\rm V}$ (mag))',)
        else:
            ax.set_xlabel(r'A$_V$ [mag]')
        ax.set_ylabel(r'Frequency',)
        ax.set_title(title)
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
        plt.savefig(filename,bbox_inches='tight')
    if show:
        fig.show()

def plot_results(results_dict):
    # unpack results_dict
    dgr_background = results_dict['params_summary']['dgr_background']
    dgr_cloud = results_dict['params_summary']['dgr_cloud']
    dgr_cloud_error = results_dict['params_summary']['dgr_cloud_error']
    av_data = results_dict['data_products']['av_data_backsub']
    av_error_data = results_dict['data']['av_error_data']
    nhi_image = results_dict['data_products']['nhi']
    nhi_image_background = \
            results_dict['data_products']['nhi_image_background']
    plot_kwargs = results_dict['plot_kwargs']
    global_args = results_dict['global_args']
    boot_result = results_dict['boot_result']

    av_cloud = create_cloud_model(av_data,
                                  nhi_image_background,
                                  dgr_background,)

    if nhi_image_background is not None:
        av_background = create_background_model(av_data,
                                     nhi_image,
                                     dgr_cloud,)
        #nhi_total = nhi_boot + nhi_back_boot
        #nhi_total = np.hstack((nhi_boot, nhi_back_boot))
        #av_boot = np.hstack((av_cloud, av_background))
        #av_images = (av_boot, av_cloud, av_background)
        av_images = (av_cloud, av_background)
        #nhi_images = (nhi_total, nhi_boot, nhi_back_boot)
        nhi_images = (nhi_image, nhi_image_background)
    else:
        nhi_total = nhi_image
        av_images = (av_data,)
        nhi_images = (nhi_total,)

    filename = plot_kwargs['figure_dir'] + \
               'av_nhi/' + plot_kwargs['filename_base'] + \
               '_av_vs_nhi.png'

    plot_av_vs_nhi(av_images,
                   nhi_images,
                   av_error=av_error_data,
                   fit_params=results_dict['params_summary'],
                   contour_plot=plot_kwargs['av_nhi_contour'],
                   limits=plot_kwargs['av_nhi_limits'],
                   filename=filename,
                   )

    # Plot distribution
    if global_args['use_intercept'] or global_args['use_background']:
        filename = plot_kwargs['figure_dir'] + \
                   'bootstrap_dists/' + plot_kwargs['filename_base'] + \
                   '_backdgr_vs_clouddgr.png'
        #print('\n\tSaving bootstrap distributions to:\n\t' + filename)
        plot_bootstrap_dist(boot_result[0], boot_result[1],
                            filename=filename,
                            axis_labels=(r'Cloud DGR [10$^{-20}$ cm$^2$ mag]',
                                    r'Background DGR [10$^{-20}$ cm$^2$ mag]'),
                            levels=4)

        filename = plot_kwargs['figure_dir'] + \
                   'bootstrap_dists/' + plot_kwargs['filename_base'] + \
                   '_int_vs_clouddgr.png'

        try:
            plot_bootstrap_dist(boot_result[0], boot_result[2],
                                filename=filename,
                                axis_labels=(r'Cloud DGR [10$^{-20}$ cm$^2$ mag]',
                                             r'Intercept [mag]'),
                                levels=4)
        except ValueError:
            plot_bootstrap_dist(boot_result[0], boot_result[2],
                                filename=filename,
                                axis_labels=(r'Cloud DGR [10$^{-20}$ cm$^2$ mag]',
                                             r'Intercept [mag]'),
                                contour_plot=False,
                                levels=4)
    else:
        filename = plot_kwargs['figure_dir'] + \
                   'bootstrap_dists/' + plot_kwargs['filename_base'] + \
                   '_clouddgr.png'
        #print('\n\tSaving bootstrap distributions to:\n\t' + filename)
        plot_bootstrap_hist(boot_result[0],
                        filename=filename,
                        axis_label=r'Cloud DGR [10$^{-20}$ cm$^2$ mag]',
                        statistics=(dgr_cloud, dgr_cloud_error))

def print_av_error_stats(av, av_error):

    from scipy.stats import nanmedian

    error_above = nanmedian(av_error[av > 5])
    error_below = nanmedian(av_error[av <= 5])
    print('\n\tMedian Av error below 5 mag = {0:.2f} mag'.format(error_below))
    print('\n\tMedian Av error above 5 mag = {0:.2f} mag'.format(error_above))

def plot_hi_vs_h_grid(hsd_list, hisd_list, core_names=None, model_results=None,
        model_analysis=None, xlimits=None, ylimits=None, model_fits=None,
        scale=('linear', 'linear'), filename=None, show_params=False, levels=5,
        ncols=2, hsd_error_list=None, hisd_error_list=None):

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Determine size of figure and number of grids
    # --------------------------------------------
    if 0:
        n = int(np.ceil(len(core_names)**0.5))
        if n**2 - n > len(core_names):
            nrows = n - 1
            ncols = ncols
            y_scaling = 1.0 - 1.0 / n
        else:
            nrows, ncols = n, n
            y_scaling = 1.0

    if 1:
        n = len(core_names)
        nrows = (n + 1) / ncols
        if n > nrows * ncols:
            nrows += 1
        y_scaling = nrows / 2.0
        x_scaling = ncols / 2.0

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Color map
    cmap = plt.cm.copper

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 4)]
    font_scale = 9

    figsize = (3.6*x_scaling, 3.6*y_scaling)

    # Create figure instance
    fig = plt.figure(figsize=figsize)

    if n == 1:
        n = 2

    c_cycle = myplt.set_color_cycle(num_colors=2, cmap_limits=[0.5, 0.8])

    if xlimits is not None and ylimits is not None:
        aspect = (xlimits[1] - xlimits[0]) / (ylimits[1] - ylimits[0])
    else:
        aspect = False

    if 1:
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.ravel(axes)
    else:
        axes = AxesGrid(fig, (1,1,1),
                        nrows_ncols=(nrows, ncols),
                        ngrids=n,
                        #axes_pad=0.1,
                        axes_pad=0.3,
                        aspect=False,
                        label_mode='all',
                        share_all=False,
                        )


    # Cycle through lists
    # -------------------
    for i, core in enumerate(core_names):
        hi_sd = hisd_list[i]
        h_sd = hsd_list[i]

        # Load parameters
        alphaG = model_analysis[core]['sternberg_results']['alphaG']
        alphaG_error = \
                model_analysis[core]['sternberg_results']['alphaG_error']
        phi_cnm = model_analysis[core]['krumholz_results']['phi_cnm']
        phi_cnm_error = \
                model_analysis[core]['krumholz_results']['phi_cnm_error']
        h_sd_fit = model_analysis[core]['krumholz_results']['hsd_fit']


        # Drop the NaNs from the images
        indices = np.where((hi_sd == hi_sd) &\
                           (h_sd == h_sd))

        hi_sd_nonans = hi_sd[indices]
        h_sd_nonans = h_sd[indices]

        # Create plot
        ax = axes[i]

        #ax.set_xticks([0, 40, 80, 120])

        if 1:
            if xlimits is None:
                xmin = np.min(h_sd_nonans)
                xmax = np.max(h_sd_nonans)
                xscalar = 0.15 * xmax
                xlimits = [xmin - xscalar, xmax + xscalar]
            if ylimits is None:
                ymin = np.min(hi_sd_nonans)
                ymax = np.max(hi_sd_nonans)
                yscalar = 0.15 * ymax
                ylimits = [ymin - yscalar, ymax + yscalar]

            cmap = myplt.truncate_colormap(plt.cm.gray_r,
                                           minval=0.2,
                                           maxval=1)

            l1 = myplt.scatter_contour(h_sd_nonans.ravel(),
                                 hi_sd_nonans.ravel(),
                                 threshold=2,
                                 log_counts=0,
                                 levels=levels,
                                 ax=ax,
                                 histogram2d_args=dict(bins=20,
                                        range=(((xlimits)[0], xlimits[1]),
                                                (ylimits[0], ylimits[1]))),
                                 plot_args=dict(marker='o',
                                                linestyle='none',
                                                markeredgewidth=0,
                                                color='black',
                                                alpha=0.4,
                                                markersize=2.5,
                                                ),
                                 contour_args=dict(
                                                   cmap=cmap,
                                                   ),
                                 )

        if xlimits is not None:
            ax.set_xlim(xlimits[0], xlimits[1])
        if ylimits is not None:
            ax.set_ylim(ylimits[0], ylimits[1])
            ylimits = None

        if model_fits is not None:
            c_cycle = myplt.set_color_cycle(num_colors=2,
                                            cmap_limits=[0.4, 0.8])

            for model in ('krumholz', 'sternberg'):
                if 'krumholz' in model:
                    label = 'K+09'
                    color = c_cycle[0]
                    alpha = 0.5
                else:
                    label = 'S+14'
                    color = c_cycle[1]
                    alpha = 0.5
                h_sd_fit = model_fits[i][model + '_results']['h_sd']
                hi_sd_fit = model_fits[i][model + '_results']['hi_sd']
                l3 = ax.plot(h_sd_fit,
                             hi_sd_fit,
                             linestyle='-',
                             label=label,
                             color=color,
                             linewidth=1,
                             zorder=1000,
                             alpha=1
                             )
        elif 0:
            # get bootstrap results
            model = 'krumholz'
            model = 'sternberg'
            core_results = model_results['cores'][core]

            params = {}
            nfits = len(core_results['krumholz_results']['phi_cnm'])
            if nfits > 500:
                nfits = 1000
            alpha = 1 / float(nfits) * 10.0
            for j in xrange(nfits):
                params['phi_cnm'] = \
                    core_results['krumholz_results']['phi_cnm'][j]
                params['Z'] = core_results['krumholz_results']['Z'][j]
                params['phi_mol'] = \
                    core_results['krumholz_results']['phi_mol'][j]
                params['alphaG'] = \
                    core_results['sternberg_results']['alphaG'][j]
                params['Z'] = core_results['sternberg_results']['Z'][j]
                params['phi_g'] = \
                    core_results['sternberg_results']['phi_g'][j]
                if 'sternberg' in model:
                    model_fits = calc_sternberg(params,
                                              h_sd_extent=(0, 300),
                                              return_fractions=False,
                                              return_hisd=True)
                    color = 'b'
                elif 'krumholz' in model:
                    model_fits = calc_krumholz(params,
                                              h_sd_extent=(0, 300),
                                              return_fractions=False,
                                              return_hisd=True)
                    color = 'r'

                ax.plot(model_fits[1], model_fits[2],
                        linestyle='-',
                        color=color,
                        alpha=alpha,
                        )
                hi_trans = core_results['krumholz_results']['hi_transition'][j]

                if 0:
                    ax.axvline(hi_trans,
                               alpha=0.5,
                               color='r')

        elif 1:
            for model in ('krumholz', 'sternberg'):
                analysis = model_analysis[core][model + '_results']
                plot_fits = calc_model_plot_fit(analysis,
                                                model=model,
                                                )
                c_cycle = myplt.set_color_cycle(num_colors=2,
                                                cmap_limits=[0.4, 0.8])

                if 'krumholz' in model:
                    label = 'K+09'
                    color = c_cycle[0]
                    alpha = 0.5
                else:
                    label = 'S+14'
                    color = c_cycle[1]
                    alpha = 0.5

                fill_between = 0
                if 0:
                    l3 = ax.plot(plot_fits[0], plot_fits[1],
                            linestyle='-',
                            label=label,
                            color=color,
                            linewidth=1,
                            zorder=1000,
                            alpha=1
                            )
                else:
                    hi_sd_fit = \
                        model_analysis[core][model + '_results']['hisd_fit']
                    l3 = ax.plot(h_sd_fit, hi_sd_fit,
                            linestyle='-',
                            label=label,
                            color=color,
                            linewidth=1,
                            zorder=1000,
                            alpha=1
                            )
                if fill_between:
                    if sum(plot_fits[2] > plot_fits[3]) > 0:
                        where = plot_fits[2] > plot_fits[3]
                    else:
                        where = plot_fits[2] < plot_fits[3]

                    if model == 'sternberg':
                        ax.fill_between(plot_fits[0],
                                        plot_fits[2],
                                        plot_fits[3],
                                        where=where,
                                        facecolor=color,
                                        edgecolor='none',
                                        alpha=0.3,
                                        interpolate=True,
                                        zorder=0,
                                        )
                elif 0:
                    l3 = ax.plot(plot_fits[0], plot_fits[1],
                            linestyle='-',
                            label=label,
                            linewidth=2,
                            color=color,
                            alpha=1,
                            )
                    if model == 'sternberg':
                        ax.plot(plot_fits[0], plot_fits[2],
                                linestyle='-',
                                color=color,
                                alpha=alpha,
                                linewidth=1.5,
                                zorder=0,
                                )
                        ax.plot(plot_fits[0], plot_fits[3],
                                linestyle='-',
                                color=color,
                                alpha=alpha,
                                linewidth=1.5,
                                zorder=0,
                                )

        if hsd_error_list[i] is not None:
            print hsd_error_list[i]
            print hisd_error_list[i]
            ax.errorbar(90,
                        np.max(hi_sd_nonans) * 1.1,
                        xerr=np.median(hsd_error_list[i]),
                        yerr=np.median(hisd_error_list[i]),
                        #markersize=1.5,
                        marker='',
                        alpha=0.3,
                        color='k',
                        )

        if i == 0:
            ax.legend(loc='upper right')

        # Annotations
        anno_xpos = 0.95

        if show_params:
            alphaG_text = r'\noindent$\alpha G$ =' + \
                           r' %.2f' % (alphaG) + \
                           r'$^{+%.2f}_{-%.2f}$ \\' % (alphaG_error[0],
                                                       alphaG_error[1])
            phi_cnm_text = r'\noindent$\phi_{\rm CNM}$ =' + \
                           r' %.2f' % (phi_cnm) + \
                           r'$^{+%.2f}_{-%.2f}$ \\' % (phi_cnm_error[0],
                                                       phi_cnm_error[1])

            ax.annotate(alphaG_text + phi_cnm_text,
                    xytext=(anno_xpos, 0.05),
                    xy=(anno_xpos, 0.05),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=8,
                    color='k',
                    bbox=dict(boxstyle='round',
                              facecolor='w',
                              alpha=1),
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    )

        ax.annotate(core,
                    xytext=(0.9, 0.1),
                    xy=(0, 0),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=10,
                    color='k',
                    bbox=dict(boxstyle='square',
                              facecolor='w',
                              alpha=1),
                    verticalalignment='bottom',
                    horizontalalignment='right',
                    zorder=10000,
                    )

        ax.set_xscale(scale[0])
        ax.set_yscale(scale[1])

        # turn labels on or off
        if i % ncols == 0:
            ylabel = True
        else:
            #ax.yaxis.set_ticklabels([])
            ylabel = False

        if i >= len(core_names) - ncols:
            #ax.set_xlabel(labels[x_i])
            xlabel = True
        else:
            xlabel = False

        if len(core_names) % ncols > 0 and i == len(core_names) - 1:
            axes[i + 1].axis('off')

        # Adjust asthetics
        if xlabel:
            ax.set_xlabel(r'$\Sigma_{\rm H\,I}$ + $\Sigma_{\rm H_2}$ ' + \
                           '[M$_\odot$ pc$^{-2}$]',)
        if ylabel:
            ax.set_ylabel(r'$\Sigma_{\rm H\,I}$ [M$_\odot$ pc$^{-2}$]',)

        if 'log' not in scale:
            ax.locator_params(nbins=5)

    if 0:
        fig.legend(l2 + l3,
                   ('S+14', 'K+09'),
                   loc='upper center',
                   #loc=3,
                   ncol=2,
                   numpoints=6,
                   bbox_transform = plt.gcf().transFigure,
                   #bbox_transform = imagegrid.transFigure,
                   #bbox_to_anchor=(0., 0.95, 1.05, -0.0),
                   mode="expand",
                   bbox_to_anchor=(0.3, 1.0, 0.5, 0.1),
                   borderaxespad=0.0)

    if filename is not None:
        if filename[-3:] == 'pdf':
            dpi = 600
        else:
            dpi = 100
        plt.savefig(filename, bbox_inches='tight', dpi=dpi)

def plot_hi_cdf_grid(hsd_list, hisd_list, core_names=None,
        model_results=None, model_analysis=None, xlimits=None, ylimits=None,
        scale=('linear', 'linear'), filename=None, show_params=False,
        levels=5, ncols=2):

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Determine size of figure and number of grids
    # --------------------------------------------
    if 0:
        n = int(np.ceil(len(core_names)**0.5))
        if n**2 - n > len(core_names):
            nrows = n - 1
            ncols = ncols
            y_scaling = 1.0 - 1.0 / n
        else:
            nrows, ncols = n, n
            y_scaling = 1.0

    if 1:
        n = len(core_names)
        nrows = (n + 1) / ncols
        if n > nrows * ncols:
            nrows += 1
        y_scaling = nrows / 2.0
        x_scaling = ncols / 2.0

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Color map
    cmap = plt.cm.copper

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 4)]
    font_scale = 9

    figsize = (3.6*x_scaling, 3.6*y_scaling)

    # Create figure instance
    fig = plt.figure(figsize=figsize)

    if n == 1:
        n = 2

    c_cycle = myplt.set_color_cycle(num_colors=2, cmap_limits=[0.5, 0.8])

    if xlimits is not None and ylimits is not None:
        aspect = (xlimits[1] - xlimits[0]) / (ylimits[1] - ylimits[0])
    else:
        aspect = False

    if 1:
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.ravel(axes)
    else:
        axes = AxesGrid(fig, (1,1,1),
                        nrows_ncols=(nrows, ncols),
                        ngrids=n,
                        #axes_pad=0.1,
                        axes_pad=0.3,
                        aspect=False,
                        label_mode='all',
                        share_all=False,
                        )


    # Cycle through lists
    # -------------------
    for i, core in enumerate(core_names):
        hi_sd = hisd_list[i]
        h_sd = hsd_list[i]

        # Load parameters
        alphaG = model_analysis[core]['sternberg_results']['alphaG']
        alphaG_error = \
                model_analysis[core]['sternberg_results']['alphaG_error']
        phi_cnm = model_analysis[core]['krumholz_results']['phi_cnm']
        phi_cnm_error = \
                model_analysis[core]['krumholz_results']['phi_cnm_error']
        hi_sd_fit_sternberg = \
                model_analysis[core]['sternberg_results']['hisd_fit']
        hi_sd_fit_krumholz = \
                model_analysis[core]['krumholz_results']['hisd_fit']
        h_sd_fit = model_analysis[core]['krumholz_results']['hsd_fit']


        # Drop the NaNs from the images
        indices = np.where((hi_sd == hi_sd) &\
                           (h_sd == h_sd))

        hi_sd_nonans = hi_sd[indices]
        h_sd_nonans = h_sd[indices]

        # Create plot
        ax = axes[i]

        #ax.set_xticks([0, 40, 80, 120])

        if 1:
            kwargs = {'linewidth': 2}
            myplt.plot_cdf(hi_sd_nonans,
                           ax=ax,
                           plot_kwargs=kwargs,
                           )

        if xlimits is not None:
            ax.set_xlim(xlimits[0], xlimits[1])
        if ylimits is not None:
            ax.set_ylim(ylimits[0], ylimits[1])
            ylimits = None

        for model in ('sternberg', ):
            analysis = model_analysis[core][model + '_results']
            hi_trans = analysis['hi_transition']
            hi_trans_error = analysis['hi_transition_error']

            ax.axvline(hi_trans,
                       color='k',
                       alpha=1,
                       linestyle='--',
                       linewidth=2,
                       )
            ax.axvspan(hi_trans - hi_trans_error[0],
                       hi_trans + hi_trans_error[1],
                       color='k',
                       alpha=0.3,
                       )

        if i == 0:
            ax.legend(loc='upper left')

        # Annotations
        anno_xpos = 0.95

        ax.annotate(core,
                    xytext=(0.9, 0.1),
                    xy=(0, 0),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=10,
                    color='k',
                    bbox=dict(boxstyle='square',
                              facecolor='w',
                              alpha=1),
                    verticalalignment='bottom',
                    horizontalalignment='right',
                    )

        ax.set_xscale(scale[0])
        ax.set_yscale(scale[1])

        # turn labels on or off
        if i % ncols == 0:
            ylabel = True
        else:
            #ax.yaxis.set_ticklabels([])
            ylabel = False

        if i >= len(core_names) - ncols:
            #ax.set_xlabel(labels[x_i])
            xlabel = True
        else:
            xlabel = False

        if len(core_names) % ncols > 0 and i == len(core_names) - 1:
            axes[i + 1].axis('off')

        # Adjust asthetics
        if xlabel:
            ax.set_xlabel(r'$\Sigma_{\rm H\,I}$ ' + \
                           '[M$_\odot$ pc$^{-2}$]',)
        if ylabel:
            ax.set_ylabel(r'CDF',)

        if 'log' not in scale:
            ax.locator_params(nbins=5)

    if filename is not None:
        if filename[-3:] == 'pdf':
            dpi = 600
        else:
            dpi = 100
        plt.savefig(filename, bbox_inches='tight', dpi=dpi)

def plot_rh2_vs_h_grid(hsd_list, hisd_list, core_names=None, model_fits=None,
        model_results=None, model_analysis=None, xlimits=None, ylimits=None,
        scale=('linear', 'linear'), filename=None, show_params=False,
        levels=5, ncols=2):

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Determine size of figure and number of grids
    # --------------------------------------------
    if 0:
        n = int(np.ceil(len(core_names)**0.5))
        if n**2 - n > len(core_names):
            nrows = n - 1
            ncols = ncols
            y_scaling = 1.0 - 1.0 / n
        else:
            nrows, ncols = n, n
            y_scaling = 1.0

    if 1:
        n = len(core_names)
        nrows = (n + 1) / ncols
        if n > nrows * ncols:
            nrows += 1
        y_scaling = nrows / 2.0
        x_scaling = ncols / 2.0

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Color map
    cmap = plt.cm.copper

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 4)]
    font_scale = 9

    figsize = (3.6*x_scaling, 3.6*y_scaling)

    # Create figure instance
    fig = plt.figure(figsize=figsize)

    if n == 1:
        n = 2

    c_cycle = myplt.set_color_cycle(num_colors=2, cmap_limits=[0.5, 0.8])

    if xlimits is not None and ylimits is not None:
        aspect = (xlimits[1] - xlimits[0]) / (ylimits[1] - ylimits[0])
    else:
        aspect = False

    if 0:
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.ravel(axes)
    else:
        axes = AxesGrid(fig, (1,1,1),
                        nrows_ncols=(nrows, ncols),
                        ngrids=n,
                        #axes_pad=0.1,
                        axes_pad=0.1,
                        aspect=False,
                        label_mode='L',
                        share_all=False,
                        )


    # Cycle through lists
    # -------------------
    for i, core in enumerate(core_names):
        hi_sd = hisd_list[i]
        h_sd = hsd_list[i]

        # Load parameters
        alphaG = model_analysis[core]['sternberg_results']['alphaG']
        alphaG_error = \
                model_analysis[core]['sternberg_results']['alphaG_error']
        phi_cnm = model_analysis[core]['krumholz_results']['phi_cnm']
        phi_cnm_error = \
                model_analysis[core]['krumholz_results']['phi_cnm_error']
        hi_sd_fit_sternberg = \
                model_analysis[core]['sternberg_results']['hisd_fit']
        hi_sd_fit_krumholz = \
                model_analysis[core]['krumholz_results']['hisd_fit']
        h_sd_fit = model_analysis[core]['krumholz_results']['hsd_fit']


        # Drop the NaNs from the images
        indices = np.where((hi_sd == hi_sd) &\
                           (h_sd == h_sd))

        hi_sd_nonans = hi_sd[indices]
        h_sd_nonans = h_sd[indices]

        rh2 = (h_sd_nonans - hi_sd_nonans) / hi_sd_nonans

        # Create plot
        ax = axes[i]

        #ax.set_xticks([0, 40, 80, 120])

        if 1:
            if xlimits is None:
                xmin = np.min(h_sd_nonans)
                xmax = np.max(h_sd_nonans)
                xscalar = 0.15 * xmax
                xlimits = [xmin - xscalar, xmax + xscalar]
            if ylimits is None:
                ymin = np.min(rh2)
                ymax = np.max(rh2)
                yscalar = 0.15 * ymax
                ylimits = [ymin - yscalar, ymax + yscalar]

            cmap = myplt.truncate_colormap(plt.cm.gray_r,
                                           minval=0.2,
                                           maxval=1)

            ax.scatter(h_sd_nonans.ravel(),
                       rh2.ravel(),
                       #markersize=1.5,
                       s=5,
                       alpha=0.4,
                       color='k',
                       )

        if xlimits is not None:
            ax.set_xlim(xlimits[0], xlimits[1])
        if ylimits is not None:
            ax.set_ylim(ylimits[0], ylimits[1])
            #ylimits = None

        if model_fits is not None:
            c_cycle = myplt.set_color_cycle(num_colors=2,
                                            cmap_limits=[0.4, 0.8])

            for model in ('krumholz', 'sternberg'):
                if 'krumholz' in model:
                    label = 'K+09'
                    color = c_cycle[0]
                    alpha = 0.5
                else:
                    label = 'S+14'
                    color = c_cycle[1]
                    alpha = 0.5
                h_sd_fit = model_fits[i][model + '_results']['h_sd']
                rh2_fit = model_fits[i][model + '_results']['rh2']
                l3 = ax.plot(h_sd_fit,
                             rh2_fit,
                             linestyle='-',
                             label=label,
                             color=color,
                             linewidth=1,
                             zorder=1000,
                             alpha=1
                             )
        elif 0:
            # get bootstrap results
            model = 'krumholz'
            core_results = model_results['cores'][core]

            params = {}
            nfits = len(core_results['krumholz_results']['phi_cnm'])
            if nfits > 500:
                nfits = 500
            alpha = 1 / float(nfits) * 10.0
            for j in xrange(nfits):
                params['phi_cnm'] = \
                    core_results['krumholz_results']['phi_cnm'][j]
                params['Z'] = core_results['krumholz_results']['Z'][j]
                params['phi_mol'] = \
                    core_results['krumholz_results']['phi_mol'][j]
                if 'sternberg' in model:
                    model_fits = calc_sternberg(params,
                                              h_sd_extent=(0, limits[1]),
                                              return_fractions=False,
                                              return_hisd=True)
                elif 'krumholz' in model:
                    model_fits = calc_krumholz(params,
                                              h_sd_extent=(0, limits[1]),
                                              return_fractions=False,
                                              return_hisd=True)

                ax.plot(model_fits[1], model_fits[2],
                        linestyle='-',
                        color=c_cycle[2],
                        alpha=alpha,
                        )
                hi_trans = core_results['krumholz_results']['hi_transition'][j]
                ax.axvline(hi_trans,
                           alpha=0.5,
                           color='r')

        else:
            if 0:
                l2 = ax.plot(h_sd_fit, hi_sd_fit_sternberg,
                        label='S+14',
                        color=c_cycle[1],
                        alpha=0.75,
                        )
            for model in ('krumholz', 'sternberg'):
                analysis = model_analysis[core][model + '_results']
                plot_fits = calc_model_plot_fit(analysis,
                                                model=model)

                rh2_fit = (plot_fits[0] - plot_fits[1]) / plot_fits[1]

                if 'krumholz' in model:
                    label = 'K+09'
                    color = c_cycle[0]
                    alpha = 0.3
                else:
                    label = 'S+14'
                    color = c_cycle[1]
                    alpha = 0.8


                l3 = ax.plot(plot_fits[0], rh2_fit,
                        linestyle='-',
                        label=label,
                        color=color,
                        linewidth=1,
                        zorder=1000,
                        alpha=1
                        )

        if i == 0:
            ax.legend(loc='upper left')

        # Annotations
        anno_xpos = 0.95

        if show_params:
            alphaG_text = r'\noindent$\alpha G$ =' + \
                           r' %.2f' % (alphaG) + \
                           r'$^{+%.2f}_{-%.2f}$ \\' % (alphaG_error[0],
                                                       alphaG_error[1])
            phi_cnm_text = r'\noindent$\phi_{\rm CNM}$ =' + \
                           r' %.2f' % (phi_cnm) + \
                           r'$^{+%.2f}_{-%.2f}$ \\' % (phi_cnm_error[0],
                                                       phi_cnm_error[1])

            ax.annotate(alphaG_text + phi_cnm_text,
                    xytext=(anno_xpos, 0.05),
                    xy=(anno_xpos, 0.05),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=8,
                    color='k',
                    bbox=dict(boxstyle='round',
                              facecolor='w',
                              alpha=1),
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    )

        ax.annotate(core,
                    xytext=(0.9, 0.1),
                    xy=(0, 0),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=10,
                    color='k',
                    bbox=dict(boxstyle='square',
                              facecolor='w',
                              alpha=1),
                    verticalalignment='bottom',
                    horizontalalignment='right',
                    )

        ax.set_xscale(scale[0])
        ax.set_yscale(scale[1])

        # turn labels on or off
        if i % ncols == 0:
            ylabel = True
        else:
            #ax.yaxis.set_ticklabels([])
            ylabel = False

        if i >= len(core_names) - ncols:
            #ax.set_xlabel(labels[x_i])
            xlabel = True
        else:
            xlabel = False

        if len(core_names) % ncols > 0 and i == len(core_names) - 1:
            axes[i + 1].axis('off')

        # Adjust asthetics
        if xlabel:
            ax.set_xlabel(r'$\Sigma_{\rm H\,I}$ + $\Sigma_{\rm H_2}$ ' + \
                           '[M$_\odot$ pc$^{-2}$]',)
        if ylabel:
            ax.set_ylabel(r'$R_{H2}$',)

        if 'log' not in scale:
            ax.locator_params(nbins=5)

    if 0:
        fig.legend(l2 + l3,
                   ('S+14', 'K+09'),
                   loc='upper center',
                   #loc=3,
                   ncol=2,
                   numpoints=6,
                   bbox_transform = plt.gcf().transFigure,
                   #bbox_transform = imagegrid.transFigure,
                   #bbox_to_anchor=(0., 0.95, 1.05, -0.0),
                   mode="expand",
                   bbox_to_anchor=(0.3, 1.0, 0.5, 0.1),
                   borderaxespad=0.0)

    if filename is not None:
        if filename[-3:] == 'pdf':
            dpi = 600
        else:
            dpi = 100
        plt.savefig(filename, bbox_inches='tight', dpi=dpi)

def plot_rh2_vs_h(hsd_list, hisd_list, core_names=None,
        model_results=None, model_analysis=None, xlimits=None, ylimits=None,
        scale=('linear', 'linear'), filename=None, show_params=False,
        levels=5, ncols=2):

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid
    import matplotlib
    matplotlib.use('Agg')

    # Determine size of figure and number of grids
    # --------------------------------------------
    if 0:
        n = int(np.ceil(len(core_names)**0.5))
        if n**2 - n > len(core_names):
            nrows = n - 1
            ncols = ncols
            y_scaling = 1.0 - 1.0 / n
        else:
            nrows, ncols = n, n
            y_scaling = 1.0

    if 1:
        n = len(core_names)
        nrows = (n + 1) / ncols
        if n > nrows * ncols:
            nrows += 1
        y_scaling = nrows / 2.0
        x_scaling = ncols / 2.0

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Color map
    cmap = plt.cm.copper

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 4)]
    font_scale = 9

    figsize = (3.6*x_scaling, 3.6*y_scaling)

    # Create figure instance
    fig = plt.figure(figsize=figsize)

    if n == 1:
        n = 2

    c_cycle = myplt.set_color_cycle(num_colors=2, cmap_limits=[0.5, 0.8])

    if xlimits is not None and ylimits is not None:
        aspect = (xlimits[1] - xlimits[0]) / (ylimits[1] - ylimits[0])
    else:
        aspect = False

    if 0:
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.ravel(axes)
    else:
        axes = AxesGrid(fig, (1,1,1),
                        nrows_ncols=(nrows, ncols),
                        ngrids=n,
                        #axes_pad=0.1,
                        axes_pad=0.1,
                        aspect=False,
                        label_mode='L',
                        share_all=False,
                        )


    # Cycle through lists
    # -------------------
    for i, core in enumerate(core_names):
        hi_sd = hisd_list[i]
        h_sd = hsd_list[i]

        # Load parameters
        alphaG = model_analysis[core]['sternberg_results']['alphaG']
        alphaG_error = \
                model_analysis[core]['sternberg_results']['alphaG_error']
        phi_cnm = model_analysis[core]['krumholz_results']['phi_cnm']
        phi_cnm_error = \
                model_analysis[core]['krumholz_results']['phi_cnm_error']
        hi_sd_fit_sternberg = \
                model_analysis[core]['sternberg_results']['hisd_fit']
        hi_sd_fit_krumholz = \
                model_analysis[core]['krumholz_results']['hisd_fit']
        h_sd_fit = model_analysis[core]['krumholz_results']['hsd_fit']


        # Drop the NaNs from the images
        indices = np.where((hi_sd == hi_sd) &\
                           (h_sd == h_sd))

        hi_sd_nonans = hi_sd[indices]
        h_sd_nonans = h_sd[indices]

        rh2 = (h_sd_nonans - hi_sd_nonans) / hi_sd_nonans

        # Create plot
        ax = axes[i]

        #ax.set_xticks([0, 40, 80, 120])

        if 1:
            if xlimits is None:
                xmin = np.min(h_sd_nonans)
                xmax = np.max(h_sd_nonans)
                xscalar = 0.15 * xmax
                xlimits = [xmin - xscalar, xmax + xscalar]
            if ylimits is None:
                ymin = np.min(rh2)
                ymax = np.max(rh2)
                yscalar = 0.15 * ymax
                ylimits = [ymin - yscalar, ymax + yscalar]

            cmap = myplt.truncate_colormap(plt.cm.gray_r,
                                           minval=0.2,
                                           maxval=1)

            ax.scatter(h_sd_nonans.ravel(),
                       rh2.ravel(),
                       #markersize=1.5,
                       marker='o',
                       s=3,
                       alpha=0.3,
                       color='k',
                       )
        elif 0:
            if xlimits is None:
                xmin = np.min(h_sd_nonans)
                xmax = np.max(h_sd_nonans)
                xscalar = 0.15 * xmax
                xlimits = [xmin - xscalar, xmax + xscalar]
            if ylimits is None:
                ymin = np.min(rh2)
                ymax = np.max(rh2)
                yscalar = 0.15 * ymax
                ylimits = [ymin - yscalar, ymax + yscalar]

            cmap = myplt.truncate_colormap(plt.cm.gray_r,
                                           minval=0.2,
                                           maxval=1)

            bins = (np.linspace(xlimits[0], xlimits[1], 20),
                    np.logspace(np.log10(ylimits[0]), np.log10(ylimits[1]), 20))

            hist_range=((xlimits[0], xlimits[1]),
                        (np.log10(ylimits[0]), np.log10(ylimits[1])))
            print bins

            l1 = myplt.scatter_contour(h_sd_nonans.ravel(),
                                 rh2.ravel(),
                                 threshold=2,
                                 log_counts=0,
                                 levels=levels,
                                 ax=ax,
                                 histogram2d_args=dict(bins=bins,
                                                       range=hist_range,),
                                 plot_args=dict(marker='o',
                                                linestyle='none',
                                                markeredgewidth=0,
                                                color='black',
                                                alpha=0.4,
                                                markersize=2.5,
                                                ),
                                 contour_args=dict(
                                                   cmap=cmap,
                                                   ),
                                 )

        if xlimits is not None:
            ax.set_xlim(xlimits[0], xlimits[1])
        if ylimits is not None:
            ax.set_ylim(ylimits[0], ylimits[1])
            ylimits = None


        if xlimits is not None:
            ax.set_xlim(xlimits[0], xlimits[1])
        if ylimits is not None:
            ax.set_ylim(ylimits[0], ylimits[1])
            #ylimits = None

        if 0:
            # get bootstrap results
            model = 'krumholz'
            core_results = model_results['cores'][core]

            params = {}
            nfits = len(core_results['krumholz_results']['phi_cnm'])
            if nfits > 500:
                nfits = 500
            alpha = 1 / float(nfits) * 10.0
            for j in xrange(nfits):
                params['phi_cnm'] = \
                    core_results['krumholz_results']['phi_cnm'][j]
                params['Z'] = core_results['krumholz_results']['Z'][j]
                params['phi_mol'] = \
                    core_results['krumholz_results']['phi_mol'][j]
                if 'sternberg' in model:
                    model_fits = calc_sternberg(params,
                                              h_sd_extent=(0, limits[1]),
                                              return_fractions=False,
                                              return_hisd=True)
                elif 'krumholz' in model:
                    model_fits = calc_krumholz(params,
                                              h_sd_extent=(0, limits[1]),
                                              return_fractions=False,
                                              return_hisd=True)

                ax.plot(model_fits[1], model_fits[2],
                        linestyle='-',
                        color=c_cycle[2],
                        alpha=alpha,
                        )
                hi_trans = core_results['krumholz_results']['hi_transition'][j]
                ax.axvline(hi_trans,
                           alpha=0.5,
                           color='r')

        else:
            if 0:
                l2 = ax.plot(h_sd_fit, hi_sd_fit_sternberg,
                        label='S+14',
                        color=c_cycle[1],
                        alpha=0.75,
                        )
            for model in ('krumholz',):
                analysis = model_analysis[core][model + '_results']
                hsd = np.logspace(np.log10(0.1), np.log10(200), 1000)
                plot_fits = calc_model_plot_fit(analysis,
                                                model=model,
                                                hsd=hsd)

                rh2_fit = (plot_fits[0] - plot_fits[1]) / plot_fits[1]

                if 'krumholz' in model:
                    label = 'K+09'
                    color = c_cycle[1]
                    alpha = 0.6
                else:
                    label = 'S+14'
                    color = c_cycle[1]
                    alpha = 0.8


                l3 = ax.plot(plot_fits[0], rh2_fit,
                        linestyle='-',
                        label=label,
                        color='b',
                        linewidth=1,
                        zorder=1000,
                        alpha=0.6
                        )

        if 0:
            if i == 0:
                ax.legend(loc='upper left')

        # Annotations
        anno_xpos = 0.95

        if show_params:
            alphaG_text = r'\noindent$\alpha G$ =' + \
                           r' %.2f' % (alphaG) + \
                           r'$^{+%.2f}_{-%.2f}$ \\' % (alphaG_error[0],
                                                       alphaG_error[1])
            phi_cnm_text = r'\noindent$\phi_{\rm CNM}$ =' + \
                           r' %.2f' % (phi_cnm) + \
                           r'$^{+%.2f}_{-%.2f}$ \\' % (phi_cnm_error[0],
                                                       phi_cnm_error[1])

            ax.annotate(alphaG_text + phi_cnm_text,
                    xytext=(anno_xpos, 0.05),
                    xy=(anno_xpos, 0.05),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=8,
                    color='k',
                    bbox=dict(boxstyle='round',
                              facecolor='w',
                              alpha=1),
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    )

        if core == 'G169.32-16.17':
            core_name = 'B213'
        else:
            core_name = core
        ax.annotate(core_name,
                    xytext=(0.9, 0.1),
                    xy=(0, 0),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=10,
                    color='k',
                    bbox=dict(boxstyle='square',
                              facecolor='w',
                              alpha=1),
                    verticalalignment='bottom',
                    horizontalalignment='right',
                    )

        ax.set_xscale(scale[0])
        ax.set_yscale(scale[1])

        # turn labels on or off
        if i % ncols == 0:
            ylabel = True
        else:
            #ax.yaxis.set_ticklabels([])
            ylabel = False

        if i >= len(core_names) - ncols:
            #ax.set_xlabel(labels[x_i])
            xlabel = True
        else:
            xlabel = False

        if len(core_names) % ncols > 0 and i == len(core_names) - 1:
            axes[i + 1].axis('off')

        # Adjust asthetics
        if xlabel:
            ax.set_xlabel(r'$\Sigma_{\rm H\,I}$ + $\Sigma_{\rm H_2}$ ' + \
                           '[M$_\odot$ pc$^{-2}$]',)
        if ylabel:
            ax.set_ylabel(r'$R_{H2}$',)

        if 'log' not in scale:
            ax.locator_params(nbins=5)

    if 0:
        fig.legend(l2 + l3,
                   ('S+14', 'K+09'),
                   loc='upper center',
                   #loc=3,
                   ncol=2,
                   numpoints=6,
                   bbox_transform = plt.gcf().transFigure,
                   #bbox_transform = imagegrid.transFigure,
                   #bbox_to_anchor=(0., 0.95, 1.05, -0.0),
                   mode="expand",
                   bbox_to_anchor=(0.3, 1.0, 0.5, 0.1),
                   borderaxespad=0.0)

    if filename is not None:
        if filename[-3:] == 'pdf':
            dpi = 600
        else:
            dpi = 600
        #plt.savefig(filename, bbox_inches='tight', dpi=dpi)
        plt.savefig(filename, dpi=dpi)

def plot_rh2_vs_h_diagnostic(h_sd, rh2, h_sd_error=None, rh2_error=None,
        xlimits=None, ylimits=None, scale=('log', 'log'), filename=None,
        show_params=False, levels=5, ncols=2, scatter=True,
        model_results=None):

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid
    import matplotlib
    matplotlib.use('Agg')

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Color map
    cmap = plt.cm.copper

    # Color cycle, grabs colors from cmap
    font_scale = 9

    figsize = (3.6, 3.6)

    # Create figure instance
    fig = plt.figure(figsize=figsize)

    c_cycle = myplt.set_color_cycle(num_colors=3, cmap_limits=[0.5, 0.8])

    if 0:
        fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
        axes = np.ravel(axes)
    else:
        axes = AxesGrid(fig, (1,1,1),
                        nrows_ncols=(1, 1),
                        ngrids=1,
                        #axes_pad=0.1,
                        axes_pad=0.1,
                        aspect=False,
                        label_mode='L',
                        share_all=False,
                        )

    # Create plot
    ax = axes[0]

    #ax.set_xticks([0, 40, 80, 120])

    if 1:
        cmap = myplt.truncate_colormap(plt.cm.gray_r,
                                       minval=0.2,
                                       maxval=1)

        ax.errorbar(h_sd,
                    rh2,
                    xerr=h_sd_error,
                    yerr=rh2_error,
                    #markersize=1.5,
                    marker='^',
                    alpha=0.3,
                    color='k',
                    )

    if xlimits is not None:
        ax.set_xlim(xlimits[0], xlimits[1])
    if ylimits is not None:
        ax.set_ylim(ylimits[0], ylimits[1])

    if 1:
        alpha = 0.6
        s14 = model_results['sternberg_results']
        s14_params = (s14['alphaG'], s14['Z'], s14['phi_g'])
        k09 = model_results['krumholz_results']
        k09_params = (k09['phi_cnm'], k09['Z'], k09['phi_mol'])
        h_sd_extent = np.logspace(-2, 3, 1000)

        # Sternberg
        model_fits = calc_sternberg(s14,
                                  h_sd=h_sd_extent,
                                  return_fractions=False,
                                  return_hisd=True)

        if len(model_fits[0]) > 1:
            ax.plot(model_fits[1], model_fits[0],
                    linestyle='-',
                    label='S+14',
                    color='r',
                    alpha=alpha,
                    )

        model_fits = calc_krumholz(k09,
                                  h_sd=h_sd_extent,
                                  return_fractions=False,
                                  return_hisd=True)

        if len(model_fits[0]) > 1:
            ax.plot(model_fits[1], model_fits[0],
                    linestyle='-',
                    label='K+09',
                    color='b',
                    alpha=alpha,
                    )

    ax.legend(loc='best')

    ax.set_xscale(scale[0], nonposx='clip')
    ax.set_yscale(scale[1], nonposy='clip')

    # Adjust asthetics
    ax.set_xlabel(r'$\Sigma_{\rm H\,I}$ + $\Sigma_{\rm H_2}$ ' + \
                   '[M$_\odot$ pc$^{-2}$]',)
    ax.set_ylabel(r'$R_{H2}$',)

    if 'log' not in scale:
        ax.locator_params(nbins=5)

    if filename is not None:
        if filename[-3:] == 'pdf':
            dpi = 100
        else:
            dpi = 100
        #plt.savefig(filename, bbox_inches='tight', dpi=dpi)
        plt.savefig(filename, dpi=dpi)

def calc_model_plot_fit(analysis, model='krumholz', hsd=None,):

    if 'sternberg' in model:
        alphaG = analysis['alphaG']
        alphaG_low = alphaG - analysis['alphaG_error'][0]
        alphaG_high = alphaG + analysis['alphaG_error'][1]

        params = {'alphaG': alphaG,
                  'phi_g': analysis['phi_g'],
                  'Z': analysis['Z'],
                  }

        h_sd, hi_sd = calc_sternberg(params,
                                  h_sd_extent=(0, 100),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[1:3]

        params = {'alphaG': alphaG_low,
                  'phi_g': analysis['phi_g'],
                  'Z': analysis['Z'],
                  }
        hi_sd_low = calc_sternberg(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[2]

        params = {'alphaG': alphaG_high,
                  'phi_g': analysis['phi_g'],
                  'Z': analysis['Z'],
                  }
        hi_sd_high = calc_sternberg(params,
                                  h_sd_extent=(0, 100),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[2]
    elif 'krumholz' in model:
        phi_cnm = analysis['phi_cnm']
        phi_cnm_low = phi_cnm - analysis['phi_cnm_error'][0]
        phi_cnm_high = phi_cnm + analysis['phi_cnm_error'][1]

        params = {'phi_cnm': phi_cnm,
                  'phi_mol': analysis['phi_mol'],
                  'Z': analysis['Z'],
                  }
        h_sd, hi_sd = calc_krumholz(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[1:3]

        params = {'phi_cnm': phi_cnm_high,
                  'phi_mol': analysis['phi_mol'],
                  'Z': analysis['Z'],
                  }
        hi_sd_low = calc_krumholz(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[2]

        params = {'phi_cnm': phi_cnm_low,
                  'phi_mol': analysis['phi_mol'],
                  'Z': analysis['Z'],
                  }
        hi_sd_high = calc_krumholz(params,
                                  h_sd_extent=(0, 200),
                                  h_sd=hsd,
                                  return_fractions=False,
                                  return_hisd=True)[2]

    return h_sd, hi_sd, hi_sd_low, hi_sd_high

def plot_multicloud_results(results):

    print('\nPlotting multicloud results...')

    # Collect Data
    # =========================================================================

    spectra_list = []
    hi_range_kwargs_list = []
    av_list = []
    av_error_list = []
    nhi_list = []
    nh2_list = []
    rh2_list = []
    hsd_list = []
    hisd_list = []
    nhi_error_list = []
    nh2_error_list = []
    rh2_error_list = []
    hsd_error_list = []
    hisd_error_list = []
    hsd_median_error_list = []
    hisd_median_error_list = []
    hisd_cores_list = []
    hsd_cores_list = []
    rh2_cores_list = []
    hisd_error_cores_list = []
    hsd_error_cores_list = []
    hisd_median_error_cores_list = []
    hsd_median_error_cores_list = []
    rh2_error_cores_list = []
    model_results_list = []
    model_analysis_list = []
    model_analysis_dict = {}
    dgr_list = []
    fit_params_list = []
    cloud_name_list = []
    core_names_list = []
    core_list = []
    cloud_model_fits_list = []
    for i, cloud_name in enumerate(results):
        results_dict = results[cloud_name]
        figure_dir = results_dict['filenames']['figure_dir']
        results_dir = results_dict['filenames']['results_dir']
        plot_kwargs = results_dict['plot_kwargs']
        data_products = results_dict['data_products']
        spectra_list.append((data_products['hi_spectrum'],
                             data_products['hi_std_spectrum'],
                             data_products['co_spectrum'],)
                             )
        cloud_name_list.append(cloud_name)
        hi_range_kwargs_list.append(data_products['hi_range_kwargs'])
        #av_list.append(data_products['av_data_backsub'])
        av_list.append(data_products['av'])
        av_error_list.append(results_dict['data']['av_error_data'])
        nhi_list.append(data_products['nhi'])
        nh2_list.append(data_products['nh2'])
        rh2_list.append(data_products['rh2'])
        hsd_list.append(data_products['h_sd'])
        hisd_list.append(data_products['hi_sd'])
        nhi_error_list.append(data_products['nhi_error'])
        nh2_error_list.append(data_products['nh2_error'])
        rh2_error_list.append(data_products['rh2_error'])
        hsd_error_list.append(data_products['h_sd_error'])
        hisd_error_list.append(data_products['hi_sd_error'])
        hsd_median_error_list.append(data_products['h_sd_median_error'])
        hisd_median_error_list.append(data_products['hi_sd_median_error'])
        dgr_list.append(results_dict['mc_analysis']['dgr'])
        #fit_params_list.append(results_dict['params_summary'])
        model_results_list.append(results_dict['mc_results']['ss_model_results'])
        model_analysis_list.append(results_dict['mc_analysis'])
        model_analysis_dict[cloud_name] = results_dict['mc_analysis']

        # Modeling and cores
        global_args = results_dict['global_args']
        cores = global_args['ss_model_kwargs']['cores']
        cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']
        model_kwargs = global_args['ss_model_kwargs']['model_kwargs']
        rh2_core_list = []
        hsd_core_list = []
        hisd_core_list = []
        rh2_error_core_list = []
        hsd_error_core_list = []
        hsd_median_error_core_list = []
        hisd_error_core_list = []
        hisd_median_error_core_list = []
        core_names = []
        model_fits_list = []
        core_list.append(cores)
        for j, core in enumerate(cores_to_plot):
            core_indices = cores[core]['indices_orig']
            core_names.append(core)
            hisd_core_list.append(hisd_list[i][core_indices])
            hsd_core_list.append(hsd_list[i][core_indices])
            rh2_core_list.append(rh2_list[i][core_indices])
            hisd_error_core_list.append(hisd_error_list[i][core_indices])
            hsd_error_core_list.append(hsd_error_list[i][core_indices])
            hisd_median_error_core_list.append(hisd_median_error_list[i])
            hsd_median_error_core_list.append(hsd_median_error_list[i])
            rh2_error_core_list.append(rh2_error_list[i][core_indices])

            model_fits_list.append(refit_data(hsd_core_list[j],
                                              rh2_core_list[j],
                                              h_sd_error=hsd_error_core_list[j],
                                              rh2_error=rh2_error_core_list[j],
                                              model_kwargs=model_kwargs,
                                              )
                                   )

            if 0:
                rh2_copy = rh2_list[i].copy()
                rh2_copy[core_indices] = 1000
                plt.imshow(rh2_copy, origin='lower')
                plt.savefig('/d/bip3/ezbc/scratch/core_' + core + '.png')

        cloud_model_fits_list.append(model_fits_list)
        core_names_list.append(core_names)
        hisd_cores_list.append(hisd_core_list)
        hsd_cores_list.append(hsd_core_list)
        hisd_error_cores_list.append(hisd_error_core_list)
        hsd_error_cores_list.append(hsd_error_core_list)
        hisd_median_error_cores_list.append(hisd_median_error_core_list)
        hsd_median_error_cores_list.append(hsd_median_error_core_list)

        if 0:
            print(cloud_name)
            print('vel range', data_products['hi_range_kwargs']['vel_range'])
            print(results_dict['global_args']['vel_range_error'])



    # Print results
    # =========================================================================
    # Write results to a
    #print_av_error_stats(av_list[0], av_error_list[0])

    filename = results_dir + 'tables/multicloud_model_params.tex'
    write_model_params_table(model_analysis_dict,
                             filename,
                             models=('krumholz','sternberg'))

    # Write param summary to dataframe for ease of use
    filename = results_dir + 'tables/multicloud_model_params.csv'
    write_param_csv(model_analysis_dict,
                    core_list,
                    cloud_name_list,
                    filename,
                    )

    # Write nhi properties
    filename = results_dir + 'tables/nhi_properties.csv'
    write_nhi_properties_csv(nhi_list,
                             cloud_name_list,
                             filename,
                             )

    # Write param summary to dataframe for ease of use
    filename = results_dir + 'tables/multicloud_model_summary.pickle'
    write_fit_summary_dict(model_analysis_dict,
                           core_list,
                           cloud_name_list,
                           filename,
                           )

    # Write table for
    filename = results_dir + 'tables/multicloud_hi_transitions.csv'
    hi_trans_dict = collect_hi_transition_results(model_analysis_list,
                                                  cloud_name_list,
                                                  filename=filename,)

    # Plot the results
    # =========================================================================

    # Plot HI spectra
    # -------------------------------------------------------------------------
    hi_vel_axis = results_dict['data']['hi_vel_axis']
    co_vel_axis = results_dict['data']['co_vel_axis']


    filetypes = ['png', 'pdf']
    for filetype in filetypes:
        filename = plot_kwargs['figure_dir'] + \
                   'spectra/multicloud_spectra.' + filetype

        plot_spectra_grid(spectra_list,
                     hi_range_kwargs_list=hi_range_kwargs_list,
                     names_list=cloud_name_list,
                     hi_vel_axis=hi_vel_axis,
                     co_vel_axis=co_vel_axis,
                     filename=filename,
                     limits=[-30, 30, -10, 59],
                     )

        filename = results_dir + 'tables/multicloud_vel_ranges.tex'
        write_hi_vel_range_table(cloud_name_list,
                                 hi_range_kwargs_list,
                                 filename)

        # Plot N(HI) vs. Av
        # ---------------------------------------------------------------------
        filename = plot_kwargs['figure_dir'] + \
                   'av_nhi/multicloud_av_vs_nhi.' + filetype
        plot_av_vs_nhi_grid(av_list,
                            nhi_list,
                            av_error_list=av_error_list,
                            fit_params_list=model_analysis_list,
                            names_list=cloud_name_list,
                            filename=filename,
                            levels=(0.99, 0.98, 0.95, 0.86, 0.59),
                            poly_fit=False,
                            limits=[-2, 19, 0, 20]
                            )

        # Plot Av PDF
        # ---------------------------------------------------------------------
        filename = plot_kwargs['figure_dir'] + \
                   'pdfs/multicloud_pdfs.' + filetype
        plot_pdf_grid(av_list,
                  nhi_list=nhi_list,
                  nh2_list=nh2_list,
                  dgr_list=dgr_list,
                  hi_trans_dict=hi_trans_dict,
                  #limits = [-4,3,1,10000],
                  names_list=cloud_name_list,
                  #limits = [0.07,14,7,6000],
                  limits=[10**-2, 3 * 10**1, 5 * 10**-4, 3 * 10**-1],
                  scale=(1,1),
                  filename=filename,
                  #core_names=core_name_list,
                  show=False)

        # Plot HI vs H
        # ---------------------------------------------------------------------
        levels = (0.9, 0.8, 0.6, 0.3)
        for i, cloud in enumerate(cloud_name_list):
            core_names = core_names_list[i]
            print('\n\tPlotting Models')
            if len(core_names) > 10:
                ncols = 5
                if 0:
                    filename = plot_kwargs['figure_dir'] + \
                               'models/' + cloud + '_hisd_vs_hsd_1.' + filetype
                    plot_hi_vs_h_grid(hsd_cores_list[i][:10],
                                      hisd_cores_list[i][:10],
                                      core_names=core_names[:10],
                                      model_results=model_results_list[i],
                                      model_analysis=\
                                              model_analysis_list[i]['cores'],
                                      model_fits=model_fits_list,
                                      limits=[-9, 159, -1.5, 15],
                                      scale=('linear', 'linear'),
                                      levels=levels,
                                      filename=filename,
                                      )
                    filename = plot_kwargs['figure_dir'] + \
                               'models/' + cloud + '_hisd_vs_hsd_2.' + filetype
                    plot_hi_vs_h_grid(hsd_cores_list[i][10:],
                                      hisd_cores_list[i][10:],
                                      core_names=core_names[10:],
                                      model_results=model_results_list[i],
                                      model_analysis=\
                                          model_analysis_list[i]['cores'],
                                      model_fits=cloud_model_fits_list,
                                      limits=[-9, 159, -1.5, 15],
                                      levels=levels,
                                      scale=('linear', 'linear'),
                                      filename=filename,
                                      )
            else:
                ncols = 2


            # RH2 vs. H SD for L1478
            # -----------------------------------------------------------------
            if cloud == 'taurus':
                single_core = 'G169.32-16.17'
                filename = plot_kwargs['figure_dir'] + \
                           'models/' + cloud + '_rh2_vs_hsd_' + \
                           single_core + '.' + filetype
                index = core_names.index(single_core)
                plot_rh2_vs_h((hsd_cores_list[i][index],),
                              (hisd_cores_list[i][index],),
                              core_names=(core_names[index],),
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              ylimits=[10**-3, 10**2],
                              levels=levels,
                              #scale=('log', 'linear'),
                              scale=('linear', 'log'),
                              filename=filename,
                              ncols=ncols
                              )

                filename = plot_kwargs['figure_dir'] + \
                           'models/' + cloud + '_rh2_vs_hsd_' + \
                           single_core + '.ps'
                plot_rh2_vs_h((hsd_cores_list[i][index],),
                              (hisd_cores_list[i][index],),
                              core_names=(core_names[index],),
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              ylimits=[10**-3, 10**2],
                              levels=levels,
                              #scale=('log', 'linear'),
                              scale=('linear', 'log'),
                              filename=filename,
                              ncols=ncols
                              )

            # HI SD vs. H SD
            # -----------------------------------------------------------------
            filename = plot_kwargs['figure_dir'] + \
                       'models/' + cloud + '_hisd_vs_hsd.' + filetype
            plot_hi_vs_h_grid(hsd_cores_list[i],
                              hisd_cores_list[i],
                              core_names=core_names,
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              model_fits=cloud_model_fits_list[i],
                              hsd_error_list=hsd_median_error_cores_list[i],
                              hisd_error_list=hisd_median_error_cores_list[i],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              levels=levels,
                              #scale=('log', 'linear'),
                              #scale=('linear', 'log'),
                              scale=('linear', 'linear'),
                              filename=filename,
                              ncols=ncols
                              )
            # HI CDF
            # -----------------------------------------------------------------
            filename = plot_kwargs['figure_dir'] + \
                       'models/' + cloud + '_hisd_cdf.' + filetype
            plot_hi_cdf_grid(hsd_cores_list[i],
                              hisd_cores_list[i],
                              core_names=core_names,
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              #xlimits=[-1, 20],
                              levels=levels,
                              #scale=('log', 'linear'),
                              #scale=('linear', 'log'),
                              scale=('linear', 'linear'),
                              filename=filename,
                              ncols=ncols
                              )

            # RH2 vs. H SD
            # -----------------------------------------------------------------
            filename = plot_kwargs['figure_dir'] + \
                       'models/' + cloud + '_rh2_vs_hsd.' + filetype
            plot_rh2_vs_h_grid(hsd_cores_list[i],
                              hisd_cores_list[i],
                              core_names=core_names,
                              model_results=model_results_list[i],
                              model_analysis=\
                                  model_analysis_list[i]['cores'],
                              model_fits=cloud_model_fits_list[i],
                              #limits=[-9, 100, 2, 14],
                              #limits=[-9, 159, 3, 14],
                              xlimits=[-9, 100],
                              ylimits=[10**-3, 10**2],
                              levels=levels,
                              #scale=('log', 'linear'),
                              scale=('linear', 'log'),
                              filename=filename,
                              ncols=ncols
                              )

        if 0:
            filename = plot_kwargs['figure_dir'] + \
                       'av_nhi/multicloud_av_vs_nhi_log.' + filetype
            plot_av_vs_nhi_grid(av_list,
                                nhi_list,
                                av_error_list=av_error_list,
                                fit_params_list=fit_params_list,
                                names_list=cloud_name_list,
                                filename=filename,
                                levels=(0.99, 0.98, 0.95, 0.86, 0.59),
                                poly_fit=True,
                                scale=['log','log'],
                                #limits=[2, 20, -13, 20]
                                limits=[0.1, 100, 0.01, 100]
                                )

            filename = plot_kwargs['figure_dir'] + \
                       'av_nhi/multicloud_av_vs_nhi_scatter.' + filetype
            plot_av_vs_nhi_grid(av_list,
                                nhi_list,
                                av_error_list=av_error_list,
                                fit_params_list=fit_params_list,
                                names_list=cloud_name_list,
                                filename=filename,
                                levels=(0.99, 0.98, 0.95, 0.86, 0.59),
                                poly_fit=True,
                                contour=False,
                                limits=[2, 20, -6, 10]
                                )

