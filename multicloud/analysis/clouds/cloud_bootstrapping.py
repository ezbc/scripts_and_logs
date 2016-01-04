#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import myimage_analysis as myia
from multiprocessing.queues import Queue
import mygeometry as myg
import scipy


global debugging
debugging = True
#debugging = False


'''
Plotting
'''

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

    myplt.set_color_cycle(num_colors=4, cmap_limits=[0, 0.6])

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
            dgr_error_text = r'$^{+%.1f}_{-%.1f}$ ' % fit_params['dgr_error']
        else:
            dgr_error_text = ''

        # Plot 1 to 1 pline
        y_fit = np.linspace(-10, 100, 1000)
        dgr_error_text = \
            r'$^{+%.1f}_{-%.1f}$ ' % (fit_params['dgr_error'][0] * 100.,
                                      fit_params['dgr_error'][1] * 100.)
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

def plot_hi_vs_h_grid(hsd_list, hisd_list, core_names=None,
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

        elif 1:
            for model in ('krumholz', 'sternberg'):
                analysis = model_analysis[core][model + '_results']
                plot_fits = calc_model_plot_fit(analysis,
                                                model=model)
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
                l3 = ax.plot(plot_fits[0], plot_fits[1],
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

def plot_rh2_vs_h_grid(hsd_list, hisd_list, core_names=None,
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
    hisd_cores_list = []
    hsd_cores_list = []
    model_results_list = []
    model_analysis_list = []
    model_analysis_dict = {}
    dgr_list = []
    fit_params_list = []
    cloud_name_list = []
    core_names_list = []
    core_list = []
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
        dgr_list.append(results_dict['mc_analysis']['dgr'])
        #fit_params_list.append(results_dict['params_summary'])
        model_results_list.append(results_dict['mc_results']['ss_model_results'])
        model_analysis_list.append(results_dict['mc_analysis'])
        model_analysis_dict[cloud_name] = results_dict['mc_analysis']

        # Modeling and cores
        global_args = results_dict['global_args']
        cores = global_args['ss_model_kwargs']['cores']
        cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']
        rh2_core_list = []
        hsd_core_list = []
        hisd_core_list = []
        core_names = []
        core_list.append(cores)
        for core in cores_to_plot:
            core_indices = cores[core]['indices_orig']
            core_names.append(core)
            hisd_core_list.append(hisd_list[i][core_indices])
            hsd_core_list.append(hsd_list[i][core_indices])
            rh2_core_list.append(rh2_list[i][core_indices])

            if 0:
                rh2_copy = rh2_list[i].copy()
                rh2_copy[core_indices] = 1000
                plt.imshow(rh2_copy, origin='lower')
                plt.savefig('/d/bip3/ezbc/scratch/core_' + core + '.png')

        core_names_list.append(core_names)
        hisd_cores_list.append(hisd_core_list)
        hsd_cores_list.append(hsd_core_list)

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

'''
Data Prep functions
'''

def collect_hi_transition_results(model_analysis_list, cloud_list,
        filename=None):

    import pandas as pd

    hi_trans_dict = {}

    for i, cloud in enumerate(cloud_list):
        for core_name in model_analysis_list[i]['cores']:
            core = model_analysis_list[i]['cores'][core_name]
            hi_trans_k09 = \
                core['krumholz_results']['hi_transition']
            #print hi_trans_k09
            hi_trans_k09_error = \
                core['krumholz_results']['hi_transition_error']
            hi_trans_s14 = \
                core['sternberg_results']['hi_transition']
            hi_trans_s14_error = \
                core['sternberg_results']['hi_transition_error']

            hi_trans_dict[core_name] = \
                {'cloud': cloud,
                 'k09_transition': hi_trans_k09,
                 'k09_transition_error': hi_trans_k09_error,
                 's14_transition': hi_trans_s14,
                 's14_transition_error': hi_trans_s14_error,
                 }

    return hi_trans_dict

def print_dict_keys(d):

    for key in d:
        print(key)
        if type(d[key]) is dict:
            print('--')
            print_dict_keys(d[key])

def write_fit_summary_dict(mc_analysis_dict, core_list, cloud_name_list,
        filename):

    import pandas as pd
    import pickle

    d = {}

    # Collect parameter names for each model for each core
    for i, cloud in enumerate(cloud_name_list):
        mc_analysis = mc_analysis_dict[cloud]
        core_names = np.sort(mc_analysis['cores'].keys())

        cores = core_list[i]

        # each core correspond to a row
        for cloud_row, core_name in enumerate(core_names):
            core = mc_analysis['cores'][core_name]

            core_props = cores[core_name]

            d[core_name] = {}
            core_new = d[core_name]
            core_new['cloud'] = cloud
            ra_deg = core_props['center_wcs'][0]
            dec_deg = core_props['center_wcs'][1]
            #ra_deg = 15*(ra[0] + ra[1] / 60. + ra[2] / 3600.)
            #dec_deg = dec[0] + dec[1] / 60. + dec[2] / 3600.
            core_new['ra'] = ra_deg
            core_new['dec'] = dec_deg
            try:
                core_new['temp'] = core_props['temp']
                core_new['temp_error'] = core_props['temp_error']
            except KeyError:
                core_new['temp'], core_new['temp_error'] = 17, 1
            core_new['region_vertices'] = core_props['poly_verts']['wcs']

            # append model params and errors to row
            for model in ('krumholz', 'sternberg'):
                if model == 'krumholz':
                    params_to_write = ['phi_cnm', 'Z', 'phi_mol',
                                       'hi_transition']
                else:
                    params_to_write = ['alphaG', 'Z', 'phi_g',
                                       'hi_transition']
                core_new[model] = {}
                for i, param_name in enumerate(params_to_write):
                    param = \
                        core[model + '_results'][param_name]
                    param_error = \
                        core[model + '_results'][param_name + '_error']

                    core_new[model][param_name] = param
                    #core_new[model][param_name + '_error_low'] = \
                    #        param_error[0]
                    #core_new[model][param_name + '_error_high'] = \
                    #        param_error[1]
                    core_new[model][param_name + '_error'] = param_error

    with open(filename, 'wb') as f:
        pickle.dump(d, f)

def write_param_csv(mc_analysis_dict, core_list, cloud_name_list, filename):

    import pandas as pd

    d = {}

    #print_dict_keys(mc_analysis_dict)
    params_to_write = ['phi_cnm', 'alphaG', 'hi_transition']

    d['cloud'] = []
    d['core'] = []
    d['ra'] = []
    d['dec'] = []
    d['region_vertices'] = []
    for param in params_to_write:
        d[param] = []
        d[param + '_error_low'] = []
        d[param + '_error_high'] = []

    # Collect parameter names for each model for each core
    for i, cloud in enumerate(cloud_name_list):
        mc_analysis = mc_analysis_dict[cloud]
        core_names = np.sort(mc_analysis['cores'].keys())

        cores = core_list[i]

        # each core correspond to a row
        for cloud_row, core_name in enumerate(core_names):
            core = mc_analysis['cores'][core_name]

            core_props = cores[core_name]

            d['cloud'].append(cloud)
            d['core'].append(core_name)
            ra_deg = core_props['center_wcs'][0]
            dec_deg = core_props['center_wcs'][1]
            #ra_deg = 15*(ra[0] + ra[1] / 60. + ra[2] / 3600.)
            #dec_deg = dec[0] + dec[1] / 60. + dec[2] / 3600.
            d['ra'].append(ra_deg)
            d['dec'].append(dec_deg)
            d['region_vertices'].append(core_props['poly_verts']['wcs'])

            # append model params and errors to row
            for model in ('krumholz', 'sternberg'):
                if model == 'krumholz':
                    params_to_write = ['phi_cnm', 'hi_transition']
                else:
                    params_to_write = ['alphaG',]
                for i, param_name in enumerate(params_to_write):
                    param = \
                        core[model + '_results'][param_name]
                    param_error = \
                        core[model + '_results'][param_name + '_error']

                    d[param_name].append(param)
                    d[param_name + '_error_low'].append(param_error[0])
                    d[param_name + '_error_high'].append(param_error[1])

            #print d[param_name]

    # Create dataframe and write it!
    df = pd.DataFrame(data=d,)
    df.to_csv(filename,
              sep=',',
              columns=('cloud',
                       'core',
                       'ra',
                       'dec',
                       'phi_cnm',
                       'phi_cnm_error_low',
                       'phi_cnm_error_high',
                       'alphaG',
                       'alphaG_error_low',
                       'alphaG_error_high',
                       'hi_trans',
                       'hi_trans_error_low',
                       'hi_trans_error_high',
                       ),
              index=False,
              )

    df.save(filename.replace('csv', 'pickle'))

def write_hi_vel_range_table(names_list, hi_range_kwargs_list, filename):

    # Open file to be appended
    f = open(filename, 'wb')

    for i in xrange(0, len(names_list)):
        cloud_name = names_list[i]
        hi_range_kwargs = hi_range_kwargs_list[i]
        vel_range = hi_range_kwargs['vel_range']
        vel_range_error = hi_range_kwargs['hi_range_error']

        row_text = cloud_name.capitalize()

        row_text = add_row_element(row_text,
                        vel_range,
                        text_format='[{0:.0f}, {1:.0f}]')
        row_text = add_row_element(row_text,
                        vel_range_error,
                        text_format='{0:.0f}')
        row_text += ' \\\\[0.1cm] \n'

        f.write(row_text)

    f.close()



def write_model_params_table(mc_analysis_dict, filename, models=('krumholz',)):

    # Open file to be appended
    f = open(filename, 'wb')

    text_param_format ='{0:.1f}$^{{+{1:.1f}}}_{{-{2:.1f}}}$'

    #print_dict_keys(mc_analysis_dict)
    params_to_write = ['phi_cnm', 'alphaG']

    # Collect parameter names for each model for each core
    for cloud in ('california', 'perseus', 'taurus'):
        mc_analysis = mc_analysis_dict[cloud]
        core_names = np.sort(mc_analysis['cores'].keys())
        # each core correspond to a row
        for cloud_row, core_name in enumerate(core_names):
            core = mc_analysis['cores'][core_name]
            if cloud_row == 0:
                row_text = cloud.capitalize()
            else:
                row_text = ''
            row_text = add_row_element(row_text,
                                       core_name)

            # append model params and errors to row
            for model in ('sternberg', 'krumholz',):
                if model == 'krumholz':
                    params_to_write = ['phi_cnm', 'hi_transition']
                else:
                    params_to_write = ['alphaG',]
                for i, param_name in enumerate(params_to_write):
                    param = \
                        core[model + '_results'][param_name]
                    param_error = \
                        core[model + '_results'][param_name + '_error']

                    param_info = (param, param_error[1], param_error[0])

                    #if param_name == 'alphaG':
                        #print core_name, param_info

                    row_text = \
                        add_row_element(row_text,
                                        param_info,
                                        text_format=text_param_format)


            row_text += ' \\\\[0.1cm] \n'
            if cloud_row == len(mc_analysis['cores']) - 1 \
                and cloud != 'taurus':
                row_text += '\hline  \\\\[-0.2cm] \n'
            elif cloud_row == len(mc_analysis['cores']) - 1 and \
                    cloud == 'taurus':
                row_text.replace(r'\\[0.1cm] \n', '')

            f.write(row_text)

    f.close()

def add_row_element(row_text, element, text_format='{0:s}'):

    if type(element) is list or type(element) is tuple:
        return row_text + ' & ' + text_format.format(*element)
    else:
        return row_text + ' & ' + text_format.format(element)


'''
Multiprocessing functions
'''


class QueueGet(Queue):
    """Queue which will retry if interrupted with EINTR."""
    def get( block=True, timeout=None):
        return retry_on_eintr(Queue.get,  block, timeout)

def retry_on_eintr(function, *args, **kw):
    from multiprocessing.queues import Queue
    import errno

    while True:
        try:
            return function(*args, **kw)
        except IOError, e:
            if e.errno == errno.EINTR:
                continue
            else:
                raise

def _my_queue_get(queue, block=True, timeout=None):
    import errno
    while True:
        try:
            return queue.get(block, timeout)
        except IOError, e:
            if e.errno != errno.EINTR:
                raise

class KeyboardInterruptError(Exception): pass

'''
Av core functions
'''
def read_ds9_region(filename):

    ''' Converts DS9 region file into format for plotting region.

    Need the following format:
        angle : degrees
        xy : pixels
        width : pixels
        height : pixels

    Region file provides following format:
        # Region file format: DS9 version 4.1
        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        fk5
        box(4:17:04.740,+29:20:31.32,5854.33",11972.7",130) # text={test}

    pyregion module reads DS9 regions:
    http://leejjoon.github.io/pyregion/users/overview.html


    '''

    # Import external modules
    import pyregion as pyr

    # Read region file
    try:
        region = pyr.open(filename)
    except IOError:
        return None

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    return region

def load_ds9_core_region(cores, filename='',
        header=None):

    from myimage_analysis import get_pix_coords

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    regions = read_ds9_region(filename)

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        core = tag[tag.find('{')+1:tag.find('}')]

        if core in cores:
            # Format vertices to be 2 x N array
            poly_verts = []
            for i in xrange(0, len(region.coord_list)/2):
                poly_verts.append((region.coord_list[2*i],
                                   region.coord_list[2*i+1]))

            poly_verts_pix = []
            for i in xrange(0, len(poly_verts)):
                poly_verts_pix.append(get_pix_coords(ra=poly_verts[i][0],
                                          dec=poly_verts[i][1],
                                          header=header)[:-1][::-1].tolist())

            cores[core]['poly_verts'] = {}
            cores[core]['poly_verts']['wcs'] = poly_verts
            cores[core]['poly_verts']['pixel'] = poly_verts_pix

    return cores

def load_wedge_core_region(cores, filename='', header=None):

    from myimage_analysis import get_pix_coords
    import pickle
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from astropy.io import fits
    from astropy.wcs import WCS

    # Create WCS object
    wcs_header = WCS(header)

    with open(filename, 'rb') as f:
        region_dict = pickle.load(f)

    for core_name in region_dict:
        if core_name in cores:
            core = region_dict[core_name]
            # Format vertices to be 2 x N array
            poly_verts = np.array((core['ra'], core['dec']))

            # Make a galactic coords object and convert to Ra/dec
            coords_fk5 = SkyCoord(core['ra'] * u.deg,
                                  core['dec'] * u.deg,
                                  frame='fk5',
                                  )
            # convert to pixel
            coords_pixel = np.array(coords_fk5.to_pixel(wcs_header))

            #print coords_pixel
            #print coords_pixel.shape

            # write data to dataframe
            poly_verts_pix = np.array((coords_pixel[1], coords_pixel[0])).T

            #print poly_verts_pix.shape

            #poly_verts_pix = np.array((core['ypix'], core['xpix'])).T

            cores[core_name]['poly_verts'] = {}
            cores[core_name]['poly_verts']['wcs'] = poly_verts
            cores[core_name]['poly_verts']['pixel'] = poly_verts_pix

    return cores

def get_cores_to_plot():

    '''

    '''

    # Which cores to include in analysis?
    cores_to_keep = [ # N(H2) cores
                     'G168.54-6.22',
                     'G168.12-6.42',
                     'G166.91-7.76',
                     'G165.36-7.51',
                     'G164.70-7.63',
                     'G164.65-8.12',
                     'G165.71-9.15',
                     'G164.99-8.60',
                     'G164.26-8.39',
                     'G164.18-8.84',
                     'G174.40-13.45',
                     'G174.70-15.47',
                     'G174.05-15.82',
                     'G171.49-14.91',
                     'G172.93-16.73',
                     'G171.00-15.80',
                     'G172.12-16.94',
                     'G171.14-17.57',
                     'G169.32-16.17',
                     'G168.10-16.38',
                     'G160.49-16.81',
                     'G160.46-17.99',
                     'G159.80-18.49',
                     'G160.14-19.08',
                     'G160.53-19.73',
                     'G159.19-20.11',
                     'G159.17-21.09',
                     'G158.39-20.72',
                     'G158.89-21.60',
                     'G158.26-21.81',
                     ]

    if 0: # Random cores
        cores_to_keep = [
                 'G166.83-8.68',
                 'G168.82-6.37',
                 'G168.05-7.01',
                 'G164.16-8.46',
                 'G165.23-8.78',
                 'G167.06-7.77',
                 'G168.12-6.42',
                 'G167.58-6.64',
                 'G164.70-7.63',
                 'G166.35-8.77',
                 'G166.73-15.06',
                 'G173.08-16.50',
                 'G172.74-14.53',
                 'G169.44-16.18',
                 'G173.86-17.65',
                 'G173.71-13.91',
                 'G171.75-14.18',
                 'G173.70-15.21',
                 'G170.28-19.48',
                 'G171.00-15.80',
                 'G158.23-20.15',
                 'G159.01-22.19',
                 'G159.19-20.11',
                 'G157.12-23.49',
                 'G160.10-19.90',
                 'G160.34-18.42',
                 'G158.40-21.86',
                 'G159.79-21.32',
                 'G158.89-21.60',
                 'G159.51-18.41',
                ]

    if 0:
        cores_to_keep = [# taur
                         'L1495',
                         'L1495A',
                         'B213',
                         'L1498',
                         'B215',
                         'B18',
                         'B217',
                         'B220-1',
                         'B220-2',
                         'L1521',
                         'L1524',
                         'L1527-1',
                         'L1527-2',
                         # Calif
                         'L1536',
                         'L1483-1',
                         'L1483-2',
                         'L1482-1',
                         'L1482-2',
                         'L1478-1',
                         'L1478-2',
                         'L1456',
                         'NGC1579',
                         #'L1545',
                         #'L1517',
                         #'L1512',
                         #'L1523',
                         #'L1512',
                         # Pers
                         'B5',
                         'IC348',
                         'B1E',
                         'B1',
                         'NGC1333',
                         'B4',
                         'B3',
                         'L1455',
                         'L1448',
                         ]

    return cores_to_keep

def get_core_properties(data_dict, cloud_name):

    from myimage_analysis import load_ds9_region, get_pix_coords
    import json

    box_method = 'ds9'
    core_dir = '/d/bip3/ezbc/multicloud/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'

    header = data_dict['av_header']
    if 0:
        # define core properties
        with open(core_dir + cloud_name + '_core_properties.txt', 'r') as f:
            cores = json.load(f)

        #cores = convert_core_coordinates(cores, header)

        # convert center WCS coordinate to pixel
        for core in cores:
            cores[core].update({'box_pixel': 0})
            cores[core].update({'center_pixel': 0})

            center_wcs = cores[core]['center_wcs']

            # convert centers to pixel coords
            center_pixel = get_pix_coords(ra=center_wcs[0],
                                          dec=center_wcs[1],
                                          header=header)[:2]
            cores[core]['center_pixel'] = center_pixel

    else:
        filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
            'planck11_coldclumps.pickle'

        import pickle
        with open(filename, 'rb') as f:
            core_sample = pickle.load(f)
        #core_sample = pd.load(filename)
        cores = {}
        df_cores = core_sample[cloud_name]
        for name in df_cores['Name'].values:
            cores[name] = {}
            ra = df_cores[df_cores['Name'] == name]['ra'].values[0]
            dec = df_cores[df_cores['Name'] == name]['dec'].values[0]
            cores[name]['center_wcs'] = (ra, dec)
            cores[name]['temp'] = \
                df_cores[df_cores['Name'] == name]['temp'].values[0]
            cores[name]['temp_error'] = \
                df_cores[df_cores['Name'] == name]['temp_error'].values[0]

            #print cores[name]['temp_error']

            # convert centers to pixel coords
            center_pixel = get_pix_coords(ra=ra,
                                          dec=dec,
                                          header=header)[:2]
            cores[name]['center_pixel'] = center_pixel[::-1]

    # load the bounding regions
    if 0:
        region_filename = region_dir + 'multicloud_coldclump_divisions.reg'
        cores = load_ds9_core_region(cores,
                                filename=region_filename,
                                header=header)
    else:
        filename = region_dir + 'multicloud_divisions_coldcore_wedges.pickle'
        cores = load_wedge_core_region(cores,
                                       filename=filename,
                                       header=header)

    # add indices of core in data
    add_core_mask(cores, data_dict['av_data'])

    return cores

def trim_cores_to_plot(cores, cores_to_plot):

    # Trim down cores to keep list to include only cores for cloud
    cores_to_keep_old = list(cores_to_plot)
    for core in cores_to_keep_old:
        if core not in cores or 'poly_verts' not in cores[core]:
            cores_to_plot.remove(core)

    return cores_to_plot

def add_core_mask(cores, data):

    for core in cores:
        try:
            vertices = cores[core]['poly_verts']['pixel']

            mask = np.logical_not(myg.get_polygon_mask(data,
                                                       vertices))

            cores[core]['indices_orig'] = np.where(mask == 0)

            cores[core]['mask'] = mask

            #print 'npix in core ' + core + ':'
            #print np.where(mask == 0)[0].size

        except KeyError:
            cores[core]['mask'] = None
            cores[core]['indices_orig'] = None

'''
Modeling Functions
'''

def fit_av_model(av, nhi, av_error=None, algebraic=False, nhi_background=None,
        plot_kwargs=None, init_guesses=[0.05, 0.05, 0], use_intercept=True,
        return_fit=False, fit_method='lbfgsb'):

    from lmfit import minimize, Parameters
    import lmfit

    if nhi_background is None:
        use_background = False
        init_guesses[1] = 0.0
        nhi_background = 0.0
    else:
        use_background = True

    # Use linear algebra?
    if algebraic:
        b = av
        A = np.array([nhi, np.ones(nhi.shape)]).T

        # weights
        if av_error is not None:
            W = 1.0 / av_error**2
        else:
            W = np.ones(av.shape)

        A = np.array([nhi, np.ones(nhi.shape)]).T

        params = np.dot(np.linalg.pinv(A), b)

    else:
        # Set parameter limits and initial guesses
        params = Parameters()
        params.add('dgr_cloud',
                   value=init_guesses[0],
                   min=-0.5,
                   max=1,
                   )
        params.add('dgr_background',
                   value=init_guesses[1],
                   min=0.0,
                   max=1,
                   vary=use_background,
                   )
        params.add('intercept',
                   value=init_guesses[2],
                   min=-5,
                   max=5,
                   vary=use_intercept,
                   )

        #bin_edges = residuals_crop
        #counts = np.ones(residuals_crop.size - 1)

        def norm(params, av, nhi, av_error=None, nhi_background=None,):
            if nhi_background is None:
                nhi_background = 0.0
            if av_error is None:
                av_error = np.ones(av.shape)


            model = params['dgr_cloud'] * nhi + \
                    params['dgr_background'] * nhi_background + \
                    params['intercept']

            #if fit_method == 'leastsq':
                #norm = np.sum((av - model)**2 * (1.0/av_error**2)) / \
                #       np.sum(1.0/av_error**2)
            norm = np.sum((av - model)**2 * (1.0/av_error**2)) / \
                    np.sum(1.0/av_error**2)
            #norm = np.sum((av - model)**2)
            #else:
            #    norm = (av - model)

            return norm

        #print('fitting')
        # Perform the fit!
        result = minimize(norm,
                          params,
                          args=(av, nhi, av_error, nhi_background),
                          #method='leastsq',
                          #method=fit_method,
                          method='nelder',
                          )

        #print lmfit.report_fit(params)
        #print lmfit.printfuncs.report_ci(lmfit.conf_interval(result))
        #print lmfit.conf_interval(result)

        dgr_cloud = params['dgr_cloud'].value
        dgr_background = params['dgr_background'].value
        intercept = params['intercept'].value

        #if debugging:
        if 0:
            print('dgr = ', dgr_cloud)
            print('dgr background = ', dgr_background)
            print('intercept = ', intercept)
        if 0:
            plt.close(); plt.clf()
            background = dgr_background * nhi_background
            plt.errorbar(nhi, av - background,
                     yerr=(av_error),
                     linestyle='',
                     marker='o',
                     alpha=0.1,
                     markersize=2)
            xfit = np.linspace(0,50)
            plt.plot(xfit, dgr_cloud * xfit + intercept)
            plt.xlim(0, 22)
            plt.ylim(0, 15)
            plt.xlabel(r'N(HI)')
            plt.ylabel(r'$A_V$')
            plt.savefig(plot_kwargs['figure_dir'] + \
                        'diagnostics/av_nhi/' + plot_kwargs['filename_base']+ \
                        '_avnhi_bootsrap' + \
                        '{0:03d}.png'.format(plot_kwargs['bootstrap_num']))

        results = {'dgr_cloud': dgr_cloud,
                  'dgr_background': dgr_background,
                  'intercept': intercept}

        if return_fit:
            return (result, params), (dgr_cloud, dgr_background, intercept)
        else:
            return results

def fit_steady_state_models(h_sd, rh2, model_kwargs, rh2_error=None,
        h_sd_error=None, bootstrap_residuals=False, nboot=100, G0=1.0):

    # Fit R_H2
    #---------
    sternberg_params = model_kwargs['sternberg_params']
    sternberg_results = {}
    krumholz_params = model_kwargs['krumholz_params']
    krumholz_results = {}

    # Fit to sternberg model
    if rh2.size > 3:
        result = \
            fit_sternberg(h_sd,
                          rh2,
                          guesses=sternberg_params['guesses'],
                          vary=sternberg_params['param_vary'],
                          radiation_type=sternberg_params['radiation_type'],
                          bootstrap_residuals=bootstrap_residuals,
                          nboot=nboot,
                          rh2_error=rh2_error,
                          h_sd_error=h_sd_error,
                          )
        if bootstrap_residuals:
            alphaG, alphaG_error, Z_s14, Z_s14_error, phi_g, phi_g_error = \
                    result
            sternberg_results['alphaG_error'] = alphaG_error
            sternberg_results['Z_error'] = Z_s14_error
            sternberg_results['phi_g_error'] = phi_g_error
        else:
            alphaG, Z_s14, phi_g = result
            alphaG_error, Z_s14_error, phi_g_error = 3*[np.nan]

        # Fit to krumholz model
        result = \
            fit_krumholz(h_sd,
                         rh2,
                         guesses=krumholz_params['guesses'],
                         vary=krumholz_params['param_vary'],
                         bootstrap_residuals=bootstrap_residuals,
                         nboot=nboot,
                         G0=G0,
                         rh2_error=rh2_error,
                         h_sd_error=h_sd_error,
                         )
        if bootstrap_residuals:
            phi_cnm, phi_cnm_error,Z_k09, Z_k09_error,phi_mol, phi_mol_error= \
                    result
            krumholz_results['phi_cnm_error'] = phi_cnm_error
            krumholz_results['Z_error'] = Z_k09_error
            krumholz_results['phi_mol_error'] = phi_mol_error
        else:
            phi_cnm, Z_k09, phi_mol = result
            phi_cnm_error, Z_k09_error, phi_mol_error = 3*[np.nan]
    else:
        alphaG, Z_s14, phi_g, phi_cnm, Z_k09, phi_mol = 6 * [np.nan]
        alphaG_error, Z_s14_error, phi_g_error, \
        phi_cnm_error, Z_k09_error, phi_mol_error = 6*[np.nan]

    # keep results
    sternberg_results['alphaG'] = alphaG
    sternberg_results['Z'] = Z_s14
    sternberg_results['phi_g'] = phi_g

    # keep results
    krumholz_results['phi_cnm'] = phi_cnm
    krumholz_results['Z'] = Z_k09
    krumholz_results['phi_mol'] = phi_mol

    # see eq 6 of sternberg+09
    # alphaG is the number density of the CNM over the minimum number
    # density required for pressure balance
    # the lower alphaG values than for taurus mean that taurus
    # has a more diffuse CNM

    # By fitting the model to the observed R_H2 vs total H, you
    # basically constrained psi in Equation (35) of sternberg+09.  This
    # means that you can calculate f_H2 for a given total hydrogen
    # surface density.  In this case, H2 surface density = f_H2 *
    # total hydrogen surface density HI surface density = (1 - f_HI) *
    # total hydrogen surface density

    results = {}
    results['sternberg_results'] = sternberg_results
    results['krumholz_results'] = krumholz_results

    return results

def add_hi_transition_calc(ss_model_result):

    h_sd_fit = np.linspace(0, 100, 1000)

    # To get HI transition, calculate model fits, then find where RH2 = 1
    for model_name in ss_model_result:
        model = ss_model_result[model_name]

        params = {}

        if 'sternberg' in model_name:
            params['phi_g'] = model['phi_g']
            params['Z'] = model['Z']
            params['alphaG'] = model['alphaG']
            model_fits = calc_sternberg(params,
                                      h_sd=h_sd_fit,
                                      return_fractions=False,
                                      return_hisd=False,
                                      )
        elif 'krumholz' in model_name:
            params['phi_cnm'] = model['phi_cnm']
            params['Z'] = model['Z']
            params['phi_mol'] = model['phi_mol']
            model_fits = calc_krumholz(params,
                                      h_sd=h_sd_fit,
                                      return_fractions=False,
                                      return_hisd=False,
                                      )

        rh2_fit = model_fits[0]

        try:
            if np.isnan(np.sum(rh2_fit)):
                hi_transition = np.nan
            else:
                # when R_H2 = 1, HI-to-H2 transition
                hi_transition = np.interp(1, rh2_fit, h_sd_fit) / 2.0

        except ValueError:
            hi_transition = np.nan

        model['hi_transition'] = hi_transition

def calc_krumholz(params, h_sd_extent=(0.001, 500), return_fractions=True,
        return_hisd=False, h_sd=None):

    '''
    Parameters
    ----------
    phi_cnm, Z : float
        Phi_cnm and Z parameters for Krumholz model.
    h_sd_extent : tuple
        Lower and upper bound of hydrogen surface densities with which to
        build the output model array.
    return_fractions : bool
        Return f_H2 and f_HI?

    Returns
    -------
    rh2_fit : array-like
        Model ratio between molecular and atomic hydrogen masses.
    h_sd_extended : list
        Model hydrogen surface density in units of solar mass per parsec**2.
    f_H2, f_HI : array-like, optional
        f_H2 = mass fraction of molecular hydrogen
        f_HI = mass fraction of atomic hydrogen

    '''

    from scipy import stats
    from myscience import krumholz09 as k09

    # Get large array of h_sd
    if h_sd is None:
        h_sd = np.linspace(h_sd_extent[0], h_sd_extent[1], 1e2)

    params = [params['phi_cnm'], params['Z'], params['phi_mol']]
    if params[0] <= 0 or np.isnan(params[0]):
        rh2_fit, f_H2, f_HI = np.empty(1), np.empty(1), np.empty(1)
    else:
        rh2_fit, f_H2, f_HI = k09.calc_rh2(h_sd,
                                           phi_cnm=params[0],
                                           Z=params[1],
                                           phi_mol=params[2],
                                           return_fractions=True)

    output = [rh2_fit, h_sd]

    if return_fractions:
        output.append(f_H2)
        output.append(f_HI)
    if return_hisd:
        #hi_sd = f_HI * h_sd
        hi_sd = (1 - f_H2) * h_sd
        output.append(hi_sd)

    return output

def fit_krumholz(h_sd, rh2, guesses=[10.0, 1.0, 10.0], h_sd_error=None,
        rh2_error=None, verbose=False, vary=[True, True, True],
        bootstrap_residuals=False, nboot=100, G0=1.0):

    '''
    Parameters
    ----------
    h_sd : array-like
        Hydrogen surface density in units of solar mass per parsec**2
    rh2 : array-like
        Ratio between molecular and atomic hydrogen masses.
    guesses : None, scalar, or M-length sequence.
        Initial guess for the parameters. See scipy.optimize.curve_fit.
    rh2_error : bool
        Error in rh2 parameter. Calculates a more accurate chi^2 statistic
    G0 : float
        radiation field

    Returns
    -------
    rh2_fit_params : array-like, optional
        Model parameter fits.

    '''

    from scipy.optimize import curve_fit
    from scipy import stats
    from lmfit import minimize, Parameters, report_fit
    from myscience import krumholz09 as k09
    import mystats

    def chisq(params, h_sd, rh2):
        phi_cnm = params['phi_cnm'].value
        phi_mol = params['phi_mol'].value
        Z = params['Z'].value

        rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol)

        chisq = np.sum(np.abs(rh2 - rh2_model))

        return chisq

    def calc_residual(params, h_sd, rh2, G0):
        phi_cnm = params['phi_cnm'].value
        phi_mol = params['phi_mol'].value
        Z = params['Z'].value

        rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol, G_0=G0)

        residual = rh2 - rh2_model

        return residual

    # Set parameter limits and initial guesses
    params = Parameters()
    params.add('phi_cnm',
               value=guesses[0],
               min=0.001,
               max=1000,
               vary=vary[0])
    params.add('phi_mol',
               value=guesses[2],
               min=1,
               max=20,
               vary=vary[2])
    params.add('Z',
               value=guesses[1],
               min=0.1,
               max=4,
               vary=vary[1])

    # Perform the fit!
    result = minimize(calc_residual,
                      params,
                      args=(h_sd, rh2, G0),
                      method='leastsq')

    if bootstrap_residuals:
        def resample_residuals(residuals):
            return np.random.choice(residuals,
                                    size=residuals.size,
                                    replace=True)
        phi_cnm = params['phi_cnm'].value
        phi_mol = params['phi_mol'].value
        Z = params['Z'].value
        rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol)
        residual = rh2 - rh2_model

        empty = np.empty(nboot)
        param_dict = {}
        param_names = ('phi_cnm', 'Z', 'phi_mol')
        for param_name in param_names:
            param_dict[param_name] = empty.copy()

        for i in xrange(nboot):
            rh2_resampled = rh2_model + resample_residuals(residual)
            result = minimize(calc_residual,
                              params,
                              args=(h_sd, rh2_resampled),
                              #method='anneal',
                              method='leastsq',
                              )

            for param_name in param_names:
                param_dict[param_name][i] = params[param_name].value

        rh2_fit_params = []
        for param_name in param_names:
            conf = mystats.calc_cdf_error(param_dict[param_name])
            rh2_fit_params.append(conf[0])
            rh2_fit_params.append(conf[1])
            #print param_name, conf
    else:
        rh2_fit_params = (params['phi_cnm'].value, params['Z'].value,
                params['phi_mol'].value)

    return rh2_fit_params

def analyze_krumholz_model(krumholz_results):

    ''' Calculates various properties of Krumholz model from fitted phi_cnm
    values.

    '''

    from myscience.krumholz09 import calc_T_cnm

    # Calculate T_cnm from Krumholz et al. (2009) Eq 19
    phi_cnm, phi_cnm_error, Z = \
        [krumholz_results[key] for key in ('phi_cnm', 'phi_cnm_error', 'Z')]
    T_cnm = calc_T_cnm(phi_cnm, Z=Z)
    T_cnm_error = []
    T_cnm_error.append(\
            T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[0], Z=Z))
    T_cnm_error.append(\
            T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[1], Z=Z))

    krumholz_results.update({'T_cnm': T_cnm,
                             'T_cnm_error': T_cnm_error})

    # Get fitted surface density ratios...
    params = [krumholz_results[param] for param in \
              krumholz_results['parameters']]

    rh2_fit, h_sd_fit, f_H2, f_HI = \
           calc_krumholz(params=params,
                          h_sd_extent=krumholz_results['h_sd_fit_range'],
                          return_fractions=True)

    krumholz_results.update({'rh2_fit': rh2_fit,
                             'h_sd_fit' : h_sd_fit,
                             'hi_sd_fit' : f_HI * h_sd_fit,
                             'f_H2' : f_H2,
                             'f_HI' : f_HI,})

    return krumholz_results

def calc_sternberg(params, h_sd_extent=(0.001, 500), return_fractions=True,
        return_hisd=False, h_sd=None):

    '''
    Parameters
    ----------
    alphaG, Z : float
        alphaG and Z parameters for sternberg model.
    h_sd_extent : tuple
        Lower and upper bound of hydrogen surface densities with which to
        build the output model array.
    return_fractions : bool
        Return f_H2 and f_HI?

    Returns
    -------
    rh2_fit : array-like
        Model ratio between molecular and atomic hydrogen masses.
    h_sd_extended : list
        Model hydrogen surface density in units of solar mass per parsec**2.
    f_H2, f_HI : array-like, optional
        f_H2 = mass fraction of molecular hydrogen
        f_HI = mass fraction of atomic hydrogen

    '''

    from scipy import stats
    from myscience import sternberg14 as s14

    # Create large array of h_sd
    if h_sd is None:
        h_sd = np.linspace(h_sd_extent[0], h_sd_extent[1], 1e2)

    params = [params['alphaG'], params['Z'], params['phi_g']]
    if params[0] <= 0 or np.isnan(params[0]):
        rh2_fit, f_H2, f_HI = np.empty(1), np.empty(1), np.empty(1)
    else:
        rh2_fit, f_H2, f_HI = s14.calc_rh2(h_sd,
                                           alphaG=params[0],
                                           Z=params[1],
                                           phi_g=params[2],
                                           return_fractions=True)

    output = [rh2_fit, h_sd]

    if return_fractions:
        output.append(f_H2)
        output.append(f_HI)
    if return_hisd:
        hi_sd = f_HI * h_sd
        output.append(hi_sd)

    return output

def fit_sternberg(h_sd, rh2, guesses=[10.0, 1.0, 10.0], rh2_error=None,
        h_sd_error=None, verbose=False, vary=[True, True, True],
        radiation_type='isotropic', bootstrap_residuals=False, nboot=100,
        odr_fit=True):

    '''
    Parameters
    ----------
    h_sd : array-like
        Hydrogen surface density in units of solar mass per parsec**2
    rh2 : array-like
        Ratio between molecular and atomic hydrogen masses.
    guesses : None, scalar, or M-length sequence.
        Initial guess for the parameters. See scipy.optimize.curve_fit.
    rh2_error : bool
        Error in rh2 parameter. Calculates a more accurate chi^2 statistic

    Returns
    -------
    rh2_fit_params : array-like, optional
        Model parameter fits.


    Residual bootstrapping:
    http://stats.stackexchange.com/questions/67519/bootstrapping-residuals-am-i-doing-it-right


    '''

    from scipy.optimize import curve_fit
    from scipy import stats
    from lmfit import minimize, Parameters, report_fit, Minimizer
    import lmfit
    from myscience import sternberg14 as s14
    import mystats

    def chisq(params, h_sd, rh2):
        alphaG = params['alphaG'].value
        phi_g = params['phi_g'].value
        Z = params['Z'].value

        rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                 return_fractions=False,
                                 radiation_type=radiation_type)

        chisq = np.sum(np.abs(rh2 - rh2_model))

        return chisq

    def calc_residual(params, h_sd, rh2):
        alphaG = params['alphaG'].value
        phi_g = params['phi_g'].value
        Z = params['Z'].value

        rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                 return_fractions=False)

        residual = rh2 - rh2_model

        return residual

    # Set parameter limits and initial guesses
    params = Parameters()
    params.add('alphaG',
               value=guesses[0],
               min=0.1,
               max=500,
               vary=vary[0])
    params.add('phi_g',
               value=guesses[2],
               min=0.01,
               max=10,
               vary=vary[2])
    params.add('Z',
               value=guesses[1],
               min=0.1,
               max=4,
               vary=vary[1])

    # Perform the fit!
    result = minimize(calc_residual,
                      params,
                      args=(h_sd, rh2),
                      #method='leastsq',
                      method='anneal',
                      )

    if bootstrap_residuals:
        def resample_residuals(residuals):
            return np.random.choice(residuals,
                                    size=residuals.size,
                                    replace=True)
        alphaG = params['alphaG'].value
        phi_g = params['phi_g'].value
        Z = params['Z'].value
        rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                 return_fractions=False)
        residual = rh2 - rh2_model

        empty = np.empty(nboot)
        param_dict = {}
        param_names = ('alphaG', 'Z', 'phi_g',)
        for param_name in param_names:
            param_dict[param_name] = empty.copy()

        for i in xrange(nboot):
            rh2_resampled = rh2_model + resample_residuals(residual)
            result = minimize(calc_residual,
                              params,
                              args=(h_sd, rh2_resampled),
                              #method='leastsq',
                              method='anneal',
                              )

            for param_name in param_names:
                param_dict[param_name][i] = params[param_name].value

        rh2_fit_params = []
        for param_name in param_names:
            conf = mystats.calc_cdf_error(param_dict[param_name])
            rh2_fit_params.append(conf[0])
            rh2_fit_params.append(conf[1])
            #print param_name, conf
    else:
        rh2_fit_params = (params['alphaG'].value, params['Z'].value,
                params['phi_g'].value)

    if odr_fit:
        import scipy.odr as odr

        alphaG, Z, phi_g = guesses
        def odr_func(alphaG, h_sd):

            return s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                     return_fractions=False)

        model = odr.Model(odr_func)
        data = odr.RealData(h_sd, rh2, sx=0.1 * h_sd, sy=rh2_error)
        odr_instance = odr.ODR(data, model, beta0=[10.0,])
        output = odr_instance.run()
        alphaG = output.beta[0]
        print 'alphaG =', alphaG

        rh2_fit_params = (alphaG, Z, phi_g)

    return rh2_fit_params

def analyze_sternberg_model(sternberg_results):

    ''' Calculates various properties of Krumholz model from fitted phi_cnm
    values.

    '''

    # Get fitted surface density ratios...
    params = [sternberg_results[param] for param in \
              sternberg_results['parameters']]

    rh2_fit, h_sd_fit, f_H2, f_HI = \
           calc_sternberg(params=params,
                          h_sd_extent=sternberg_results['h_sd_fit_range'],
                          return_fractions=True)

    sternberg_results.update({'rh2_fit': rh2_fit,
                              'h_sd_fit' : h_sd_fit,
                              'hi_sd_fit' : (1 - f_H2) * h_sd_fit,
                              'f_H2' : f_H2,
                              'f_HI' : f_HI,})

    return sternberg_results

'''
Bootstrapping functions
'''

def create_cloud_model(av, nhi_background, dgr_background,):

    if nhi_background is None:
        return av

    return av - dgr_background * nhi_background

def create_background_model(av, nhi_cloud, dgr_cloud):

    return av - dgr_cloud * nhi_cloud

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

    filename_extension = global_args['cloud_name'] + '_' + global_args['data_type'] + \
            background_name + \
            bin_name + weights_name + \
            region_name + width_name + avthres_name + \
            intercept_name + error_name + compsub_name + backdgr_name + \
            '_' + hi_range_name + '_' + global_args['radiation_type'] + \
            rotate_cores_name + vary_phi_g_name

    return filename_extension, global_args

def mask_nans(arrays, return_mask=False):

    """ Masks any positions where any array in the list has a NaN.

    Parameters
    ----------
    arrays : tuple
        Tuple of arrays. The mask will be the shape of the first array. The
        last axes of the rest of the arrays will be masked.

    """

    mask = np.zeros(arrays[0].shape, dtype=bool)

    for array in arrays:
        mask[array != array] = 1

    masked_arrays = []
    for array in arrays:
        if isinstance(array, np.ndarray):
            masked_arrays.append(array[~mask])
        else:
            masked_arrays.append(array)

    if return_mask:
        return masked_arrays, mask
    else:
        return masked_arrays

def simulate_noise(av, av_error):

    ''' Simulates noise of Av data

    Possible noise contributions:
        + uncertainty in dust params, tau_353, Beta, T

        + CIB background - section 4.2, Planck 2013

        + background subtraction

        + variation in dust temperature, e.g. difference between radiance and
        tau_353

    '''

    #av_error = (av_error**2 + (0.07 * av)**2)**0.5

    # get bootstrap indices and apply them to each dataset
    np.random.seed()

    # empirical uncertainty from comparison with Schlegel 98
    #sigma_ = 0.003*3.1

    av_noise_sim = np.random.normal(0, scale=av_error)

    # empirical uncertainty from comparison with Schlegel 98
    #av_noise_sim += np.random.normal(0, scale=0.003*3.1)

    return av_noise_sim

def simulate_rescaling(av, scalar=1.0):

    denominator = np.random.uniform(low=1, high=scalar)
    av_rescale = av / denominator

    return av_rescale, denominator

def simulate_background_error(av, scale=1.0):

    # bias from 2MASS image, see lombardi et al. (2009), end of section 2
    av_bias = 0.2
    scale = (scale**2 + (av_bias)**2)**0.5

    background = np.random.normal(0, scale=scale)

    av_background_sim = av + background

    return av_background_sim, background

def simulate_nhi(hi_data, vel_axis, vel_range, vel_range_error,
        hi_data_error=None):

    vel_range_sim = [vel_range[0] + np.random.normal(0, scale=vel_range_error),
                     vel_range[1] + np.random.normal(0, scale=vel_range_error)]

    if hi_data_error is None:
        nhi_sim = myia.calculate_nhi(cube=hi_data,
                                  velocity_axis=vel_axis,
                                  velocity_range=vel_range_sim,
                                  )
        return nhi_sim.ravel(), vel_range_sim
    else:
        nhi_sim, nhi_sim_error = \
            myia.calculate_nhi(cube=hi_data,
                               velocity_axis=vel_axis,
                               velocity_range=vel_range_sim,
                               noise_cube=hi_data_error,
                               return_nhi_error=True,
                               )
        return nhi_sim.ravel(), nhi_sim_error.ravel(), vel_range_sim

def get_rotated_core_indices(core, mask, corename=None, iteration=None):

    import mygeometry as myg

    # Get random angle
    angle = np.random.uniform(low=-1, high=1) * 90.0 # deg

    # rotate vertices
    # ---------------
    vertices = np.array(core['poly_verts']['pixel'])
    center = core['center_pixel']
    vertices_rotated = myg.rotate_wedge(vertices,
                                        center,
                                        angle)

    # Get mask
    # --------
    core_mask = np.logical_not(myg.get_polygon_mask(mask,
                                                    vertices_rotated))

    if 0:
        import matplotlib.pyplot as plt
        plt.close(); plt.clf()
        plt.imshow(core_mask, origin='lower')
        plt.savefig('/d/bip3/ezbc/scratch/' + corename + '_mask' + \
                    str(iteration) + '.png')

    # Map the mask to the unraveled mask and get the indices
    # ------------------------------------------------------
    core_indices = np.where(core_mask == 0)[0]

    return core_indices

def bootstrap_worker(global_args, i):

    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    av = global_args['av']
    av_error = global_args['av_error']
    nhi = global_args['nhi']
    nhi_back = global_args['nhi_back']
    hi_data = global_args['hi_data']
    hi_data_error = global_args['hi_data_error']
    mask = global_args['mask']
    vel_axis = global_args['vel_axis']
    vel_range = global_args['vel_range']
    vel_range_error = global_args['vel_range_error']
    init_guesses = global_args['init_guesses']
    plot_kwargs = global_args['plot_kwargs']
    use_intercept = global_args['use_intercept']
    probabilities = global_args['probabilities']
    av_scalar = global_args['scale_kwargs']['av_scalar']
    intercept_error = global_args['scale_kwargs']['intercept_error']
    model_kwargs = global_args['ss_model_kwargs']
    rotate_cores = global_args['rotate_cores']

    #i = global_args['i']
    #np.random.seed()

    #queue = global_args['queue']

    # Create simulated data
    # -------------------------------------------------------------------------
    # add random noise
    av_sim = av + simulate_noise(av, av_error)

    # rescale the data somewhere between Planck and 2MASS:
    # rescaled Planck = Planck / beta where beta is between 1.0 and 1.4
    av_sim, av_scalar_sim = simulate_rescaling(av_sim, scalar=av_scalar)

    # remove background
    if 0:
        av_sim, av_background_sim = \
                simulate_background_error(av_sim,
                                          scale=intercept_error)
    else:
        av_background_sim = 0.0

    # calculate N(HI)
    if global_args['sim_hi_error']:
        nhi_sim, nhi_sim_error, vel_range_sim = \
            simulate_nhi(hi_data,
                         vel_axis,
                         vel_range,
                         vel_range_error,
                         hi_data_error=hi_data_error,
                         )
    else:
        nhi_sim = nhi

    # Bootstrap data
    # -------------------------------------------------------------------------
    boot_indices = np.random.choice(av.size, size=av.size,)# p=probabilities)

    av_boot = av_sim[boot_indices]
    av_error_boot = av_error[boot_indices]
    nhi_boot = nhi_sim[boot_indices]
    nhi_error_boot = nhi_sim_error[boot_indices]

    if nhi_back is not None:
        nhi_back_boot = nhi_back[boot_indices]
    else:
        nhi_back_boot = None

    # for plotting
    plot_kwargs['bootstrap_num'] = i

    # Fit the bootstrapped data
    # -------------------------------------------------------------------------
    av_model_results = fit_av_model(av_boot,
                            nhi_boot,
                            av_error=av_error_boot,
                            nhi_error=nhi_error_boot,
                            nhi_background=nhi_back_boot,
                            init_guesses=init_guesses,
                            plot_kwargs=plot_kwargs,
                            use_intercept=use_intercept)

    # Calculate N(H2), then HI + H2 surf dens, fit relationship
    # -------------------------------------------------------------------------
    # calculate N(H2) maps
    nh2_image = calculate_nh2(nhi_image=nhi_sim,
                              av_image=av,
                              dgr=av_model_results['dgr_cloud'])
    nh2_image_error = calculate_nh2(nhi_image=nhi_sim_error,
                              av_image=av_error,
                              dgr=av_model_results['dgr_cloud'])

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi_sim,
                               sd_factor=1/1.25)
    hi_sd_image_error = calculate_sd(nhi_sim_error,
                                     sd_factor=1/1.25)

    h2_sd_image = calculate_sd(nh2_image,
                               sd_factor=1/0.625)
    h2_sd_image_error = calculate_sd(nh2_image_error,
                               sd_factor=1/0.625)

    h_sd_image = hi_sd_image + h2_sd_image
    h_sd_image_error = (hi_sd_image_error**2 + h2_sd_image_error**2)**0.5

    # Write ratio between H2 and HI
    rh2_image = h2_sd_image / hi_sd_image
    rh2_image_error = rh2_image * (h2_sd_image_error**2 / h2_sd_image**2 + \
                      h_sd_image_error**2 / h_sd_image**2)**0.5

    model_kwargs = global_args['ss_model_kwargs']['model_kwargs']
    cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']
    cores = global_args['ss_model_kwargs']['cores']

    ss_model_results = {}

    # scale phi_g by relative to the galactic DGR
    # Galactic DGR = 5.3 x 10^-22 cm^2 mag
    # our DGR in units of 10^-20 cm^2 mag
    # Galactic DGR in our units = 5.3 x 10^-22 / 10^-20 = 5.3 x 10^-2 = 0.053
    if 1:
        DGR = av_model_results['dgr_cloud']
        #print 'DGR = ', DGR
        phi_g = DGR / 0.053
        Z = DGR / 0.053
        #print 'phi_g', phi_g
        new_model_kwargs = dict(model_kwargs)
        new_model_kwargs['sternberg_params']['guesses'][2] = phi_g
        #new_model_kwargs['sternberg_params']['guesses'][1] = Z

    # cycle through each core, bootstrapping the pixels
    for core in cores_to_plot:
        if rotate_cores:
            core_indices = get_rotated_core_indices(cores[core],
                                                    mask,
                                                    corename=core,
                                                    iteration=i,
                                                    )
            #if core == 'G168.54-6.22':
            #    print np.sum(core_indices)
        else:
            # grab the indices of the core in the unraveled array
            core_indices = cores[core]['indices']

        if 0:
            assert av[core_indices] == cores[core]['test_pix']

        # Bootstrap core indices
        #core_boot_indices = core_indices[index_ints]
        np.random.seed()
        core_boot_indices = np.random.choice(core_indices.size,
                                             size=core_indices.size,)
        # get bootstrapped pixels
        #h_sd_core = h_sd_image[core_boot_indices]
        #rh2_core = rh2_image[core_boot_indices]
        h_sd_core = h_sd_image[core_indices]
        h_sd_core_error = h_sd_image_error[core_indices]
        rh2_core = rh2_image[core_indices]
        rh2_core_error = rh2_image_error[core_indices]

        # mask negative ratios
        if 1:
            mask_rh2 = (rh2_core < 0) | (np.isnan(rh2_core))
            rh2_core = rh2_core[~mask_rh2]
            rh2_core_error = rh2_core_error[~mask_rh2]
            h_sd_core = h_sd_core[~mask_rh2]
            h_sd_core_error = h_sd_core_error[~mask_rh2]

        G0 = model_kwargs['krumholz_params']['G0'] + \
             np.random.normal(model_kwargs['krumholz_params']['G0_error'])
        ss_model_result = \
            fit_steady_state_models(h_sd_core,
                                    rh2_core,
                                    rh2_error=rh2_core_error,
                                    h_sd_error=h_sd_core_error,
                                    model_kwargs=model_kwargs,
                                    G0=G0,
                                    )

        #print '\npost fit phi_g:'
        #print ss_model_result['sternberg_results']['phi_g']

        # get HI transition result
        add_hi_transition_calc(ss_model_result)
        ss_model_results[core] = ss_model_result

        phi_cnm = ss_model_result['krumholz_results']['phi_cnm']
        #if phi_cnm < 0:
            #print 'core = ', core, ' phi_cnm =', phi_cnm

    # Write results
    # -------------------------------------------------------------------------
    #global_args['init_guesses'] = av_model_result

    # Write results
    mc_results = {}
    mc_results['data_params'] = {'av_background_sim': av_background_sim,
                                 'vel_range_sim': vel_range_sim,
                                 'av_scalar_sim': av_scalar_sim}
    mc_results['av_model_results'] = av_model_results
    mc_results['ss_model_results'] = ss_model_results


    # Plot distribution and fit
    if plot_kwargs['plot_diagnostics']:
    #if 1:
        dgr_cloud = av_model_results['dgr_cloud']
        dgr_background = av_model_results['dgr_background']
        intercept = av_model_results['intercept']

        filename = plot_kwargs['figure_dir'] + \
                   'diagnostics/av_nhi/' + plot_kwargs['filename_base'] + \
                   '_av_vs_nhi_bootstrap' + \
                   '{0:03d}.png'.format(plot_kwargs['bootstrap_num'])
        av_cloud = create_cloud_model(av_boot,
                                     nhi_back_boot,
                                     dgr_background,)
        av_background = create_background_model(av_boot,
                                     nhi_boot,
                                     dgr_cloud,)
        if nhi_back_boot is not None:
            #nhi_total = nhi_boot + nhi_back_boot
            #nhi_total = np.hstack((nhi_boot, nhi_back_boot))
            #av_boot = np.hstack((av_cloud, av_background))
            #av_images = (av_boot, av_cloud, av_background)
            av_images = (av_cloud, av_background)
            #nhi_images = (nhi_total, nhi_boot, nhi_back_boot)
            nhi_images = (nhi_boot, nhi_back_boot)
        else:
            nhi_total = nhi_boot
            av_images = (av_boot,)
            nhi_images = (nhi_total,)

        fit_params = {
                      'dgr_cloud': dgr_cloud,
                      'dgr_cloud_error': (0, 0),
                      'dgr_background': dgr_background,
                      'dgr_background_error': (0, 0),
                      'intercept': intercept,
                      'intercept_error': (0,0),
                      }

        #print('plotting')
        plot_av_vs_nhi(av_images,
                       nhi_images,
                       av_error=av_error_boot,
                       fit_params=fit_params,
                       contour_plot=plot_kwargs['av_nhi_contour'],
                       limits=plot_kwargs['av_nhi_limits'],
                       filename=filename,
                       )
    #queue.put(result)
    return mc_results

def bootstrap_worker_wrapper(args, i):

    import sys
    import traceback

    try:
        output = bootstrap_worker(args, i)
        return output
    except Exception as error:
        # capture the exception and bundle the traceback
        # in a string, then raise new exception with the string traceback
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

def bootstrap_fits(av_data, nhi_image=None, hi_data=None,
        nhi_image_error=None, hi_data_error=None, vel_axis=None,
        vel_range=None, vel_range_error=1, av_error_data=None,
        av_reference=None, nhi_image_background=None, num_bootstraps=100,
        plot_kwargs=None, scale_kwargs=None, use_intercept=True,
        sim_hi_error=False, ss_model_kwargs=None, multiprocess=True,
        rotate_cores=False,):

    import multiprocessing as mp
    import sys
    import traceback

    if av_error_data is None:
        av_error_data = np.ones(av_data.size)

    # mask for nans, arrays will be 1D
    (av, av_error, nhi, nhi_error, nhi_back, ), mask = \
        mask_nans((av_data, av_error_data, nhi_image, nhi_image_error, nhi_image_background),
                   return_mask=True)
    hi_data = hi_data[:, ~mask]
    hi_data_error = hi_data_error[:, ~mask]
    cores = ss_model_kwargs['cores']
    for core in cores:
        if cores[core]['mask'] is not None:
            cores[core]['mask_raveled'] = cores[core]['mask'][~mask]
            cores[core]['indices'] = \
                np.where(cores[core]['mask_raveled'] == 0)[0]
        else:
            cores[core]['indices'] = None

    probabilities = 1.0 / av_error**2
    probabilities /= np.nansum(probabilities)

    # for plotting
    plot_kwargs['num_bootstraps'] = num_bootstraps

    # initialize array for storing output
    boot_results = np.empty((3, num_bootstraps))
    init_guesses = [0.05, 0.05, 0.0] # dgr_cloud, dgr_background, intercept

    # Prep arguments
    global_args = {}
    global_args['av'] = av
    global_args['av_unmasked'] = av_data
    global_args['av_error'] = av_error
    global_args['nhi'] = nhi
    global_args['rotate_cores'] = rotate_cores
    global_args['mask'] = mask
    global_args['nhi_back'] = nhi_back
    global_args['init_guesses'] = init_guesses
    global_args['plot_kwargs'] = plot_kwargs
    global_args['use_intercept'] = use_intercept
    global_args['probabilities'] = probabilities
    global_args['scale_kwargs'] = scale_kwargs
    global_args['sim_hi_error'] = sim_hi_error
    global_args['ss_model_kwargs'] = ss_model_kwargs
    if sim_hi_error:
        global_args['hi_data'] = hi_data
        global_args['vel_axis'] = vel_axis
        global_args['vel_range'] = vel_range
        global_args['vel_range_error'] = vel_range_error
    else:
        global_args['hi_data'] = None
        global_args['vel_axis'] = None
        global_args['vel_range'] = None
        global_args['vel_range_error'] = None

    # Prep multiprocessing
    queue = mp.Queue(10)
    pool = mp.Pool()
    processes = []

    # bootstrap
    if multiprocess:
        for i in xrange(num_bootstraps):
            processes.append(pool.apply_async(bootstrap_worker_wrapper,
                                              args=(global_args,i,)))
        pool.close()
        pool.join()

        # Get the results
        mc_results = collect_bootstrap_results(processes, ss_model_kwargs)

        for i in xrange(len(processes)):
            result = processes[i].get()
            boot_results[:, i] = result['av_model_results'].values()
    else:
        for i in xrange(num_bootstraps):
            processes.append(bootstrap_worker(global_args, i))

        mc_results = collect_bootstrap_results(processes, ss_model_kwargs,
                                               multiprocess=False)

        for i in xrange(len(processes)):
            result = processes[i]
            boot_results[:, i] = result['av_model_results'].values()

    return boot_results, mc_results

def residual_worker(global_args, core):

    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    rh2 = global_args['rh2']
    rh2_error = global_args['rh2_error']
    h_sd = global_args['h_sd']
    nboot = global_args['nboot']
    av = global_args['av']
    av_error = global_args['av_error']
    nhi = global_args['nhi']
    nhi_back = global_args['nhi_back']
    hi_data = global_args['hi_data']
    mask = global_args['mask']
    vel_axis = global_args['vel_axis']
    vel_range = global_args['vel_range']
    vel_range_error = global_args['vel_range_error']
    init_guesses = global_args['init_guesses']
    plot_kwargs = global_args['plot_kwargs']
    use_intercept = global_args['use_intercept']
    probabilities = global_args['probabilities']
    av_scalar = global_args['scale_kwargs']['av_scalar']
    intercept_error = global_args['scale_kwargs']['intercept_error']
    model_kwargs = global_args['ss_model_kwargs']
    rotate_cores = global_args['rotate_cores']

    model_kwargs = global_args['ss_model_kwargs']['model_kwargs']
    cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']
    cores = global_args['ss_model_kwargs']['cores']

    # cycle through each core, bootstrapping the pixels
    if rotate_cores:
        core_indices = get_rotated_core_indices(cores[core],
                                                mask,
                                                corename=core,
                                                iteration=i,
                                                )
        #if core == 'G168.54-6.22':
        #    print np.sum(core_indices)
    else:
        # grab the indices of the core in the unraveled array
        core_indices = cores[core]['indices']

    h_sd_core = h_sd[core_indices]
    rh2_core = rh2[core_indices]
    rh2_core_error = rh2_error[core_indices]

    # mask negative ratios
    mask_rh2 = (rh2_core < 0) | (np.isnan(rh2_core))
    rh2_core = rh2_core[~mask_rh2]
    rh2_core_error = rh2_core_error[~mask_rh2]
    h_sd_core = h_sd_core[~mask_rh2]

    ss_model_result = \
        fit_steady_state_models(h_sd_core,
                                rh2_core,
                                rh2_error=rh2_core_error,
                                model_kwargs=model_kwargs,
                                bootstrap_residuals=True,
                                nboot=nboot,
                                )

    # get HI transition result
    add_hi_transition_calc(ss_model_result)

    return ss_model_result, core

def residual_worker_wrapper(args, core):

    import sys
    import traceback

    try:
        output = residual_worker(args, core)
        return output
    except Exception as error:
        # capture the exception and bundle the traceback
        # in a string, then raise new exception with the string traceback
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

def bootstrap_residuals(av_data, nhi_image=None, hi_data=None, vel_axis=None,
        vel_range=None, vel_range_error=1, av_error_data=None,
        av_reference=None, nhi_image_background=None, num_bootstraps=100,
        plot_kwargs=None, scale_kwargs=None, use_intercept=True,
        sim_hi_error=False, ss_model_kwargs=None, multiprocess=True,
        rotate_cores=False,):

    import multiprocessing as mp
    import sys
    import traceback
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    if av_error_data is None:
        av_error_data = np.ones(av_data.size)

    # mask for nans, arrays will be 1D
    (av, av_error, nhi, nhi_back, ), mask = \
        mask_nans((av_data, av_error_data, nhi_image, nhi_image_background),
                   return_mask=True)
    hi_data = hi_data[:, ~mask]
    cores = ss_model_kwargs['cores']
    for core in cores:
        if cores[core]['mask'] is not None:
            cores[core]['mask_raveled'] = cores[core]['mask'][~mask]
            cores[core]['indices'] = \
                np.where(cores[core]['mask_raveled'] == 0)[0]
        else:
            cores[core]['indices'] = None

    probabilities = 1.0 / av_error**2
    probabilities /= np.nansum(probabilities)

    # for plotting
    plot_kwargs['num_bootstraps'] = num_bootstraps

    # initialize array for storing output
    boot_results = np.empty((3, num_bootstraps))
    init_guesses = [0.05, 0.05, 0.0] # dgr_cloud, dgr_background, intercept

    # calculate N(HI)
    nhi = myia.calculate_nhi(cube=hi_data,
                             velocity_axis=vel_axis,
                             velocity_range=vel_range,
                             )

    # Fit the data
    # -------------------------------------------------------------------------
    av_model_results = fit_av_model(av,
                            nhi,
                            av_error=av_error,
                            nhi_background=nhi_back,
                            init_guesses=init_guesses,
                            plot_kwargs=plot_kwargs,
                            use_intercept=use_intercept)

    # Calculate N(H2), then HI + H2 surf dens, fit relationship
    # -------------------------------------------------------------------------
    # calculate N(H2) maps
    nh2_image = calculate_nh2(nhi_image=nhi,
                              av_image=av,
                              dgr=av_model_results['dgr_cloud'])
    nh2_image_error = calculate_nh2(nhi_image=nhi,
                              av_image=av_error,
                              dgr=av_model_results['dgr_cloud'])

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi,
                               sd_factor=1/1.25)
    hi_sd_image_error = 0.02 * hi_sd_image

    h2_sd_image = calculate_sd(nh2_image,
                               sd_factor=1/0.625)
    h2_sd_image_error = calculate_sd(nh2_image_error,
                               sd_factor=1/0.625)

    h_sd_image = hi_sd_image + h2_sd_image
    h_sd_image_error = (hi_sd_image_error**2 + h2_sd_image_error**2)**0.5

    # Write ratio between H2 and HI
    rh2_image = h2_sd_image / hi_sd_image
    rh2_image_error = rh2_image * (h2_sd_image_error**2 / h2_sd_image**2 + \
                      h_sd_image_error**2 / h_sd_image**2)**0.5

    # Prep arguments
    global_args = {}
    global_args['av'] = av
    global_args['av_unmasked'] = av_data
    global_args['av_error'] = av_error
    global_args['rh2'] = rh2_image
    global_args['rh2_error'] = rh2_image_error
    global_args['h_sd'] = h_sd_image
    global_args['nhi'] = nhi
    global_args['rotate_cores'] = rotate_cores
    global_args['mask'] = mask
    global_args['nboot'] = num_bootstraps
    global_args['nhi_back'] = nhi_back
    global_args['init_guesses'] = init_guesses
    global_args['plot_kwargs'] = plot_kwargs
    global_args['use_intercept'] = use_intercept
    global_args['probabilities'] = probabilities
    global_args['scale_kwargs'] = scale_kwargs
    global_args['sim_hi_error'] = sim_hi_error
    global_args['ss_model_kwargs'] = ss_model_kwargs
    if sim_hi_error:
        global_args['hi_data'] = hi_data
        global_args['vel_axis'] = vel_axis
        global_args['vel_range'] = vel_range
        global_args['vel_range_error'] = vel_range_error
        print 'vel_range_error', vel_range_error
    else:
        global_args['hi_data'] = None
        global_args['vel_axis'] = None
        global_args['vel_range'] = None
        global_args['vel_range_error'] = None
    cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']

    # Prep multiprocessing
    queue = mp.Queue(10)
    pool = mp.Pool()
    processes = []

    # bootstrap
    if multiprocess:
        for core in cores_to_plot:
            processes.append(pool.apply_async(residual_worker_wrapper,
                                              args=(global_args,core)))
        pool.close()
        pool.join()

        # Get the results
        resid_results = collect_residual_results(processes, ss_model_kwargs)

        for i in xrange(len(processes)):
            result = processes[i].get()
    else:
        for i in xrange(num_bootstraps):
            processes.append(residual_worker(global_args, i))

        resid_results = collect_residual_results(processes, ss_model_kwargs,
                                               multiprocess=False)

    return resid_results

def collect_residual_results(processes, ss_model_kwargs, multiprocess=True):

    empty = lambda: np.empty(len(processes))
    mc_analysis = {}
    mc_analysis['cores'] = {}

    for i in xrange(len(processes)):
        #result = queue.get())

        if multiprocess:
            result = processes[i].get()
        else:
            result = processes[i]

        model_result, core = result
        mc_analysis[core] = {}

        for model in model_result:
            mc_analysis[core][model] = {}
            model_dict = \
                mc_analysis[core][model]
            for param in model_result[model]:
                param_result = \
                    model_result[model][param]
                model_dict[param] = param_result

    return mc_analysis

def collect_bootstrap_results(processes, ss_model_kwargs, multiprocess=True):

    empty = lambda: np.empty(len(processes))
    mc_results = {}
    mc_results['av_model_results'] = {}
    mc_results['ss_model_results'] = {}
    mc_results['av_model_results']['dgr'] = empty()
    mc_results['ss_model_results']['cores'] = {}
    for core in ss_model_kwargs['cores_to_plot']:
        mc_results['ss_model_results']['cores'][core] = {}
        core_dict = mc_results['ss_model_results']['cores'][core]
        core_dict['krumholz_results'] = \
                {'phi_cnm': empty(),
                 'Z': empty(),
                 'phi_mol': empty(),
                 'hi_transition': empty(),
                 }
        core_dict['sternberg_results'] = \
                {'alphaG': empty(),
                 'Z': empty(),
                 'phi_g': empty(),
                 'hi_transition': empty(),
                 }
    mc_results['data_params'] = \
            {'av_background_sim': empty(),
             'vel_range_sim': np.empty((len(processes), 2)),
             'av_scalar_sim': empty()}

    for i in xrange(len(processes)):
        #result = queue.get())

        if multiprocess:
            result = processes[i].get()
        else:
            result = processes[i]

        for data_param in mc_results['data_params']:
            mc_results['data_params'][data_param][i] = \
                result['data_params'][data_param]
        mc_results['av_model_results']['dgr'][i] = \
            result['av_model_results']['dgr_cloud']

        for core in result['ss_model_results']:
            for model in result['ss_model_results'][core]:
                for param in result['ss_model_results'][core][model]:
                    param_result = \
                        result['ss_model_results'][core][model][param]
                    model_dict = \
                        mc_results['ss_model_results']['cores'][core][model]
                    model_dict[param][i] = param_result

    return mc_results

def scale_av_with_refav(av_data, av_reference, av_error_data):

    import scipy as sp


    if av_reference is not None:
        # Mask AV for Av > 10 mag. The 2MASS Av will saturate at high Av thus is
        # unreliable.

        nan_mask = (np.isnan(av_reference) | \
                    np.isnan(av_data) | \
                    np.isnan(av_error_data) | \
                    (av_reference > 10.0))
        p, V = np.polyfit(av_reference[~nan_mask], av_data[~nan_mask], deg=1,
                       #w=1.0/av_error_data[~nan_mask]**2,
                       cov=True,
                       )
        av_scalar, intercept = p
        intercept_error = V[1, 1]

        # Perform residual bootstrapping to get errors on intercept
        bootindex = sp.random.random_integers
        nboot = 1000
        x_fit = av_reference[~nan_mask]
        y_fit = p[0] * x_fit  + p[1]
        residuals = av_data[~nan_mask] - y_fit
        fits = np.empty((2, nboot))

        weights = np.abs(1.0 / av_error_data[~nan_mask])
        weights /= np.nansum(weights)

        for i in xrange(nboot): # loop over n bootstrap samples from the resids
            if 0:
                boot_indices = bootindex(0,
                                         len(residuals)-1,
                                         len(residuals))
                residuals_bootstrapped = residuals[boot_indices]
                fits[:, i] = sp.polyfit(av_reference[~nan_mask],
                                        y_fit + residuals_bootstrapped,
                                        deg=1)
            else:
                boot_indices = bootindex(0,
                                         len(x_fit)-1,
                                         len(x_fit))
                x = x_fit[boot_indices] + np.random.normal(0, 0.4,
                                                           size=x_fit.size)
                y = av_data[~nan_mask][boot_indices] + \
                        np.random.normal(0,
                                av_error_data[~nan_mask][boot_indices])

                fits[:, i] = np.polyfit(x,
                                        y,
                                        deg=1,)

        intercept_error = np.std(fits[1])
        av_scalar_error = np.std(fits[0])

    kwargs = {}
    kwargs['av_scalar'] = av_scalar
    kwargs['av_scalar_error'] = av_scalar_error
    kwargs['intercept'] = intercept
    kwargs['intercept_error'] = intercept_error

    print 'Reference Av characteristcs:'
    print '\t', av_scalar, av_scalar_error
    print '\t', intercept, intercept_error
    print ''

    return kwargs

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

        print 'min error', hi_width_error_min, 'km/s'
        print 'max error', hi_width_error_max, 'km/s'
        print 'avg error', hi_width_error, 'km/s'

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

def calc_param_errors(results_dict):

    import mystats

    boot_result = results_dict['boot_result']
    global_args = results_dict['global_args']

    dgr_cloud, dgr_background, intercept = np.mean(boot_result, axis=1)
    dgr_cloud_error, dgr_background_error, intercept_error = \
                np.std(boot_result, axis=1)

    # Calculate conf interval
    dgrs = boot_result[0]
    dgr_cloud, dgr_cloud_error = mystats.calc_cdf_error(dgrs,
                                                        alpha=0.32)

    if global_args['use_background']:
        dgrs = boot_result[1]
        dgr_background_error = (mean - conf_int_a[0], conf_int_a[1] - mean)
    else:
        dgr_background_error = (0.0, 0.0)
        dgr_background = 0.0

    if global_args['use_intercept']:
        intercepts = boot_result[2]
        intercept, intercept_error = mystats.calc_cdf_error(intercepts,
                                                            alpha=0.32)
    else:
        intercept_error = (0.0, 0.0)
        intercept = 0.0

    fit_params = {
                  'dgr_cloud': dgr_cloud,
                  'dgr_cloud_error': dgr_cloud_error,#, dgr_cloud_error),
                  'dgr_background': dgr_background,
                  'dgr_background_error': dgr_background_error,
                  'intercept': intercept,
                  'intercept_error': intercept_error,
                  }

    return fit_params

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
        results['data']['co_data'], results['data']['co_header'] = \
                fits.getdata(results['filenames']['co_filename'], header=True)

    return results

def save_results(results_dict, filename, write_fits=False):

    import pickle
    from astropy.io import fits

    if not write_fits:
        results_dict['data']['av_data'] = None
        results_dict['data']['av_error_data'] = None
        results_dict['data']['av_data_ref'] = None
        results_dict['data']['hi_data'] = None
        results_dict['data']['co_data'] = None

    with open(filename, 'wb') as output:
        pickle.dump(results_dict, output)
    output.close()

def get_model_fit_kwargs(cloud_name, vary_phi_g=False):

    '''

    '''
    vary_alphaG = True # Vary alphaG in S+14 fit?
    vary_Z = False # Vary metallicity in S+14 fit?
    vary_phi_g = vary_phi_g # Vary phi_g in S+14 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    guesses=[10.0, 1.0, 1] # Guesses for (alphaG, Z, phi_g)
    h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

    # Monte carlo results file bases
    results_filename = '/d/bip3/ezbc/multicloud/' + \
                       '/data/python_output/' + \
                       'monte_carlo_results/' + \
                       cloud_name + '_mc_results_' + \
                       'planck' + '_'

    sternberg_params = {}
    sternberg_params['param_vary'] = [vary_alphaG, vary_Z, vary_phi_g]
    sternberg_params['error_method'] = error_method
    sternberg_params['alpha'] = alpha
    sternberg_params['guesses'] = guesses
    sternberg_params['h_sd_fit_range'] = h_sd_fit_range
    sternberg_params['results_filename'] = results_filename
    sternberg_params['parameters'] = ['alphaG', 'Z', 'phi_g']

    # Krumholz Parameters
    # --------------------
    vary_phi_cnm = True # Vary phi_cnm in K+09 fit?
    vary_Z = False # Vary metallicity in K+09 fit?
    vary_phi_mol = False # Vary phi_mol in K+09 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    guesses=[8.0, 1.0, 10.0] # Guesses for (phi_cnm, Z, phi_mol)
    h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

    krumholz_params = {}
    krumholz_params['param_vary'] = [vary_phi_cnm, vary_Z, vary_phi_mol]
    krumholz_params['error_method'] = error_method
    krumholz_params['alpha'] = alpha
    krumholz_params['guesses'] = guesses
    krumholz_params['h_sd_fit_range'] = h_sd_fit_range
    krumholz_params['results_filename'] = results_filename
    krumholz_params['parameters'] = ['phi_cnm', 'Z', 'phi_mol']
    if cloud_name == 'taurus':
        G0 = 0.6
        G0_error = 0.1
    elif cloud_name == 'california':
        G0 = 1.0
        G0_error = 0.2
    elif cloud_name == 'perseus':
        G0 = 0.7
        G0_error = 0.2
    krumholz_params['G0'] = G0
    krumholz_params['G0_error'] = G0_error

    results = {}
    results['results_filename'] = results_filename
    results['krumholz_params'] = krumholz_params
    results['sternberg_params'] = sternberg_params

    return results

def add_coldens_images(data_products, mc_analysis):

    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    nhi = data_products['nhi']
    av = data_products['av']

    # calculate N(H2) maps
    nh2 = calculate_nh2(nhi_image=nhi,
                        av_image=av,
                        dgr=mc_analysis['dgr'])

    # convert to column density to surface density
    hi_sd = calculate_sd(nhi,
                         sd_factor=1/1.25)

    h2_sd = calculate_sd(nh2,
                         sd_factor=1/0.625)

    h_sd = hi_sd + h2_sd

    # Write ratio between H2 and HI
    rh2 = h2_sd / hi_sd

    data_products['nh2'] = nh2
    data_products['h2_sd'] = h2_sd
    data_products['hi_sd'] = hi_sd
    data_products['h_sd'] = h_sd
    data_products['rh2'] = rh2

def calc_mc_analysis(mc_results, resid_results):

    import mystats

    mc_analysis = {}
    core_analysis = {}

    # Calculate conf intervals on parameters
    # -------------------------------------------------------------------------
    # DGR
    dgrs = mc_results['av_model_results']['dgr']
    dgr, dgr_error = mystats.calc_cdf_error(dgrs,
                                            alpha=0.32)

    # Model params fit for each core:
    for core in mc_results['ss_model_results']['cores']:
        core_analysis[core] = {}
        for model in mc_results['ss_model_results']['cores'][core]:
            params = {}
            core_analysis[core][model] = {}
            for param in mc_results['ss_model_results']['cores'][core][model]:
                param_results = \
                    mc_results['ss_model_results']['cores'][core][model][param]
                param_value, param_error = \
                    mystats.calc_cdf_error(param_results,
                                           alpha=0.32)

                # add fit error
                if param != 'hi_transition' and resid_results is not None:
                    fit_error = \
                        np.array(resid_results[core][model][param + '_error'])

                    param_error = np.array(param_error)
                    #print ''
                    #print 'before', param, param_error / param_value
                    #print 'fit_error', fit_error / param_value
                    #print 'mc error', param_error / param_value
                    if 0:
                        param_error = np.sqrt(fit_error**2 + \
                                              np.array(param_error)**2)
                    #print 'after', param, param_error / param_value


                core_analysis[core][model][param] = param_value
                core_analysis[core][model][param + '_error'] = param_error

                params[param] = param_value

                #if param == 'phi_g':
                    #print '\nanalysis'
                    #print param_value
                    #print param_error


                if 0:
                    if param == 'phi_g' or param == 'alphaG':
                        import matplotlib.pyplot as plt
                        import mystats
                        plt.close(); plt.clf()
                        cdf = mystats.calc_cdf(param_results)
                        plt.plot(np.sort(param_results), cdf)
                        plt.xlabel(param)
                        plt.savefig('/d/bip3/ezbc/scratch/'+core+\
                                    '_' + param + '_cdf.png')

                    print('core ' + core + ': ' + param  + \
                          ' = {0:.2f}'.format(param_value) + \
                          '+{0:.2f} / -{1:.2f}'.format(*param_error))

            if 'sternberg' in model:
                model_fits = calc_sternberg(params,
                                          h_sd_extent=(0.001, 1000),
                                          return_fractions=False,
                                          return_hisd=True)
            elif 'krumholz' in model:
                model_fits = calc_krumholz(params,
                                          h_sd_extent=(0.001, 1000),
                                          return_fractions=False,
                                          return_hisd=True)

            core_analysis[core][model]['rh2_fit'] = model_fits[0]
            core_analysis[core][model]['hsd_fit'] = model_fits[1]
            core_analysis[core][model]['hisd_fit'] = model_fits[2]


    mc_analysis = {'dgr': dgr,
                   'dgr_error': dgr_error,
                   'cores': core_analysis,
                   }

    return mc_analysis

def add_results_analysis(results_dict):

    # calculate statistics of bootstrapped model values
    results_dict['mc_analysis'] = calc_mc_analysis(results_dict['mc_results'],
                                            results_dict['resid_mc_results'])

    if 0:
        for core in results_dict['mc_analysis']['cores']:
            for key in results_dict['mc_analysis']['cores'][core]:
                print results_dict['mc_analysis']['cores'][core][key]

    # derive N(H2), N(H) etc...
    add_coldens_images(results_dict['data_products'],
                       results_dict['mc_analysis'])

def run_cloud_analysis(global_args,):

    from astropy.io import fits
    from myimage_analysis import calculate_nhi, calc_region_mask
    import myimage_analysis as myia
    from mycoords import make_velocity_axis
    from mystats import calc_symmetric_error, calc_logL
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

    # define filenames
    prop_filename = property_dir + \
       cloud_name + '_global_properties.txt'
    hi_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres.fits'
    hi_error_filename = hi_dir + \
       cloud_name + '_hi_galfa_cube_regrid_planckres_noise.fits'
    co_filename = co_dir + \
       cloud_name + '_co_cfa_cube_regrid_planckres.fits'

    if cloud_name == 'perseus' and data_type == 'lee12':
        av_filename = av_dir + \
           cloud_name + '_av_lee12_iris_regrid_planckres.fits'
        av_error_filename = None
        av_error = 0.1
        if background_subtract:
            av_background = 0.5
        else:
            av_background = None
    if data_type == 'planck':
        av_filename = av_dir + \
           cloud_name + '_av_planck_tau353_5arcmin.fits'
        av_error_filename = av_dir + \
           cloud_name + '_av_error_planck_tau353_5arcmin.fits'
        av_error = None
        if 0:
            av_error_filename = None
            av_error = 1
        av_background = None
        av_ref_filename = av_dir + \
           cloud_name + '_av_k09_regrid_planckres.fits'
    if cloud_name == 'perseus' and data_type == 'planck_lee12mask':
        av_filename = av_dir + \
           cloud_name + '_av_planck_tau353_5arcmin_lee12mask.fits'
        av_error_filename = av_dir + \
           cloud_name + '_av_error_planck_tau353_5arcmin.fits'
        av_error = None
        av_background = None
    if data_type == 'k09':
        av_filename = av_dir + \
           cloud_name + '_av_k09_regrid_planckres.fits'

        av_error_filename = None
        av_error = 0.4

        av_background = 0.0

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

    # Load data
    if global_args['bin_image']:
        av_filename = av_filename.replace('.fits', '_bin.fits')
        if av_error_filename is not None:
            av_error_filename = av_error_filename.replace('.fits', '_bin.fits')
        hi_filename = hi_filename.replace('.fits', '_bin.fits')
        av_nhi_contour = False
    else:
        av_nhi_contour = True

    av_data, av_header = fits.getdata(av_filename, header=True)
    av_data_ref, av_header = fits.getdata(av_ref_filename, header=True)
    if av_error_filename is not None:
        av_error_data, av_error_header = fits.getdata(av_error_filename,
                                                      header=True)
    else:
        av_error_data = av_error * np.ones(av_data.shape)

    # mask data
    region_filename = region_dir + 'multicloud_divisions.reg'
    region_mask = calc_region_mask(region_filename,
                                   av_data,
                                   av_header,
                                   region_name=global_args['region_name'])

    av_data[region_mask] = np.nan
    av_data_ref[region_mask] = np.nan

    # Scale the data to the 2MASS K+09 data
    scale_kwargs = scale_av_with_refav(av_data, av_data_ref, av_error_data)
    av_data_backsub = av_data - scale_kwargs['intercept']
    avg_scalar = (scale_kwargs['av_scalar'] + 1) / 2.0
    av_data_backsub_scaled = av_data_backsub / avg_scalar

    if 0:
        print('\n\tSubtracting background of ' + \
              str(scale_kwargs['intercept']) + \
              ' from ' + cloud_name)

    # Load HI and CO cubes
    hi_data, hi_header = fits.getdata(hi_filename, header=True)
    hi_data_error = fits.getdata(hi_error_filename)
    co_data, co_header = fits.getdata(co_filename, header=True)

    hi_data[:, region_mask] = np.nan
    co_data[:, region_mask] = np.nan

    hi_vel_axis = make_velocity_axis(hi_header)
    co_vel_axis = make_velocity_axis(co_header)

    # Derive N(HI)
    # -------------------------------------------------------------------------
    # get fit kwargs
    gauss_fit_kwargs, ncomps_in_cloud = get_gauss_fit_kwargs(global_args)

    # derive spectra or load
    spectra_filename = results_dir + 'spectra/' + global_args['cloud_name'] + \
            '_spectra.pickle'
    load_spectra = myio.check_file(spectra_filename,
                                   clobber=global_args['clobber_spectra'])
    if load_spectra:
        hi_spectrum, hi_std_spectrum, co_spectrum = \
                myio.load_pickle(spectra_filename)
    else:
        hi_spectrum = myia.calc_spectrum(hi_data)
        hi_std_spectrum = myia.calc_spectrum(hi_data, statistic=np.nanstd)
        co_spectrum = myia.calc_spectrum(co_data)
        myio.save_pickle(spectra_filename,
                         (hi_spectrum, hi_std_spectrum, co_spectrum))

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
    filename = plot_kwargs['figure_dir'] + \
               'spectra/' + plot_kwargs['filename_base'] + \
               '_spectra.png'
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

    if 0:
        print('\n\tVelocity range = ' + \
              '{0:.1f} to {1:.1f}'.format(*velocity_range))

    # use the vel range to derive N(HI)
    nhi_image, nhi_image_error = \
        calculate_nhi(cube=hi_data,
                      velocity_axis=hi_vel_axis,
                      velocity_range=velocity_range,
                      noise_cube=hi_data_error,
                      return_nhi_error=True,
                      )
    nhi_image_background = calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=(-100,velocity_range[0]),
                              )
    nhi_image_background += calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=(velocity_range[1],100),
                              )

    # mask for erroneous pixels
    nhi_image[nhi_image < 0] = np.nan
    nhi_image_error[nhi_image_error < 0] = np.nan
    nhi_image_background[nhi_image_background < 0] = np.nan

    if not global_args['use_background']:
        nhi_image_background = None


    # Write filenames
    filenames = {
                 'region_filename': region_filename,
                 'av_filename': av_filename,
                 'av_error_filename': av_error_filename,
                 'av_ref_filename': av_ref_filename,
                 'hi_filename': hi_filename,
                 'hi_error_filename': hi_error_filename,
                 'co_filename': co_filename,
                 'results_dir': results_dir,
                 'figure_dir': figure_dir,
                 }

    # Collect data
    data = {
            'av_data': av_data,
            'av_error_data': av_error_data,
            'av_data_ref': av_data_ref,
            'hi_data': hi_data,
            'hi_data_error': hi_data_error,
            'hi_vel_axis': hi_vel_axis,
            'co_data': co_data,
            'co_vel_axis': co_vel_axis,
            'av_header': av_header,
            'av_error_header': av_error_header,
            'hi_header': hi_header,
            'co_header': co_header,
            }

    # Collect data products
    data_products = {
                     'av_data_backsub': av_data_backsub,
                     'av': av_data_backsub_scaled,
                     'nhi': nhi_image,
                     'nhi_error': nhi_image_error,
                     'nhi_image_background': nhi_image_background,
                     'region_mask': region_mask,
                     'scale_kwargs': scale_kwargs,
                     'hi_spectrum': hi_spectrum,
                     'hi_std_spectrum': hi_std_spectrum,
                     'co_spectrum': co_spectrum,
                     'hi_range_kwargs': hi_range_kwargs,
                     }

    # Get model fitting params
    model_fitting = get_model_fit_kwargs(cloud_name,
                                         vary_phi_g=global_args['vary_phi_g'])
    model_fitting['sternberg_params']['radiation_type'] = \
            global_args['radiation_type']

    # Get cores params
    cores = get_core_properties(data, cloud_name)

    # Calculate average dust temperature of each core
    #cores = add_core_temps(data, )

    # get the cores in the cloud
    cores_to_plot = get_cores_to_plot()
    cores_to_plot = trim_cores_to_plot(cores, cores_to_plot)

    if 0:
        import sys
        sys.exit()

    global_args['ss_model_kwargs'] = {}
    global_args['ss_model_kwargs']['cores'] = cores
    global_args['ss_model_kwargs']['cores_to_plot'] = cores_to_plot
    global_args['ss_model_kwargs']['model_kwargs'] = model_fitting

    # Bootstrap data
    # -------------------------------------------------------------------------

    # crop hi_data to be a reasonable size
    hi_data_crop, hi_vel_axis_crop = myia.crop_cube(hi_data,
                                                    hi_vel_axis,
                                                    [-20, 30])
    hi_data_error_crop, hi_vel_axis_crop = myia.crop_cube(hi_data_error,
                                                    hi_vel_axis,
                                                    [-20, 30])

    bootstrap_filename = results_dir + filename_base + '_bootresults.npy'
    results_filename = results_dir + \
               'bootstrap_results/' + filename_base + \
               '_bootstrap_results.pickle'

    # Bootstrap residuals of best fitting models
    if global_args['bootstrap_fit_residuals']:
        print('\n\tBeginning residual bootstrapping...')
        resid_mc_results = \
            bootstrap_residuals(av_data_backsub,
                           nhi_image=nhi_image,
                           nhi_image_error=nhi_image_error,
                           av_error_data=av_error_data,
                           nhi_image_background=nhi_image_background,
                           plot_kwargs=plot_kwargs,
                           hi_data=hi_data_crop,
                           hi_data_error=hi_data_error_crop,
                           vel_axis=hi_vel_axis_crop,
                           vel_range=velocity_range,
                           vel_range_error=2,
                           av_reference=av_data_ref,
                           use_intercept=global_args['use_intercept'],
                           num_bootstraps=global_args['num_resid_bootstraps'],
                           scale_kwargs=scale_kwargs,
                           sim_hi_error=global_args['sim_hi_error'],
                           ss_model_kwargs=global_args['ss_model_kwargs'],
                           multiprocess=global_args['multiprocess'],
                           rotate_cores=global_args['rotate_cores'],
                           )
    else:
        resid_mc_results = None

    print('\n\tBeginning bootstrap monte carlo...')
    # Perform bootsrapping
    boot_result, mc_results = \
        bootstrap_fits(av_data_backsub,
                       nhi_image=nhi_image,
                       nhi_image_error=nhi_image_error,
                       av_error_data=av_error_data,
                       nhi_image_background=nhi_image_background,
                       plot_kwargs=plot_kwargs,
                       hi_data=hi_data_crop,
                       hi_data_error=hi_data_error_crop,
                       vel_axis=hi_vel_axis_crop,
                       vel_range=velocity_range,
                       vel_range_error=hi_range_error,
                       av_reference=av_data_ref,
                       use_intercept=global_args['use_intercept'],
                       num_bootstraps=global_args['num_bootstraps'],
                       scale_kwargs=scale_kwargs,
                       sim_hi_error=global_args['sim_hi_error'],
                       ss_model_kwargs=global_args['ss_model_kwargs'],
                       multiprocess=global_args['multiprocess'],
                       rotate_cores=global_args['rotate_cores'],
                       )
    np.save(bootstrap_filename, boot_result)

    results_dict = {'boot_result': boot_result,
                    'data': data,
                    'data_products': data_products,
                    'global_args': global_args,
                    'plot_kwargs': plot_kwargs,
                    'filenames': filenames,
                    'mc_results': mc_results,
                    'resid_mc_results': resid_mc_results,
                    }

    print('\n\tSaving results...')
    save_results(results_dict, global_args['results_filename'])
    #results_dict = load_results(global_args['results_filename'])

    return results_dict

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
    add_results_analysis(results_dict)

    # calculate errors on dgrs and intercept
    results_dict['params_summary'] = calc_param_errors(results_dict)

    return results_dict

def main():

    import itertools

    results = {}

    clouds = (
              'california',
              'taurus',
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
                'clobber_spectra': False,
                'use_background': permutation[10],
                #'num_bootstraps': 10000,
                'num_bootstraps': 10,
                'num_resid_bootstraps': 100,
                'bootstrap_fit_residuals': False,
                'hi_range_calc': permutation[11],
                'sim_hi_error': True,
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

            print('\n\tPlotting')
            plot_results(results[global_args['cloud_name']])

    plot_multicloud_results(results)

if __name__ == '__main__':
    main()



