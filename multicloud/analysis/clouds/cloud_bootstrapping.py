#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import myimage_analysis as myia
from multiprocessing.queues import Queue
import mygeometry as myg


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

    myplt.set_color_cycle(num_colors=4)

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
                        color='r',
                        marker='s',
                        linestyle='None',
                        label=label,
                        alpha=0.5,
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

        myplt.set_color_cycle(num_colors=2, cmap_limits=[0.5, 0.8])
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
                    #color='r',
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
                color='#6B47B2',
                linestyle='--',
                linewidth=2,
                alpha=0.8,
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
        title='', base=10.0, normalize=False):

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


    c_cycle = myplt.set_color_cycle(num_colors=3, cmap_limits=[0, 0.8])

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

        ax.errorbar(
            bin_centers,
            n_av,
            #yerr = n**0.5,
            marker = '',
            label=r'A$_V$',
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
            color=c_cycle[1],
            drawstyle = 'steps-mid'
        )
        ax.errorbar(
            bin_centers,
            n_av_nh2,
            #yerr = n**0.5,
            label=r'$2\,N($H$_2) \times$ DGR',
            marker = '',
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
        model_analysis=None, limits=None, scale=('linear', 'linear'),
        filename=None, show_params=False, levels=5):

    # Import external modules
    import numpy as np
    import matplotlib.pyplot as plt
    import myplotting as myplt
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Determine size of figure and number of grids
    # --------------------------------------------
    n = int(np.ceil(len(core_names)**0.5))
    if n**2 - n > len(core_names):
        nrows = n - 1
        ncols = 2
        y_scaling = 1.0 - 1.0 / n
    else:
        nrows, ncols = n, n
        y_scaling = 1.0

    n = len(core_names)
    ncols = 2
    nrows = (n + 1) / ncols
    y_scaling = nrows / float(ncols)

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 3)]
    font_scale = 9

    figsize = (3.6, 3.6*y_scaling)

    # Create figure instance
    fig = plt.figure(figsize=figsize)

    if n == 1:
        n = 2

    c_cycle = myplt.set_color_cycle(num_colors=3, cmap_limits=[0, 0.8])

    axes = AxesGrid(fig, (1,1,1),
                    nrows_ncols=(nrows, ncols),
                    ngrids=n,
                    axes_pad=0.1,
                    aspect=False,
                    label_mode='L',
                    share_all=True,
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
            l1 = myplt.scatter_contour(h_sd_nonans.ravel(),
                                 hi_sd_nonans.ravel(),
                                 threshold=2,
                                 log_counts=0,
                                 levels=3,
                                 ax=ax,
                                 histogram2d_args=dict(bins=20,
                                        range=((limits[0], limits[1]),
                                               (limits[2], limits[3]))),
                                 plot_args=dict(marker='o',
                                                linestyle='none',
                                                color='black',
                                                alpha=0.7,
                                                markersize=2),
                                 contour_args=dict(
                                                   cmap=plt.cm.gray_r,
                                                   #cmap=cmap,
                                                   ),
                                 )



        if 1:
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
                                              h_sd_extent=(0.001, 1000),
                                              return_fractions=False,
                                              return_hisd=True)
                elif 'krumholz' in model:
                    model_fits = calc_krumholz(params,
                                              h_sd_extent=(0.001, 1000),
                                              return_fractions=False,
                                              return_hisd=True)

                ax.plot(model_fits[1], model_fits[2],
                        linestyle='-',
                        color=c_cycle[2],
                        alpha=alpha,
                        )
        else:
            if 0:
                l2 = ax.plot(h_sd_fit, hi_sd_fit_sternberg,
                        label='S+14',
                        color=c_cycle[1],
                        alpha=0.75,
                        )

            l3 = ax.plot(h_sd_fit, hi_sd_fit_krumholz,
                    linestyle='--',
                    label='K+09',
                    color=c_cycle[2],
                    alpha=0.75
                    )

        if i == 0:
            ax.legend(loc='lower right')

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
                    xytext=(0.95, 0.95),
                    xy=(0.95, 0.95),
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

        ax.set_xscale(scale[0])
        ax.set_yscale(scale[1])

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'$\Sigma_{\rm H\,I}$ + $\Sigma_{\rm H_2}$ ' + \
                       '[M$_\odot$ pc$^{-2}$]',)
        ax.set_ylabel(r'$\Sigma_{\rm H\,I}$ [M$_\odot$ pc$^{-2}$]',)

        if 'log' not in scale:
            ax.locator_params(nbins=6)

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
        plt.savefig(filename, bbox_inches='tight')

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
                plt.savefig('/usr/users/ezbc/scratch/core_' + core + '.png')

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
    print_av_error_stats(av_list[0], av_error_list[0])

    filename = results_dir + 'tables/multicloud_model_params.tex'
    write_model_params_table(model_analysis_dict,
                             filename,
                             models=('krumholz',))

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
                  #limits = [-4,3,1,10000],
                  names_list=cloud_name_list,
                  limits = [0.07,14,7,6000],
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
                filename = plot_kwargs['figure_dir'] + \
                           'models/' + cloud + '_hisd_vs_hsd.' + filetype
                plot_hi_vs_h_grid(hsd_cores_list[i],
                                  hisd_cores_list[i],
                                  core_names=core_names,
                                  model_results=model_results_list[i],
                                  model_analysis=\
                                      model_analysis_list[i]['cores'],
                                  limits=[-9, 159, -1.5, 15],
                                  levels=levels,
                                  scale=('linear', 'linear'),
                                  filename=filename,
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
def print_dict_keys(d):

    for key in d:
        print key
        if type(d[key]) is dict:
            print '--'
            print_dict_keys(d[key])

def write_model_params_table(mc_analysis_dict, filename, models=('krumholz',)):

    # Open file to be appended
    f = open(filename, 'wb')

    text_param_format ='{0:.1f}$^{{+{1:.1f}}}_{{-{2:.1f}}}$'

    #print_dict_keys(mc_analysis_dict)
    params_to_write = ['phi_cnm', 'alphaG']

    # Collect parameter names for each model for each core
    for cloud in ('california', 'perseus', 'taurus'):
        mc_analysis = mc_analysis_dict[cloud]
        for model in models:
            # each core correspond to a row
            for cloud_row, core_name in enumerate(mc_analysis['cores']):
                core = mc_analysis['cores'][core_name]
                if cloud_row == 0:
                    row_text = cloud.capitalize() + ' & '
                else:
                    row_text = ' & '
                row_text = add_row_element(row_text,
                                           core_name)
                if model == 'krumholz':
                    params_to_write = ['phi_cnm',]
                else:
                    params_to_write = ['alphaG',]
                for i, param_name in enumerate(params_to_write):
                    param = \
                        core[model + '_results'][param_name]
                    param_error = \
                        core[model + '_results'][param_name + '_error']

                    param_info = (param, param_error[0], param_error[1])

                    row_text = \
                        add_row_element(row_text,
                                        param_info,
                                        text_format=text_param_format)

                row_text += ' \\\\ \n'

                if cloud_row == len(mc_analysis['cores']) - 1:
                    row_text += '\hrule\n'

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

def get_cores_to_plot():

    '''

    '''

    # Which cores to include in analysis?
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
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/core_properties/'

    # define core properties
    with open(core_dir + cloud_name + '_core_properties.txt', 'r') as f:
        cores = json.load(f)

    header = data_dict['av_header']
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

    # load the bounding regions
    cores = load_ds9_core_region(cores,
                            filename=region_dir + \
                                     cloud_name + '_av_poly_cores.reg',
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

            if 0:
                print np.array(cores[core]['indices_orig']).shape
                print mask.shape
                print data.shape

                plt.clf(); plt.close()
                rh2_copy = data.copy()
                rh2_copy[cores[core]['indices_orig']] = 1000
                plt.imshow(rh2_copy, origin='lower')
                plt.title(core)
                plt.savefig('/usr/users/ezbc/scratch/core_' + core + '.png')

            cores[core]['mask'] = mask
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

def fit_steady_state_models(h_sd, rh2, model_kwargs):

    # Fit R_H2
    #---------
    sternberg_params = model_kwargs['sternberg_params']
    sternberg_results = {}
    krumholz_params = model_kwargs['krumholz_params']
    krumholz_results = {}

    # Fit to sternberg model
    if rh2.size > 3:
        alphaG, Z, phi_g = \
            fit_sternberg(h_sd,
                          rh2,
                          guesses=sternberg_params['guesses'],
                          vary=sternberg_params['param_vary'],
                          )

        # Fit to krumholz model
        phi_cnm, Z, phi_mol = \
            fit_krumholz(h_sd,
                         rh2,
                         guesses=krumholz_params['guesses'],
                         vary=krumholz_params['param_vary'],
                         )
    else:
        alphaG, Z, phi_g, phi_cnm, Z, phi_mol = 6 * [np.nan]

    # keep results
    sternberg_results['alphaG'] = alphaG
    sternberg_results['Z'] = Z
    sternberg_results['phi_g'] = phi_g


    # keep results
    krumholz_results['phi_cnm'] = phi_cnm
    krumholz_results['Z'] = Z
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

def calc_krumholz(params, h_sd_extent=(0.001, 500), return_fractions=True,
        return_hisd=False):

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

    # Create large array of h_sd
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

def fit_krumholz(h_sd, rh2, guesses=[10.0, 1.0, 10.0], rh2_error=None,
        verbose=False, vary=[True, True, True]):

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

    '''

    from scipy.optimize import curve_fit
    from scipy import stats
    from lmfit import minimize, Parameters, report_fit
    from myscience import krumholz09 as k09

    def chisq(params, h_sd, rh2):
        phi_cnm = params['phi_cnm'].value
        phi_mol = params['phi_mol'].value
        Z = params['Z'].value

        rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol)

        chisq = np.sum(np.abs(rh2 - rh2_model))

        return chisq

    def residual(params, h_sd, rh2):
        phi_cnm = params['phi_cnm'].value
        phi_mol = params['phi_mol'].value
        Z = params['Z'].value

        rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol)

        residual = rh2 - rh2_model

        return residual

    # Set parameter limits and initial guesses
    params = Parameters()
    params.add('phi_cnm',
               value=guesses[0],
               min=0.001,
               max=10000,
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
    result = minimize(residual,
                      params,
                      args=(h_sd, rh2),
                      method='leastsq')

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

def calc_sternberg(params, h_sd_extent=(0.001, 500),
        return_fractions=True, return_hisd=False):

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
        verbose=False, vary=[True, True, True]):

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

    '''

    from scipy.optimize import curve_fit
    from scipy import stats
    from lmfit import minimize, Parameters, report_fit
    from myscience import sternberg14 as s14

    def chisq(params, h_sd, rh2):
        alphaG = params['alphaG'].value
        phi_g = params['phi_g'].value
        Z = params['Z'].value

        rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                 return_fractions=False)

        chisq = np.sum(np.abs(rh2 - rh2_model))

        return chisq

    def residual(params, h_sd, rh2):
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
               min=0.001,
               max=1000,
               vary=vary[0])
    params.add('phi_g',
               value=guesses[2],
               min=0.5,
               max=2,
               vary=vary[2])
    params.add('Z',
               value=guesses[1],
               min=0.1,
               max=4,
               vary=vary[1])

    # Perform the fit!
    result = minimize(residual,
                      params,
                      args=(h_sd, rh2),
                      method='leastsq')

    rh2_fit_params = (params['alphaG'].value, params['Z'].value,
            params['phi_g'].value)

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

    filename_extension = global_args['cloud_name'] + '_' + global_args['data_type'] + \
            background_name + \
            bin_name + weights_name + \
            region_name + width_name + avthres_name + \
            intercept_name + error_name + compsub_name + backdgr_name + \
            '_' + hi_range_name

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

def simulate_nhi(hi_data, vel_axis, vel_range, vel_range_error):

    vel_range_sim = [vel_range[0] + np.random.normal(0, scale=vel_range_error),
                     vel_range[1] + np.random.normal(0, scale=vel_range_error)]

    nhi_sim = myia.calculate_nhi(cube=hi_data,
                              velocity_axis=vel_axis,
                              velocity_range=vel_range_sim,
                              )

    return nhi_sim.ravel(), vel_range_sim

def bootstrap_worker(global_args, i):

    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    av = global_args['av']
    av_error = global_args['av_error']
    nhi = global_args['nhi']
    nhi_back = global_args['nhi_back']
    hi_data = global_args['hi_data']
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
    av_sim, av_background_sim = \
            simulate_background_error(av_sim,
                                      scale=intercept_error)

    # calculate N(HI)
    if global_args['sim_hi_error']:
        nhi_sim, vel_range_sim = simulate_nhi(hi_data,
                                              vel_axis,
                                              vel_range,
                                              vel_range_error)
    else:
        nhi_sim = nhi

    # Bootstrap data
    # -------------------------------------------------------------------------
    boot_indices = np.random.choice(av.size, size=av.size,)# p=probabilities)

    av_boot = av_sim[boot_indices]
    av_error_boot = av_error[boot_indices]
    nhi_boot = nhi_sim[boot_indices]

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
                            nhi_background=nhi_back_boot,
                            init_guesses=init_guesses,
                            plot_kwargs=plot_kwargs,
                            use_intercept=use_intercept)

    # Calculate N(H2), then HI + H2 surf dens, fit relationship
    # -------------------------------------------------------------------------
    # calculate N(H2) maps
    nh2_image = calculate_nh2(nhi_image=nhi_sim,
                              av_image=av_sim,
                              dgr=av_model_results['dgr_cloud'])

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi_sim,
                               sd_factor=1/1.25)

    h2_sd_image = calculate_sd(nh2_image,
                               sd_factor=1/0.625)

    h_sd_image = hi_sd_image + h2_sd_image

    # Write ratio between H2 and HI
    rh2_image = h2_sd_image / hi_sd_image


    model_kwargs = global_args['ss_model_kwargs']['model_kwargs']
    cores_to_plot = global_args['ss_model_kwargs']['cores_to_plot']
    cores = global_args['ss_model_kwargs']['cores']

    ss_model_results = {}

    # cycle through each core, bootstrapping the pixels
    for core in cores_to_plot:
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
        rh2_core = rh2_image[core_indices]

        # mask negative ratios
        if 1:
            mask = (rh2_core < 0) | (np.isnan(rh2_core))
            #mask = (rh2_core < 0)
            rh2_core = rh2_core[~mask]
            h_sd_core = h_sd_core[~mask]

        ss_model_result = \
            fit_steady_state_models(h_sd_core,
                                    rh2_core,
                                    model_kwargs=model_kwargs,
                                    )
        ss_model_results[core] = ss_model_result

        if 0:
            print core, 'phi_cnm', ss_model_result['krumholz_results']['phi_cnm']
        if 0:
            import matplotlib.pyplot as plt
            plt.close(); plt.clf();
            rh2_fit, hsd_fit = \
                calc_sternberg(ss_model_result['sternberg_results'],
                              h_sd_extent=[1, 100],
                              return_fractions=False)
            rh2_fit, hsd_fit = \
                calc_krumholz(ss_model_result['krumholz_results'],
                              h_sd_extent=[1, 100],
                              return_fractions=False)
            #print ss_model_result['sternberg_results']['alphaG']

            plt.scatter(h_sd_core, rh2_core, color='k', alpha=0.1)
            plt.plot(hsd_fit, rh2_fit, color='r', alpha=1)
            plt.yscale('log')
            plt.xlim([0,80]); plt.ylim([0.001, 10])
            plt.savefig('/usr/users/ezbc/scratch/rh2_vs_hsd' + \
                        core + '_' + str(i) + '.png')


        phi_cnm = ss_model_result['krumholz_results']['phi_cnm']
        if phi_cnm < 0:
            print 'core = ', core, ' phi_cnm =', phi_cnm

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

    try:
        output = bootstrap_worker(args, i)
    except KeyboardInterrupt:
        raise KeyboardInterruptError()

    return output

def bootstrap_fits(av_data, nhi_image=None, hi_data=None, vel_axis=None,
        vel_range=None, vel_range_error=1, av_error_data=None,
        av_reference=None, nhi_image_background=None, num_bootstraps=100,
        plot_kwargs=None, scale_kwargs=None, use_intercept=True,
        sim_hi_error=False, ss_model_kwargs=None, multiprocess=True):

    import multiprocessing as mp
    import sys

    if av_error_data is None:
        av_error_data = np.ones(av_data.size)

    # mask for nans, arrays will be 1D
    (av, av_error, nhi, nhi_back, ), mask = \
        mask_nans((av_data, av_error_data, nhi_image, nhi_image_background),
                   return_mask=True)
    hi_data = hi_data[:, ~mask]
    cores = ss_model_kwargs['cores']
    for core in cores:
        if 0:
            import matplotlib.pyplot as plt
            plt.close(); plt.clf()
            nh2_copy = np.copy(av_data)
            nh2_copy[cores[core]['mask']] = np.nan
            plt.imshow(nh2_copy, origin='lower', interpolation='nearest')
            plt.savefig('/usr/users/ezbc/scratch/core.png')
        if cores[core]['mask'] is not None:
            if 0:
                av_copy = np.copy(av_data)
                av_copy[cores[core]['mask']] = np.nan
                av_copy[mask] = np.nan
                cores[core]['test_pix'] = av_copy[~np.isnan(av_copy)]

            cores[core]['mask'] = cores[core]['mask'][~mask]
            cores[core]['indices'] = np.where(cores[core]['mask'] == 0)[0]
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

    #global_args['queue'] = queue
    if 0:
        args_list = []
        for i in xrange(num_bootstraps):
            global_args['i'] = i
            args_list.append(global_args.copy())

    # Prep multiprocessing
    queue = mp.Queue(10)
    pool = mp.Pool()
    processes = []

    # bootstrap
    if multiprocess:
        try:
            for i in xrange(num_bootstraps):
                processes.append(pool.apply_async(bootstrap_worker_wrapper,
                                                  args=(global_args,i,)))
        except KeyboardInterruptError:
            pool.terminate()
            pool.join()
        pool.close()
        pool.join()

        # Get the results
        if 0:
            for p in processes:
                result = p.get()
                #result = queue.get()
                boot_results[:, result[0]] = result[1]

        mc_results = collect_bootstrap_results(processes, ss_model_kwargs)

        for i in xrange(len(processes)):
            #result = queue.get()
            result = processes[i].get()
            boot_results[:, i] = result['av_model_results'].values()

    else:
        for i in xrange(num_bootstraps):
            processes.append(bootstrap_worker(global_args, i))

        mc_results = collect_bootstrap_results(processes, ss_model_kwargs)

        for i in xrange(len(processes)):
            #result = queue.get()
            result = processes[i].get()
            boot_results[:, i] = result['av_model_results'].values()

    #for process in processes:
    #    process.set()

    # Start the processes
    #for process in processes:
    # Join so that we wait until all processes are finished
    #for process in processes:
    #    process.join()
    #    process.terminate()

    #if 1:

    return boot_results, mc_results

def collect_bootstrap_results(processes, ss_model_kwargs):

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
                 }
        core_dict['sternberg_results'] = \
                {'alphaG': empty(),
                 'Z': empty(),
                 'phi_g': empty(),
                 }
    mc_results['data_params'] = \
            {'av_background_sim': empty(),
             'vel_range_sim': np.empty((len(processes), 2)),
             'av_scalar_sim': empty()}

    for i in xrange(len(processes)):
        #result = queue.get()
        result = processes[i].get()

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
        nan_mask = (np.isnan(av_reference) | \
                    np.isnan(av_data) | \
                    np.isnan(av_error_data))
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

    return kwargs

def calc_hi_vel_range(hi_spectrum, hi_vel_axis, gauss_fit_kwargs,
        width_scale=2, co_spectrum=None, co_vel_axis=None, ncomps=1):

    from scipy.stats import nanmedian
    from myfitting import fit_gaussians

    hi_fits = fit_gaussians(hi_vel_axis,
            hi_spectrum, **gauss_fit_kwargs)

    # use either the gaussian closest to the CO peak, or the tallest gaussian
    if co_spectrum is not None:
        co_peak_vel = co_vel_axis[co_spectrum == np.nanmax(co_spectrum)][0]
        vel_diffs = []
        for i, param in enumerate(hi_fits[2]):
            vel_center = param[1]
            width = param[2]
            vel_diffs.append(np.abs(param[1] - co_peak_vel))

        # get the velocity range
        vel_diffs = np.asarray(vel_diffs)
        sort_indices = np.argsort(vel_diffs)
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
        hi_width_error = np.max(vel_diffs[cloud_comp_num])
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
    elif global_args['cloud_name'] == 'taurus':
        guesses = (35, 3, 5,
                   5, -15, 20,
                   #3, -2, 2,
                   )
        ncomps = 2
        ncomps_in_cloud = 1
    elif global_args['cloud_name'] == 'california':
        guesses = (50, 3, 5,
                   10, -3, 3,
                   12, -10, 10,
                   3, -30, 10,
                   #2, -20, 20,
                   )
        ncomps = 4
        ncomps_in_cloud = 2
    gauss_fit_kwargs = {'guesses': guesses,
                        'ncomps': ncomps,
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

def get_model_fit_kwargs(cloud_name):

    '''

    '''
    vary_alphaG = True # Vary alphaG in S+14 fit?
    vary_Z = False # Vary metallicity in S+14 fit?
    vary_phi_g = False # Vary phi_g in S+14 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    guesses=(1.0, 1.0, 1.0) # Guesses for (alphaG, Z, phi_g)
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
    guesses=(8.0, 1.0, 10.0) # Guesses for (phi_cnm, Z, phi_mol)
    h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

    krumholz_params = {}
    krumholz_params['param_vary'] = [vary_phi_cnm, vary_Z, vary_phi_mol]
    krumholz_params['error_method'] = error_method
    krumholz_params['alpha'] = alpha
    krumholz_params['guesses'] = guesses
    krumholz_params['h_sd_fit_range'] = h_sd_fit_range
    krumholz_params['results_filename'] = results_filename
    krumholz_params['parameters'] = ['phi_cnm', 'Z', 'phi_mol']

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

def calc_mc_analysis(mc_results):

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

                core_analysis[core][model][param] = param_value
                core_analysis[core][model][param + '_error'] = param_error

                params[param] = param_value


                if 0:
                    if param == 'phi_cnm':
                        import matplotlib.pyplot as plt
                        plt.close(); plt.clf()
                        plt.hist(param_results, bins=50)
                        plt.savefig('/usr/users/ezbc/scratch/'+core+\
                                    '_hist.png')

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
    results_dict['mc_analysis'] = calc_mc_analysis(results_dict['mc_results'])

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
    hi_dir = '/d/bip3/ezbc/' + cloud_name + '/data/hi/'
    co_dir = '/d/bip3/ezbc/' + cloud_name + '/data/co/'
    core_dir = \
       '/d/bip3/ezbc/' + cloud_name + '/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/' + cloud_name + '/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
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
    nhi_image = calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=velocity_range,
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
                     'nhi_image_background': nhi_image_background,
                     'region_mask': region_mask,
                     'scale_kwargs': scale_kwargs,
                     'hi_spectrum': hi_spectrum,
                     'hi_std_spectrum': hi_std_spectrum,
                     'co_spectrum': co_spectrum,
                     'hi_range_kwargs': hi_range_kwargs,
                     }

    # Get model fitting params
    model_fitting = get_model_fit_kwargs(cloud_name)

    # Get cores params
    cores = get_core_properties(data, cloud_name)

    # get the cores in the cloud
    cores_to_plot = get_cores_to_plot()
    cores_to_plot = trim_cores_to_plot(cores, cores_to_plot)

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

    bootstrap_filename = results_dir + filename_base + '_bootresults.npy'
    results_filename = results_dir + \
               'bootstrap_results/' + filename_base + \
               '_bootstrap_results.pickle'

    print('\n\tBeginning bootstrap monte carlo...')

    # Perform bootsrapping
    boot_result, mc_results = \
        bootstrap_fits(av_data_backsub,
                       nhi_image=nhi_image,
                       av_error_data=av_error_data,
                       nhi_image_background=nhi_image_background,
                       plot_kwargs=plot_kwargs,
                       hi_data=hi_data_crop,
                       vel_axis=hi_vel_axis_crop,
                       vel_range=velocity_range,
                       vel_range_error=2,
                       av_reference=av_data_ref,
                       use_intercept=global_args['use_intercept'],
                       num_bootstraps=global_args['num_bootstraps'],
                       scale_kwargs=scale_kwargs,
                       sim_hi_error=global_args['sim_hi_error'],
                       ss_model_kwargs=global_args['ss_model_kwargs'],
                       multiprocess=global_args['multiprocess'],
                       )
    np.save(bootstrap_filename, boot_result)

    results_dict = {'boot_result': boot_result,
                    'data': data,
                    'data_products': data_products,
                    'global_args': global_args,
                    'plot_kwargs': plot_kwargs,
                    'filenames': filenames,
                    'mc_results': mc_results,
                    }

    print('\n\tSaving results...')
    save_results(results_dict, global_args['results_filename'])
    results_dict = load_results(global_args['results_filename'])

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

    elements = (clouds, data_types, recalculate_likelihoods, bin_image,
            init_vel_width, fixed_width, use_intercept, av_mask_threshold,
            regions, subtract_comps, use_background, hi_range_calc)

    permutations = list(itertools.product(*elements))

    print('Number of permutations to run: ' + str(len(permutations)))

    #for cloud in clouds:
    for permutation in permutations:
        global_args = {
                'cloud_name':permutation[0],
                'load': 1,
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
                'num_bootstraps': 10000,
                'hi_range_calc': permutation[11],
                'sim_hi_error': True,
                'multiprocess': True,
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



