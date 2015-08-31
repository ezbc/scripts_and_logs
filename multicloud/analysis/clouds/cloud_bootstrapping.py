#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import cloudpy
import matplotlib.pyplot as plt

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

def plot_av_vs_nhi_grid(av_grid, nhi_grid, fit_params=None, filename=None,
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
                '\nIntercept = {0:.2f}'.format(fit_params['intercept']) + \
                intercept_error_text
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

'''

'''
from multiprocessing.queues import Queue

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

def create_cloud_model(av, nhi_background, dgr_background,):

    if nhi_background is None:
        return av

    return av - dgr_background * nhi_background

def create_background_model(av, nhi_cloud, dgr_cloud):

    return av - dgr_cloud * nhi_cloud

def create_filename_base(args):

    # Name of diagnostic files
    if args['background_subtract']:
        background_name = '_backsub'
    else:
        background_name = ''

    if args['bin_image']:
        bin_name = '_binned'
        args['bin_procedure'] = 'all'
    else:
        bin_name = ''
        args['bin_procedure'] = 'none'
    if args['fixed_width'] is None:
        width_name = ''
        init_vel_width = args['init_vel_width']
        vel_center_gauss_fit_kwargs = None
    else:
        if args['fixed_width'] == 'gaussfit':
            if args['cloud_name'] == 'perseus':
                guesses = (28, 3, 5,
                           2, -20, 20)
                ncomps = 2
            elif args['cloud_name'] == 'taurus':
                guesses = (28, 3, 5,
                           5, -30, 20,
                           3, -15, 5,
                           )
                ncomps = 3
            elif args['cloud_name'] == 'california':
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
        init_vel_width = args['fixed_width']
    if args['use_weights']:
        weights_name = '_weights'
        weights_filename = av_dir + \
           args['cloud_name'] + '_binweights.fits'
    else:
        weights_name = ''
        weights_filename = None
    if args['region'] is None:
        region_name = ''
        args['region_name'] = args['cloud_name']
    else:
        region_name = '_region' + args['region']
        args['region_name'] = args['cloud_name'] + args['region']
    if args['av_mask_threshold'] is not None:
        avthres_name = '_avthres'
    else:
        avthres_name = ''
    if not args['use_intercept']:
        intercept_name = '_noint'
    else:
        intercept_name = ''
    if args['recalculate_likelihoods']:
        error_name = '_errorrecalc'
    else:
        error_name = ''
    if args['subtract_comps']:
        compsub_name = '_compsub'
    else:
        compsub_name = ''
    if args['use_background']:
        backdgr_name = '_backdgr'
    else:
        backdgr_name = ''
    if args['hi_range_calc'] == 'gaussian':
        hi_range_name = 'gaussrange'
    else:
        hi_range_name = 'stdrange'

    filename_extension = args['cloud_name'] + '_' + args['data_type'] + \
            background_name + \
            bin_name + weights_name + \
            region_name + width_name + avthres_name + \
            intercept_name + error_name + compsub_name + backdgr_name

    return filename_extension, args

def mask_nans(arrays, return_mask=False):

    """ Masks any positions where any array in the list has a NaN.

    Parameters
    ----------
    arrays : tuple
        Tuple of arrays with same dimension.

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

def fit_model(av, nhi, av_error=None, algebraic=False, nhi_background=None,
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
                   min=0.0,
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

        if return_fit:
            return (result, params), (dgr_cloud, dgr_background, intercept)
        else:
            return (dgr_cloud, dgr_background, intercept)

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

    av_rescale = av / np.random.uniform(low=1, high=scalar)

    return av_rescale

def simulate_background_error(av, scale=1.0):

    # bias from 2MASS image, see lombardi et al. (2009), end of section 2
    av_bias = 0.1
    scale = (scale**2 + (av_bias)**2)**0.5

    av_background_sim = av + np.random.normal(0, scale=scale)

    return av_background_sim

def bootstrap_worker(args, i):

    av = args['av']
    av_error = args['av_error']
    nhi = args['nhi']
    nhi_back = args['nhi_back']
    init_guesses = args['init_guesses']
    plot_kwargs = args['plot_kwargs']
    use_intercept = args['use_intercept']
    probabilities = args['probabilities']
    av_scalar = args['scale_kwargs']['av_scalar']
    intercept_error = args['scale_kwargs']['intercept_error']
    #i = args['i']

    #queue = args['queue']

    # Create simulated data
    # -------------------------------------------------------------------------
    # add random noise
    av_sim = av + simulate_noise(av, av_error)

    # rescale the data somewhere between Planck and 2MASS:
    # rescaled Planck = Planck / beta where beta is between 1.0 and 1.4
    av_sim = simulate_rescaling(av_sim, scalar=av_scalar)

    # remove background
    av_sim = simulate_background_error(av_sim, scale=intercept_error)

    # Bootstrap data
    # -------------------------------------------------------------------------
    boot_indices = np.random.choice(av.size, size=av.size, p=probabilities)

    av_boot = av_sim[boot_indices]
    av_error_boot = av_error[boot_indices]
    nhi_boot = nhi[boot_indices]

    if nhi_back is not None:
        nhi_back_boot = nhi_back[boot_indices]
    else:
        nhi_back_boot = None

    # for plotting
    plot_kwargs['bootstrap_num'] = i

    # Fit the bootstrapped data
    # -------------------------------------------------------------------------
    boot_result = fit_model(av_boot,
                            nhi_boot,
                            av_error=av_error_boot,
                            nhi_background=nhi_back_boot,
                            init_guesses=init_guesses,
                            plot_kwargs=plot_kwargs,
                            use_intercept=use_intercept)

    args['init_guesses'] = boot_result

    # Plot distribution and fit
    if plot_kwargs['plot_diagnostics']:
    #if 1:
        dgr_cloud, dgr_background, intercept = boot_result
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
        plot_av_vs_nhi_grid(av_images,
                       nhi_images,
                       av_error=av_error_boot,
                       fit_params=fit_params,
                       contour_plot=plot_kwargs['av_nhi_contour'],
                       limits=plot_kwargs['av_nhi_limits'],
                       filename=filename,
                       )
        if 0:
        #if nhi_back_boot is not None:
            filename = plot_kwargs['figure_dir'] + \
                       'diagnostics/av_nhi/' + \
                       plot_kwargs['filename_base'] + \
                       '_av_vs_nhibackground_bootstrap' + \
                       '{0:03d}.png'.format(plot_kwargs['bootstrap_num'])
            fit_params = {'dgr': dgr_background,
                          'intercept': intercept}
            av_plot = create_background_model(av_boot,
                                         nhi_boot,
                                         dgr_cloud,)
            plot_av_vs_nhi(av_plot,
                           nhi_back_boot,
                           av_error=av_error_boot,
                           fit_params=fit_params,
                           contour_plot=plot_kwargs['av_nhi_contour'],
                           filename=filename,
                           )
    result = [i, boot_result]


    #queue.put(result)
    return result

def bootstrap_worker_wrapper(args, i):

    try:
        output = bootstrap_worker(args, i)
    except KeyboardInterrupt:
        raise KeyboardInterruptError()

    return output

def bootstrap_fits(av_data, nhi_image, av_error_data=None, av_reference=None,
        nhi_image_background=None, num_bootstraps=100, plot_kwargs=None,
        scale_kwargs=None, use_intercept=True):

    import multiprocessing as mp
    import sys

    if av_error_data is None:
        av_error_data = np.ones(av_data.size)

    # mask for nans, arrays will be 1D
    av, av_error, nhi, nhi_back = \
            mask_nans((av_data, av_error_data, nhi_image, nhi_image_background))

    probabilities = 1.0 / av_error**2
    probabilities /= np.nansum(probabilities)

    # for plotting
    plot_kwargs['num_bootstraps'] = num_bootstraps

    # initialize array for storing output
    boot_results = np.empty((3, num_bootstraps))
    init_guesses = [0.05, 0.05, 0.0] # dgr_cloud, dgr_background, intercept


    # Prep arguments
    args = {}
    args['av'] = av
    args['av_error'] = av_error
    args['nhi'] = nhi
    args['nhi_back'] = nhi_back
    args['init_guesses'] = init_guesses
    args['plot_kwargs'] = plot_kwargs
    args['use_intercept'] = use_intercept
    args['probabilities'] = probabilities
    args['scale_kwargs'] = scale_kwargs

    #args['queue'] = queue
    args_list = []
    for i in xrange(num_bootstraps):
        args['i'] = i
        args_list.append(args.copy())

    # Prep multiprocessing
    queue = mp.Queue(10)
    pool = mp.Pool()
    processes = []

    # bootstrap
    if 1:
        try:
            for i in xrange(num_bootstraps):
                processes.append(pool.apply_async(bootstrap_worker_wrapper,
                                                  args=(args,i,)))
        except KeyboardInterruptError:
            pool.terminate()
            pool.join()
        pool.close()
        pool.join()

        # Get the results
        for p in processes:
            result = p.get()
            #result = queue.get()
            boot_results[:, result[0]] = result[1]
    else:
        for i in xrange(num_bootstraps):
            processes.append(bootstrap_worker(args, i))

        for result in processes:
            #result = queue.get()
            boot_results[:, result[0]] = result[1]

    #for process in processes:
    #    process.set()

    # Start the processes
    #for process in processes:
    # Join so that we wait until all processes are finished
    #for process in processes:
    #    process.join()
    #    process.terminate()

    #if 1:

    return boot_results

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

                fits[:, i] = sp.polyfit(x,
                                        y,
                                        deg=1)

        intercept_error = np.std(fits[1])

    print intercept_error

    kwargs = {}
    kwargs['av_scalar'] = av_scalar
    kwargs['intercept'] = intercept
    kwargs['intercept_error'] = intercept_error

    return kwargs

def calc_gaussian_hi_range(hi_data, hi_vel_axis, gauss_fit_kwargs, mask=None,
        width_scale=2, co_data=None, co_vel_axis=None):

    from scipy.stats import nanmedian
    from myfitting import fit_gaussians

    if mask is None:
        mask = np.zeros(hi_data.shape[1:], dtype=bool)

    hi_spectrum = nanmedian(hi_data[:, ~mask], axis=1)

    vel_center_fits = fit_gaussians(hi_vel_axis,
            hi_spectrum, **gauss_fit_kwargs)

    # use either the gaussian closest to the CO peak, or the tallest gaussian
    if co_data is not None:
        co_spectrum = nanmedian(co_data[:, ~mask], axis=1)
        co_peak_vel = co_vel_axis[co_spectrum == np.nanmax(co_spectrum)]
        vel_diff = np.Inf
        for i, param in enumerate(vel_center_fits[2]):
            print np.abs(param[1] - co_peak_vel)
            if np.abs(param[1] - co_peak_vel) < vel_diff:
                vel_center = param[1]
                width = param[2]
                vel_diff = np.abs(param[1] - co_peak_vel)
                print('width', width)
                print('center', param[1])
    else:
        amp_max = -np.Inf
        for i, param in enumerate(vel_center_fits[2]):
            if param[0] > amp_max:
                amp_max = param[0]
                vel_center = param[1]
                width = param[2] * 4
                cloud_comp_num = i

    velocity_range = [vel_center - width * width_scale,
                      vel_center + width * width_scale]
    if 1:
        import matplotlib.pyplot as plt
        plt.close(); plt.clf()
        plt.plot(co_vel_axis, co_spectrum * 100)
        plt.plot(hi_vel_axis, hi_spectrum)
        plt.axvline(co_peak_vel)
        plt.savefig('/usr/users/ezbc/Desktop/spectrum.png')

    return velocity_range

def get_gauss_fit_kwargs(args):
    if args['cloud_name'] == 'perseus':
        guesses = (28, 3, 5,
                   2, -20, 20)
        ncomps = 2
    elif args['cloud_name'] == 'taurus':
        guesses = (28, 3, 5,
                   5, -30, 20,
                   3, -15, 5,
                   )
        ncomps = 3
    elif args['cloud_name'] == 'california':
        guesses = (50, 3, 5,
                   20, -10, 10,
                   3, -45, 10,
                   #2, -20, 20,
                   )
        ncomps = 3
    gauss_fit_kwargs = {'guesses': guesses,
                                   'ncomps': ncomps,
                                   #'width_scale': 2,
                                   }

    return gauss_fit_kwargs

def run_cloud_analysis(args,):

    from astropy.io import fits
    from myimage_analysis import calculate_nhi, calc_region_mask
    from mycoords import make_velocity_axis
    from mystats import calc_symmetric_error, calc_logL
    import myio
    import mystats

    if 1:
        cloud_name = args['cloud_name']
        region = args['region']
        load = args['load']
        data_type = args['data_type']
        background_subtract = args['background_subtract']

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
    filename_base, args = create_filename_base(args)

    # Load data
    if args['bin_image']:
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
                                   region_name=args['region_name'])

    av_data[region_mask] = np.nan
    av_data_ref[region_mask] = np.nan


    # Scale the data to the 2MASS K+09 data
    scale_kwargs = scale_av_with_refav(av_data, av_data_ref, av_error_data)
    av_data -= scale_kwargs['intercept']

    print('Subtracting background of ' + str(scale_kwargs['intercept']) + \
          ' from ' + cloud_name)

    #if debugging:
    if 1:
        import matplotlib.pyplot as plt
        plt.close(); plt.clf()
        plt.imshow(av_data, origin='lower')
        plt.savefig('/usr/users/ezbc/Desktop/avmap.png')

    hi_data, hi_header = fits.getdata(hi_filename, header=True)
    co_data, co_header = fits.getdata(co_filename, header=True)

    hi_data[:, region_mask] = np.nan
    co_data[:, region_mask] = np.nan

    hi_vel_axis = make_velocity_axis(hi_header)
    co_vel_axis = make_velocity_axis(co_header)



    # derive HI range
    gauss_fit_kwargs = get_gauss_fit_kwargs(args)
    if args['hi_range_calc'] == 'gaussian':
        velocity_range = \
                calc_gaussian_hi_range(hi_data,
                                      hi_vel_axis,
                                      gauss_fit_kwargs,
                                      mask=region_mask,
                                      co_data=co_data,
                                      co_vel_axis=co_vel_axis)
    else:
        velocity_range = [-5, 15]

    print('\nVelocity range = ', velocity_range)

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
    nhi_image[nhi_image < 0] = np.nan
    nhi_image_background[nhi_image_background < 0] = np.nan

    if not args['use_background']:
        nhi_image_background = None

    # set up plotting variables
    plot_kwargs = {
                   'figure_dir': figure_dir,
                   'cloud_name': cloud_name,
                   'filename_base': filename_base,
                   'plot_diagnostics': args['plot_diagnostics'],
                   #'av_nhi_contour': av_nhi_contour,
                   'av_nhi_contour': True,
                   'av_nhi_limits': [0, 20, -1, 9],
                   #'av_nhi_limits': None,
                    }

    print('Filename base = \n' + filename_base + '\n')

    bootstrap_filename = results_dir + filename_base + '_bootresults.npy'
    run_analysis = True
    if args['load']:
        exists = myio.check_file(bootstrap_filename)
        if exists:
            boot_result = np.load(bootstrap_filename)
            run_analysis = False
    if run_analysis:
        # Perform bootsrapping
        boot_result = bootstrap_fits(av_data,
                                     nhi_image,
                                     av_error_data=av_error_data,
                                     nhi_image_background=nhi_image_background,
                                     plot_kwargs=plot_kwargs,
                                     av_reference=av_data_ref,
                                     use_intercept=args['use_intercept'],
                                     num_bootstraps=args['num_bootstraps'],
                                     scale_kwargs=scale_kwargs,
                                     )
        np.save(results_dir + filename_base + '_bootresults.npy', boot_result)

    if 0:
        # fit model to data, get confidence intervals of fit
        mask = (np.isnan(av_data) | np.isnan(nhi_image))
        p, V = np.polyfit(av_data[~mask],
                          nhi_image[~mask], cov=True, deg=1)

        print(V[0,0], V[1,1])

        from scipy.optimize import curve_fit

        def func(x, dgr, intercept):
            return dgr * x + intercept

        popt, pcov = curve_fit(func, av_data[~mask], nhi_image[~mask])


    if 1:
        dgr_cloud, dgr_background, intercept = np.mean(boot_result, axis=1)
        dgr_cloud_error, dgr_background_error, intercept_error = \
                    np.std(boot_result, axis=1)

        filename = plot_kwargs['figure_dir'] + \
                   'av_nhi/' + plot_kwargs['filename_base'] + \
                   '_av_vs_nhi.png'
        av_cloud = create_cloud_model(av_data,
                                     nhi_image_background,
                                     dgr_background,)

        # Calculate conf interval
        from scipy import stats
        dgrs = boot_result[0]
        N = dgrs.size
        mean, sigma = np.nanmean(dgrs), stats.nanstd(dgrs)
        conf_int_a = stats.norm.interval(0.68, loc=mean, scale=sigma)
        dgr_cloud_error = (mean - conf_int_a[0], conf_int_a[1] - mean)
        dgr_cloud, dgr_cloud_error = mystats.calc_cdf_error(dgrs,
                                                            alpha=0.32)

        if args['use_background']:
            dgrs = boot_result[1]
            N = dgrs.size
            mean, sigma = np.nanmean(dgrs), stats.nanstd(dgrs)
            conf_int_a = stats.norm.interval(0.68, loc=mean, scale=sigma)
            dgr_background_error = (mean - conf_int_a[0], conf_int_a[1] - mean)
        else:
            dgr_background_error = (0.0, 0.0)

        if args['use_intercept']:
            intercepts = boot_result[2]
            N = intercepts.size
            mean, sigma = np.nanmean(intercepts), stats.nanstd(intercepts)
            conf_int_a = stats.norm.interval(0.68, loc=mean, scale=sigma)
            intercept_error = (mean - conf_int_a[0], conf_int_a[1] - mean)
            intercept, intercept_error = mystats.calc_cdf_error(intercepts,
                                                                alpha=0.32)
        else:
            intercept_error = (0.0, 0.0)

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

        fit_params = {
                      'dgr_cloud': dgr_cloud,
                      'dgr_cloud_error': dgr_cloud_error,#, dgr_cloud_error),
                      'dgr_background': dgr_background,
                      'dgr_background_error': dgr_background_error,
                      'intercept': intercept,
                      'intercept_error': intercept_error,
                      }


        print('plotting')
        plot_av_vs_nhi_grid(av_images,
                       nhi_images,
                       av_error=av_error_data,
                       fit_params=fit_params,
                       contour_plot=plot_kwargs['av_nhi_contour'],
                       limits=plot_kwargs['av_nhi_limits'],
                       filename=filename,
                       )

        # Plot distribution
        if args['use_intercept'] or args['use_background']:
            filename = plot_kwargs['figure_dir'] + \
                       'bootstrap_dists/' + plot_kwargs['filename_base'] + \
                       '_backdgr_vs_clouddgr.png'
            print('\nSaving bootstrap distributions to:\n' + filename)
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
            print('\nSaving bootstrap distributions to:\n' + filename)
            plot_bootstrap_hist(boot_result[0],
                            filename=filename,
                            axis_label=r'Cloud DGR [10$^{-20}$ cm$^2$ mag]',
                            statistics=(dgr_cloud, dgr_cloud_error)
                            )

    return boot_result

def main():

    import itertools

    results = {}

    clouds = (
              'perseus',
              'taurus',
              'california',
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
                     'std',
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
        args = {'cloud_name':permutation[0],
                'load': 0,
                'load_props': 0,
                #'data_type': 'planck',
                #'data_type': 'k09',
                #'data_type': 'planck_lee12mask',
                #'data_type': 'lee12',
                'data_type' : permutation[1],
                'background_subtract': 0,
                'recalculate_likelihoods': permutation[2],
                'bin_image': permutation[3],
                'use_weights': 0,
                'init_vel_width': permutation[4],
                #'fixed_width': 20,
                'fixed_width': permutation[5],
                'use_intercept': permutation[6],
                'av_mask_threshold': permutation[7],
                #'av_mask_threshold': 1.2,
                'binned_data_filename_ext': '_bin',
                #'likelihood_resolution': 'fine',
                'likelihood_resolution': 'coarse',
                'region': permutation[8],
                'subtract_comps': permutation[9],
                'plot_diagnostics': 1,
                'use_background': permutation[10],
                'num_bootstraps': 100,
                'hi_range_calc': permutation[11],
                }
        run_analysis = False
        if args['data_type'] in ('planck_lee12mask', 'lee12'):
            if args['cloud_name'] == 'perseus':
                run_analysis = True
        else:
            if args['cloud_name'] == 'california':
                if args['region'] is None:
                    run_analysis = True
            else:
                run_analysis = True

        if run_analysis:
            results[args['cloud_name']] = run_cloud_analysis(args)

if __name__ == '__main__':
    main()



