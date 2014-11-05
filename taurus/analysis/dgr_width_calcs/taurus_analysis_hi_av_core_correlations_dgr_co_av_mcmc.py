#!/usr/bin/python

''' Calculates the N(HI) / Av likelihoodelation for the taurus molecular cloud.
'''

import pyfits as pf
import numpy as np

class KeyboardInterruptError(Exception): pass

''' Plotting Functions
'''
def plot_likelihoods(likelihoods,velocity_centers,vel_widths,
        filename=None,show=True, returnimage=False):
    ''' Plots a heat map of likelihoodelation values as a function of velocity width
    and velocity center.
    '''

    # Import external modules
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 8
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale * 3 / 4.0,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              #'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)


    fig = plt.figure(figsize=(3,2))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="3%",
                 cbar_size='6%',
                 axes_pad=0,
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    # Unravel the likelihoods if raveled
    if len(likelihoods.shape) == 1:
        likelihoods = np.empty((velocity_centers.shape[0],
                                       vel_widths.shape[0]))
        likelihoods[:,:] = np.NaN
        count = 0
        try:
            for i, center in enumerate(velocity_centers):
                for j, width in enumerate(vel_widths):
                    likelihoods[i,j] = likelihoods[count]
                    count += 1
        except IndexError:
            print(' plot_likelihoods: O-d array input, cannot proceed')
    else:
           likelihoods = likelihoods

    image = np.ma.array(likelihoods, mask=np.isnan(likelihoods))

    ax = imagegrid[0]

    ax.set_xlabel('Velocity Width (km/s)')
    ax.set_ylabel('Velocity Center (km/s)')

    #ax.set_xticks(np.arange(0,vel_widths.shape[0],1)[::5],
    #        velocity_centers[::5])

    plt.rc('text', usetex=False)
    im = ax.imshow(image, interpolation='nearest', origin='lower',
            extent=[vel_widths[0],vel_widths[-1],
                    velocity_centers[0],velocity_centers[-1]],
            cmap=plt.cm.gist_stern,
            #cmap=plt.cm.gray,
            #norm=matplotlib.colors.LogNorm(),
            )
    cb = ax.cax.colorbar(im)

    #cb.set_clim(vmin=0.)
    # Write label to colorbar
    cb.set_label_text(r'log L')

    fractions = np.array([0.95, 0.68])
    levels = (1 + fractions * image.min())

    cs = ax.contour(image, levels=levels, origin='lower',
            extent=[vel_widths[0],vel_widths[-1],
                    velocity_centers[0],velocity_centers[-1]],
            colors='k'
            )

    # Define a class that forces representation of float to look a certain way
    # This remove trailing zero so '1.0' becomes '1'
    class nf(float):
         def __repr__(self):
             str = '%.1f' % (self.__float__(),)
             if str[-1]=='0':
                 return '%.0f' % self.__float__()
             else:
                 return '%.1f' % self.__float__()

    # Recast levels to new class
    cs.levels = [nf(val) for val in fractions*100.0]

    #fmt = {}
    #for level, fraction in zip(cs.levels, fractions):
    #    fmt[level] = fraction
    fmt = '%r %%'

    ax.clabel(cs, cs.levels, fmt=fmt, fontsize=9, inline=1)

    if filename is not None:
        plt.savefig(filename,bbox_inches='tight')
    if show:
        plt.draw()
        plt.show()
    if returnimage:
        return likelihoods

def plot_likelihoods_hist(likelihoods, x_grid, y_grid, y_pdf=None,
        x_pdf=None, x_confint=None, y_confint=None, filename=None, show=True,
        returnimage=False, plot_axes=('centers', 'widths'),
        contour_confs=None, npix=None, av_thres=None):

    ''' Plots a heat map of likelihoodelation values as a function of velocity width
    and velocity center.
    '''

    # Import external modules
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale * 3 / 4.0,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              #'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    fig, ax_image = plt.subplots(figsize=(6,6))

    # Mask NaNs
    image = np.ma.array(likelihoods, mask=np.isnan(likelihoods))

    if plot_axes[0] == 'centers':
        x_extent = x_grid[0], x_grid[-1]
        ax_image.set_xlabel(r'Velocity Center (km/s)')
        x_sum_axes = (1, 2)
        y_pdf_label = r'Centers PDF'
    if plot_axes[1] == 'centers':
        y_extent = y_grid[0], y_grid[-1]
        ax_image.set_ylabel(r'Velocity Center (km/s)')
        y_sum_axes = (1, 2)
        x_pdf_label = r'Centers PDF'
    if plot_axes[0] == 'widths':
        x_extent = x_grid[0], x_grid[-1]
        ax_image.set_xlabel(r'Velocity Width (km/s)')
        x_sum_axes = (0, 2)
        y_pdf_label = r'Width PDF'
        x_limits = (0, x_grid[-1])
    if plot_axes[1] == 'widths':
        y_extent = y_grid[0], y_grid[-1]
        ax_image.set_ylabel(r'Velocity Width (km/s)')
        y_sum_axes = (0, 2)
        x_pdf_label = r'Width PDF'
    if plot_axes[0] == 'dgrs':
        x_extent = x_grid[0], x_grid[-1]
        ax_image.set_xlabel(r'DGR (10$^{-20}$ cm$^2$ mag$^1$)')
        x_sum_axes = (0, 1)
        y_pdf_label = r'DGR PDF'
    if plot_axes[1] == 'dgrs':
        y_extent = y_grid[0], y_grid[-1]
        ax_image.set_ylabel(r'DGR (10$^{-20}$ cm$^2$ mag$^1$)')
        y_sum_axes = (0, 1)
        x_pdf_label = r'DGR PDF'
        y_limits = (0.0, y_grid[-1])

    sum_axes = np.array((x_sum_axes, y_sum_axes))
    sum_axis = np.argmax(np.bincount(np.ravel(sum_axes)))

    # Create likelihood image
    image = np.sum(likelihoods, axis=sum_axis) / np.sum(likelihoods)

    # Derive marginal distributions of both centers and widths

    x_sum = np.sum(likelihoods, axis=x_sum_axes)
    x_pdf = x_sum / np.sum(x_sum)
    y_sum = np.sum(likelihoods, axis=y_sum_axes)
    y_pdf = y_sum / np.sum(y_sum)

    extent = np.ravel(np.array((x_extent, y_extent)))

    #plt.rc('text', usetex=False)
    im = ax_image.imshow(image.T, interpolation='nearest', origin='lower',
            extent=extent,
            #cmap=plt.cm.gist_stern,
            #cmap=plt.cm.gray,
            cmap=plt.cm.binary,
            #norm=matplotlib.colors.LogNorm(),
            aspect='auto',
            )

    show_pdfs = 1

    if show_pdfs:
        divider = make_axes_locatable(ax_image)
        ax_pdf_x = divider.append_axes("top", 1, pad=0.1, sharex=ax_image)
        ax_pdf_y  = divider.append_axes("right", 1, pad=0.1,
                sharey=ax_image)

        # make some labels invisible
        plt.setp(ax_pdf_x.get_xticklabels() + \
                 ax_pdf_y.get_yticklabels(),
                 visible=False)

        ax_pdf_x.plot(x_grid,
                      x_pdf,
                      color='k',
                      drawstyle='steps-post',
                      linewidth=2,
                      )

        ax_pdf_y.plot(y_pdf,
                      y_grid,
                      color='k',
                      drawstyle='steps-post',
                      linewidth=2,
                      )

        #axHistx.axis["bottom"].major_ticklabels.set_visible(False)

        # Tick marks on the pdf?
        pdf_ticks = False

        for tl in ax_pdf_x.get_xticklabels():
            tl.set_visible(False)

        if pdf_ticks:
            wmax = x_pdf.max()
            ticks = [0, 0.5*wmax, 1.0*wmax]
            tick_labels = ['{0:.1f}'.format(ticks[0]),
                           '{0:.1f}'.format(ticks[1]),
                           '{0:.1f}'.format(ticks[2]),
                            ]
            ax_pdf_x.set_yticks(ticks)
            ax_pdf_x.set_yticklabels(tick_labels)
        else:
            for tl in ax_pdf_x.get_yticklabels():
                tl.set_visible(False)

        ax_pdf_x.set_ylabel(y_pdf_label)

        for tl in ax_pdf_y.get_yticklabels():
            tl.set_visible(False)
        if pdf_ticks:
            cmax = y_pdf.max()
            ticks = [0, 0.5*cmax, 1.0*cmax]
            tick_labels = ['{0:.1f}'.format(ticks[0]),
                           '{0:.1f}'.format(ticks[1]),
                           '{0:.1f}'.format(ticks[2]),
                            ]
            ax_pdf_y.set_xticks(ticks)
            ax_pdf_y.set_xticklabels(tick_labels)
        else:
            for tl in ax_pdf_y.get_xticklabels():
                tl.set_visible(False)

        ax_pdf_y.set_xlabel(x_pdf_label)

        # Show confidence limits
        if y_confint is not None:
            ax_pdf_y.axhspan(y_confint[0] - y_confint[1],
                             y_confint[0] + y_confint[2],
                             color='k',
                             linewidth=1,
                             alpha=0.2)
            ax_pdf_y.axhline(y_confint[0],
                             color='k',
                             linestyle='--',
                             linewidth=3,
                             alpha=1)
        if x_confint is not None:
            ax_pdf_x.axvspan(x_confint[0] - x_confint[1],
                                 x_confint[0] + x_confint[2],
                                  color='k',
                                 linewidth=1,
                                  alpha=0.2)
            ax_pdf_x.axvline(x_confint[0],
                                 color='k',
                                 linestyle='--',
                                 linewidth=3,
                                 alpha=1)

    #cb.set_clim(vmin=0.)
    # Write label to colorbar
    #cb.set_label_text(r'log L')

    # Plot contours
    if contour_confs is not None:

        fractions = (1.0 - np.asarray(contour_confs))
        levels = (fractions * image.max())

        cs = ax_image.contour(image.T, levels=levels, origin='lower',
                extent=extent,
                colors='k'
                )

        # Define a class that forces representation of float to look a certain
        # way This remove trailing zero so '1.0' becomes '1'
        class nf(float):
             def __repr__(self):
                 str = '%.1f' % (self.__float__(),)
                 if str[-1]=='0':
                     return '%.0f' % self.__float__()
                 else:
                     return '%.1f' % self.__float__()

        # Recast levels to new class
        cs.levels = [nf(val) for val in np.asarray(contour_confs)*100.0]

        #fmt = {}
        #for level, fraction in zip(cs.levels, fractions):
        #    fmt[level] = fraction
        fmt = '%r %%'

        ax_image.clabel(cs, cs.levels, fmt=fmt, fontsize=9, inline=1)

    try:
        ax_image.set_xlim(x_limits)
        ax_image.set_ylim(y_limits)
    except UnboundLocalError:
        pass

    if npix is not None or av_thres is not None:
        text = ''
        if npix is not None:
            text += r'N$_{\rm pix}$ = ' + \
                     '{0:.0f}'.format(npix)
            if av_thres is not None:
                text += '\n'
        if av_thres is not None:
            text += r'$A_V$ threshold = {0:.1f} mag'.format(av_thres)
            text += '\n'
        text += r'DGR = {0:.2f} '.format(y_confint[0]) + \
                r'$\times$ 10$^{-20}$ (cm$^2$ mag$^1$)'
        text += '\n'
        text += r'Velocity width = {0:.2f} '.format(x_confint[0]) + \
                r'km/s'
        ax_image.annotate(text,
                xytext=(0.95, 0.95),
                xy=(0.95, 0.95),
                textcoords='axes fraction',
                xycoords='axes fraction',
                color='k',
                fontsize=font_scale*0.75,
                bbox=dict(boxstyle='round',
                          facecolor='w',
                          alpha=0.3),
                horizontalalignment='right',
                verticalalignment='top',
                )

    if filename is not None:
        plt.draw()
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.draw()
        plt.show()
    if returnimage:
        return likelihoods

def plot_av_image(av_image=None, header=None, title=None,
        limits=None, savedir='./', filename=None, show=True):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 15
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (8, 7),
              'figure.titlesize': font_scale
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    nrows_ncols=(1,1)
    ngrids=1

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=1,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # create axes
    ax = imagegrid[0]
    cmap = cm.jet # colormap
    # show the image
    im = ax.imshow(av_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            #norm=matplotlib.colors.LogNorm()
            vmin=0,
            vmax=1.4
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',)
    ax.set_ylabel('Declination (J2000)',)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Write label to colorbar
    cb.set_label_text(r'A$_V$ (Mag)',)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        fig.show()

''' Calculations
'''

def calc_likelihood(return_likelihoods=True, plot_results=True,
        results_filename='', likelihood_filename=None, clobber=False,
        conf=0.68, contour_confs=None, npix=None):

    '''
    Parameters
    ----------

    '''

    import numpy as np
    from myimage_analysis import calculate_nhi
    from scipy import signal
    from os import path
    from astropy.io import fits
    import multiprocessing
    import itertools
    import emcee
    import triangle
    import matplotlib.pyplot as plt

    # Check if likelihood grid should be derived
    if likelihood_filename is not None:
        if not path.isfile(likelihood_filename):
            perform_mle = True
            write_mle = True
        elif clobber:
            perform_mle = True
            write_mle = True
        else:
            perform_mle = False
            write_mle = False
    # If no filename provided, do not read file and do not write file
    else:
        write_mle = False
        perform_mle = True

    # Perform MCMC estimation of likelihood space
    # Else if estimation already performed, then read in previous estimation
    if perform_mle:
        # initialize the walkers in a tiny Gaussian ball around the maximum
        # likelihood result
        walker_pos = [_init_guesses + _init_spread * np.random.randn(_ndim) \
                for i in range(_nwalkers)]

        sampler = emcee.EnsembleSampler(_nwalkers, _ndim, logL_prob,
                threads=_mc_threads)

        sampler.run_mcmc(walker_pos, _niter)

        # Run the MCMC saving the results incrementally
        '''
        f = open(_likelihood_dir + _progress_filename, "w")
        f.close()

        for result in sampler.sample(walker_pos,
        	                         iterations=_niter,
        	                         storechain=False):
            position = result[0]
            f = open(_likelihood_dir + _progress_filename, "a")
            for k in range(position.shape[0]):
                f.write("{0:4d\t {1:s}\n".format(k, " ".join(position[k])))
            f.close()
        '''

        # Flatten the chain
        samples = sampler.chain[:, 50:, :].reshape((-1, _ndim))

        # Save the samples
        np.save(likelihood_filename, samples)

        # plot
        fig = triangle.corner(samples, labels=["Width (km/s)",
            "DGR 10$^{20}$ cm$^2$ mag$^1$)", "$A_V$ (mag)"])
        plt.savefig(results_filename)

    elif not perform_mle:
        print('Reading likelihood grid file:')
        print(likelihood_filename)

        hdu = fits.open(likelihood_filename)
        likelihoods = hdu[0].data

        if len(velocity_centers) != likelihoods.shape[0] or \
            len(vel_widths) != likelihoods.shape[1]:
            raise ValueError('Specified parameter grid not the same as in' + \
                    'loaded data likelihoods.')

        likelihoods = np.ma.array(likelihoods,
                mask=(likelihoods != likelihoods))

    # Define parameter resolutions
    #delta_center = velocity_centers[1] - velocity_centers[0]
    #delta_width = vel_widths[1] - vel_widths[0]

    # Derive marginal distributions of both centers and widths
    center_likelihood = np.sum(likelihoods, axis=(1,2)) / \
            np.sum(likelihoods)
    width_likelihood = np.sum(likelihoods, axis=(0,2)) / \
            np.sum(likelihoods)
    dgr_likelihood = np.sum(likelihoods, axis=(0,1)) / \
            np.sum(likelihoods)

    # Derive confidence intervals of parameters
    center_confint = threshold_area(velocity_centers,
                                    center_likelihood,
                                    area_fraction=conf)
    width_confint = threshold_area(vel_widths,
                                   width_likelihood,
                                   area_fraction=conf)
    dgr_confint = threshold_area(dgrs,
                                 dgr_likelihood,
                                 area_fraction=conf)

    print('Velocity widths = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(width_confint[0],
                                                    width_confint[2],
                                                    np.abs(width_confint[1])))
    print('Velocity centers = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(center_confint[0],
                                                    center_confint[2],
                                                    np.abs(center_confint[1])))
    print('DGRs = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(dgr_confint[0],
                                                    dgr_confint[2],
                                                    np.abs(dgr_confint[1])))

    # Write PDF
    center = center_confint[0]
    upper_lim = (center_confint[0] + width_confint[0]/2.)
    lower_lim = (center_confint[0] - width_confint[0]/2.)
    upper_lim_error = (center_confint[2]**2 + width_confint[2]**2)**0.5
    lower_lim_error = (center_confint[1]**2 + width_confint[1]**2)**0.5

    vel_range_confint = (lower_lim, upper_lim, lower_lim_error,
            upper_lim_error)

    if plot_results:
        pass
    if not return_likelihoods:
        return vel_range_confint, dgr_confint
    else:
        return (vel_range_confint, dgr_confint, likelihoods,
            center_likelihood, width_likelihood, dgr_likelihood)

def calc_logL(theta):

    '''
    Calculates log likelihood

    http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node2.html

    '''

    from myimage_analysis import calculate_nhi

    # Unpack parameters
    vel_width, dgr, av_thres = theta


    # Define velocity range
    vel_range = (_velocity_center - vel_width / 2.,
                 _velocity_center + vel_width / 2.)

    # Derive N(HI) maps
    nhi_image_temp, nhi_image_error = \
            calculate_nhi(cube=_hi_cube,
                velocity_axis=_hi_velocity_axis,
                velocity_range=vel_range,
                noise_cube=_hi_noise_cube)

    # Avoid NaNs, and mask images with Av threshold
    indices = np.where((nhi_image_temp == nhi_image_temp) & \
                       (_av_image == _av_image) & \
                       (_av_image <= av_thres))

    nhi_image_sub = nhi_image_temp[indices]
    nhi_image_sub_error = nhi_image_error[indices]
    av_image_sub = _av_image[indices]
    av_image_sub_error = np.average(_av_image_error[indices])

    # Create model of Av with N(HI) and DGR
    av_image_model = nhi_image_sub * dgr

    data = av_image_sub
    model = av_image_model
    error = av_image_sub_error
    N = av_image_sub.size

    if error > 0 and N > 0:
        #logL = -(np.sum((av_image_sub - av_image_model)**2 / \
        #         (2 * av_image_sub_error**2)) - np.log(av_image_sub.size))
        logL = - np.sum((data - model)**2 / (2 * error**2)) \
               - 0.5 * N * np.log(np.sum(2 * np.pi * error**2))
    else:
    	logL = -np.inf

    print theta, np.median(av_image_model) - np.median(av_image_sub), logL

    if np.isnan(logL):
     	return -np.inf

    return logL

def logL_prior(theta):

    ''' Function to define range of parameter space.

    '''

    vel_width, dgr, av_thres = theta
    if _vel_width_range[0] < vel_width < _vel_width_range[1] and \
       _dgr_range[0] < dgr < _dgr_range[1] and \
       _av_thres_range[0] < av_thres < _av_thres_range[1]:
        return 0.0
    return -np.inf

def logL_prob(theta):

    ''' Log likelihood function.

    '''

    lp = logL_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + calc_logL(theta)

def threshold_area(x, y, area_fraction=0.68):

    '''
    Finds the limits of a 1D array which includes a given fraction of the
    integrated data.

    Parameters
    ----------
    data : array-like
        1D array.
    area_fraction : float
        Fraction of area.

    Returns
    -------
    limits : tuple
        Lower and upper bound including fraction of area.

    '''

    import numpy as np
    from scipy.integrate import simps as integrate

    # Check if size of data
    if x.size == 1:
        return x[0], 0, 0

    # Step for lowering threshold
    step = (np.max(y) - np.median(y)) / 10000.0

    # initial threshold
    threshold = np.max(y) - step
    threshold_area = 0.0

    # area under whole function
    area = integrate(y, x)

    # Stop when the area below the threshold is greater than the max area
    while threshold_area < area * area_fraction and threshold > 0:

        threshold_indices = np.where(y > threshold)[0]

        try:
            bounds_indices = (threshold_indices[0], threshold_indices[-1])
        except IndexError:
            bounds_indices = ()

        try:
            threshold_area = integrate(y[bounds_indices[0]:bounds_indices[1]],
                                       x[bounds_indices[0]:bounds_indices[1]])
            threshold_area += threshold * (x[bounds_indices[1]] - \
                                           x[bounds_indices[0]])
        except IndexError:
            threshold_area = 0

        threshold -= step

    if threshold < 0:
        bounds_indices = (0, len(x) - 1)

    x_peak = x[y == y.max()][0]
    low_error, up_error = x_peak - x[bounds_indices[0]], \
                          x[bounds_indices[1]] - x_peak

    return (x_peak, low_error, up_error)

def write_mle_tofits(filename='', velocity_centers=None,
        vel_widths=None, dgrs=None, likelihoods=None, clobber=False):

    from astropy.io import fits

    print('Writing likelihood grid to file:')
    print(filename)

    header = fits.Header()
    header['NAXIS'] = 3
    header['CTYPE1'] = 'CENTERS'
    header['CTYPE2'] = 'WIDTHS'
    header['CTYPE3'] = 'DGR'
    header['CRPIX1'] = 0
    header['CRPIX2'] = 0
    header['CRPIX3'] = 0
    header['CRVAL1'] = velocity_centers[0]
    header['CRVAL2'] = vel_widths[0]
    header['CRVAL3'] = dgrs[0]
    try:
        header['CDELT1'] = velocity_centers[1] - velocity_centers[0]
    except IndexError:
        header['CDELT1'] = 1
    try:
        header['CDELT2'] = vel_widths[1] - vel_widths[0]
    except IndexError:
        header['CDELT2'] = 1
    try:
        header['CDELT3'] = dgrs[1] - dgrs[0]
    except IndexError:
        header['CDELT3'] = 1

    fits.writeto(filename,
                 likelihoods,
                 header,
                 clobber=clobber)

''' DS9 Region and Coordinate Functions
'''
def convert_core_coordinates(cores, header):

    for core in cores:
        cores[core].update({'box_pixel': 0})
        cores[core].update({'center_pixel': 0})

        #box_wcs = cores[core]['box_wcs']
        #box_pixel = len(box_wcs) * [0,]
        center_wcs = cores[core]['center_wcs']

        # convert centers to pixel coords
        cores[core]['center_pixel'] = get_pix_coords(ra=center_wcs[0],
                                                     dec=center_wcs[1],
                                                     header=header)[:2].tolist()
        # convert box corners to pixel coords
        #for i in range(len(box_wcs)/2):
        #    pixels = get_pix_coords(ra=box_wcs[2*i], dec=box_wcs[2*i + 1],
        #            header=header)
        #    box_pixel[2*i], box_pixel[2*i + 1] = int(pixels[0]), int(pixels[1])
        #cores[core]['box_pixel'] = box_pixel

    return cores

def convert_limit_coordinates(prop_dict):

    prop_dict.update({'limit_pixels': []})

    header = prop_dict['av_header']

    limit_wcs = prop_dict['limit_wcs']

    for limits in limit_wcs:
        # convert centers to pixel coords
        limit_pixels = get_pix_coords(ra=limits[0],
                                     dec=limits[1],
                                     header=header)[:2].tolist()

        prop_dict['limit_pixels'].append(limit_pixels[0])
        prop_dict['limit_pixels'].append(limit_pixels[1])

    return prop_dict

def load_fits(filename,return_header=False):
    ''' Loads a fits file.
    '''

    import pyfits as pf

    f = pf.open(filename)
    if return_header:
        return f[0].data,f[0].header
    else:
        return f[0].data

def get_sub_image(image, indices):

    return image[indices[1]:indices[3],
            indices[0]:indices[2]]

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec), or Ra in degrees
    and dec in degrees.
    '''

    import pywcsgrid2 as wcs
    import pywcs

    # convert to degrees if ra and dec are array-like
    try:
        if len(ra) == 3 and len(dec) == 3:
            ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)
        else:
            raise ValueError('RA and Dec must be in (hrs,min,sec) and' + \
                    ' (deg,arcmin,arcsec) or in degrees.')
    except TypeError:
        ra_deg, dec_deg = ra, dec

    wcs_header = pywcs.WCS(header)
    pix_coords = wcs_header.wcs_sky2pix([[ra_deg, dec_deg, 0]], 0)[0]

    return pix_coords

def hrs2degs(ra=None, dec=None):
    ''' Ra and dec tuples in hrs min sec and deg arcmin arcsec.
    '''

    ra_deg = 15*(ra[0] + ra[1]/60. + ra[2]/3600.)
    dec_deg = dec[0] + dec[1]/60. + dec[2]/3600.

    return (ra_deg, dec_deg)

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

    return region[0].coord_list

def load_ds9_region(cores, filename_base = 'taurus_av_boxes_', header=None):

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    for core in cores:
        region = read_ds9_region(filename_base + core + '.reg')
        if region is not None:
            box_center_pixel = get_pix_coords(ra = region[0],
                                              dec = region[1],
                                              header = header)
            box_center_pixel = (int(box_center_pixel[1]),
                    int(box_center_pixel[0]))
            wedge_radius = region[2] / header['CDELT1']
            wedge_angle = region[3] / header['CDELT2']
            cores[core].update({'box_center_pix': box_center_pixel})
            cores[core].update({'wedge_angle': wedge_angle})
            cores[core].update({'wedge_radius': wedge_radius})
            cores[core].update({'box_angle': region[4]})

    return cores

'''
The main script
'''

def main():

    import grid
    import numpy as np
    import numpy
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    from multiprocessing import Pool

    global _hi_cube
    global _hi_velocity_axis
    global _hi_noise_cube
    global _av_image
    global _av_image_error

    # parameters used in script
    # -------------------------
    # HI velocity integration range
    # Determine HI integration velocity by CO or likelihoodelation with Av?
    hi_av_likelihoodelation = True

    center_vary = False
    width_vary = True
    dgr_vary = True

    # Check if likelihood file already written, rewrite?
    clobber = 1

    # Include only pixels within core regions for analysis?
    core_mask = 0

    # Confidence of parameter errors
    conf = 0.68
    # Confidence of contour levels
    contour_confs = (0.68, 0.95)

    # Results and fits filenames
    likelihood_filename = 'taurus_nhi_av_likelihoods_mcmc_co_av'
    results_filename = 'taurus_likelihood_mcmc_co_av'
    global _progress_filename
    _progress_filename = 'taurus_mcmc_samples.dat'

    # Define ranges of parameters
    global _av_thres_range
    _av_thres_range = (1.0, 1.1)
    _av_thres_range = (0.1, 2.0)
    global _vel_width_range
    _vel_width_range = (0.0, 80.0)
    global _dgr_range
    _dgr_range = (0.01, 0.4)
    global _velocity_center
    _velocity_center = 5.0 # km/s

    # MCMC parameters
    global _ndim
    _ndim = 3
    global _nwalkers
    _nwalkers = 30
    global _niter
    _niter = 100
    global _init_guesses
    _init_guesses = np.array((10, 0.10, 1.0))
    global _init_spread
    _init_spread = np.array((0.1, 0.01, 0.01))
    global _mc_threads
    _mc_threads = 10

    # Name of property files results are written to
    global_property_file = 'taurus_global_properties.txt'
    core_property_file = 'taurus_core_properties.txt'

    # Name of noise cube
    noise_cube_filename = 'taurus_hi_galfa_cube_regrid_planckres_noise.fits'

    # Define limits for plotting the map
    prop_dict = {}
    prop_dict['limit_wcs'] = (((3, 58, 0), (27, 6, 0)),
                              ((3, 20, 0), (35, 0, 0)))
    prop_dict['limit_wcs'] = (((3, 58, 0), (26, 6, 0)),
                              ((3, 0, 0), (35, 0, 0)))

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/hi_velocity_range/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/hi/'
    co_dir = '/d/bip3/ezbc/taurus/data/co/'
    core_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/taurus/data/python_output/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'
    likelihood_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    global _likelihood_dir
    _likelihood_dir = likelihood_dir

    # load Planck Av and GALFA HI images, on same grid
    av_data_planck, av_header = load_fits(av_dir + \
                'taurus_av_planck_5arcmin.fits',
            return_header=True)
    prop_dict['av_header'] = av_header

    av_error_data_planck, av_error_header = load_fits(av_dir + \
                'taurus_av_error_planck_5arcmin.fits',
            return_header=True)

    hi_data, h = load_fits(hi_dir + \
                'taurus_hi_galfa_cube_regrid_planckres.fits',
            return_header=True)

    co_data, co_header = load_fits(co_dir + \
                'taurus_co_cfa_cube_regrid_planckres.fits',
            return_header=True)

    # make the velocity axis
    velocity_axis = make_velocity_axis(h)

    # Plot NHI vs. Av for a given velocity range
    if not path.isfile(hi_dir + noise_cube_filename):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=velocity_axis,
                velocity_noise_range=[90,110], header=h, Tsys=30.,
                filename=hi_dir + noise_cube_filename)
    else:
        noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    # define core properties
    with open(core_dir + core_property_file, 'r') as f:
        cores = json.load(f)
    with open(property_dir + global_property_file, 'r') as f:
        global_props = json.load(f)

    cores = convert_core_coordinates(cores, h)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'taurus_av_boxes_',
            header = h)

    print('\nCalculating likelihoods globally')

    mask = np.zeros(av_data_planck.shape)
    for core in cores:
        # Grab the mask
        mask += myg.get_polygon_mask(av_data_planck,
                cores[core]['wedge_vertices_rotated'])

    co_mom0 = np.sum(co_data, axis=0)

    # Mask images
    core_mask = 0
    if core_mask:
        indices = ((mask == 1) & \
                   (co_mom0 < np.std(co_mom0[~np.isnan(co_mom0)])*2.0))
        mask_type = '_core_mask'
    else:
        indices = (co_mom0 < np.std(co_mom0[~np.isnan(co_mom0)])*2.0)
        mask_type = ''

    hi_data_sub = np.copy(hi_data[:, indices])
    noise_cube_sub = np.copy(noise_cube[:, indices])
    av_data_sub = np.copy(av_data_planck[indices])
    av_error_data_sub = np.copy(av_error_data_planck[indices])

    # Set global variables
    _hi_cube = hi_data_sub
    _hi_velocity_axis = velocity_axis
    _hi_noise_cube = noise_cube_sub
    _av_image = av_data_sub
    _av_image_error = av_error_data_sub

    # Define filename for plotting results
    results_filename = figure_dir + results_filename

    # likelihoodelate each core region Av and N(HI) for velocity ranges
    vel_range_confint, dgr_confint, likelihoods, center_likelihood,\
        width_likelihood, dgr_likelihood = \
            calc_likelihood(return_likelihoods=True,
                            plot_results=True,
                            results_filename=results_filename + mask_type,
                            likelihood_filename=likelihood_dir + \
                                    likelihood_filename + \
                                    mask_type + '.npy',
                            clobber=clobber,
                            conf=conf,
                            contour_confs=contour_confs)

    '''
    print('HI velocity integration range:')
    print('%.1f to %.1f km/s' % (vel_range_confint[0],
                                 vel_range_confint[1]))
    print('DGR:')
    print('%.1f to %.1f km/s' % (dgr_confint[0],
                                 dgr_confint[1]))

    global_props['dust2gas_ratio'] = {}
    global_props['dust2gas_ratio_error'] = {}

    global_props['hi_vel_range'] = vel_range_confint[0:2]
    global_props['hi_vel_range_error'] = vel_range_confint[2:]
    global_props['dust2gas_ratio']['value'] = dgr_confint[0]
    global_props['dust2gas_ratio_error']['value'] = dgr_confint[1:]
    global_props['hi_vel_range_conf'] = conf
    global_props['center_likelihood'] = center_likelihood.tolist()
    global_props['width_likelihood'] = width_likelihood.tolist()
    global_props['dgr_likelihood'] = dgr_likelihood.tolist()
    global_props['vel_centers'] = velocity_centers.tolist()
    global_props['vel_widths'] = vel_widths.tolist()
    global_props['dgrs'] = dgrs.tolist()
    global_props['likelihoods'] = likelihoods.tolist()

    with open(property_dir + global_property_file, 'w') as f:
        json.dump(global_props, f)
    '''

if __name__ == '__main__':
    main()



