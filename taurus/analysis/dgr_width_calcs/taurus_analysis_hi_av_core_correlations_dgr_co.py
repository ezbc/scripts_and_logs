#!/usr/bin/python

''' Calculates the N(HI) / Av likelihoodelation for the taurus molecular cloud.
'''

import pyfits as pf
import numpy as np

class KeyboardInterruptError(Exception): pass

''' Plotting Functions
'''
def plot_likelihoods(likelihoods,velocity_centers,velocity_widths,
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
                                       velocity_widths.shape[0]))
        likelihoods[:,:] = np.NaN
        count = 0
        try:
            for i, center in enumerate(velocity_centers):
                for j, width in enumerate(velocity_widths):
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

    #ax.set_xticks(np.arange(0,velocity_widths.shape[0],1)[::5],
    #        velocity_centers[::5])

    plt.rc('text', usetex=False)
    im = ax.imshow(image, interpolation='nearest', origin='lower',
            extent=[velocity_widths[0],velocity_widths[-1],
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
            extent=[velocity_widths[0],velocity_widths[-1],
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
        contour_confs=None):

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
        x_limits = (0, 20)
        x_limits = (0, 20)
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
        y_limits = (0.0, 0.8)
        y_limits = (0.1, 0.4)

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

def search_likelihoods(mesh):

    from myimage_analysis import calculate_nhi

    # calculate the likelihoodelation coefficient for each velocity range

    # Progress bar parameters
    #total = float(likelihoods.size)
    #count = 0

    try:
        velocity_center = mesh[0]
        velocity_width = mesh[1]
        dgr = mesh[2]

        velocity_range = (velocity_center - velocity_width / 2.,
                          velocity_center + velocity_width / 2.)

        nhi_image_temp, nhi_image_error = \
                calculate_nhi(cube=hi_cube,
                    velocity_axis=hi_velocity_axis,
                    velocity_range=velocity_range,
                    noise_cube=hi_noise_cube)

        # Avoid NaNs
        indices = np.where((nhi_image_temp == nhi_image_temp) & \
                           (av_image == av_image) & \
                           (nhi_image_temp > 0))

        nhi_image_likelihood = nhi_image_temp[indices]
        nhi_image_error_likelihood = nhi_image_error[indices]
        av_image_likelihood = av_image[indices]
        if type(av_image_error) != float:
            av_image_error_likelihood = np.median(av_image_error[indices])
        else:
            av_image_error_likelihood = np.median(av_image_error)

        # Create model of Av with N(HI) and DGR
        av_image_model = nhi_image_likelihood * dgr

        logL = calc_logL(av_image_model,
                         av_image_likelihood,
                         data_error=av_image_error_likelihood)

        likelihood = logL

        return likelihood
    except KeyboardInterrupt:
        raise KeyboardInterruptError

def setup_likelihood_mesh(velocity_centers, velocity_widths, dgrs):
    mesh = np.meshgrid(velocity_centers, velocity_widths, dgrs)
    mesh = np.vstack(mesh).reshape(3, -1).T

    return mesh

def reshape_likelihoods(likelihoods, velocity_centers=None,
        velocity_widths=None, dgrs=None):

    likelihoods = np.asarray(likelihoods).T.reshape(len(velocity_centers),
                                                    len(velocity_widths),
                                                    len(dgrs))

    return likelihoods

def calc_likelihood_hi_av(#hi_cube=None, hi_velocity_axis=None,
        #hi_noise_cube=None, av_image=None, av_image_error=None,
        velocity_centers=None, velocity_widths=None, return_likelihoods=True,
        dgrs=None, plot_results=True, results_filename='',
        likelihood_filename=None, clobber=False, conf=0.68,
        contour_confs=None):

    '''
    Parameters
    ----------

    Returns
    -------
    hi_vel_range : tuple
        Lower and upper bound of HI velocity range in km/s which provides the
        best likelihoodelated N(HI) distribution with Av.
    likelihoods : array-like, optional
        Array of Pearson likelihoodelation coefficients likelihoodesponding to each
        permutation through the velocity centers and velocity widths.

    '''

    import numpy as np
    from scipy.stats import pearsonr
    from scipy.stats import kendalltau
    from myimage_analysis import calculate_nhi
    from scipy import signal
    from os import path
    from astropy.io import fits
    import multiprocessing
    import itertools

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

    if perform_mle:
        # calculate the velocity ranges given a set of centers and widths
        velocity_ranges = np.zeros(shape=[len(velocity_centers) * \
                len(velocity_widths),2])
        count = 0
        for i, center in enumerate(velocity_centers):
            for j, width in enumerate(velocity_widths):
                velocity_ranges[count, 0] = center - width/2.
                velocity_ranges[count, 1] = center + width/2.
                count += 1

        # Set up iterable whereby each row contains the parameter values
        mesh = setup_likelihood_mesh(velocity_centers, velocity_widths, dgrs)

        # Use multiple processors to iterate through parameter combinations
        p = multiprocessing.Pool()
        likelihoods = p.map(search_likelihoods, mesh)
        p.close()

        # reshape likelihoods
        likelihoods = reshape_likelihoods(likelihoods,
                                          velocity_centers=velocity_centers,
                                          velocity_widths=velocity_widths,
                                          dgrs=dgrs)

        # Normalize the log likelihoods
        likelihoods -= likelihoods.max()

        # Convert to likelihoods
        likelihoods = np.exp(likelihoods)

        # Normalize the likelihoods
        likelihoods = likelihoods / \
            np.sum(likelihoods[~np.isnan(likelihoods)])

        # Write out fits file of likelihoods
        if write_mle:
            write_mle_tofits(filename=likelihood_filename,
                             velocity_centers=velocity_centers,
                             velocity_widths=velocity_widths,
                             dgrs=dgrs,
                             likelihoods=likelihoods,
                             clobber=clobber)

        # Avoid nans
        likelihoods = np.ma.array(likelihoods,
                mask=(likelihoods != likelihoods))

    # Load file of likelihoods
    elif not perform_mle:
        print('Reading likelihood grid file:')
        print(likelihood_filename)

        hdu = fits.open(likelihood_filename)
        likelihoods = hdu[0].data

        if len(velocity_centers) != likelihoods.shape[0] or \
            len(velocity_widths) != likelihoods.shape[1]:
            raise ValueError('Specified parameter grid not the same as in' + \
                    'loaded data likelihoods.')

        likelihoods = np.ma.array(likelihoods,
                mask=(likelihoods != likelihoods))

    # Define parameter resolutions
    #delta_center = velocity_centers[1] - velocity_centers[0]
    #delta_width = velocity_widths[1] - velocity_widths[0]

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
    width_confint = threshold_area(velocity_widths,
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
        plot_likelihoods_hist(likelihoods,
                              velocity_widths,
                              dgrs,
                              x_confint=width_confint,
                              y_confint=dgr_confint,
                              plot_axes=('widths', 'dgrs'),
                              show=0,
                              returnimage=False,
                              filename=results_filename + '_wd.png',
                              contour_confs=contour_confs)

    if not return_likelihoods:
        return vel_range_confint, dgr_confint
    else:
        return (vel_range_confint, dgr_confint, likelihoods,
            center_likelihood, width_likelihood, dgr_likelihood)

def calc_logL(model, data, data_error=None):

    '''
    Calculates log likelihood

    http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node2.html

    '''

    if data_error is None:
        data_error = np.std(data)
    N = data.size

    logL = -np.sum((data - model)**2 / (2 * data_error**2)) / N

    return logL

def calc_model_chisq(params, av_image=None, av_image_error=None, hi_cube=None,
        hi_velocity_axis=None, hi_noise_cube=None, dgr=None):

    from myimage_analysis import calculate_nhi

    velocity_range = params['low_vel'].value, params['high_vel'].value

    nhi_image_temp, nhi_image_error = calculate_nhi(cube=hi_cube,
            velocity_axis=hi_velocity_axis,
            velocity_range=velocity_range,
            noise_cube=hi_noise_cube)

    # Select pixels with Av > 1.0 mag and Av_SNR > 5.0.
    # Av > 1.0 mag is used to avoid too low Av.
    # 1.0 mag likelihoodesponds to SNR = 1 / 0.2 ~ 5
    # (see Table 2 of Ridge et al. 2006).
    indices = np.where((nhi_image_temp == nhi_image_temp) & \
                       (av_image == av_image) & \
                       (av_image_error == av_image_error))

    nhi_image_likelihood = nhi_image_temp[indices]
    nhi_image_error_likelihood = nhi_image_error[indices]
    av_image_data = av_image[indices]
    av_image_error_data = av_image_error[indices]

    # Create model image
    av_image_model = nhi_image_likelihood * dgr
    av_image_error_model = nhi_image_error_likelihood * dgr

    chisq = np.sum((av_image_data - av_image_model)**2 / \
                   (av_image_error_data**2))

    return chisq

def fit_hi_vel_range(guesses=None, av_image=None, av_image_error=None,
        hi_cube=None, hi_velocity_axis=None, hi_noise_cube=None, dgr=None):

    from scipy.optimize import curve_fit
    from scipy import stats
    from lmfit import minimize, Parameters, report_fit, report_errors
    from pprint import pprint

    params = Parameters()
    params.add('low_vel',
               value=guesses[0],
               min=-100,
               max=100,
               vary=True)
    params.add('high_vel',
               value=guesses[1],
               min=-100,
               max=100,
               vary=True)

    result = minimize(calc_model_chisq,
                      params,
                      kws={'av_image': av_image,
                           'av_image_error': av_image_error,
                           'hi_cube': hi_cube,
                           'hi_velocity_axis': hi_velocity_axis,
                           'hi_noise_cube': hi_noise_cube,
                           'dgr': dgr},
                      #method='BFGS',
                      method='anneal',
                      #method='powell',
                      #method='SLSQP',
                      options={'gtol': 1e-6,
                               'disp': True,
                               'maxiter' : 1e9}
                      )

    report_fit(params)
    report_errors(params)
    print result.values
    print result.errorbars
    #print(result.__dict__)
    #print(dir(result))
    #print result.vars
    #print help(result)
    print result.view_items()
    print result.Bopt


    #report_fit(result)

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
        velocity_widths=None, dgrs=None, likelihoods=None, clobber=False):

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
    header['CRVAL2'] = velocity_widths[0]
    header['CRVAL3'] = dgrs[0]
    try:
        header['CDELT1'] = velocity_centers[1] - velocity_centers[0]
    except IndexError:
        header['CDELT1'] = 1
    try:
        header['CDELT2'] = velocity_widths[1] - velocity_widths[0]
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

def calc_co_noise(co_mom0, prop_dict):

    co_noise_region = []

    # Append pixels from each region to CO region map
    for region in prop_dict['co_noise_limits']['pixel']:
        co_noise_region.append(co_mom0[region[0][1]:region[1][1],
                                       region[0][0]:region[1][0]])

    # Calc noise
    noise = 0.0
    for region in co_noise_region:
    	std = np.std(np.array(region)[~np.isnan(region)])
    	noise += std

    # Take average of stds
    noise = noise / len(co_noise_region)

    return noise

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

def convert_limit_coordinates(prop_dict,
        coords=('region_limit', 'co_noise_limits'), header=None):

    # Initialize pixel keys
    for coord in coords:
        prop_dict[coord].update({'pixel': []})

        if coord == 'region_limit':
            limit_wcs = prop_dict[coord]['wcs']

            for limits in limit_wcs:
                # convert centers to pixel coords
                limit_pixels = get_pix_coords(ra=limits[0],
                                             dec=limits[1],
                                             header=header)[:2].tolist()

                prop_dict[coord]['pixel'].append(limit_pixels[0])
                prop_dict[coord]['pixel'].append(limit_pixels[1])
        elif coord == 'co_noise_limits':
            region_limits = prop_dict[coord]['wcs']

            # Cycle through each region, convert WCS limits to pixels
            for region in region_limits:
                region_pixels = []
                for limits in region:
                    # convert centers to pixel coords
                    limit_pixels = get_pix_coords(ra=limits[0],
                                                  dec=limits[1],
                                                  header=header)[:2].tolist()
                    region_pixels.append(limit_pixels)

                # Append individual regions back to CO noise
                prop_dict[coord]['pixel'].append(region_pixels)

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
            box_height = region[2] / header['CDELT1']
            box_width = region[3] / header['CDELT2']
            cores[core].update({'box_center_pix': box_center_pixel})
            cores[core].update({'box_width': box_width})
            cores[core].update({'box_height': box_height})
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

    global hi_cube
    global hi_velocity_axis
    global hi_noise_cube
    global av_image
    global av_image_error

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

    # Confidence of parameter errors
    conf = 0.68
    # Confidence of contour levels
    contour_confs = (0.68, 0.95)

    # Course, large grid or fine, small grid?
    grid_res = 'course'
    grid_res = 'fine'

    # Results and fits filenames
    likelihood_filename = 'taurus_nhi_av_likelihoods_co_only'
    results_filename = 'taurus_likelihood_co_only'

    # Name of property files results are written to
    global_property_file = 'taurus_global_properties.txt'
    core_property_file = 'taurus_core_properties.txt'

    # Threshold of Av below which we expect only atomic gas, in mag
    av_thres = 20

    # Name of noise cube
    noise_cube_filename = 'taurus_hi_galfa_cube_regrid_planckres_noise.fits'

    # Define ranges of parameters
    if center_vary and width_vary and dgr_vary:
        likelihood_filename += '_width_dgr_center'
        results_filename += '_width_dgr_center'

        velocity_centers = np.arange(-15, 30, 1)
        velocity_widths = np.arange(1, 80, 1)
        dgrs = np.arange(1e-2, 1, 2e-2)
    elif not center_vary and width_vary and dgr_vary:
        if grid_res == 'course':
            likelihood_filename += '_dgr_width_lowres'
            results_filename += '_dgr_width_lowres'
            velocity_centers = np.arange(5, 6, 1)
            velocity_widths = np.arange(1, 15, 0.16667)
            dgrs = np.arange(0.05, 0.7, 0.5e-3)
        elif grid_res == 'fine':
            likelihood_filename += '_dgr_width_highres'
            results_filename += '_dgr_width_highres'
            velocity_centers = np.arange(5, 6, 1)
            velocity_widths = np.arange(1, 100, 0.16667)
            dgrs = np.arange(0.15, 0.4, 1e-3)
            velocity_widths = np.arange(1, 20, 0.16667)
            dgrs = np.arange(0.1, 0.9, 2e-3)
            #velocity_widths = np.arange(1, 40, 1)
            #dgrs = np.arange(0.15, 0.4, 1e-1)
    elif center_vary and width_vary and not dgr_vary:
        likelihood_filename += '_width_center'
        results_filename += '_width_center'

        velocity_centers = np.arange(-15, 30, 1)
        velocity_widths = np.arange(1, 80, 1)
        dgrs = np.arange(1.1e-1, 1.2e-1, 0.1e-1)
    elif not center_vary and width_vary and not dgr_vary:
        likelihood_filename += '_width'
        results_filename += '_width'

        velocity_centers = np.arange(5, 6, 1)
        velocity_widths = np.arange(1, 80, 1)
        dgrs = np.arange(1.1e-1, 1.2e-1, 0.1e-1)

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

    # load Planck Av and GALFA HI images, on same grid
    av_data_planck, av_header = load_fits(av_dir + \
                'taurus_av_planck_5arcmin.fits',
            return_header=True)

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
    co_velocity_axis = make_velocity_axis(co_header)

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

    # Change WCS coords to pixel coords of images
    cores = convert_core_coordinates(cores, h)
    cores = load_ds9_region(cores,
            filename_base = region_dir + 'taurus_av_boxes_',
            header = h)
    global_props = convert_limit_coordinates(global_props, header=av_header)

    print('\nCalculating likelihoods globally')

    # Set velocity center as CO peak
    co_data_nonans = np.copy(co_data)
    co_data_nonans[np.isnan(co_data_nonans)] = 0.0
    co_spectrum = np.sum(co_data_nonans, axis=(1,2))
    co_avg_vel = np.average(co_velocity_axis, weights=co_spectrum)
    co_peak_vel = co_velocity_axis[co_spectrum == np.max(co_spectrum)]
    #velocity_centers = np.arange(co_peak_vel, co_peak_vel + 1, 1)
    velocity_centers = np.arange(co_avg_vel, co_avg_vel + 1, 1)

    print('\nVelocity center from CO = ' +\
            '{0:.2f} km/s'.format(velocity_centers[0]))

    # Create mask where CO is present
    core_mask = np.zeros(av_data_planck.shape)
    #for core in cores:
    #    # Grab the mask
    #    core_mask += myg.get_polygon_mask(av_data_planck,
    #            cores[core]['box_vertices_rotated'])

    # Calc moment 0 map of CO
    co_mom0 = np.sum(co_data_nonans, axis=0)

    # calc noise without any emission
    co_noise = calc_co_noise(co_mom0, global_props)
    co_thres = 2.0 * co_noise
    co_thres = 6.0 # K km/s

    pix = global_props['region_limit']['pixel']
    region_vertices = ((pix[0], pix[1]),
                       (pix[1], pix[2]),
                       (pix[2], pix[3]),
                       (pix[3], pix[0])
                       )

    # block offregion
    region_mask = myg.get_polygon_mask(av_data_planck, region_vertices)

    # Get indices which trace only atomic gas, i.e., no CO emission
    indices = ((co_mom0 < co_thres) & \
               (av_data_planck < av_thres) & \
               (region_mask == 1))

    print('\nTotal number of pixels in analysis = ' + \
            '{0:.0f}'.format(indices[indices].size)) + \
            '\ngiven a CO threshold of {0:.2f} K km/s'.format(co_thres)

    # Mask global data with CO indices
    hi_data_sub = np.copy(hi_data[:, indices])
    noise_cube_sub = np.copy(noise_cube[:, indices])
    av_data_sub = np.copy(av_data_planck[indices])
    av_error_data_sub = np.copy(av_error_data_planck[indices])

    # Plot the masked image
    av_data_masked = np.copy(av_data_planck)
    av_data_masked[~indices] = np.nan
    figure_types = ['png',]
    for figure_type in figure_types:
        plot_av_image(av_image=av_data_masked, header=av_header,
                savedir=figure_dir + '../maps/',
                limits=global_props['region_limit']['pixel'],
                filename='taurus_dgr_co_masked_map.' + figure_type,
                show=0)

    # Set global variables
    hi_cube = hi_data_sub
    hi_velocity_axis = velocity_axis
    hi_noise_cube = noise_cube_sub
    av_image = av_data_sub
    av_image_error = av_error_data_sub

    # Define filename for plotting results
    results_filename = figure_dir + results_filename

    # likelihoodelate each core region Av and N(HI) for velocity ranges
    pool = Pool()
    vel_range_confint, dgr_confint, likelihoods, center_likelihood,\
        width_likelihood, dgr_likelihood = \
            calc_likelihood_hi_av(
                            dgrs=dgrs,
                            velocity_centers=velocity_centers,
                            velocity_widths=velocity_widths,
                            return_likelihoods=True,
                            plot_results=True,
                            results_filename=results_filename,
                            likelihood_filename=likelihood_dir + \
                                    likelihood_filename + \
                                    '_global.fits',
                            clobber=clobber,
                            conf=conf,
                            contour_confs=contour_confs)

    print('HI velocity integration range:')
    print('%.1f to %.1f km/s' % (vel_range_confint[0],
                                 vel_range_confint[1]))
    print('DGR:')
    print('%.1f to %.1f km/s' % (dgr_confint[0],
                                 dgr_confint[1]))

    global_props['dust2gas_ratio'] = {}
    global_props['dust2gas_ratio_error'] = {}

    global_props['hi_velocity_range'] = vel_range_confint[0:2]
    global_props['hi_velocity_range_error'] = vel_range_confint[2:]
    global_props['dust2gas_ratio']['value'] = dgr_confint[0]
    global_props['dust2gas_ratio_error']['value'] = dgr_confint[1:]
    global_props['hi_velocity_range_conf'] = conf
    global_props['center_likelihood'] = center_likelihood.tolist()
    global_props['width_likelihood'] = width_likelihood.tolist()
    global_props['dgr_likelihood'] = dgr_likelihood.tolist()
    global_props['vel_centers'] = velocity_centers.tolist()
    global_props['vel_widths'] = velocity_widths.tolist()
    global_props['dgrs'] = dgrs.tolist()
    global_props['likelihoods'] = likelihoods.tolist()

    with open(property_dir + global_property_file, 'w') as f:
        json.dump(global_props, f)

if __name__ == '__main__':
    main()








