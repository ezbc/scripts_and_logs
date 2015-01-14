import matplotlib
matplotlib.use('Agg')

import warnings
warnings.filterwarnings('ignore')

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

def plot_likelihoods_hist(global_props, filename=None, show=True,
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

    if plot_axes[0] == 'widths':
    	x_grid = global_props['vel_widths']
    	x_confint = (global_props['hi_velocity_width']['value'],
    	             global_props['hi_velocity_width_error']['value'][0],
    	             global_props['hi_velocity_width_error']['value'][1],
    	             )
        x_extent = x_grid[0], x_grid[-1]
        ax_image.set_xlabel(r'Velocity Width (km/s)')
        x_sum_axes = 1
        y_pdf_label = r'Width PDF'
        x_limits = (x_grid[0], x_grid[-1])
    if plot_axes[1] == 'dgrs':
    	y_grid = global_props['dgrs']
    	y_confint = (global_props['dust2gas_ratio']['value'],
    	             global_props['dust2gas_ratio_error']['value'][0],
    	             global_props['dust2gas_ratio_error']['value'][1],
    	             )
        y_extent = y_grid[0], y_grid[-1]
        ax_image.set_ylabel(r'DGR (10$^{-20}$ cm$^2$ mag$^1$)')
        y_sum_axes = 0
        x_pdf_label = r'DGR PDF'
        y_limits = (y_grid[0], y_grid[-1])

    # Create axes
    sum_axes = np.array((x_sum_axes, y_sum_axes))
    sum_axis = np.argmax(np.bincount(np.ravel(sum_axes)))

    # Mask NaNs
    likelihoods = global_props['likelihoods']
    image = np.ma.array(likelihoods, mask=np.isnan(likelihoods))

    # Create likelihood image
    #image = np.sum(likelihoods, axis=sum_axis) / np.sum(likelihoods)
    image = likelihoods / np.sum(likelihoods)

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
                      drawstyle='steps-mid',
                      linewidth=2,
                      )

        ax_pdf_y.plot(y_pdf,
                      y_grid,
                      color='k',
                      drawstyle='steps-mid',
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

    if 0:
    #if npix is not None or av_threshold is not None:
    	text = ''
        if npix is not None:
            text += r'N$_{\rm pix}$ = ' + \
                     '{0:.0f}'.format(npix)
            if av_threshold is not None:
            	text += '\n'
        if av_threshold is not None:
            text += r'$A_V$ threshold = {0:.1f} mag'.format(av_threshold)
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

def calc_logL(model, data, data_error=None):

    '''
    Calculates log likelihood

    http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node2.html

    '''

    import numpy as np

    if data_error is None:
        data_error = np.std(data)

    logL = -np.sum((data - model)**2 / (2 * (data_error)**2))

    return logL

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

def write_mle_tofits(filename='', vel_widths=None, dgrs=None,
        likelihoods=None, clobber=False):

    from astropy.io import fits

    print('Writing likelihood grid to file:')
    print(filename)

    header = fits.Header()
    header['NAXIS'] = 2
    header['CTYPE1'] = 'WIDTHS'
    header['CTYPE2'] = 'DGR'
    header['CRPIX1'] = 0
    header['CRPIX2'] = 0
    header['CRVAL1'] = vel_widths[0]
    header['CRVAL2'] = dgrs[0]

    try:
        header['CDELT1'] = vel_widths[1] - vel_widths[0]
    except IndexError:
        header['CDELT1'] = 1
    try:
        header['CDELT2'] = dgrs[1] - dgrs[0]
    except IndexError:
        header['CDELT2'] = 1

    fits.writeto(filename,
                 likelihoods,
                 header,
                 clobber=clobber)

def gauss(x, width, amp, x0):
    import numpy as np

    return amp * np.exp(-(x - x0)**2 / (2 * width**2))

def get_residual_mask(residuals, resid_width_scale=3.0, plot_progress=False):

    '''

    '''

    import numpy as np
    from scipy.optimize import curve_fit

    # Fit the rising portion of the residuals
    residuals_crop = residuals[(residuals < 0) & ~np.isnan(residuals)]

    counts, bin_edges = np.histogram(np.ravel(residuals_crop),
                                     bins=100,
                                     )
    fit_params = curve_fit(gauss,
                           bin_edges[:-1],
                           counts,
                           p0=(2, np.nanmax(counts), 0),
                           maxfev=10000,
                           )[0]

    # Include only residuals within 3 sigma
    mask = residuals > resid_width_scale * np.abs(fit_params[0])

    if plot_progress:
        import matplotlib.pyplot as plt
        x_fit = np.linspace(-10,
                            10,
                            1000)

        y_fit = gauss(x_fit, *fit_params)

        counts, bin_edges = \
            np.histogram(np.ravel(residuals[~np.isnan(residuals)]),
                                     bins=100,
                                     )

        bin_edges_ext = np.zeros(len(counts) + 1)
        counts_ext = np.zeros(len(counts) + 1)

        bin_edges_ext[0] = bin_edges[0] - (bin_edges[1] - bin_edges[0])
        bin_edges_ext[1:] = bin_edges[:-1]
        counts_ext[0] = 0
        counts_ext[1:] = counts

        counts_ext /= np.nanmax(counts_ext)
        y_fit /= np.max(y_fit)

        plt.close;plt.clf()
        plt.plot(bin_edges_ext, counts_ext, drawstyle='steps-mid')
        plt.plot(x_fit, y_fit, color='r')
        plt.xlim([np.nanmin(bin_edges_ext) - \
                  np.abs(0.8 * np.nanmin(bin_edges_ext)),4])
        plt.ylim([-0.1, 1.1])
        plt.axvline(resid_width_scale * np.abs(fit_params[0]),
                    color='k',
                    linestyle='--',
                    linewidth=3)
        plt.xlabel(r'Residual $A_V$ [mag]')
        plt.ylabel('Normalized PDF')
        plt.show()

    return mask

def iterate_residual_masking(
                             nhi_image=None,
                             av_data=None,
                             av_data_error=None,
                             init_mask=None,
                             vel_range=None,
                             threshold_delta_dgr=None,
                             resid_width_scale=3.0,
                             plot_progress=False,
                             verbose=False,
                             ):

    '''

    Returns
    -------
    av_model : numpy array
    mask : numpy array
    dgr : float

    '''


    import numpy as np

    # Mask out nans
    mask = (np.isnan(av_data) | \
            np.isnan(av_data_error) | \
            (av_data_error == 0) | \
            np.isnan(nhi_image))

    if init_mask is not None:
        mask += init_mask

    # solve for DGR using linear least squares
    print('\nBeginning iterative DGR calculations + masking...')

    # Iterate masking pixels which are correlated and rederiving a linear least
    # squares solution for the DGR
    # -------------------------------------------------------------------------
    delta_dgr = 1e10
    dgr = 1e10
    while delta_dgr > threshold_delta_dgr:
        A = np.array((np.ravel(nhi_image[~mask] / av_data_error[~mask]),))
        b = np.array((np.ravel(av_data[~mask] / av_data_error[~mask]),))
        A = np.matrix(A).T
        b = np.matrix(b).T
        dgr_new = (np.linalg.pinv(A) * b)[0, 0]

        # Create model with the DGR
        if verbose:
            print('\tDGR = {0:.2} 10^20 cm^2 mag'.format(dgr))
        av_image_model = nhi_image * dgr_new

        residuals = av_data - av_image_model

        # Include only residuals which are white noise
        mask_new = get_residual_mask(residuals,
                                     resid_width_scale=resid_width_scale,
                                     plot_progress=plot_progress)

        # Mask non-white noise, i.e. correlated residuals.
        mask[mask_new] = 1

        if verbose:
            npix = mask.size - np.sum(mask)
            print('\tNumber of non-masked pixels = {0:.0f}'.format(npix))

        # Reset while loop conditions
        delta_dgr = np.abs(dgr - dgr_new)
        dgr = dgr_new

    # Create model of Av
    av_model = dgr * nhi_image
    av_model[mask] = np.nan

    return (av_model, mask, dgr)

def calc_likelihoods(
        hi_cube=None,
        hi_vel_axis=None,
        av_image=None,
        av_image_error=None,
        vel_center=None,
        vel_widths=None,
        dgrs=None,
        plot_results=False,
        results_filename='',
        return_likelihoods=True,
        likelihood_filename=None,
        clobber=False,
        conf=0.68,
        threshold_delta_dgr=0.0005,
        ):

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
    from myimage_analysis import calculate_nhi
    from os import path
    from astropy.io import fits

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
        # calculate the likelihoodelation coefficient for each velocity
        # range
        likelihoods = np.zeros((len(vel_widths),
                                len(dgrs)))

        # Progress bar parameters
        total = float(likelihoods.size)
        count = 0

        for j, vel_width in enumerate(vel_widths):
            # Construct N(HI) image outside of DGR loop, then apply
            # DGRs in loop
            vel_range = (vel_center - vel_width / 2.,
                         vel_center + vel_width / 2.)

            nhi_image = calculate_nhi(cube=hi_cube,
                                      velocity_axis=hi_vel_axis,
                                      velocity_range=vel_range,
                                      return_nhi_error=False)

            # Cycle through DGR to estimate error
            for k, dgr in enumerate(dgrs):
                # Create model of Av with N(HI) and DGR
                av_image_model = nhi_image * dgr
                av_image_model_error = nhi_image * dgr

                logL = calc_logL(av_image_model,
                                 av_image,
                                 data_error=av_image_error)

                likelihoods[j, k] = logL

                # Shows progress each 10%
                count += 1
                abs_step = int((total * 1)/10) or 10
                if count and not count % abs_step:
                    print "\t{0:.0%} processed".format(count/total)

    # Load file of likelihoods
    elif not perform_mle:
        print('Reading likelihood grid file:')
        print(likelihood_filename)

        hdu = fits.open(likelihood_filename)
        likelihoods = hdu[0].data

        if len(vel_widths) != likelihoods.shape[0] or \
           len(dgrs) != likelihoods.shape[1]:
            raise ValueError('Specified parameter grid not the same as in' + \
                    'loaded data likelihoods.')

        likelihoods = np.ma.array(likelihoods,
                mask=(likelihoods != likelihoods))

    # Normalize the log likelihoods
    likelihoods -= likelihoods.max()

    # Convert to likelihoods
    likelihoods = np.exp(likelihoods)

    # Normalize the likelihoods
    likelihoods = likelihoods / np.nansum(likelihoods)

    # Derive marginal distributions of both centers and widths
    width_likelihood = np.sum(likelihoods, axis=1) / \
            np.sum(likelihoods)
    dgr_likelihood = np.sum(likelihoods, axis=0) / \
            np.sum(likelihoods)

    # Derive confidence intervals of parameters
    width_confint = threshold_area(vel_widths,
                                   width_likelihood,
                                   area_fraction=conf)
    dgr_confint = threshold_area(dgrs,
                                 dgr_likelihood,
                                 area_fraction=conf)
    from mystats import calc_symmetric_error

    width_confint = calc_symmetric_error(vel_widths,
                                   width_likelihood,
                                   alpha=1.0 - conf)
    dgr_confint = calc_symmetric_error(dgrs,
                                 dgr_likelihood,
                                 alpha=1.0 - conf)

    # Get values of best-fit model parameters
    max_loc = np.where(likelihoods == np.max(likelihoods))
    width_max = vel_widths[max_loc[0]][0]
    dgr_max = dgrs[max_loc[1]][0]

    print('\nVelocity widths = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(width_confint[0],
                                                    width_confint[2],
                                                    np.abs(width_confint[1])))
    print('\nDGRs = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} 10^20 cm^2 mag'.format(dgr_confint[0],
                                                    dgr_confint[2],
                                                    np.abs(dgr_confint[1])))

    # Write PDF
    upper_lim = (vel_center + width_confint[0]/2.)
    lower_lim = (vel_center - width_confint[0]/2.)
    upper_lim_error = width_confint[2]**2
    lower_lim_error = width_confint[1]**2

    vel_range_confint = (lower_lim, upper_lim, lower_lim_error, upper_lim_error)
    vel_range_max = (vel_center - width_max/2.0, vel_center + width_max/2.0)

    if not return_likelihoods:
        return vel_range_confint, dgr_confint
    else:
        return (vel_range_confint, width_confint, dgr_confint, likelihoods,
                width_likelihood, dgr_likelihood, width_max, dgr_max,
                vel_range_max)

def rebin_image(image, bin_size):

    ''' From stack overflow

    http://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array

    '''

    import numpy as np

    shape = (image.shape[0] / bin_size, image.shape[1] / bin_size)

    # Crop image for binning
    image = image[:-image.shape[0] % shape[0], :-image.shape[1] % shape[1]]

    # Rebin image
    sh = shape[0],image.shape[0]//shape[0],shape[1],image.shape[1]//shape[1]
    return image.reshape(sh).mean(-1).mean(1)

def bin_image(image, in_ext='', out_ext='_bin', width=1, clobber=True,
        bin_dim=(0, 1), image_dir='./', verbose=False):

    ''' Bins and smooths image with MIRIAD imbin and imsmooth.

    File must end in '.fits', leave out from in_ext and out_ext

    width is in degrees

    '''

    from mirpy import fits, smooth, imbin
    import os

    # Change directories so that string is not too long for miriad
    current_dir = os.getcwd()
    os.chdir(image_dir)

    if not check_file(image + in_ext + '.mir',
                      clobber=clobber,
                      verbose=verbose):
        fits(image + in_ext + '.fits',
             out=image + in_ext + '.mir',
             op='xyin')

    # Determine size of beam to smooth with
    conv_beam = (width**2 - (5.0/60.0)**2)

    # Smooth the image
    if not check_file(image + in_ext + '_smooth.mir',
                      clobber=clobber,
                      verbose=verbose):
        smooth(image + in_ext + '.mir',
               out=image + in_ext + '_smooth.mir',
               fwhm=conv_beam,
               pa=0,
               scale=0.0)

    if not check_file(image + in_ext + '_smooth.fits',
                      clobber=clobber,
                      verbose=verbose):
        fits(image + in_ext + '_smooth.mir',
             out=image + in_ext + '_smooth.fits',
             op='xyout')

    # Determine number of pixels to bin
    binsize = width * 60.0 / 5.0

    if bin_dim == (0, 1):
        bins = 4*(binsize,)
    elif bin_dim == (1, 2):
        bins = (binsize, binsize, binsize, binsize, 1, 1)

    if not check_file(image + out_ext + '.mir',
                      clobber=clobber,
                      verbose=verbose):
        imbin(image + in_ext + '.mir',
              out=image + out_ext + '.mir',
              bin=bins,
              options='sum')

    if not check_file(image + out_ext + '.fits',
                      clobber=clobber,
                      verbose=verbose):
        fits(image + out_ext + '.mir',
             out=image + out_ext + '.fits',
             op='xyout')

    # Go back to original working directory
    os.chdir(current_dir)

    return binsize

def check_file(filename, clobber=False, verbose=False):

    import os

    exists = False

    if os.path.isfile(filename) or os.path.isdir(filename):
        exists = True
        if verbose:
            print('\tImage {:s} exists'.format(filename))
        if clobber:
            if verbose:
                print('\tDeleting image {:s}'.format(filename))
            os.system('rm -rf {:s}'.format(filename))
            exists = False

    return exists

''' DS9 Region and Coordinate Functions
'''

def convert_limit_coordinates(prop_dict,
        coords=('region_limit', 'co_noise_limits', 'plot_limit'), header=None):

    # Initialize pixel keys
    for coord in coords:
        prop_dict[coord].update({'pixel': []})

        if coord == 'region_limit' or coord == 'plot_limit':
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

def load_ds9_region(cores, filename_base = 'perseus_av_boxes_', header=None):

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
Main Script
'''

def main(av_data_type='planck'):

    # Import external modules
    # -----------------------
    import numpy as np
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json
    from myimage_analysis import calculate_nhi, calculate_noise_cube, bin_image
    #from astropy.io import fits
    import pyfits as fits
    import matplotlib.pyplot as plt

    # Set parameters
    # --------------
    # Check if likelihood file already written, rewrite?
    clobber = 1

    # Confidence of parameter errors
    conf = 0.68
    # Confidence of contour levels
    contour_confs = (0.95,)

    # Name of HI noise cube
    noise_cube_filename = 'perseus_hi_galfa_cube_regrid_planckres_noise'

    # Threshold for converging DGR
    threshold_delta_dgr = 0.00005

    # Number of white noise standard deviations with which to fit the
    # residuals in iterative masking
    resid_width_scale = 3.0

    # Name of property files results are written to
    global_property_file = 'perseus_global_properties.txt'

    # Likelihood axis resolutions
    vel_widths = np.arange(1, 20, 2*0.16667)
    dgrs = np.arange(0.01, 0.2, 3e-4)
    #vel_widths = np.arange(1, 50, 8*0.16667)
    #dgrs = np.arange(0.01, 0.2, 1e-2)

    # Velocity range over which to integrate HI for deriving the mask
    vel_range = (-2.3, 9)

    # Bin width in degrees
    bin_width_deg = 1.0

    # Clobber the binned images and remake them?
    clobber_bin_images = 0

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'
    figure_dir = \
        '/d/bip3/ezbc/perseus/figures/'
    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    hi_dir = '/d/bip3/ezbc/perseus/data/hi/'
    co_dir = '/d/bip3/ezbc/perseus/data/co/'
    core_dir = '/d/bip3/ezbc/perseus/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/perseus/data/python_output/'
    region_dir = '/d/bip3/ezbc/perseus/data/python_output/ds9_regions/'
    likelihood_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'

    # Load data
    # ---------
    # Adjust filenames
    #noise_cube_filename += bin_string
    likelihood_filename = 'perseus_likelihood_{0:s}_bin'.format(av_data_type)
    results_filename = 'perseus_likelihood_{0:s}_bin'.format(av_data_type)

    av_data, av_header = fits.getdata(av_dir + \
                            'perseus_av_planck_5arcmin.fits',
                                      header=True)

    av_data_error, av_error_header = fits.getdata(av_dir + \
                'perseus_av_error_planck_5arcmin.fits',
            header=True)

    hi_data, hi_header = fits.getdata(hi_dir + \
                'perseus_hi_galfa_cube_regrid_planckres.fits',
            header=True)

    # Load global properties
    with open(property_dir + global_property_file, 'r') as f:
        global_props = json.load(f)

    # Prepare data products
    # ---------------------
    # Change WCS coords to pixel coords of images
    global_props = convert_limit_coordinates(global_props, header=av_header)

    # make the velocity axes
    hi_vel_axis = make_velocity_axis(hi_header)

    # Derive relevant region
    pix = global_props['region_limit']['pixel']
    region_vertices = ((pix[1], pix[0]),
                       (pix[1], pix[2]),
                       (pix[3], pix[2]),
                       (pix[3], pix[0])
                       )

    # block off region
    region_mask = np.logical_not(myg.get_polygon_mask(av_data, region_vertices))

    print('\nRegion size = ' + \
          '{0:.0f} pix'.format(region_mask[region_mask == 1].size))

    # Derive mask by excluding correlated residuals
    # ---------------------------------------------
    nhi_image = calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=vel_range,
                              return_nhi_error=False,
                              )

    print('\nDeriving mask for correlated residuals...')

    av_model, mask, dgr = iterate_residual_masking(
                             nhi_image=nhi_image,
                             av_data=av_data,
                             av_data_error=av_data_error,
                             vel_range=vel_range,
                             threshold_delta_dgr=threshold_delta_dgr,
                             resid_width_scale=resid_width_scale,
                             init_mask=region_mask,
                             verbose=1,
                             plot_progress=0,
                             )

    # Combine region mask with new mask
    #mask += np.logical_not(region_mask)
    mask += region_mask

    if 0:
        import matplotlib.pyplot as plt
        plt.imshow(np.ma.array(av_data, mask=mask), origin='lower')
        plt.show()

    # Bin the masked images
    # ---------------------
    print('\nBinning masked images...')

    av_data[mask] = np.nan
    av_data_error[mask] = np.nan
    hi_data[:, mask] = np.nan
    #av_data[mask] = 0
    #av_data_error[mask] = 0
    #hi_data[:, mask] = 0

    if not check_file(av_dir + 'perseus_av_planck_5arcmin_masked.fits',
                      clobber=clobber_bin_images):
        fits.writeto(av_dir + 'perseus_av_planck_5arcmin_masked.fits',
                     av_data,
                     av_header)

    if not check_file(av_dir + 'perseus_av_error_planck_5arcmin_masked.fits',
                      clobber=clobber_bin_images):
        fits.writeto(av_dir + 'perseus_av_error_planck_5arcmin_masked.fits',
                     av_data_error,
                     av_header)

    if not check_file(hi_dir + \
                      'perseus_hi_galfa_cube_regrid_planckres_masked.fits',
                      clobber=clobber_bin_images):
        fits.writeto(hi_dir + \
                     'perseus_hi_galfa_cube_regrid_planckres_masked.fits',
                     hi_data,
                     hi_header)

    if clobber_bin_images:
        # Define number of pixels in each bin
        binsize = bin_width_deg * 60.0 / 5.0

        # Bin the images
        # Av image
        av_data_bin, av_header_bin = bin_image(av_data,
                                       binsize=(binsize, binsize),
                                       header=av_header,
                                       func=np.nanmean)

        if not check_file(av_dir + 'perseus_av_planck_5arcmin_bin.fits',
                          clobber=clobber_bin_images):
            fits.writeto(av_dir + 'perseus_av_planck_5arcmin_bin.fits',
                         av_data_bin,
                         av_header_bin)

        # Av image error
        # Errors add in square
        # mean = sum(a_i) / n
        # error on mean = sqrt(sum(a_i**2 / n**2))
        noise_func = lambda x: np.nansum(x**2)**0.5 / x[~np.isnan(x)].size

        av_data_error_bin, av_header_bin = bin_image(av_data_error,
                                             binsize=(binsize, binsize),
                                             header=av_header,
                                             func=noise_func,)

        if not check_file(av_dir + 'perseus_av_error_planck_5arcmin_bin.fits',
                          clobber=clobber_bin_images):
            fits.writeto(av_dir + 'perseus_av_error_planck_5arcmin_bin.fits',
                         av_data_error_bin,
                         av_header_bin)

        # Hi image
        hi_data_bin, hi_header_bin = bin_image(hi_data,
                                       binsize=(binsize, binsize),
                                       header=hi_header,
                                       func=np.nanmean)

        if not check_file(hi_dir + \
                          'perseus_hi_galfa_cube_regrid_planckres_bin.fits',
                          clobber=clobber_bin_images):
            fits.writeto(hi_dir + \
                         'perseus_hi_galfa_cube_regrid_planckres_bin.fits',
                         hi_data_bin,
                         hi_header_bin)

    # Load data
    # ---------
    bin_string = '_bin'

    # Adjust filenames
    noise_cube_filename += bin_string
    likelihood_filename = 'perseus_likelihood_{0:s}_bin'.format(av_data_type)
    results_filename = 'perseus_likelihood_{0:s}_bin'.format(av_data_type)

    av_data, av_header = fits.getdata(av_dir + \
                            'perseus_av_planck_5arcmin' + bin_string + '.fits',
                                      header=True)

    av_data_error, av_error_header = fits.getdata(av_dir + \
                'perseus_av_error_planck_5arcmin' + bin_string + '.fits',
            header=True)

    hi_data, hi_header = fits.getdata(hi_dir + \
                            'perseus_hi_galfa_cube_regrid_planckres' + \
                            bin_string + '.fits',
                            header=True)

    # Load global properties
    with open(property_dir + global_property_file, 'r') as f:
        global_props = json.load(f)

    # Prepare data products
    # ---------------------
    # Change WCS coords to pixel coords of images
    global_props = convert_limit_coordinates(global_props, header=av_header)

    # make the velocity axes
    hi_vel_axis = make_velocity_axis(hi_header)

    # Load the HI noise cube if it exists, else make it
    if not path.isfile(hi_dir + noise_cube_filename + '.fits'):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=hi_vel_axis,
                velocity_noise_range=[90,110], header=hi_header, Tsys=30.,
                filename=hi_dir + noise_cube_filename + '.fits')
    else:
        noise_cube, noise_header = fits.getdata(hi_dir +
                noise_cube_filename + '.fits',
            header=True)

    # Derive relevant region
    pix = global_props['region_limit']['pixel']
    region_vertices = ((pix[1], pix[0]),
                       (pix[1], pix[2]),
                       (pix[3], pix[2]),
                       (pix[3], pix[0])
                       )

    # block off region
    region_mask = np.logical_not(myg.get_polygon_mask(av_data,
                                                      region_vertices))

    print('\nRegion size = ' + \
          '{0:.0f} pix'.format(region_mask[region_mask == 1].size))

    # Mask out the NaNs
    mask = (np.isnan(av_data) & \
            np.isnan(av_data_error) & \
            np.isnan(np.sum(hi_data, axis=0)))
    mask += region_mask
    mask = mask.astype(bool)

    # Derive center velocity from hi
    # ------------------------------
    hi_spectrum = np.sum(hi_data[:, ~mask], axis=(1))
    vel_center = np.array((np.average(hi_vel_axis,
                           weights=hi_spectrum**2),))[0]
    print('\nVelocity center from HI = ' +\
            '{0:.2f} km/s'.format(vel_center))

    # Perform likelihood calculation of masked images
    # -----------------------------------------------
    # Define filename for plotting results
    results_filename = figure_dir + 'likelihood/'+ results_filename

    print('\nPerforming likelihood calculations with initial error ' + \
          'estimate...')
    results = calc_likelihoods(
                     hi_cube=hi_data[:, ~mask],
                     hi_vel_axis=hi_vel_axis,
                     av_image=av_data[~mask],
                     av_image_error=av_data_error[~mask],
                     vel_center=vel_center,
                     vel_widths=vel_widths,
                     dgrs=dgrs,
                     results_filename='',
                     return_likelihoods=True,
                     likelihood_filename=None,
                     clobber=False,
                     conf=conf,
                     )

    # Unpack output of likelihood calculation
    (vel_range_confint, width_confint, dgr_confint, likelihoods,
            width_likelihood, dgr_likelihood, width_max, dgr_max,
            vel_range_max) = results

    print('\nHI velocity integration range:')
    print('%.1f to %.1f km/s' % (vel_range_confint[0],
                                 vel_range_confint[1]))

    # Calulate chi^2 for best fit models
    # ----------------------------------
    nhi_image_temp, nhi_image_error = \
            calculate_nhi(cube=hi_data,
                velocity_axis=hi_vel_axis,
                velocity_range=vel_range_max,
                noise_cube=noise_cube,
                return_nhi_error=True)
    av_image_model = nhi_image_temp * dgr_max
    # avoid NaNs
    indices = ((av_image_model == av_image_model) & \
               (av_data == av_data))
    # add nan locations to the mask
    mask[~indices] = 1

    # count number of pixels used in analysis
    npix = mask[~mask].size

    # finally calculate chi^2
    chisq = np.sum((av_data[~mask] - av_image_model[~mask])**2 / \
            av_data_error[~mask]**2) / av_data[~mask].size

    print('\nTotal number of pixels in analysis, after masking = ' + \
            '{0:.0f}'.format(npix))

    print('\nReduced chi^2 = {0:.1f}'.format(chisq))

    # Write results to global properties
    global_props['dust2gas_ratio'] = {}
    global_props['dust2gas_ratio_error'] = {}
    global_props['hi_velocity_width'] = {}
    global_props['hi_velocity_width_error'] = {}
    global_props['dust2gas_ratio_max'] = {}
    global_props['hi_velocity_center_max'] = {}
    global_props['hi_velocity_width_max'] = {}
    global_props['hi_velocity_range_max'] =  {}
    global_props['av_threshold'] = {}
    global_props['co_threshold'] = {}
    global_props['hi_velocity_width']['value'] = width_confint[0]
    global_props['hi_velocity_width']['unit'] = 'km/s'
    global_props['hi_velocity_width_error']['value'] = width_confint[1:]
    global_props['hi_velocity_width_error']['unit'] = 'km/s'
    global_props['hi_velocity_range'] = vel_range_confint[0:2]
    global_props['hi_velocity_range_error'] = vel_range_confint[2:]
    global_props['dust2gas_ratio']['value'] = dgr_confint[0]
    global_props['dust2gas_ratio_error']['value'] = dgr_confint[1:]
    global_props['dust2gas_ratio_max']['value'] = dgr_max
    global_props['hi_velocity_center_max']['value'] = vel_center
    global_props['hi_velocity_width_max']['value'] = width_max
    global_props['hi_velocity_range_max']['value'] = vel_range_max
    global_props['hi_velocity_range_conf'] = conf
    global_props['width_likelihood'] = width_likelihood.tolist()
    global_props['dgr_likelihood'] = dgr_likelihood.tolist()
    global_props['vel_centers'] = [vel_center,]
    global_props['vel_widths'] = vel_widths.tolist()
    global_props['dgrs'] = dgrs.tolist()
    global_props['likelihoods'] = likelihoods.tolist()
    global_props['av_threshold']['value'] = None
    global_props['av_threshold']['unit'] = 'mag'
    global_props['co_threshold']['value'] = None
    global_props['co_threshold']['unit'] = 'K km/s'
    global_props['chisq'] = chisq
    global_props['npix'] = npix
    global_props['mask'] = mask.tolist()
    global_props['use_binned_image'] = True

    # Write the file
    print('\nWriting results to\n' + global_property_file + '_init.txt')

    with open(property_dir + global_property_file + '_init.txt', 'w') as f:
        json.dump(global_props, f)

    # Plot likelihood space
    print('\nWriting likelihood image to\n' + results_filename + \
            '_init_wd.png')
    plot_likelihoods_hist(global_props,
                          plot_axes=('widths', 'dgrs'),
                          show=0,
                          returnimage=False,
                          filename=results_filename + '_init_wd.png',
                          contour_confs=contour_confs)

    # Rerun analysis with new error calculated
    # Error should be calculated across entire image, not just atomic regions,
    # in order to understand variation in DGR
    # -------------------------------------------------------------------------
    # Calculate new standard deviation, set global variable
    # npix - 2 is the number of degrees of freedom
    # see equation 15.1.6 in Numerical Recipes
    std = np.sqrt(np.sum((av_data[~mask] - av_image_model[~mask])**2 \
                         / (av_data[~mask].size - 2)))

    av_data_error = std * np.ones(av_data_error.shape)
    #av_image_error += np.std(av_data[~mask] - av_image_model[~mask])

    print('\nSystematic error between model and data Av images:')
    print('\tstd(model - data) = {0:.3f} mag'.format(av_data_error[0, 0]))

    # Perform likelihood calculation of masked images
    # -----------------------------------------------
    print('\nPerforming likelihood calculations with scaled error ' + \
          'estimate...')
    results = calc_likelihoods(
                     hi_cube=hi_data[:, ~mask],
                     hi_vel_axis=hi_vel_axis,
                     av_image=av_data[~mask],
                     av_image_error=av_data_error[~mask],
                     vel_center=vel_center,
                     vel_widths=vel_widths,
                     dgrs=dgrs,
                     results_filename='',
                     return_likelihoods=True,
                     likelihood_filename=None,
                     clobber=False,
                     conf=conf,
                     )

    # Unpack output of likelihood calculation
    (vel_range_confint, width_confint, dgr_confint, likelihoods,
            width_likelihood, dgr_likelihood, width_max, dgr_max,
            vel_range_max) = results

    print('\nHI velocity integration range:')
    print('%.1f to %.1f km/s' % (vel_range_confint[0],
                                 vel_range_confint[1]))

    # Calulate chi^2 for best fit models
    # ----------------------------------
    nhi_image_temp, nhi_image_error = \
            calculate_nhi(cube=hi_data,
                velocity_axis=hi_vel_axis,
                velocity_range=vel_range_max,
                noise_cube=noise_cube,
                return_nhi_error=True)
    av_image_model = nhi_image_temp * dgr_max
    # avoid NaNs
    indices = ((av_image_model == av_image_model) & \
               (av_data == av_data))
    # add nan locations to the mask
    mask[~indices] = 1

    # count number of pixels used in analysis
    npix = mask[~mask].size

    # finally calculate chi^2
    chisq = np.sum((av_data[~mask] - av_image_model[~mask])**2 / \
            av_data_error[~mask]**2) / av_data[~mask].size

    print('\nTotal number of pixels in analysis, after masking = ' + \
            '{0:.0f}'.format(npix))

    print('\nReduced chi^2 = {0:.1f}'.format(chisq))

    # Write results to global properties
    global_props['dust2gas_ratio'] = {}
    global_props['dust2gas_ratio_error'] = {}
    global_props['hi_velocity_width'] = {}
    global_props['hi_velocity_width_error'] = {}
    global_props['dust2gas_ratio_max'] = {}
    global_props['hi_velocity_center_max'] = {}
    global_props['hi_velocity_width_max'] = {}
    global_props['hi_velocity_range_max'] =  {}
    global_props['av_threshold'] = {}
    global_props['co_threshold'] = {}
    global_props['hi_velocity_width']['value'] = width_confint[0]
    global_props['hi_velocity_width']['unit'] = 'km/s'
    global_props['hi_velocity_width_error']['value'] = width_confint[1:]
    global_props['hi_velocity_width_error']['unit'] = 'km/s'
    global_props['hi_velocity_range'] = vel_range_confint[0:2]
    global_props['hi_velocity_range_error'] = vel_range_confint[2:]
    global_props['dust2gas_ratio']['value'] = dgr_confint[0]
    global_props['dust2gas_ratio_error']['value'] = dgr_confint[1:]
    global_props['dust2gas_ratio_max']['value'] = dgr_max
    global_props['hi_velocity_center_max']['value'] = vel_center
    global_props['hi_velocity_width_max']['value'] = width_max
    global_props['hi_velocity_range_max']['value'] = vel_range_max
    global_props['hi_velocity_range_conf'] = conf
    global_props['width_likelihood'] = width_likelihood.tolist()
    global_props['dgr_likelihood'] = dgr_likelihood.tolist()
    global_props['vel_centers'] = [vel_center,]
    global_props['vel_widths'] = vel_widths.tolist()
    global_props['dgrs'] = dgrs.tolist()
    global_props['likelihoods'] = likelihoods.tolist()
    global_props['av_threshold']['value'] = None
    global_props['av_threshold']['unit'] = 'mag'
    global_props['co_threshold']['value'] = None
    global_props['co_threshold']['unit'] = 'K km/s'
    global_props['chisq'] = chisq
    global_props['npix'] = npix
    global_props['mask'] = mask.tolist()
    global_props['use_binned_image'] = True

    # Write the file
    print('\nWriting results to\n' + global_property_file + '_scaled.txt')

    with open(property_dir + global_property_file + '_scaled.txt', 'w') as f:
        json.dump(global_props, f)

    # Plot likelihood space
    print('\nWriting likelihood image to\n' + results_filename + \
            '_scaled_wd.png')
    plot_likelihoods_hist(global_props,
                          plot_axes=('widths', 'dgrs'),
                          show=0,
                          returnimage=False,
                          filename=results_filename + '_scaled_wd.png',
                          contour_confs=contour_confs)

if __name__ == '__main__':
    main()









