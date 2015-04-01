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
              #'figure.figsize': (8, 8 ),
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
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    font_scale = 9
    line_weight = 600
    font_weight = 600
    params = {
              'axes.color_cycle': color_cycle, # colors of different plots
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              #'axes.weight': line_weight,
              'axes.linewidth': 1.2,
              'axes.labelweight': font_weight,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': font_weight,
              'font.serif': 'computer modern roman',
              'text.fontsize': font_scale,
              'text.usetex': True,
              'text.latex.preamble': r'\usepackage[T1]{fontenc}',
              #'font.family': 'sans-serif',
              'figure.figsize': (3.6, 3.6),
              'figure.dpi': 600,
              'backend' : 'pdf',
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    pgf_with_pdflatex = {
        "pgf.texsystem": "pdflatex",
        "pgf.preamble": [
             r"\usepackage[utf8x]{inputenc}",
             r"\usepackage[T1]{fontenc}",
             r"\usepackage{cmbright}",
             ]
    }
    plt.rcParams.update(pgf_with_pdflatex)


    fig, ax_image = plt.subplots()

    if plot_axes[0] == 'widths':
    	x_grid = global_props['vel_widths']
    	x_confint = (global_props['hi_velocity_width']['value'],
    	             global_props['hi_velocity_width_error']['value'][0],
    	             global_props['hi_velocity_width_error']['value'][1],
    	             )
        x_extent = x_grid[0], x_grid[-1]
        ax_image.set_xlabel(r'Velocity Width [km/s]')
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
        ax_image.set_ylabel(r'DGR [10$^{-20}$ cm$^2$ mag]')
        y_sum_axes = 0
        x_pdf_label = r'DGR PDF'
        y_limits = (y_grid[0], y_grid[-1])

    # Create axes
    sum_axes = np.array((x_sum_axes, y_sum_axes))
    sum_axis = np.argmax(np.bincount(np.ravel(sum_axes)))

    # Mask NaNs
    likelihoods = np.sum(global_props['likelihoods'], axis=(2))
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
        ax_pdf_x = divider.append_axes("top", 0.6, pad=0.1, sharex=ax_image)
        ax_pdf_y  = divider.append_axes("right", 0.6, pad=0.1,
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
    params = {
              'axes.color_cycle': color_cycle, # colors of different plots
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              #'axes.weight': line_weight,
              'axes.linewidth': 1.2,
              'axes.labelweight': font_weight,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': font_weight,
              'font.serif': 'computer modern roman',
              'text.fontsize': font_scale,
              'text.usetex': True,
              'text.latex.preamble': r'\usepackage[T1]{fontenc}',
              #'font.family': 'sans-serif',
              'figure.figsize': (3.6, 3.6),
              'figure.dpi': 600,
              'backend' : 'pdf',
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    pgf_with_pdflatex = {
        "pgf.texsystem": "pdflatex",
        "pgf.preamble": [
             r"\usepackage[utf8x]{inputenc}",
             r"\usepackage[T1]{fontenc}",
             r"\usepackage{cmbright}",
             ]
    }
    plt.rcParams.update(pgf_with_pdflatex)


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

def plot_mask_residuals(residuals=None, x_fit=None, y_fit=None,
        residual_thres=None, filename=None, show=True, title=None):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from scipy.integrate import simps as integrate

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    font_scale = 9
    line_weight = 600
    font_weight = 600

    params = {
              'axes.color_cycle': color_cycle, # colors of different plots
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              #'axes.weight': line_weight,
              'axes.linewidth': 1.2,
              'axes.labelweight': font_weight,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': font_weight,
              'font.serif': 'computer modern roman',
              'text.fontsize': font_scale,
              'text.usetex': True,
              'text.latex.preamble': r'\usepackage[T1]{fontenc}',
              #'font.family': 'sans-serif',
              'figure.figsize': (3.6, 3.6),
              'figure.dpi': 600,
              'backend' : 'pdf',
              #'figure.titlesize': font_scale,
             }
    plt.rcParams.update(params)

    pgf_with_pdflatex = {
        "pgf.texsystem": "pdflatex",
        "pgf.preamble": [
             r"\usepackage[utf8x]{inputenc}",
             r"\usepackage[T1]{fontenc}",
             r"\usepackage{cmbright}",
             ]
    }
    plt.rcParams.update(pgf_with_pdflatex)



    # Create figure instance
    fig = plt.figure()

    ax = fig.add_subplot(111)

    counts, bin_edges = \
        np.histogram(np.ravel(residuals[~np.isnan(residuals)]),
                                 bins=1000,
                                 )

    bin_edges_ext = np.zeros(len(counts) + 1)
    counts_ext = np.zeros(len(counts) + 1)

    bin_edges_ext[0] = bin_edges[0] - (bin_edges[1] - bin_edges[0])
    bin_edges_ext[1:] = bin_edges[:-1]
    counts_ext[0] = 0
    counts_ext[1:] = counts

    # Normalize so area = 1
    #counts_ext /= np.nansum(counts_ext) * (bin_edges_ext[2] - bin_edges_ext[1])
    counts_ext = counts_ext / integrate(counts_ext, x=bin_edges_ext)
    y_fit /= np.max(y_fit)
    y_fit *= np.max(counts_ext)
    print('max counts', np.max(counts_ext))

    ax.plot(bin_edges_ext, counts_ext, drawstyle='steps-mid',
            linewidth=1.5)
    ax.plot(x_fit, y_fit,
            linewidth=2,
            alpha=0.6)
    ax.set_xlim([np.nanmin(bin_edges_ext) - \
                 np.abs(0.8 * np.nanmin(bin_edges_ext)),4])
    ax.set_ylim([-0.1, 1.1])
    ax.axvline(residual_thres,
               color='k',
               linestyle='--',
               linewidth=1.5)
    ax.set_xlabel(r'Residual $A_V$ [mag]')
    ax.set_ylabel('Normalized PDF')

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=600)
    if show:
        plt.show()

''' Calculations
'''

def calc_logL(model, data, data_error=None, weights=None):

    '''
    Calculates log likelihood

    http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node2.html

    '''

    import numpy as np

    if data_error is None:
        data_error = np.std(data)

    if weights is None:
        weights = 1.0

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

def get_residual_mask(residuals, resid_width_scale=3.0, plot_progress=False,
        results_filename=None):

    '''

    '''

    import numpy as np
    from scipy.optimize import curve_fit, minimize


    # Fit the rising portion of the residuals
    residuals_crop = residuals[(residuals < 0) & \
                               #(residuals > -1.5) & \
                                ~np.isnan(residuals)]

    counts, bin_edges = np.histogram(np.ravel(residuals_crop),
                                     bins=100,
                                     )

    p0=(2, np.nanmax(counts), 0)

    if 0:
        fit_params = curve_fit(gauss,
                               bin_edges[:-1],
                               counts,
                               p0=p0,
                               maxfev=1000000,
                               )[0]
    elif 1:
        from lmfit import minimize, Parameters

        # Set parameter limits and initial guesses
        params = Parameters()
        params.add('width',
                   value=p0[0],
                   min=0.1,
                   max=10,
                   )
        params.add('amp',
                   value=p0[1],
                   min=0,
                   max=2 * np.nanmax(counts),
                   )
        params.add('x0',
                   value=p0[2],
                   min=-4,
                   max=4,
                   )

        def norm(params, bin_edges, counts):
            width = params['width'].value
            amp = params['amp'].value
            x0 = params['x0'].value

            model = gauss(bin_edges, width, amp, x0)

            norm = np.sum((counts - model)**2)

            return norm

        # Perform the fit!
        result = minimize(norm,
                          params,
                          args=(bin_edges[:-1], counts),
                          method='lbfgsb')

        fit_params = (params['width'].value, params['amp'].value,
                params['x0'].value)
    else:
        bounds = ((0, 10), (0, 5 * np.nanmax(counts)), (-10, 10))
        fit_params = minimize(gauss,
                              counts,
                              method='L-BFGS-B',
                              bounds=bounds,)

    # Include only residuals within 3 sigma
    residual_thres = resid_width_scale * np.abs(fit_params[0]) + fit_params[2]
    mask = residuals > residual_thres

    import matplotlib.pyplot as plt
    plt.clf(); plt.close();
    x_fit = np.linspace(np.nanmin(residuals),
                        np.nanmax(residuals),
                        1000)

    y_fit = gauss(x_fit, *fit_params)
    plt.plot(bin_edges[:-1], counts)
    plt.plot(x_fit, y_fit)
    plt.savefig('/usr/users/ezbc/Desktop/residuals.png')

    if results_filename is not None:
        x_fit = np.linspace(np.nanmin(residuals),
                            np.nanmax(residuals),
                            1000)

        y_fit = gauss(x_fit, *fit_params)
        y_fit / np.nanmax(residuals)

        print('\nSaving residual mask PDF figure to\n' + results_filename)
        plot_mask_residuals(residuals=residuals,
                            x_fit=x_fit,
                            y_fit=y_fit,
                            residual_thres=residual_thres,
                            filename=results_filename,
                            show=plot_progress)
        plot_mask_residuals(residuals=residuals,
                            x_fit=x_fit,
                            y_fit=y_fit,
                            residual_thres=residual_thres,
                            filename=results_filename.replace('.pdf', '.png'),
                            show=plot_progress)

    return mask

def iterate_residual_masking(
                             nhi_image=None,
                             nhi_image_error=None,
                             av_data=None,
                             av_data_error=None,
                             init_mask=None,
                             vel_range=None,
                             dgrs=None,
                             intercepts=None,
                             threshold_delta_dgr=None,
                             resid_width_scale=3.0,
                             plot_progress=False,
                             results_filename=None,
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
            np.isnan(nhi_image) | \
            np.isnan(nhi_image_error) | \
            (nhi_image_error == 0))

    # Apply initial mask to exclude throughout process
    if init_mask is not None:
        mask += init_mask

    # solve for DGR using linear least squares
    print('\nBeginning iterative DGR calculations + masking...')

    # Iterate masking pixels which are correlated and rederiving a linear least
    # squares solution for the DGR
    # -------------------------------------------------------------------------
    use_intercept = True
    delta_dgr = 1e10
    dgr = 1e10
    iteration = 0
    while delta_dgr > threshold_delta_dgr:
        if 0:
            N = len(np.ravel(nhi_image[~mask]))
            if use_intercept:
                A = np.array((np.ones(N),
                              np.ravel(nhi_image[~mask] / \
                                       nhi_image_error[~mask]),))
            else:
                A = np.array((np.ravel(nhi_image[~mask] / \
                              nhi_image_error[~mask]),))
            b = np.array((np.ravel(av_data[~mask] / av_data_error[~mask]),))
            A = np.matrix(A).T
            b = np.matrix(b).T

            a = (np.linalg.pinv(A) * b)
            if use_intercept:
                intercept = a[0, 0]
                dgr_new = a[1, 0]
            else:
                dgr_new = a[0, 0]
                intercept = 0
        else:
            results = calc_likelihoods(
                             nhi_image=nhi_image[~mask],
                             av_image=av_data[~mask],
                             av_image_error=av_data_error[~mask],
                             #image_weights=bin_weights[~mask],
                             #vel_center=vel_center_masked,
                             vel_widths=np.arange(0,1,1),
                             dgrs=dgrs,
                             intercepts=intercepts,
                             results_filename='',
                             return_likelihoods=True,
                             likelihood_filename=None,
                             clobber=False,
                             verbose=False
                             )

        # Unpack output of likelihood calculation
        (vel_range_confint, width_confint, dgr_confint, intercepts_confint,
                likelihoods, width_likelihood, dgr_likelihood,
                intercept_likelihood, width_max, dgr_max, intercept_max,
                vel_range_max) = results

        dgr_new = dgr_max
        intercept = intercept_max

        # Create model with the DGR
        if verbose:
            print('Iteration {0:.0f} results:'.format(iteration))
            print('\tDGR = {0:.2} 10^20 cm^2 mag'.format(dgr_new))
            print('\tIntercept = {0:.2f} mag'.format(intercept))
            print('')

        av_image_model = nhi_image * dgr_new + intercept

        #if dgr == 1e10:
        #    residuals = av_data - av_image_model
        #else:
        #    residuals = av_data - av_image_model + intercept
        residuals = av_data - av_image_model

        residuals[mask] = np.nan

        if 0:
            import matplotlib.pyplot as plt
            plt.imshow(residuals, origin='lower')
            plt.colorbar(cmap=plt.cm.gnuplot)
            plt.show()

        # Include only residuals which are white noise
        if iteration == 0:
            plot_filename = results_filename
        else:
            plot_filename = None
        mask_new = get_residual_mask(residuals,
                                     resid_width_scale=resid_width_scale,
                                     plot_progress=plot_progress,
                                     results_filename=plot_filename)

        # Mask non-white noise, i.e. correlated residuals.
        mask[mask_new] = 1

        if verbose:
            npix = mask.size - np.sum(mask)
            print('\tNumber of non-masked pixels = {0:.0f}'.format(npix))

        # Reset while loop conditions
        delta_dgr = np.abs(dgr - dgr_new)
        dgr = dgr_new
        iteration += 1

    # Plot results
    if 0:
        mask_new = get_residual_mask(residuals,
                                     resid_width_scale=resid_width_scale,
                                     plot_progress=plot_progress,
                                     results_filename=results_filename)

    # Create model of Av
    av_model = dgr * nhi_image
    av_model[mask] = np.nan

    return (av_model, mask, dgr)

def calc_likelihoods(
        hi_cube=None,
        hi_vel_axis=None,
        nhi_image=None,
        av_image=None,
        av_image_error=None,
        image_weights=None,
        vel_center=None,
        vel_widths=None,
        dgrs=None,
        intercepts=None,
        plot_results=False,
        results_filename='',
        return_likelihoods=True,
        likelihood_filename=None,
        clobber=False,
        conf=0.68,
        threshold_delta_dgr=0.0005,
        verbose=False,
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
    from mystats import calc_symmetric_error

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
                                len(dgrs),
                                len(intercepts)))

        # Progress bar parameters
        total = float(likelihoods.size)
        count = 0

        for j, vel_width in enumerate(vel_widths):
            # Construct N(HI) image outside of DGR loop, then apply
            # DGRs in loop

            # use the hi cube and vel range if no nhi image provided
            if nhi_image is None:
                vel_range = np.array((vel_center - vel_width / 2.,
                                      vel_center + vel_width / 2.))
                nhi_image = calculate_nhi(cube=hi_cube,
                                          velocity_axis=hi_vel_axis,
                                          velocity_range=vel_range,
                                          return_nhi_error=False)

            # Cycle through DGR to estimate error
            for k, dgr in enumerate(dgrs):
                for m, intercept in enumerate(intercepts):
                    # Create model of Av with N(HI) and DGR
                    av_image_model = nhi_image * dgr + intercept

                    logL = calc_logL(av_image_model,
                                     av_image,
                                     data_error=av_image_error,
                                     weights=image_weights)

                    likelihoods[j, k, m] = logL
                    #print 'logL =', logL

                    # Shows progress each 10%
                    count += 1
                    abs_step = int((total * 1)/10) or 10
                    if count and not count % abs_step:
                        print "\t{0:.0%} processed".format(count/total)

            nhi_image = None

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
    intercept_likelihood = np.sum(likelihoods, axis=(0, 1)) / \
                                  np.sum(likelihoods)
    width_likelihood = np.sum(likelihoods, axis=(1, 2)) / \
            np.sum(likelihoods)
    dgr_likelihood = np.sum(likelihoods, axis=(0, 2)) / \
            np.sum(likelihoods)

    # Derive confidence intervals of parameters
    width_confint = calc_symmetric_error(vel_widths,
                                   width_likelihood,
                                   alpha=1.0 - conf)
    dgr_confint = calc_symmetric_error(dgrs,
                                 dgr_likelihood,
                                 alpha=1.0 - conf)
    intercept_confint = calc_symmetric_error(intercepts,
                                 intercept_likelihood,
                                 alpha=1.0 - conf)

    # Get values of best-fit model parameters
    max_loc = np.where(likelihoods == np.max(likelihoods))
    width_max = vel_widths[max_loc[0][0]]
    dgr_max = dgrs[max_loc[1][0]]
    intercept_max = intercepts[max_loc[2][0]]

    if verbose:
        print('\nVelocity widths = ' + \
                '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(width_confint[0],
                                                    width_confint[2],
                                                    np.abs(width_confint[1])))
        print('\nDGRs = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} 10^20 cm^2 mag'.format(dgr_confint[0],
                                                    dgr_confint[2],
                                                    np.abs(dgr_confint[1])))
        print('\nIntercepts = ' + \
        '{0:.2f} +{1:.2f}/-{2:.2f} 10^20 cm^2 mag'.format(intercept_confint[0],
                                                intercept_confint[2],
                                                np.abs(intercept_confint[1])))

    # Write PDF
    if vel_center is None:
        vel_center = 0.0

    upper_lim = (np.nanmean(vel_center) + width_confint[0]/2.)
    lower_lim = (np.nanmean(vel_center) - width_confint[0]/2.)
    upper_lim_error = width_confint[2]**2
    lower_lim_error = width_confint[1]**2

    vel_range_confint = (lower_lim, upper_lim, lower_lim_error,
                         upper_lim_error)
    vel_range_max = (vel_center - width_max/2.0, vel_center + width_max/2.0)

    if not return_likelihoods:
        return vel_range_confint, dgr_confint
    else:
        return (vel_range_confint, width_confint, dgr_confint,
                intercept_confint, likelihoods,
                width_likelihood, dgr_likelihood,
                intercept_likelihood, width_max, dgr_max, intercept_max,
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

        if coord in ('region_limit',
                     'plot_limit',
                     'region_limit_bin',
                     'plot_limit_bin'):
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

def load_ds9_region(props, filename=None, header=None):

    import pyregion as pyr

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    props['regions'] = {}


    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        region_name = tag[tag.find('text={')+6:tag.find('}')].lower()

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

        props['regions'][region_name] = {}
        props['regions'][region_name]['poly_verts'] = {}
        props['regions'][region_name]['poly_verts']['wcs'] = poly_verts
        props['regions'][region_name]['poly_verts']['pixel'] = poly_verts_pix

    return props

def run_likelihood_analysis(av_data_type='planck', region=None,
        vel_range=None, resid_width_scale=3.0, single_vel_center=True):

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
    noise_cube_filename = 'california_hi_galfa_cube_regrid_planckres_noise'

    # Threshold for converging DGR
    threshold_delta_dgr = 0.00005

    # Name of property files results are written to
    global_property_file = 'california_global_properties'

    # Likelihood axis resolutions
    vel_widths = np.arange(1, 75, 2*0.16667)
    dgrs = np.arange(0.001, 0.8, 1e-3)
    intercepts = np.arange(0, 1, 1)

    # Velocity range over which to integrate HI for deriving the mask
    if vel_range is None:
        vel_range = (2.2,7.6)
        vel_range = (1.7,7.7)
        vel_range = (-5, 15)

    # Bin width in degrees
    bin_width_deg = 1.0

    # Clobber the binned images and remake them?
    clobber_bin_images = True

    # Filetype extensions for figures
    figure_types = ('png', 'pdf')

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'
    figure_dir = \
        '/d/bip3/ezbc/california/figures/'
    av_dir = '/d/bip3/ezbc/california/data/av/'
    hi_dir = '/d/bip3/ezbc/california/data/hi/'
    co_dir = '/d/bip3/ezbc/california/data/co/'
    core_dir = '/d/bip3/ezbc/california/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/california/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    likelihood_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'

    # Load data
    # ---------
    # Adjust filenames
    #noise_cube_filename += bin_string
    likelihood_filename = 'california_likelihood_{0:s}_bin'.format(av_data_type)
    results_filename = 'california_likelihood_{0:s}_bin'.format(av_data_type)
    # load Planck Av and GALFA HI images, on same grid
    if av_data_type == 'k09':
        print('\nLoading K+09 2MASS data...')
        av_data, av_header = fits.getdata(av_dir + \
                                  'california_av_k09_regrid_planckres.fits',
                                  header=True)
        av_data_error = 0.1 * np.ones(av_data.shape)
    else:
    	print('\nLoading Planck data...')
        av_data, av_header = fits.getdata(av_dir + \
                                          'california_av_planck_5arcmin.fits',
                                          header=True)

        av_data_error, av_error_header = fits.getdata(av_dir + \
                                    'california_av_error_planck_5arcmin.fits',
                                    header=True)

    hi_data, hi_header = fits.getdata(hi_dir + \
                'california_hi_galfa_cube_regrid_planckres.fits',
            header=True)

    # Load global properties
    with open(property_dir + global_property_file + '.txt', 'r') as f:
        global_props = json.load(f)

    # Prepare data products
    # ---------------------
    # Name correct region of cloud
    if region == 1:
        region_name = 'california1'
    elif region == 2:
        region_name = 'california2'
    else:
        region_name = 'california'
    global_property_file = global_property_file.replace('california', region_name)

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

    # Load cloud division regions from ds9
    global_props = load_ds9_region(global_props,
                            filename=region_dir + 'multicloud_divisions.reg',
                            header=av_header)

    # Derive relevant region
    region_vertices = \
        global_props['regions'][region_name]['poly_verts']['pixel']

    # block off region
    region_mask = np.logical_not(myg.get_polygon_mask(av_data,
                                                      region_vertices))

    if 0:
        import matplotlib.pyplot as plt
        plt.imshow(np.ma.array(av_data, mask=region_mask), origin='lower')
        plt.colorbar()
        plt.show()

    print('\nRegion size = ' + \
          '{0:.0f} pix'.format(region_mask[region_mask == 1].size))

    # Derive mask by excluding correlated residuals
    # ---------------------------------------------
    nhi_image, nhi_image_error = calculate_nhi(cube=hi_data,
                              velocity_axis=hi_vel_axis,
                              velocity_range=vel_range,
                              noise_cube=noise_cube,
                              velocity_noise_range=[90, 110],
                              Tsys=30.0,
                              return_nhi_error=True,
                              )
    if 0:
        vel_center = np.zeros(hi_data.shape[1:])
        for i in xrange(0, hi_data.shape[1]):
            for j in xrange(0, hi_data.shape[2]):
                hi_spectrum = hi_data[:, i, j]
                hi_spectrum[np.isnan(hi_spectrum)] = 0.0
                if np.nansum(hi_spectrum) > 0:
                    vel_center[i,j] = \
                            np.array((np.average(hi_vel_axis,
                                                 weights=hi_spectrum**2),))[0]
                else:
                    vel_center[i,j] = np.nan
        vel_range = (vel_center - 5, vel_center + 5)
        nhi_image1, nhi_image_error = calculate_nhi(cube=hi_data,
                                  velocity_axis=hi_vel_axis,
                                  velocity_range=vel_range,
                                  noise_cube=noise_cube,
                                  velocity_noise_range=[90, 110],
                                  Tsys=30.0,
                                  return_nhi_error=True,
                                  )
        hi_spectrum = np.sum(hi_data, axis=(1,2))
        vel_center = np.array((np.average(hi_vel_axis,
                               weights=hi_spectrum**2),))[0]
        vel_range = (vel_center - 5, vel_center + 5)
        nhi_image2, nhi_image_error = calculate_nhi(cube=hi_data,
                                  velocity_axis=hi_vel_axis,
                                  velocity_range=vel_range,
                                  noise_cube=noise_cube,
                                  velocity_noise_range=[90, 110],
                                  Tsys=30.0,
                                  return_nhi_error=True,
                                  )
        import matplotlib.pyplot as plt
        plt.close(); plt.clf()
        plt.imshow(nhi_image1 - nhi_image2, origin='lower left')
        plt.colorbar()
        plt.show()

    print('\nDeriving mask for correlated residuals...')

    av_model, mask, dgr = iterate_residual_masking(
                             nhi_image=nhi_image,
                             nhi_image_error=nhi_image_error,
                             av_data=av_data,
                             av_data_error=av_data_error,
                             vel_range=vel_range,
                             dgrs=dgrs,
                             intercepts=intercepts,
                             threshold_delta_dgr=threshold_delta_dgr,
                             resid_width_scale=resid_width_scale,
                             init_mask=region_mask,
                             verbose=1,
                             plot_progress=0,
                             results_filename=figure_dir + 'likelihood/'\
                                              'california_residual_pdf.pdf'
                             )

    # Combine region mask with new mask
    #mask += np.logical_not(region_mask)
    mask += region_mask
    mask = mask.astype('bool')

    # Write full resolution mask to parameters
    global_props['mask'] = mask.tolist()

    if 1:
        import matplotlib.pyplot as plt
        plt.imshow(np.ma.array(av_data, mask=mask), origin='lower')
        plt.show()

    # Bin the masked images
    # ---------------------
    print('\nBinning masked images...')

    # Mask the data with nans
    av_data[mask] = np.nan
    av_data_error[mask] = np.nan
    hi_data[:, mask] = np.nan

    if not check_file(av_dir + 'california_av_planck_5arcmin_masked.fits',
                      clobber=clobber_bin_images):
        fits.writeto(av_dir + 'california_av_planck_5arcmin_masked.fits',
                     av_data,
                     av_header)

    if not check_file(av_dir + 'california_av_error_planck_5arcmin_masked.fits',
                      clobber=clobber_bin_images):
        fits.writeto(av_dir + 'california_av_error_planck_5arcmin_masked.fits',
                     av_data_error,
                     av_header)

    if not check_file(hi_dir + \
                      'california_hi_galfa_cube_regrid_planckres_masked.fits',
                      clobber=clobber_bin_images):
        fits.writeto(hi_dir + \
                     'california_hi_galfa_cube_regrid_planckres_masked.fits',
                     hi_data,
                     hi_header)

    if clobber_bin_images:
        # Define number of pixels in each bin
        binsize = bin_width_deg * 60.0 / 5.0

        # Bin the images, retain only one bin_weight image since they are all
        # the same
        # -------------------------------------------------------------------
        # Av image
        av_data_bin, av_header_bin, bin_weights = \
                bin_image(av_data,
                          binsize=(binsize, binsize),
                          header=av_header,
                          func=np.nanmean,
                          return_weights=True)

        if not check_file(av_dir + 'california_av_planck_5arcmin_bin.fits',
                          clobber=clobber_bin_images):
            fits.writeto(av_dir + 'california_av_planck_5arcmin_bin.fits',
                         av_data_bin,
                         av_header_bin)
        if not check_file(av_dir + 'california_av_planck_5arcmin_bin_weights.fits',
                          clobber=clobber_bin_images):
            fits.writeto(av_dir + 'california_av_planck_5arcmin_bin_weights.fits',
                         bin_weights,
                         av_header_bin)

        # Av image error
        # Errors add in square
        # mean = sum(a_i) / n
        # error on mean = sqrt(sum(a_i**2 / n**2))
        noise_func = lambda x: np.nansum(x**2)**0.5 / x[~np.isnan(x)].size

        av_data_error_bin, av_header_bin = \
                bin_image(av_data_error,
                          binsize=(binsize, binsize),
                          header=av_header,
                          func=noise_func,)

        if not check_file(av_dir + 'california_av_error_planck_5arcmin_bin.fits',
                          clobber=clobber_bin_images):
            fits.writeto(av_dir + 'california_av_error_planck_5arcmin_bin.fits',
                         av_data_error_bin,
                         av_header_bin)

        # Hi image
        hi_data_bin, hi_header_bin = \
                bin_image(hi_data,
                          binsize=(binsize, binsize),
                          header=hi_header,
                          func=np.nanmean)

        if not check_file(hi_dir + \
                          'california_hi_galfa_cube_regrid_planckres_bin.fits',
                          clobber=clobber_bin_images):
            fits.writeto(hi_dir + \
                         'california_hi_galfa_cube_regrid_planckres_bin.fits',
                         hi_data_bin,
                         hi_header_bin)

    # Load data
    # ---------
    bin_string = '_bin'

    # Adjust filenames
    noise_cube_filename += bin_string
    likelihood_filename = '{0:s}_likelihood_{1:s}_bin'.format(region_name,
                                                              av_data_type)
    results_filename = '{0:s}_likelihood_{1:s}_bin'.format(region_name,
                                                           av_data_type)

    av_data, av_header = fits.getdata(av_dir + \
                            'california_av_planck_5arcmin' + bin_string + '.fits',
                                      header=True)

    av_data_error, av_error_header = fits.getdata(av_dir + \
                'california_av_error_planck_5arcmin' + bin_string + '.fits',
            header=True)

    bin_weights = fits.getdata(av_dir + \
                               'california_av_planck_5arcmin' + bin_string + \
                               '_weights.fits',)

    hi_data, hi_header = fits.getdata(hi_dir + \
                            'california_hi_galfa_cube_regrid_planckres' + \
                            bin_string + '.fits',
                            header=True)

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

    # Prepare data products
    # ---------------------
    # Change WCS coords to pixel coords of images
    global_props['region_limit_bin'] = global_props['region_limit'].copy()
    global_props['plot_limit_bin'] = global_props['plot_limit'].copy()
    global_props = convert_limit_coordinates(global_props,
                                             header=av_header,
                                             coords=('region_limit_bin',
                                                     'plot_limit_bin'))

    # Derive relevant region
    pix = global_props['region_limit_bin']['pixel']
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
    if single_vel_center:
        hi_spectrum = np.sum(hi_data[:, ~mask], axis=(1))
        vel_center = np.array((np.average(hi_vel_axis,
                               weights=hi_spectrum**2),))[0]
        print('\nVelocity center from HI = ' +\
                '{0:.2f} km/s'.format(vel_center))
        vel_center_masked = vel_center
    else:
        vel_center = np.zeros(hi_data.shape[1:])
        for i in xrange(0, hi_data.shape[1]):
            for j in xrange(0, hi_data.shape[2]):
                hi_spectrum = hi_data[:, i, j]
                hi_spectrum[np.isnan(hi_spectrum)] = 0.0
                if np.nansum(hi_spectrum) > 0:
                    vel_center[i,j] = \
                            np.array((np.average(hi_vel_axis,
                                                 weights=hi_spectrum**2),))[0]
                else:
                    vel_center[i,j] = np.nan

        vel_center_masked = vel_center[~mask]

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
                     #image_weights=bin_weights[~mask],
                     vel_center=vel_center_masked,
                     vel_widths=vel_widths,
                     dgrs=dgrs,
                     intercepts=intercepts,
                     results_filename='',
                     return_likelihoods=True,
                     likelihood_filename=None,
                     clobber=False,
                     conf=conf,
                     )

    # Unpack output of likelihood calculation
    (vel_range_confint, width_confint, dgr_confint, intercepts_confint,
            likelihoods, width_likelihood, dgr_likelihood,
            intercept_likelihood, width_max, dgr_max, intercept_max,
            vel_range_max) = results

    print('\nHI velocity integration range:')
    print('%.1f to %.1f km/s' % (vel_range_confint[0],
                                 vel_range_confint[1]))

    vel_range_max = (vel_center - width_max / 2.0,
                     vel_center + width_max / 2.0)

    # Calulate chi^2 for best fit models
    # ----------------------------------
    nhi_image_temp = calculate_nhi(cube=hi_data,
                                   velocity_axis=hi_vel_axis,
                                   velocity_range=vel_range_max,
                                   noise_cube=noise_cube,
                                   return_nhi_error=False)

    av_image_model = nhi_image_temp * dgr_max + intercept_max

    # count number of pixels used in analysis
    npix = mask[~mask].size

    # finally calculate chi^2
    #chisq = np.sum((av_data[~mask] - av_image_model[~mask])**2 / \
    #        av_data_error[~mask]**2) / av_data[~mask].size

    print('\nTotal number of pixels in analysis, after masking = ' + \
            '{0:.0f}'.format(npix))

    #print('\nReduced chi^2 = {0:.1f}'.format(chisq))

    # Write results to global properties
    global_props['dust2gas_ratio'] = {}
    global_props['dust2gas_ratio_error'] = {}
    global_props['intercept'] = {}
    global_props['intercept_error'] = {}
    global_props['hi_velocity_width'] = {}
    global_props['hi_velocity_width_error'] = {}
    global_props['dust2gas_ratio_max'] = {}
    global_props['intercept_max'] = {}
    global_props['hi_velocity_center'] = {}
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
    global_props['intercept_max']['value'] = intercept_max
    global_props['intercept']['value'] = intercepts_confint[0]
    global_props['intercept_error']['value'] = intercepts_confint[1:]
    global_props['hi_velocity_center']['value'] = vel_center.tolist()
    #global_props['hi_velocity_width_max']['value'] = width_max
    #global_props['hi_velocity_range_max']['value'] = vel_range_max
    global_props['hi_velocity_range_conf'] = conf
    global_props['width_likelihood'] = width_likelihood.tolist()
    global_props['dgr_likelihood'] = dgr_likelihood.tolist()
    global_props['vel_centers'] = vel_center.tolist()
    global_props['vel_widths'] = vel_widths.tolist()
    global_props['dgrs'] = dgrs.tolist()
    global_props['likelihoods'] = likelihoods.tolist()
    global_props['av_threshold']['value'] = None
    global_props['av_threshold']['unit'] = 'mag'
    global_props['co_threshold']['value'] = None
    global_props['co_threshold']['unit'] = 'K km/s'
    #global_props['chisq'] = chisq
    global_props['npix'] = npix
    global_props['mask_bin'] = mask.tolist()
    global_props['use_binned_image'] = True
    global_props['residual_width_scale'] = resid_width_scale
    global_props['threshold_delta_dgr'] = threshold_delta_dgr

    # Write the file
    print('\nWriting results to\n' + global_property_file + '_' + \
          av_data_type + '_init.txt')

    with open(property_dir + global_property_file + '_' + av_data_type + \
              '_init.txt', 'w') as f:
        json.dump(global_props, f, allow_nan=True)

    # Plot likelihood space
    for figure_type in figure_types:
        print('\nWriting likelihood image to\n' + results_filename + \
              '_init_wd.{0:s}'.format(figure_type))
        plot_likelihoods_hist(global_props,
                              plot_axes=('widths', 'dgrs'),
                              show=0,
                              returnimage=False,
                              filename=results_filename + \
                                       '_init_wd.{0:s}'.format(figure_type),
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
    #std = np.sqrt(np.sum((av_data - av_image_model)**2 \
    #                     / (av_data.size - 2)))

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
                     image_weights=bin_weights[~mask],
                     vel_center=vel_center_masked,
                     vel_widths=vel_widths,
                     dgrs=dgrs,
                     intercepts=intercepts,
                     results_filename='',
                     return_likelihoods=True,
                     likelihood_filename=None,
                     clobber=False,
                     conf=conf,
                     )

    # Unpack output of likelihood calculation
    (vel_range_confint, width_confint, dgr_confint, intercept_confint,
            likelihoods, width_likelihood, dgr_likelihood,
            intercept_likelihood, width_max, dgr_max, intercept_max,
            vel_range_max) = results


    print('\nHI velocity integration range:')
    print('%.1f to %.1f km/s' % (vel_range_confint[0],
                                 vel_range_confint[1]))

    # Calulate chi^2 for best fit models
    # ----------------------------------
    nhi_image_temp = \
            calculate_nhi(cube=hi_data,
                velocity_axis=hi_vel_axis,
                velocity_range=vel_range_max,
                noise_cube=noise_cube)

    av_image_model = nhi_image_temp * dgr_max
    # avoid NaNs
    indices = ((av_image_model == av_image_model) & \
               (av_data == av_data))
    # add nan locations to the mask
    mask[~indices] = 1

    # count number of pixels used in analysis
    npix = mask[~mask].size

    # finally calculate chi^2
    #chisq = np.sum((av_data[~mask] - av_image_model[~mask])**2 / \
    #        av_data_error[~mask]**2) / av_data[~mask].size

    print('\nTotal number of pixels in analysis, after masking = ' + \
            '{0:.0f}'.format(npix))

    #print('\nReduced chi^2 = {0:.1f}'.format(chisq))

    # Write results to global properties
    global_props['dust2gas_ratio'] = {}
    global_props['dust2gas_ratio_error'] = {}
    global_props['intercept'] = {}
    global_props['intercept_error'] = {}
    global_props['hi_velocity_width'] = {}
    global_props['hi_velocity_width_error'] = {}
    global_props['dust2gas_ratio_max'] = {}
    global_props['intercept_max'] = {}
    global_props['hi_velocity_center'] = {}
    global_props['hi_velocity_width_max'] = {}
    global_props['hi_velocity_range_max'] =  {}
    global_props['av_threshold'] = {}
    global_props['co_threshold'] = {}
    global_props['hi_velocity_width']['value'] = width_confint[0]
    global_props['hi_velocity_width']['unit'] = 'km/s'
    global_props['hi_velocity_width_max']['value'] = width_max
    global_props['hi_velocity_width_max']['unit'] = 'km/s'
    global_props['hi_velocity_width_error']['value'] = width_confint[1:]
    global_props['hi_velocity_width_error']['unit'] = 'km/s'
    global_props['hi_velocity_range'] = vel_range_confint[0:2]
    global_props['hi_velocity_range_error'] = vel_range_confint[2:]
    global_props['dust2gas_ratio']['value'] = dgr_confint[0]
    global_props['dust2gas_ratio_error']['value'] = dgr_confint[1:]
    global_props['dust2gas_ratio_max']['value'] = dgr_max
    global_props['intercept']['value'] = intercept_confint[0]
    global_props['intercept_error']['value'] = intercept_confint[1:]
    global_props['intercept_max']['value'] = intercept_max
    global_props['hi_velocity_center']['value'] = vel_center.tolist()
    #global_props['hi_velocity_width_max']['value'] = width_max
    #global_props['hi_velocity_range_max']['value'] = vel_range_max
    global_props['hi_velocity_range_conf'] = conf
    global_props['width_likelihood'] = width_likelihood.tolist()
    global_props['dgr_likelihood'] = dgr_likelihood.tolist()
    global_props['intercept_likelihood'] = intercept_likelihood.tolist()
    global_props['vel_centers'] = vel_center.tolist()
    global_props['vel_widths'] = vel_widths.tolist()
    global_props['dgrs'] = dgrs.tolist()
    global_props['intercepts'] = intercepts.tolist()
    global_props['likelihoods'] = likelihoods.tolist()
    global_props['av_threshold']['value'] = None
    global_props['av_threshold']['unit'] = 'mag'
    global_props['co_threshold']['value'] = None
    global_props['co_threshold']['unit'] = 'K km/s'
    #global_props['chisq'] = chisq
    global_props['npix'] = npix
    global_props['mask_bin'] = mask.tolist()
    global_props['use_binned_image'] = True
    global_props['residual_width_scale'] = resid_width_scale
    global_props['threshold_delta_dgr'] = threshold_delta_dgr

    # Write the file
    print('\nWriting results to\n' + global_property_file + '_' + \
          av_data_type + '_scaled.txt')

    with open(property_dir + global_property_file + '_' + av_data_type + \
              '_scaled.txt', 'w') as f:
        json.dump(global_props, f)

    # Plot likelihood space
    for figure_type in figure_types:
        print('\nWriting likelihood image to\n' + results_filename + \
              '_scaled_wd.{0:s}'.format(figure_type))
        plot_likelihoods_hist(global_props,
                              plot_axes=('widths', 'dgrs'),
                              show=0,
                              returnimage=False,
                              filename=results_filename + \
                                       '_scaled_wd.{0:s}'.format(figure_type),
                              contour_confs=contour_confs)

    #return global_props['hi_velocity_range']
    return global_props

'''
Main Script
'''

def main():

    import numpy as np
    from os import path
    import json
    from pandas import DataFrame

    av_data_type = 'planck'

    # threshold in velocity range difference
    vel_range_diff_thres = 3.0 # km/s

    property_dir = \
        '/d/bip3/ezbc/california/data/python_output/residual_parameter_results/'

    final_property_dir = '/d/bip3/ezbc/california/data/python_output/'

    property_filename = 'california_global_properties_planck'

    # Number of white noise standard deviations with which to fit the
    # residuals in iterative masking
    residual_width_scales = [3.0,]

    regions = [None, ]

    clobber_results = True

    # Use single velocity center for entire image?
    single_vel_center = False

    table_cols = ('dust2gas_ratio', 'hi_velocity_width',
                  'hi_velocity_width', 'intercept', 'residual_width_scale')
    n = len(residual_width_scales)
    table_df = DataFrame({col:np.empty(n) for col in table_cols})

    for region in regions:
        # Grab correct region
        if region == 1:
            region_name = 'california1'
        elif region == 2:
            region_name = 'california2'
        else:
            region_name = 'california'

        property_filename = 'california_global_properties_planck'
        property_filename = property_filename.replace('california', region_name)

        print('\nPerforming likelihood derivations for ' + region_name)

        for i, residual_width_scale in enumerate(residual_width_scales):
            iteration = 0
            vel_range = (-20.0, 30.0)
            vel_range_new = (-1.0, 1.0)
            vel_range_diff = np.sum(np.abs(np.array(vel_range) - \
                                           np.array(vel_range_new)))

            while vel_range_diff > vel_range_diff_thres:
                json_filename = property_dir + property_filename + '_' + \
                            av_data_type + \
                            '_residscale{0:.1f}'.format(residual_width_scale)\
                            + '_iter{0:.0f}'.format(iteration) + \
                            '_centervary.txt'

                exists = path.isfile(json_filename)

                print('Writing iteration data file to ' + json_filename)

                if exists and not clobber_results:
                    with open(json_filename, 'r') as f:
                        global_props = json.load(f)
                else:
                    global_props = \
                        run_likelihood_analysis(av_data_type=av_data_type,
                                        vel_range=vel_range,
                                        region=region,
                                        single_vel_center=single_vel_center,
                                        resid_width_scale=residual_width_scale)

                vel_range_new = global_props['hi_velocity_range']

                vel_range_diff = np.sum(np.abs(np.array(vel_range) - \
                                               np.array(vel_range_new)))

                if clobber_results:
                    with open(json_filename, 'w') as f:
                        json.dump(global_props, f)

                print('\n\n\n Next iteration \n-------------------\n\n\n')
                print('Velocity range difference =' + \
                      ' {0:.1f}'.format(vel_range_diff))

                vel_range = vel_range_new

                iteration += 1

            # Write important results to table
            for col in table_df:
                if col == 'residual_width_scale':
                    table_df[col][i] = global_props[col]
                else:
                    table_df[col][i] = global_props[col]['value']

            # Write the file
            print('\nWriting results to\n' + property_filename + \
                    '_' + av_data_type + '_scaled_centervary.txt')

            with open(final_property_dir + property_filename +\
                    '_' + av_data_type + '_scaled_centervary.txt', 'w') as f:
                json.dump(global_props, f)

if __name__ == '__main__':
    main()









