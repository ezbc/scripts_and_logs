#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the perseus molecular cloud.
'''

import pyfits as pf
import numpy as np


''' Plotting Functions
'''
def plot_correlations(correlations,velocity_centers,velocity_widths,
        filename=None,show=True, returnimage=False):
    ''' Plots a heat map of correlation values as a function of velocity width
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

    # Unravel the correlations if raveled
    if len(correlations.shape) == 1:
        correlations_image = np.empty((velocity_centers.shape[0],
                                       velocity_widths.shape[0]))
        correlations_image[:,:] = np.NaN
        count = 0
        try:
            for i, center in enumerate(velocity_centers):
                for j, width in enumerate(velocity_widths):
                    correlations_image[i,j] = correlations[count]
                    count += 1
        except IndexError:
            print(' plot_correlations: O-d array input, cannot proceed')
    else:
       	correlations_image = correlations

    image = np.ma.array(correlations_image, mask=np.isnan(correlations_image))

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
        return correlations_image

def plot_correlations_hist(correlations, velocity_centers, velocity_widths,
        center_pdf=None, width_pdf=None, center_confint=None,
        width_confint=None, filename=None, show=True, returnimage=False):

    ''' Plots a heat map of correlation values as a function of velocity width
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

    fig, ax_image = plt.subplots(figsize=(8,8))

    # Unravel the correlations if raveled
    if len(correlations.shape) == 1:
        correlations_image = np.empty((velocity_centers.shape[0],
                                       velocity_widths.shape[0]))
        correlations_image[:,:] = np.NaN
        count = 0
        try:
            for i, center in enumerate(velocity_centers):
                for j, width in enumerate(velocity_widths):
                    correlations_image[i,j] = correlations[count]
                    count += 1
        except IndexError:
            print(' plot_correlations: O-d array input, cannot proceed')
    else:
    	correlations_image = correlations

    # Mask NaNs
    image = np.ma.array(correlations_image, mask=np.isnan(correlations_image))

    #plt.rc('text', usetex=False)
    im = ax_image.imshow(image, interpolation='nearest', origin='lower',
            extent=[velocity_widths[0], velocity_widths[-1],
                    velocity_centers[0], velocity_centers[-1]],
            #cmap=plt.cm.gist_stern,
            cmap=plt.cm.gray,
            norm=matplotlib.colors.LogNorm(),
            )

    ax_image.set_aspect(1.)

    ax_image.set_xlabel(r'Velocity Width (km/s)')
    ax_image.set_ylabel(r'Velocity Center (km/s)')

    # Ticks only every 5 km/s
    #ax_image.set_xticks(np.arange(0, velocity_widths.shape[0], 1),#[::],
    #                    velocity_centers)#[::])
    #ax_image.set_yticks(np.arange(0, velocity_centers.shape[0], 1)[::5],
    #                    velocity_centers[::5])

    show_pdfs = 1

    if show_pdfs:
        divider = make_axes_locatable(ax_image)
        ax_pdf_width = divider.append_axes("top", 1, pad=0.1, sharex=ax_image)
        ax_pdf_center = divider.append_axes("right", 1, pad=0.1,
                sharey=ax_image)

        # make some labels invisible
        plt.setp(ax_pdf_width.get_xticklabels() + \
                 ax_pdf_center.get_yticklabels(),
                 visible=False)

        ax_pdf_width.plot(velocity_widths,
                          width_pdf,
                          color='k',
                          drawstyle='steps-pre',
                          linewidth=2,
                          )
        ax_pdf_center.plot(center_pdf,
                           velocity_centers,
                           color='k',
                           drawstyle='steps-pre',
                           linewidth=2,
                           )

        #axHistx.axis["bottom"].major_ticklabels.set_visible(False)
        for tl in ax_pdf_width.get_xticklabels():
            tl.set_visible(False)
        wmax = width_pdf.max()
        ticks = [0, 0.5*wmax, 1.0*wmax]
        tick_labels = ['{0:.1f}'.format(ticks[0]),
                       '{0:.1f}'.format(ticks[1]),
                       '{0:.1f}'.format(ticks[2]),
                        ]
        ax_pdf_width.set_yticks(ticks)
        ax_pdf_width.set_yticklabels(tick_labels)

        for tl in ax_pdf_center.get_yticklabels():
            tl.set_visible(False)
        cmax = center_pdf.max()
        ticks = [0, 0.5*cmax, 1.0*cmax]
        tick_labels = ['{0:.1f}'.format(ticks[0]),
                       '{0:.1f}'.format(ticks[1]),
                       '{0:.1f}'.format(ticks[2]),
                        ]
        ax_pdf_center.set_xticks(ticks)
        ax_pdf_center.set_xticklabels(tick_labels)

        # Show confidence limits
        if center_confint is not None:
            ax_pdf_center.axhspan(center_confint[0] - center_confint[1],
                                  center_confint[0] + center_confint[2],
                                  color='k',
                                  linewidth=1,
                                  alpha=0.2)
            ax_pdf_center.axhline(center_confint[0],
                                  color='k',
                                  linestyle='--',
                                  linewidth=3,
                                  alpha=1)
        if width_confint is not None:
            ax_pdf_width.axvspan(width_confint[0] - width_confint[1],
                                 width_confint[0] + width_confint[2],
                                  color='k',
                                 linewidth=1,
                                  alpha=0.2)
            ax_pdf_width.axvline(width_confint[0],
                                 color='k',
                                 linestyle='--',
                                 linewidth=3,
                                 alpha=1)

    #cb.set_clim(vmin=0.)
    # Write label to colorbar
    #cb.set_label_text(r'log L')

    fractions = np.array([0.95, 0.68])
    levels = (1 + fractions * image.min())

    cs = ax_image.contour(image, levels=levels, origin='lower',
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

    ax_image.clabel(cs, cs.levels, fmt=fmt, fontsize=9, inline=1)

    if filename is not None:
        plt.draw()
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.draw()
        plt.show()
    if returnimage:
        return correlations_image

''' Calculations
'''

def correlate_hi_av(hi_cube=None, hi_velocity_axis=None, hi_noise_cube=None,
        av_image=None, av_image_error=None, velocity_centers=None,
        velocity_widths=None, return_correlations=True, dgr=None,
        plot_results=True, results_filename='', likelihood_filename=None,
        clobber=False, hi_vel_range_conf=0.68):

    '''
    Parameters
    ----------

    Returns
    -------
    hi_vel_range : tuple
        Lower and upper bound of HI velocity range in km/s which provides the
        best correlated N(HI) distribution with Av.
    correlations : array-like, optional
        Array of Pearson correlation coefficients corresponding to each
        permutation through the velocity centers and velocity widths.

    '''

    import numpy as np
    from scipy.stats import pearsonr
    from scipy.stats import kendalltau
    from myimage_analysis import calculate_nhi
    from scipy import signal
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
        # calculate the velocity ranges given a set of centers and widths
        velocity_ranges = np.zeros(shape=[len(velocity_centers) * \
                len(velocity_widths),2])
        count = 0
        for i, center in enumerate(velocity_centers):
            for j, width in enumerate(velocity_widths):
                velocity_ranges[count, 0] = center - width/2.
                velocity_ranges[count, 1] = center + width/2.
                count += 1

        # calculate the correlation coefficient for each velocity range
        correlations = np.zeros(velocity_ranges.shape[0])
        pvalues = np.zeros(velocity_ranges.shape[0])

        for i, velocity_range in enumerate(velocity_ranges):
            nhi_image_temp, nhi_image_error = calculate_nhi(cube=hi_cube,
                    velocity_axis=hi_velocity_axis,
                    velocity_range=velocity_range,
                    noise_cube=hi_noise_cube)

            nhi_image = np.ma.array(nhi_image_temp,
                                    mask=np.isnan(nhi_image_temp))

            # Avoid NaNs
            indices = np.where((nhi_image_temp == nhi_image_temp) & \
                               (av_image == av_image))

            nhi_image_corr = nhi_image_temp[indices]
            nhi_image_error_corr = nhi_image_error[indices]
            av_image_corr = av_image[indices]
            if type(av_image_error) != float:
                av_image_error_corr = av_image_error[indices]
            else:
                av_image_error_corr = av_image_error


            # Create model of Av with N(HI) and DGR
            av_image_model = nhi_image_corr * dgr
            av_image_model_error = nhi_image_error_corr * dgr

            logL = calc_logL(av_image_model,
                             av_image_corr,
                             data_error=av_image_error_corr)

            correlations[i] = -logL

            # Shows progress each 10%
            total = float(correlations.shape[0])
            abs_step = int((total * 1)/10) or 10
            if i and not i % abs_step:
                print "\t{0:.0%} processed".format(i/total)

        # Normalize the log likelihoods
        correlations -= correlations.max()

        # Convert to likelihoods
        correlations = np.exp(correlations)

        # Normalize the likelihoods
        correlations = correlations / np.sum(correlations)

        # Avoid nans
        correlations = np.ma.array(correlations,
                mask=(correlations != correlations))

        # Reshape array
        correlations_image = np.empty((velocity_centers.shape[0],
                                       velocity_widths.shape[0]))
        correlations_image[:,:] = np.NaN
        count = 0
        for i, center in enumerate(velocity_centers):
            for j, width in enumerate(velocity_widths):
                correlations_image[i,j] = correlations[count]
                count += 1

        # Write out fits file of likelihoods
        if write_mle:
            print('Writing likelihood grid to file:')
            print(likelihood_filename)
            header = fits.Header()
            header['NAXIS'] = 2
            header['CTYPE1'] = 'CENTERS'
            header['CTYPE2'] = 'WIDTHS'
            header['CRPIX1'] = 0
            header['CRPIX2'] = 0
            header['CRVAL1'] = velocity_centers[0]
            header['CRVAL2'] = velocity_widths[0]
            header['CDELT1'] = velocity_centers[1] - velocity_centers[0]
            header['CDELT2'] = velocity_widths[1] - velocity_widths[0]

            hdu = fits.PrimaryHDU(correlations_image, header=header)

            hdu.writeto(likelihood_filename, clobber=clobber)

    # Load file of likelihoods
    elif not perform_mle:
        print('Reading likelihood grid file:')
        print(likelihood_filename)

        hdu = fits.open(likelihood_filename)
        correlations_image = hdu[0].data

        if len(velocity_centers) != correlations_image.shape[0] or \
            len(velocity_widths) != correlations_image.shape[1]:
            raise ValueError('Specified parameter grid not the same as in' + \
                    'loaded data likelihoods.')

    # Define parameter resolutions
    delta_center = velocity_centers[1] - velocity_centers[0]
    delta_width = velocity_widths[1] - velocity_widths[0]

    # Derive marginal distributions of both centers and widths
    center_corr = np.sum(correlations_image, axis=1) / \
            np.sum(correlations_image)
    width_corr = np.sum(correlations_image, axis=0) / \
            np.sum(correlations_image)

    # Derive confidence intervals of parameters
    center_confint = threshold_area(velocity_centers,
                                    center_corr,
                                    area_fraction=hi_vel_range_conf)
    width_confint = threshold_area(velocity_widths,
                                   width_corr,
                                   area_fraction=hi_vel_range_conf)

    print('Velocity widths = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(width_confint[0],
                                                    width_confint[2],
                                                    np.abs(width_confint[1])))
    print('Velocity centers = ' + \
            '{0:.2f} +{1:.2f}/-{2:.2f} km/s'.format(center_confint[0],
                                                    center_confint[2],
                                                    np.abs(center_confint[1])))

    # Write PDF
    center = center_confint[0]
    upper_lim = (center_confint[0] + width_confint[0]/2.)
    lower_lim = (center_confint[0] - width_confint[0]/2.)
    upper_lim_error = (center_confint[2]**2 + width_confint[2]**2)**0.5
    lower_lim_error = (center_confint[1]**2 + width_confint[1]**2)**0.5

    vel_range_confint = (lower_lim, upper_lim, lower_lim_error,
            upper_lim_error)

    if plot_results:
        plot_correlations(correlations_image,
                          velocity_centers,
                          velocity_widths,
                          show=0,
                          returnimage=False,
                          filename=results_filename)
        plot_correlations_hist(correlations_image,
                              velocity_centers,
                              velocity_widths,
                              center_pdf=center_corr,
                              width_pdf=width_corr,
                              center_confint=center_confint,
                              width_confint=width_confint,
                              show=0,
                              returnimage=False,
                              filename=results_filename)

    if not return_correlations:
        return vel_range_confint
    else:
        return vel_range_confint, correlations_image, center_corr, width_corr

def calc_logL(model, data, data_error=None):

    '''
    Calculates log likelihood

    http://www.physics.utah.edu/~detar/phys6720/handouts/curve_fit/curve_fit/node2.html

    '''

    if data_error is None:
        data_error = np.std(data)

    logL = -np.sum(-(data - model)**2 / (2 * data_error**2)) / data.size

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
    # 1.0 mag corresponds to SNR = 1 / 0.2 ~ 5
    # (see Table 2 of Ridge et al. 2006).
    indices = np.where((nhi_image_temp == nhi_image_temp) & \
                       (av_image == av_image) & \
                       (av_image_error == av_image_error))

    nhi_image_corr = nhi_image_temp[indices]
    nhi_image_error_corr = nhi_image_error[indices]
    av_image_data = av_image[indices]
    av_image_error_data = av_image_error[indices]

    # Create model image
    av_image_model = nhi_image_corr * dgr
    av_image_error_model = nhi_image_error_corr * dgr

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

    # Step for lowering threshold
    step = (np.max(y) - np.median(y)) / 1000.0

    # initial threshold
    threshold = np.max(y) - step
    threshold_area = 0.0

    # area under whole function
    area = integrate(y, x)

    # Stop when the area below the threshold is greater than the max area
    while threshold_area < area * area_fraction:

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

    x_peak = x[y == y.max()][0]
    low_error, up_error = x_peak - x[bounds_indices[0]], \
                          x[bounds_indices[1]] - x_peak

    return (x_peak, low_error, up_error)

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

    # parameters used in script
    # -------------------------
    # HI velocity integration range
    # Determine HI integration velocity by CO or correlation with Av?
    hi_av_correlation = True
    velocity_centers = np.arange(0, 20, 0.5)
    velocity_widths = np.arange(1, 50, 0.5)

    # Which likelihood fits should be performed?
    core_correlation = 0
    global_correlation = 1

    # Name of property files results are written to
    global_property_file = 'perseus_global_properties.txt'
    core_property_file = 'perseus_core_properties.txt'

    # Threshold of Av below which we expect only atomic gas, in mag
    av_threshold = 1

    # Check if likelihood file already written, rewrite?>
    likelihood_filename = 'perseus_nhi_av_likelihoods'
    clobber=True
    hi_vel_range_conf = 0.95

    # Name of noise cube
    noise_cube_filename = 'perseus_hi_galfa_cube_regrid_planckres_noise.fits'

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/perseus/figures/hi_velocity_range/'
    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    hi_dir = '/d/bip3/ezbc/perseus/data/hi/'
    co_dir = '/d/bip3/ezbc/perseus/data/co/'
    core_dir = '/d/bip3/ezbc/perseus/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/perseus/data/python_output/'
    region_dir = '/d/bip3/ezbc/perseus/data/python_output/ds9_regions/'
    likelihood_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'

    # load Planck Av and GALFA HI images, on same grid
    av_data_planck, av_header = load_fits(av_dir + \
                'perseus_av_planck_5arcmin.fits',
            return_header=True)

    av_error_data_planck, av_error_header = load_fits(av_dir + \
                'perseus_av_error_planck_5arcmin.fits',
            return_header=True)

    hi_data, h = load_fits(hi_dir + \
                'perseus_hi_galfa_cube_regrid_planckres.fits',
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

    dgr = global_props['dust2gas_ratio']['value']

    cores = convert_core_coordinates(cores, h)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'perseus_av_boxes_',
            header = h)

    if core_correlation:
        for core in cores:
            print('\nCalculating for core %s' % core)

            # Grab the mask
            mask = myg.get_polygon_mask(av_data_planck,
                    cores[core]['box_vertices_rotated'])

            indices = ((mask == 0) &\
                       (av_data_planck < av_threshold))

            hi_data_sub = np.copy(hi_data[:, indices])
            noise_cube_sub = np.copy(noise_cube[:, indices])
            av_data_sub = np.copy(av_data_planck[indices])
            av_error_data_sub = np.copy(av_error_data_planck[indices])

            # Define filename for plotting results
            results_filename = figure_dir + 'perseus_logL_%s.png' % core

            # Correlate each core region Av and N(HI) for velocity ranges
            vel_range_confint, correlations, center_corr, width_corr = \
                    correlate_hi_av(hi_cube=hi_data_sub,
                                    hi_velocity_axis=velocity_axis,
                                    hi_noise_cube=noise_cube_sub,
                                    av_image=av_data_sub,
                                    av_image_error=av_error_data_sub,
                                    dgr=dgr,
                                    velocity_centers=velocity_centers,
                                    velocity_widths=velocity_widths,
                                    return_correlations=True,
                                    plot_results=True,
                                    results_filename=results_filename,
                                    likelihood_filename=likelihood_dir + \
                                            likelihood_filename + \
                                            '{0:s}.fits'.format(core),
                                    clobber=clobber,
                                    hi_vel_range_conf=hi_vel_range_conf)

            print('HI velocity integration range:')
            print('%.1f to %.1f km/s' % (vel_range_confint[0],
                                         vel_range_confint[1]))

            cores[core]['hi_velocity_range'] = vel_range_confint[0:2]
            cores[core]['hi_velocity_range_error'] = vel_range_confint[2:]
            cores[core]['center_corr'] = center_corr.tolist()
            cores[core]['width_corr'] = width_corr.tolist()
            cores[core]['vel_centers'] = velocity_centers.tolist()
            cores[core]['vel_widths'] = velocity_widths.tolist()

        with open(core_dir + core_property_file, 'w') as f:
            json.dump(cores, f)

    if global_correlation:
        print('\nCalculating correlations globally')

        indices = ((av_data_planck < av_threshold))

        hi_data_sub = np.copy(hi_data[:, indices])
        noise_cube_sub = np.copy(noise_cube[:, indices])
        av_data_sub = np.copy(av_data_planck[indices])
        av_error_data_sub = np.copy(av_error_data_planck[indices])

        # Define filename for plotting results
        results_filename = figure_dir + 'perseus_logL_global.png'

        # Correlate each core region Av and N(HI) for velocity ranges
        vel_range_confint, correlations, center_corr, width_corr = \
                correlate_hi_av(hi_cube=hi_data_sub,
                                hi_velocity_axis=velocity_axis,
                                hi_noise_cube=noise_cube_sub,
                                av_image=av_data_sub,
                                av_image_error=av_error_data_sub,
                                dgr=dgr,
                                velocity_centers=velocity_centers,
                                velocity_widths=velocity_widths,
                                return_correlations=True,
                                plot_results=True,
                                results_filename=results_filename,
                                likelihood_filename=likelihood_dir + \
                                        likelihood_filename + '_global.fits',
                                clobber=clobber,
                                hi_vel_range_conf=hi_vel_range_conf)

        '''
        fit_hi_vel_range(guesses=(0, 30),
                         av_image=av_data_sub,
                         av_image_error=av_error_data_sub,
                         hi_cube=hi_data_sub,
                         hi_velocity_axis=velocity_axis,
                         hi_noise_cube=noise_cube_sub,
                         dgr=dgr)
        '''

        print('HI velocity integration range:')
        print('%.1f to %.1f km/s' % (vel_range_confint[0],
                                     vel_range_confint[1]))

        global_props['hi_velocity_range'] = vel_range_confint[0:2]
        global_props['hi_velocity_range_error'] = vel_range_confint[2:]
        global_props['hi_velocity_range_conf'] = hi_vel_range_conf
        global_props['center_corr'] = center_corr.tolist()
        global_props['width_corr'] = width_corr.tolist()
        global_props['vel_centers'] = velocity_centers.tolist()
        global_props['vel_widths'] = velocity_widths.tolist()

        with open(property_dir + global_property_file, 'w') as f:
            json.dump(global_props, f)

if __name__ == '__main__':
    main()








