#!/usr/bin/python

''' Calculates the N(HI) map for multicloud

'''

import numpy as np
import matplotlib
matplotlib.use('Agg')

from astropy.io import fits
import numpy as np
import warnings
warnings.filterwarnings('ignore')


''' Plotting Functions
'''

def plot_nhi_image(nhi_image=None, header=None, contour_image=None,
        av_image=None, cores=None, props=None, regions=None, title=None,
        limits=None, contours=None, boxes=False, filename=None,
        show=True, hi_vlimits=None, av_vlimits=None,):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    import pyfits as fits
    import matplotlib.pyplot as plt
    import myplotting as myplt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Set up plot aesthetics
    plt.clf(); plt.close()

    # Color map
    cmap = plt.cm.gnuplot
    #cmap = myplt.reverse_colormap(plt.cm.copper)
    cmap = plt.cm.copper

    # Create figure instance
    fig = plt.figure(figsize=(7.5, 3.6*2))

    if av_image is not None:
        nrows_ncols=(2,1)
        ngrids=2
    else:
        nrows_ncols=(1,1)
        ngrids=1

    axes = AxesGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0.1,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # ------------------
    # NHI image
    # ------------------
    # create axes
    ax = axes[0]

    # show the image
    im = ax.imshow(nhi_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=hi_vlimits[0],
            vmax=hi_vlimits[1],
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)

    ax.locator_params(nbins=6)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')

    # plot limits
    if limits is not None:
        limits_pix = myplt.convert_wcs_limits(limits, header, frame='fk5')
        ax.set_xlim(limits_pix[0],limits_pix[1])
        ax.set_ylim(limits_pix[2],limits_pix[3])

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Write label to colorbar
    cb.set_label_text(r'$N$(H\textsc{i}) [10$^{20}$ cm$^{-2}$]',)

    # Convert sky to pix coordinates
    #wcs_header = pywcs.WCS(header)
    if type(cores) is dict:
        for core in cores:
            pix_coords = cores[core]['center_pixel']

            anno_color = (0.3, 0.5, 1)

            ax.scatter(pix_coords[0],pix_coords[1],
                    color=anno_color,
                    s=200,
                    marker='+',
                    linewidths=2)

            ax.annotate(core,
                    xy=[pix_coords[0], pix_coords[1]],
                    xytext=(5,5),
                    textcoords='offset points',
                    color=anno_color)

            if boxes:
                rect = ax.add_patch(Polygon(
                    cores[core]['box_vertices'][:, ::-1],
                        facecolor='none',
                        edgecolor=anno_color))

    if regions is not None:
        for region in props['region_name_pos']:
            vertices = np.copy(regions[region]['poly_verts']['pixel'])
            rect = ax.add_patch(Polygon(
                    vertices[:, ::-1],
                    facecolor='none',
                    edgecolor='w'))

    # ------------------
    # Av image
    # ------------------
    if av_image is not None:
        # create axes
        ax = axes[1]
        # show the image
        im = ax.imshow(av_image,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                vmin=av_vlimits[0],
                vmax=av_vlimits[1],
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension [J2000]',)
        ax.set_ylabel('Declination [J2000]',)

        ax.locator_params(nbins=6)

        # colorbar
        cb = ax.cax.colorbar(im)
        cmap.set_bad(color='w')

        # plot limits
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits, header, frame='fk5')
            ax.set_xlim(limits_pix[0],limits_pix[1])
            ax.set_ylim(limits_pix[2],limits_pix[3])

        ax.tick_params(axis='xy', which='major', colors='w')

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='r')

        # Write label to colorbar
        cb.set_label_text(r'$A_V$ [mag]',)

        # Convert sky to pix coordinates
        #wcs_header = pywcs.WCS(header)
        if type(cores) is dict:
            for core in cores:
                pix_coords = cores[core]['center_pixel']

                anno_color = (0.3, 0.5, 1)

                if boxes:
                    rect = ax.add_patch(Polygon(
                        cores[core]['box_vertices'][:, ::-1],
                            facecolor='none',
                            edgecolor='w'))

        if regions is not None:
            for region in props['region_name_pos']:
                vertices = np.copy(regions[region]['poly_verts']['pixel'])
                rect = ax.add_patch(Polygon(
                        vertices[:, ::-1],
                        facecolor='none',
                        edgecolor='w'))

            if props is not None:
                for region in props['region_name_pos']:
                    if region == 'taurus1':
                        region = 'taurus 1'
                    if region == 'taurus2':
                        region = 'taurus 2'
                    ax.annotate(region.capitalize(),
                                xy=props['region_name_pos'][region]['pixel'],
                                xytext=(0,0),
                                textcoords='offset points',
                                color='w',
                                fontsize=10,
                                zorder=10)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        fig.show()

def plot_av_model(av_image=None, header=None, contour_image=None,
        av_model=None, hi_velocity_axis=None, vel_range=None, hi_spectrum=None,
        co_spectrum=None, co_velocity_axis=None, co_scaling=50,
        hi_limits=None, cores=None, results=None, title=None, limits=None,
        contours=None, boxes=False, savedir='./', filename=None, show=True,
        plot_residuals=True):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import axes
    import pyfits as fits
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon

    # Set up plot aesthetics
    plt.clf()
    plt.close()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 12
    if plot_residuals:
        figsize = (17,8)
    else:
        figsize = (13,10)

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
              'figure.figsize': figsize,
              'figure.titlesize': font_scale
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    #==========================================================================
    # Av maps
    #==========================================================================

    if plot_residuals:
        nrows_ncols=(1,3)
        ngrids=3
    else:
        nrows_ncols=(1,2)
        ngrids=2

    axes = axes(fig, (2,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode='each',
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0.4,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # ------------------
    # Av model
    # ------------------
    # create axes
    ax = axes[0]
    cmap = cm.Greys # colormap
    cmap = cm.gnuplot # colormap
    cmap.set_bad(color='w')
    # show the image
    vmax = np.max((np.max(av_image[av_image == av_image]),
                   np.max(av_model[av_model == av_model])))
    vmax = 1.4
    im = ax.imshow(av_model,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',)
    ax.set_ylabel('Declination (J2000)',)

    ax.set_title(r'Model $A_V$')

    # colorbar
    cb = ax.cax.colorbar(im)
    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    if type(cores) is dict:
        for core in cores:
            pix_coords = cores[core]['center_pixel']

            anno_color = (0.3, 0.5, 1)

            ax.scatter(pix_coords[0],pix_coords[1],
                    color=anno_color,
                    s=200,
                    marker='+',
                    linewidths=2)

            ax.annotate(core,
                    xy=[pix_coords[0], pix_coords[1]],
                    xytext=(5,5),
                    textcoords='offset points',
                    color=anno_color)

            if boxes:
                rect = ax.add_patch(Polygon(
                    cores[core]['box_vertices'][:, ::-1],
                        facecolor='none',
                        edgecolor=anno_color))

    # ------------------
    # Av image
    # ------------------
    if av_image is not None:
        # create axes
        ax = axes[1]
        # show the image
        im = ax.imshow(av_image,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                vmin=0,
                vmax=vmax,
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension (J2000)',)
        ax.set_ylabel('Declination (J2000)',)

        ax.set_title(r'Observed $A_V$')

        # colorbar
        cb = ax.cax.colorbar(im)
        #cmap.set_bad(color='c')

        # plot limits
        if limits is not None:
            ax.set_xlim(limits[0],limits[2])
            ax.set_ylim(limits[1],limits[3])

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='r')

        # Write label to colorbar
        #cb.set_label_text(r'$A_V$ (mag)',)

    # ------------------
    # Av residuals
    # ------------------
    if plot_residuals:
        ax = axes[2]
        #cmap = cm.Greys # colormap
        # show the image
        vmax = np.max((np.max(av_image[av_image == av_image]),
                       np.max(av_model[av_model == av_model])))
        vmax = 1.4
        im = ax.imshow(av_image - av_model,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                vmin=-0.25,
                vmax=0.5,
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension (J2000)',)
        ax.set_ylabel('Declination (J2000)',)

        ax.set_title(r'Residual $A_V$')

        # colorbar
        cb = ax.cax.colorbar(im)
        #cmap.set_bad(color='w')
        # plot limits
        if limits is not None:
            ax.set_xlim(limits[0],limits[2])
            ax.set_ylim(limits[1],limits[3])


    # ==========================================================================
    # HI Spectrum
    # ==========================================================================

    vel_range = np.asarray(results['hi_velocity_range'][:2])

    # create axes
    ax = fig.add_subplot(2,1,2)

    ax.plot(hi_velocity_axis,
            hi_spectrum,
            color='k',
            drawstyle = 'steps-mid',
            label='HI',
            )
    if co_spectrum is not None:
        ax.plot(co_velocity_axis,
                co_spectrum*co_scaling,
                color='r',
                drawstyle='steps-mid',
                label=r'CO x ' + str(co_scaling),
                )
        ax.legend(loc='upper right')

    # Plot velocity range
    if vel_range.ndim == 1:
        ax.axvspan(vel_range[0], vel_range[1], color='k', alpha=0.3)
    elif vel_range.ndim == 2:
        for i in xrange(0, vel_range.shape[0]):
            ax.axvspan(vel_range[i, 0], vel_range[i, 1], color='k', alpha=0.3)

    # Plot center
    ax.axvline(results['co_center']['value'], color='k', linestyle='--', )

    if hi_limits is not None:
        ax.set_xlim(hi_limits[0], hi_limits[1])
        ax.set_ylim(hi_limits[2], hi_limits[3])

    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel(r'T$_b$ (K)')

    # Plot results
    if results is not None:
        av_thres = results['av_threshold']['value']
        co_thres = results['co_threshold']['value']
        if av_thres > 15:
            av_thres = None
        text = ''
        text += r'N$_{\rm pix}$ = ' + \
                 '{0:.0f}'.format(results['npix'])
        text += '\n'
        if av_thres is not None:
            text += r'$A_V$ threshold = ' + \
                    '{0:.1f} mag'.format(av_thres)
        else:
            text += r'$A_V$ threshold = ' + \
                    '{0:s}'.format(av_thres)
        text += '\n'
        text += r'CO threshold = ' + \
                '{0:.1f} K km/s'.format(co_thres)
        text += '\n'
        text += r'DGR = {0:.2f} '.format(results['dust2gas_ratio']['value']) + \
                r'$\times$ 10$^{-20}$ (cm$^2$ mag$^1$)'
        text += '\n'
        if vel_range.ndim == 1:
            text += r'Velocity range = ' + \
                    '{0:.1f} to {1:.1f} km/s'.format(vel_range[0], vel_range[1])
            text += '\n'
        text += r'$\chi^2$ / $\nu$ = {0:.1f}'.format(results['chisq'])

        ax.annotate(text,
                xytext=(0.03, 0.95),
                xy=(0.03, 0.95),
                textcoords='axes fraction',
                xycoords='axes fraction',
                color='k',
                fontsize=font_scale*0.75,
                bbox=dict(boxstyle='round',
                          facecolor='w',
                          alpha=0.8),
                horizontalalignment='left',
                verticalalignment='top',
                )

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        plt.show()

def plot_avmod_vs_av(avmod_images, av_images, av_errors=None, limits=None,
        fit=True, savedir='./', filename=None, show=True,
        scale=('linear','linear'), title = '', gridsize=(100,100), std=None):

    # Import external modules
    import numpy as np
    import math
    import pyfits as fits
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import axes
    from myscience.krumholz09 import calc_T_cnm
    from matplotlib import cm

    n = int(np.ceil(len(av_images)**0.5))
    if n**2 - n > len(av_images):
        nrows = n - 1
        ncols = n
        y_scaling = 1.0 - 1.0 / n
    else:
        nrows, ncols = n, n
        y_scaling = 1.0

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
              'text.usetex': True,
              'figure.figsize': (5, 5 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    axes = axes(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=len(av_images),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 #cbar_mode='single',
                 cbar_pad=0.1,
                 cbar_size=0.2,
                 )

    # Cycle through lists
    for i in xrange(len(av_images)):
        av = av_images[i]
        avmod = avmod_images[i]
        av_error = av_errors[i]

        # Drop the NaNs from the images
        if type(av_error) is float:
            indices = np.where((av == av) &\
                               (avmod == avmod)
                               )

        if type(av_error) is np.ndarray or \
                type(av_error) is np.ma.core.MaskedArray or \
                type(avmod_error) is np.ndarray or \
                type(avmod_error) is np.ma.core.MaskedArray:
            indices = np.where((av == av) &\
                               (avmod == avmod) &\
                               (av_error == av_error)
                               )

        av_nonans = av[indices]
        avmod_nonans = avmod[indices]

        # Fix error data types
        if type(av_error) is np.ndarray:
            av_error_nonans = av_error[indices]
        else:
            av_error_nonans = np.array(av_error[indices])

        # Create plot
        ax = axes[i]

        if 1:
            image = ax.errorbar(av_nonans.ravel(),
                    avmod_nonans.ravel(),
                    xerr=(av_error_nonans.ravel()),
                    alpha=0.2,
                    color='k',
                    marker='^',
                    ecolor='k',
                    linestyle='None',
                    markersize=3
                    )
        if 0:
            image = ax.hexbin(av_nonans.ravel(),
                avmod_nonans.ravel(),
                #norm=matplotlib.colors.LogNorm(),
                mincnt=1,
                yscale=scale[1],
                xscale=scale[0],
                gridsize=gridsize,
                cmap=cm.Greys,
                #cmap=cm.gist_stern,
                )
            cb = ax.cax.colorbar(image,)
            # Write label to colorbar
            cb.set_label_text('Bin Counts',)

        # Plot sensitivies
        av_limit = np.median(av_errors[0])
        ax.axvline(av_limit, color='k', linestyle='--')

        # Plot 1 to 1 pline
        ax.plot((0, 10),
                (0, 10),
                color='0.5',
                linewidth=3,
                alpha=0.5)
        if std is not None:
            ax.fill_between((-std, 10),
                            (0, 10 + std),
                            (-2*std, 10 - std),
                            color='0.2',
                            alpha=0.2,
                            edgecolor='none'
                            )

        # Annotations
        anno_xpos = 0.95

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'$A_V$ Data (mag)')
        ax.set_ylabel(r'$A_V$ Model (mag)')
        #ax.set_title(core_names[i])

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        fig.show()

def plot_power_spectrum(image, title=None, filename_prefix=None,
        filename_suffix='.png', show=False, savedir='./'):

    '''
    Plots power spectrum derived from a fourier transform of the image.

    '''

    # import external modules
    import numpy as np
    from agpy import azimuthalAverage as radial_average
    from scipy import fftpack
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import axes
    from matplotlib import cm

    if 0:
        plt.close(); plt.clf()
        plt.imshow(image)
        plt.show()

    image[np.isnan(image)] = 1e10

    # Determine power spectrum
    # -------------------------------------------------------------------------
    # Take the fourier transform of the image.
    #F1 = fftpack.fft2(np.ma.array(image, mask=np.isnan(image)))
    F1 = fftpack.fft2(image)

    # Now shift the quadrants around so that low spatial frequencies are in
    # the center of the 2D fourier transformed image.
    F2 = fftpack.fftshift(F1)

    # Calculate a 2D power spectrum
    psd2D = np.abs(F2)**2

    if 0:
        plt.close(); plt.clf()
        plt.imshow(psd2D)
        plt.show()

    power_spectrum = radial_average(psd2D, interpnan=True)

    # Write frequency in arcmin
    freq = fftpack.fftfreq(len(power_spectrum))
    freq *= 5.0

    # Simulate power spectrum for white noise
    noise_image = np.random.normal(scale=0.1, size=image.shape)

    F1 = fftpack.fft2(noise_image)

    # Now shift the quadrants around so that low spatial frequencies are in
    # the center of the 2D fourier transformed image.
    F2 = fftpack.fftshift(F1)

    # Calculate a 2D power spectrum
    psd2D_noise = np.abs(F2)**2

    power_spectrum_noise = radial_average(psd2D_noise, interpnan=True)

    # Plot power spectrum 1D
    # -------------------------------------------------------------------------
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
              'figure.figsize': (5, 5),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    nrows = 1; ncols = 1; ngrids = 1;

    axes = axes(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=ngrids,
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 #cbar_mode='single',
                 cbar_pad=0.1,
                 cbar_size=0.2,
                 )

    ax = axes[0]

    ax.plot(freq, power_spectrum / np.nanmax(power_spectrum),
            color='k',
            linestyle='-',
            linewidth=1.5,
            drawstyle='steps-mid',
            label='Data Residuals')

    ax.plot(freq, power_spectrum_noise / np.nanmax(power_spectrum_noise),
            color='r',
            linestyle='-',
            linewidth=0.4,
            drawstyle='steps-mid',
            label='White Noise Residuals')

    #ax.set_xscale('log')
    ax.legend(loc='best')
    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Spatial Frequency [1/arcmin]')
    ax.set_ylabel('Normalized Power Spectrum')
    ax.set_xlim(0, 0.4)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename_prefix is not None:
        plt.savefig(savedir + filename_prefix + '_1D' + filename_suffix,
                bbox_inches='tight')
    if show:
        plt.show()

    # Plot power spectrum image
    # -------------------------
    # Create figure instance
    fig = plt.figure()

    nrows = 1; ncols = 1; ngrids = 1;

    axes = axes(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=ngrids,
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 #cbar_mode='single',
                 cbar_pad=0.1,
                 cbar_size=0.2,
                 )

    ax = axes[0]

    extent = [- image.shape[0] / 2.0, + image.shape[0] / 2.0,
              - image.shape[1] / 2.0, + image.shape[1] / 2.0]

    ax.imshow(psd2D,
              origin='lower',
              cmap=cm.gist_heat,
              norm=matplotlib.colors.LogNorm(),
              extent=extent
              )

    #ax.set_xscale('log')
    #ax.legend(loc='center right')
    #ax.set_yscale('log')
    ax.set_xlabel('Spatial Frequency in Right Ascension')
    ax.set_ylabel('Spatial Frequency in Declination')
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename_prefix is not None:
        plt.savefig(savedir + filename_prefix + '_2D' + filename_suffix,
                bbox_inches='tight')
    if show:
        plt.show()

''' DS9 Region and Coordinate Functions
'''

def convert_limit_coordinates(prop_dict,
        coords=('region_limit', 'co_noise_limits', 'plot_limit'), header=None):

    # Initialize pixel keys
    for coord in coords:

        if coord == 'region_limit' or coord == 'plot_limit':
            prop_dict[coord].update({'pixel': []})
            limit_wcs = prop_dict[coord]['wcs']

            for limits in limit_wcs:
                # convert centers to pixel coords
                limit_pixels = get_pix_coords(ra=limits[0],
                                             dec=limits[1],
                                             header=header)[:2].tolist()

                prop_dict[coord]['pixel'].append(limit_pixels[0])
                prop_dict[coord]['pixel'].append(limit_pixels[1])
        elif coord == 'co_noise_limits':
            prop_dict[coord].update({'pixel': []})
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
        elif coord == 'region_name_pos':

            # convert centers to pixel coords
            for region in prop_dict[coord]:
                prop_dict[coord][region].update({'pixel': []})

                coord_wcs = prop_dict[coord][region]['wcs']

                coord_pixel = get_pix_coords(ra=coord_wcs[0],
                                             dec=coord_wcs[1],
                                             header=header)[:2].tolist()

                prop_dict[coord][region]['pixel'].append(coord_pixel[0])
                prop_dict[coord][region]['pixel'].append(coord_pixel[1])

    return prop_dict

def convert_core_coordinates(cores, header):

    for core in cores:
        cores[core].update({'box_pixel': 0})
        cores[core].update({'center_pixel': 0})

        #box_wcs = cores[core]['box_wcs']
        #box_pixel = len(box_wcs) * [0,]
        center_wcs = cores[core]['center_wcs']

        # convert centers to pixel coords
        center_pixel = get_pix_coords(ra=center_wcs[0],
                                      dec=center_wcs[1],
                                      header=header)[:2]
        cores[core]['center_pixel'] = center_pixel.tolist()

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

    import pyfits as fits
    from astropy.io import fits

    f = fits.open(filename)
    if return_header:
        return f[0].data,f[0].header
    else:
        return f[0].data

def get_sub_image(image, indices):

    return image[indices[1]:indices[3],
            indices[0]:indices[2]]

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

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec), or Ra in degrees
    and dec in degrees.
    '''

    import pywcsgrid2 as wcs
    import pywcs
    from astropy.wcs import WCS

    # convert to degrees if ra and dec are array-like
    try:
        if len(ra) == 3 and len(dec) == 3:
            ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)
        else:
            raise ValueError('RA and Dec must be in (hrs,min,sec) and' + \
                    ' (deg,arcmin,arcsec) or in degrees.')
    except TypeError:
        ra_deg, dec_deg = ra, dec

    #wcs_header = pywcs.WCS(header)
    wcs_header = WCS(header)
    #pix_coords = wcs_header.wcs_sky2pix([[ra_deg, dec_deg, 0]], 0)[0]
    pix_coords = wcs_header.wcs_world2pix([[ra_deg, dec_deg],], 0)[0]

    return np.hstack((pix_coords, -1))

def hrs2degs(ra=None, dec=None):
    ''' Ra and dec tuples in hrs min sec and deg arcmin arcsec.
    '''

    ra_deg = 15*(ra[0] + ra[1]/60. + ra[2]/3600.)
    dec_deg = dec[0] + dec[1]/60. + dec[2]/3600.

    return (ra_deg, dec_deg)


'''
The main script
'''

def main(dgr=None, vel_range=None, vel_range_type='single', region=None,
        av_data_type='planck'):
    ''' Executes script.

    Parameters
    ----------
    dgr : float
        If None, pulls best-fit value from properties.
    vel_range : tuple
        If None, pulls best-fit value from properties.
    '''

    # import external modules
    #import pyfits as fits
    from astropy.io import fits
    import numpy as np
    from mycoords import make_velocity_axis
    import mygeometry as myg
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    import json
    from os import system,path

    # Script parameters
    # -----------------
    # Name of noise cube
    noise_cube_filename = 'multicloud_hi_galfa_cube_regrid_planckres_noise.fits'

    # Use Planck dust Av map or Kainulainen 2009 optical extinction Av map?
    # options are 'planck' or 'lee12'
    #av_data_type = 'lee12'
    #av_data_type = 'planck'

    # Global parameter file
    prop_file = 'multicloud_global_properties'

    # Regions, regions to edit the global properties with
    if region == 1:
        region_limit = {'wcs' : (((5, 10, 0), (19, 0, 0)),
                                 ((4, 30, 0), (27, 0, 0))),
                          'pixel' : ()
                         }
    elif region == 2:
        region_limit = {'wcs' : (((4, 30, 0), (19, 0, 0)),
                                 ((3, 50, 0), (29, 0, 0))),
                          'pixel' : ()
                        }
    elif region == 3:
        region_limit = {'wcs' : (((4, 30, 0), (29, 0, 0)),
                                 ((3, 50, 0), (33, 0, 0))),
                          'pixel' : ()
                        }
    else:
        region_limit = None

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/multicloud/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/multicloud/figures/'
    av_dir = '/d/bip3/ezbc/multicloud/data/av/'
    hi_dir = '/d/bip3/ezbc/multicloud/data/hi/'
    co_dir = '/d/bip3/ezbc/multicloud/data/co/'
    core_dir = '/d/bip3/ezbc/multicloud/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'

    # load Planck Av and GALFA HI images, on same grid
    nhi, nhi_header = load_fits(nhi_dir + \
                '.fits',
            return_header=True)

    props['plot_limit']['wcs'] = (((5, 20, 0), (19, 0 ,0)),
                                  ((2, 30, 0, (37, 0, 0))
                                  )
    props['region_name_pos'] = {
             #'taurus 1' : {'wcs' : ((3, 50,  0),
             #                       (21.5, 0, 0)),
             #             },
             #'taurus 2' : {'wcs' : ((5, 10,  0),
             #                       (21.5, 0, 0)),
             #             },
             'taurus' : {'wcs' : ((4, 40,  0),
                                  (21, 0, 0)),
                          },
             'perseus' : {'wcs' : ((3, 30,  0),
                                   (26, 0, 0)),
                          },
             #'perseus 1' : {'wcs' : ((3, 0,  0),
             #                      (34, 0, 0)),
             #             },
             #'perseus 2' : {'wcs' : ((3, 10,  0),
             #                      (22.5, 0, 0)),
             #             },
             'california' : {'wcs' : ((4, 28,  0),
                                      (34, 0, 0)),
                             },
             }

    # Change WCS coords to pixel coords of images
    props = convert_limit_coordinates(props,
                                      header=av_header,
                                      coords=('region_limit',
                                              'co_noise_limits',
                                              'plot_limit',
                                              'region_name_pos'))



    # Plot
    figure_types = ['png', 'pdf']
    for figure_type in figure_types:
        plot_nhi_images(nhi_image=nhi_image,
                        header=av_header,
                        av_image=av_image,
                        #limits=props['plot_limit']['pixel'],
                        #limits=[76, 43.5, 19.5, 38,],
                        regions=props['regions'],
                        props=props,
                        #av_vlimits=(0,15.5),
                        #av_vlimits=(0.1,30),
                        filename=filename,
                        show=False)

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
        results['data']['hi_data_error'] = \
                fits.getdata(results['filenames']['hi_error_filename'])
        results['data']['co_data'], results['data']['co_header'] = \
                fits.getdata(results['filenames']['co_filename'], header=True)

    return results

def main()

    # define constants
    DIR_FIGURES = '/d/bip3/ezbc/multicloud/figures/'
    FILENAME_EXT = '_planck_noint_gaussrange_isotropic.pickle'
    FILENAME_PLOT_BASE = DIR_FIGURES + 'maps/h2_maps'
    PLOT_FILETYPES = ['png', 'pdf']

    CLOUD_NAMES = ['california', 'perseus', 'taurus']

    plot_dict = {}
    for cloud_name in CLOUD_NAMES:
        plot_dict[cloud_name] = {}

        results = load_results(cloud_name + FILENAME_EXT)

        hi_sd = results['data_products']['hi_sd']
        h2_sd = results['data_products']['h2_sd']

        plot_dict['hi_sd'] = hi_sd
        plot_dict['h2_sd'] = h2_sd

    for filetype in PLOT_FILETYPES:
        plot_h2_sd(plot_dict
                   filename=FILENAME_PLOT_BASE + filetype,
                   )


if __name__ == '__main__':
    main()

