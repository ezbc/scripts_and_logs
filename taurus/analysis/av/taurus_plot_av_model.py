#!/usr/bin/python

''' Calculates the N(HI) map for Taurus

'''

import numpy as np

''' Plotting Functions
'''

def plot_nhi_image(nhi_image=None, header=None, contour_image=None,
        av_image=None,
        cores=None, title=None, limits=None,
        contours=None, boxes=False, savedir='./', filename=None, show=True):

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
              'figure.figsize': (15, 7),
              'figure.titlesize': font_scale
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    if av_image is not None:
        nrows_ncols=(1,2)
        ngrids=2
    else:
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

    # ------------------
    # NHI image
    # ------------------
    # create axes
    ax = imagegrid[0]
    cmap = cm.Greys # colormap
    # show the image
    im = ax.imshow(nhi_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=0,
            #norm=matplotlib.colors.LogNorm()
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

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Write label to colorbar
    cb.set_label_text(r'N(HI) $\times$ 10$^{20}$ cm$^{-2}$',)

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
        ax = imagegrid[1]
        # show the image
        im = ax.imshow(av_image,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                #norm=matplotlib.colors.LogNorm()
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

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='r')

        # Write label to colorbar
        cb.set_label_text(r'$A_V$ (mag)',)

        # Convert sky to pix coordinates
        wcs_header = pywcs.WCS(header)
        if type(cores) is dict:
            for core in cores:
                pix_coords = cores[core]['center_pixel']

                anno_color = (0.3, 0.5, 1)

                if boxes:
                    rect = ax.add_patch(Polygon(
                        cores[core]['box_vertices'][:, ::-1],
                            facecolor='none',
                            edgecolor=anno_color))

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
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
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
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

    imagegrid = ImageGrid(fig, (2,1,1),
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
    ax = imagegrid[0]
    cmap = cm.Greys # colormap
    cmap = cm.gist_stern # colormap
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
        ax = imagegrid[1]
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
    # Av image
    # ------------------
    if plot_residuals:
        ax = imagegrid[2]
        #cmap = cm.Greys # colormap
        # show the image
        vmax = np.max((np.max(av_image[av_image == av_image]),
                       np.max(av_model[av_model == av_model])))
        vmax = 1.4
        im = ax.imshow(av_image - av_model,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                #vmin=0,
                #vmax=vmax,
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
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
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

    imagegrid = ImageGrid(fig, (1,1,1),
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
        ax = imagegrid[i]

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
    from mpl_toolkits.axes_grid1 import ImageGrid
    from matplotlib import cm

    if 0:
        plt.close(); plt.clf()
        plt.imshow(image)
        plt.show()

    image[np.isnan(image)] = 1e8

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

    imagegrid = ImageGrid(fig, (1,1,1),
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

    ax = imagegrid[0]

    ax.plot(power_spectrum / np.nanmax(power_spectrum),
            color='k',
            linestyle='-',
            linewidth=1.5,
            drawstyle='steps-mid',
            label='Data Residuals')

    ax.plot(power_spectrum_noise / np.nanmax(power_spectrum_noise),
            color='r',
            linestyle='-',
            linewidth=0.4,
            drawstyle='steps-mid',
            label='White Noise Residuals')

    #ax.set_xscale('log')
    ax.legend(loc='center right')
    ax.set_yscale('log')
    ax.set_xlabel('Spatial Frequency')
    ax.set_ylabel('Normalized Power Spectrum')

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

    imagegrid = ImageGrid(fig, (1,1,1),
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

    ax = imagegrid[0]

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
    import pyfits as pf
    import numpy as np
    from mycoords import make_velocity_axis
    import mygeometry as myg
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    import json

    # Script parameters
    # -----------------
    # Name of noise cube
    noise_cube_filename = 'taurus_hi_galfa_cube_regrid_planckres_noise.fits'

    # Use Planck dust Av map or Kainulainen 2009 optical extinction Av map?
    # options are 'planck' or 'k09'
    #av_data_type = 'k09'
    #av_data_type = 'planck'

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
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/hi/'
    co_dir = '/d/bip3/ezbc/taurus/data/co/'
    core_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/taurus/data/python_output/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'

    # load Planck Av and GALFA HI images, on same grid
    if av_data_type == 'k09':
        av_image, av_header = load_fits(av_dir + \
                    'taurus_av_k09_regrid_planckres.fits',
                return_header=True)
        av_image_error = 0.1 * np.ones(av_image.shape)
    else:
        av_image, av_header = load_fits(av_dir + \
                    'taurus_av_planck_5arcmin.fits',
                return_header=True)

        av_image_error, av_error_header = load_fits(av_dir + \
                    'taurus_av_error_planck_5arcmin.fits',
                return_header=True)


    hi_cube, hi_header = load_fits(hi_dir + \
                'taurus_hi_galfa_cube_regrid_planckres.fits',
            return_header=True)

    hi_noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    co_data, co_header = load_fits(co_dir + \
                'taurus_co_cfa_cube_regrid_planckres.fits',
            return_header=True)

    # Load global properties of cloud
    # global properties written from script
    # 'av/taurus_analysis_global_properties.txt'
    prop_file = 'taurus_global_properties_' + av_data_type
    if region is not None:
        likelihood_filename += '_region{0:.0f}'.format(region)
        results_filename += '_region{0:.0f}'.format(region)
    with open(property_dir + prop_file + '.txt', 'r') as f:
        props = json.load(f)
        '''
        dgr = props['dust2gas_ratio_max']['value']
        vel_center = props['hi_velocity_center_max']['value']
        vel_width = props['hi_velocity_center_max']['value']
        vel_range = np.asarray(props['hi_velocity_range_max']['value'])
        dgr = dgr/3
        scale = 1.
        vel_range = (vel_center - scale*vel_width/2.0,
                     vel_center + scale*vel_width/2.0)
        vel_range = (0, 15)
        '''

    if vel_range is not None:
        props['hi_velocity_range'] = vel_range
    else:
        vel_range = props['hi_velocity_range']
    if dgr is not None:
        props['dust2gas_ratio']['value'] = dgr
    else:
        dgr = props['dust2gas_ratio']['value']

    # define core properties
    with open(core_dir + 'taurus_core_properties.txt', 'r') as f:
        cores = json.load(f)

    # make velocity axis for hi cube
    velocity_axis = make_velocity_axis(hi_header)
    # make velocity axis for co cube
    co_velocity_axis = make_velocity_axis(co_header)

    # Write core coordinates in pixels
    cores = convert_core_coordinates(cores, hi_header)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'taurus_av_boxes_',
            header = hi_header)

    # create nhi image
    nhi_image, nhi_image_error = calculate_nhi(cube=hi_cube,
            velocity_axis=velocity_axis,
            velocity_range=vel_range,
            header=hi_header,
            noise_cube=hi_noise_cube)

    # create model av map
    av_model = nhi_image * dgr

    # Mask the images based on av trheshol
    co_data_nonans = np.copy(co_data)
    co_data_nonans[np.isnan(co_data_nonans)] = 0.0
    co_mom0 = np.sum(co_data_nonans, axis=0)
    mask = ((av_image > props['av_threshold']['value']) & \
            (co_mom0 > props['co_threshold']['value']))

    # Derive relevant region
    pix = props['region_limit']['pixel']
    region_vertices = ((pix[1], pix[0]),
                       (pix[1], pix[2]),
                       (pix[3], pix[2]),
                       (pix[3], pix[0])
                       )

    # block offregion
    region_mask = myg.get_polygon_mask(av_image, region_vertices)

    print('\nRegion size = ' + \
          '{0:.0f} pix'.format(region_mask[region_mask == 1].size))

    if vel_range_type == 'single':
        print('\nHI velocity integration range:')
        print('%.1f to %.1f km/s' % (vel_range[0],
                                     vel_range[1]))
    elif vel_range_type == 'multiple':
        print('\nHI velocity integration ranges:')
        for i in xrange(0, vel_range.shape[0]):
            print('%.1f to %.1f km/s' % (vel_range[i, 0],
                                         vel_range[i, 1]))

    print('\nDGR:')
    print('%.2f x 10^-20 cm^2 mag' % (dgr))

    # Get mask and mask images
    mask = np.asarray(props['mask'])

    av_image_masked = np.copy(av_image)
    #av_image_masked[(mask == 1) & (region_mask == 1)] = np.nan
    av_image_masked[mask == 1] = np.nan

    av_error_masked = np.copy(av_image_error)
    #av_image_masked[(mask == 1) & (region_mask == 1)] = np.nan
    av_error_masked[mask == 1] = np.nan

    av_model_masked = np.copy(av_model)
    #av_model_masked[(mask == 1) & (region_mask == 1)] = np.nan
    av_model_masked[mask == 1] = np.nan

    indices = ((np.isnan(av_model_masked)) & \
               (np.isnan(av_image_masked)) & \
               (np.isnan(av_image_error)))

    print('\nTotal number of pixels after masking = ' + str(props['npix']))

    # import matplotlib.pyplot as plt
    # av_plot_data = np.copy(av_image)
    # av_plot_data[~indices] = np.nan
    # plt.imshow(av_plot_data, origin='lower')
    # plt.show()

    # Calc chi^2
    chisq = np.sum((av_image[~mask] - av_model[~mask])**2 / \
            av_image_error[~mask]**2) / props['npix']
    props['chisq'] = chisq

    # Create HI spectrum
    hi_cube[hi_cube != hi_cube] = 0
    hi_cube[:, mask==1] = 0
    hi_spectrum = np.mean(hi_cube, axis=(1,2))

    # Derive CO spectrum
    co_data[:, region_mask == 1] = 0
    co_data[np.isnan(co_data)] = 0
    co_spectrum = np.mean(co_data, axis=(1,2))

    # Plot
    figure_types = ['png', 'pdf']
    for figure_type in figure_types:
        if region is None:
            if vel_range_type == 'single':
                filename = 'single_vel_range/taurus_av_model_map_' + \
                    'dgr{0:.3f}_'.format(dgr) + \
                    '{0:.1f}to{1:.1f}kms'.format(vel_range[0], vel_range[1]) + \
                    '_' + av_data_type + '.%s' % figure_type
            elif vel_range_type == 'multiple':
                filename = 'multiple_vel_range/taurus_av_model_map_' + \
                           'dgr{0:.3f}'.format(dgr)
                for i in xrange(0, vel_range.shape[0]):
                    filename += '_{0:.1f}to{1:.1f}kms'.format(vel_range[i, 0],
                                                              vel_range[i, 1])
                filename += '.%s' % figure_type
        else:
            filename = 'taurus_av_model_map_region{0:.0f}'.format(region) + \
                       '.{0:s}'.format(figure_type)

        print('\nSaving Av model image to \n' + filename)

        plot_av_model(av_image=av_image_masked,
                      av_model=av_model_masked,
                      header=av_header,
                      results=props,
                      hi_velocity_axis=velocity_axis,
                      vel_range=vel_range,
                      hi_spectrum=hi_spectrum,
                      #hi_limits=[-15, 25, -1, 10],
                      hi_limits=[-15, 25, None, None],
                      co_spectrum=co_spectrum,
                      co_velocity_axis=co_velocity_axis,
                      limits=props['plot_limit']['pixel'],
                      savedir=figure_dir + 'maps/av_models/',
                      filename=filename,
                      show=False)

        plot_avmod_vs_av((av_model_masked,),
                (av_image_masked,),
                av_errors=(av_error_masked,),
                #limits=[10**-1, 10**1.9, 10**0, 10**1.7],
                limits=[0,1.5,0,1.5],
                savedir=figure_dir + 'av/',
                gridsize=(10,10),
                #scale=('log', 'log'),
                #scale=('linear', 'linear'),
                filename='taurus_avmod_vs_av.%s' % figure_type,
                show = False,
                std=0.22,
                )

        plot_power_spectrum(av_image_masked - av_model_masked,
            filename_prefix='taurus_av_resid_power_spectrum_' + \
                            '{0:s}'.format(av_data_type),
            filename_suffix='.{0:s}'.format(figure_type),
            savedir=figure_dir + 'power_spectra/',
            show=False)

if __name__ == '__main__':
    dgrs = np.arange(0.05, 0.4, 0.025)
    vel_widths = (4, 6, 8, 10, 12, 18, 34, 50)
    vel_widths = (50,)
    vel_center = 6.15

    regions = (1,2,3)

    #

    # Use Planck dust Av map or Kainulainen 2009 optical extinction Av map?

    main(dgr=None, vel_range=None, vel_range_type='single', region=None,
            av_data_type='planck')
    main(dgr=None, vel_range=None, vel_range_type='single', region=None,
            av_data_type='k09')

    '''
    for region in regions:
        main(region=region, vel_range=None)
    '''

    '''
    # Single velocity ranges
    for i in xrange(0, len(dgrs)):
        for j in xrange(0, len(vel_widths)):
            vel_range = (vel_center - vel_widths[j]/2.0,
                         vel_center + vel_widths[j]/2.0)
            vel_range = np.asarray(vel_range)
            main(dgr=dgrs[i], vel_range=vel_range,
                    vel_range_type='single')

    # Multiple velocity ranges
    for i in xrange(0, len(dgrs)):
        for j in xrange(0, len(vel_widths)):
            vel_range = ((-10, vel_center - vel_widths[j]/2.0),
                         (vel_center + vel_widths[j]/2.0, 20),)
            vel_range = np.asarray(vel_range)
            main(dgr=dgrs[i], vel_range=vel_range,
                    vel_range_type='multiple')


    '''



