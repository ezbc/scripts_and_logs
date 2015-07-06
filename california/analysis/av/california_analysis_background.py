#!/usr/bin/python

''' Calculates the N(HI) map for california

'''

import numpy as np
import warnings
warnings.filterwarnings('ignore')

#import matplotlib
#matplotlib.use('Agg')

''' Plotting Functions
'''

def plot_nhi_image(nhi_image=None, header=None, contour_image=None,
        av_data=None,
        cores=None, title=None, limits=None,
        contours=None, boxes=False, savedir='./', filename=None, show=True):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
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

    if av_data is not None:
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

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)

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
    if av_data is not None:
        # create axes
        ax = imagegrid[1]
        # show the image
        im = ax.imshow(av_data,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension [J2000]',)
        ax.set_ylabel('Declination [J2000]',)

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

def plot_av_model(av_data=None, header=None, contour_image=None,
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
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon

    # Set up plot aesthetics
    plt.clf()
    plt.close()
    plt.rcdefaults()
    cmap = plt.cm.gnuplot
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 4)]
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 9
    if plot_residuals:
        if hi_spectrum is not None:
            figsize = (7, 8.5)
        else:
            figsize = (3.3, 7)
    else:
    	figsize = (13,10)

    font_scale = 9
    line_weight = 600
    font_weight = 600
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'axes.weight': line_weight,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'xtick.weight': line_weight,
              'ytick.labelsize': font_scale,
              'ytick.weight': line_weight,
              'font.weight': font_weight,
              'axes.labelweight': font_weight,
              'text.usetex': True,
              #'font.family': 'sans-serif',
              'figure.figsize': figsize,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    #==========================================================================
    # Av maps
    #==========================================================================

    if plot_residuals:
        nrows_ncols=(3,1)
        ngrids=3
        if hi_spectrum is not None:
            subplots = 210
        else:
            subplots = 110
    else:
        nrows_ncols=(1,2)
        ngrids=2

    imagegrid = ImageGrid(fig, subplots + 1,
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode='each',
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='5%',
                 axes_pad=0.0,
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
    cmap = cm.gnuplot # colormap
    cmap.set_bad(color='w')
    # show the image
    vmax = np.max((np.max(av_data[av_data == av_data]),
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

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)

    #ax.set_title(r'Model $A_V$')
    ax.annotate(r'Model $A_V$ [mag]',
                xytext=(0.05, 0.9),
                xy=(0.05, 0.9),
                textcoords='axes fraction',
                xycoords='axes fraction',
                size=font_scale,
                color='k',
                bbox=dict(boxstyle='square',
                          facecolor='w',
                          alpha=1),
                horizontalalignment='left',
                verticalalignment='top',
                )

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
    if av_data is not None:
        # create axes
        ax = imagegrid[1]
        # show the image
        im = ax.imshow(av_data,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                vmin=0,
                vmax=vmax,
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension [J2000]',)
        ax.set_ylabel('Declination [J2000]',)

        #ax.set_title(r'Observed $A_V$')
        ax.annotate(r'Observed $A_V$ [mag]',
                    xytext=(0.05, 0.9),
                    xy=(0.05, 0.9),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=font_scale,
                    color='k',
                    bbox=dict(boxstyle='square',
                              facecolor='w',
                              alpha=1),
                    horizontalalignment='left',
                    verticalalignment='top',
                    )
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
        ax = imagegrid[2]
        #cmap = cm.Greys # colormap
        # show the image
        vmax = np.max((np.max(av_data[av_data == av_data]),
                       np.max(av_model[av_model == av_model])))
        vmax = 1.4
        im = ax.imshow(av_data - av_model,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                vmin=-0.35,
                vmax=0.35,
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension [J2000]',)
        ax.set_ylabel('Declination [J2000]',)

        #ax.set_title(r'Residual $A_V$')
        ax.annotate(r'Residual $A_V$ [mag]',
                    xytext=(0.05, 0.9),
                    xy=(0.05, 0.9),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=font_scale,
                    color='k',
                    bbox=dict(boxstyle='square',
                              facecolor='w',
                              alpha=1),
                    horizontalalignment='left',
                    verticalalignment='top',
                    )
        # colorbar
        cb = ax.cax.colorbar(im)
        #cmap.set_bad(color='w')
        # plot limits
        if limits is not None:
            ax.set_xlim(limits[0],limits[2])
            ax.set_ylim(limits[1],limits[3])


    # =========================================================================
    # HI Spectrum
    # =========================================================================

    if hi_spectrum is not None:

        print('\nIncluding HI spectrum in Av model plot...')

        vel_range = np.asarray(results['hi_velocity_range'][:2])

        # create axes
        ax = fig.add_subplot(subplots + 1)

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
        ax.axvline(results['vel_center'], color='k', linestyle='--', )

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
                text += '\n'
            if co_thres is not None:
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
        plt.savefig(savedir + filename, bbox_inches='tight', dpi=400)
    if show:
        plt.show()

def plot_avmod_vs_av(avmod_images, av_datas, av_error_datas=None, limits=None,
        fit=True, savedir='./', filename=None, show=True,
        scale=('linear','linear'), title = '', gridsize=(100,100), std=None):

    # import external modules
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from matplotlib import cm
    from astroML.plotting import scatter_contour

    n = int(np.ceil(len(av_datas)**0.5))
    if n**2 - n > len(av_datas):
        nrows = n - 1
        ncols = n
        y_scaling = 1.0 - 1.0 / n
    else:
        nrows, ncols = n, n
        y_scaling = 1.0

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

    # Create figure instance
    fig = plt.figure(figsize=(3.6, 3.6))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=len(av_datas),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 #cbar_mode='single',
                 cbar_pad=0.1,
                 cbar_size=0.2,
                 )

    # Cycle through lists
    for i in xrange(len(av_datas)):
        av = av_datas[i]
        avmod = avmod_images[i]
        av_error_data = av_error_datas[i]

        # Drop the NaNs from the images
        if type(av_error_data) is float:
            indices = np.where((av == av) &\
                               (avmod == avmod)
                               )

        if type(av_error_data) is np.ndarray or \
                type(av_error_data) is np.ma.core.MaskedArray or \
                type(avmod_error) is np.ndarray or \
                type(avmod_error) is np.ma.core.MaskedArray:
            indices = np.where((av == av) &\
                               (avmod == avmod) &\
                               (av_error_data == av_error_data)
                               )

        av_nonans = av[indices]
        avmod_nonans = avmod[indices]

        # Fix error data types
        if type(av_error_data) is np.ndarray:
            av_error_data_nonans = av_error_data[indices]
        else:
            av_error_data_nonans = np.array(av_error_data[indices])

        # Create plot
        ax = imagegrid[i]

        if 0:
            image = ax.errorbar(av_nonans.ravel(),
                    avmod_nonans.ravel(),
                    xerr=(av_error_data_nonans.ravel()),
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
        if 1:
            l1 = scatter_contour(av_nonans.ravel(),
                                 avmod_nonans.ravel(),
                                 threshold=2,
                                 log_counts=0,
                                 levels=5,
                                 ax=ax,
                                 histogram2d_args=dict(bins=30,
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



        # Plot sensitivies
        av_limit = np.median(av_error_datas[0])
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
        ax.set_xlabel(r'$A_V$ Data [mag]')
        ax.set_ylabel(r'$A_V$ Model [mag]')
        #ax.set_title(core_names[i])

    if filename is not None:
        plt.savefig(savedir + filename)
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

def plot_av_images(av_image=None, header=None, contour_image=None,
        av_image_backsub=None, av_background=None, cores=None, props=None,
        regions=None, title=None, limits=None, contours=None, boxes=False,
        filename=None, show=True, hi_vlimits=None, av_vlimits=None,
        av_back_vlimits=None,):

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

    # Color map
    cmap = plt.cm.gnuplot

    params = {
              'figure.figsize': (7.3, 9),
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(3,1),
                 ngrids=3,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0.3,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # ------------------
    # 2MASS image
    # ------------------
    # create axes
    ax = imagegrid[0]

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
    ax.set_title('Original $A_V$', fontsize=9)

    ax.locator_params(nbins=6)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')

    # Write label to colorbar
    cb.set_label_text(r'$A_V$ [mag]',)

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

    if regions is not None:
        for region in regions:
            vertices = np.copy(regions[region]['poly_verts']['pixel'])
            rect = ax.add_patch(Polygon(
                    vertices[:, ::-1],
                    facecolor='none',
                    edgecolor='w'))

    # ------------------
    # Planck image
    # ------------------
    # create axes
    ax = imagegrid[1]
    # show the image
    im = ax.imshow(av_image_backsub,
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
    ax.set_title('Background-subtracted $A_V$', fontsize=9)

    ax.locator_params(nbins=6)

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
    cb.set_label_text(r'$A_V$ [mag]',)

    # ------------------
    # Division image
    # ------------------
    # create axes
    ax = imagegrid[2]

    # show the image
    if type(av_background) is not np.ndarray:
        av_background = np.ones(av_image.shape) * av_background

    im = ax.imshow(av_background,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=av_back_vlimits[0],
            vmax=av_back_vlimits[1],
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)
    ax.set_title('Fitted Background $A_V$', fontsize=9)

    ax.locator_params(nbins=6)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Write label to colorbar
    cb.set_label_text(r'$A_V$ [mag]',)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()

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

def load_ds9_region(props, filename=None, header=None, key='regions'):

    import pyregion as pyr

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    props[key] = {}

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        if region.comment is not None:
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

            props[key][region_name] = {}
            props[key][region_name]['poly_verts'] = {}
            props[key][region_name]['poly_verts']['wcs'] = poly_verts
            props[key][region_name]['poly_verts']['pixel'] = poly_verts_pix

    return props

def fit_background(av_data, background_mask=None, background_dim=1):

    from scipy.interpolate import interp2d
    from scipy.interpolate import SmoothBivariateSpline as spline

    if background_mask is None:
        background_mask = np.zeros(av_data.shape)

    if background_dim == 1:
        background = np.nanmean(av_data[~background_mask])

    if background_dim == 2:
        #av_data = np.ma.array(av_data, mask=background_mask)

        loc = np.where(~background_mask)

        x = loc[0]
        y = loc[1]

        z = av_data[~background_mask]

        #print av_data.size, z.shape
        assert z.shape == x.shape

        bbox = [0, av_data.shape[0], 0, av_data.shape[1]]

        result_interp = spline(x, y, z, bbox=bbox, kx=1, ky=1)

        x_grid, y_grid = np.where(av_data)

        background_flat = result_interp.ev(x_grid, y_grid)

        background = np.reshape(background_flat, av_data.shape)

    return background

'''
The main script
'''

def main(dgr=None, vel_range=None, vel_range_type='single', region=None,
        av_data_type='planck', use_binned_images=False, background_dim=1):
    ''' Executes script.

    Parameters
    ----------
    dgr : float
        If None, pulls best-fit value from properties.
    vel_range : tuple
        If None, pulls best-fit value from properties.
    '''

    # import external modules
    import pyfits as fits
    import numpy as np
    from mycoords import make_velocity_axis
    import mygeometry as myg
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    import json

    # Script parameters
    # -----------------
    if use_binned_images:
        bin_string = '_bin'
    else:
        bin_string = ''

    # Name of noise cube
    noise_cube_filename = \
            'california_hi_galfa_cube_regrid_planckres_noise' + bin_string + \
            '.fits'

    # Name of property files results are written to
    prop_file = 'california_global_properties_' + av_data_type + '_scaled'

    # Name of property files results are written to
    background_file = 'california_background_' + av_data_type

    # Regions, regions to edit the global properties with
    region_limit = None

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/california/figures/'
    av_dir = '/d/bip3/ezbc/california/data/av/'
    hi_dir = '/d/bip3/ezbc/california/data/hi/'
    co_dir = '/d/bip3/ezbc/california/data/co/'
    core_dir = '/d/bip3/ezbc/california/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/california/data/python_output/'
    region_dir = '/d/bip3/ezbc/california/data/python_output/ds9_regions/'

    # load Planck Av and GALFA HI images, on same grid
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
                                          'california_av_planck_tau353_5arcmin.fits',
                                          header=True)

        av_error_data, av_error_data_header = fits.getdata(av_dir + \
                                    'california_av_error_planck_tau353_5arcmin.fits',
                                    header=True)

        #av_data -= 0.9 # background

    # Load global properties of cloud
    # global properties written from script
    # 'av/california_analysis_global_properties.txt'
    if region is not None:
        likelihood_filename += '_region{0:.0f}'.format(region)
        results_filename += '_region{0:.0f}'.format(region)

    print('\nReading global parameter file\n' + prop_file + '.txt')
    if 1:
        with open(property_dir + prop_file + '.txt', 'r') as f:
            props = json.load(f)

    # Load background regions from ds9
    props = load_ds9_region(props,
                            filename=region_dir + 'california_background.reg',
                            header=av_header,
                            key='background_regions')

    # Convert plot limits
    props['plot_limit'] = {}
    props['plot_limit']['wcs'] = (((4, 50, 0), (33, 0, 0)),
                                  ((4, 10, 0), (39, 0, 0)))

    props = convert_limit_coordinates(props, coords=('plot_limit',),
            header=av_header)

    # Derive relevant region
    background_mask = np.ones(av_data.shape)
    for background_region in props['background_regions']:
        background_vertices = \
          props['background_regions'][background_region]['poly_verts']['pixel']

        # block off region
        background_mask_temp = ~np.logical_not(myg.get_polygon_mask(av_data,
                                            background_vertices))

        background_mask[background_mask_temp] = 0

    background_mask = ~np.logical_not(background_mask)

    if 0:
        import matplotlib.pyplot as plt
        av_plot_data = np.copy(av_data)
        av_plot_data[background_mask] = np.nan
        plt.imshow(av_plot_data, origin='lower')
        #plt.xlim(props['plot_limit_bin']['pixel'][0:3:2])
        #plt.ylim(props['plot_limit_bin']['pixel'][1:4:2])
        plt.show()

    background = fit_background(av_data, background_mask,
            background_dim=background_dim)

    if background_dim == 1:
        print('\nBackground A_V = {0:.1f} mag'.format(background))
        props['background_1D'] = float(background)

    if background_dim == 2:
        print('\nBackground A_V is 2D')
        props['background_2D'] = background.tolist()

    print('\nWriting global parameter file\n' + prop_file + '.txt')
    with open(property_dir + background_file + '.txt', 'w') as f:
        json.dump(props, f)

    # Plot
    figure_types = ['png',]
    for figure_type in figure_types:
        filename = figure_dir + 'maps/california_av_background_maps_' + \
                   '{0:d}D.'.format(background_dim) + figure_type

        print('\nSaving maps to \n' + filename)

        plot_av_images(av_image=av_data,
                       av_image_backsub=av_data - background,
                       av_background=background,
                       header=av_header,
                       regions=props['background_regions'],
                       av_vlimits=(-1,16),
                       av_back_vlimits=(0,3),
                       limits=props['plot_limit']['pixel'],
                       filename=filename,
                       show=False)

if __name__ == '__main__':
    main(av_data_type='planck', background_dim=1)
    main(av_data_type='planck', background_dim=2)

