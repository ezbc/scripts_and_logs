#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the taurus molecular cloud.
'''

import pyfits as pf
import numpy as np


''' Plotting Functions
'''

def plot_nhi_vs_av(nhi_image, av_image,
        nhi_image_error=None, av_image_error=None, limits=None,
        savedir='./', filename=None, show=True, scale='linear',
        returnimage=False, hess_binsize=None, title='', plot_type='hexbin'):
    ''' Plots N(HI) as a function of Av for individual pixels in an N(HI) image
    and an Av image.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    # Drop the NaNs from the images
    indices = np.where((nhi_image == nhi_image) &\
                       (av_image == av_image)&\
                       (av_image > 0) &\
                       (nhi_image > -5))

    nhi_image_nonans = nhi_image[indices]
    av_image_nonans = av_image[indices]

    if type(nhi_image_error) is np.ndarray:
        nhi_image_error_nonans = nhi_image_error[indices]
    else:
        nhi_image_error_nonans = np.array(nhi_image_error[indices])

    if type(av_image_error) is np.ndarray:
        av_image_error_nonans = av_image_error[indices]
    else:
        av_image_error_nonans = av_image_error * \
                np.ones(av_image[indices].shape)

    # Create figure
    plt.clf()
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    if nhi_image_error is None:
        if plot_type is 'hexbin':
            image = ax.hexbin(av_image_nonans.ravel(),
                    nhi_image_nonans.ravel(),
                    norm=matplotlib.colors.LogNorm(),
                    mincnt=1,
                    yscale='log',
                    xscale='log')
            # Adjust color bar of density plot
            cb = plt.colorbar(image)
            cb.set_label('Counts')
        elif plot_type is 'scatter':
            image = ax.scatter(av_image_nonans.ravel(),
                    nhi_image_nonans.ravel(),
                    alpha=0.3,
                    color='k'
                    )
            ax.set_xscale('log')
            ax.set_yscale('log')
    else:
        image = ax.errorbar(av_image_nonans.ravel(),
                nhi_image_nonans.ravel(),
                xerr=(av_image_error_nonans.ravel()),
                yerr=(nhi_image_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=2
                )

        ax.set_xscale(scale)
        ax.set_yscale(scale)

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel('A$_v$ (mag)',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'N(HI) (1 $\times 10^{20}$ cm$^{-2}$)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_hisd_vs_hsd(hi_sd_image, h_sd_image,
        hi_sd_image_error=None, h_sd_image_error=None, limits=None,
        savedir='./', filename=None, show=True, scale='linear',
        returnimage=False, hess_binsize=None, title='', plot_type='hexbin'):

    ''' Plots N(HI) as a function of Av for individual pixels in an N(HI) image
    and an Av image.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    # Drop the NaNs from the images
    indices = np.where((hi_sd_image == hi_sd_image) &\
                       (h_sd_image == h_sd_image)&\
                       (h_sd_image > 0) &\
                       (hi_sd_image > -5))

    hi_sd_image_nonans = hi_sd_image[indices]
    h_sd_image_nonans = h_sd_image[indices]

    if type(hi_sd_image_error) is np.ndarray:
        hi_sd_image_error_nonans = hi_sd_image_error[indices]
    else:
        hi_sd_image_error_nonans = np.array(hi_sd_image_error[indices])

    if type(h_sd_image_error) is np.ndarray:
        h_sd_image_error_nonans = h_sd_image_error[indices]
    elif type(h_sd_image_error) is np.ma.core.MaskedArray:
        #h_sd_image_error_nonans = np.copy(h_sd_image_error[indices])
        h_sd_image_error_nonans = h_sd_image_error[indices]
    else:
        h_sd_image_error_nonans = h_sd_image_error * \
                np.ones(h_sd_image[indices].shape)

    # Create figure
    plt.clf()
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    if hi_sd_image_error is None:
        if plot_type is 'hexbin':
            image = ax.hexbin(h_sd_image_nonans.ravel(),
                    hi_sd_image_nonans.ravel(),
                    norm=matplotlib.colors.LogNorm(),
                    mincnt=1,
                    yscale='log',
                    xscale='log')
            # Adjust color bar of density plot
            cb = plt.colorbar(image)
            cb.set_label('Counts')
        elif plot_type is 'scatter':
            image = ax.scatter(h_sd_image_nonans.ravel(),
                    hi_sd_image_nonans.ravel(),
                    alpha=0.3,
                    color='k'
                    )
            ax.set_xscale('log')
            ax.set_yscale('log')
    else:
        image = ax.errorbar(h_sd_image_nonans.ravel(),
                hi_sd_image_nonans.ravel(),
                xerr=(h_sd_image_error_nonans.ravel()),
                yerr=(hi_sd_image_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=2
                )

        ax.set_xscale(scale)
        ax.set_yscale(scale)

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{H2}$ (M$_\odot$ / pc$^2$)',
              size = 'small',
              family='serif')
    ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_sd_vs_av(sd_image, av_image,
        sd_image_error=None, av_image_error=None, limits=None,
        savedir='./', filename=None, show=True, scale='linear',
        returnimage=False, hess_binsize=None, title='', plot_type='hexbin'):
    ''' Plots N(HI) as a function of Av for individual pixels in an N(HI) image
    and an Av image.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    # Drop the NaNs from the images
    indices = np.where((sd_image == sd_image) &\
                       (av_image == av_image)&\
                       (av_image > 0) &\
                       (sd_image > -5))

    sd_image_nonans = sd_image[indices]
    av_image_nonans = av_image[indices]

    if type(sd_image_error) is np.ndarray:
        sd_image_error_nonans = sd_image_error[indices]
    else:
        sd_image_error_nonans = np.array(sd_image_error[indices])

    if type(av_image_error) is np.ndarray:
        av_image_error_nonans = av_image_error[indices]
    else:
        av_image_error_nonans = av_image_error * \
                np.ones(av_image[indices].shape)

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
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (6, 6),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)


    # Create figure
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if sd_image_error is None:
        if plot_type is 'hexbin':
            image = ax.hexbin(av_image_nonans.ravel(),
                    sd_image_nonans.ravel(),
                    norm=matplotlib.colors.LogNorm(),
                    mincnt=1,
                    yscale='log',
                    xscale='log')
            # Adjust color bar of density plot
            cb = plt.colorbar(image)
            cb.set_label('Counts')
        elif plot_type is 'scatter':
            image = ax.scatter(av_image_nonans.ravel(),
                    sd_image_nonans.ravel(),
                    alpha=0.3,
                    color='k'
                    )
            ax.set_xscale('log')
            ax.set_yscale('log')
    else:
        image = ax.errorbar(av_image_nonans.ravel(),
                sd_image_nonans.ravel(),
                xerr=(av_image_error_nonans.ravel()),
                yerr=(sd_image_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=2
                )

        ax.set_xscale('linear')
        ax.set_yscale(scale)

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'A$_{\rm V}$ (mag)',)
    ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_sd_vs_av_grid(sd_images, av_images,
        sd_image_errors=None, av_image_errors=None, limits=None,
        savedir='./', filename=None, show=True, scale=('linear', 'linear'),
        returnimage=False, title=None, core_names=''):

    ''' Plots N(HI) as a function of Av for individual pixels in an N(HI) image
    and an Av image.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid

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
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (10, 10),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    n = int(np.ceil(len(av_images)**0.5))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(n, n),
                 ngrids=len(av_images),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    for i in xrange(len(av_images)):
        sd_image = sd_images[i]
        av_image = av_images[i]
        sd_image_error = sd_image_errors[i]
        av_image_error = av_image_errors[i]

        # Drop the NaNs from the images
        indices = np.where((sd_image == sd_image) &\
                           (av_image == av_image) &\
                           (av_image > 0) &\
                           (sd_image > -5))

        sd_image_nonans = sd_image[indices]
        av_image_nonans = av_image[indices]

        if type(sd_image_error) is np.ndarray:
            sd_image_error_nonans = sd_image_error[indices]
        else:
            sd_image_error_nonans = np.array(sd_image_error[indices])

        if type(av_image_error) is np.ndarray:
            av_image_error_nonans = av_image_error[indices]
        else:
            av_image_error_nonans = av_image_error * \
                    np.ones(av_image[indices].shape)

        # Create plot
        ax = imagegrid[i]

        image = ax.errorbar(av_image_nonans.ravel(),
                sd_image_nonans.ravel(),
                xerr=(av_image_error_nonans.ravel()),
                yerr=(sd_image_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=2
                )

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'A$_{\rm V}$ (mag)',)
        ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
        ax.set_title(core_names[i])
        ax.grid(True)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_hisd_vs_hsd_grid(hi_sd_images, h_sd_images,
        hi_sd_image_errors=None, h_sd_image_errors=None, limits=None,
        savedir='./', filename=None, show=True, scale=('linear', 'linear'),
        returnimage=False, title=None, core_names=''):
    ''' Plots N(HI) as a function of Av for individual pixels in an N(HI) image
    and an Av image.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid

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
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (10, 10),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    n = int(np.ceil(len(h_sd_images)**0.5))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(n, n),
                 ngrids=len(h_sd_images),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    for i in xrange(len(h_sd_images)):
        hi_sd_image = hi_sd_images[i]
        h_sd_image = h_sd_images[i]
        hi_sd_image_error = hi_sd_image_errors[i]
        h_sd_image_error = h_sd_image_errors[i]

        # Drop the NaNs from the images
        indices = np.where((hi_sd_image == hi_sd_image) &\
                           (hi_sd_image_error == hi_sd_image_error)&\
                           (h_sd_image == h_sd_image)&\
                           (h_sd_image_error == h_sd_image_error)&\
                           (h_sd_image > 0) &\
                           (hi_sd_image > 0))

        hi_sd_image_nonans = hi_sd_image[indices]
        h_sd_image_nonans = h_sd_image[indices]

        if type(hi_sd_image_error) is np.ndarray:
            hi_sd_image_error_nonans = hi_sd_image_error[indices]
        else:
            hi_sd_image_error_nonans = np.array(hi_sd_image_error[indices])

        if type(h_sd_image_error) is np.ndarray:
            h_sd_image_error_nonans = h_sd_image_error[indices]
        else:
            h_sd_image_error_nonans = h_sd_image_error * \
                    np.ones(h_sd_image[indices].shape)

        # Create plot
        ax = imagegrid[i]

        image = ax.errorbar(h_sd_image_nonans.ravel(),
                hi_sd_image_nonans.ravel(),
                xerr=(h_sd_image_error_nonans.ravel()),
                yerr=(hi_sd_image_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=2
                )

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'$\Sigma_{HI}$ + $\Sigma_{H_2}$ (M$_\odot$ / pc$^2$)',)
        ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
        ax.set_title(core_names[i])
        ax.grid(True)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_rh2_vs_h_grid(rh2_images, h_sd_images, rh2_error_images=None,
        h_sd_error_images=None, rh2_fits = None, h_sd_fits = None, limits =
        None, fit = True, savedir = './', filename = None, show = True, scale =
        'linear', title = '', core_names='', phi_cnm_list=None,
        phi_cnm_error_list=None, Z_list=None, Z_error_list=None,
        chisq_list=None, p_value_list=None,):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid

    n = int(np.ceil(len(rh2_images)**0.5))
    if n**2 - n > len(rh2_images):
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
              'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=len(rh2_images),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    # Cycle through lists

    if len(phi_cnm_error_list) == 0:
        phi_cnm_error_list = None
    if len(Z_error_list) == 0:
        Z_error_list = None

    for i in xrange(len(rh2_images)):
        rh2 = rh2_images[i]
        h_sd = h_sd_images[i]
        rh2_error = rh2_error_images[i]
        h_sd_error = h_sd_error_images[i]
        rh2_fit = rh2_fits[i]
        h_sd_fit = h_sd_fits[i]
        if phi_cnm_list is not None:
            phi_cnm = phi_cnm_list[i]
        if phi_cnm_error_list is not None:
            phi_cnm_error = phi_cnm_error_list[i]
        if Z_list is not None:
            Z = Z_list[i]
        if Z_error_list is not None:
            Z_error = Z_error_list[i]
        if chisq_list is not None:
            chisq = chisq_list[i]
        if p_value_list is not None:
            p_value = p_value_list[i]

        # Drop the NaNs from the images
        if type(rh2_error) is float:
            indices = np.where((rh2 == rh2) &\
                               (h_sd == h_sd)&\
                               (h_sd > 0) &\
                               (rh2 > 0))

        if type(rh2_error) is np.ndarray or \
                type(rh2_error) is np.ma.core.MaskedArray or \
                type(h_sd_error) is np.ndarray or \
                type(h_sd_error) is np.ma.core.MaskedArray:
            indices = np.where((rh2 == rh2) &\
                               (h_sd == h_sd) &\
                               (h_sd_error == h_sd_error) &\
                               (rh2_error == rh2_error) &\
                               (h_sd > 0) &\
                               (rh2 > 0))

        rh2_nonans = rh2[indices]
        h_sd_nonans = h_sd[indices]

        if type(rh2_error) is np.ndarray:
            rh2_error_nonans = rh2_error[indices]
        else:
            rh2_error_nonans = np.array(rh2_error[indices])

        if type(h_sd_error) is np.ndarray or \
                type(h_sd_error) is np.ma.core.MaskedArray:
            h_sd_error_nonans = h_sd_error[indices]
        else:
            h_sd_error_nonans = h_sd_error * \
                    np.ones(h_sd[indices].shape)

                # Create plot
        ax = imagegrid[i]

        image = ax.errorbar(h_sd_nonans.ravel(),
                rh2_nonans.ravel(),
                xerr=(h_sd_error_nonans.ravel()),
                yerr=(rh2_error_nonans.ravel()),
                alpha=0.5,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=4
                )

        if rh2_fit is not None:
            ax.plot(h_sd_fit, rh2_fit,
                    color = 'r')

        # Annotations
        anno_xpos = 0.95

        if phi_cnm_list is not None and Z_list is not None:
            if phi_cnm_error_list is None and Z_error_list is not None:
                ax.annotate(r'$\phi_{\rm CNM}$ = {0:.2f}\n'.format(phi_cnm) + \
                            r'Z = {0:.2f} Z$_\odot$'.format(Z),
                        xytext=(anno_xpos, 0.05),
                        xy=(anno_xpos, 0.05),
                        textcoords='axes fraction',
                        xycoords='axes fraction',
                        color='k',
                        bbox=dict(boxstyle='round',
                                  facecolor='w',
                                  alpha=0.5),
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        )
            else:
                ax.annotate(r'\noindent$\phi_{\rm CNM}$ =' + \
                            r' %.2f' % (phi_cnm) + \
                            r'$^{+%.2f}_{-%.2f}$ \\' % (phi_cnm_error[0],
                                                     phi_cnm_error[1]) + \
                            r'Z = %.2f' % (Z) + \
                            r'$^{+%.2f}_{-%.2f}$ Z$_\odot$' % (Z_error[0],
                                                               Z_error[1]) + \
                            r'',
                        xytext=(anno_xpos, 0.05),
                        xy=(anno_xpos, 0.05),
                        textcoords='axes fraction',
                        xycoords='axes fraction',
                        size=font_scale*3/4.0,
                        color='k',
                        bbox=dict(boxstyle='round',
                                  facecolor='w',
                                  alpha=1),
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        )
        if chisq_list is not None:
            ax.annotate(r'$\chi^2/\nu$ = %.2f' % \
                    chisq,
                    xytext=(0.48, 0.3),
                    xy=(0.48, 0.3),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    color='k'
                    )

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{H2}$ (M$_\odot$ / pc$^2$)',)
        ax.set_ylabel(r'R$_{H2}$',)
        ax.set_title(core_names[i])
        ax.grid(True)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        fig.show()

def plot_hi_vs_av_grid(hi_images, av_images, hi_error_images=None,
        av_error_images=None, limits=None, savedir='./',
        filename=None, show=True, scale='linear', title='', core_names=''):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid

    n = int(np.ceil(len(hi_images)**0.5))
    if n**2  - n > len(hi_images):
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
              'figure.figsize': (8, 8 * y_scaling),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=len(hi_images),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    # Cycle through lists
    for i in xrange(len(hi_images)):
        hi = hi_images[i]
        av = av_images[i]
        hi_error = hi_error_images[i]
        av_error = av_error_images[i]

        # Drop the NaNs from the images
        if type(hi_error) is float:
            indices = np.where((hi == hi) &\
                               (av == av)&\
                               (av > 0) &\
                               (hi > 0))

        if type(hi_error) is np.ndarray or \
                type(hi_error) is np.ma.core.MaskedArray or \
                type(av_error) is np.ndarray or \
                type(av_error) is np.ma.core.MaskedArray:
            indices = np.where((hi == hi) &\
                               (av == av) &\
                               (av_error == av_error) &\
                               (hi_error == hi_error) &\
                               (av > 0) &\
                               (hi > 0))

        hi_nonans = hi[indices]
        av_nonans = av[indices]

        if type(hi_error) is np.ndarray:
            hi_error_nonans = hi_error[indices]
        else:
            hi_error_nonans = np.array(hi_error[indices])

        if type(av_error) is np.ndarray or \
                type(av_error) is np.ma.core.MaskedArray:
            av_error_nonans = av_error[indices]
        else:
            av_error_nonans = av_error * \
                    np.ones(av[indices].shape)

        # Create plot
        ax = imagegrid[i]

        image = ax.errorbar(av_nonans.ravel(),
                hi_nonans.ravel(),
                xerr=(av_error_nonans.ravel()),
                yerr=(hi_error_nonans.ravel()),
                alpha=0.5,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=4
                )

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'$A_V$ (Mag)',)
        ax.set_ylabel(r'$\Sigma_{HI}$ ($M_\odot$ pc$^{-2}$)',)
        ax.set_title(core_names[i])
        ax.grid(True)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()

def plot_hi_vs_h_grid(hi_images, h_sd_images, hi_sd_error_images=None,
        h_sd_error_images=None, hi_fits = None, h_sd_fits = None, limits =
        None, fit = True, savedir = './', filename = None, show = True, scale =
        'linear', title = '', core_names='', phi_cnm_list=None):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid

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
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (10, 10),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    n = int(np.ceil(len(hi_images)**0.5))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(n, n),
                 ngrids=len(hi_images),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    # Cycle through lists
    for i in xrange(len(hi_images)):
        hi_sd = hi_images[i]
        h_sd = h_sd_images[i]
        hi_sd_error = hi_sd_error_images[i]
        h_sd_error = h_sd_error_images[i]
        hi_sd_fit = hi_fits[i]
        h_sd_fit = h_sd_fits[i]
        if phi_cnm_list is not None:
            phi_cnm = phi_cnm_list[i]


        # Drop the NaNs from the images
        if type(hi_sd_error) is float:
            indices = np.where((hi_sd == hi_sd) &\
                               (h_sd == h_sd)&\
                               (h_sd > 0) &\
                               (hi_sd > 0))

        if type(hi_sd_error) is np.ndarray or \
                type(hi_sd_error) is np.ma.core.MaskedArray or \
                type(h_sd_error) is np.ndarray or \
                type(h_sd_error) is np.ma.core.MaskedArray:
            indices = np.where((hi_sd == hi_sd) &\
                               (h_sd == h_sd) &\
                               (h_sd_error == h_sd_error) &\
                               (hi_sd_error == hi_sd_error) &\
                               (h_sd > 0) &\
                               (hi_sd > 0))

        hi_sd_nonans = hi_sd[indices]
        h_sd_nonans = h_sd[indices]

        if type(hi_sd_error) is np.ndarray:
            hi_sd_error_nonans = hi_sd_error[indices]
        else:
            hi_sd_error_nonans = np.array(hi_sd_error[indices])

        if type(h_sd_error) is np.ndarray or \
                type(h_sd_error) is np.ma.core.MaskedArray:
            h_sd_error_nonans = h_sd_error[indices]
        else:
            h_sd_error_nonans = h_sd_error * \
                    np.ones(h_sd[indices].shape)

                # Create plot
        ax = imagegrid[i]

        image = ax.errorbar(h_sd_nonans.ravel(),
                hi_sd_nonans.ravel(),
                xerr=(h_sd_error_nonans.ravel()),
                yerr=(hi_sd_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=3
                )
        if hi_sd_fit is not None:
            ax.plot(h_sd_fit, hi_sd_fit,
                    color = 'r')

        if phi_cnm_list is not None:
            ax.annotate(r'$\phi$ = %.2f' % phi_cnm,
                    xytext=(0.6, 0.1),
                    xy=(0.6, 0.1),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    color='k'
                    )

        ax.set_xscale(scale[0])
        ax.set_yscale(scale[1])

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{H2}$ (M$_\odot$ / pc$^2$)',)
        ax.set_ylabel(r'$\Sigma_{HI}$',)
        ax.set_title(core_names[i])
        ax.grid(True)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()

def plot_co_spectrum_grid(vel_axis, co_spectrum_list,
        vel_range_list=None,
        vel_range_hiav_list=None,
        limits = None, savedir = './', filename = None, show = True,
        scale = 'linear', title = '', core_names='',):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid

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
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (10, 10),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    n = int(np.ceil(len(co_spectrum_list)**0.5))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(n, n),
                 ngrids=len(co_spectrum_list),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    # Cycle through lists
    for i in xrange(len(co_spectrum_list)):
        co_spectrum = co_spectrum_list[i]
        try:
            hi_velocity_range = vel_range_list[i]
        except IndexError:
            hi_velocity_range = None
        try:
            hi_velocity_range_corr = vel_range_hiav_list[i]
        except IndexError:
            hi_velocity_range_corr = None

        # Create plot
        ax = imagegrid[i]

        image = ax.plot(vel_axis, co_spectrum,
                color='k',
                marker=None,
                drawstyle='steps-mid',
                markersize=3,
                )

        if hi_velocity_range is not None:
            ax.axvspan(hi_velocity_range[0], hi_velocity_range[1],
                    color = 'r', alpha=0.2, label=r'$^{12}$CO Spectrum Width')
        if hi_velocity_range_corr is not None:
            ax.axvspan(hi_velocity_range_corr[0], hi_velocity_range_corr[1],
                    color = 'b', alpha=0.2, label=r'HI/A$_V$ Correlation')

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel('Velocity (km/s)')
        ax.set_ylabel(r'<I$_{\rm CO}$> (K)',)
        ax.set_title(core_names[i])
        ax.grid(True)

    if vel_range_list is not None or vel_range_hiav_list is not None:
        # Single legend
        ax.legend(bbox_to_anchor=(3.1, 0.2),
                loc='lower right',
                borderaxespad=0.)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()

''' Calculations
'''

def correlate_hi_av(hi_cube=None, hi_velocity_axis=None,
        hi_noise_cube=None, av_image=None, av_image_error=None,
        vel_centers=None, vel_widths=None, return_correlations=True):

    '''
            hi_vel_range, av_correlations = correlate_hi_av(hi_cube=hi_data,
                    hi_velocity_axis=velocity_axis,
                    hi_noise_cube=noise_cube,
                    av_image=av_data_planck,
                    av_image_error=av_error_data_planck,
                    vel_centers=vel_centers,
                    vel_widths=vel_widths,
                    return_correlations=True)

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

    # calculate the velocity ranges given a set of centers and widths
    velocity_ranges = np.zeros(shape=[len(vel_centers) * \
            len(vel_widths),2])
    count = 0
    for i, center in enumerate(vel_centers):
        for j, width in enumerate(vel_widths):
            velocity_ranges[count,0] = center - width/2.
            velocity_ranges[count,1] = center + width/2.
            count += 1

    # calculate the correlation coefficient for each velocity range
    correlations = np.zeros(velocity_ranges.shape[0])
    pvalues = np.zeros(velocity_ranges.shape[0])

    for i, velocity_range in enumerate(velocity_ranges):
        nhi_image_temp, nhi_image_error = calculate_nhi(cube=hi_cube,
                velocity_axis=hi_velocity_axis,
                velocity_range=velocity_range,
                noise_cube=hi_noise_cube)

        #nhi_image = np.ma.array(nhi_image_temp,
        #                        mask=np.isnan(nhi_image_temp))

        # Select pixels with Av > 1.0 mag and Av_SNR > 5.0.
        # Av > 1.0 mag is used to avoid too low Av.
        # 1.0 mag corresponds to SNR = 1 / 0.2 ~ 5
        # (see Table 2 of Ridge et al. 2006).
        indices = np.where((nhi_image_temp == nhi_image_temp) & \
                (av_image == av_image))

        nhi_image_corr = nhi_image_temp[indices]
        av_image_corr = av_image[indices] #- 0.8 # subtract background of 0.8
        # Use Pearson's correlation test to compare images
        correlations[i] = pearsonr(nhi_image_corr.ravel(),
                av_image_corr.ravel())[0]

        # Shows progress each 10%
        #total = float(velocity_ranges.shape[0])
        #abs_step = int((total * 1)/100) or 1
        #if i and not i % abs_step:
        #    print "{0:.2%} processed".format(i/total)

    # Find best-correlating velocity range
    best_corr = correlations.max()
    best_corr_index = np.where(correlations == best_corr)
    best_corr_vel_range = velocity_ranges[best_corr_index][0]

    if not return_correlations:
        return best_corr_vel_range, best_corr
    else:
        return best_corr_vel_range, best_corr, correlations

def select_hi_vel_range(co_data, co_header, flux_threshold=0.80,
        width_scale=1.):

    from mycoords import make_velocity_axis

    vel_axis = make_velocity_axis(co_header)

    # create average spectrum
    n_pix = len(co_data[co_data == co_data])
    co_data_nonans = np.ma.array(co_data, mask=(co_data != co_data))
    spectrum = np.sum(co_data_nonans, axis=(1)) / n_pix

    # determine velocity where 25% of the flux
    peak_vel = np.interp(spectrum.max(), spectrum, vel_axis,)

    peak_vel = vel_axis[np.argmin(np.abs(spectrum - spectrum.max()))]

    max_vel_hw = np.min((peak_vel - vel_axis[0], vel_axis[1] - peak_vel))

    vel_hw = 0. # velocity half-width
    vel_delta = vel_axis[1] - vel_axis[0]

    found_threshold = False

    total_sum = np.sum(spectrum)

    # Start with small region, grow the region until integrated
    # flux in region is less than flux_threshold * spectrum total sum
    while vel_hw < max_vel_hw or not found_threshold:
        low_vel = peak_vel - vel_hw
        high_vel = peak_vel + vel_hw

        spectrum_crop = spectrum[(vel_axis > low_vel) & (vel_axis < high_vel)]

        partial_sum = np.sum(spectrum_crop)

        if partial_sum >= flux_threshold * total_sum:
            found_threshold = True
        else:
            vel_hw += vel_delta

    return peak_vel - width_scale * vel_hw, peak_vel + width_scale * vel_hw

def derive_images(hi_cube=None, hi_velocity_axis=None, hi_noise_cube=None,
        hi_vel_range=None, hi_header=None, dgr=None, av_image=None,
        av_image_error=None, sub_image_indices=None, nhi_error=None):

    '''

    Derives N(HI), Sigma_HI, Sigma_H, Sigma_H2 and RH2, plus errors on each
    image.

    Parameters
    ----------
    nhi_error : float, optional
        If not None, then this error is used in the N(HI) calculation instead
        of the cube noise.

    '''

    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error

    # Calculate N(HI) and the error
    nhi_image, nhi_image_error = calculate_nhi(cube=hi_cube,
            velocity_axis=hi_velocity_axis,
            noise_cube=hi_noise_cube,
            velocity_range=hi_vel_range,
            return_nhi_error=True,
            header=hi_header)

    if nhi_error is not None:
        nhi_image_error = nhi_error

    # mask the image for NaNs
    #nhi_image = np.ma.array(nhi_image,
    #        mask=(nhi_image != nhi_image))
    #nhi_image_error = np.ma.array(nhi_image_error,
    #        mask=(nhi_image_error != nhi_image_error))

    # calculate N(H2) maps
    nh2_image = calculate_nh2(nhi_image = nhi_image,
            av_image = av_image, dgr = dgr)

    nh2_image_error = \
        calculate_nh2_error(nhi_image_error=nhi_image_error,
            av_image_error=av_image_error, dgr = dgr)

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi_image, sd_factor=1/1.25)
    hi_sd_image_error = calculate_sd(nhi_image_error, sd_factor=1/1.25)

    h2_sd_image = calculate_sd(nh2_image, sd_factor=1/0.625)
    h2_sd_image_error = calculate_sd(nh2_image_error,
            sd_factor=1/0.625)

    h_sd_image = hi_sd_image + h2_sd_image
    h_sd_image_error = (hi_sd_image_error**2 + \
            h2_sd_image_error**2)**0.5

    # h2 surface density
    #h2_sd_image = h_sd_image - hi_sd_image
    #h2_sd_image_error=(h_sd_image_error**2-hi_sd_image_error**2)**0.5

    # Write ratio between H2 and HI
    rh2_image = h2_sd_image / hi_sd_image
    rh2_image_error = rh2_image * (hi_sd_image_error**2 / \
            hi_sd_image**2 + h2_sd_image_error**2 / \
            h2_sd_image**2)**0.5

    if sub_image_indices is not None:
        av_image = av_image[sub_image_indices]
        av_image_error = av_image_error[sub_image_indices]
        hi_sd_image = hi_sd_image[sub_image_indices]
        hi_sd_image_error = hi_sd_image_error[sub_image_indices]
        h_sd_image = h_sd_image[sub_image_indices]
        h_sd_image_error = h_sd_image_error[sub_image_indices]
        rh2_image = rh2_image[sub_image_indices]
        rh2_image_error = rh2_image_error[sub_image_indices]
        nhi_image = nhi_image[sub_image_indices]
        nhi_image_error = nhi_image_error[sub_image_indices]

    images = {'av' : av_image,
              'av_error' : av_image_error,
              'hi_sd' : hi_sd_image,
              'hi_sd_error' : hi_sd_image_error,
              'h_sd' : h_sd_image,
              'h_sd_error' : h_sd_image_error,
              'rh2' : rh2_image,
              'rh2_error' : rh2_image_error,
              'nhi' : nhi_image,
              'nhi_error' : nhi_image_error
              }

    return images

def run_analysis(hi_cube=None, hi_noise_cube=None, hi_velocity_axis=None,
        hi_header=None, dgr=None, dgr_error=None, av_image=None,
        av_image_error=None, hi_vel_range=None, hi_vel_range_error=None,
        N_runs=1, verbose=False, guesses=(10.0,1.0),
        parameter_vary=[True,True], core_dict=None, results_figure_name=None,
        error_method='bootstrap', alpha=0.05):

    '''

    Runs Monte Carlo to vary measure the error induced in RH2 from the HI
    integration velocity range. If N_runs = 1 then a simple fit without adding
    any random errors is performed.

    Parameters
    ----------
    guesses : tuple, list
        Initial guesses for phi_cnm and Z.
    paramter_vary : tuple, bool
        Vary phi_cnm and Z?
    hi_vel_range : tuple
        Tuple of floats, lower and upper bound to HI velocity range.
    hi_vel_range_error : float
        1 sigma Gaussian error on hi_vel_range limits. If error_method ==
        'gaussian', and hi_vel_range_error is None then the confidence limit
        from the HI velocity range error in the core dict will be used.
    N_runs : int
        Number of runs in the Monte Carlo simulation
    core_dict : dict
        Dictionary of core parameters.
    results_figure_name : str
        Base name for results of bootstrapping error analysis. '.png' will be
        added.
    error_method : str
        Method with which to vary the hi_velocity_range. Options are
        'bootstrap' and 'gaussian' and 'threshold'.
    alpha : float
        Significance level of confidence interval.

    Returns
    -------
    images : dict
        Dictionary of output images.
    params : dict
        Parameters in Krumholz fit. If N_runs > 1, then errors on the
        parameters are returned as well.


    '''

    from numpy.random import normal
    from scipy.stats import rv_discrete
    import mystats
    import matplotlib.pyplot as plt
    from scikits.bootstrap import ci

    if N_runs < 1:
        raise ValueError('N_runs must be >= 1')

    verbose = False

    # Results of monte carlo will be stored here
    phi_cnm_fits = np.empty((3))
    Z_fits = np.empty((3))
    hi_vel_range_list = np.empty((3, 2))

    # Get smallest and largest HI velocity intervals
    hi_vel_range_list[0, 0] = hi_vel_range[0] - hi_vel_range_error[0]
    hi_vel_range_list[0, 1] = hi_vel_range[1] + hi_vel_range_error[1]
    hi_vel_range_list[1, 0] = hi_vel_range[0]
    hi_vel_range_list[1, 1] = hi_vel_range[1]
    hi_vel_range_list[2, 0] = hi_vel_range[0] + hi_vel_range_error[0]
    hi_vel_range_list[2, 1] = hi_vel_range[1] - hi_vel_range_error[1]

    print 'hivel ranges', hi_vel_range_list

    # Run error analysis
    # -------------------
    for i in xrange(3):

        hi_vel_range = hi_vel_range_list[i]

        # Derive new images
        images = derive_images(hi_cube=hi_cube,
                               hi_velocity_axis=hi_velocity_axis,
                               hi_noise_cube=hi_noise_cube,
                               hi_vel_range=hi_vel_range,
                               hi_header=hi_header,
                               dgr=dgr,
                               av_image=av_image,
                               av_image_error=av_image_error,
                               )

        # Fit R_H2
        #---------
        # Unravel images to single dimension
        rh2_ravel = images['rh2'].ravel()
        rh2_error_ravel = images['rh2_error'].ravel()
        h_sd_ravel = images['h_sd'].ravel()
        h_sd_error_ravel = images['h_sd_error'].ravel()

        # write indices for only ratios > 0
        indices = np.where(rh2_ravel > 0)

        rh2_ravel = rh2_ravel[indices]
        rh2_error_ravel = rh2_error_ravel[indices]
        h_sd_ravel = h_sd_ravel[indices]
        h_sd_error_ravel = h_sd_error_ravel[indices]

        # Fit to krumholz model, init guess of phi_CNM = 10
        phi_cnm, Z = fit_krumholz(h_sd_ravel,
                                  rh2_ravel,
                                  guesses=guesses, # phi_cnm, Z
                                  vary=parameter_vary,
                                  rh2_error=rh2_error_ravel,
                                  verbose=verbose)

        phi_cnm_fits[i] = phi_cnm
        Z_fits[i] = Z

        if verbose:
            print('phi = %.2f' % phi_cnm)
            print('Z = %.2f' % Z)

        # see eq 6 of krumholz+09
        # phi_cnm is the number density of the CNM over the minimum number
        # density required for pressure balance
        # the lower phi_cnm values than for taurus mean that taurus
        # has a more diffuse CNM

        # By fitting the model to the observed R_H2 vs total H, you basically
        # constrained psi in Equation (35) of Krumholz+09.  This means that
        # you can calculate f_H2 for a given total hydrogen surface density.
        # In this case, H2 surface density = f_H2 * total hydrogen surface
        # density HI surface density = (1 - f_HI) * total hydrogen surface
        # density

    # Derive images
    # -------------
    images = derive_images(hi_cube=hi_cube,
                           hi_velocity_axis=hi_velocity_axis,
                           hi_noise_cube=hi_noise_cube,
                           hi_vel_range=hi_vel_range_list[1],
                           hi_header=hi_header,
                           dgr=dgr,
                           av_image=av_image,
                           av_image_error=av_image_error,
                           )

    np.sort(phi_cnm_fits)
    np.sort(Z_fits)

    phi_cnm_confint = (phi_cnm_fits[1],
                       phi_cnm_fits[2] - phi_cnm_fits[1],
                       phi_cnm_fits[1] - phi_cnm_fits[0],
                       )
    Z_confint = (Z_fits[1],
                 Z_fits[1] - Z_fits[0],
                 Z_fits[2] - Z_fits[1],
                 )

    phi_cnm = phi_cnm_confint[0]
    phi_cnm_error = phi_cnm_confint[1:]
    Z = Z_confint[0]
    Z_error = Z_confint[1:]

    # Print results
    print('Results:')
    print('phi_cnm = {0:.2f} +{1:.2f}/-{2:.2f}'.format(phi_cnm_confint[0],
                                                      phi_cnm_confint[1],
                                                      phi_cnm_confint[2]))
    print('Z = {0:.2f} +{1:.2f}/-{2:.2f}'.format(Z_confint[0],
                                                  Z_confint[1],
                                                  Z_confint[2]))

    print('HI velocity range:' + \
          '{0:.1f} to {1:.1f} km/s'.format(hi_vel_range[0],
                                           hi_vel_range[1],))

    return images, (phi_cnm, Z, phi_cnm_error, Z_error)

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
    import numpy
    from scipy.integrate import simps as integrate

    # Step for lowering threshold
    step = (np.max(y) - np.median(y)) / 100.0

    # initial threshold
    threshold = np.max(y) - step
    threshold_area = 0.0

    # area under whole function
    area = integrate(y, x)

    print area
    print y

    if type(area) != numpy.float64 or np.isnan(area):
        raise ValueError('Integration of y and x failed. Check types.')

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

    x_peak = x[y == y.max()][0]
    print x_peak, x[bounds_indices[0]]
    low_error, up_error = x_peak - x[bounds_indices[0]], \
                          x[bounds_indices[1]] - x_peak

    return (x_peak, low_error, up_error)

''' Fitting Functions
'''

def calc_krumholz(params=[10.0, 1.0], h_sd_extent=(0.001, 500),
        return_fractions=True):

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
    h_sd_extent.append(1e4)
    h_sd = np.linspace(h_sd_extent[0], h_sd_extent[1], h_sd_extent[2])

    if not return_fractions:
        rh2_fit = k09.calc_rh2(h_sd, params)
    elif return_fractions:
        rh2_fit, f_H2, f_HI = k09.calc_rh2(h_sd,
                                           *params,
                                           return_fractions=True)

    output = [rh2_fit, h_sd]

    if return_fractions:
        output.append(f_H2)
        output.append(f_HI)

    return output

def fit_krumholz(h_sd, rh2, guesses=[10.0, 1.0], rh2_error=None,
        verbose=False, vary=[True, True]):

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

    def chisq(params, h_sd, rh2, rh2_error):
        phi_cnm = params['phi_cnm'].value
        Z = params['Z'].value

        rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z)

        chisq = (rh2 - rh2_model)**2 / rh2_error**2

        return chisq

    # Set parameter limits and initial guesses
    params = Parameters()
    params.add('phi_cnm',
               value=guesses[0],
               min=0.001,
               max=1000,
               vary=vary[0])
    params.add('Z',
               value=guesses[1],
               min=0.01,
               max=100,
               vary=vary[1])

    # Perform the fit!
    result = minimize(chisq,
                      params,
                      args=(h_sd, rh2, rh2_error),
                      method='L-BFGS-B')

    # Print fit results?
    #if verbose:
    #    report_fit(params)

    rh2_fit_params = (params['phi_cnm'].value, params['Z'].value)

    return_chisq = False
    if return_chisq:
        dof = len(rh2) - 1
        rh2_fit_grid = k09.calc_rh2(h_sd, *rh2_fit_params)
        #chisq = np.sum((rh2 - rh2_fit_grid)**2 / rh2_fit_grid)
        chisq = np.sum((rh2 - rh2_fit_grid)**2 / rh2_error**2)
        #if rh2_error is not None:
        #    chisq /= np.sum(rh2_error**2)
        p_value = 1.0 - stats.chi2.cdf(chisq, dof)
        chisq_reduced = chisq / dof

    if verbose:
        #print('pvalue = %.4f' % p_value)
        #print('dof = %.2f' % dof)
        #print('chi2 = %.2f' % chisq)
        pass

    return rh2_fit_params

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
                                                     header=header)[:2]
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

def main(verbose=True):

    '''

    This script requires analysis output from three other scripts.
    Script                                          Purpose
    ---------------------------------------------------------------------------
    hi/taurus_analysis_core_properties.py           Defines core positions
    av/taurus_analysis_av_derive_core_boxes.py      Calculates box sizes
    hi/taurus_analysis_hi_av_core_correlations.py   Calculates HI velocity
                                                    range

    '''

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
    hi_co_width = True
    hi_av_correlation = True
    reproduce_lee12 = True

    # Error analysis
    calc_errors = True # Run monte carlo error analysis?
    vary_phi_cnm = True # Vary phi_cnm in K+09 fit?
    vary_Z = False # Vary metallicity in K+09 fit?

    # Regions
    # Options are 'ds9' or 'av_gradient'
    box_method = 'av_gradient'

    #dgr = 5.33e-2 # dust to gas ratio [10^-22 mag / 10^20 cm^-2
    h_sd_fit_range = [0.001, 1000] # range of fitted values for krumholz model

    # Figures
    write_pdf_figures = True

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/cores/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/hi/'
    co_dir = '/d/bip3/ezbc/taurus/data/co/'
    core_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/taurus/data/python_output/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'

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
    co_vel_axis = make_velocity_axis(co_header)

    # Load global properties of cloud
    # global properties written from script
    # 'av/taurus_analysis_global_properties.txt'
    with open(property_dir + 'taurus_global_properties.txt', 'r') as f:
        properties = json.load(f)
        dgr = properties['dust2gas_ratio']['value']
        dgr_error = properties['dust2gas_ratio_error']['value']
        Z = properties['metallicity']['value']
        hi_vel_range_conf = properties['hi_velocity_range_conf']
        hi_vel_range = properties['hi_velocity_range']
        hi_vel_range_error = properties['hi_velocity_range_error']

    print dgr

    # Plot NHI vs. Av for a given velocity range
    noise_cube_filename = 'taurus_hi_galfa_cube_regrid_planckres_noise.fits'
    if not path.isfile(hi_dir + noise_cube_filename):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=velocity_axis,
                velocity_noise_range=[90,110], header=h, Tsys=30.,
                filename=hi_dir + noise_cube_filename)
    else:
        noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    # define core properties
    with open(core_dir + 'taurus_core_properties.txt', 'r') as f:
        cores = json.load(f)

    cores = convert_core_coordinates(cores, h)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'taurus_av_boxes_',
            header = h)

    # Set up lists
    hi_image_list = []
    hi_sd_image_list = []
    hi_sd_image_error_list = []
    h_sd_image_list = []
    h_sd_image_error_list = []
    av_image_list = []
    av_image_error_list = []
    rh2_image_list = []
    rh2_image_error_list = []
    rh2_fit_list = []
    h_sd_fit_list = []
    phi_cnm_list = []
    phi_cnm_error_list = []
    Z_list = []
    Z_error_list = []
    chisq_list = []
    p_value_list = []
    core_name_list = []
    co_image_list = []
    hi_vel_range_list = []
    hi_vel_range_corr_list = []
    hi_sd_fit_list = []

    hi_data_orig = np.copy(hi_data)
    noise_cube_orig = np.copy(noise_cube)
    av_data_planck_orig = np.copy(av_data_planck)
    av_error_data_planck_orig = np.copy(av_error_data_planck)

    for core in cores:
        print('\nCalculating for core %s' % core)

        if box_method == 'ds9':
            # Grab the mask from the DS9 regions
            xy = cores[core]['box_center_pix']
            box_width = cores[core]['box_width']
            box_height = cores[core]['box_height']
            box_angle = cores[core]['box_angle']
            mask = myg.get_rectangular_mask(av_data_planck_orig,
                    xy[0], xy[1],
                    width = box_width,
                    height = box_height,
                    angle = box_angle)
        elif box_method == 'av_gradient':
            mask = myg.get_polygon_mask(av_data_planck_orig,
                    cores[core]['box_vertices_rotated'])
        else:
            raise ValueError('Method for boxes is either ds9 or av_gradient')

        indices = mask == 1

        # Get only the relevant pixels to decrease computation time
        hi_data = np.copy(hi_data_orig[:, indices])
        noise_cube = np.copy(noise_cube_orig[:, indices])
        av_data_planck = np.copy(av_data_planck_orig[indices])
        av_error_data_planck = np.copy(av_error_data_planck_orig[indices])

        # ---------------------------------------------------------------------
        # Perform analysis on cores, including fitting the Krumholz model.
        # ---------------------------------------------------------------------
        images, params = \
                run_analysis(hi_cube=hi_data,
                             hi_noise_cube=noise_cube,
                             hi_velocity_axis=velocity_axis,
                             hi_header=h,
                             dgr=dgr,
                             dgr_error=dgr_error,
                             av_image=av_data_planck,
                             av_image_error=av_error_data_planck,
                             hi_vel_range=hi_vel_range,
                             hi_vel_range_error=hi_vel_range_error,
                             guesses=[10.0, 1.0],
                             parameter_vary=[vary_phi_cnm, vary_Z],
                             core_dict=cores[core],
                             results_figure_name=figure_dir + \
                                     'monte_carlo_results/' + \
                                     'taurus_%s' % core,
                             )

        rh2_fit, h_sd_fit, f_H2, f_HI = calc_krumholz(params=params[:2],
                                            h_sd_extent=h_sd_fit_range,
                                            return_fractions=True)

        hi_sd_fit = f_HI * h_sd_fit
        phi_cnm = params[0]
        phi_cnm_error = params[2]
        Z = params[1]
        Z_error = params[3]

        # append to the lists
        hi_sd_image_list.append(images['hi_sd'])
        hi_sd_image_error_list.append(images['hi_sd_error'])
        h_sd_image_list.append(images['h_sd'])
        h_sd_image_error_list.append(images['h_sd_error'])
        av_image_list.append(images['av'])
        av_image_error_list.append(images['av_error'])
        rh2_image_list.append(images['rh2'])
        rh2_image_error_list.append(images['rh2_error'])
        rh2_fit_list.append(rh2_fit)
        h_sd_fit_list.append(h_sd_fit)
        hi_sd_fit_list.append(hi_sd_fit)
        phi_cnm_list.append(phi_cnm)
        phi_cnm_error_list.append(phi_cnm_error)
        Z_list.append(Z)
        Z_error_list.append(Z_error)
        #chisq_list.append(chisq)
        #p_value_list.append(p_value)
        core_name_list.append(core)

    # Create the figures!
    # -------------------
    print('\nCreating figures...')

    figure_types = ['png',]
    if write_pdf_figures:
        figure_types.append('pdf')

    for figure_type in figure_types:
        plot_rh2_vs_h_grid(rh2_image_list,
                h_sd_image_list,
                rh2_error_images = rh2_image_error_list,
                h_sd_error_images = h_sd_image_error_list,
                rh2_fits = rh2_fit_list,
                h_sd_fits = h_sd_fit_list,
                limits = [0, 80, 10**-3, 10**2],
                savedir = figure_dir + 'panel_cores/',
                scale = ('linear', 'log'),
                filename = 'taurus_rh2_vs_hsd_panels_planck.%s' % figure_type,
                #title = r'$R_{\rm H2}$ vs. $\Sigma_{\rm HI}$'\
                #        + ' of taurus Cores',
                core_names=core_name_list,
                phi_cnm_list=phi_cnm_list,
                phi_cnm_error_list=phi_cnm_error_list,
                Z_list=Z_list,
                Z_error_list=Z_error_list,
                show = False)

        plot_hi_vs_h_grid(hi_sd_image_list,
                h_sd_image_list,
                hi_sd_error_images = hi_sd_image_error_list,
                h_sd_error_images = h_sd_image_error_list,
                hi_fits = hi_sd_fit_list,
                h_sd_fits = h_sd_fit_list,
                #limits = [1, 100, 1, 100],
                savedir = figure_dir + 'panel_cores/',
                scale = ('linear', 'linear'),
                filename = 'taurus_hi_vs_h_panels_planck_linear.%s' % \
                        figure_type,
                #title = r'$\Sigma_{\rm HI}$ vs. $\Sigma_{\rm H}$'\
                #        + ' of taurus Cores',
                core_names=core_name_list,
                phi_cnm_list=phi_cnm_list,
                show = False)

        plot_hi_vs_h_grid(hi_sd_image_list,
                h_sd_image_list,
                hi_sd_error_images = hi_sd_image_error_list,
                h_sd_error_images = h_sd_image_error_list,
                hi_fits = hi_sd_fit_list,
                h_sd_fits = h_sd_fit_list,
                limits = [1, 100, 1, 100],
                savedir = figure_dir + 'panel_cores/',
                scale = ('log', 'log'),
                filename = 'taurus_hi_vs_h_panels_planck_log.%s' % \
                        figure_type,
                #title = r'$\Sigma_{\rm HI}$ vs. $\Sigma_{\rm H}$'\
                #        + ' of taurus Cores',
                core_names=core_name_list,
                phi_cnm_list=phi_cnm_list,
                show = False)

        plot_hi_vs_av_grid(hi_sd_image_list,
                av_image_list,
                hi_error_images = hi_sd_image_error_list,
                av_error_images = h_sd_image_error_list,
                #limits = [10**-1, 10**2, 10**0, 10**2],
                limits = [1e-1, 70, 2, 10],
                savedir = figure_dir + 'panel_cores/',
                scale = ('log', 'linear'),
                filename = 'taurus_hi_vs_av_panels_planck_log.%s' % \
                        figure_type,
                core_names=core_name_list,
                #title = r'$\Sigma_{\rm HI}$ vs. $\Sigma_{\rm H}$'\
                #        + ' of taurus Cores',
                show = False)


if __name__ == '__main__':
    main()



