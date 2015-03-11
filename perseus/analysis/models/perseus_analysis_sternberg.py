#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the perseus molecular cloud.
'''
import matplotlib
matplotlib.use('Agg')

import pyfits as pf
import numpy as np
import warnings
warnings.filterwarnings('ignore')

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
                marker='^',ecolor='k',linestyle='None',
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
    ax.set_ylabel('N(HI) (1 $\times 10^{20}$ cm$^{-2}$)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return likelihoods_image

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
        elif plot_type is 'scatte':
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
                marker='^',ecolor='k',linestyle='None',
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
    ax.set_ylabel('$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return likelihoods_image

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
                marker='^',ecolor='k',linestyle='None',
                markersize=2
                )

        ax.set_xscale('linear')
        ax.set_yscale(scale)

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel('A$_{\rm V}$ (mag)',)
    ax.set_ylabel('$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return likelihoods_image

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
                marker='^',ecolor='k',linestyle='None',
                markersize=2
                )

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel('A$_{\rm V}$ (mag)',)
        ax.set_ylabel('$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
        ax.set_title(core_names[i])
        ax.grid(True)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return likelihoods_image

def plot_hisd_vs_hsd_grid(hi_sd_images, h_sd_images,
        hi_sd_image_errors=None, h_sd_image_errors=None, limits=None,
        savedir='./', filename=None, show=True, scale=('linear', 'linear'),
        returnimage=False, title=None, core_names='', alphaG_list=None,
        alphaG_error_list=None, Z_list=None,
        Z_error_list=None, phi_g_list=None,
        phi_g_error_list=None,):
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
              'figure.figsize': (7.3, 7.3 * y_scaling),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update()


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
                marker='^',ecolor='k',linestyle='None',
                markersize=2
                )

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{H_2}$ (M$_\odot$ / pc$^2$)',)
        ax.set_ylabel('$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
        ax.set_title(core_names[i])
        ax.grid(True)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return likelihoods_image

def plot_rh2_vs_h_grid(rh2_images, h_sd_images, rh2_error_images=None,
        h_sd_error_images=None, rh2_fits = None, h_sd_fits = None, limits =
        None, fit = True, savedir = './', filename = None, show = True, scale =
        'linear', title = '', core_names='', alphaG_list=None,
        alphaG_error_list=None, Z_list=None,
        Z_error_list=None, phi_g_list=None,
        phi_g_error_list=None,
        chisq_list=None, p_value_list=None,):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from myscience.sternberg14 import calc_T_cnm

    n = int(np.ceil(len(rh2_images)**0.5))
    if n**2 - n > len(rh2_images):
        nrows = n - 1
        ncols = n
        y_scaling = 1.0 - 1.0 / n
    else:
        nrows, ncols = n, n
        y_scaling = 1.0

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
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'axes.weight': line_weight,
              'axes.labelweight': font_weight,
              'axes.color_cycle': color_cycle, # colors of different plots
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'xtick.weight': line_weight,
              'ytick.labelsize': font_scale,
              'ytick.weight': line_weight,
              'font.weight': font_weight,
              #'font.family': 'sans-serif',
              'text.fontsize': font_scale,
              'text.usetex': True,
              'figure.figsize': (7.3, 7.3 * y_scaling),
              'figure.titlesize': font_scale,
             }
    plt.rcParams.update()

    # Create figure instance
    fig = plt.figure()

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=len(rh2_images),
                 axes_pad=0.0,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    # Cycle through lists

    if len(alphaG_error_list) == 0:
        alphaG_error_list = None
    if len(Z_error_list) == 0:
        Z_error_list = None

    for i in xrange(len(rh2_images)):
        rh2 = rh2_images[i]
        h_sd = h_sd_images[i]
        rh2_error = rh2_error_images[i]
        h_sd_error = h_sd_error_images[i]
        rh2_fit = rh2_fits[i]
        h_sd_fit = h_sd_fits[i]
        if alphaG_list is not None:
            alphaG = alphaG_list[i]
        if alphaG_error_list is not None:
            alphaG_error = alphaG_error_list[i]
        if Z_list is not None:
            Z = Z_list[i]
        if Z_error_list is not None:
            Z_error = Z_error_list[i]
        if phi_g_list is not None:
            phi_g = phi_g_list[i]
        if phi_g_error_list is not None:
            phi_g_error = phi_g_error_list[i]

        # Drop the NaNs from the images
        if type(rh2_error) is float:
            indices = np.where((rh2 == rh2) &\
                               (h_sd == h_sd)&\
                               (h_sd > 0) &\
                               (rh2 > 0))

        if type(rh2_error) is np.ndarray or \
                type(rh2_error) is np.ma.core.MaskedArray or \
                type(h_sd_error) is np.ndarray or \
                type(h_sd_error) is np.ma.core.MaskedArray or \
                type(rh2_error) is list:

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
                alpha=0.3,
                color='k',
                marker='^',
                ecolor='k',
                linestyle='None',
                markersize=4
                )

        if rh2_fit is not None:
            ax.plot(h_sd_fit, rh2_fit,
                    color = 'r', alpha=1)

        # Annotations
        anno_xpos = 0.95

        if alphaG_list is not None and Z_list is not None:
            if alphaG_error_list is None and Z_error_list is not None:
                ax.annotate(r'$\alpha G$ = {0:.2f}\n'.format(alphaG) + \
                            r'Z = {0:.2f} Z$_{\odot}$'.format(Z),
                        xytext=(anno_xpos, 0.05),
                        xy=(anno_xpos, 0.05),
                        textcoords='axes fraction',
                        xycoords='axes fraction',
                        color='k',
                        bbox=dict(boxstyle='round',
                                  facecolor='w',
                                  alpha=0.3),
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        )

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'$\Sigma_{HI}$ + $\Sigma_{H2}$ (M$_\odot$ / pc$^2$)',)
        ax.set_ylabel(r'R$_{H2}$ = $\Sigma_{H2}$ / $\Sigma_{HI}$',)
        ax.set_title(core_names[i])
        ax.grid(False)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        plt.show()

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
              'figure.figsize': (7.3, 7.3 * y_scaling),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update()

    # Create figure instance
    fig = plt.figure()

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=len(hi_images),
                 axes_pad=0,
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

        if 0:
            image = ax.errorbar(av_nonans.ravel(),
                    hi_nonans.ravel(),
                    xerr=(av_error_nonans.ravel()),
                    yerr=(hi_error_nonans.ravel()),
                    alpha=0.5,
                    color='k',
                    marker='^',ecolor='k',linestyle='None',
                    markersize=4
                    )
        else:
            h, x, y = np.histogram2d(av_nonans.ravel(),
                                hi_nonans.ravel(),
                                bins=100,)
            ax.imshow(h,
                      origin='lower',
                      interpolation='bicubic',
                      extent=[x[0], x[-1], y[0], y[-1]])

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel('$A_V$ (Mag)',)
        ax.set_ylabel('$\Sigma_{HI}$ ($M_\odot$ pc$^{-2}$)',)
        #ax.set_title(core_names[i])
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
        'linear', title = '', core_names='', alphaG_list=None,
        alphaG_error_list=None, Z_list=None, Z_error_list=None,
        phi_g_list=None, phi_g_error_list=None,):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from myscience.sternberg14 import calc_T_cnm

    n = int(np.ceil(len(hi_images)**0.5))
    if n**2 - n > len(hi_images):
        nrows = n - 1
        ncols = 2
        y_scaling = 1.0 - 1.0 / n
    else:
        nrows, ncols = n, n
        y_scaling = 1.0

    n = len(hi_images)
    ncols = 2
    nrows = (n + 1) / ncols
    y_scaling = nrows / float(ncols)

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
              'figure.figsize': (3.6, 3.6*y_scaling),
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)


    # Create figure instance
    fig = plt.figure()

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=n,
                 axes_pad=0.,
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
        if alphaG_list is not None:
            alphaG = alphaG_list[i]
        if alphaG_error_list is not None:
            alphaG_error = alphaG_error_list[i]
        if Z_list is not None:
            Z = Z_list[i]
        if Z_error_list is not None:
            Z_error = Z_error_list[i]
        if phi_g_list is not None:
            phi_g = phi_g_list[i]
        if phi_g_error_list is not None:
            phi_g_error = phi_g_error_list[i]

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

        if 1:
            image = ax.errorbar(h_sd_nonans.ravel(),
                    hi_sd_nonans.ravel(),
                    xerr=(h_sd_error_nonans.ravel()),
                    yerr=(hi_sd_error_nonans.ravel()),
                    alpha=0.1,
                    color='k',
                    marker='^',ecolor='k',linestyle='None',
                    markersize=3
                    )
        elif 0:
            ax.hexbin(h_sd_nonans.ravel(),
                      hi_sd_nonans.ravel(),
                      cmap=cmap,
                      mincnt=1)
        else:
            h, x, y = np.histogram2d(h_sd_nonans.ravel(),
                                hi_sd_nonans.ravel(),
                                bins=1000,)
            ax.imshow(h,
                      origin='lower',
                      interpolation='bicubic',
                      extent=[x[0], x[-1], y[0], y[-1]])


        if hi_sd_fit is not None:
            ax.plot(h_sd_fit, hi_sd_fit,
                    color = 'r')
        # Annotations
        anno_xpos = 0.95

        if alphaG_list is not None and Z_list is not None:
            T_cnm = calc_T_cnm(alphaG, Z=Z)
            T_cnm_error = []
            T_cnm_error.append(\
                    T_cnm - calc_T_cnm(alphaG + alphaG_error[0], Z=Z))
            T_cnm_error.append(\
                    T_cnm - calc_T_cnm(alphaG + alphaG_error[1], Z=Z))

            alphaG_text = r'\noindent$\alpha G$ =' + \
                           r' %.2f' % (alphaG) + \
                           r'$^{+%.2f}_{-%.2f}$ \\' % (alphaG_error[0],
                                                     alphaG_error[1])
            T_cnm_text = r'\noindent T$_{\rm CNM}$ =' + \
                         r' %.2f' % (T_cnm) + \
                         r'$^{+%.2f}_{-%.2f}$ \\' % (T_cnm_error[0],
                                                     T_cnm_error[1])
            if Z != 1:
                if Z_error == (0.0, 0.0):
                    Z_text = r'Z = %.1f Z$_\odot$ \\' % (Z)
                else:
                    Z_text = r'Z = %.2f' % (Z) + \
                    r'$^{+%.2f}_{-%.2f}$ Z$_\odot$ \\' % (Z_error[0],
                                                          Z_error[1])
            else:
                Z_text = ''

            if phi_g_error == (0.0, 0.0):
                phi_g_text = r'\noindent$\phi_{\rm mol}$ = ' + \
                                 '%.1f' % (phi_g)
            else:
                phi_g_text = r'\noindent$\phi_{\rm mol}$ =' + \
                            r' %.2f' % (phi_g) + \
                            r'$^{+%.2f}_{-%.2f}$ \\' % (phi_g_error[0],
                                                     phi_g_error[1])

            ax.annotate(alphaG_text + T_cnm_text + Z_text,
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

            ax.annotate(core_names[i],
                        xytext=(0.95, 0.95),
                        xy=(0.95, 0.95),
                        textcoords='axes fraction',
                        xycoords='axes fraction',
                        size=font_scale,
                        color='k',
                        bbox=dict(boxstyle='square',
                                  facecolor='w',
                                  alpha=1),
                        horizontalalignment='right',
                        verticalalignment='top',
                        )

        ax.set_xscale(scale[0])
        ax.set_yscale(scale[1])

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel(r'$\Sigma_{HI}$ + $\Sigma_{H2}$ ' + \
                       '(M$_\odot$ / pc$^2$)',)
        ax.set_ylabel(r'$\Sigma_{HI}$',)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
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
            hi_velocity_range_likelihood = vel_range_hiav_list[i]
        except IndexError:
            hi_velocity_range_likelihood = None

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
                    color = 'r', alpha=0.2, label='$^{12}$CO Spectrum Width')
        if hi_velocity_range_likelihood is not None:
            ax.axvspan(hi_velocity_range_likelihood[0], hi_velocity_range_likelihood[1],
                    color = 'b', alpha=0.2, label='HI/A$_V$ Correlation')

        ax.set_xscale(scale[0], nonposx = 'clip')
        ax.set_yscale(scale[1], nonposy = 'clip')

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        # Adjust asthetics
        ax.set_xlabel('Velocity (km/s)')
        ax.set_ylabel('<I$_{\rm CO}$> (K)',)
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

def plot_parameter_hist(results_dict, parameter='Z fits',
        results_figure_name=None, show=False):

    import matplotlib.pyplot as plt
    import numpy as np

    if parameter == 'Z fits':
        bins = np.linspace(0, 10, 100)
    elif parameter == 'alphaG fits':
        bins = 100

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    ax.hist(results_dict[parameter],
            bins=bins)
    ax.set_ylabel('Counts')
    if parameter == 'Z fits':
        ax.set_xlabel(r'$Z$ $(Z_\odot)$')
        plt.savefig(results_figure_name + '_Z_hist.png')
    elif parameter == 'alphaG fits':
        ax.set_xlabel(r'$\alpha G$')
        plt.savefig(results_figure_name + '_alphaG_hist.png')

    if show:
        plt.show()

    plt.close()

def recreate_PDFs(vel_centers=None, vel_widths=None, dgrs=None,
        center_likelihoods=None, width_likelihoods=None, dgr_likelihoods=None,
        center_rv=None, width_rv=None, dgr_rv=None, results_figure_name=None):

    import matplotlib.pyplot as plt

    # Recreate the distribution of likelihoods
    center_likelihoods_recreate = np.zeros(10000)
    width_likelihoods_recreate = np.zeros(10000)
    dgr_likelihoods_recreate = np.zeros(10000)
    for i in range(len(center_likelihoods_recreate)):
        center_likelihoods_recreate[i] = center_rv.rvs()
        width_likelihoods_recreate[i] = width_rv.rvs()
        dgr_likelihoods_recreate[i] = dgr_rv.rvs() * 1000.0

    center_likelihood_normed = center_likelihoods / np.sum(center_likelihoods)
    width_likelihood_normed = width_likelihoods / np.sum(width_likelihoods)
    dgr_likelihood_normed = dgr_likelihoods / np.sum(dgr_likelihoods)

    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
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

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)
    center_bins = np.arange(vel_centers[0], vel_centers[-1] + 2, 1)
    width_bins = np.arange(vel_widths[0], vel_widths[-1] + 2, 1)
    ax.hist(center_likelihoods_recreate, bins=center_bins, alpha=0.5,
            label='Centers Reproduced', color='b', normed=True)
    ax.hist(width_likelihoods_recreate, bins=width_bins, alpha=0.5,
            label='Widths Reproduced', color='r', normed=True)
    ax.plot(vel_centers, center_likelihood_normed, alpha=0.5,
            label='Centers', color='k')
    ax.plot(vel_widths, width_likelihood_normed, alpha=0.5,
            label='Widths', color='g')
    ax.legend(fontsize=font_scale * 3/4.0)
    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel('Normalized value')
    plt.savefig(results_figure_name + '_PDF_hist.png')

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)
    dgr_bins = np.arange(dgrs[0], dgrs[-1] + 2, 0.02)
    ax.hist(dgr_likelihoods_recreate, bins=dgr_bins, alpha=0.5,
            label='DGRs Reproduced', color='k', normed=True)
    ax.plot(dgrs, dgr_likelihood_normed, alpha=0.5,
            label='Widths', color='k')
    ax.legend(fontsize=font_scale * 3/4.0)
    ax.set_xlabel('DGR')
    ax.set_ylabel('Normalized value')
    plt.savefig(results_figure_name + '_PDF_dgr_hist.png')

''' Calculations
'''

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
        hi_vel_range=None, hi_header=None, dgr=None, intercept=None,
        av_image=None, av_image_error=None, sub_image_indices=None,
        nhi_error=None):

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
    nh2_image = calculate_nh2(nhi_image=nhi_image,
                              av_image=av_image + intercept,
                              dgr=dgr)

    nh2_image_error = calculate_nh2_error(nhi_image_error=nhi_image_error,
                                          av_image_error=av_image_error,
                                          dgr=dgr)

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi_image,
                               sd_factor=1/1.25)
    hi_sd_image_error = calculate_sd(nhi_image_error,
                                     sd_factor=1/1.25)


    h2_sd_image = calculate_sd(nh2_image,
                               sd_factor=1/0.625)
    h2_sd_image_error = calculate_sd(nh2_image_error,
                                     sd_factor=1/0.625)


    h_sd_image = hi_sd_image + h2_sd_image
    h_sd_image_error = (hi_sd_image_error**2 + h2_sd_image_error**2)**0.5

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
        hi_header=None, dgr=None, dgr_error=None, intercept=None,
        av_image=None, av_image_error=None, vel_center=None,
        hi_vel_range=None, hi_vel_range_error=None, verbose=False,
        core_dict=None, results_figure_name=None, properties=None,
        results_filename=None):

    '''

    Runs Monte Carlo to vary measure the error induced in RH2 from the HI
    integration velocity range. If N_monte_carlo_runs = 1 then a simple fit
    without adding any random errors is performed.

    Parameters
    ----------
    guesses : tuple, list
        Initial guesses for alphaG and Z.
    paramter_vary : tuple, bool
        Vary alphaG and Z?
    hi_vel_range : tuple
        Tuple of floats, lower and upper bound to HI velocity range.
    hi_vel_range_error : float
        1 sigma Gaussian error on hi_vel_range limits. If error_method ==
        'gaussian', and hi_vel_range_error is None then the confidence limit
        from the HI velocity range error in the core dict will be used.
    N_monte_carlo_runs : int
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
        Parameters in sternberg fit. If N_monte_carlo_runs > 1, then errors on the
        parameters are returned as well.


    '''

    from numpy.random import normal
    from scipy.stats import rv_discrete
    import mystats
    from mystats import rv2d_discrete, rv3d_discrete
    import matplotlib.pyplot as plt
    import json
    from os import path

    if N_monte_carlo_runs < 1:
        raise ValueError('N_monte_carlo_runs must be >= 1')

    verbose = False

    if results_filename is not None:
        if not path.isfile(results_filename):
            perform_mc = True
            write_mc = True
        elif clobber:
            perform_mc = True
            write_mc = True
        else:
            perform_mc = False
            write_mc = False
    # If no filename provided, do not read file and do not write file
    else:
        write_mc = False
        perform_mc = True

    # Results of monte carlo will be stored here
    results_dict = {'alphaG fits' : np.empty((N_monte_carlo_runs)),
                    'Z fits' : np.empty((N_monte_carlo_runs)),
                    'phi_g fits' : np.empty((N_monte_carlo_runs)),
                    'nhi errors' : np.empty((N_monte_carlo_runs))}
    hi_vel_range_list = np.empty((N_monte_carlo_runs, 2))
    dgr_list = np.empty((N_monte_carlo_runs))
    intercept_list = np.empty((N_monte_carlo_runs))

    # Get standard errors on images + the HI velocity range
    hi_error = np.median(hi_noise_cube)
    av_error = np.median(av_image_error)

    if likelihood_derivation == 'cores':
        likelihood_param_dict = core_dict
    elif likelihood_derivation == 'global':
        likelihood_param_dict = properties
    else:
    	raise ValueError('likelihood_derivation must be either cores or ' + \
    	        'global')

    # Load the calculate marginalized PDFs for each parameter
    width_likelihoods = np.asarray(likelihood_param_dict['width_likelihood'])
    dgr_likelihoods = np.asarray(likelihood_param_dict['dgr_likelihood'])
    intercept_likelihoods = \
        np.asarray(likelihood_param_dict['intercept_likelihood'])

    # Load the full 3D likelihood space
    likelihoods = np.asarray(likelihood_param_dict['likelihoods'])

    # Load the grids
    vel_widths = np.asarray(likelihood_param_dict['vel_widths'])
    dgrs = np.asarray(likelihood_param_dict['dgrs'])
    intercepts = np.asarray(likelihood_param_dict['intercepts'])

    # Derive PDF of the velocity centers and widths.
    # The Monte Carlo will draw from this PDF randomly.
    # rv_discrete requires the PDF be normalized to have area of 1
    width_likelihood_normed = width_likelihoods / np.sum(width_likelihoods)
    dgr_likelihood_normed = dgr_likelihoods / np.sum(dgr_likelihoods)
    intercept_likelihood_normed = intercept_likelihoods / \
                                  np.sum(intercept_likelihoods)
    likelihoods_normed = likelihoods / np.sum(likelihoods)

    # Create 2-D PDF of velocity widths + DGRs
    likelihood_pdf = rv3d_discrete(likelihoods=likelihoods_normed,
                                   param_grid1=vel_widths,
                                   param_grid2=dgrs,
                                   param_grid3=intercepts,
                                   param_name1='widths',
                                   param_name2='dgrs',
                                   param_name3='intercepts',
                                   #L_scalar=N_monte_carlo_runs,
                                   )

    if perform_mc:
        # Run the Monte Carlo
        # -------------------

        # Progress bar parameters
        total = float(N_monte_carlo_runs)
        count = 0

        for i in xrange(N_monte_carlo_runs):
            if N_monte_carlo_runs > 1:
                # Randomly sample from Gaussian distribution
                hi_random_error = normal(scale=hi_error,\
                        size=hi_noise_cube.shape)
                av_random_error = normal(scale=av_error, size=av_image.shape)

                # Get a random sample from the DGR / width /
                # intercept likelihood space
                vel_width_random, dgr_random, intercept_random = \
                        likelihood_pdf.rvs()

                # Create the velocity range
                hi_vel_range_noise = (vel_center - \
                                        vel_width_random / 2.,
                                      vel_center + \
                                        vel_width_random / 2.)

                hi_vel_range_list[i, 0] = hi_vel_range_noise[0]
                hi_vel_range_list[i, 1] = hi_vel_range_noise[1]
                dgr_list[i] = dgr_random
                intercept_list[i] = intercept_random

                # Add random error to images
                hi_cube_noise = np.copy(hi_cube) + hi_random_error
                av_image_noise = np.copy(av_image) + av_random_error
            elif N_monte_carlo_runs == 1:
                hi_cube_noise = np.copy(hi_cube)
                av_image_noise = np.copy(av_image)
                dgr_noise = dgr
                hi_vel_range_noise = np.asarray(hi_vel_range)

            # Derive new images
            images = derive_images(hi_cube=hi_cube_noise,
                                   hi_velocity_axis=hi_velocity_axis,
                                   hi_noise_cube=hi_noise_cube,
                                   hi_vel_range=hi_vel_range_noise,
                                   hi_header=hi_header,
                                   dgr=dgr_random,
                                   intercept=intercept_random,
                                   av_image=av_image_noise,
                                   av_image_error=av_image_error,
                                   )
            # grab the N(HI) error
            results_dict['nhi errors'][i] = \
                    np.mean(images['nhi'][images['nhi'] == images['nhi']])

            # Fit R_H2
            #---------
            # Unravel images to single dimension
            rh2_ravel = images['rh2'].ravel()
            rh2_error_ravel = images['rh2_error'].ravel()
            h_sd_ravel = images['h_sd'].ravel()
            h_sd_error_ravel = images['h_sd_error'].ravel()

            # write indices for only ratios > 0
            indices = np.where((rh2_ravel > 1) & \
                               (rh2_ravel == rh2_ravel) & \
                               (rh2_error_ravel == rh2_error_ravel))
            rh2_ravel = rh2_ravel[indices]
            rh2_error_ravel = rh2_error_ravel[indices]
            h_sd_ravel = h_sd_ravel[indices]
            h_sd_error_ravel = h_sd_error_ravel[indices]

            # Fit to sternberg model, init guess of alphaG = 10
            alphaG, Z, phi_g = fit_sternberg(h_sd_ravel,
                                      rh2_ravel,
                                      guesses=guesses, # alphaG, Z
                                      vary=[vary_alphaG,
                                            vary_Z, vary_phi_g],
                                      rh2_error=rh2_error_ravel,
                                      verbose=verbose)

            # keep results
            results_dict['alphaG fits'][i] = alphaG
            results_dict['phi_g fits'][i] = phi_g
            results_dict['Z fits'][i] = Z

            if verbose:
                print('phi = %.2f' % alphaG)
                print('phi = %.2f' % phi_g)
                print('Z = %.2f' % Z)

            # see eq 6 of sternberg+09
            # alphaG is the number density of the CNM over the minimum number
            # density required for pressure balance
            # the lower alphaG values than for perseus mean that perseus
            # has a more diffuse CNM

            # By fitting the model to the observed R_H2 vs total H, you
            # basically constrained psi in Equation (35) of sternberg+09.  This
            # means that you can calculate f_H2 for a given total hydrogen
            # surface density.  In this case, H2 surface density = f_H2 *
            # total hydrogen surface density HI surface density = (1 - f_HI) *
            # total hydrogen surface density

            # Shows progress each 10%
            count += 1
            abs_step = int((total * 1)/10) or 10
            if count and not count % abs_step:
                print "\t{0:.0%} processed".format(count/total)

        # Write the MC results in JSON human readable txt format
        if write_mc:
            alphaGs = results_dict['alphaG fits'].tolist()
            phi_gs = results_dict['phi_g fits'].tolist()
            Zs = results_dict['Z fits'].tolist()
            vel_ranges = hi_vel_range_list.tolist()
            dgr_list = dgr_list.tolist()
            intercept_list = intercept_list.tolist()

            results = {'alphaGs': alphaGs,
                       'phi_gs': phi_gs,
                       'Zs': Zs,
                       'dgrs': dgr_list,
                       'intercepts': intercept_list,
                       'vel_ranges': hi_vel_range_list.tolist()}

            with open(results_filename, 'w') as f:
                json.dump(results, f)

    # If not performing the MC, read a previously written MC file
    elif not perform_mc:
        print('Reading monte carlo results file:')
        print(results_filename)

        with open(results_filename, 'r') as f:
            results = json.load(f)

        results_dict['alphaG fits'] = np.asarray(results['alphaGs'])
        results_dict['phi_g fits'] = np.asarray(results['phi_gs'])
        results_dict['Z fits'] = np.asarray(results['Zs'])
        dgr_list = np.asarray(results['dgrs'])
        intercept_list = np.asarray(results['intercepts'])
        hi_vel_range_list = np.asarray(results['vel_ranges'])

    # Remove failed fits
    results_dict['alphaG fits'] = \
        results_dict['alphaG fits']\
            [~np.isnan(results_dict['alphaG fits'])]
    results_dict['phi_g fits'] = \
        results_dict['phi_g fits']\
            [~np.isnan(results_dict['phi_g fits'])]
    results_dict['Z fits'] = \
        results_dict['Z fits'] \
            [~np.isnan(results_dict['Z fits'])]

    # Plot the distributions of parameters from the MC
    if results_figure_name is not None and perform_mc:
    	plot_parameter_hist(results_dict,
    	                    parameter='alphaG fits',
    	                    results_figure_name=results_figure_name,
                            show=False)
    	plot_parameter_hist(results_dict,
    	                    parameter='Z fits',
    	                    results_figure_name=results_figure_name)

    # Derive images
    # -------------
    hi_vel_range_sample = (np.median(hi_vel_range_list[:, 0]),
                           np.median(hi_vel_range_list[:, 1]))

    nhi_error = np.std(results_dict['nhi errors'])

    images = derive_images(hi_cube=hi_cube,
                           hi_velocity_axis=hi_velocity_axis,
                           hi_noise_cube=hi_noise_cube,
                           hi_vel_range=hi_vel_range_sample,
                           hi_header=hi_header,
                           dgr=dgr,
                           intercept=intercept,
                           av_image=av_image,
                           av_image_error=av_image_error,
                           )

    alphaG_confint, Z_confint, phi_g_confint = calc_MC_errors(results_dict,
                                                error_method=error_method,
                                                alpha=alpha,
                                                parameter_vary=[vary_alphaG,
                                                    vary_Z, vary_phi_g])

    alphaG = alphaG_confint[0]
    alphaG_error = alphaG_confint[1:]
    phi_g = phi_g_confint[0]
    phi_g_error = phi_g_confint[1:]
    Z = Z_confint[0]
    Z_error = Z_confint[1:]

    # Print results
    print('results are:')
    print('alphaG = {0:.2f} +{1:.2f}/-{2:.2f}'.format(alphaG_confint[0],
                                                      alphaG_confint[1],
                                                      alphaG_confint[2]))
    print('phi_g = {0:.2f} +{1:.2f}/-{2:.2f}'.format(phi_g_confint[0],
                                                       phi_g_confint[1],
                                                       phi_g_confint[2]))
    print('Z = {0:.2f} +{1:.2f}/-{2:.2f}'.format(Z_confint[0],
                                                  Z_confint[1],
                                                  Z_confint[2]))

    print('Median HI velocity range:' + \
          '{0:.1f} to {1:.1f} km/s'.format(hi_vel_range_sample[0],
                                           hi_vel_range_sample[1],))

    if N_monte_carlo_runs > 1:
        return images, hi_vel_range_sample, (alphaG, Z, alphaG_error,
                Z_error, phi_g, phi_g_error)
    if N_monte_carlo_runs == 1:
        return images, hi_vel_range_sample, (alphaG, Z, phi_g)

def calc_MC_errors(results_dict, error_method='edges', alpha=0.05,
        parameter_vary=[True, False, True]):

    from mystats import calc_symmetric_error
    from scikits.bootstrap import ci
    import numpy as np

    if error_method == 'bootstrap':
        # Bootstrap for errors
        # --------------------
        # Returns errors of a bootstrap simulation at the 100.*(1 - alpha)
        # confidence interval. Errors are computed by deriving a cumulative
        # distribution function of the medians of the sampled data and
        # determining the distance between the median and the value including
        # alpha/2 % of the data, and the value including alpha / 2 % of the
        # data.
        # samples = mystats.bootstrap(results_dict['alphaG fits'], 1000)
        # alphaG_confint = mystats.calc_bootstrap_error(samples, alpha)
        # samples = mystats.bootstrap(results_dict['Z fits'], 1000)
        # Z_confint = mystats.calc_bootstrap_error(samples, alpha)
        if parameter_vary[0]:
            alphaG_confint = ci(results_dict['alphaG fits'],
                                 statfunction=np.median,
                                 alpha=alpha,
                                 method='bca')
            alphaG = np.median(results_dict['alphaG fits'])
            alphaG_confint = (alphaG,
                               alphaG + alphaG_confint[1],
                               alphaG - alphaG_confint[0],
                               )
        else:
            alphaG_confint = (results_dict['alphaG fits'][0], 0.0, 0.0)

        if parameter_vary[1]:
            Z_confint = ci(results_dict['Z fits'],
                           statfunction=np.mean,
                           alpha=alpha)
            Z = np.median(results_dict['Z fits'])
            Z_confint = (Z,
                         Z + Z_confint[0],
                         Z - Z_confint[1],
                         )
        else:
            Z_confint = (results_dict['Z fits'][0], 0.0, 0.0)
    elif error_method == 'edges':
        # If there is a distribution of the parameter, then find the
        # confidence interval
        if parameter_vary[0]:
            alphaG = results_dict['alphaG fits']

            # Histogram will act as distribution of parameter values
            counts, bins = np.histogram(results_dict['alphaG fits'],
                bins=np.linspace(alphaG.min(), alphaG.max(), 500))

            # Calculate confidence interval by taking the PDF weighted mean as
            # value, where errors are calculated by moving vertical threshold
            # from edges towards mean until alpha/2 *100 % of the area of the
            # PDF is included within the threshold
            #alphaG_confint = calc_symmetric_error(bins[:-1], counts,
            #        alpha=alpha)
            alphaG_confint = calc_symmetric_error(alphaG,
                                                   alpha=alpha)

        else:
            alphaG_confint = (results_dict['alphaG fits'][0], 0.0, 0.0)

        if parameter_vary[1]:
            Z_upper_lim = np.log10(np.max(results_dict['Z fits']))
            Z_lower_lim = np.log10(np.min(results_dict['Z fits']))
            Z = results_dict['Z fits']
            counts, bins = np.histogram(results_dict['Z fits'],
                bins=np.linspace(Z.min(), Z.max(), 1000))
            Z_confint = calc_symmetric_error(Z,
                    alpha=alpha)
        else:
            Z_confint = (results_dict['Z fits'][0], 0.0, 0.0)

        if parameter_vary[2]:
            phi_g_upper_lim = np.log10(np.max(results_dict['phi_g fits']))
            phi_g_lower_lim = np.log10(np.min(results_dict['phi_g fits']))
            phi_g = results_dict['phi_g fits']
            phi_g_confint = calc_symmetric_error(phi_g,
                    alpha=alpha)
        else:
            phi_g_confint = (results_dict['phi_g fits'][0], 0.0, 0.0)

    else:
        raise ValueError('Error method must be "bootstrap" or "threshold"')

    return alphaG_confint, Z_confint, phi_g_confint

''' Fitting Functions
'''


def calc_sternberg(params=[10.0, 1.0], h_sd_extent=(0.001, 500),
        return_fractions=True):

    '''
    Parameters
    ----------
    alphaG, Z : float
        alphaG and Z parameters for sternberg model.
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
    from myscience import sternberg14 as s14

    # Create large array of h_sd
    h_sd_extent.append(1e4)
    h_sd = np.linspace(h_sd_extent[0], h_sd_extent[1], h_sd_extent[2])

    if not return_fractions:
        rh2_fit = s14.calc_rh2(h_sd, params)
    elif return_fractions:
        rh2_fit, f_H2, f_HI = s14.calc_rh2(h_sd,
                                           *params,
                                           return_fractions=True)

    output = [rh2_fit, h_sd]

    if return_fractions:
        output.append(f_H2)
        output.append(f_HI)

    return output

def fit_sternberg(h_sd, rh2, guesses=[10.0, 1.0, 10.0], rh2_error=None,
        verbose=False, vary=[True, True, True]):

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
    from myscience import sternberg14 as s14

    def chisq(params, h_sd, rh2, rh2_error):
        alphaG = params['alphaG'].value
        phi_g = params['phi_g'].value
        Z = params['Z'].value

        rh2_model = s14.calc_rh2(h_sd, alphaG, Z, phi_g=phi_g,
                                 return_fractions=False)

        chisq = np.sum((rh2 - rh2_model)**2 / rh2_error**2)

        print('chisq = {0:.1e}, alphaG = {1:.1f}'.format(chisq, alphaG))
        if 0:
            import matplotlib.pyplot as plt
            plt.plot(rh2 - rh2_model, linestyle='', marker='o')
            plt.show()
        if 0:
            import matplotlib.pyplot as plt
            plt.plot(h_sd, rh2, linestyle='', marker='o')
            plt.plot(h_sd, rh2_model, linestyle='-', marker='')
            plt.show()
        print np.asarray(rh2_model).shape

        return chisq

    # Set parameter limits and initial guesses
    params = Parameters()
    params.add('alphaG',
               value=guesses[0],
               min=0.1,
               max=10,
               vary=vary[0])
    params.add('phi_g',
               value=guesses[2],
               min=0.5,
               max=2,
               vary=vary[2])
    params.add('Z',
               value=guesses[1],
               min=0.1,
               max=4,
               vary=vary[1])

    # Perform the fit!
    result = minimize(chisq,
                      params,
                      args=(h_sd, rh2, rh2_error),
                      method='nelder-mead')

    rh2_fit_params = (params['alphaG'].value, params['Z'].value,
            params['phi_g'].value)

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

    return region

def load_ds9_core_region(cores, filename_base = 'perseus_av_boxes_',
        header=None):

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    regions = read_ds9_region(filename_base + '.reg')


    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        core = tag[tag.find('{')+1:tag.find('}')]

        if core in cores:

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

            cores[core]['poly_verts'] = {}
            cores[core]['poly_verts']['wcs'] = poly_verts
            cores[core]['poly_verts']['pixel'] = poly_verts_pix

    return cores

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
        region_name = tag[tag.find('{')+1:tag.find('}')].lower()

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

'''
The main script
'''

def main(verbose=True, av_data_type='planck', region=None):

    '''

    This script requires analysis output from three other scripts.
    Script                                          Purpose
    ---------------------------------------------------------------------------
    hi/perseus_analysis_core_properties.py           Defines core positions
    av/perseus_analysis_av_derive_core_boxes.py      Calculates box sizes
    hi/perseus_analysis_hi_av_core_likelihoods.py    Calculates HI velocity
                                                        range
    av/perseus_analysis_av_load_regions.py

    '''

    import numpy as np
    import numpy
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    from myscience.sternberg14 import calc_T_cnm

    # parameters used in script
    # -------------------------
    # HI velocity integration range
    # Determine HI integration velocity by CO or correlation with Av?
    hi_co_width = True
    hi_av_likelihoodelation = True
    reproduce_lee12 = True

    # Error analysis
    global clobber
    global vary_alphaG
    global vary_Z
    global vary_phi_g
    global error_method
    global alpha
    global N_monte_carlo_runs
    global calc_errors
    global guesses

    calc_errors = True # Run monte carlo error analysis?
    N_monte_carlo_runs = 50 # Number of monte carlo runs
    vary_alphaG = True # Vary alphaG in K+09 fit?
    vary_Z = False # Vary metallicity in K+09 fit?
    vary_phi_g = False # Vary phi_g in K+09 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence

    # Monte carlo results file bases
    results_filename = '/d/bip3/ezbc/perseus/data/python_output/' + \
            'monte_carlo_results/perseus_mc_results_' + av_data_type + '_'

    # global property filename
    global_property_filename = 'perseus_global_properties'

    # Name correct region of cloud
    if region == 1:
        region_name = 'perseus1'
    elif region == 2:
        region_name = 'perseus2'
    else:
        region_name = 'perseus'

    global_property_filename = \
            global_property_filename.replace('perseus', region_name)

    # Which cores to include in analysis?
    cores_to_keep = ('L1495', 'L1495A', 'B213', 'L1498', 'B215',
                     'L1483', 'L1478', 'L1456', 'NGC1579',
                     'B5', 'IC348', 'B1E', 'B1', 'NGC1333', 'L1482')

    # perform MC and write over current results?
    clobber = 1

    # Guesses for (alphaG, Z, phi_g)
    guesses=(1.0, 1.0, 1.0)

    # Use core-derived or global-derived likelihoods forr DGR - vel width
    # combinations. Options are 'cores' and 'global'
    global likelihood_derivation
    likelihood_derivation = 'global'

    # Regions
    # Options are 'ds9' or 'av_gradient'
    global box_method
    global region_type
    box_method = 'av_gradient'
    box_method = 'ds9'
    region_type = 'wedge' # Shape of core region, box or wedge

    #dgr = 5.33e-2 # dust to gas ratio [10^-22 mag / 10^20 cm^-2
    global h_sd_fit_range
    h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

    # Figures
    write_pdf_figures = True

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/perseus/figures/cores/'
    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    hi_dir = '/d/bip3/ezbc/perseus/data/hi/'
    co_dir = '/d/bip3/ezbc/perseus/data/co/'
    core_dir = '/d/bip3/ezbc/perseus/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/perseus/data/python_output/'
    region_dir = '/d/bip3/ezbc/perseus/data/python_output/ds9_regions/'
    multicloud_region_dir = \
            '/d/bip3/ezbc/multicloud/data/python_output/'

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

    co_data, co_header = load_fits(co_dir + \
                'perseus_co_cfa_cube_regrid_planckres.fits',
            return_header=True)

    # make the velocity axis
    velocity_axis = make_velocity_axis(h)
    co_vel_axis = make_velocity_axis(co_header)

    # Load global properties of cloud
    # global properties written from script
    # 'av/perseus_analysis_global_properties.txt'
    with open(property_dir + \
              global_property_filename + '_' + av_data_type + \
              '_scaled.txt', 'r') as f:
        properties = json.load(f)
        dgr = properties['dust2gas_ratio']['value']
        intercept = properties['intercept']['value']
        intercept_error = properties['intercept_error']['value']
        dgr_error = properties['dust2gas_ratio_error']['value']
        vel_center = properties['hi_velocity_center']['value']
        Z = properties['metallicity']['value']

    # Plot NHI vs. Av for a given velocity range
    noise_cube_filename = 'perseus_hi_galfa_cube_regrid_planckres_noise.fits'
    if not path.isfile(hi_dir + noise_cube_filename):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=velocity_axis,
                velocity_noise_range=[90,110], header=h, Tsys=30.,
                filename=hi_dir + noise_cube_filename)
    else:
        noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    # define core properties
    with open(core_dir + 'perseus_core_properties.txt', 'r') as f:
        cores = json.load(f)

    # Remove cores
    cores_to_remove = []
    for core in cores:
        if core not in cores_to_keep:
            cores_to_remove.append(core)
    for core_to_remove in cores_to_remove:
        del cores[core_to_remove]

    cores = convert_core_coordinates(cores, h)
    cores = load_ds9_core_region(cores,
                            filename_base = region_dir + \
                                            'perseus_av_poly_cores',
                            header = av_header)

    # Load cloud division regions from ds9
    properties = load_ds9_region(properties,
                            filename=multicloud_region_dir + \
                                     'multicloud_divisions.reg',
                            header=av_header)

    region_vertices = properties['regions']['perseus']['poly_verts']['pixel']

    # block off region
    region_mask = np.logical_not(myg.get_polygon_mask(av_data_planck,
                                                      region_vertices))

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
    alphaG_list = []
    alphaG_error_list = []
    phi_g_list = []
    phi_g_error_list = []
    Z_list = []
    Z_error_list = []
    chisq_list = []
    p_value_list = []
    core_name_list = []
    co_image_list = []
    hi_vel_range_list = []
    hi_vel_range_likelihood_list = []
    hi_sd_fit_list = []

    if clobber:
    #if 0:
        hi_data_orig = np.copy(hi_data)
        noise_cube_orig = np.copy(noise_cube)
        av_data_planck_orig = np.copy(av_data_planck)
        av_error_data_planck_orig = np.copy(av_error_data_planck)

        for core in np.sort(cores.keys()):
            print('\nCalculating for core %s' % core)

            if box_method == 'ds9':
                indices = myg.get_polygon_mask(av_data_planck_orig,
                        cores[core]['poly_verts']['pixel'])
            elif box_method == 'av_gradient':
                indices = myg.get_polygon_mask(av_data_planck_orig,
                    cores[core]['{0:s}_vertices_rotated'.format(region_type)])
            else:
                raise ValueError('Method for boxes is either ds9 or ' + \
                                 'av_gradient')

            indices = indices.astype('bool')
            mask = ~indices

            mask += region_mask

            if 0:
                import matplotlib.pyplot as plt
                plt.imshow(np.ma.array(av_data_planck_orig, mask=mask),
                           origin='lower left')
                plt.show()

            # Get only the relevant pixels to decrease computation time
            hi_data = np.copy(hi_data_orig[:, ~mask])
            noise_cube = np.copy(noise_cube_orig[:, ~mask])
            av_data_planck = np.copy(av_data_planck_orig[~mask])
            av_error_data_planck = np.copy(av_error_data_planck_orig[~mask])

            # -----------------------------------------------------------------
            # Perform analysis on cores, including fitting the sternberg
            # model.  If calc_errors is True then a monte carlo is run by
            # adding noise to AV and HI and
            # refitting.
            # -----------------------------------------------------------------

            if calc_errors:
                images, hi_vel_range, params = \
                        run_analysis(hi_cube=hi_data,
                                     hi_noise_cube=noise_cube,
                                     hi_velocity_axis=velocity_axis,
                                     hi_header=h,
                                     dgr=dgr,
                                     dgr_error=dgr_error,
                                     intercept=intercept,
                                     vel_center=vel_center,
                                     av_image=av_data_planck,
                                     av_image_error=av_error_data_planck,
                                     core_dict=cores[core],
                                     results_figure_name=figure_dir + \
                                             'monte_carlo_results/' + \
                                             'perseus_%s' % core,
                                     properties=properties,
                                     results_filename=results_filename + core,
                                     )
            else:
                images, params = run_analysis(hi_cube=hi_data,
                                     hi_noise_cube=noise_cube,
                                     hi_velocity_axis=velocity_axis,
                                     hi_header=h,
                                     dgr=dgr,
                                     dgr_error=dgr_error,
                                     intercept=intercept,
                                     av_image=av_data_planck,
                                     av_image_error=av_error_data_planck,
                                     hi_vel_range=hi_vel_range,
                                     parameter_vary=[vary_alphaG, vary_Z])

            rh2_fit, h_sd_fit, f_H2, f_HI = calc_sternberg(params=params[:2],
                                                h_sd_extent=h_sd_fit_range,
                                                return_fractions=True)

            hi_sd_fit = f_HI * h_sd_fit
            alphaG = params[0]
            alphaG_error = params[2]
            Z = params[1]
            Z_error = params[3]
            phi_g = params[4]
            phi_g_error = params[5]

            # Calculate T_cnm from sternberg et al. (2009) Eq 19
            T_cnm = calc_T_cnm(alphaG, Z=Z)
            T_cnm_error = []
            T_cnm_error.append(\
                    T_cnm - calc_T_cnm(alphaG + alphaG_error[0], Z=Z))
            T_cnm_error.append(\
                    T_cnm - calc_T_cnm(alphaG + alphaG_error[1], Z=Z))

            cores[core]['hi_sd_fit'] = hi_sd_fit.tolist()
            cores[core]['rh2'] = images['rh2'].tolist()
            cores[core]['rh2_error'] = images['rh2_error'].tolist()
            cores[core]['hi_sd'] = images['hi_sd'].tolist()
            cores[core]['hi_sd_error'] = images['hi_sd_error'].tolist()
            cores[core]['av'] = images['av'].tolist()
            cores[core]['av_error'] = images['av_error'].tolist()
            cores[core]['h_sd'] = images['h_sd'].tolist()
            cores[core]['h_sd_error'] = images['h_sd_error'].tolist()
            cores[core]['alphaG'] = alphaG
            cores[core]['alphaG_error'] = alphaG_error
            cores[core]['T_cnm'] = T_cnm
            cores[core]['T_cnm_error'] = T_cnm_error
            cores[core]['Z'] = Z
            cores[core]['Z_error'] = Z_error
            cores[core]['phi_g'] = phi_g
            cores[core]['phi_g_error'] = phi_g_error
            cores[core]['rh2_fit'] = rh2_fit.tolist()
            cores[core]['h_sd_fit'] = h_sd_fit.tolist()
            cores[core]['f_H2'] = f_H2.tolist()
            cores[core]['f_HI_fit'] = f_HI.tolist()

        with open(core_dir + 'perseus_core_properties.txt', 'w') as f:
            json.dump(cores, f)

    for core in cores:
        if core in cores_to_keep and \
           len(cores[core]['hi_sd']) > 0:

            # append to the lists
            hi_sd_image_list.append(np.array(cores[core]['hi_sd']))
            hi_sd_image_error_list.append(np.array(cores[core]['hi_sd_error']))
            h_sd_image_list.append(np.array(cores[core]['h_sd']))
            h_sd_image_error_list.append(np.array(cores[core]['h_sd_error']))
            av_image_list.append(np.array(cores[core]['av']))
            av_image_error_list.append(np.array(cores[core]['av_error']))
            rh2_image_list.append(np.array(cores[core]['rh2']))
            rh2_image_error_list.append(np.array(cores[core]['rh2_error']))
            rh2_fit_list.append(np.array(cores[core]['rh2_fit']))
            h_sd_fit_list.append(np.array(cores[core]['h_sd_fit']))
            hi_sd_fit_list.append(np.array(cores[core]['hi_sd_fit']))
            alphaG_list.append(cores[core]['alphaG'])
            alphaG_error_list.append(cores[core]['alphaG_error'])
            Z_list.append(cores[core]['Z'])
            Z_error_list.append(cores[core]['Z_error'])
            phi_g_list.append(cores[core]['phi_g'])
            phi_g_error_list.append(cores[core]['phi_g_error'])
            core_name_list.append(core)

    # Create the figures!
    # -------------------
    print('\nCreating figures...')

    figure_types = ['png',]
    if write_pdf_figures:
        figure_types.append('pdf')

    for figure_type in figure_types:
        fig_name_rh2 = 'perseus_rh2_vs_hsd_panels_{0:s}'.format(av_data_type)
        print('\nWriting Rh2 figures to\n' + fig_name_rh2)

        plot_rh2_vs_h_grid(rh2_image_list,
                h_sd_image_list,
                rh2_error_images = rh2_image_error_list,
                h_sd_error_images = h_sd_image_error_list,
                rh2_fits = rh2_fit_list,
                h_sd_fits = h_sd_fit_list,
                limits = [0, 80, 10**-3, 10**2],
                savedir = figure_dir + 'panel_cores/',
                scale = ('linear', 'log'),
                filename=fig_name_rh2 + '_linear.{0:s}'.format(figure_type),
                core_names=core_name_list,
                alphaG_list=alphaG_list,
                alphaG_error_list=alphaG_error_list,
                phi_g_list=phi_g_list,
                phi_g_error_list=phi_g_error_list,
                Z_list=Z_list,
                Z_error_list=Z_error_list,
                show = False)

        plot_rh2_vs_h_grid(rh2_image_list,
                h_sd_image_list,
                rh2_error_images = rh2_image_error_list,
                h_sd_error_images = h_sd_image_error_list,
                rh2_fits = rh2_fit_list,
                h_sd_fits = h_sd_fit_list,
                limits = [1, 200, 10**-3, 10**2],
                savedir = figure_dir + 'panel_cores/',
                scale = ('log', 'log'),
                filename=fig_name_rh2 + '_log.{0:s}'.format(figure_type),
                core_names=core_name_list,
                alphaG_list=alphaG_list,
                alphaG_error_list=alphaG_error_list,
                phi_g_list=phi_g_list,
                phi_g_error_list=phi_g_error_list,
                Z_list=Z_list,
                Z_error_list=Z_error_list,
                show = False)

        # Calif limits = [0, 80, 0, 8]
        # taur limits = [0, 80, -1.5, 6.5]
        # Pers limits = [0, 80, 1, 8],
        fig_name_hivsh = 'perseus_hi_vs_h_panels_{0:s}'.format(av_data_type)
        print('\nWriting HI vs H figures to\n' + fig_name_hivsh)

        plot_hi_vs_h_grid(hi_sd_image_list,
                h_sd_image_list,
                hi_sd_error_images=hi_sd_image_error_list,
                h_sd_error_images=h_sd_image_error_list,
                hi_fits=hi_sd_fit_list,
                h_sd_fits=h_sd_fit_list,
                limits=[-9, 80, -1.5, 6.5],
                savedir=figure_dir + 'panel_cores/',
                scale=('linear', 'linear'),
                filename=fig_name_hivsh + '_linear.{0:s}'.format(figure_type),
                core_names=core_name_list,
                alphaG_list=alphaG_list,
                alphaG_error_list=alphaG_error_list,
                phi_g_list=phi_g_list,
                phi_g_error_list=phi_g_error_list,
                Z_list=Z_list,
                Z_error_list=Z_error_list,
                show = False)

        plot_hi_vs_h_grid(hi_sd_image_list,
                h_sd_image_list,
                hi_sd_error_images=hi_sd_image_error_list,
                h_sd_error_images=h_sd_image_error_list,
                hi_fits=hi_sd_fit_list,
                h_sd_fits=h_sd_fit_list,
                limits=[1, 100, 1, 100],
                savedir=figure_dir + 'panel_cores/',
                scale=('log', 'log'),
                filename=fig_name_hivsh + '_log.{0:s}'.format(figure_type),
                core_names=core_name_list,
                alphaG_list=alphaG_list,
                alphaG_error_list=alphaG_error_list,
                phi_g_list=phi_g_list,
                phi_g_error_list=phi_g_error_list,
                Z_list=Z_list,
                Z_error_list=Z_error_list,
                show = False)

        fig_name_hivsav = 'perseus_hi_vs_av_panels_{0:s}'.format(av_data_type)
        print('\nWriting HI vs H figures to\n' + fig_name_hivsav)

        plot_hi_vs_av_grid(hi_sd_image_list,
                av_image_list,
                hi_error_images = hi_sd_image_error_list,
                av_error_images = h_sd_image_error_list,
                #limits = [10**-1, 10**2, 10**0, 10**2],
                limits = [0, 50, 1, 8],
                savedir = figure_dir + 'panel_cores/',
                scale = ('linear', 'linear'),
                filename=fig_name_hivsav + '_linear.{0:s}'.format(figure_type),
                core_names=core_name_list,
                #title = r'$\Sigma_{\rm HI}$ vs. $\Sigma_{\rm H}$'\
                #        + ' of perseus Cores',
                show = False)

if __name__ == '__main__':
    main(av_data_type='planck', region=None)
    main(av_data_type='k09')



