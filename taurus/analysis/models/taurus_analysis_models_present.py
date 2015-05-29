#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the taurus molecular cloud.
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
              'figure.dpi': 600,
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

def plot_hi_vs_h_grid(cores, limits=None, fit=True, savedir='./',
        cores_to_plot=None, filename=None, show=True, scale=('linear',
            'linear'), title='', plot_fit=True, plot_errors=False):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from astroML.plotting import scatter_contour

    if len(cores_to_plot) > 1:
        n = int(np.ceil(len(cores_to_plot)**0.5))
        if n**2 - n > len(cores_to_plot):
            nrows = n - 1
            ncols = 2
            y_scaling = 1.0 - 1.0 / n
        else:
            nrows, ncols = n, n
            y_scaling = 1.0

        n = len(cores_to_plot)
        ncols = 4
        nrows = n / ncols + 1
        y_scaling = nrows / float(ncols)

        figsize = (7.5, 7.5*y_scaling)
        font_scale = 9
    else:
        nrows, ncols = 1, 1
        n = 1
        figsize = (4, 4)
        font_scale = 18

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 3)]
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
              #'figure.figsize': (3.6, 3.6*y_scaling),
              'figure.figsize': figsize,
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

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(nrows, ncols),
                 ngrids=n,
                 axes_pad=0.,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    # Cycle through lists
    for i, core in enumerate(cores_to_plot):
        for key in cores[core]:
            if type(cores[core][key]) is list:
                cores[core][key] = np.array(cores[core][key])

        hi_sd = cores[core]['hi_sd']
        hi_sd_error = cores[core]['hi_sd_error']
        h_sd = cores[core]['h_sd']
        h_sd_error = cores[core]['h_sd_error']

        # Load parameters
        alphaG = cores[core]['sternberg_results']['alphaG']
        alphaG_error = cores[core]['sternberg_results']['alphaG_error']
        phi_cnm = cores[core]['krumholz_results']['phi_cnm']
        phi_cnm_error = cores[core]['krumholz_results']['phi_cnm_error']
        hi_sd_fit_sternberg = cores[core]['sternberg_results']['hi_sd_fit']
        hi_sd_fit_krumholz = cores[core]['krumholz_results']['hi_sd_fit']
        h_sd_fit = cores[core]['krumholz_results']['h_sd_fit']

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

        #ax.set_xticks([0, 40, 80, 120])


        if plot_errors:
            ax.axvline(np.nanmean(h_sd_error_nonans.ravel()),
                       color='k',
                       linestyle='--',
                       linewidth=2,
                       alpha=0.5,
                       zorder=0)
            ax.axhline(np.nanmean(hi_sd_error_nonans.ravel()),
                       color='k',
                       linestyle='--',
                       linewidth=2,
                       alpha=0.5,
                       zorder=1)

        if 0:
            l1 = ax.errorbar(h_sd_nonans.ravel(),
                    hi_sd_nonans.ravel(),
                    xerr=(h_sd_error_nonans.ravel()),
                    yerr=(hi_sd_error_nonans.ravel()),
                    alpha=0.1,
                    #color='k',
                    marker='^',ecolor='k',linestyle='None',
                    markersize=3
                    )
        elif 1:
            l1 = scatter_contour(h_sd_nonans.ravel(),
                                 hi_sd_nonans.ravel(),
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

        if plot_fit:
            if 0:
                l2 = ax.plot(h_sd_fit, hi_sd_fit_sternberg,
                        label='S+14',
                        color=color_cycle[1],
                        alpha=0.75,
                        )
                if i == 0:
                    ax.legend(loc='upper left')

            l3 = ax.plot(h_sd_fit, hi_sd_fit_krumholz,
                    linestyle='-',
                    linewidth=3,
                    label='K+09',
                    color=color_cycle[2],
                    alpha=0.75
                    )


        # Annotations
        anno_xpos = 0.95

        alphaG_text = r'\noindent$\alpha G$ =' + \
                       r' %.2f' % (alphaG) + \
                       r'$^{+%.2f}_{-%.2f}$ \\' % (alphaG_error[0],
                                                   alphaG_error[1])
        phi_cnm_text = r'\noindent$\phi_{\rm CNM}$ =' + \
                       r' %.2f' % (phi_cnm) + \
                       r'$^{+%.2f}_{-%.2f}$ \\' % (phi_cnm_error[0],
                                                   phi_cnm_error[1])


        if 0:
            ax.annotate(alphaG_text + phi_cnm_text,
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

        ax.annotate(core,
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
        ax.set_xlabel(r'$\Sigma_{\rm H\,I}$ + $\Sigma_{\rm H_2}$ ' + \
                       '[M$_\odot$ pc$^{-2}$]',)
        ax.set_ylabel(r'$\Sigma_{\rm H\,I}$ [M$_\odot$ pc$^{-2}$]',)

        if 'log' not in scale:
            ax.locator_params(nbins=6)

    if 0:
        fig.legend(l2 + l3,
                   ('S+14', 'K+09'),
                   loc='upper center',
                   #loc=3,
                   ncol=2,
                   numpoints=6,
                   bbox_transform = plt.gcf().transFigure,
                   #bbox_transform = imagegrid.transFigure,
                   #bbox_to_anchor=(0., 0.95, 1.05, -0.0),
                   mode="expand",
                   bbox_to_anchor=(0.3, 1.0, 0.5, 0.1),
                   borderaxespad=0.0)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight', dpi=800)
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

def run_analysis(hi_cube=None,
                 hi_noise_cube=None,
                 hi_velocity_axis=None,
                 hi_header=None,
                 dgr=None,
                 dgr_error=None,
                 intercept=None,
                 av_image=None,
                 av_image_error=None,
                 vel_center=None,
                 hi_vel_range=None,
                 hi_vel_range_error=None,
                 verbose=False,
                 core_dict=None,
                 results_figure_name=None,
                 properties=None,
                 results_filename=None,
                 sternberg_params=None,
                 krumholz_params=None):

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


    verbose = False

    clobber = sternberg_params['clobber']
    N_monte_carlo_runs = sternberg_params['N_monte_carlo_runs']
    if N_monte_carlo_runs < 1:
        raise ValueError('N_monte_carlo_runs must be >= 1')

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
    sternberg_results = sternberg_params.copy()
    krumholz_results = krumholz_params.copy()
    sternberg_results.update({'alphaG fits' : np.empty((N_monte_carlo_runs)),
                         'Z fits' : np.empty((N_monte_carlo_runs)),
                         'phi_g fits' : np.empty((N_monte_carlo_runs)),})
    krumholz_results.update({'phi_cnm fits' : np.empty((N_monte_carlo_runs)),
                        'Z fits' : np.empty((N_monte_carlo_runs)),
                        'phi_mol fits' : np.empty((N_monte_carlo_runs)),})
    results_dict = {'nhi errors' : np.empty((N_monte_carlo_runs))}

    hi_vel_range_list = np.empty((N_monte_carlo_runs, 2))
    width_list = np.empty((N_monte_carlo_runs))
    dgr_list = np.empty((N_monte_carlo_runs))
    intercept_list = np.empty((N_monte_carlo_runs))

    # Get standard errors on images + the HI velocity range
    hi_error = np.median(hi_noise_cube)
    av_error = np.median(av_image_error)

    # Load the calculate marginalized PDFs for each parameter
    width_likelihoods = np.asarray(properties['width_likelihood'])
    dgr_likelihoods = np.asarray(properties['dgr_likelihood'])
    intercept_likelihoods = np.asarray(properties['intercept_likelihood'])

    # Load center velocity
    vel_center = properties['hi_velocity_center']['value']

    # Load the full 3D likelihood space
    likelihoods = np.asarray(properties['likelihoods'])

    # Load the grids
    vel_widths = np.asarray(properties['vel_widths'])
    dgrs = np.asarray(properties['dgrs'])
    intercepts = np.asarray(properties['intercepts'])

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

                width_list[i] = vel_width_random
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

            # Fit to sternberg model
            alphaG, Z, phi_g = fit_sternberg(h_sd_ravel,
                                      rh2_ravel,
                                      guesses=sternberg_params['guesses'],
                                      vary=sternberg_params['param_vary'],
                                      rh2_error=rh2_error_ravel,
                                      verbose=verbose)

            # keep results
            sternberg_results['alphaG fits'][i] = alphaG
            sternberg_results['phi_g fits'][i] = phi_g
            sternberg_results['Z fits'][i] = Z

            # Fit to krumholz model
            phi_cnm, Z, phi_mol = fit_krumholz(h_sd_ravel,
                                      rh2_ravel,
                                      guesses=krumholz_params['guesses'],
                                      vary=krumholz_params['param_vary'],
                                      rh2_error=rh2_error_ravel,
                                      verbose=verbose)

            # keep results
            krumholz_results['phi_cnm fits'][i] = phi_cnm
            krumholz_results['phi_mol fits'][i] = phi_mol
            krumholz_results['Z fits'][i] = Z

            # see eq 6 of sternberg+09
            # alphaG is the number density of the CNM over the minimum number
            # density required for pressure balance
            # the lower alphaG values than for taurus mean that taurus
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


        sternberg_results = \
                calc_MC_errors(sternberg_results)

        sternberg_results = \
                analyze_sternberg_model(sternberg_results)

        krumholz_results = \
                calc_MC_errors(krumholz_results)

        krumholz_results = \
                analyze_krumholz_model(krumholz_results)

        # Update results dictionary to include fitted results
        for key in sternberg_results:
            if type(sternberg_results[key]) is np.ndarray:
                sternberg_results[key] = sternberg_results[key].tolist()
        for key in krumholz_results:
            if type(krumholz_results[key]) is np.ndarray:
                    krumholz_results[key] = krumholz_results[key].tolist()

        results_dict.update({'sternberg_results': sternberg_results,
                             'krumholz_results': krumholz_results,
                             'dgrs': dgr_list,
                             'intercepts': intercept_list,
                             'widths': width_list})

    # If not performing the MC, read a previously written MC file
    elif not perform_mc:
        print('Reading monte carlo results file:')
        print(results_filename)

        with open(results_filename, 'r') as f:
            results_dict = json.load(f)

    # Derive images
    # -------------
    hi_vel_range_sample = (np.median(hi_vel_range_list[:, 0]),
                           np.median(hi_vel_range_list[:, 1]))

    nhi_error = np.std(results_dict['nhi errors'])

    dgr = properties['dust2gas_ratio_max']['value']
    intercept = properties['intercept_max']['value']
    #intercept_error = properties['intercept_error']['value']
    #dgr_error = properties['dust2gas_ratio_error']['value']
    vel_center = properties['hi_velocity_center']['value']
    vel_width = properties['hi_velocity_width_max']['value']
    hi_vel_range = (vel_center - vel_width/2.0, vel_center + vel_width/2.0)
    #hi_vel_range = properties['hi_velocity_range']


    # Get values of best-fit model parameters
    #max_loc = np.where(likelihoods == np.max(likelihoods))
    #vel_width = vel_widths[max_loc[0][0]]
    #hi_velocity_range =


    images = derive_images(hi_cube=hi_cube,
                           hi_velocity_axis=hi_velocity_axis,
                           hi_noise_cube=hi_noise_cube,
                           hi_vel_range=hi_vel_range,
                           hi_header=hi_header,
                           dgr=dgr,
                           intercept=intercept,
                           av_image=av_image,
                           av_image_error=av_image_error,
                           )
    if 0:
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

        # Fit to sternberg model
        alphaG, Z, phi_g = fit_sternberg(h_sd_ravel,
                                  rh2_ravel,
                                  guesses=sternberg_params['guesses'],
                                  vary=sternberg_params['param_vary'],
                                  rh2_error=rh2_error_ravel,
                                  verbose=verbose)

        # keep results
        sternberg_results['alphaG fits'][i] = alphaG
        sternberg_results['phi_g fits'][i] = phi_g
        sternberg_results['Z fits'][i] = Z

        # Fit to krumholz model
        phi_cnm, Z, phi_mol = fit_krumholz(h_sd_ravel,
                                  rh2_ravel,
                                  guesses=krumholz_params['guesses'],
                                  vary=krumholz_params['param_vary'],
                                  rh2_error=rh2_error_ravel,
                                  verbose=verbose)

        # keep results
        krumholz_results['phi_cnm fits'][i] = phi_cnm
        krumholz_results['phi_mol fits'][i] = phi_mol
        krumholz_results['Z fits'][i] = Z

    results_dict.update(core_dict)
    results_dict.update({'rh2' : images['rh2'],
                         'rh2_error' : images['rh2_error'],
                         'hi_sd' : images['hi_sd'],
                         'hi_sd_error' : images['hi_sd_error'],
                         'av' : images['av'],
                         'av_error' : images['av_error'],
                         'h_sd' : images['h_sd'],
                         'h_sd_error' : images['h_sd_error'],
                         'sternberg_results' : sternberg_results,
                         'krumholz_results' : krumholz_results,
                         })

    # Write the MC results in JSON human readable txt format
    if write_mc:
        for key in results_dict:
            if type(results_dict[key]) is np.ndarray:
                    results_dict[key] = results_dict[key].tolist()

        with open(results_filename, 'w') as f:
            json.dump(results_dict, f)

    # Print results
    print('\n\tMonte Carlo Results:')
    print('\tSternberg Fit Parameters')
    for parameter in results_dict['sternberg_results']['parameters']:
        confint = results_dict['sternberg_results'][parameter + '_confint']
        print('\t\t{0:s} = {1:.2f} +{2:.2f}/-{3:.2f}'.format(parameter,
                                                           *confint))
    print('\tKrumholz Fit Parameters')
    for parameter in results_dict['krumholz_results']['parameters']:
        confint = results_dict['krumholz_results'][parameter + '_confint']
        print('\t\t{0:s} = {1:.2f} +{2:.2f}/-{3:.2f}'.format(parameter,
                                                           *confint))

    return results_dict

def calc_MC_errors(results_dict, error_method='edges', alpha=0.05,):

    from mystats import calc_symmetric_error
    from scikits.bootstrap import ci
    import numpy as np

    for i, parameter in enumerate(results_dict['parameters']):
        # If there is a distribution of the parameter, then find the
        # confidence interval
        if results_dict['param_vary'][i]:
            param_fits = results_dict[parameter + ' fits']

            # Histogram will act as distribution of parameter values
            counts, bins = np.histogram(param_fits,
                                        bins=np.linspace(np.nanmin(param_fits),
                                                         np.nanmax(param_fits),
                                                         1000))

            # Calculate confidence interval by taking the PDF weighted mean as
            # value, where errors are calculated by moving vertical threshold
            # from edges towards mean until alpha/2 *100 % of the area of the
            # PDF is included within the threshold
            #alphaG_confint = calc_symmetric_error(bins[:-1], counts,
            #        alpha=alpha)
            results_dict[parameter + '_confint'] = \
                 calc_symmetric_error(param_fits,
                                      alpha=results_dict['alpha'])

        else:
            results_dict[parameter + '_confint'] = \
                    (results_dict[parameter + ' fits'][0], 0.0, 0.0)

        # Split up errors and best estimate
        results_dict[parameter] = results_dict[parameter + '_confint'][0]
        results_dict[parameter + '_error'] = \
                results_dict[parameter + '_confint'][1:]

    return results_dict


''' Fitting Functions
'''

def calc_krumholz(params=[10.0, 1.0, 10.0], h_sd_extent=(0.001, 500),
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
                                           phi_cnm=params[0],
                                           Z=params[1],
                                           phi_mol=params[2],
                                           return_fractions=True)

    output = [rh2_fit, h_sd]

    if return_fractions:
        output.append(f_H2)
        output.append(f_HI)

    return output

def fit_krumholz(h_sd, rh2, guesses=[10.0, 1.0, 10.0], rh2_error=None,
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
    from myscience import krumholz09 as k09

    def chisq(params, h_sd, rh2, rh2_error):
        phi_cnm = params['phi_cnm'].value
        phi_mol = params['phi_mol'].value
        Z = params['Z'].value

        rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol)

        chisq = np.sum((rh2 - rh2_model)**2 / rh2_error**2)

        return chisq

    # Set parameter limits and initial guesses
    params = Parameters()
    params.add('phi_cnm',
               value=guesses[0],
               min=0.5,
               max=100,
               vary=vary[0])
    params.add('phi_mol',
               value=guesses[2],
               min=1,
               max=20,
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
                      method='lbfgsb')

    rh2_fit_params = (params['phi_cnm'].value, params['Z'].value,
            params['phi_mol'].value)

    return rh2_fit_params

def analyze_krumholz_model(krumholz_results):

    ''' Calculates various properties of Krumholz model from fitted phi_cnm
    values.

    '''

    from myscience.krumholz09 import calc_T_cnm

    # Calculate T_cnm from Krumholz et al. (2009) Eq 19
    phi_cnm, phi_cnm_error, Z = \
        [krumholz_results[key] for key in ('phi_cnm', 'phi_cnm_error', 'Z')]
    T_cnm = calc_T_cnm(phi_cnm, Z=Z)
    T_cnm_error = []
    T_cnm_error.append(\
            T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[0], Z=Z))
    T_cnm_error.append(\
            T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[1], Z=Z))

    krumholz_results.update({'T_cnm': T_cnm,
                             'T_cnm_error': T_cnm_error})

    # Get fitted surface density ratios...
    params = [krumholz_results[param] for param in \
              krumholz_results['parameters']]

    rh2_fit, h_sd_fit, f_H2, f_HI = \
           calc_krumholz(params=params,
                          h_sd_extent=krumholz_results['h_sd_fit_range'],
                          return_fractions=True)

    krumholz_results.update({'rh2_fit': rh2_fit,
                             'h_sd_fit' : h_sd_fit,
                             'hi_sd_fit' : f_HI * h_sd_fit,
                             'f_H2' : f_H2,
                             'f_HI' : f_HI,})

    return krumholz_results

def calc_sternberg(params=[1.5, 1.0, 1.0], h_sd_extent=(0.001, 500),
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
                                           alphaG=params[0],
                                           Z=params[1],
                                           phi_g=params[2],
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

def analyze_sternberg_model(sternberg_results):

    ''' Calculates various properties of Krumholz model from fitted phi_cnm
    values.

    '''

    # Get fitted surface density ratios...
    params = [sternberg_results[param] for param in \
              sternberg_results['parameters']]

    rh2_fit, h_sd_fit, f_H2, f_HI = \
           calc_sternberg(params=params,
                          h_sd_extent=sternberg_results['h_sd_fit_range'],
                          return_fractions=True)

    sternberg_results.update({'rh2_fit': rh2_fit,
                              'h_sd_fit' : h_sd_fit,
                              'hi_sd_fit' : (1 - f_H2) * h_sd_fit,
                              'f_H2' : f_H2,
                              'f_HI' : f_HI,})

    return sternberg_results

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

def load_ds9_core_region(cores, filename_base = 'taurus_av_boxes_',
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

def main(verbose=True, av_data_type='planck', regions=None):

    '''

    This script requires analysis output from three other scripts.
    Script                                          Purpose
    ---------------------------------------------------------------------------
    hi/taurus_analysis_core_properties.py           Defines core positions
    av/taurus_analysis_av_derive_core_boxes.py      Calculates box sizes
    hi/taurus_analysis_hi_av_core_likelihoods.py    Calculates HI velocity
                                                        range
    av/taurus_analysis_av_load_regions.py

    '''

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
    hi_av_likelihoodelation = True
    reproduce_lee12 = True

    # Sternberg Parameters
    # --------------------
    N_monte_carlo_runs = 1000 # Number of monte carlo runs
    vary_alphaG = True # Vary alphaG in S+14 fit?
    vary_Z = False # Vary metallicity in S+14 fit?
    vary_phi_g = False # Vary phi_g in S+14 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    guesses=(1.0, 1.0, 1.0) # Guesses for (alphaG, Z, phi_g)
    h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

    clobber = 0 # perform MC and write over current results?

    # Monte carlo results file bases
    results_filename = '/d/bip3/ezbc/taurus/data/python_output/' + \
                       'monte_carlo_results/taurus_mc_results_sternberg' + \
                       av_data_type + '_'

    sternberg_params = {}
    sternberg_params['N_monte_carlo_runs'] = N_monte_carlo_runs
    sternberg_params['param_vary'] = [vary_alphaG, vary_Z, vary_phi_g]
    sternberg_params['error_method'] = error_method
    sternberg_params['alpha'] = alpha
    sternberg_params['guesses'] = guesses
    sternberg_params['clobber'] = clobber
    sternberg_params['h_sd_fit_range'] = h_sd_fit_range
    sternberg_params['results_filename'] = results_filename
    sternberg_params['parameters'] = ['alphaG', 'Z', 'phi_g']

    # Krumholz Parameters
    # --------------------
    vary_phi_cnm = True # Vary phi_cnm in K+09 fit?
    vary_Z = False # Vary metallicity in K+09 fit?
    vary_phi_mol = False # Vary phi_mol in K+09 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    guesses=(10.0, 1.0, 10.0) # Guesses for (phi_cnm, Z, phi_mol)
    h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

    # Monte carlo results file bases
    results_filename = '/d/bip3/ezbc/taurus/data/python_output/' + \
                       'monte_carlo_results/taurus_mc_results_krumholz' + \
                       av_data_type + '_'

    krumholz_params = {}
    krumholz_params['N_monte_carlo_runs'] = N_monte_carlo_runs
    krumholz_params['param_vary'] = [vary_phi_cnm, vary_Z, vary_phi_mol]
    krumholz_params['error_method'] = error_method
    krumholz_params['alpha'] = alpha
    krumholz_params['guesses'] = guesses
    krumholz_params['clobber'] = clobber
    krumholz_params['h_sd_fit_range'] = h_sd_fit_range
    krumholz_params['results_filename'] = results_filename
    krumholz_params['parameters'] = ['phi_cnm', 'Z', 'phi_mol']

    # Universal properties
    # --------------------
    # global property filename
    global_property_filename = 'taurus_global_properties'

    if 0:
        # Name correct region of cloud
        if region == 1:
            region_name = 'taurus1'
        elif region == 2:
            region_name = 'taurus2'
        else:
            region_name = 'taurus'


    # Which cores to include in analysis?
    cores_to_keep = [# taur
                     'L1495',
                     #'L1495A',
                     'B213',
                     'L1498',
                     'B215',
                     'B18',
                     #'B217',
                     #'B220-1',
                     'B220-2',
                     #'L1521',
                     'L1524',
                     #'L1527-1',
                     'L1527-2',
                     # Calif
                     #'L1536',
                     'L1483',
                     'L1478',
                     'L1456',
                     'NGC1579',
                     'L1545',
                     'L1517',
                     'L1512',
                     'L1523',
                     'L1512',
                     # Pers
                     'B5',
                     'IC348',
                     'B1E',
                     'B1',
                     'NGC1333',
                     'L1482'
                     ]

    # Regions
    # Options are 'ds9' or 'av_gradient'
    box_method = 'ds9'
    region_type = 'wedge' # Shape of core region, box or wedge

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
    multicloud_region_dir = \
            '/d/bip3/ezbc/multicloud/data/python_output/'

    # =========================================================================
    # Begin Analysis
    # =========================================================================

    # load Planck Av and GALFA HI images, on same grid
    # ------------------------------------------------
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
    global_properties_list = []
    for region_name in regions:
        filename = \
                global_property_filename.replace('taurus', region_name)

        with open(property_dir + \
                  filename + '_' + av_data_type + \
                  '_scaled.txt', 'r') as f:
            properties = json.load(f)

        global_properties_list.append(properties)

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
    cores = load_ds9_core_region(cores,
                            filename_base = region_dir + \
                                            'taurus_av_poly_cores',
                            header = av_header)

    # Load cloud division regions from ds9
    properties = load_ds9_region(properties,
                            filename=multicloud_region_dir + \
                                     'multicloud_divisions.reg',
                            header=av_header)

    region_vertices = properties['regions']['taurus']['poly_verts']['pixel']

    # block off region
    region_mask = np.logical_not(myg.get_polygon_mask(av_data_planck,
                                                      region_vertices))

    # Trim down cores to keep list to include only cores for cloud
    cores_to_keep_old = list(cores_to_keep)
    for core in cores_to_keep_old:
        if core not in cores:
            cores_to_keep.remove(core)

    # Set up lists
    sternberg_results = {}
    krumholz_results = {}

    if clobber:
        hi_data_orig = np.copy(hi_data)
        noise_cube_orig = np.copy(noise_cube)
        av_data_planck_orig = np.copy(av_data_planck)
        av_error_data_planck_orig = np.copy(av_error_data_planck)

        # Initialize dictionary to store MC results for each core
        #results = {}

        for core in cores_to_keep:
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

            #mask += region_mask

            # find which region the core is in
            for global_properties in global_properties_list:
                center = cores[core]['center_pixel']
                if myg.point_in_polygon(center, region_vertices):
                    properties = global_properties

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

            cores[core] = \
                    run_analysis(hi_cube=hi_data,
                                 hi_noise_cube=noise_cube,
                                 hi_velocity_axis=velocity_axis,
                                 hi_header=h,
                                 av_image=av_data_planck,
                                 av_image_error=av_error_data_planck,
                                 core_dict=cores[core],
                                 sternberg_params=sternberg_params,
                                 krumholz_params=krumholz_params,
                                 results_figure_name=figure_dir + \
                                         'monte_carlo_results/' + \
                                         'taurus_%s' % core,
                                 properties=properties,
                                 results_filename=results_filename + core,
                                 )

        # Write results to cores
        for key in cores:
            if type(cores[key]) is np.ndarray:
                    cores[key] = cores[key].tolist()
        with open(core_dir + 'taurus_core_properties.txt', 'w') as f:
            json.dump(cores, f)

    # Create the figures!
    # -------------------
    print('\nCreating figures...')

    figure_types = ['png',]
    if write_pdf_figures:
        figure_types.append('pdf')



    for figure_type in figure_types:

        # Calif limits = [0, 80, 0, 8]
        # taur limits = [0, 80, -1.5, 6.5]
        # Pers limits = [0, 80, 1, 8],
        fig_name_hivsh = 'taurus_hi_vs_h_panels_{0:s}'.format(av_data_type)
        print('\nWriting HI vs H figures to\n' + fig_name_hivsh)

        plot_hi_vs_h_grid(cores,
                #limits=[-9, 80, -1.5, 6.5],
                limits=[-9, 159, -1.5, 30],
                savedir=figure_dir + 'panel_cores/',
                scale=('linear', 'linear'),
                filename=fig_name_hivsh + '_linear.{0:s}'.format(figure_type),
                cores_to_plot=cores_to_keep,
                show = False)

        plot_hi_vs_h_grid(cores,
                #limits=[-9, 80, -1.5, 6.5],
                limits=[-9, 159, -1.5, 20],
                savedir=figure_dir + 'panel_cores/',
                scale=('linear', 'linear'),
                filename=fig_name_hivsh + '_linear_B213_fit.{0:s}'.format(figure_type),
                cores_to_plot=('B213',),
                show = False)

        plot_hi_vs_h_grid(cores,
                #limits=[-9, 80, -1.5, 6.5],
                limits=[-9, 159, -1.5, 20],
                savedir=figure_dir + 'panel_cores/',
                scale=('linear', 'linear'),
                plot_errors=True,
                filename=fig_name_hivsh + '_linear_B213_fit_errors.{0:s}'.format(figure_type),
                cores_to_plot=('B213',),
                show = False)

        plot_hi_vs_h_grid(cores,
                #limits=[-9, 80, -1.5, 6.5],
                limits=[-9, 159, -1.5, 20],
                savedir=figure_dir + 'panel_cores/',
                scale=('linear', 'linear'),
                filename=fig_name_hivsh + '_linear_B213_nofit.{0:s}'.format(figure_type),
                cores_to_plot=('B213',),
                plot_fit=False,
                show = False)

if __name__ == '__main__':
    main(av_data_type='planck', regions=('taurus1', 'taurus2',))

