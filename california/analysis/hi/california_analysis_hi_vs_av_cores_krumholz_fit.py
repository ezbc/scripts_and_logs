#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the california molecular cloud.
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
    fontScale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
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
    fontScale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
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
        fig.suptitle(title, fontsize=fontScale*1.5)
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
    fontScale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
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
        fig.suptitle(title, fontsize=fontScale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_rh2_vs_h(rh2, h_sd, rh2_error=None, h_sd_error=None, rh2_fit = None,
        h_sd_fit = None, limits = None, fit = True, savedir = './', filename =
        None, show = True, scale = 'linear', title = ''):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib


    # Drop the NaNs from the images
    if type(rh2_error) is float:
        indices = np.where((rh2 == rh2) &\
                           (h_sd == h_sd)&\
                           (h_sd > 0) &\
                           (rh2 > -5))

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

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    fontScale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
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
    if 1:
        image = ax.errorbar(h_sd_nonans.ravel(),
                rh2_nonans.ravel(),
                xerr=(h_sd_error_nonans.ravel()),
                yerr=(rh2_error_nonans.ravel()),
                alpha=0.3,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=3
                )
    if rh2_fit is not None:
        ax.plot(h_sd_fit, rh2_fit,
                color = 'r')

    ax.set_xscale(scale)
    ax.set_yscale(scale)

    if limits is not None:
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)
    ax.set_ylabel(r'R$_{H2}$',)
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()

def plot_rh2_vs_h_grid(rh2_images, h_sd_images, rh2_error_images=None,
        h_sd_error_images=None, rh2_fits = None, h_sd_fits = None, limits =
        None, fit = True, savedir = './', filename = None, show = True, scale =
        'linear', title = '', core_names='', phi_cnm_list=None,
        chisq_list=None, p_value_list=None):

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
    fontScale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (10, 10),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    n = int(np.ceil(len(rh2_images)**0.5))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(n, n),
                 ngrids=len(rh2_images),
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True)

    # Cycle through lists
    for i in xrange(len(rh2_images)):
        rh2 = rh2_images[i]
        h_sd = h_sd_images[i]
        rh2_error = rh2_error_images[i]
        h_sd_error = h_sd_error_images[i]
        rh2_fit = rh2_fits[i]
        h_sd_fit = h_sd_fits[i]
        if phi_cnm_list is not None:
            phi_cnm = phi_cnm_list[i]
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

        if phi_cnm_list is not None:
            ax.annotate(r'$\phi$ = %.2f' % phi_cnm,
                    xytext=(0.6, 0.1),
                    xy=(0.6, 0.1),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    color='k'
                    )
        if chisq_list is not None:
            conf = (1.0 - p_value) * 100.
            ax.annotate(r'$\chi^2$ = %.2f' % \
                    chisq,
                    xytext=(0.48, 0.2),
                    xy=(0.48, 0.2),
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
        fig.suptitle(title, fontsize=fontScale*1.5)
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
    fontScale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
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
        fig.suptitle(title, fontsize=fontScale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()

def plot_co_spectrum_grid(vel_axis, co_spectrum_list,
        vel_range_list=None,
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
    fontScale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': fontScale,
              'axes.titlesize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
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
        hi_velocity_range = vel_range_list[i]

        # Create plot
        ax = imagegrid[i]

        image = ax.plot(vel_axis, co_spectrum,
                color='k',
                marker=None,
                drawstyle='steps-mid',
                markersize=3
                )

        if hi_velocity_range is not None:
            ax.axvspan(hi_velocity_range[0], hi_velocity_range[1],
                    color = 'r', alpha=0.3)

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

    if title is not None:
        fig.suptitle(title, fontsize=fontScale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
    if show:
        fig.show()

''' Calculations
'''

def calculate_nhi(cube=None, velocity_axis=None, velocity_range=[],
        return_nhi_error=True, noise_cube=None,
        velocity_noise_range=[90,100], Tsys=30., header=None,
        fits_filename=None, fits_error_filename=None, verbose=True):

    ''' Calculates an N(HI) image given a velocity range within which to
    include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to
    fits_filename : str
        If specified, and a header is provided, the nhi image will be written.
    header : pyfits.Header
        Header from cube.

    '''

    import numpy as np

    # Calculate NHI from cube if set
    if cube is not None and velocity_axis is not None:
        image = np.empty((cube.shape[1],
                          cube.shape[2]))
        image[:,:] = np.NaN
        indices = np.where((velocity_axis > velocity_range[0]) & \
                (velocity_axis < velocity_range[1]))[0]
        image[:,:] = cube[indices,:,:].sum(axis=0)
        # Calculate image error
        if return_nhi_error:
            image_error = np.empty((cube.shape[1],
                              cube.shape[2]))
            image_error[:,:] = np.NaN
            image_error[:,:] = (noise_cube[indices,:,:]**2).sum(axis=0)**0.5

    # NHI in units of 1e20 cm^-2
    nhi_image = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2
    nhi_image = image * 1.823e-2

    if fits_filename is not None and header is not None:
        if verbose:
            print('Writing N(HI) image to FITS file %s' % fits_filename)
        header['BUNIT'] = '1e20 cm^-2'
        header.remove('CDELT3')
        header.remove('CRVAL3')
        header.remove('CRPIX3')
        header.remove('CTYPE3')
        header.remove('NAXIS3')
        header['NAXIS'] = 2

        pf.writeto(fits_filename, image*1.823e-2, header = header, clobber =
                True, output_verify = 'fix')

    if fits_error_filename is not None and header is not None:
        if verbose:
            print('Writing N(HI) error image to FITS file %s' % fits_filename)

        pf.writeto(fits_error_filename, image_error * 1.823e-2, header =
                header, clobber = True, output_verify = 'fix')

    if return_nhi_error:
        nhi_image_error = np.ma.array(image_error,
                mask=np.isnan(image_error)) * 1.823e-2
        nhi_image_error = image_error * 1.823e-2
        return nhi_image, nhi_image_error
    else:
        return nhi_image

def calculate_noise_cube(cube=None, velocity_axis=None,
            velocity_noise_range=[-110,-90,90,110], header=None, Tsys=30.,
            filename=None):

    """ Calcuates noise envelopes for each pixel in a cube
    """

    import numpy as np
    import pyfits as pf

    noise_cube = np.zeros(cube.shape)
    for i in range(cube.shape[1]):
        for j in range(cube.shape[2]):
            profile = cube[:,i,j]
            noise = calculate_noise(profile, velocity_axis,
                    velocity_noise_range)
            #noise = 0.1 # measured in line free region
            noise_cube[:,i,j] = calculate_noise_scale(Tsys,
                    profile, noise=noise)

    if filename is not None:
        pf.writeto(filename, noise_cube, header=header)

    return noise_cube

def calculate_noise(profile, velocity_axis, velocity_range):
    """ Calculates rms noise of Tile profile given velocity ranges.
    """
    import numpy as np

    std = 0

    # calculate noises for each individual region
    for i in range(len(velocity_range) / 2):
        velMin = velocity_range[2*i + 0]
        velMax = velocity_range[2*i + 1]

        noise_region = np.where((velocity_axis >= velMin) & \
                        (velocity_axis <= velMax))

        std += np.std(profile[noise_region])

    std /= len(velocity_range) / 2
    return std

def calculate_noise_scale(Tsys, profile, noise=None):
    """ Creates an array for scaling the noise by (Tsys + Tb) / Tsys
    """
    import numpy as np
    n = np.zeros(len(profile))
    n = (Tsys + profile) / Tsys * noise

    return n

def calculate_sd(image, sd_factor=1/1.25):

    ''' Calculates a surface density image given a velocity range within which
    to include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to '''

    import numpy as np

    # NHI in units of 1e20 cm^-2
    sd_image = image * sd_factor

    return sd_image

def calculate_nh2(nhi_image = None, av_image = None, dgr = 1.1e-1):

    ''' Calculates the total gas column density given N(HI), A_v and a
    dust-to-gas ratio.

    Parameters
    ----------
    '''

    import numpy as np

    nh_image = 0.5 * (av_image / dgr - nhi_image)

    return nh_image

def calculate_nh2_error(nhi_image_error=None, av_image_error=None, dgr=1.1e-1):

    ''' Calculates the total gas column density given N(HI), A_v and a
    dust-to-gas ratio.

    Parameters
    ----------
    '''

    import numpy as np

    nh2_image_error = 0.5 * ((av_image_error / dgr)**2 - \
            nhi_image_error**2)**0.5

    return nh2_image_error

def correlate_hi_av(hi_cube=None, hi_velocity_axis=None,
        hi_noise_cube=None, av_image=None, av_image_error=None,
        velocity_centers=None, velocity_widths=None, return_correlations=True):

    '''
            hi_vel_range, av_correlations = correlate_hi_av(hi_cube=hi_data,
                    hi_velocity_axis=velocity_axis,
                    hi_noise_cube=noise_cube,
                    av_image=av_data_planck,
                    av_image_error=av_error_data_planck,
                    velocity_centers=velocity_centers,
                    velocity_widths=velocity_widths,
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
    velocity_ranges = np.zeros(shape=[len(velocity_centers) * \
            len(velocity_widths),2])
    count = 0
    for i, center in enumerate(velocity_centers):
        for j, width in enumerate(velocity_widths):
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

        nhi_image = np.ma.array(nhi_image_temp,
                                mask=np.isnan(nhi_image_temp))

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

''' Fitting Functions
'''

def krumholz_eq(h_sd, phi_CNM = None,
        Z = 1., # metallicity
        a = 0.2, # ?
        f_diss = 0.1, # fraction of absorbing H2 which disociates
        phi_mol = 10.0, # molecular gas fraction
        mu_H = 2.3e-24, # molecular weight of H, g
        return_fractions=False
        ):
    '''
    Krumholz et al. 2008
    '''

    # Constants
    c = 3.0e10 # speed of light, cm/s

    # solar values
    sigma_d_solar = 1e-21 # solar dust grain cross section, cm^2
    R_d_solar = 10**-16.5 # solar cloud radius, cm
    E_0_solar = 7.5e-4 # Radiation field, erg/s

    # cloud values
    sigma_d = sigma_d_solar * Z # dust grain cross section, cm^2
    R_d = R_d_solar * Z # cloud radius, cm

    # normalized radiation field strength, EQ 7
    chi = ((f_diss * sigma_d_solar * c * E_0_solar) \
            * (1.0 + (3.1 * Z**0.365))) \
            / (31.0 * phi_CNM * R_d_solar)

    # dust-adjusted radiation field, EQ 10
    psi = chi * (2.5 + chi) / (2.5 + (chi * np.e))

    # cloud optical depth, EQ 21
    tau_c = (3.0 * h_sd * sigma_d) / (4.0 * (3.1 * Z**0.365) * mu_H)

    # cloud optical depth, EQ 21
    tau_c = 0.067 * Z * h_sd

    f_H2_sub1 = (3.0 * psi) / (3.0 * tau_c)
    f_H2_sub2 = (4.0 * a * psi * phi_mol) / ((4.0 * tau_c) + (3.0 * (phi_mol \
            - 1.0) * psi))
    f_H2 = 1.0 - (f_H2_sub1 / (1.0 + f_H2_sub2))
    f_HI = 1.0 - f_H2

    # Keep fractions within real fractional value range
    f_HI[f_HI > 1] = 1.0
    f_HI[f_HI < 0] = 0.0
    f_H2[f_H2 > 1] = 1.0
    f_H2[f_H2 < 0] = 0.0

    # ratio of molecular to atomic fraction, EQ 17 Lee et al. 2012
    R_H2 = 4 * tau_c / (3 * psi) \
            * (1+ 0.8 * psi * phi_mol \
                / (4 * tau_c + 3 * (phi_mol - 1) * psi)) -1

    R_H2 = f_H2 / f_HI

    if not return_fractions:
        return R_H2
    elif return_fractions:
        return R_H2, f_H2, f_HI

def krumholz_eq_simple(h_sd, phi_CNM, return_fractions=False):

   return krumholz_eq(h_sd, phi_CNM = phi_CNM,
        Z=1.0, # metallicity
        a=0.2, # ?
        f_diss=0.1, # fraction of absorbing H2 which disociates
        phi_mol=10.0, # molecular gas fraction
        mu_H=2.3e-24, # molecular weight of H, g
        return_fractions=return_fractions,
        )

def fit_krumholz(h_sd, rh2, h_sd_extent, p0 = 10, return_params = False,
        return_fractions=False, return_chisq=False, rh2_error=None):

    '''
    Parameters
    ----------
    h_sd : array-like
        Hydrogen surface density in units of solar mass per parsec**2
    rh2 : array-like
        Ratio between molecular and atomic hydrogen masses.
    h_sd_extent : tuple
        Lower and upper bound of hydrogen surface densities with which to
        build the output model array.
    p0 : None, scalar, or M-length sequence.
        Initial guess for the parameters. See scipy.optimize.curve_fit.
    return_params : bool
        Return parameters from fit?
    return_fractions : bool
        Return f_H2 and f_HI?
    return_chisq : bool
        Return the chi^2 statistic and p-value?
    rh2_error : bool
        Error in rh2 parameter. Calculates a more accurate chi^2 statistic

    Returns
    -------
    rh2_fit : array-like
        Model ratio between molecular and atomic hydrogen masses.
    h_sd_extended : array-like
        Model hydrogen surface density in units of solar mass per parsec**2.
    rh2_fit_params : array-like, optional
        Model parameter fits.
    f_H2, f_HI : array-like, optional
        f_H2 = mass fraction of molecular hydrogen
        f_HI = mass fraction of atomic hydrogen
    chisq, p_value : float, optional
        Chi squared statistic and p-value.

    '''

    from scipy.optimize import curve_fit
    from scipy import stats

    # fit the krumholz model, choose first element of tuple = parameters
    rh2_fit_params = curve_fit(krumholz_eq_simple, h_sd, rh2,
            maxfev=10000000, p0 = p0)[0]

    # Create large array of h_sd
    h_sd_extended = np.linspace(h_sd_extent[0], h_sd_extent[1], h_sd_extent[2])

    if not return_fractions:
        rh2_fit = krumholz_eq_simple(h_sd_extended, rh2_fit_params)
    elif return_fractions:
        rh2_fit, f_H2, f_HI = krumholz_eq_simple(h_sd_extended,
                                                 rh2_fit_params,
                                                 return_fractions=True)

    if return_chisq:
        dof = len(rh2) - 1
        rh2_fit_grid = krumholz_eq_simple(h_sd, rh2_fit_params)
        chisq = np.sum((rh2 - rh2_fit_grid)**2)
        if rh2_error is not None:
            chisq /= np.sum(rh2_error**2)
        p_value = 1 - stats.chi2.cdf(chisq, 1)
        chisq_reduced = chisq / dof

    print('pvalue = %.2f' % p_value)
    print('dof = %.2f' % dof)
    print('chi2 = %.2f' % chisq)

    output = [rh2_fit, h_sd_extended]

    if return_params:
        output.append(rh2_fit_params)
    if return_fractions:
        output.append(f_H2)
        output.append(f_HI)
    if return_chisq:
        output.append(chisq)
        output.append(p_value)

    return output

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

def load_ds9_region(cores, filename_base = 'california_av_boxes_', header=None):

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
    import myclumpfinder as clump_finder
    reload(clump_finder)
    import mygeometry as myg
    reload(myg)
    from mycoords import make_velocity_axis
    import json

    # parameters used in script
    # -------------------------
    # HI velocity integration range
    # Determine HI integration velocity by CO or correlation with Av?
    hi_co_width = False
    hi_av_correlation = True
    co_width_scale = 5.0 # for determining N(HI) vel range
    hi_vel_range_scale = 2.0 # scale hi velocity range for Av/N(HI) correlation

    dgr = 6.33e-2 # dust to gas ratio [10^-22 mag / 10^20 cm^-2
    rh2_fit_range = [0.001, 1000] # range of fitted values for krumholz model


    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/california/figures/cores/'
    av_dir = '/d/bip3/ezbc/california/data/av/'
    hi_dir = '/d/bip3/ezbc/california/data/hi/'
    co_dir = '/d/bip3/ezbc/california/data/co/'
    core_dir = '/d/bip3/ezbc/california/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/california/data/python_output/ds9_regions/'

    # load Planck Av and GALFA HI images, on same grid
    av_data_planck, av_header = load_fits(av_dir + \
                'california_av_planck_5arcmin.fits',
            return_header=True)

    av_error_data_planck, av_error_header = load_fits(av_dir + \
                'california_av_error_planck_5arcmin.fits',
            return_header=True)

    hi_data, h = load_fits(hi_dir + \
                'california_hi_galfa_cube_regrid_planckres.fits',
            return_header=True)

    co_data, co_header = load_fits(co_dir + \
                'california_co_cfa_cube_regrid_planckres.fits',
            return_header=True)

    # make the velocity axis
    velocity_axis = make_velocity_axis(h)
    co_vel_axis = make_velocity_axis(co_header)

    # Plot NHI vs. Av for a given velocity range
    noise_cube_filename = 'california_hi_galfa_cube_regrid_planckres_noise.fits'
    if not path.isfile(hi_dir + noise_cube_filename):
        noise_cube = calculate_noise_cube(cube=hi_data,
                velocity_axis=velocity_axis,
                velocity_noise_range=[90,110], header=h, Tsys=30.,
                filename=hi_dir + noise_cube_filename)
    else:
        noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    # calculate nhi and error maps, write nhi map to fits file
    nhi_image, nhi_image_error = calculate_nhi(cube=hi_data,
        velocity_axis=velocity_axis,
        noise_cube = noise_cube,
        velocity_range=[0,15],
        return_nhi_error=True,
        fits_filename = hi_dir + 'california_nhi_galfa_5arcmin.fits',
        fits_error_filename = hi_dir + 'california_nhi_galfa_5arcmin_error.fits',
        header = h)

    # calculate N(H2) maps
    nh2_image = calculate_nh2(nhi_image = nhi_image,
            av_image = av_data_planck, dgr = dgr)
    nh2_image_error = calculate_nh2(nhi_image = nhi_image_error,
            av_image = av_error_data_planck, dgr = dgr)

    # convert to column density to surface density
    hi_sd_image = calculate_sd(nhi_image, sd_factor=1/1.25)
    hi_sd_image_error = calculate_sd(nhi_image_error, sd_factor=1/1.25)

    # Paradis et al. (2012) gives DGR for california
    h_sd_image = av_data_planck / (1.25 * dgr) # DGR = 1.1e-22 mag / cm^-2
    h_sd_image_error = av_error_data_planck / (1.25 * dgr)

    # h2 surface density
    h2_sd_image = h_sd_image - hi_sd_image
    h2_sd_image_error = (h_sd_image_error**2 - hi_sd_image_error**2)**0.5

    # Write ratio between H2 and HI
    rh2_image = h2_sd_image / hi_sd_image
    rh2_image_error = rh2_image * (hi_sd_image_error**2 / hi_sd_image**2 \
                 + h2_sd_image_error**2 / h2_sd_image**2)**0.5

    # define core properties
    with open(core_dir + 'california_core_properties.txt', 'r') as f:
        cores = json.load(f)

    cores = convert_core_coordinates(cores, h)

    hsd_limits =[0.1,300]
    hisd_limits = [2,20]
    av_limits =[0.01,100]
    nhi_limits = [2,20]

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'california_av_boxes_',
            header = h)

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
    chisq_list = []
    p_value_list = []
    core_name_list = []
    co_image_list = []
    hi_vel_range_list = []
    hi_sd_fit_list = []

    for core in cores:
        print('\nCalculating for core %s' % core)

        # Grab the mask from the DS9 regions
        '''
        xy = cores[core]['box_center_pix']
        box_width = cores[core]['box_width']
        box_height = cores[core]['box_height']
        box_angle = cores[core]['box_angle']
        mask = myg.get_rectangular_mask(nhi_image,
                xy[0], xy[1],
                width = box_width,
                height = box_height,
                angle = box_angle)

        # Get indices where there is no mask, and extract those pixels
        indices = np.where(mask == 1)
        indices = mask == 1
        '''

        mask = myg.get_polygon_mask(av_data_planck,
                cores[core]['box_vertices_rotated'])

        indices = mask == 1

        if hi_co_width:
            co_data_sub = co_data[:, indices]
            co_image_list.append(np.sum(co_data_sub, axis=1) / \
                    co_data_sub.shape[0])

            # fit gaussians to CO data
            hi_vel_range = select_hi_vel_range(co_data_sub, co_header,
                    flux_threshold=0.85, width_scale=co_width_scale)
            hi_vel_range_list.append(hi_vel_range)

        if hi_av_correlation:
            hi_vel_range = cores[core]['hi_velocity_range']
            correlation_coeff = cores[core]['correlation_coeff']

            hi_vel_range = (hi_vel_range[0] * hi_vel_range_scale,
                    hi_vel_range[1] * hi_vel_range_scale)


        if 1:
            print('HI velocity integration range:')
            print('%.0f to %.0f km/s' % (hi_vel_range[0], hi_vel_range[1]))
            print('With correlation of:')
            print('%.2f' % correlation_coeff)

            nhi_image, nhi_image_error = calculate_nhi(cube=hi_data,
                    velocity_axis=velocity_axis,
                    noise_cube=noise_cube,
                    velocity_range=hi_vel_range,
                    return_nhi_error=True,
                    header=h)

            nhi_image = np.ma.array(nhi_image,
                    mask=(nhi_image != nhi_image))
            nhi_image_error = np.ma.array(nhi_image_error,
                    mask=(nhi_image_error != nhi_image_error))

            # calculate N(H2) maps
            nh2_image = calculate_nh2(nhi_image = nhi_image,
                    av_image = av_data_planck, dgr = dgr)

            nh2_image_error = \
                calculate_nh2_error(nhi_image_error=nhi_image_error,
                    av_image_error=av_error_data_planck, dgr = dgr)

            # convert to column density to surface density
            hi_sd_image = calculate_sd(nhi_image, sd_factor=1/1.25)
            hi_sd_image_error = calculate_sd(nhi_image_error, sd_factor=1/1.25)

            h2_sd_image = calculate_sd(nh2_image, sd_factor=1/0.625)
            h2_sd_image_error = calculate_sd(nh2_image_error,
                    sd_factor=1/0.625)

            # Paradis et al. (2012) gives DGR for california
            #h_sd_image = av_data_planck / (1.25 * dgr)
            #h_sd_image_error = av_error_data_planck / (1.25 * dgr)
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

        # Write velocity range properties to a file
        with open(core_dir + 'california_hi_vel_properties.txt', 'w') as f:
            f.write('core\tco_range\thi_range')


        av_data_planck_sub = av_data_planck[indices]
        av_error_data_planck_sub = av_error_data_planck[indices]
        hi_sd_image_sub = hi_sd_image[indices]
        hi_sd_image_error_sub = hi_sd_image_error[indices]
        h_sd_image_sub = h_sd_image[indices]
        h_sd_image_error_sub = h_sd_image_error[indices]
        rh2_image_sub = rh2_image[indices]
        rh2_image_error_sub = rh2_image_error[indices]

        # Fit R_H2
        #---------
        # Unravel images to single dimension
        rh2_ravel = rh2_image_sub.ravel()
        rh2_error_ravel = rh2_image_error_sub.ravel()
        h_sd_ravel = h_sd_image_sub.ravel()
        h_sd_error_ravel = h_sd_image_error_sub.ravel()

        # write indices for only ratios > 0
        indices = np.where(rh2_ravel > 0)

        # Fit to krumholz model, init guess of phi_CNM = 10
        rh2_fit_range.append(1e6)
        rh2_fit, h_sd_fit, phi_cnm, f_H2, f_HI, chisq, p_value= \
                fit_krumholz(h_sd_ravel[indices],
                             rh2_ravel[indices],
                             rh2_fit_range,
                             p0=10.0,
                             return_params=True,
                             return_fractions=True,
                             return_chisq=True,
                             rh2_error=rh2_image_error_sub)

        print('phi = %.2f' % phi_cnm)

        # see eq 6 of krumholz+09
        # phi_cnm is the number density of the CNM over the minimum number
        # density required for pressure balance
        # the lower phi_cnm values than for california mean that california
        # has a more diffuse CNM

        # By fitting the model to the observed R_H2 vs total H, you basically
        # constrained psi in Equation (35) of Krumholz+09.  This means that
        # you can calculate f_H2 for a given total hydrogen surface density.
        # In this case, H2 surface density = f_H2 * total hydrogen surface
        # density HI surface density = (1 - f_HI) * total hydrogen surface
        # density

        hi_sd_fit = f_HI * h_sd_fit

        # append to the lists
        hi_sd_image_list.append(hi_sd_image_sub)
        hi_sd_image_error_list.append(hi_sd_image_error_sub)
        h_sd_image_list.append(h_sd_image_sub)
        h_sd_image_error_list.append(h_sd_image_error_sub)
        av_image_list.append(av_data_planck_sub)
        av_image_error_list.append(av_error_data_planck_sub)
        rh2_image_list.append(rh2_image_sub)
        rh2_image_error_list.append(rh2_image_error_sub)
        rh2_fit_list.append(rh2_fit)
        h_sd_fit_list.append(h_sd_fit)
        hi_sd_fit_list.append(hi_sd_fit)
        phi_cnm_list.append(phi_cnm)
        chisq_list.append(chisq)
        p_value_list.append(p_value)
        core_name_list.append(core)

    figure_types = ['png', 'pdf']
    for figure_type in figure_types:
        if hi_co_width:
            plot_co_spectrum_grid(co_vel_axis,
                    co_image_list,
                    vel_range_list=hi_vel_range_list,
                    #limits = [0, 80, 10**-3, 10**2],
                    savedir = figure_dir + 'panel_cores/',
                    scale = ('linear', 'linear'),
                    filename = 'california_core_co_spectra_grid.%s' % figure_type,
                    title = r'Average $^{12}$CO spectra of california Cores',
                    core_names=core_name_list,
                    show = False)

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
                filename = 'california_rh2_vs_hsd_panels_planck.%s' % figure_type,
                title = r'$R_{\rm H2}$ vs. $\Sigma_{\rm HI}$'\
                        + ' of california Cores',
                core_names=core_name_list,
                phi_cnm_list=phi_cnm_list,
                chisq_list=chisq_list,
                p_value_list=p_value_list,
                show = False)

        plot_hi_vs_h_grid(hi_sd_image_list,
                h_sd_image_list,
                hi_sd_error_images = hi_sd_image_error_list,
                h_sd_error_images = h_sd_image_error_list,
                hi_fits = hi_sd_fit_list,
                h_sd_fits = h_sd_fit_list,
                #limits = [10**-1, 10**2, 10**0, 10**2],
                #limits = [-5, 50, 1, 8],
                savedir = figure_dir + 'panel_cores/',
                scale = ('linear', 'linear'),
                filename = 'california_hi_vs_h_panels_planck.%s' % figure_type,
                title = r'$\Sigma_{\rm HI}$ vs. $\Sigma_{\rm H}$'\
                        + ' of california Cores',
                core_names=core_name_list,
                phi_cnm_list=phi_cnm_list,
                show = False)

        plot_hi_vs_h_grid(hi_sd_image_list,
                h_sd_image_list,
                hi_sd_error_images = hi_sd_image_error_list,
                h_sd_error_images = h_sd_image_error_list,
                hi_fits = hi_sd_fit_list,
                h_sd_fits = h_sd_fit_list,
                limits = [10**0, 10**2, 10**0, 10**2],
                savedir = figure_dir + 'panel_cores/',
                scale = ('log', 'log'),
                filename = 'california_hi_vs_h_panels_planck_log.%s' %\
                        figure_type,
                title = r'$\Sigma_{\rm HI}$ vs. $\Sigma_{\rm H}$'\
                        + ' of california Cores',
                core_names=core_name_list,
                phi_cnm_list=phi_cnm_list,
                show = False)

if __name__ == '__main__':
    main()



