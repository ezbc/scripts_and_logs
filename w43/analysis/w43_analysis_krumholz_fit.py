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
        return likelihoods_image

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
        return likelihoods_image

def plot_rh2_vs_h_grid(rh2_images, h_sd_images, rh2_error_images=None,
        h_sd_error_images=None, rh2_fits = None, h_sd_fits = None, limits =
        None, fit = True, savedir = './', filename = None, show = True, scale =
        'linear', title = '', core_names='', phi_cnm_list=None,
        phi_cnm_error_list=None, Z_list=None,
        Z_error_list=None, phi_mol_list=None,
        phi_mol_error_list=None,
        chisq_list=None, p_value_list=None,):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from myscience.krumholz09 import calc_T_cnm

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
        if phi_mol_list is not None:
            phi_mol = phi_mol_list[i]
        if phi_mol_error_list is not None:
            phi_mol_error = phi_mol_error_list[i]
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
                alpha=0.75,
                color='k',
                marker='^',ecolor='k',linestyle='none',
                markersize=4
                )

        if rh2_fit is not None:
            ax.plot(h_sd_fit, rh2_fit,
                    color = 'r', alpha=0.5)

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
                                  alpha=0.3),
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        )
            else:
            	T_cnm = calc_T_cnm(phi_cnm, Z=Z)
            	T_cnm_error = []
            	T_cnm_error.append(\
            	        T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[0], Z=Z))
            	T_cnm_error.append(\
            	        T_cnm - calc_T_cnm(phi_cnm + phi_cnm_error[1], Z=Z))

                phi_cnm_text = r'\noindent$\phi_{\rm CNM}$ =' + \
                               r' %.2f' % (phi_cnm) + \
                               r'$^{+%.2f}_{-%.2f}$ \\' % (phi_cnm_error[0],
                                                         phi_cnm_error[1])
                T_cnm_text = r'\noindent T$_{\rm CNM}$ =' + \
                             r' %.2f' % (T_cnm) + \
                             r'$^{+%.2f}_{-%.2f}$ \\' % (T_cnm_error[0],
                                                         T_cnm_error[1])
                if Z_error == (0.0, 0.0):
                	Z_text = r'Z = %.1f Z$_\odot$ \\' % (Z)
                else:
                	Z_text = r'Z = %.2f' % (Z) + \
                    r'$^{+%.2f}_{-%.2f}$ Z$_\odot$ \\' % (Z_error[0],
                                                          Z_error[1])
                if phi_mol_error == (0.0, 0.0):
                    phi_mol_text = r'\noindent$\phi_{\rm mol}$ = ' + \
                                     '%.1f \\' % (phi_mol)
                else:
                    phi_mol_text = r'\noindent$\phi_{\rm mol}$ =' + \
                                r' %.2f' % (phi_mol) + \
                                r'$^{+%.2f}_{-%.2f}$ \\' % (phi_mol_error[0],
                                                         phi_mol_error[1]) + \
                                r''

                ax.annotate(phi_cnm_text + T_cnm_text + Z_text + phi_mol_text,
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
        ax.grid(False)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale*1.5)
    if filename is not None:
        plt.savefig(savedir + filename) #, bbox_inches='tight')
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
                    color = 'r', alpha=0.2, label=r'$^{12}$CO Spectrum Width')
        if hi_velocity_range_likelihood is not None:
            ax.axvspan(hi_velocity_range_likelihood[0], hi_velocity_range_likelihood[1],
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

def plot_parameter_hist(results_dict, parameter='Z fits',
        results_figure_name=None):

    import matplotlib.pyplot as plt
    import numpy as np

    if parameter == 'Z fits':
        bins = np.linspace(0, 10, 100)
    elif parameter == 'phi_cnm fits':
        bins = 100

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    ax.hist(results_dict[parameter],
            bins=bins)
    ax.set_ylabel('Counts')
    if parameter == 'Z fits':
        ax.set_xlabel(r'$Z$ $(Z_\odot)$')
        plt.savefig(results_figure_name + '_Z_hist.png')
    elif parameter == 'phi_cnm fits':
        ax.set_xlabel(r'$\phi_{\rm CNM}$')
        plt.savefig(results_figure_name + '_phi_cnm_hist.png')

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
    ax.set_xlabel(r'Velocity (km/s)')
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
    ax.set_xlabel(r'DGR')
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
        parameter_vary=[True,True, True], core_dict=None,
        results_figure_name=None, error_method='bootstrap', alpha=0.05,
        properties=None, results_filename=None, clobber=True,
        likelihood_derivation='global'):

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
    from mystats import rv2d_discrete
    import matplotlib.pyplot as plt
    import json
    from os import path

    if N_runs < 1:
        raise ValueError('N_runs must be >= 1')

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
    results_dict = {'phi_cnm fits' : np.empty((N_runs)),
                    'Z fits' : np.empty((N_runs)),
                    'phi_mol fits' : np.empty((N_runs)),
                    'nhi errors' : np.empty((N_runs))}
    hi_vel_range_list = np.empty((N_runs, 2))
    dgr_list = np.empty((N_runs))

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

    center_likelihoods = np.asarray(likelihood_param_dict['center_likelihood'])
    width_likelihoods = np.asarray(likelihood_param_dict['width_likelihood'])
    dgr_likelihoods = np.asarray(likelihood_param_dict['dgr_likelihood'])
    likelihoods = np.asarray(likelihood_param_dict['likelihoods'])
    vel_centers = np.asarray(likelihood_param_dict['vel_centers'])
    vel_widths = np.asarray(likelihood_param_dict['vel_widths'])
    dgrs = np.asarray(likelihood_param_dict['dgrs'])

    # Derive PDF of the velocity centers and widths.
    # The Monte Carlo will draw from this PDF randomly.
    # rv_discrete requires the PDF be normalized to have area of 1
    center_likelihood_normed = center_likelihoods / np.sum(center_likelihoods)
    width_likelihood_normed = width_likelihoods / np.sum(width_likelihoods)
    dgr_likelihood_normed = dgr_likelihoods / np.sum(dgr_likelihoods)
    likelihoods_normed = likelihoods / np.sum(likelihoods)

    # Create 2-D PDF of velocity widths + DGRs
    likelihood_pdf = rv2d_discrete(likelihoods=likelihoods_normed,
                           param_grid1=vel_widths,
                           param_grid2=dgrs,
                           param_name1='widths',
                           param_name2='dgrs',
                           L_scalar=N_runs)

    if perform_mc:
        # Run the Monte Carlo
        # -------------------

        # Progress bar parameters
        total = float(N_runs)
        count = 0

        for i in xrange(N_runs):
            if N_runs > 1:
                # Randomly sample from Gaussian distribution
                hi_random_error = normal(scale=hi_error,\
                        size=hi_noise_cube.shape)
                av_random_error = normal(scale=av_error, size=av_image.shape)

                # Get a random sample from the DGR / width likelihood space
                vel_width_random, dgr_random = likelihood_pdf.rvs()

                # Create the velocity range
                hi_vel_range_noise = (vel_centers[0] - \
                                        vel_width_random / 2.,
                                      vel_centers[0] + \
                                        vel_width_random / 2.)

                hi_vel_range_list[i, 0] = hi_vel_range_noise[0]
                hi_vel_range_list[i, 1] = hi_vel_range_noise[1]
                dgr_list[i] = dgr_random

                # Add random error to images
                hi_cube_noise = np.copy(hi_cube) + hi_random_error
                av_image_noise = np.copy(av_image) + av_random_error
            elif N_runs == 1:
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

            # Fit to krumholz model, init guess of phi_CNM = 10
            phi_cnm, Z, phi_mol = fit_krumholz(h_sd_ravel,
                                      rh2_ravel,
                                      guesses=guesses, # phi_cnm, Z
                                      vary=parameter_vary,
                                      rh2_error=rh2_error_ravel,
                                      verbose=verbose)

            # keep results
            results_dict['phi_cnm fits'][i] = phi_cnm
            results_dict['phi_mol fits'][i] = phi_mol
            results_dict['Z fits'][i] = Z

            if verbose:
                print('phi = %.2f' % phi_cnm)
                print('phi = %.2f' % phi_mol)
                print('Z = %.2f' % Z)

            # see eq 6 of krumholz+09
            # phi_cnm is the number density of the CNM over the minimum number
            # density required for pressure balance
            # the lower phi_cnm values than for taurus mean that taurus
            # has a more diffuse CNM

            # By fitting the model to the observed R_H2 vs total H, you
            # basically constrained psi in Equation (35) of Krumholz+09.  This
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
            phi_cnms = results_dict['phi_cnm fits'].tolist()
            phi_mols = results_dict['phi_mol fits'].tolist()
            Zs = results_dict['Z fits'].tolist()
            vel_ranges = hi_vel_range_list.tolist()
            dgr_list = dgr_list.tolist()

            results = {'phi_cnms': phi_cnms,
                       'phi_mols': phi_mols,
                       'Zs': Zs,
                       'dgrs': dgr_list,
                       'vel_ranges': hi_vel_range_list.tolist()}

            with open(results_filename, 'w') as f:
                json.dump(results, f)

    # If not performing the MC, read a previously written MC file
    elif not perform_mc:
        print('Reading monte carlo results file:')
        print(results_filename)

        with open(results_filename, 'r') as f:
            results = json.load(f)

        results_dict['phi_cnm fits'] = np.asarray(results['phi_cnms'])
        results_dict['phi_mol fits'] = np.asarray(results['phi_mols'])
        results_dict['Z fits'] = np.asarray(results['Zs'])
        dgr_list = np.asarray(results['dgrs'])
        hi_vel_range_list = np.asarray(results['vel_ranges'])

    # Remove failed fits
    results_dict['phi_cnm fits'] = \
        results_dict['phi_cnm fits']\
            [~np.isnan(results_dict['phi_cnm fits'])]
    results_dict['phi_mol fits'] = \
        results_dict['phi_mol fits']\
            [~np.isnan(results_dict['phi_mol fits'])]
    results_dict['Z fits'] = \
        results_dict['Z fits'] \
            [~np.isnan(results_dict['Z fits'])]

    # Plot the distributions of parameters from the MC
    if results_figure_name is not None:
    	plot_parameter_hist(results_dict,
    	                    parameter='phi_cnm fits',
    	                    results_figure_name=results_figure_name)
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
                           av_image=av_image,
                           av_image_error=av_image_error,
                           )

    phi_cnm_confint, Z_confint, phi_mol_confint = calc_MC_errors(results_dict,
                                                error_method=error_method,
                                                alpha=alpha,
                                                parameter_vary=parameter_vary)

    phi_cnm = phi_cnm_confint[0]
    phi_cnm_error = phi_cnm_confint[1:]
    phi_mol = phi_mol_confint[0]
    phi_mol_error = phi_mol_confint[1:]
    Z = Z_confint[0]
    Z_error = Z_confint[1:]

    # Print results
    print('results are:')
    print('phi_cnm = {0:.2f} +{1:.2f}/-{2:.2f}'.format(phi_cnm_confint[0],
                                                      phi_cnm_confint[1],
                                                      phi_cnm_confint[2]))
    print('phi_mol = {0:.2f} +{1:.2f}/-{2:.2f}'.format(phi_mol_confint[0],
                                                       phi_mol_confint[1],
                                                       phi_mol_confint[2]))
    print('Z = {0:.2f} +{1:.2f}/-{2:.2f}'.format(Z_confint[0],
                                                  Z_confint[1],
                                                  Z_confint[2]))

    print('Median HI velocity range:' + \
          '{0:.1f} to {1:.1f} km/s'.format(hi_vel_range_sample[0],
                                           hi_vel_range_sample[1],))

    if N_runs > 1:
        return images, hi_vel_range_sample, (phi_cnm, Z, phi_cnm_error,
                Z_error, phi_mol, phi_mol_error)
    if N_runs == 1:
        return images, hi_vel_range_sample, (phi_cnm, Z, phi_mol)

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
        # samples = mystats.bootstrap(results_dict['phi_cnm fits'], 1000)
        # phi_cnm_confint = mystats.calc_bootstrap_error(samples, alpha)
        # samples = mystats.bootstrap(results_dict['Z fits'], 1000)
        # Z_confint = mystats.calc_bootstrap_error(samples, alpha)
        if parameter_vary[0]:
            phi_cnm_confint = ci(results_dict['phi_cnm fits'],
                                 statfunction=np.median,
                                 alpha=alpha,
                                 method='bca')
            phi_cnm = np.median(results_dict['phi_cnm fits'])
            phi_cnm_confint = (phi_cnm,
                               phi_cnm + phi_cnm_confint[1],
                               phi_cnm - phi_cnm_confint[0],
                               )
        else:
            phi_cnm_confint = (results_dict['phi_cnm fits'][0], 0.0, 0.0)

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
            phi_cnm = results_dict['phi_cnm fits']

            # Histogram will act as distribution of parameter values
            counts, bins = np.histogram(results_dict['phi_cnm fits'],
                bins=np.linspace(phi_cnm.min(), phi_cnm.max(), 500))

            # Calculate confidence interval by taking the PDF weighted mean as
            # value, where errors are calculated by moving vertical threshold
            # from edges towards mean until alpha/2 *100 % of the area of the
            # PDF is included within the threshold
            #phi_cnm_confint = calc_symmetric_error(bins[:-1], counts,
            #        alpha=alpha)
            phi_cnm_confint = calc_symmetric_error(phi_cnm,
                                                   alpha=alpha)

        else:
            phi_cnm_confint = (results_dict['phi_cnm fits'][0], 0.0, 0.0)

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
            phi_mol_upper_lim = np.log10(np.max(results_dict['phi_mol fits']))
            phi_mol_lower_lim = np.log10(np.min(results_dict['phi_mol fits']))
            phi_mol = results_dict['phi_mol fits']
            phi_mol_confint = calc_symmetric_error(phi_mol,
                    alpha=alpha)
        else:
            phi_mol_confint = (results_dict['phi_mol fits'][0], 0.0, 0.0)

    else:
        raise ValueError('Error method must be "bootstrap" or "threshold"')

    return phi_cnm_confint, Z_confint, phi_mol_confint

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

        rh2_model = k09.calc_rh2(h_sd, phi_cnm, Z, phi_mol=phi_mol, G_0=100)

        chisq = (rh2 - rh2_model)**2 / rh2_error**2

        return chisq

    # Set parameter limits and initial guesses
    params = Parameters()
    params.add('phi_cnm',
               value=guesses[0],
               min=0.1,
               max=100,
               vary=vary[0])
    params.add('phi_mol',
               value=guesses[2],
               min=0.1,
               max=100,
               vary=vary[2])
    params.add('Z',
               value=guesses[1],
               min=0.01,
               max=100,
               vary=vary[1])

    # Perform the fit!
    result = minimize(chisq,
                      params,
                      args=(h_sd, rh2, rh2_error),
                      method='lbfgsb')

    rh2_fit_params = (params['phi_cnm'].value, params['Z'].value,
            params['phi_mol'].value)

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
    hi/taurus_analysis_hi_av_core_likelihoods.py   Calculates HI velocity
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
    from myscience.krumholz09 import calc_T_cnm

    # parameters used in script
    # -------------------------
    # HI velocity integration range
    # Determine HI integration velocity by CO or correlation with Av?
    hi_co_width = True
    hi_av_likelihoodelation = True
    reproduce_lee12 = True

    # Error analysis
    calc_errors = True # Run monte carlo error analysis?
    N_monte_carlo_runs = 1000 # Number of monte carlo runs
    vary_phi_cnm = True # Vary phi_cnm in K+09 fit?
    vary_Z = False # Vary metallicity in K+09 fit?
    vary_phi_mol = False # Vary phi_mol in K+09 fit?
    # Error method:
    # options are 'edges', 'bootstrap'
    error_method = 'edges'
    alpha = 0.32 # 1 - alpha = confidence
    results_filename = '/d/bip3/ezbc/taurus/data/python_output/' + \
            'monte_carlo_results/taurus_mc_results_'
    clobber = 0 # perform MC and write over current results?

    # Use core-derived or global-derived likelihoods for DGR - vel width
    # combinations. Options are 'cores' and 'global'
    likelihood_derivation = 'global'

    # Regions
    # Options are 'ds9' or 'av_gradient'
    box_method = 'av_gradient'

    #dgr = 5.33e-2 # dust to gas ratio [10^-22 mag / 10^20 cm^-2
    h_sd_fit_range = [0.001, 1000] # range of fitted values for krumholz model

    # Figures
    write_pdf_figures = False

    rh2 = np.array((0.1, 0.4, 1, 1.5, 3, 4, 10,))
    hsd = np.array((25, 50, 60, 100, 200, 275, 800,))

    from scipy.optimize import curve_fit
    from myscience import krumholz09 as k09

    #params = fit_krumholz(hsd, rh2, guesses=[10.0, 1.0, 10.0],
    #    verbose=False, vary=[vary_phi_cnm, vary_Z, vary_phi_mol])

    Z = 1.4 # solar masses
    G_0 = 100.0 # solar values

    #Z = 1.0 # solar masses
    #G_0 = 1.0 # solar values

    def krumholz_fitting(hsd, phi_cnm=10):
        rh2 = k09.calc_rh2(hsd, phi_CNM=phi_cnm,
                    Z = Z, # metallicity
                    a = 0.2, # ?
                    f_diss = 0.1, # fraction of absorbing H2 which disociates
                    phi_mol = 10.0, # molecular gas fraction
                    mu_H = 2.3e-24, # molecular weight of H, g
                    return_fractions=False)
        return rh2

    coeff, var_matrix = curve_fit(krumholz_fitting,
                                  hsd,
                                  rh2,
                                  p0=(10.0))

    phi_cnm = coeff[0]
    #phi_mol = coeff[1]
    #phi_cnm = 3

    T_cnm = k09.calc_T_cnm(phi_cnm, Z)
    n_min = k09.calc_n_min(G_0=G_0, Z=Z)
    n_cnm = k09.calc_n_cnm(G_0=G_0, T_cnm=T_cnm, Z=Z)

    print '\n'
    print 'G_0 =', G_0
    print '--------'
    print 'phi_cnm =', phi_cnm,
    #print 'phi_mol =', phi_mol
    print 'T_cnm =', T_cnm, 'K'
    print 'n_min =', n_min, 'cm^-3'
    print 'n_cnm =', n_cnm, 'cm^-3'
    print 'P/k =', n_cnm * T_cnm, 'K cm^-3'

    G_0 = 20
    n_min = k09.calc_n_min(G_0=G_0, Z=Z)
    n_cnm = k09.calc_n_cnm(G_0=G_0, T_cnm=T_cnm, Z=Z)

    print '\n'
    print 'G_0 =', G_0
    print '--------'
    print 'phi_cnm =', phi_cnm,
    #print 'phi_mol =', phi_mol
    print 'T_cnm =', T_cnm, 'K'
    print 'n_min =', n_min, 'cm^-3'
    print 'n_cnm =', n_cnm, 'cm^-3'
    print 'P/k =', n_cnm * T_cnm, 'K cm^-3'


    import matplotlib.pyplot as plt

    h_sd_fit = np.linspace(0.001, 1000, 400)
    rh2_fit = krumholz_fitting(h_sd_fit, phi_cnm=coeff)

    plt.close()
    plt.plot(hsd, rh2, linestyle='', marker='^', label='Motte et al. (2014)')
    plt.plot(h_sd_fit, rh2_fit, label='Krumholz+09 fit')
    plt.xlabel(r'$\Sigma_{H}$')
    plt.ylabel(r'$R_{\rm H2}$')
    plt.yscale('log')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()



