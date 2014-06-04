#!/usr/bin/python

from astropy.io import fits as pf
import numpy as np

def plot_av(av_x, av_y,
        av_x_error=None, av_y_error=None, limits=None,
        savedir='./', filename=None, show=True, scale='linear',
        returnimage=False, hess_binsize=None, title='', plot_type='scatter',
        color_scale = 'linear'):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    # Drop the NaNs from the images
    indices = np.where((av_x == av_x) &\
                       (av_y == av_y) &\
                       (av_x > 0) &\
                       (av_y > 0))

    try:
        av_x_nonans = av_x[indices]
        av_y_nonans = av_y[indices]

        if type(av_y_error) is float:
            av_y_error_nonans = sd_image_error * \
                    np.ones(av_y[indices].shape)
        else:
            av_y_error_nonans = sd_image_error[indices]

        if type(av_x_error) is np.ndarray:
            av_x_error_nonans = av_x_error[indices]
        else:
            av_x_error_nonans = av_x_error * \
                    np.ones(av_x[indices].shape)
    except NameError:
        no_errors = True

    # Create figure
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
    ax = fig.add_subplot(111, aspect = 'equal')

    if av_y_error is None:
        if plot_type is 'hexbin':
            if color_scale == 'linear':
                image = ax.hexbin(av_x_nonans.ravel(),
                        av_y_nonans.ravel(),
                        mincnt=1,
                        yscale=scale,
                        xscale=scale)
            elif color_scale == 'log':
                image = ax.hexbin(av_x_nonans.ravel(),
                    av_y_nonans.ravel(),
                    norm=matplotlib.colors.LogNorm(),
                    mincnt=1,
                    yscale=scale,
                    xscale=scale)
            # Adjust color bar of density plot
            cb = plt.colorbar(image)
            cb.set_label('Bin Counts')
        elif plot_type is 'scatter':
            image = ax.scatter(av_x_nonans.ravel(),
                    av_y_nonans.ravel(),
                    alpha=0.3,
                    color='k'
                    )
            ax.set_xscale(scale)
            ax.set_yscale(scale)
    else:
        image = ax.errorbar(av_x_nonans.ravel(),
                av_y_nonans.ravel(),
                xerr=(av_x_error_nonans.ravel()),
                yerr=(av_y_error_nonans.ravel()),
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
    ax.set_xlabel(r'Planck A$_{\rm V}$ (mag)')
    ax.set_ylabel(r'K+09 A$_{\rm V}$ (mag)')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def main():
    ''' Executes script.

    For a log of which files Min used to establish the correlation see:
    HI is copied from /d/leffe2/lee/taurus_cloud/Min/R_H2/H2_102311/
        FIR_based/mask/HI_cdensity.sub.conv4.3.congrid4.3.sub.mask.tmask.fits
    to
        /d/bip3/ezbc/taurus/data/galfa/taurus_galfa_lee12_masked.fits

    Av is copied from /d/leffe2/lee/taurus_cloud/Min/R_H2/H2_102311/
        FIR_based/T_dust/Av_add.fits
    to
        /d/bip3/ezbc/taurus/data/2mass/taurus_av_lee12_masked.fits

    '''

    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)
    import mygeometry as myg
    reload(myg)

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    core_dir = output_dir + 'core_arrays/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'

    # load 2mass Av and GALFA HI images, on same grid
    av_data_2mass, av_header = pf.getdata(av_dir + \
                'taurus_av_k09_regrid.fits',
            header = True)

    av_data_planck, av_header = pf.getdata(av_dir + \
                '../av/taurus_planck_av_regrid.fits',
            header=True)

    figure_types = ['pdf', 'png']
    for figure_type in figure_types:
        plot_av(av_data_planck, av_data_2mass,
                plot_type = 'hexbin',
                color_scale = 'log',
                scale = 'log',
                limits = (10**-2, 30, 10**-2, 30),
                title = 'Taurus: K+09 / Planck Comparison',
                savedir = figure_dir,
                filename = 'taurus_av_k09_planck_compare.%s' % figure_type,
                show = False,
                )

if __name__ == '__main__':
    main()



