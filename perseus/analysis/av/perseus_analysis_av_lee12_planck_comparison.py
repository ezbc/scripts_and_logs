#!/usr/bin/python

from astropy.io import fits as pf
import numpy as np
import pyfits as pf

def plot_av_residuals(av_image1=None, av_image2=None, title=None, limits=None,
        savedir='./', filename=None, show=True, header=None):

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
              'figure.figsize': (8, 6),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    print(type(header))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    av_image1 = np.ma.array(av_image1, mask=(av_image1 != av_image1))
    av_image2 = np.ma.array(av_image2, mask=(av_image2 != av_image2))

    # create axes
    ax = imagegrid[0]
    cmap = cm.gist_stern # colormap

    # show the image
    im = ax.imshow(av_image1 - av_image2,
            interpolation='nearest',
            origin='lower',
            #norm=matplotlib.colors.LogNorm(),#vmin=0.01, vmax=10,),
            cmap=cmap,)

    # Asthetics
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',)
    ax.set_ylabel('Declination (J2000)',)
    if title is not None:
    	ax.set_title(title)

    # colorbar
    cb = ax.cax.colorbar(im,) #ticks=[0.01, 0.1, 1, 10])
    #cb.ax.set_xticklabels(['%.2f'%0.01, '0.1', '1', '10'])# horizontal colorbar

    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Write label to colorbar
    cb.set_label_text('Residuals (Mag)',)

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        fig.show()

def plot_av(av_x, av_y,
        av_x_error=None, av_y_error=None, limits=None,
        savedir='./', filename=None, show=True, scale='linear',
        returnimage=False, hess_binsize=None, title='', plot_type='scatter',
        color_scale = 'linear', errorbar = None, errorbar_pos = None):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib import cm

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
              'figure.figsize': (8, 8),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)


    # Create figure
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

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
                    xscale=scale,
                    cmap = cm.gist_stern)
                for i, pos in enumerate(errorbar_pos):
                	error_bar = ax.errorbar(pos[0], pos[1],
                            xerr = (errorbar[0]),
                            yerr = (errorbar[1]),
                            color = 'c')
                ax.set_xscale(scale, nonposx = 'clip')
                ax.set_yscale(scale, nonposy = 'clip')
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
    ax.set_ylabel(r'Lee+12 A$_{\rm V}$ (mag)')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def print_stats(av_image = None, filename = ''):

    #av_image = np.ma.array(av_image, mask=(av_image != av_image))

    av_image = av_image[av_image == av_image]

    mean = np.mean(av_image)
    median = np.median(av_image)
    maximum = np.max(av_image)
    minimum = np.min(av_image)
    std = np.std(av_image)

    f = open(filename, 'w')
    f.write('mean\t median\t max\t min\t std\t\n')
    f.write('%.2f\t %.2f\t %.2f\t %.2f\t %.2f\t' % \
            (mean, median, maximum, minimum, std))
    f.close

    print(filename)
    print('mean\t median\t max\t min\t std\t')
    print('%.2f\t %.2f\t %.2f\t %.2f\t %.2f\t' % \
            (mean, median, maximum, minimum, std))

def main():
    ''' Executes script.

    For a log of which files Min used to establish the correlation see:
    HI is copied from /d/leffe2/lee/perseus_cloud/Min/R_H2/H2_102311/
        FIR_based/mask/HI_cdensity.sub.conv4.3.congrid4.3.sub.mask.tmask.fits
    to
        /d/bip3/ezbc/perseus/data/galfa/perseus_galfa_lee12_masked.fits

    Av is copied from /d/leffe2/lee/perseus_cloud/Min/R_H2/H2_102311/
        FIR_based/T_dust/Av_add.fits
    to
        /d/bip3/ezbc/perseus/data/2mass/perseus_av_lee12_masked.fits

    '''

    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)
    import mygeometry as myg
    reload(myg)

    # define directory locations
    output_dir = '/d/bip3/ezbc/perseus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/perseus/figures/av_comparison/'
    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    hi_dir = '/d/bip3/ezbc/perseus/data/galfa/'
    core_dir = output_dir + 'core_arrays/'
    region_dir = '/d/bip3/ezbc/perseus/data/python_output/ds9_regions/'

    # load 2mass Av and GALFA HI images, on same grid
    av_data_2mass, av_header = pf.getdata(av_dir + \
                'perseus_av_lee12_masked_regrid_planckres.fits',
            header = True)

    av_data_planck, av_header = pf.getdata(av_dir + \
                'perseus_av_planck_5arcmin.fits',
            header=True)

    figure_types = ['pdf', 'png']
    if 1:
        for figure_type in figure_types:
            plot_av(av_data_planck, av_data_2mass,
                    plot_type = 'hexbin',
                    color_scale = 'log',
                    scale = 'log',
                    limits = (0.2, 30, 0.2, 30),
                    errorbar = (0.1, 0.3),
                    errorbar_pos = ((2., 2.),
                                    (0.5, 0.5),
                                    (9., 9.)),
                    title = 'Perseus: 2MASS / Planck Comparison',
                    savedir = figure_dir,
                    filename = 'perseus_av_lee12_planck_compare.%s' % \
                            figure_type,
                    show = False,
                    )

    # --------------------------------------------------------------------------
    # Print the stats of the images
    # --------------------------------------------------------------------------

    print_stats(av_image = av_data_2mass, filename = output_dir +
            'perseus_av_stats_2mass.txt')

    print_stats(av_image = av_data_planck, filename = output_dir +
            'perseus_av_stats_planck.txt')

    # --------------------------------------------------------------------------
    # Create map of the residuals
    # --------------------------------------------------------------------------
    for figure_type in figure_types:
        plot_av_residuals(av_image1 = av_data_2mass,
                    av_image2 = av_data_planck,
                    #limits=[70, 50, 250, 200],
                    header = av_header,
                    title = 'Perseus: 2MASS - Planck Residuals',
                    savedir = figure_dir,
                    filename = 'perseus_av_2mass_planck_residual_map.%s' % \
                            figure_type,
                    show = False,)



if __name__ == '__main__':
    main()



