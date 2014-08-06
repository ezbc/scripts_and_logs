#!/usr/bin/python

''' Calculates the N(H) / Av correlation for the Taurus molecular cloud. Uses
Xco factor from Paradis et al. (2012) A&A, 543, 103 to calculate N(H2).
Integrates GALFA HI image to determine N(HI).

'''


def plot_av_vs_nhi(nhi_image, av_image, limits=None,
        savedir='./', filename=None, show=False, scale='linear',
        returnimage=False, hess_binsize=None, title='', plot_type='hexbin',
        color_scale='linear'):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from matplotlib import cm
    from mpl_toolkits.axes_grid1 import ImageGrid

    # Drop the NaNs from the images
    indices = np.where((nhi_image == nhi_image) &\
                       (av_image == av_image) &\
                       (nhi_image > 0) &\
                       (av_image > 0))

    try:
        nhi_image_nonans = nhi_image[indices]
        av_image_nonans = av_image[indices]

        if type(av_image_error) is float:
            av_image_error_nonans = sd_image_error * \
                    np.ones(av_image[indices].shape)
        else:
            av_image_error_nonans = sd_image_error[indices]

        if type(nhi_image_error) is np.ndarray:
            nhi_image_error_nonans = nhi_image_error[indices]
        else:
            nhi_image_error_nonans = nhi_image_error * \
                    np.ones(nhi_image[indices].shape)
    except NameError:
        no_errors = True

    # Create figure
    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    fig_size = (4,4)
    font_scale = 10
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
              'figure.figsize': fig_size,
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)


    # Create figure
    plt.clf()
    fig = plt.figure()
    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1, 1),
                 ngrids=1,
                 axes_pad=0.25,
                 aspect=False,
                 label_mode='L',
                 share_all=True,
                 cbar_mode='single',
                 cbar_pad=0.1,
                 cbar_size=0.2)

    ax = imagegrid[0]

    if plot_type is 'hexbin':
        if color_scale == 'linear':
            image = ax.hexbin(nhi_image_nonans.ravel(),
                    av_image_nonans.ravel(),
                    mincnt=1,
                    yscale=scale,
                    xscale=scale,
                    cmap = cm.gist_stern)
            ax.set_xscale(scale, nonposx = 'clip')
            ax.set_yscale(scale, nonposy = 'clip')
            cb = ax.cax.colorbar(image,)
            # Write label to colorbar
            cb.set_label_text('Bin Counts',)
        elif color_scale == 'log':
            image = ax.hexbin(nhi_image_nonans.ravel(),
                av_image_nonans.ravel(),
                norm=matplotlib.colors.LogNorm(),
                mincnt=1,
                yscale=scale,
                xscale=scale,
                gridsize=(100,200),
                cmap = cm.gist_stern)
            ax.set_xscale(scale, nonposx = 'clip')
            ax.set_yscale(scale, nonposy = 'clip')
            cb = ax.cax.colorbar(image,)
            # Write label to colorbar
            cb.set_label_text('Bin Counts',)
        # Adjust color bar of density plot
        #cb = image.colorbar(image)
        #cb.set_label('Bin Counts')
        elif plot_type is 'scatter':
            image = ax.scatter(nhi_image_nonans.ravel(),
                    av_image_nonans.ravel(),
                    alpha=0.3,
                    color='k'
                    )
            ax.set_xscale(scale)
            ax.set_yscale(scale)

    if limits is not None:
    	ax.set_xlim(limits[0],limits[1])
    	ax.set_ylim(limits[2],limits[3])

    # Adjust asthetics
    ax.set_xlabel(r'$N(HI)$ (10$^{20}$ cm$^{-2}$)')
    ax.set_ylabel(r'A$_{\rm V}$ (mag)')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight', dpi=600)
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def main():

    import grid
    import numpy as np
    from myimage_analysis import calculate_nhi
    from mycoords import make_velocity_axis
    import pyfits as pf

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/dgr/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/hi/'
    co_dir = '/d/bip3/ezbc/taurus/data/cfa/'

    av_data_planck, planck_header = pf.getdata(av_dir + \
                'taurus_av_planck_5arcmin.fits',
            header=True)

    # load GALFA HI
    hi_data, hi_header = pf.getdata(hi_dir + \
            'taurus_hi_galfa_cube_regrid_planckres.fits',
            header=True)
    velocity_axis = make_velocity_axis(hi_header)

    noise_cube, noise_header = pf.getdata(hi_dir + \
            'taurus_hi_galfa_cube_regrid_planckres_noise.fits', header=True)

    # Derive N(HI) image
    nhi_image, nhi_image_error = calculate_nhi(cube=hi_data,
            velocity_axis=velocity_axis,
            noise_cube=noise_cube,
            velocity_range=[0,15])

    # Plot correlation, similar to Figure 3 of Paradis et al. (2012)
    plot_av_vs_nhi(nhi_image,
            av_data_planck,
            #limits=[0,14, 0,10],
            savedir=figure_dir,
            filename='taurus_av_vs_nhi.png',
            color_scale='linear')

if __name__ == '__main__':
    main()


