#!/usr/bin/python

''' Calculates the N(HI) map for multicloud

'''

import numpy as np
import matplotlib
matplotlib.use('Agg')

from astropy.io import fits
import numpy as np
import warnings
warnings.filterwarnings('ignore')


''' Plotting Functions
'''

def plot_h2_sd(plot_dict, contour_image=None, limits=None,
        filename=None, show=True, vlimits=None, av_vlimits=None,
        cloud_names=('california', 'perseus', 'taurus'),):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    import pyfits as fits
    import matplotlib.pyplot as plt
    import myplotting as myplt
    import pywcsgrid2
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Set up plot aesthetics
    plt.clf(); plt.close()

    cmap = plt.cm.copper
    norm = matplotlib.colors.LogNorm()


    # Create figure instance
    fig = plt.figure(figsize=(3.6, 9))

    if 1:
        nrows = 3
        ncols = 1
        nrows_ncols=(3, 1)
        ngrids=3
        figsize = (3.6, 8)
        fig = pywcsgrid2.plt.figure(figsize=figsize)
    else:
        nrows_ncols=(1,1)
        ngrids=1

        axes = AxesGrid(fig, (1,1,1),
                     nrows_ncols=nrows_ncols,
                     ngrids=ngrids,
                     cbar_mode="each",
                     cbar_location='right',
                     cbar_pad="2%",
                     cbar_size='3%',
                     axes_pad=0.1,
                     axes_class=(wcs.Axes,
                                 dict(header=plot_dict['header'])),
                     aspect=True,
                     label_mode='L',
                     share_all=True)

    colorbar_axes = [0.05, 0.97, 0.95, 0.02,]
    map_axes = np.array([0.05, 0.69, 0.95, 0.25])

    for i, cloud_name in enumerate(cloud_names):
        header = plot_dict[cloud_name]['header']
        limits = plot_dict[cloud_name]['limits']
        contour_image = plot_dict[cloud_name]['contour_image']
        contours = plot_dict[cloud_name]['contours']

        #ax = pywcsgrid2.subplot(311+i, header=header)
        ax = pywcsgrid2.axes(map_axes, header=header)

        #ax = fig.add_subplot(311 + i, header=)

        h2_sd = plot_dict[cloud_name]['h2_sd']
        # create axes

        # show the image
        im = ax.imshow(h2_sd,
                interpolation='nearest',
                origin='lower',
                cmap=cmap,
                vmin=vlimits[0],
                vmax=vlimits[1],
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        if i == 2:
            ax.set_xlabel('Right Ascension [J2000]',)
        else:
            ax.set_xlabel('')
        ax.set_xlabel('Right Ascension [J2000]',)

        ax.set_ylabel('Declination [J2000]',)

        ax.locator_params(nbins=6)

        #cb_axes.colorbar(im)
        #cb_axes.axis["right"].toggle(ticklabels=False)

        # colorbar
        #cb = ax.cax.colorbar(im)
        cmap.set_bad(color='w')

        # plot limits
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits, header, frame='fk5')
            ax.set_xlim(limits_pix[0],limits_pix[1])
            ax.set_ylim(limits_pix[2],limits_pix[3])

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='w')

        # Write label to colorbar
        #cb.set_label_text(r'$\Sigma_{\rm H_2} [M_\odot$\,pc$^{-2}$]',)

        if 0:
            if regions is not None:
                for region in props['region_name_pos']:
                    vertices = np.copy(regions[region]['poly_verts']['pixel'])
                    rect = ax.add_patch(Polygon(
                            vertices[:, ::-1],
                            facecolor='none',
                            edgecolor='w'))

        ax.annotate(cloud_name.capitalize(),
                    xytext=(0.96, 0.94),
                    xy=(0.96, 0.94),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=10,
                    color='k',
                    bbox=dict(boxstyle='square',
                              facecolor='w',
                              alpha=1),
                    horizontalalignment='right',
                    verticalalignment='top',
                    )

        # create new box for plot
        map_axes[1] -= 0.3
        #map_axes[3] -= 0.3

    ax = plt.axes(colorbar_axes,
                  )
    cb = plt.colorbar(im,
                      cax=ax,
                      orientation='horizontal',
                      )
    cb.ax.set_xlabel(r'$\Sigma_{\rm H_2} [{\rm M_\odot}$\,pc$^{-2}$]',)
    cb.ax.xaxis.set_label_position('top')

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        fig.show()

def plot_hi_h2_sd(plot_dict, contour_image=None, limits=None, filename=None,
        show=True, hi_vlimits=None, h2_vlimits=None, av_vlimits=None,
        cloud_names=('california', 'perseus', 'taurus'),):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    import pyfits as fits
    import matplotlib.pyplot as plt
    import myplotting as myplt
    import pywcsgrid2
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    from mpl_toolkits.axes_grid1.axes_grid import AxesGrid

    # Set up plot aesthetics
    plt.clf(); plt.close()

    # Colormap settings
    cmap = plt.cm.copper
    cmap.set_bad(color='w')
    norm = matplotlib.colors.LogNorm()

    # Create figure instance
    figsize = (7.5, 9)
    fig = pywcsgrid2.plt.figure(figsize=figsize)

    colorbar_axes = [0.05, 0.97, 0.95, 0.02,]
    map_axes = np.array([0.05, 0.69, 0.95, 0.25])
    map_axes = np.array([0.05, 0.31, 0.95, 0.25])

    # plot HI and H2 maps for each cloud
    # ---------------------------------------------------------------------------
    for i, cloud_name in enumerate(cloud_names):
        # get image properties
        # -----------------------------------------------------------------------
        header = plot_dict[cloud_name]['header']
        limits = plot_dict[cloud_name]['limits']
        contour_image = plot_dict[cloud_name]['contour_image']
        contours = plot_dict[cloud_name]['contours']

        # Prep axes instance
        # -----------------------------------------------------------------------
        nrows_ncols=(1, 2)
        ngrids=2

        # only show cbar for first row
        if i == 0:
            cbar_mode = "each"
        else:
            cbar_mode = None

        axes = AxesGrid(fig, 311 + i, #(3,1,i),
                     nrows_ncols=nrows_ncols,
                     ngrids=ngrids,
                     cbar_mode=cbar_mode,
                     cbar_location='top',
                     cbar_pad="5%",
                     cbar_size='10%',
                     axes_pad=0.1,
                     axes_class=(pywcsgrid2.Axes,
                                 dict(header=header)),
                     aspect=True,
                     label_mode='L',
                     share_all=True,
                     )

        # get axes from row
        ax1 = axes[0]
        ax2 = axes[1]

        # Plot the maps
        # -----------------------------------------------------------------------
        # HI surface density map
        hi_sd = plot_dict[cloud_name]['hi_sd']

        # H2 surface density map
        h2_sd = plot_dict[cloud_name]['h2_sd']

        # show the images
        im_hi = ax1.imshow(hi_sd,
                interpolation='nearest',
                origin='lower',
                cmap=cmap,
                vmin=hi_vlimits[0],
                vmax=hi_vlimits[1],
                #norm=matplotlib.colors.LogNorm()
                )
        im_h2 = ax2.imshow(h2_sd,
                interpolation='nearest',
                origin='lower',
                cmap=cmap,
                vmin=h2_vlimits[0],
                vmax=h2_vlimits[1],
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        # -----------------------------------------------------------------------
        for ax in (ax1, ax2):
            ax.set_display_coord_system("fk5")
            ax.set_ticklabel_type("hms", "dms")
            ax.set_xlabel('Right Ascension [J2000]',)
            ax.set_ylabel('Declination [J2000]',)
            ax.locator_params(nbins=5)

        # colorbar
        # -----------------------------------------------------------------------
        cb_hi = ax1.cax.colorbar(im_hi)
        cb_h2 = ax2.cax.colorbar(im_h2)

        # Write label to colorbar
        cb_hi.set_label_text(r'$\Sigma_{\rm HI}$ [M$_\odot$\,pc$^{-2}$]',)
        cb_h2.set_label_text(r'$\Sigma_{\rm H_2}$ [M$_\odot$\,pc$^{-2}$]',)

        # plot limits
        # -----------------------------------------------------------------------
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits, header, frame='fk5')
            ax1.set_xlim(limits_pix[0],limits_pix[1])
            ax1.set_ylim(limits_pix[2],limits_pix[3])
            ax2.set_xlim(limits_pix[0],limits_pix[1])
            ax2.set_ylim(limits_pix[2],limits_pix[3])

        # Plot Av contours
        # -----------------------------------------------------------------------
        if contour_image is not None:
            ax1.contour(contour_image, levels=contours, colors='w')
            ax2.contour(contour_image, levels=contours, colors='w')

        # Cloud names
        # -----------------------------------------------------------------------
        ax1.annotate(cloud_name.capitalize(),
                    xytext=(0.96, 0.94),
                    xy=(0.96, 0.94),
                    textcoords='axes fraction',
                    xycoords='axes fraction',
                    size=10,
                    color='k',
                    bbox=dict(boxstyle='square',
                              facecolor='w',
                              alpha=1),
                    horizontalalignment='right',
                    verticalalignment='top',
                    )

        # create new box for plot
        # -----------------------------------------------------------------------
        #map_axes[1] -= 0.3
        map_axes[1] += 0.3

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        fig.show()

'''
The main script
'''

def load_results(filename, load_fits=True):

    import pickle
    from astropy.io import fits

    with open(filename, 'rb') as input:
        results = pickle.load(input)
    input.close()

    if load_fits:
        results['data']['av_data'], results['data']['av_header'] = \
                fits.getdata(results['filenames']['av_filename'],
                             header=True)
        results['data']['av_error_data'], results['data']['av_error_header'] = \
                fits.getdata(results['filenames']['av_error_filename'],
                             header=True)
        results['data']['av_data_ref'] = \
                fits.getdata(results['filenames']['av_ref_filename'])
        results['data']['hi_data'], results['data']['hi_header'] = \
                fits.getdata(results['filenames']['hi_filename'], header=True)
        results['data']['hi_data_error'] = \
                fits.getdata(results['filenames']['hi_error_filename'])
        results['data']['co_data'], results['data']['co_header'] = \
                fits.getdata(results['filenames']['co_filename'], header=True)

    return results

def main():

    # define constants
    # ---------------------------------------------------------------------------
    DIR_FIGURES = '/d/bip3/ezbc/multicloud/figures/'
    DIR_RESULTS = '/d/bip3/ezbc/multicloud/data/python_output/bootstrap_results/'
    FILENAME_EXT = '_planck_noint_gaussrange_isotropic_bootstrap_results.pickle'
    FILENAME_PLOT_BASE_H2 = DIR_FIGURES + 'maps/h2_maps'
    FILENAME_PLOT_BASE = DIR_FIGURES + 'maps/hi_h2_maps'
    PLOT_FILETYPES = ['png', 'pdf']
    CLOUD_NAMES = ['california', 'perseus', 'taurus']

    # get the data needed to plot
    # ---------------------------------------------------------------------------
    plot_dict = {}
    for cloud_name in CLOUD_NAMES:
        plot_dict[cloud_name] = {}
        cloud_dict = plot_dict[cloud_name]

        # load the analysis
        results = load_results(DIR_RESULTS + cloud_name + FILENAME_EXT)

        # Write data to dictionary for the plotting function
        # -----------------------------------------------------------------------
        hi_sd = results['data_products']['hi_sd']
        h2_sd = results['data_products']['h2_sd']

        cloud_dict['contour_image'] = None# hi_sd
        cloud_dict['contours'] = [4, 8]
        cloud_dict['hi_sd'] = hi_sd
        cloud_dict['h2_sd'] = h2_sd
        cloud_dict['header'] = results['data']['av_header']

        # set cloud limits
        # -----------------------------------------------------------------------
        if cloud_name == 'california':
            plot_dict[cloud_name]['limits'] = [74, 60, 30, 38]
        if cloud_name == 'perseus':
            plot_dict[cloud_name]['limits'] = [61, 46, 24, 36]
        if cloud_name == 'taurus':
            plot_dict[cloud_name]['limits'] = [74, 57, 20, 33]

    # plot the 3-panel H2 surface density map
    # ---------------------------------------------------------------------------
    for filetype in PLOT_FILETYPES:
        plot_hi_h2_sd(plot_dict,
                      filename=FILENAME_PLOT_BASE + '.' + filetype,
                      hi_vlimits=[5, 13],
                      h2_vlimits=[-0.1, 50],
                      )
        plot_h2_sd(plot_dict,
                   filename=FILENAME_PLOT_BASE_H2 + '.' + filetype,
                   vlimits=[-0.1, 50],
                   )

if __name__ == '__main__':
    main()

