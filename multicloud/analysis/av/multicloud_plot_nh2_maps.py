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
    DIR_FIGURES = '/d/bip3/ezbc/multicloud/figures/'
    DIR_RESULTS = '/d/bip3/ezbc/multicloud/data/python_output/bootstrap_results/'
    FILENAME_EXT = '_planck_noint_gaussrange_isotropic_bootstrap_results.pickle'
    FILENAME_PLOT_BASE = DIR_FIGURES + 'maps/h2_maps'
    PLOT_FILETYPES = ['png', 'pdf']
    CLOUD_NAMES = ['california', 'perseus', 'taurus']

    # get the data needed to plot
    plot_dict = {}
    for cloud_name in CLOUD_NAMES:
        plot_dict[cloud_name] = {}
        cloud_dict = plot_dict[cloud_name]

        # load the analysis
        results = load_results(DIR_RESULTS + cloud_name + FILENAME_EXT)

        hi_sd = results['data_products']['hi_sd']
        h2_sd = results['data_products']['h2_sd']

        cloud_dict['contour_image'] = None# hi_sd
        cloud_dict['contours'] = [4, 8]
        cloud_dict['h2_sd'] = h2_sd
        cloud_dict['header'] = results['data']['av_header']

        if cloud_name == 'california':
            plot_dict[cloud_name]['limits'] = [75, 60, 30, 38]
        if cloud_name == 'perseus':
            plot_dict[cloud_name]['limits'] = [61, 45, 24, 36]
        if cloud_name == 'taurus':
            plot_dict[cloud_name]['limits'] = [75, 57, 20, 33]

    # plot the 3-panel H2 surface density map
    for filetype in PLOT_FILETYPES:
        plot_h2_sd(plot_dict,
                   filename=FILENAME_PLOT_BASE + '.' + filetype,
                   vlimits=[-0.1, 50],
                   )


if __name__ == '__main__':
    main()

