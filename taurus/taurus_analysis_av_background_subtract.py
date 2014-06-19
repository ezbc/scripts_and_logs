#!/usr/bin/python

''' Subtracts a polynomial fit from the Av image of Taurus.
'''


import numpy as np
import matplotlib.pyplot as plt


def plot_nhi_image(av_image=None, header=None, contour_image=None, title=None,
        limits=None, contours=None, boxes=False, savedir='./', filename=None,
        show=True):

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
              'figure.figsize': (6, 6),
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

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

    # create axes
    ax = imagegrid[0]
    cmap = cm.winter # colormap
    # show the image
    im = ax.imshow(av_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,)

    # Asthetics
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',)
    ax.set_ylabel('Declination (J2000)',)
    if title is not None:
    	ax.set_title(title)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Write label to colorbar
    cb.set_label_text(r'A$_V$ (mag)',)

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        fig.show()

def main():
    ''' Executes script.
    '''

    # import external modules
    import pyfits as pf
    import numpy as np
    from mycoords import make_velocity_axis
    import mygeometry as myg
    reload(myg)
    import mymath
    import myimage_analysis as myimage

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/maps/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    core_dir = output_dir + 'core_arrays/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'

    # Load av fits file
    av_image, av_header = pf.getdata(av_dir + 'taurus_planck_av_regrid.fits',
            header=True)

    fit_regions = ((51, 50, 82, 81),
                   (167, 52, 186, 79),
                   (229, 257, 263, 289),
                   (138, 332, 175, 372),
                   (34, 309, 69, 344))
    av_image_backsub = myimage.subtract_background(av_image,
                                                  degree=1,
                                                  fit_regions=fit_regions)

    figure_types = ['pdf', 'png']
    for figure_type in figure_types:
        plot_av_image(av_image=av_image_backsub,
                header=av_header,
                #contour_image=av_image,
                #contours=[5,10,15],
                #limits=[128,37,308,206],
                title=r'Taurus: A$_V$ map with background sub.',
                savedir=figure_dir,
                filename='taurus_av_back_sub_map.%s' % figure_type,
                show=False,
                )

if __name__ == '__main__':
    main()


