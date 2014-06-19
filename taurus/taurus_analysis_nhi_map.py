#!/usr/bin/python

''' Calculates the N(HI) map for Taurus

'''

import numpy as np

def calculate_nhi(hi_cube=None, velocity_axis=None, velocity_range=[]):
    ''' Calculates an N(HI) image given a velocity range within which to
    include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to '''

    import numpy as np

    # Calculate NHI from cube if set
    if hi_cube is not None and velocity_axis is not None:
        image = np.empty((hi_cube.shape[1],
                          hi_cube.shape[2]))
        image[:,:] = np.NaN
        indices = np.where((velocity_axis > velocity_range[0]) & \
                (velocity_axis < velocity_range[1]))[0]
        image[:,:] = hi_cube[indices,:,:].sum(axis=0)

    # NHI in units of 1e20 cm^-2
    nhi_image = np.ma.array(image, mask=np.isnan(image)) * 1.823e-2

    return nhi_image

def convert_core_coordinates(cores, header):

    for core in cores:
        cores[core].update({'box_pixel': 0})
        cores[core].update({'center_pixel': 0})

    	box_wcs = cores[core]['box_wcs']
    	box_pixel = len(box_wcs) * [0,]
    	center_wcs = cores[core]['center_wcs']

        # convert centers to pixel coords
        cores[core]['center_pixel'] = get_pix_coords(ra=center_wcs[0],
                                                     dec=center_wcs[1],
                                                     header=header)
        # convert box corners to pixel coords
        for i in range(len(box_wcs)/2):
        	pixels = get_pix_coords(ra=box_wcs[2*i], dec=box_wcs[2*i + 1],
        	        header=header)
        	box_pixel[2*i], box_pixel[2*i + 1] = int(pixels[0]), int(pixels[1])
        cores[core]['box_pixel'] = box_pixel

    return cores

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec).
    '''

    import pywcsgrid2 as wcs
    import pywcs

    # convert to degrees
    if type(ra) is tuple and type(dec) is tuple:
        ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)
    else:
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

def plot_nhi_image(nhi_image=None, header=None, contour_image=None,
        cores=None, title=None, limits=None,
        contours=None, boxes=False, savedir='./', filename=None, show=True):

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
    im = ax.imshow(nhi_image,
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
    cb.set_label_text(r'N(HI) $\times$ 10$^{20}$ cm$^{-2}$',)

    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    for core in cores:
        pix_coords = cores[core]['center_pixel']

        ax.scatter(pix_coords[0],pix_coords[1], color='k', s=200, marker='+',
                linewidths=2)

        ax.annotate(core,
                xy=[pix_coords[0], pix_coords[1]],
                xytext=(5,5),
                textcoords='offset points',
                color='k')

        if boxes:
            rect = ax.add_patch(Polygon(
                cores[core]['box_vertices'][:, ::-1],
                    facecolor='none',
                    edgecolor='k' ))

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        fig.show()

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
    region = pyr.open(filename)

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
        box_center_pixel = get_pix_coords(ra = region[0],
                                          dec = region[1],
                                          header = header)
        box_center_pixel = (int(box_center_pixel[1]), int(box_center_pixel[0]))
        box_height = region[2] / header['CDELT1']
        box_width = region[3] / header['CDELT2']
        cores[core].update({'box_center_pix': box_center_pixel})
        cores[core].update({'box_width': box_width})
        cores[core].update({'box_height': box_height})
        cores[core].update({'box_angle': region[4]})

    return cores

def main():
    ''' Executes script.
    '''

    # import external modules
    import pyfits as pf
    import numpy as np
    from mycoords import make_velocity_axis
    import mygeometry as myg
    reload(myg)

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/maps/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    core_dir = output_dir + 'core_arrays/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'

    # Load hi fits file
    hi_image, hi_header = pf.getdata(hi_dir + \
            'taurus_galfa_cube_bin_3.7arcmin.fits', header=True)
    h = hi_header

    # Load av fits file
    av_image, av_header = pf.getdata(av_dir + 'taurus_av_k09_regrid.fits',
            header=True)

    # make velocity axis for hi cube
    velocity_axis = make_velocity_axis(hi_header)

    # create nhi image
    nhi_image = calculate_nhi(hi_cube=hi_image,
            velocity_axis=velocity_axis, velocity_range=[-100,100])

    if False:
        # trim hi_image to av_image size
        nhi_image_trim = np.ma.array(nhi_image, mask=av_image != av_image)

        plot_nhi_image(nhi_image=nhi_image_trim, header=hi_header,
                contour_image=av_image, contours=[5,10,15],
                savedir=figure_dir, filename='taurus_nhi_cores_map.png',
                show=True)

    cores = {'L1495':
                {'center_wcs': [(4,14,0), (28, 11, 0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(4,16,30), (27,44,30), (4,5,20), (28,28,33)]
                 },
             'L1495A':
                {'center_wcs': [(4,18,0), (28,23., 0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(4,28,23),(28,12,50),(4,16,23),(29,46,5)],
                 },
             'B213':
                {'center_wcs': [(4, 19, 0), (27, 15,0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(4,22,27), (26,45,47),(4,5,25),(27,18,48)],
                },
             'B220':
                {'center_wcs': [(4, 41, 0.), (26,7,0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,47,49),(25,31,13),(4,40,37),(27,31,17)],
                 },
             'L1527':
                {'center_wcs': [(4, 39, 0.), (25,47, 0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,40,13), (24,46,38), (4,34,35), (25,56,7)],
                 },
             'B215':
                {'center_wcs': [(4, 23, 0), (25, 3, 0)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,24,51), (22,36,7), (4,20,54), (25,26,31)],
                 },
             'L1524':
                {'center_wcs': [(4,29,0.), (24,31.,0)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,31,0), (22,4,6), (4,25,33), (25,0,55)],
                 }
                }

    cores = convert_core_coordinates(cores, h)

    if False:
        nhi_image = np.zeros(nhi_image.shape)

        for core in cores:
        	core_image = np.load(core_dir + core + '.npy')
        	core_indices = np.where(core_image == core_image)
        	nhi_image[core_indices] += core_image[core_indices]

        nhi_image_trim =np.ma.array(nhi_image, mask=((av_image != av_image) &\
                (nhi_image == 0)))

        nhi_image_trim[nhi_image_trim == 0] = np.NaN

        read_ds9_region(av_dir + 'taurus_av_boxes.reg')

        plot_nhi_image(nhi_image=nhi_image_trim, header=hi_header,
            savedir=figure_dir,
            cores=cores,
            filename='taurus_nhi_core_regions_map.png',
            show=True)

    if True:
        cores = load_ds9_region(cores,
                filename_base = region_dir + 'taurus_av_boxes_',
                header = h)

        # Grab the mask
        mask = np.zeros((nhi_image.shape))
        for core in cores:
        	xy = cores[core]['box_center_pix']
        	box_width = cores[core]['box_width']
        	box_height = cores[core]['box_height']
        	box_angle = cores[core]['box_angle']
        	mask += myg.get_rectangular_mask(nhi_image,
        	        xy[0], xy[1],
                    width = box_width,
                    height = box_height,
                    angle = box_angle)

        	cores[core]['box_vertices'] = myg.get_rect(
                        xy[0], xy[1],
                        width = box_width,
                        height = box_height,
                        angle = box_angle,)

        #print(cores[core]['box_vertices'])
        #print core, xy, box_width, box_height, box_angle

        mask[mask > 1] = 1

        #nhi_image[mask == 0] = np.nan

        # trim hi_image to av_image size
        nhi_image_trim = np.ma.array(nhi_image,
                mask = (av_image != av_image))

        # Plot
        figure_types = ['pdf', 'png']
        for figure_type in figure_types:
            plot_nhi_image(nhi_image=nhi_image_trim, header=hi_header,
                    contour_image=av_image, contours=[5,10,15],
                    boxes=True, cores = cores, limits=[128,37,308,206],
                    title='Taurus: N(HI) map with core boxed-regions.',
                    savedir=figure_dir, filename='taurus_nhi_cores_map.%s' % \
                            figure_type,
                    show=False)

if __name__ == '__main__':
    main()


