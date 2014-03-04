#!/usr/bin/python

''' Calculates the N(HI) map for Taurus

'''

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

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec).
    '''

    import pywcsgrid2 as wcs
    import pywcs

    # convert to degrees
    ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)

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
    from matplotlib.patches import Rectangle

    # Create figure instance
    fig = plt.figure(figsize=(8,8))

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

    ax.set_xlabel('Right Ascension (J2000)',
              size = 'small',
              family='serif')
    ax.set_ylabel('Declination (J2000)',
              size = 'small',
              family='serif')
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
    cb.set_label_text(r'N(HI) $\times$ 10$^{20}$ cm$^{-2}$',
                   size='small',
                   family='serif')

    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    for key in cores:
        pix_coords = wcs_header.wcs_sky2pix([cores[key]['wcs_position']], 0)[0]

        ax.scatter(pix_coords[0],pix_coords[1], color='w', s=200, marker='+',
                linewidths=2)

        ax.annotate(key,
                xy=[pix_coords[0], pix_coords[1]],
                xytext=(5,5),
                textcoords='offset points',
                color='w')

        if boxes:
        	box = cores[key]['box']
        	width = box[2] - box[0]
        	height = box[3] - box[1]
        	print width,height
        	rect = ax.add_patch(Rectangle(box[0:2],
                width, height, facecolor='none', edgecolor='k' ))

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

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    core_dir = output_dir + 'core_arrays/'

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
            velocity_axis=velocity_axis, velocity_range=[-5,15])

    if False:
        # trim hi_image to av_image size
        nhi_image_trim = np.ma.array(nhi_image, mask=av_image != av_image)

        plot_nhi_image(nhi_image=nhi_image_trim, header=hi_header,
                contour_image=av_image, contours=[5,10,15],
                savedir=figure_dir, filename='taurus_nhi_cores_map.png',
                show=True)

    cores = {'L1495':
                {'wcs_position': [15*(4+14/60.), 28+11/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [get_pix_coords(ra=(4,16,30.031),
                                        dec=(27,44,30),
                                        header=h),
                         get_pix_coords(ra=(4,5,20),
                                        dec=(28,28,33),
                                        header=h)]
                 },
             'L1495A':
                {'wcs_position': [15*(4+18/60.), 28+23/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [get_pix_coords(ra=(4,28,23),
                                        dec=(28,12,50),
                                        header=h),
                         get_pix_coords(ra=(4,16,23),
                                        dec=(29,46,5),
                                        header=h)],
                 },
             'B213':
                {'wcs_position': [15*(4+19/60.), 27+15/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [get_pix_coords(ra=(4,22,27),
                                        dec=(26,45,47),
                                        header=h),
                         get_pix_coords(ra=(4,5,25),
                                        dec=(27,18,48),
                                        header=h)],
                },
             'B220':
                {'wcs_position': [15*(4+41/60.), 26+7/60., 0],
                 'map': None,
                 'threshold': 7,
                 'box': [get_pix_coords(ra=(4,47,49),
                                        dec=(25,31,13),
                                        header=h),
                         get_pix_coords(ra=(4,40,37),
                                        dec=(27,31,17),
                                        header=h)],
                 },
             'L1527':
                {'wcs_position': [15*(4+39/60.), 25+47/60., 0],
                 'map': None,
                 'threshold': 7,
                 'box': [get_pix_coords(ra=(4,40,13),
                                        dec=(24,46,38),
                                        header=h),
                         get_pix_coords(ra=(4,34,35),
                                        dec=(25,56,7),
                                        header=h)],
                 },
             'B215':
                {'wcs_position': [15*(4+23/60.), 25+3/60., 0],
                 'map': None,
                 'threshold': 3,
                 'box': [get_pix_coords(ra=(4,24,51),
                                        dec=(22,36,7),
                                        header=h),
                         get_pix_coords(ra=(4,20,54),
                                        dec=(25,26,31),
                                        header=h)],
                 },
             'L1524':
                {'wcs_position': [15*(4+29/60.), 24+31/60., 0],
                 'map': None,
                 'threshold': 3,
                 'box': [get_pix_coords(ra=(4,31,0),
                                        dec=(22,4,6),
                                        header=h),
                         get_pix_coords(ra=(4,25,33),
                                        dec=(25,0,55),
                                        header=h)],
                 }
                }

    # write out box parameter into single list
    for core in cores:
    	box = cores[core]['box']
    	cores[core]['box'] = (int(box[0][0]),int(box[0][1]),
    	        int(box[1][0]),int(box[1][1]))

    if False:
        nhi_image = np.zeros(nhi_image.shape)

        for core in cores:
        	core_image = np.load(core_dir + core + '.npy')
        	core_indices = np.where(core_image == core_image)
        	nhi_image[core_indices] += core_image[core_indices]

        nhi_image_trim =np.ma.array(nhi_image, mask=((av_image != av_image) &\
                (nhi_image == 0)))

        nhi_image_trim[nhi_image_trim == 0] = np.NaN

        plot_nhi_image(nhi_image=nhi_image_trim, header=hi_header,
            savedir=figure_dir,
            cores=cores,
            filename='taurus_nhi_core_regions_map.png',
            show=True)

    if True:
        # trim hi_image to av_image size
        nhi_image_trim = np.ma.array(nhi_image, mask = av_image != av_image)

        plot_nhi_image(nhi_image=nhi_image_trim, header=hi_header,
                contour_image=av_image, contours=[5,10,15],
                boxes=True, cores=cores, limits=[128,37,308,206],
                title='Taurus: N(HI) map with core boxed-regions.',
                savedir=figure_dir, filename='taurus_nhi_cores_map.png',
                show=True)

if __name__ == '__main__':
    main()


