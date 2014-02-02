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

def plot_nhi_image(nhi_image=None, header=None, contour_image=None,
        contours=None, savedir='./', filename=None, show=True):

    # Import external modules
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps

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

    ax = imagegrid[0]
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)', 
              size = 'small', 
              family='serif')
    ax.set_ylabel('Declination (J2000)', 
              size = 'small',
              family='serif')

    #plt.rc('text', usetex=True)
    im = ax.imshow(nhi_image, interpolation='nearest',origin='lower',
            cmap=cm.gray)
    cb = ax.cax.colorbar(im)

    # Write label to colorbar
    cb.set_label_text(r'N(H\,\textsc{i}) $\times$ 10$^{20}$ cm$^{-2}$',
                   size='small',
                   family='serif')

    ax.contour(contour_image,colors='r')

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

    # Load hi fits file
    hi_image, hi_header = pf.getdata(hi_dir + \
            'taurus.galfa.cube.bin.4arcmin.fits', header=True)

    # Load av fits file
    av_image, av_header = pf.getdata(av_dir + 'taurus_av_k09_regrid.fits',
            header=True)

    # make velocity axis for hi cube
    velocity_axis = make_velocity_axis(hi_header)

    # create nhi image
    nhi_image = calculate_nhi(hi_cube=hi_image,
            velocity_axis=velocity_axis, velocity_range=[-10,20])

    # trim hi_image to av_image size
    nhi_image_trim = np.ma.array(nhi_image, mask=av_image != av_image)

    plot_nhi_image(nhi_image=nhi_image_trim, header=hi_header,
            contour_image=av_image, contours=[2,4,6,8,10],
            savedir=figure_dir, filename='taurus_nhi_map.png',
            show=True)

if __name__ == '__main__':
    main()


