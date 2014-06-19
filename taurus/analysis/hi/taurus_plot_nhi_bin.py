#!/usr/bin/python

''' Plots histogram of GALFA N(HI) pixels of the Taurus molecular cloud.

'''

def calculate_nhi(cube=None,velocity_axis=None,SpectralGrid=None,
        velocity_range=[]):
    ''' Calculates an N(HI) image given a velocity range within which to
    include a SpectralGrid's components.

    Parameters
    ----------
    cube : array-like, optional
        Three dimensional array with velocity axis as 0th axis. Must specify
        a velocity_axis if cube is used.
    velocity_axis : array-like, optional
        One dimensional array containing velocities corresponding to

    Returns
    -------
    nhi_image : array-like
        N(HI) array in units of 1e20 cm^-2

    '''

    import numpy as np

    # Calculate NHI from cube if set
    if cube is not None and velocity_axis is not None:
        image = np.empty((cube.shape[1],
                          cube.shape[2]))
        image[:,:] = np.NaN
        indices = np.where((velocity_axis > velocity_range[0]) & \
                (velocity_axis < velocity_range[1]))[0]
        image[:,:] = cube[indices,:,:].sum(axis=0)

    # Calculate NHI from SpectralGrid if set
    elif SpectralGrid is not None:
        if len(velocity_range) != 2:
            sp = SpectralGrid.spectralAxis
            velocity_range = [min(sp),max(sp)]
        image = np.empty((SpectralGrid.get_imagexsize(),
                          SpectralGrid.get_imageysize()))
        image[:,:] = np.NaN
        region = SpectralGrid.region
        tilecount=0
        for i in xrange(region[0],region[2]):
            for j in xrange(region[1],region[3]):
                tile = SpectralGrid.get_tile(i,j,coords='pixel')
                if tile.ncomps > 0:
                    for k in range(tile.ncomps):
                        lowVel = tile.params[k][1] - tile.params[k][2]
                        highVel = tile.params[k][1] + tile.params[k][2]
                        if highVel < velocity_range[1] and \
                            lowVel > velocity_range[0]:
                            image[i,j] = 0.
                            image[i,j] += tile.compAreaList[k]
                else:
                    image[i,j] = np.NaN
                tilecount = tilecount + 1

    # NHI in units of 1e20 cm^-2
    NHI_image = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2

    return NHI_image

def make_velocity_axis(header):
    ''' Makes velocity axis given fits header. Assumes spectral axis is axis 3.
    Returns in units of km/s.
    '''

    # External modules
    import numpy as np

    # Make axis
    velocity_axis = (np.arange(header['NAXIS3']) - header['CRPIX3'] + 1) * + \
            header['CDELT3'] + header['CRVAL3']
    return velocity_axis / 1000.

def plot_nhi_histogram(image, range=None, title='',  savedir='./',
        filename=None, show=True, returnimage=False):

    ''' Plots histogram of pixel values.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    image_nonans = image[image == image]

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.hist(image_nonans,color='g',alpha=1,range=range,bins=50, normed=True)

    # Adjust asthetics
    ax.set_ylabel('Probability Density',
              size = 'small',
              family='serif')
    ax.set_xlabel(r'N(HI) (1 $\times\ 10^{20}$ cm$^{-2}$)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return None

def main():

    import grid
    import numpy as np

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'

    # load GALFA HI
    hi_data,hi_header = load_fits(hi_dir + 'taurus_galfa_paradis_bin.fits',
            return_header=True)
    # make HI velocity axis
    velocity_axis_hi = make_velocity_axis(hi_header)

    # define the parameters to derive NHI from the GALFA cube
    nhi_image = calculate_nhi(cube=hi_data, velocity_axis=velocity_axis_hi,
            velocity_range=[4,6])

    # Adjust scale for ease of plot
    #nhi_image /= 1e20

    # Plot correlation, similar to Figure 3 of Paradis et al. (2012)
    plot_nhi_histogram(nhi_image, savedir=figure_dir,
            filename='taurus_nhi_histogram.png',)

if __name__ == '__main__':
    main()


