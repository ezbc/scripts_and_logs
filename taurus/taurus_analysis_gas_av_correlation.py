#!/usr/bin/python

''' Calculates the N(H) / Av correlation for the Taurus molecular cloud. Uses
Xco factor from Paradis et al. (2012) A&A, 543, 103 to calculate N(H2).
Integrates GALFA HI image to determine N(HI).

'''

def plot_correlations(correlations,velocity_centers,velocity_widths,
        savedir='./', filename=None,show=True, returnimage=False):
    ''' Plots a heat map of correlation values as a function of velocity width
    and velocity center.
    '''

    # Import external modules
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8,2.5))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="0.3%",
                 cbar_size='6%',
                 axes_pad=0,
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    correlations_image = np.empty((velocity_centers.shape[0],
                                   velocity_widths.shape[0]))
    correlations_image[:,:] = np.NaN
    count = 0
    try:
        for i, center in enumerate(velocity_centers):
            for j, width in enumerate(velocity_widths):
                correlations_image[i,j] = correlations[count]
                count += 1
    except IndexError:
        print(' plot_correlations: O-d array input, cannot proceed')

    image = np.ma.array(correlations_image,mask=np.isnan(correlations_image))

    ax = imagegrid[0]
    ax.set_xlabel('Velocity Width (km/s)',
              size = 'small',
              family='serif')
    ax.set_ylabel('Velocity Center (km/s)',
              size = 'small',
              family='serif')

    #ax.set_xticks(np.arange(0,velocity_widths.shape[0],1)[::5],
    #        velocity_centers[::5])

    plt.rc('text', usetex=True)
    im = ax.imshow(image, interpolation='nearest',origin='lower',
            extent=[velocity_widths[0],velocity_widths[-1],
                velocity_centers[0],velocity_centers[-1]],
            cmap=plt.jet())
    cb = ax.cax.colorbar(im)

    #cb.set_clim(vmin=0.)
    # Write label to colorbar
    cb.set_label_text(r'Correlation coefficient',
                   size='small',
                   family='serif')

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_center_velocity(correlations, velocity_centers, velocity_widths,
        velocity_center=None, savedir='./', filename=None, show=True,
        returnimage=False):
    ''' Plots correlation vs. velocity width for a single center velocity.
    '''

    # Import external modules
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8,8))

    correlations_image = np.empty((velocity_centers.shape[0],
                                   velocity_widths.shape[0]))
    correlations_image[:,:] = np.NaN
    count = 0
    for i, center in enumerate(velocity_centers):
        for j, width in enumerate(velocity_widths):
            correlations_image[i,j] = correlations[count]
            count += 1

    ax = fig.add_subplot(111)
    ax.set_xlabel('Velocity Width (km/s)',
              size = 'small',
              family='serif')
    ax.set_ylabel('Correlation Coefficient',
              size = 'small',
              family='serif')
    ax.set_title('Center velocity ' + str(velocity_center) + ' km/s')

    center_index = np.where(velocity_centers == velocity_center)[0][0]

    ax.plot(velocity_widths,
            correlations_image[center_index,:])
    # fig.show()

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def plot_nh_vs_av(nh_image, av_image, savedir='./', filename=None, show=True,
        returnimage=False, hess_binsize=None):
    ''' Plots N(HI) as a function of Av for individual pixels in an N(HI) image
    and an Av image.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    # Drop the NaNs from the images
    nh_image_nonans = nh_image[(nh_image == nh_image) &\
                                 (av_image == av_image)&\
                                 (av_image > 0) &\
                                 (nh_image > 0)]

    av_image_nonans = av_image[(nh_image == nh_image) &\
                               (av_image == av_image) &\
                               (av_image > 0) &\
                               (nh_image > 0)]

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    image = ax.hexbin(nh_image_nonans.ravel(),
            av_image_nonans.ravel(),#bins='log'
            norm=matplotlib.colors.LogNorm(),
            mincnt=1,
            yscale='log',
            xscale='log',)
            #extent=[-2,1.1,-1,2.8])

    # Adjust asthetics
    ax.set_ylabel('A$_v$ (mag)',
              size = 'small',
              family='serif')
    ax.set_xlabel(r'N(H) ($1 \times\ 10^{21}$ cm$^{-2}$)',
              size = 'small',
              family='serif')
    ax.set_title(r'A$_{\rm v}$ vs. N(H) for Taurus')
    ax.grid(True)
    cb = plt.colorbar(image)
    cb.set_label('Counts')

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def calculate_nh(image_nh2=None, image_nhi=None):
    ''' Calculates nh by image_nh2*2 + image_nhi'''

    # Perform calculation
    return image_nh2 * 2. + image_nhi

def calculate_nh2(cube=None, Xco=None, velocity_axis=None, velocity_range=[],
        lower_limit=0.2):
    ''' Calculates an NH2 image from a CO cube given a Xco value.
    '''

    # External modules
    import numpy as np

    # Integrate CO cube if set
    if cube is not None and velocity_axis is not None:
        image = np.empty((cube.shape[1],
                          cube.shape[2]))
        image[:,:] = np.NaN
        indices = np.where((velocity_axis > velocity_range[0]) & \
                (velocity_axis < velocity_range[1]))[0]
        image[:,:] = cube[indices,:,:].sum(axis=0)

    # Derive N(H2) image with X-factor
    if type(Xco) is int or type(Xco) is float:
        nh2_image = np.ma.array(image,mask=((np.isnan(image)) & \
                (image < lower_limit))) * 2. * Xco
        return nh2_image
    else:
        return None

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
        One dimensional array containing velocities corresponding to '''

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
        if len(velocityrange) != 2:
            sp = SpectralGrid.spectralAxis
            velocityrange = [min(sp),max(sp)]
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

def calculate_correlation(SpectralGrid=None,cube=None,velocity_axis=None,
        av_image=None, velocity_centers=[], velocity_widths=[],av_noise=0.5,
        av_SNR=None):
    ''' Calculates the correlation coefficient between either a SpectralGrid or
    a cube and an Av image.

    Parameters
    ----------
    cube : array-like, optional

    '''

    # Import external modules
    import grid
    import numpy as np
    from scipy.stats import pearsonr

    # calculate the velocity ranges given a set of centers and widths
    velocity_ranges = np.zeros(shape=[len(velocity_centers) * \
            len(velocity_widths),2])
    count = 0
    for i, center in enumerate(velocity_centers):
        for j, width in enumerate(velocity_widths):
            velocity_ranges[count,0] = center - width/2.
            velocity_ranges[count,1] = center + width/2.
            count += 1

    # calculate the correlation coefficient for each velocity range
    correlations = np.zeros(velocity_ranges.shape[0])
    pvalues = np.zeros(velocity_ranges.shape[0])

    # calc
    for i in range(velocity_ranges.shape[0]):
        nhi_image_temp = calculate_NHI(SpectralGrid=SpectralGrid,
                cube=cube,
                velocity_axis=velocity_axis,
                velocityrange=velocity_ranges[i])
        nhi_image = np.ma.array(nhi_image_temp,
                                mask=np.isnan(nhi_image_temp))

        # Select pixels with Av > 1.0 mag and Av_SNR > 5.0.
        # Av > 1.0 mag is used to avoid too low Av.
        # 1.0 mag corresponds to SNR = 1 / 0.2 ~ 5
        # (see Table 2 of Ridge et al. 2006).
        indices = np.where((nhi_image == nhi_image) & \
                (av_image == av_image) & \
                (av_image > 1))# & \
                #(av_image > 5*av_noise))

        nhi_image_corr = nhi_image[indices]
        av_image_corr = av_image[indices] #- 0.8 # subtract background of 0.8
        # Use Pearson's correlation test to compare images
        correlations[i] = pearsonr(nhi_image_corr.ravel(),
                av_image_corr.ravel())[0]

        # Shows progress each 10%
        total = float(velocity_ranges.shape[0])
        abs_step = int((total * 1)/100) or 1
        if i and not i % abs_step:
                 print "{0:.2%} processed".format(i/total)
    return correlations

def load_fits(filename,return_header=False):
    ''' Loads a fits file.
    '''

    import pyfits as pf

    f = pf.open(filename)
    if return_header:
        return f[0].data,f[0].header
    else:
        return f[0].data

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

def main():

    import grid
    import numpy as np

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/hi_velocity_range/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    co_dir = '/d/bip3/ezbc/taurus/data/cfa/'

    # load 2mass Av from Kainulainen et al. (2009)
    av_image_jouni = load_fits(av_dir + 'taurus_av_k09_paradis_bin.fits')

    # load GALFA HI
    hi_data,hi_header = load_fits(hi_dir + 'taurus_galfa_paradis_bin.fits',
            return_header=True)
    # make HI velocity axis
    velocity_axis_hi = make_velocity_axis(hi_header)


    # load CfA CO cube
    co_data,co_header = load_fits(co_dir + 'taurus_cfa_paradis_bin.fits',
            return_header=True)
    # make CO velocity axis
    velocity_axis_co = make_velocity_axis(co_header)

    # define the parameters to derive NHI from the GALFA cube
    velocity_centers = np.arange(-20,20,0.5)
    velocity_widths = np.arange(1,120,5)
    #velocity_centers = np.arange(-40,40,5)
    #velocity_centers = np.array([5])
    #velocity_widths = np.arange(1,100,20)

    # Derive N(H2) image
    # Xco in units of cm^-2 / (K km/s)
    # Paradis+10 imposed low CO limit of 0.2 K km/s
    nh2_image = calculate_nh2(cube=co_data, Xco=1.67e20,
            velocity_axis=velocity_axis_co, velocity_range=[-50,50],
            lower_limit=0.2)

    # Derive N(HI) image
    nhi_image = calculate_nhi(cube=hi_data, velocity_axis=velocity_axis_hi,
            velocity_range=[-10,20])

    # Derive N(H) image
    nh_image = calculate_nh(image_nh2=nh2_image, image_nhi=nhi_image)

    # Adjust scale for ease of plot
    nh_image /= 1e21


    # Plot correlation, similar to Figure 3 of Paradis et al. (2012)
    plot_nh_vs_av(nh_image,av_image_jouni,
            savedir=figure_dir,
            filename='taurus_av_vs_nh_kainulainen2009.png',)

if __name__ == '__main__':
    main()


