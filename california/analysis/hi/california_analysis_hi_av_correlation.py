#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the california molecular cloud.

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

def plot_nhi_vs_av(nhi_image, av_image, savedir='./', filename=None, show=True,
        returnimage=False, hess_binsize=None, title='', plot_type='hexbin'):
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
    nhi_image_nonans = nhi_image[(nhi_image == nhi_image) &\
                                 (av_image == av_image)&\
                                 (av_image > 0) &\
                                 (nhi_image > -5)]

    av_image_nonans = av_image[(nhi_image == nhi_image) &\
                               (av_image == av_image) &\
                               (av_image > 0) &\
                               (nhi_image > -5)]

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    if plot_type is 'hexbin':
        image = ax.hexbin(
                nhi_image_nonans.ravel(),
                av_image_nonans.ravel(),
                norm=matplotlib.colors.LogNorm(),
                mincnt=1,
                yscale='log',
                xscale='log')
        # Adjust color bar of density plot
        cb = plt.colorbar(image)
        cb.set_label('Counts')
    elif plot_type is 'scatter':
        image = ax.scatter(av_image_nonans.ravel(),
                nhi_image_nonans.ravel(),
                alpha=0.5,
                color='k'
                )
        ax.set_xscale('log')
        ax.set_yscale('log')

    # Adjust asthetics
    ax.set_ylabel('A$_v$ (mag)',
              size = 'small',
              family='serif')
    ax.set_xlabel(r'N(HI) (1 $\times 10^{20}$ cm$^{-2}$)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return correlations_image

def calculate_NHI(cube=None,velocity_axis=None,SpectralGrid=None,
        velocityrange=[]):
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
        indices = np.where((velocity_axis > velocityrange[0]) & \
                (velocity_axis < velocityrange[1]))[0]
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
                        if highVel < velocityrange[1] and \
                            lowVel > velocityrange[0]:
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
                (av_image < 2))# & \
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

def main(redo_cube_correlation_calculation = True):

    import grid
    import numpy as np

    # define directory locations
    output_dir = '/d/bip3/ezbc/california/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/california/figures/hi_velocity_range/'
    av_dir = '/d/bip3/ezbc/california/data/av/'
    hi_dir = '/d/bip3/ezbc/california/data/galfa/'

    # load Planck Av and GALFA HI images, on same grid
    av_data = load_fits(av_dir + \
                'california_planck_av_regrid.fits')

    x = np.linspace(0, av_data.shape[1] - 1, av_data.shape[1])
    y = np.linspace(0, av_data.shape[0] - 1, av_data.shape[0])

    x, y = np.meshgrid(x, y)

    mask = np.ones(av_data.shape)

    mask[(x > 0) & (x < 100) & (y > 100) & (y < 1000)] = 0

    av_data = np.ma.array(av_data,
            mask=((av_data < 2) & \
                  (av_data != av_data) & \
                  (mask == 1)))

    #av_data += - 0.4 # subtracts background of 0.4 mags
    hi_data,h = load_fits(hi_dir + 'california_galfa_cube_bin_3.7arcmin.fits',
            return_header=True)

    hi_data = np.ma.array(hi_data, mask=(hi_data != hi_data))

    # make the velocity axis
    velocity_axis = (np.arange(h['NAXIS3']) - h['CRPIX3'] + 1) * h['CDELT3'] + \
            h['CRVAL3']
    velocity_axis /= 1000.

    # define the parameters to derive NHI from the GALFA cube
    velocity_centers = np.arange(-20,20,0.5)
    velocity_widths = np.arange(1,80,5)

    if redo_cube_correlation_calculation:
        cube_correlations = calculate_correlation(cube=hi_data,
                velocity_axis=velocity_axis,
                av_image=av_data,
                velocity_centers=velocity_centers,
                velocity_widths=velocity_widths)
    else:
        try:
            cube_correlations = np.load(output_dir + 'cube_correlations.npy')
        except IOError:
            cube_correlations = None
    # save arrays
    if cube_correlations is not None:
        np.save(output_dir + 'cube_correlations', cube_correlations)
    np.save(output_dir + 'velocity_centers', velocity_centers)
    np.save(output_dir + 'velocity_widths', velocity_widths)

    if True:
        # Plot heat map of correlations


        cube_correlations_array = plot_correlations(cube_correlations,
                velocity_centers, velocity_widths, returnimage=True,
                savedir=figure_dir,
                filename='california.nhi_av_correlation.png',
                show=True)

        cube_correlations_array = np.ma.array(cube_correlations_array,
                mask=(cube_correlations_array != cube_correlations_array))
        print cube_correlations_array.max()
        print np.where(cube_correlations_array == cube_correlations_array.max())
        print cube_correlations_array.shape

        plot_center_velocity(cube_correlations,
                velocity_centers, velocity_widths,
                velocity_center=5, returnimage=True,
                savedir=figure_dir,
                filename='california.nhi_av_5kms_correlation_bkgrdSub.png',
                show=False)

        plot_center_velocity(cube_correlations,
                velocity_centers, velocity_widths,
                velocity_center=10, returnimage=True,
                savedir=figure_dir,
                filename='california.nhi_av_10kms_correlation_bkgrdSub.png',
                show=False)

    # Plot NHI vs. Av for a given velocity range
    #hi_data_corrected = np.ma.array(hi_data, mask=np.where(hi_data > -5))
    nhi_image = calculate_NHI(cube=hi_data,
            velocity_axis=velocity_axis,
        velocityrange=[4,8.5])

    if True:
        plot_nhi_vs_av(nhi_image,av_data,
                savedir=figure_dir,
                filename='california.av_nhi_2dDensity.png',)

    # calculate correlation for cores
    if False:
    	indices = [np.arange(130,160),np.arange(175,200)]

        nhi_image_sub = nhi_image[206:227,227:268]
        av_data_sub = av_data[206:227,227:268]

        plot_nhi_vs_av(nhi_image_sub,av_data_sub,
                savedir=figure_dir,
                plot_type='scatter',
                filename='california_av_vs_nhi_core1.png',
                title=r'N(HI) vs. A$_v$ of california Core L1495')

        nhi_image_sub = nhi_image[140:170,160:210]
        av_data_sub = av_data[140:170,160:210]

        plot_nhi_vs_av(nhi_image_sub,av_data_sub,
                plot_type='scatter',
                savedir=figure_dir,
                filename='california_av_vs_nhi_core2.png',)

    # Plot heat map of correlations
    if correlations is not None:
        correlations_array = plot_correlations(correlations,
                velocity_centers, velocity_widths, returnimage=True,show=False)
    if False:
        # Print best-fit characteristics
        indices = np.where(cube_correlations_array == \
                cube_correlations_array.max())
        print 'Maximum correlation values: '
        print str(velocity_centers[indices[0]][0]) + ' km/s center'
        print str(velocity_widths[indices[1]][0]) + ' km/s width'

if __name__ == '__main__':
    main()



