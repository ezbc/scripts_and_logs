#!/usr/bin/python

import numpy as np
import gausspy.gp as gp
import traingp as train
import pickle

def create_synthetic_data(data_dict, filename=None):

    print('\nCreating synthetic data...')

    RMS = data_dict['RMS']
    temp_range = [30, 9000]
    amp_range = [3 * RMS, 25 * RMS]
    nspectra = 200
    ncomponents = 5

    agd_data = \
        train.create_train_data(velocity_axis=data_dict['velocity_axis'],
                                temp_range=temp_range,
                                amp_range=amp_range,
                                mean_range=[10, 400],
                                ncomponents=ncomponents,
                                nspectra=nspectra,
                                rms=RMS)

    if 0:
        import matplotlib.pyplot as plt
        plt.close(); plt.clf()
        plt.plot(agd_data['x_values'][1],
                 agd_data['data_list'][1])
        plt.savefig('/d/bip3/ezbc/scratch/spectrum_synthetic.png')

    # Save data for batch decomposition
    if filename is not None:
        pickle.dump(agd_data, open(filename, 'w'))

    return agd_data

def train_decomposer(train_dict, filename=None, one_phase=0):

    import gausspy.gp as gp
    import time
    import pickle

    print('\nTraining on synthetic data...')

    g = gp.GaussianDecomposer()
    g.load_training_data(filename)

    if one_phase:

        #One phase training
        g.set('phase', 'one')
        g.set('alpha1', 5)
        g.set('SNR_thresh', 3.)
        g.set('SNR2_thresh', 3.)

        g.train(alpha1_initial=5.0,
                verbose=False,
                mode='conv',
                learning_rate=1.0,
                eps=1.0,
                MAD = 0.1)
    else:
        #Two phase training
        g.set('phase', 'two')
        g.set('alpha1', 3)
        g.set('alpha2', 8)
        g.set('SNR_thresh', [3.,3.])
        g.set('SNR2_thresh', [3.,0.])
        #g.set('SNR_thresh', [5.,5.])
        #g.set('SNR2_thresh', [5.,0.])

        g.train(alpha1_initial=5.0,
                alpha2_initial=7,
                plot=False,
                verbose=False,
                mode='conv',
                learning_rate=1.0,
                eps=1.0,
                MAD=0.1,
                )

    return g

def get_train_results(train_dict, filename_train=None,
        filename_train_decomposed=None, load=False):

    import os

    # load decomposed train data if exists, else perform decomposition
    if load:
        if os.path.isfile(filename_train_decomposed):
            g_train = pickle.load(open(filename_train_decomposed, 'r'))
            perform_decomposition = False
        else:
            perform_decomposition = True
    else:
        perform_decomposition = True

    # Run AGD on train data?
    if perform_decomposition:
        g_train = train_decomposer(train_dict, filename=filename_train)

        if filename_train_decomposed is not None:
            pickle.dump(g_train, open(filename_train_decomposed, 'w'))

    return g_train

def add_data(data_dict, filename=None):

    ''' Adds the cube, header, velocity axis, and positions of each pixel to
    the data_dict.

    '''

    from astropy.io import fits
    from mycoords import make_velocity_axis

    cube, header = fits.getdata(filename, header=True)
    velocity_axis = make_velocity_axis(header)

    data_dict['cube'] = cube
    data_dict['header'] = header
    data_dict['velocity_axis'] = velocity_axis

    # Add positions in cube
    add_positions_to_data(data_dict)

def crop_data(data_dict):

    cube = data_dict['cube']
    vel_axis = data_dict['velocity_axis']

    # Crop channels below 75 km/s
    mask_crop = (vel_axis < 50)
    data_dict['cube'] = cube[~mask_crop, :, :]
    data_dict['velocity_axis'] = vel_axis[~mask_crop]

def format_data(data_dict, filename=None):

    data_dict['data_list'] = []
    data_dict['x_values'] = []
    data_dict['errors'] = []
    channels = np.arange(0, data_dict['cube'].shape[0], 1)

    cube = data_dict['cube']
    vel_axis = data_dict['velocity_axis']

    # Calculate the RMS as the median RMS of all spectra above 350 km/s
    rms_indices = np.where(vel_axis > 350)[0]
    RMS = np.median(np.sqrt(np.nanmean(cube[rms_indices, :, :]**2)))
    if 0:
        print 'RMS cube shape', cube[rms_indices, :, :].shape
        print 'RMS = ', RMS

    data_dict['RMS'] = RMS
    data_dict['errors'] = RMS * np.ones(len(channels))

    #print '\nRMS'
    #print RMS

    #for i in xrange(cube.shape[1]):
    #    for j in xrange(cube.shape[2]):
    #data_dict['x_values'].append(channels)
    #data_dict['x_values'].append(vel_axis)
    #for i in xrange(100,104):
    for i in xrange(cube.shape[1]):
        for j in xrange(cube.shape[2]):
        #for j in xrange(100,104):
            data_dict['data_list'].append(data_dict['cube'][:, i, j])
            #data_dict['errors'].append(RMS * np.ones(len(channels)))
            #data_dict['errors'].append(RMS)

    data_dict['cube'] = None
    if filename is not None:
        pickle.dump(data_dict, open(filename, 'w'))

def get_data(load=True, filename_dict=None, filename_cube=None):

    if not load:
        print('\nFormatting data...')
        # Load the data
        data_dict = {}
        add_data(data_dict, filename=filename_cube)

        # Crop the cube velocities
        #crop_data(data_dict)

        # Format data for decomposition
        format_data(data_dict, filename=filename_dict)
    else:
        print('\nLoading data...')
        data_dict = pickle.load(open(filename_dict, 'r'))

    return data_dict

def add_positions_to_data(data_dict):

    ''' Adds WCS and pixel coordinates to data dictionary.
    '''

    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy import wcs
    from astropy.io import fits

    print('\nAdding positions to data...')

    # Get the data
    cube, header = data_dict['cube'], data_dict['header']

    # Adjust header to make astropy happy
    header['CUNIT3'] = 'm/s'
    header['CUNIT2'] = 'deg'
    header['CUNIT1'] = 'deg'
    header['CTYPE3'] = 'VOPT'
    header['SPECSYS'] = 'LSRK'

    # initialize positions array
    data_dict['positions'] = {}
    shape = (cube.shape[1] * cube.shape[2], 2)
    data_dict['positions']['wcs'] = np.empty(shape)
    data_dict['positions']['pix'] = np.empty(shape)

    # Create header object
    w = wcs.WCS(header)

    # add position for each spectrum
    count = 0
    for i in xrange(cube.shape[1]):
        for j in xrange(cube.shape[2]):
            coords_pix = np.array([i, j])
            data_dict['positions']['pix'][count] = coords_pix
            coords_wcs = w.wcs_pix2world(([j, i, 0],), 0)[0]
            data_dict['positions']['wcs'][count] = np.array(coords_wcs[:2])
            count += 1

    return data_dict

def main():

    ''' Script to decompose GASS HI near LMC and SMC using Autonomous
    Gaussian Decomposition. The steps to decomposition are:

        1) Create synthetic data with known Gaussian parameters expected to be
        in the data.

        2) Train AGD to fit for Gaussians with similar parameters as in the
        training data.

        3) Decompose the HI cube.

    '''

    import os
    from mydecomposition import get_decomposed_data, decompose_data

    #os.chdir('/d/bip3/ezbc/magellanic_stream/scripts/gausspy_decomp/')
    os.chdir('/home/ezbc/research/magellanic_stream/scripts/gausspy_decomp/')

    DIR_HI = '/d/bip3/ezbc/multicloud/data/hi/'
    FILENAME_DATA = '/d/bip3/ezbc/multicloud/data/decomposition/agd_multicloud_data.pickle'
    FILENAME_TRAIN = '/d/bip3/ezbc/multicloud/data/decomposition/agd_multicloud_train.pickle'
    FILENAME_TRAIN_DECOMPOSED = \
            '/d/bip3/ezbc/multicloud/data/decomposition/agd_multicloud_train_decomp.pickle'
    FILENAME_DECOMPOSED = '/d/bip3/ezbc/multicloud/data/decomposition/agd_multicloud_decomp.pickle'
    FILENAME_CUBE = '/d/bip3/ezbc/multicloud/data/hi/multicloud_hi_galfa_cube.fits'

    # Which steps to load?
    LOAD_DATA = 0
    LOAD_TRAIN = 0
    LOAD_DECOMP = 0

    # Load the HI data
    if not LOAD_DECOMP or not LOAD_TRAIN:
        data_dict = get_data(load=LOAD_DATA,
                             filename_dict=FILENAME_DATA,
                             filename_cube=FILENAME_CUBE)

    print data_dict['x_values']

    # 1) Create synthetic data with known Gaussian parameters expected to be
    # in the data.
    if not LOAD_DECOMP and not LOAD_TRAIN:
        train_dict = create_synthetic_data(data_dict,
                                           filename=FILENAME_TRAIN)
    else:
        train_dict = None

    # 2) Train AGD to fit for Gaussians with similar parameters as in the
    # training data.
    if not LOAD_DECOMP:
        g_train = \
            get_train_results(train_dict,
                              filename_train=FILENAME_TRAIN,
                              filename_train_decomposed=\
                                  FILENAME_TRAIN_DECOMPOSED,
                              load=LOAD_TRAIN
                              )

    # 3) Decompose the HI cube.
    # -----------------------------------------------------------------------
    print '\nGetting decomposition results...'
    if LOAD_DECOMP:
        results_dict = \
            get_decomposed_data(FILENAME_DATA,
                                filename_decomposed=FILENAME_DECOMPOSED,
                                load=LOAD_DECOMP,
                                )
    else:
        results_dict = \
            get_decomposed_data(FILENAME_DATA,
                                g_train=g_train,
                                data_dict=data_dict,
                                filename_decomposed=FILENAME_DECOMPOSED,
                                load=LOAD_DECOMP,
                                )

if __name__ == '__main__':
    main()

