#!/usr/bin/python

def calculate_nhi(cube=None, velocity_axis=None, SpectralGrid=None,
        velocity_range=[], return_nhi_error=True, noise_cube=None,
        velocity_noise_range=[-110,-90,90,100],
        Tsys=30.):
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

        # Calculate image error
        if return_nhi_error:
            image_error = np.empty((cube.shape[1],
                              cube.shape[2]))
            image_error[:,:] = np.NaN
            image_error[:,:] = (noise_cube[indices,:,:]**2).sum(axis=0)**0.5

    # NHI in units of 1e20 cm^-2
    nhi_image = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2

    if return_nhi_error:
        nhi_image_error = np.ma.array(image_error,
                mask=np.isnan(image_error)) * 1.823e-2
        return nhi_image, nhi_image_error
    else:
        return nhi_image

def calculate_noise_cube(cube=None, velocity_axis=None,
            velocity_noise_range=[-110,-90,90,110], header=None, Tsys=30.,
            filename=None):

    """ Calcuates noise envelopes for each pixel in a cube
    """

    import numpy as np
    import pyfits as pf

    noise_cube = np.zeros(cube.shape)
    for i in range(cube.shape[1]):
        for j in range(cube.shape[2]):
            profile = cube[:,i,j]
            #noise = calculate_noise(profile, velocity_axis,
            #        velocity_noise_range)
            noise = 0.1 # measured in line free region
            noise_cube[:,i,j] = calculate_noise_scale(Tsys,
                    profile, noise=noise)

    if filename is not None:
        pf.writeto(filename, noise_cube, header=header)

    return noise_cube

def calculate_noise(profile, velocity_axis, velocity_range):
    """ Calculates rms noise of Tile profile given velocity ranges.
    """
    import numpy as np

    velMin = [velocity_range[0],velocity_range[2]]
    velMax = [velocity_range[1],velocity_range[3]]

    std = 0
    for k in xrange(len(velMin)):
        noise_region = np.where((velocity_axis >= velMin[k]) & \
                        (velocity_axis <= velMax[k]))
        std += np.std(profile[noise_region])

    std /= 2.
    return std

def calculate_noise_scale(Tsys, profile, noise=None):
    """ Creates an array for scaling the noise by (Tsys + Tb) / Tsys
    """
    import numpy as np
    n = np.zeros(len(profile))
    for i, Tb in enumerate(profile):
        n[i] = (Tsys + Tb) / Tsys * noise

    return n

def get_noise(image,region):
    ''' Calculates std of region.
    '''

    import numpy as np

    sub_image = image[region[0]:region[2], region[1]:region[3]]

    return sub_image[sub_image == sub_image].std()

def get_sub_image(image,indices):
    ''' Returns sub-region of image.
    '''

    return image[indices[0]:indices[2],indices[1]:indices[3]]

def hrs2degs(ra=None, dec=None):
    ''' Ra and dec tuples in hrs min sec and deg arcmin arcsec.
    '''

    ra_deg = 15*(ra[0], ra[1]/60., ra[2]/3600.)
    dec_deg = (dec[0], dec[1]/60., dec[2]/3600.)

    return (ra_deg, dec_deg)

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in degrees.
    '''

    import pywcsgrid2 as wcs
    import pywcs

    ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)

    wcs_header = pywcs.WCS(header)
    pix_coords = wcs_header.wcs_sky2pix(ra_deg, deg_deg, 0)[0]

    return pix_coords

def print_properties(data_dict):

    ''' Prints statistics about each core.
    '''

    import numpy as np

    for data in data_dict:
        av_data = data_dict[data]['av_data']
        av_image = np.zeros(av_data.shape)
        indices = np.where((av_data == av_data) & (av_data > -1))
        av_image[indices] = av_data[indices]
        nhi_image = data_dict[data]['nhi_image']
        nhi_error_image = data_dict[data]['nhi_error_image']

        av_min = av_image.min()
        av_max = av_image.max()
        av_error = get_noise(av_image, data_dict[data]['av_noise_region'])

        nhi_min = nhi_image[nhi_image > 0].min()
        nhi_max = nhi_image.max()
        nhi_error = np.median(nhi_error_image)
        nhi_units = '10^20 cm^-2'

        print(data + ':')
        print('Image\t Units\t \tMin\t Max\t 1sigma error')
        print('Av\t mag\t\t %.2f \t %.2f \t %.2f' % (av_min, av_max, av_error))
        print('N(HI)\t %s\t %.2f \t %.2f \t %.2f \n' % \
                (nhi_units, nhi_min, nhi_max, nhi_error))

def print_core_properties(data_dict, cores):

    ''' Prints statistics about each core.
    '''

    import numpy as np

    for data in data_dict:
        if data == 'perseus':
            av_data = data_dict[data]['av_data']
            av_image = np.zeros(av_data.shape)
            indices = np.where((av_data == av_data) & (av_data > -1))
            av_image[indices] = av_data[indices]
            nhi_image = data_dict[data]['nhi_image']
            nhi_error_image = data_dict[data]['nhi_error_image']

            av_min = av_image.min()
            av_max = av_image.max()
            av_error = get_noise(av_image, data_dict[data]['av_noise_region'])

            nhi_min = nhi_image[nhi_image > 0].min()
            nhi_max = nhi_image.max()
            nhi_error = np.median(nhi_error_image)
            nhi_units = '10^20 cm^-2'

            print(data + ':')
            print('Image\t Units\t \tMin\t Max\t 1sigma error')
            print('Av\t mag\t\t %.2f \t %.2f \t %.2f' % \
                    (av_min, av_max, av_error))
            print('N(HI)\t %s\t %.2f \t %.2f \t %.2f \n' % \
                    (nhi_units, nhi_min, nhi_max, nhi_error))

        elif data is 'taurus':
            av_data = data_dict[data]['av_data']
            av_image = np.zeros(av_data.shape)
            indices = np.where((av_data == av_data) & (av_data > -1))
            av_image[indices] = av_data[indices]
            nhi_image = data_dict[data]['nhi_image']
            nhi_error_image = data_dict[data]['nhi_error_image']

            # caclulate traits for each core
            for core in cores:
            	region = cores[core]['box']
            	av_image_sub = get_sub_image(av_image, region)
                av_min = av_image_sub.min()
                av_max = av_image_sub.max()
                av_error = get_noise(av_image,
                        data_dict[data]['av_noise_region'])

            	nhi_image_sub = get_sub_image(nhi_image, region)
            	nhi_error_image_sub = get_sub_image(nhi_error_image, region)
                nhi_min = nhi_image_sub[nhi_image_sub > 0].min()
                nhi_max = nhi_image_sub.max()
                nhi_error = np.median(nhi_error_image_sub)
                nhi_units = '10^20 cm^-2'

                print(core + ':')
                print('Image\t Units\t \tMin\t Max\t 1sigma error')
                print('Av\t mag\t\t %.2f \t %.2f \t %.2f' % \
                        (av_min, av_max, av_error))
                print('N(HI)\t %s\t %.2f \t %.2f \t %.2f \n' % \
                        (nhi_units, nhi_min, nhi_max, nhi_error))

def main():
    ''' Executes script.
    '''

    # import external modules
    import pyfits as pf
    import numpy as np
    from mycoords import make_velocity_axis
    from os import system,path

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/' + \
            'taurus_perseus_comparison/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    perseus_dir = '/d/bip3/ezbc/perseus/'
    taurus_dir = '/d/bip3/ezbc/taurus/'
    av_dir = 'data/av/'
    hi_dir = 'data/galfa/'
    core_dir = output_dir + 'core_arrays/'

    # Load hi fits file
    data_dict = {'perseus':{
    	            'hi_file': perseus_dir + hi_dir + \
    	                    'perseus.galfa.cube.bin.4arcmin.fits',
    	            'hi_noise_file': perseus_dir + hi_dir + \
    	                    'perseus.galfa.cube.bin.4arcmin.noise.fits',
    	            'av_file': perseus_dir + 'data/2mass/' + \
    	                    '2mass_av_lee12_nocal_regrid.fits',
    	            'hi_data': None,
                    'hi_header': None,
                    'hi_vel_axis': None,
                    'av_data': None,
                    'av_header': None,
                    'nhi_image': None,
                    'nhi_error_image': None,
                    'av_noise_region': [54,66,83,103]},

                'taurus':{
    	            'hi_file': taurus_dir + hi_dir + \
    	                    'taurus.galfa.cube.bin.4arcmin.fits',
    	            'hi_noise_file': taurus_dir + hi_dir + \
    	                    'taurus.galfa.cube.bin.4arcmin.noise.fits',
    	            'av_file': taurus_dir + av_dir + \
    	                    'taurus_av_k09_regrid.fits',
    	            'hi_data': None,
                    'hi_header': None,
                    'hi_vel_axis': None,
                    'av_data': None,
                    'av_header': None,
                    'nhi_image': None,
                    'nhi_error_image': None,
                    'av_noise_region': [179,233,195,260]}}

    for data in data_dict:
        # Get hi data
    	data_dict[data]['hi_data'], data_dict[data]['hi_header'] = \
            pf.getdata(data_dict[data]['hi_file'], header=True)

        # Create hi velocity axis
        data_dict[data]['hi_vel_axis'] = \
                make_velocity_axis(data_dict[data]['hi_header'])

        # Calculate or load nhi noise cube
        noise_cube_filename = data_dict[data]['hi_noise_file']
        if not path.isfile(noise_cube_filename):
            noise_cube = calculate_noise_cube(cube = data_dict[data]['hi_data'],
                velocity_axis = data_dict[data]['hi_vel_axis'],
                velocity_noise_range=[-110,-90,90,110], Tsys=30.,
                filename = noise_cube_filename)
        else:
            noise_cube, h = load_fits(noise_cube_filename,
                return_header=True)

        # calculate N(HI) image
    	data_dict[data]['nhi_image'], data_dict[data]['nhi_error_image'] = \
    	    calculate_nhi(cube = data_dict[data]['hi_data'],
                velocity_axis = data_dict[data]['hi_vel_axis'],
                velocity_range=[-5,15], Tsys=30., noise_cube = noise_cube)

        # get av data
    	data_dict[data]['av_data'], data_dict[data]['av_header'] = \
            pf.getdata(data_dict[data]['av_file'], header=True)

    cores = {'L1495':
                {'wcs_position': [15*(4+14/60.), 28+11/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [206,242,244,287]},
             'L1495A':
                {'wcs_position': [15*(4+18/60.), 28+23/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [206,212,236,242]},
             'B213':
                {'wcs_position': [15*(4+19/60.), 27+15/60., 0],
                 'map': None,
                 'threshold': 4.75,
                 'box': [177,206,206,242]},
             'B220':
                {'wcs_position': [15*(4+41/60.), 26+7/60., 0],
                 'map': None,
                 'threshold': 7,
                 'box': [179,131,199,157]},
             'L1527':
                {'wcs_position': [15*(4+39/60.), 25+47/60., 0],
                 'map': None,
                 'threshold': 7,
                 'box': [165,152,179,172]},
             'B215':
                {'wcs_position': [15*(4+23/60.), 25+3/60., 0],
                 'map': None,
                 'threshold': 3,
                 'box': [143,209,177,243]},
             'L1524':
                {'wcs_position': [15*(4+29/60.), 24+31/60., 0],
                 'map': None,
                 'threshold': 3,
                 'box': [138,177,167,209]}}


    print_properties(data_dict)
    print_core_properties(data_dict, cores)

    return data_dict[data]['nhi_image']

if __name__ == '__main__':
    main()


