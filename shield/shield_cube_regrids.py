#!/usr/bin/python

# Planck resolution is 5'. We will smooth and regrid all maps to the Planck map.

# Map                           Resolution      Convol size
# Lee et al. (2012) Av          5'              5
# GALFA HI                      3.7'            3.363
# CfA CO Dame et al. (2001)     8'              No smoothing

import os
import numpy as np

def main():

    from myimage_analysis import bin_image, calculate_nhi
    from mycoords import make_velocity_axis
    from astropy.io import fits
    import numpy as np

    os.chdir('/d/bip3/ezbc/shield/')

    # If true, deletes files to be written
    clobber = 1

    # First, change zeros in lee image to nans
    in_images = ('749237_cube.fits',)

    # Load the images into miriad
    out_images = []
    for in_image in in_images:

        print('Binning cube:\n' + in_image)

        if 1:
            cube, header = fits.getdata(in_image, header=True)

            # set freq0 setting
            header['FREQ0'] = 1.4204058E+09
            header['RESTFREQ'] = 1.4204058E+09
            header['CTYPE3'] = 'VELO'

            beamsize = header['BMAJ']
            cdelt = np.abs(header['CDELT1'])
            binsize = int(beamsize / cdelt)

            print('\tBinsize = ' + str(binsize))

            cube_bin, header_bin = bin_image(cube,
                                             binsize=(1, binsize, binsize),
                                             header=header,
                                             statistic=np.nanmean,
                                             )

            if 0:
                print header_bin['CDELT1']
                print header_bin['CDELT2']
                print header_bin['CRVAL1']
                print header_bin['CRVAL2']

            fits.writeto(in_image.replace('cube.fits', 'cube_regrid.fits'),
                         cube_bin,
                         header_bin,
                         clobber=clobber)
        else:
            cube_bin, header_bin = \
                fits.getdata(in_image.replace('cube.fits', 'cube_regrid.fits'),
                             clobber=clobber, header=True)

        # make nhi_image
        velocity_axis = make_velocity_axis(header_bin)

        # convert to T_B
        cube_bin_tb = 1.36 * 21 * cube_bin * 1000.0 / \
                      (header_bin['BMAJ'] * 3600.) / \
                      (3600. * header_bin['BMIN'])

        nhi_image = calculate_nhi(cube_bin_tb,
                                  velocity_axis=velocity_axis,
                                  header=header_bin,
                                  fits_filename=\
                            in_image.replace('cube.fits', 'nhi_regrid.fits')  )



if __name__ == '__main__':
    main()
