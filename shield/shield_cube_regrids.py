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

    os.chdir('/d/bip3/ezbc/shield/749237_lowres/')

    # If true, deletes files to be written
    clobber = 1

    # First, change zeros in lee image to nans
    in_images = ('749237_rebin_cube.fits',)

    # Load the images into miriad
    out_images = []
    for in_image in in_images:

        print('Binning cube:\n' + in_image)

        cube, header = fits.getdata(in_image, header=True)

        # set freq0 setting
        header['FREQ0'] = 1.4204058E+09
        header['RESTFREQ'] = 1.4204058E+09
        header['CTYPE3'] = 'VELO'

        beamsize = header['BMAJ']
        cdelt = np.abs(header['CDELT1'])
        binsize = int(beamsize / cdelt)

        print('\tBinsize = ' + str(binsize))

        # cube measurement error = 700 uJy/Beam = 0.7 mJy/Beam
        # cube flux calibration error = 10%
        # add errors quadratically
        cube_std = np.nanstd(cube[0, :, :])
        cube_error = ((0.1 * cube)**2 + cube_std**2)**0.5

        cube_bin, header_bin = bin_image(cube,
                                         binsize=(1, binsize, binsize),
                                         header=header,
                                         statistic=np.nanmean,
                                         )

        noise_func = lambda x: (1 / np.nansum(x**-2))**0.5
        cube_error_bin = bin_image(cube_error,
                                   binsize=(1, binsize, binsize),
                                   statistic=noise_func,
                                   )

        fits.writeto(in_image,
                     cube,
                     header,
                     clobber=clobber)

        fits.writeto(in_image.replace('cube', 'cube_error'),
                     cube_error,
                     header,
                     clobber=clobber)

        fits.writeto(in_image.replace('cube.fits', 'cube_regrid.fits'),
                     cube_bin,
                     header_bin,
                     clobber=clobber)

        fits.writeto(in_image.replace('cube', 'cube_error_regrid'),
                     cube_error_bin,
                     header_bin,
                     clobber=clobber)
        #else:
        #    cube_bin, header_bin = \
        #        fits.getdata(in_image.replace('cube.fits', 'cube_regrid.fits'),
        #                     clobber=clobber, header=True)

        # make nhi_image
        velocity_axis = make_velocity_axis(header_bin)

        # convert to T_B
        cube_bin_tb = 1.36 * 21 * cube_bin * 1000.0 / \
                      (header_bin['BMAJ'] * 3600.) / \
                      (3600. * header_bin['BMIN'])

        cube_tb = 1.36 * 21 * cube * 1000.0 / \
                      (header['BMAJ'] * 3600.) / \
                      (3600. * header['BMIN'])

        nhi_image = calculate_nhi(cube_bin_tb,
                                  velocity_axis=velocity_axis,
                                  header=header_bin,
                                  fits_filename=\
                            in_image.replace('cube.fits', 'nhi_regrid.fits')  )
        nhi_image = calculate_nhi(cube_tb,
                                  velocity_axis=velocity_axis,
                                  header=header,
                                  fits_filename=\
                            in_image.replace('cube.fits', 'nhi.fits')  )


if __name__ == '__main__':
    main()
