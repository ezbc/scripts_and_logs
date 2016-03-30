
#!/usr/bin/python

import pickle
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import myimage_analysis as myia
import mygeometry as myg
import mystats


def write_rad_field():

    from myscience import calc_radiation_field

    # Define filenames
    # ------------
    DIR_DUST = '/d/bip3/ezbc/multicloud/data/dust_temp/'
    DIR_RAD = '/d/bip3/ezbc/multicloud/data/radiation_field/'
    FILENAME_TEMPERATURE = DIR_DUST + 'multicloud_dust_temp_5arcmin.fits'
    FILENAME_TEMPERATURE_ERROR = DIR_DUST + \
            'multicloud_dust_temp_error_5arcmin.fits'
    FILENAME_BETA = DIR_DUST + 'multicloud_dust_beta_5arcmin.fits'
    FILENAME_BETA_ERROR = DIR_DUST + 'multicloud_dust_beta_error_5arcmin.fits'
    FILENAME_RAD = DIR_RAD + 'multicloud_rad_field_5arcmin.fits'
    FILENAME_RAD_ERROR = DIR_RAD + 'multicloud_rad_field_error_5arcmin.fits'

    # Get the data
    # ------------
    temp_data, temp_header = fits.getdata(FILENAME_TEMPERATURE, header=True)
    temp_error_data = fits.getdata(FILENAME_TEMPERATURE_ERROR)
    beta_data, beta_header = fits.getdata(FILENAME_BETA, header=True)
    beta_error_data = fits.getdata(FILENAME_BETA_ERROR)

    # Calculate the radiation field
    # -----------------------------
    rad_field, rad_field_error = calc_radiation_field(temp_data,
                                                      T_dust_error=temp_error_data,
                                                      beta=beta_data,
                                                      beta_error=beta_error_data)

    # convert Mathis field to Draine field
    rad_field /= 1.48

    # create rad field header
    rad_field_header = beta_header.copy()
    rad_field_header['BUNIT'] = ''

    # Write the radiation field to fits file
    fits.writeto(FILENAME_RAD, rad_field, header=rad_field_header, clobber=True)
    fits.writeto(FILENAME_RAD_ERROR, rad_field_error, header=rad_field_header,
                 clobber=True,)

def main():

    write_rad_field()

if __name__ == '__main__':
    main()

