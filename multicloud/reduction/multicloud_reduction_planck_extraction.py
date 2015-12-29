#!/usr/bin/python

import planckpy as pl
reload(pl)
from astropy.io import fits as pf

def extract_data(datatype = 'ebv'):
    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

    dec_range = (15, 48)
    ra_range = (32, 85)

    if datatype == 'ebv':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 2,
                dr_version = 2,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'ebv_err':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 3,
                dr_version = 2,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'tau353':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 0,
                dr_version = 2,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'tau353err':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 1,
                dr_version = 2,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'temp':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 4,
                dr_version = 2,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'temp_error':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 5,
                dr_version = 2,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'beta':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 6,
                dr_version = 2,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'beta_error':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 7,
                dr_version = 2,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == '857':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = '857',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 1,
                dr_version = 1,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == '545':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = '545',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 1,
                dr_version = 1,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == '353':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = '353',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 1,
                dr_version = 1,
                resolution = 0.05,
                cut_last_pixel = False,
                verbose = True)

    return (data, header)

def write_data(data, header, filename):

    pf.writeto(filename,
            data,
            header = header,
            clobber = True,
            output_verify = 'fix')

def ebv2av(data, header):
    # Av=E(B-V) x 3.05263, 152 in Tielens' ISM book.

    data *= 3.05363
    header['BUNIT'] = 'mag'

    return (data, header)

def tau353_to_ebv(data, header):

    '''

    The dust emission maps are derived from:
    Planck 2013 results. XI. All-sky model of thermal dust emission

    Whereby they fit for the dust optical depth at 353 GHz, and then calibrated
    the conversion of the optical depth to color excess using background
    quasars.

    E(B-V) can be calculated from the dust optical depth at 353 GHz by
    (E(B - V)/tau 353 = 1.49 x 10**4

    '''

    data *= 1.49e4
    header['BUNIT'] = 'mag'

    return (data, header)

def main():
    av_dir = '/d/bip3/ezbc/multicloud/data/av/'
    temp_dir = '/d/bip3/ezbc/multicloud/data/dust_temp/'
    planck_dir = '/d/bip3/ezbc/multicloud/data/planck/'

    if 0:
    	# Color excess maps
        # -----------------
        (data, header) = extract_data(datatype = 'tau353')
        write_data(data, header, planck_dir + 'multicloud_planck_tau353.fits')
        (data, header) = tau353_to_ebv(data, header)
        write_data(data, header, planck_dir + 'multicloud_planck_ebv.fits')
        (data, header) = ebv2av(data, header)
        write_data(data, header, av_dir + 'multicloud_av_planck.fits')

        (data, header) = extract_data(datatype = 'tau353err')
        write_data(data, header, planck_dir + 'multicloud_planck_tau353_error.fits')
        (data, header) = tau353_to_ebv(data, header)
        write_data(data, header, planck_dir + 'multicloud_planck_ebv_error.fits')
        (data, header) = ebv2av(data, header)
        write_data(data, header, av_dir + 'multicloud_av_error_planck.fits')

    if 0:
        (data, header) = extract_data(datatype = 'ebv')
        write_data(data, header, planck_dir + 'multicloud_planck_radiance_ebv.fits')
        (data, header) = ebv2av(data, header)
        write_data(data, header, av_dir + 'multicloud_av_planck_radiance.fits')

        (data, header) = extract_data(datatype = 'ebv_err')
        write_data(data, header,
                planck_dir + 'multicloud_planck_radiance_ebv_error.fits')
        (data, header) = ebv2av(data, header)
        write_data(data, header,
                av_dir + 'multicloud_av_error_planck_radiance.fits')
    if 0:
        # Band maps
        # ---------
        (data, header) = extract_data(datatype = '857')
        write_data(data, header, planck_dir + 'multicloud_planck_857ghz.fits')

        (data, header) = extract_data(datatype = '545')
        write_data(data, header, planck_dir + 'multicloud_planck_545ghz.fits')

        (data, header) = extract_data(datatype = '353')
        write_data(data, header, planck_dir + 'multicloud_planck_353ghz.fits')

    if 0:
        # Temp maps
        # ---------
        (data, header) = extract_data(datatype = 'temp')
        write_data(data, header, temp_dir + 'multicloud_dust_temp.fits')

        (data, header) = extract_data(datatype = 'temp_error')
        write_data(data, header, temp_dir + 'multicloud_dust_temp_error.fits')

    if 1:
        # Beta maps
        # ---------
        (data, header) = extract_data(datatype = 'beta')
        write_data(data, header, temp_dir + 'multicloud_dust_beta.fits')

        (data, header) = extract_data(datatype = 'beta_error')
        write_data(data, header, temp_dir + 'multicloud_dust_beta_error.fits')

if __name__ == '__main__':
	main()


