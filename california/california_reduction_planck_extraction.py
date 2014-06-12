#!/usr/bin/python

import planckpy as pl
reload(pl)
from astropy.io import fits as pf

def extract_data(datatype = 'ebv'):
    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

    ra_range = (55, 75)
    dec_range = (20, 40)

    if datatype == 'ebv':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 2,
                dr_version = 2,
                resolution = 0.01,
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
                resolution = 0.01,
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
                resolution = 0.01,
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
    av_dir = '/d/bip3/ezbc/california/data/av/'
    planck_dir = '/d/bip3/ezbc/california/data/planck/'

    (data, header) = extract_data(datatype = 'tau353')
    write_data(data, header, planck_dir + 'california_planck_tau353.fits')
    (data, header) = tau353_to_ebv(data, header)
    write_data(data, header, planck_dir + 'california_planck_ebv.fits')
    (data, header) = ebv2av(data, header)
    write_data(data, header, av_dir + 'california_planck_av.fits')

    (data, header) = extract_data(datatype = 'tau353err')
    write_data(data, header, planck_dir + 'california_planck_tau353_error.fits')
    (data, header) = tau353_to_ebv(data, header)
    write_data(data, header, planck_dir + 'california_planck_ebv_error.fits')
    (data, header) = ebv2av(data, header)
    write_data(data, header, av_dir + 'california_planck_av_error.fits')


if __name__ == '__main__':
	main()


