#!/usr/bin/python

import planckpy as pl
reload(pl)
from astropy.io import fits as pf

def extract_data(datatype = 'ebv'):
    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

    ra_range = (22.5, 82.5)
    dec_range = (20, 38)

    if datatype == 'co-type3':
        (data, header) = pl.get_data(data_location = data_location,
        data_type = 'CO-Type3',
        x_range = ra_range,
        y_range = dec_range,
        coord_type = 'equatorial',
        field = 0,
        resolution = 0.01,
        cut_last_pixel = False,
        verbose = True)

    if datatype == 'co-type3_error':
        (data, header) = pl.get_data(data_location = data_location,
        data_type = 'CO-Type3',
        x_range = ra_range,
        y_range = dec_range,
        coord_type = 'equatorial',
        field = 1,
        resolution = 0.01,
        cut_last_pixel = False,
        verbose = True)

    if datatype == 'co-type1-j10':
        (data, header) = pl.get_data(data_location = data_location,
            data_type = 'CO-Type1',
            x_range = ra_range,
            y_range = dec_range,
            coord_type = 'equatorial',
            field = 0,
            resolution = 0.01,
            cut_last_pixel = False,
            verbose = True)

    if datatype == 'co-type1-j10_error':
        (data, header) = pl.get_data(data_location = data_location,
            data_type = 'CO-Type1',
            x_range = ra_range,
            y_range = dec_range,
            coord_type = 'equatorial',
            field = 1,
            resolution = 0.01,
            cut_last_pixel = False,
            verbose = True)

    if datatype == 'ebv':
        (data, header) = pl.get_data(data_location = data_location,
            data_type = 'Dust Opacity',
            x_range = ra_range,
            y_range = dec_range,
            coord_type = 'equatorial',
            field = 2,
            resolution = 0.01,
            cut_last_pixel = False,
            dr_version = 2,
            verbose = True)

    if datatype == 'ebv_error':
        (data, header) = pl.get_data(data_location = data_location,
            data_type = 'Dust Opacity',
            x_range = ra_range,
            y_range = dec_range,
            coord_type = 'equatorial',
            field = 3,
            resolution = 0.01,
            cut_last_pixel = False,
            dr_version = 2,
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
    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    planck_dir = '/d/bip3/ezbc/perseus/data/planck/'

    (data, header) = extract_data(datatype = 'co-type3')
    write_data(data, header, planck_dir + 'perseus_planck_cotype3_snez.fits')

    (data, header) = extract_data(datatype = 'co-type3_error')
    write_data(data, header, planck_dir + \
            'perseus_planck_cotype3_error_snez.fits')

    (data, header) = extract_data(datatype = 'co-type1-j10')
    write_data(data, header, planck_dir + \
            'perseus_planck_cotype1_j1-0_snez.fits')

    (data, header) = extract_data(datatype = 'co-type1-j10_error')
    write_data(data, header, planck_dir + \
            'perseus_planck_cotype1_j1-0_error_snez.fits')

    (data, header) = extract_data(datatype = 'ebv')
    write_data(data, header, planck_dir + \
        'perseus_planck_ebv_snez.fits')

    (data, header) = extract_data(datatype = 'ebv_error')
    write_data(data, header, planck_dir + \
        'perseus_planck_ebv_error_snez.fits')

if __name__ == '__main__':
	main()


