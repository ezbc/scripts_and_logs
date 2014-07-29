#!/usr/bin/python

import planckpy as pl
reload(pl)
from astropy.io import fits as pf

def extract_data(datatype = 'ebv'):
    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

    # Goal: RA from 1:30 to 5:30hrs, Dec = 20 to 38deg

    dec_range = (15, 35)
    ra_range = (52.5, 82.5)

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
    elif datatype == 'co_1to0':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type2',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 0,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_1to0_error':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type2',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 1,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_2to1':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type2',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 4,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_2to1_error':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type2',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 5,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_3to2':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type1',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 8,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_3to2_error':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type1',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 9,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type3',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 0,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_error':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type3',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 1,
                dr_version = 1,
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

    '''

    Information on the CO data products can be found at:

        http://adsabs.harvard.edu/abs/2013arXiv1303.5073P

        Planck Collaboration, Ade, P.~A.~R., Aghanim, N., et al.\ 2013,
        arXiv:1303.5073


    '''

    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    co_dir = '/d/bip3/ezbc/taurus/data/co/'
    planck_dir = '/d/bip3/ezbc/taurus/data/planck/'

    if 0:
        (data, header) = extract_data(datatype = 'tau353')
        write_data(data, header, planck_dir + 'taurus_planck_tau353.fits')
        (data, header) = tau353_to_ebv(data, header)
        write_data(data, header, planck_dir + 'taurus_planck_ebv.fits')
        (data, header) = ebv2av(data, header)
        write_data(data, header, av_dir + 'taurus_av_planck.fits')

        (data, header) = extract_data(datatype = 'tau353err')
        write_data(data, header, planck_dir + 'taurus_planck_tau353_error.fits')
        (data, header) = tau353_to_ebv(data, header)
        write_data(data, header, planck_dir + 'taurus_planck_ebv_error.fits')
        (data, header) = ebv2av(data, header)
        write_data(data, header, av_dir + 'taurus_av_error_planck.fits')

        (data, header) = extract_data(datatype = 'co_1to0')
        write_data(data, header, co_dir + 'taurus_co_1-0_planck.fits')

        (data, header) = extract_data(datatype = 'co_1to0_error')
        write_data(data, header, co_dir + 'taurus_co_1-0_error_planck.fits')

        (data, header) = extract_data(datatype = 'co_2to1')
        write_data(data, header, co_dir + 'taurus_co_2-1_planck.fits')

        (data, header) = extract_data(datatype = 'co_2to1_error')
        write_data(data, header, co_dir + 'taurus_co_2-1_error_planck.fits')

        (data, header) = extract_data(datatype = 'co_3to2')
        write_data(data, header, co_dir + 'taurus_co_3-2_planck.fits')

        (data, header) = extract_data(datatype = 'co_3to2_error')
        write_data(data, header, co_dir + 'taurus_co_3-2_error_planck.fits')

    (data, header) = extract_data(datatype = 'co')
    write_data(data, header, co_dir + 'taurus_co_planck.fits')

    (data, header) = extract_data(datatype = 'co_error')
    write_data(data, header, co_dir + 'taurus_co_error_planck.fits')

if __name__ == '__main__':
	main()


