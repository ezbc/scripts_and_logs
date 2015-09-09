#!/usr/bin/python

import numpy as np
import planckpy as pl
reload(pl)
from astropy.io import fits as pf

def extract_data(datatype = 'ebv'):
    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

    ra_range = (22.5, 82.5)
    dec_range = (20, 38)

    if datatype == 'co_1to0_type1':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type1',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 0,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_1to0_error_type1':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type1',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 1,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_2to1_type1':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type1',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 4,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_2to1_error_type1':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'CO-Type1',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 5,
                dr_version = 1,
                resolution = 0.01,
                cut_last_pixel = False,
                verbose = True)
    elif datatype == 'co_3to2_type1':
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
    elif datatype == 'co_3to2_error_type1':
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
    elif datatype == 'co_1to0_type2':
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
    elif datatype == 'co_1to0_error_type2':
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
    elif datatype == 'co_2to1_type2':
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
    elif datatype == 'co_2to1_error_type2':
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
    elif datatype == 'co_type3':
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
    elif datatype == 'co_error_type3':
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
    elif datatype == 'ebv_err':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = 'Dust Opacity',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 3,
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

    if datatype == 'temp':
        ra_range = (40, 70)
        dec_range = (20, 38)

        (data, header) = pl.get_data(data_location = data_location,
            data_type = 'Dust Opacity',
            x_range = ra_range,
            y_range = dec_range,
            coord_type = 'equatorial',
            field = 4,
            resolution = 0.01,
            cut_last_pixel = False,
            dr_version = 2,
            verbose = True)

    if datatype == 'temp_error':
        ra_range = (40, 70)
        dec_range = (20, 38)
        (data, header) = pl.get_data(data_location = data_location,
            data_type = 'Dust Opacity',
            x_range = ra_range,
            y_range = dec_range,
            coord_type = 'equatorial',
            field = 5,
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

def ebv2av(data, header, error=None):
    # Av=E(B-V) x 3.05263, 152 in Tielens' ISM book.

    header['BUNIT'] = 'mag'

    rv = 3.05363
    av = data * rv

    if error is None:
        return (av, header)
    else:
        rv_data = np.loadtxt('/usr/users/ezbc/research/scripts/multicloud/' + \
                             'analysis/clouds/cloud_errors.csv',
                             skiprows=1, delimiter=',')

        rv_rms = np.mean(rv_data[:,0]**2)**0.5
        rv_std = np.std(rv_data[:,0])

        av_error = av * ((rv_std / rv)**2 + (error / data)**2)**0.5

        return (av_error, header)

def tau353_to_ebv(data, header, error=None):

    '''

    The dust emission maps are derived from:
    Planck 2013 results. XI. All-sky model of thermal dust emission

    Whereby they fit for the dust optical depth at 353 GHz, and then calibrated
    the conversion of the optical depth to color excess using background
    quasars.

    E(B-V) can be calculated from the dust optical depth at 353 GHz by
    (E(B - V)/tau 353 = 1.49 x 10**4

    '''

    header['BUNIT'] = 'mag'

    ebv = data * 1.49e4

    if error is None:
        return (ebv, header)
    else:
        ebv_error = ebv * ((0.03e4 / 1.49e4)**2 + (error / data)**2)**0.5
        return (ebv_error, header)

def main():
    av_dir = '/d/bip3/ezbc/perseus/data/min/'
    co_dir = '/d/bip3/ezbc/perseus/data/min/'
    planck_dir = '/d/bip3/ezbc/perseus/data/min/'

    if 0:


        (data, header) = extract_data(datatype = 'ebv')
        write_data(data, header, planck_dir + \
            'perseus_ebv_planck_radiance_min.fits')

        (data, header) = extract_data(datatype = 'ebv_error')
        write_data(data, header, planck_dir + \
            'perseus_ebv_planck_radiance_error_min.fits')

    if 0:
        (data, header) = extract_data(datatype = 'temp_error')
        write_data(data, header, planck_dir + \
            'perseus_planck_temp_error_min.fits')

        (data, header) = extract_data(datatype = 'temp')
        write_data(data, header, planck_dir + \
            'perseus_planck_temp_min.fits')
    if 0:
        # Color excess maps
        # -----------------
        # tau_353
        (data_tau353, header) = extract_data(datatype = 'tau353')
        write_data(data_tau353, header,
                   planck_dir + 'perseus_planck_tau353_min.fits')

        (data_ebv_tau353, header) = tau353_to_ebv(data_tau353, header)
        write_data(data_ebv_tau353, header,
                   planck_dir + 'perseus_ebv_planck_tau353_min.fits')

        (data_av_tau353, header) = ebv2av(data_ebv_tau353, header)
        write_data(data_av_tau353, header,
                   av_dir + 'perseus_av_planck_tau353_min.fits')
    if 0:

        # tau 353 error
        (data_tau353_error, header) = extract_data(datatype = 'tau353err')
        write_data(data_tau353_error, header,
                   planck_dir + 'perseus_planck_tau353_error_min.fits')

        (data_ebv_tau353_error, header) = tau353_to_ebv(data_tau353,
                                                        header,
                                                        error=data_tau353_error)
        write_data(data_ebv_tau353, header,
                   planck_dir + 'perseus_ebv_error_planck_tau353_min.fits')

        (data_av_tau353_error, header) = ebv2av(data_ebv_tau353, header,
                                                error=data_ebv_tau353_error)
        write_data(data_av_tau353_error, header,
                   av_dir + 'perseus_av_error_planck_tau353_min.fits')
    if 1:
        # CO maps
        # -----------------
        # type 1
        (data, header) = extract_data(datatype = 'co_1to0_type1')
        write_data(data, header, co_dir + 'perseus_co_type1_1-0_planck.fits')

        (data, header) = extract_data(datatype = 'co_1to0_error_type1')
        write_data(data, header, co_dir + 'perseus_co_type1_1-0_error_planck.fits')

        (data, header) = extract_data(datatype = 'co_2to1_type1')
        write_data(data, header, co_dir + 'perseus_co_type1_2-1_planck.fits')

        (data, header) = extract_data(datatype = 'co_2to1_error_type1')
        write_data(data, header, co_dir + 'perseus_co_type1_2-1_error_planck.fits')

        (data, header) = extract_data(datatype = 'co_3to2_type1')
        write_data(data, header, co_dir + 'perseus_co_type1_3-2_planck.fits')

        (data, header) = extract_data(datatype = 'co_3to2_error_type1')
        write_data(data, header, co_dir + 'perseus_co_type1_3-2_error_planck.fits')
        # type 2
        (data, header) = extract_data(datatype = 'co_1to0_type2')
        write_data(data, header, co_dir + 'perseus_co_type2_1-0_planck.fits')

        (data, header) = extract_data(datatype = 'co_1to0_error_type2')
        write_data(data, header, co_dir + 'perseus_co_type2_1-0_error_planck.fits')

        (data, header) = extract_data(datatype = 'co_2to1_type2')
        write_data(data, header, co_dir + 'perseus_co_type2_2-1_planck.fits')

        (data, header) = extract_data(datatype = 'co_2to1_error_type2')
        write_data(data, header, co_dir + 'perseus_co_type2_2-1_error_planck.fits')


        # type 3
        (data, header) = extract_data(datatype = 'co_type3')
        write_data(data, header, co_dir + 'perseus_co_type3_planck.fits')

        (data, header) = extract_data(datatype = 'co_error_type3')
        write_data(data, header, co_dir + 'perseus_co_type3_error_planck.fits')


if __name__ == '__main__':
	main()


