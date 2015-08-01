#!/usr/bin/python

import planckpy as pl
reload(pl)
from astropy.io import fits as pf

def extract_data(datatype = 'ebv'):
    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

    dec_range = (21, 38)
    ra_range = (40, 70)

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
    elif datatype == '857':
        (data, header) = pl.get_data(data_location = data_location,
                data_type = '857',
                x_range = ra_range,
                y_range = dec_range,
                coord_type = 'equatorial',
                field = 1,
                dr_version = 1,
                resolution = 0.01,
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
                resolution = 0.01,
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

def radiance_to_ebv(data, header):

    '''

    The dust emission maps are derived from:
    Planck 2013 results. XI. All-sky model of thermal dust emission

    Whereby they fit for the radiance, and then calibrated the conversion of
    the radiance to color excess using background quasars.

    E(B-V) can be calculated from the dust radiance by
    (E(B - V)/radiance = 5.40 x 10**5

    '''

    data *= 5.40e5
    header['BUNIT'] = 'mag'

    return (data, header)

def combine_ebv(ebv_radiance, ebv_tau353, Rv=3.1):

    ebv_combined = np.copy(ebv_radiance)
    ebv_combined[ebv_combined > 0.3] = ebv_tau353

    av = ebv_combined * Rv

    return av

def main():

    '''

    Information on the CO data products can be found at:

        http://adsabs.harvard.edu/abs/2013arXiv1303.5073P

        Planck Collaboration, Ade, P.~A.~R., Aghanim, N., et al.\ 2013,
        arXiv:1303.5073


    '''

    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    co_dir = '/d/bip3/ezbc/perseus/data/co/'
    planck_dir = '/d/bip3/ezbc/perseus/data/planck/'

    if 1:
        # Color excess maps
        # -----------------
        # tau_353
        (data, header) = extract_data(datatype = 'tau353')
        write_data(data, header, planck_dir + 'perseus_planck_tau353.fits')

        (data_ebv_tau353, header) = tau353_to_ebv(data, header)
        write_data(data_ebv_tau353, header,
                   planck_dir + 'perseus_ebv_planck_tau353.fits')

        (data_av_tau353, header) = ebv2av(data_ebv_tau353, header)
        write_data(data_av_tau353, header,
                   av_dir + 'perseus_av_planck_tau353.fits')

        # tau 353 error
        (data, header) = extract_data(datatype = 'tau353err')
        write_data(data, header, planck_dir + 'perseus_planck_tau353_error.fits')

        (data_ebv_tau353_error, header) = tau353_to_ebv(data, header)
        write_data(data_ebv_tau353, header,
                   planck_dir + 'perseus_ebv_error_planck_tau353.fits')

        (data_av_tau353_error, header) = ebv2av(data_ebv_tau353_error, header)
        write_data(data_av_tau353_error, header,
                   av_dir + 'perseus_av_error_planck_tau353.fits')
    if 0:
        # radiance
        (data, header) = extract_data(datatype = 'radiance')
        write_data(data, header, planck_dir + 'perseus_planck_radiance.fits')

        (data_ebv_radiance, header) = radiance_to_ebv(data, header)
        write_data(data_ebv_radiance, header,
                   planck_dir + 'perseus_ebv_planck_radiance.fits')

        (data_av_radiance, header) = ebv2av(data_ebv_radiance, header)
        write_data(data_av_radiance, header,
                   av_dir + 'perseus_av_planck_radiance.fits')

        # combine radiance and tau353 color excess maps
        data_av = combine_ebv(data_ebv_radiance, data_ebv_tau353)
        write_data(data_av, header,
                   av_dir + 'perseus_av_planck.fits')

        write_data(data_av / 10.0, header,
                   av_dir + 'perseus_av_error_planck.fits')

        # Band maps
        # ---------
    if 0:
        (data, header) = extract_data(datatype = '857')
        write_data(data, header, planck_dir + 'perseus_planck_857ghz.fits')

        (data, header) = extract_data(datatype = '545')
        write_data(data, header, planck_dir + 'perseus_planck_545ghz.fits')

        (data, header) = extract_data(datatype = '353')
        write_data(data, header, planck_dir + 'perseus_planck_353ghz.fits')

    if 0:
        # CO maps
        # -----------------
        (data, header) = extract_data(datatype = 'co_1to0')
        write_data(data, header, co_dir + 'perseus_co_1-0_planck.fits')

        (data, header) = extract_data(datatype = 'co_1to0_error')
        write_data(data, header, co_dir + 'perseus_co_1-0_error_planck.fits')

        (data, header) = extract_data(datatype = 'co_2to1')
        write_data(data, header, co_dir + 'perseus_co_2-1_planck.fits')

        (data, header) = extract_data(datatype = 'co_2to1_error')
        write_data(data, header, co_dir + 'perseus_co_2-1_error_planck.fits')

        (data, header) = extract_data(datatype = 'co_3to2')
        write_data(data, header, co_dir + 'perseus_co_3-2_planck.fits')

        (data, header) = extract_data(datatype = 'co_3to2_error')
        write_data(data, header, co_dir + 'perseus_co_3-2_error_planck.fits')

        (data, header) = extract_data(datatype = 'co')
        write_data(data, header, co_dir + 'perseus_co_planck.fits')

        (data, header) = extract_data(datatype = 'co_error')
        write_data(data, header, co_dir + 'perseus_co_error_planck.fits')


if __name__ == '__main__':
	main()


