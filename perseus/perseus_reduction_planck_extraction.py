#!/usr/bin/python

import planckpy as pl
reload(pl)
from astropy.io import fits as pf


def extract_data(field = 2):
    data_location = '/d/bip3/ezbc/planck/planck_raw_data/'

    # Goal: RA from 1:30 to 5:30hrs, Dec = 20 to 38deg

    dec_range = (25, 35)
    ra_range = (40, 65)

    (data, header) = pl.get_data(data_location = data_location,
            data_type = 'Dust Opacity',
            x_range = ra_range,
            y_range = dec_range,
            coord_type = 'equatorial',
            field = field,
            dr_version = 2,
            resolution = 0.01,
            cut_last_pixel = False,
            verbose = True)

    return (data, header)

def write_data(data, header, filename):

    pf.writeto(filename,
            data,
            header = header,
            clobber = True,)
            #output_verify = 'fix')

def ebv2av(data, header):
    # Av=E(B-V) x 3.05263, 152 in Tielens' ISM book.

    data *= 3.05363
    header['BUNIT'] = 'mag'

    return (data, header)

def main():
    av_dir = '/d/bip3/ezbc/perseus/data/av/'
    planck_dir = '/d/bip3/ezbc/perseus/data/planck/'

    (data, header) = extract_data()
    write_data(data, header, planck_dir + 'perseus_planck_ebv.fits')

    ebv2av(data, header)
    write_data(data, header, av_dir + 'perseus_planck_av.fits')

    (data, header) = extract_data(field = 3)
    write_data(data, header, planck_dir + 'perseus_planck_ebv_error.fits')

    ebv2av(data, header)
    write_data(data, header, av_dir + 'perseus_planck_av_error.fits')

if __name__ == '__main__':
	main()


