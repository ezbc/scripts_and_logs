#!/usr/bin/python

from astropy.io import fits
import numpy as np

def calc_neg_fraction(rh2):

    rh2 = rh2[~np.isnan(rh2)]

    neg_frac = float(np.nansum(rh2 < 0)) / np.size(rh2)

    return neg_frac

def main():

    DIR_DATA = '/d/bip2/lee/perseus_cloud/paper1/'
    FILENAME_RH2 = DIR_DATA + 'R_H2.fits'

    # load the data
    rh2 = fits.getdata(FILENAME_RH2)

    # calculate fraction of negative pixels
    neg_frac = calc_neg_fraction(rh2)

    print('Fraction of negative RH2 pixels: {0:.2f}'.format(neg_frac))

if __name__ == '__main__':
    main()
