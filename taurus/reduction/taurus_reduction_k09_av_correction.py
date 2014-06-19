#!/usr/bin/python

def main():
    from astropy.io import fits
    import numpy as np

    av_dir = '/d/bip3/ezbc/taurus/data/av/'

    data, header = fits.getdata(av_dir + 'taurus_av_kainulainen2009.fits',
                                header=True)

    data[data == -1] = np.NaN

    fits.writeto(av_dir + 'taurus_av_kainulainen2009_nan.fits',
                 data,
                 header=header)

if __name__ == '__main__':
	main()
