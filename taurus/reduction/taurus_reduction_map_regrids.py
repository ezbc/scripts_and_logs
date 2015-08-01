#!/usr/bin/python

# Planck resolution is 5'. We will smooth and regrid all maps to the Planck map.

# Map                           Resolution      Convol size
# GALFA HI                      3.7'            3.363
# Kainulainen et al. (2009) Av  2.4'            4.386
# Pineda et al. (2010) Av       3.3'            3.75
# CfA CO Dame et al. (2001)     8'              No smoothing

import os
import mirpy
import numpy as np

def check_file(filename, clobber=False):

    exists = False

    if os.path.isfile(filename) or os.path.isdir(filename):
        exists = True
        if clobber:
            print('\tDeleting image {:s}'.format(filename))
            os.system('rm -rf {:s}'.format(filename))
            exists = False

    return exists

def regrid_images(image, template):
    import mirpy

    # Avoiding the miriad errors avoids errors from existing files
    try:
        mirpy.fits(image+'.fits',
                out=image+'.mir',
                op='xyin')
    except mirpy.wrapper.MiriadError:
        pass

    try:
        mirpy.fits(template+'.fits',
                out=template+'.mir',
                op='xyin')
    except mirpy.wrapper.MiriadError:
        pass

    try:
        mirpy.regrid(image+'.mir', tin=template+'.mir',
                out=image+'_regrid.mir')
    except mirpy.wrapper.MiriadError:
        pass

    try:
        mirpy.fits(image+'_regrid.mir', out=image+'_regrid.fits',
                op='xyout')
    except mirpy.wrapper.MiriadError:
        pass

def main():

    ''' Planck resolution is 5'. We will smooth and regrid all maps to the
    Planck map.

    Map                           Resolution      Convol size
    GALFA HI                      3.7'            3.363
    Kainulainen et al. (2009) Av  2.4'            4.386
    Pineda et al. (2010) Av       3.3'            3.75
    CfA CO Dame et al. (2001)     8'              No smoothing
    '''

    from mirpy import fits, regrid, smooth
    from mycoords import calc_image_origin

    os.chdir('/d/bip3/ezbc/taurus/data')

    # If true, deletes files to be written
    clobber = 1

    in_images = (
              #'av/taurus_av_kainulainen2009_nan',
              '/d/bip3/ezbc/multicloud/data/av/multicloud_av_k09_nan',
              'av/taurus_av_pineda2010',
              'hi/taurus_hi_galfa_cube',
              'av/taurus_av_planck_tau353',
              'av/taurus_av_error_planck_tau353',
              'av/taurus_av_planck_radiance',
              'av/taurus_av_error_planck_radiance',
              'co/taurus_co_cfa_cube',
              'co/taurus_co_planck',
              'co/taurus_co_error_planck',
              'co/taurus_co_1-0_planck',
              'co/taurus_co_1-0_error_planck',
              'co/taurus_co_2-1_planck',
              'co/taurus_co_2-1_error_planck',
              'co/taurus_co_3-2_planck',
              'co/taurus_co_3-2_error_planck',
              )

    im_pl = 'av/taurus_av_planck_tau353'
    im_pl_err = 'av/taurus_av_error_planck_tau353'
    im_pl2 = 'av/taurus_av_planck_radiance'
    im_pl2_err = 'av/taurus_av_error_planck_radiance'
    im_hi = 'hi/taurus_hi_galfa_cube'
    im_co = 'co/taurus_co_cfa_cube'
    im_p10 = 'av/taurus_av_p10'
    im_k09 = 'av/taurus_av_k09'
    im_pl_co = 'co/taurus_co_planck'
    im_pl_co_err = 'co/taurus_co_error_planck'
    im_pl_co_10 = 'co/taurus_co_1-0_planck'
    im_pl_co_10_err = 'co/taurus_co_1-0_error_planck'
    im_pl_co_21 = 'co/taurus_co_2-1_planck'
    im_pl_co_21_err = 'co/taurus_co_2-1_error_planck'
    im_pl_co_32 = 'co/taurus_co_3-2_planck'
    im_pl_co_32_err = 'co/taurus_co_3-2_error_planck'

    # Load the images into miriad
    print('\nLoading images into miriad')
    out_images = (im_k09,
                  im_p10,
                  im_hi,
                  im_pl,
                  im_pl_err,
                  im_pl2,
                  im_pl2_err,
                  im_co,
                  im_pl_co,
                  im_pl_co_err,
                  im_pl_co_10,
                  im_pl_co_10_err,
                  im_pl_co_21,
                  im_pl_co_21_err,
                  im_pl_co_32,
                  im_pl_co_32_err,
                  )

    for i in xrange(len(in_images)):
        exists = check_file(out_images[i] + '.mir', clobber=clobber)
        print('\tReading {:s}.fits\n'.format(in_images[i]))
        if not exists:
            fits(in_images[i] + '.fits',
                    out=out_images[i] + '.mir',
                    op='xyin')

            print('\tout = ' + out_images[i] + '.mir')

    # Regrid Planck images and HI image to have one beam/pixel
    print('\nRegridding Planck images')

    images = (
              im_pl,
              im_pl_err,
              im_pl2,
              im_pl2_err,
              im_hi,
              #im_pl_co,
              #im_pl_co_err,
              #im_pl_co_10,
              #im_pl_co_10_err,
              #im_pl_co_21,
              #im_pl_co_21_err,
              #im_pl_co_32,
              #im_pl_co_32_err,
              )

    delta_ra = -0.083333333
    delta_dec = 0.083333333
    lim_ra_wcs, lim_dec_wcs = 77.5, 18.0
    ref_ra_wcs, ref_dec_wcs = 0.0, 0.0

    ref_pix, npix = calc_image_origin(x_limits=(83, 53.0),
                                      y_limits=(15.0, 35.0),
                                      delta_x=delta_ra,
                                      delta_y=delta_dec)

    desc_av = [0, ref_pix[0], delta_ra, npix[0], \
               0, ref_pix[1], delta_dec, npix[1]]

    low_vel = -100.0
    high_vel = 100.0
    vel_res = 0.16667
    vel_npix = int((high_vel - low_vel) / vel_res)
    ref_pix_vel = int(vel_npix / 2.0) * vel_res
    ref_pix_vel = vel_npix / 2.0

    desc_hi = [0, ref_pix[0], delta_ra, npix[0], \
               0, ref_pix[1], delta_dec, npix[1], \
               0, ref_pix_vel, vel_res, vel_npix]

    for image in images:

        # If HI, regrid the velocity axis as well
        if image == im_hi:
            desc = desc_hi
        else:
        	desc = desc_av

        exists = check_file(image + '_5arcmin.mir', clobber=clobber)

        print('\tReading image {:s}.mir'.format(image))
        print('\tWriting image {:s}_5arcmin.mir\n'.format(image))

        if not exists:
            regrid(image + '.mir',
                    out=image + '_5arcmin.mir',
                    desc=desc)

    # Determine convolution beam sizes to smooth to Planck res
    print('\nSmoothing images to Planck resolution')

    planck_beam = 5.0 # arcsec
    im_beams = np.array([2.4, 3.3, 3.7]) # arcsec
    conv_beams = (planck_beam**2 - im_beams**2)**0.5

    images = [im_k09, im_p10, im_hi]

    for i in xrange(len(images)):

        exists = check_file(images[i] + '_smooth_planckres.mir',
                clobber=clobber)

        print('\t{:s}.mir\n'.format(image))

        if not exists:
            if images[i] == im_hi:
                image = images[i] + '_5arcmin'
            else:
            	image = images[i]

            smooth(image + '.mir',
                   out=images[i] + '_smooth_planckres.mir',
                   fwhm=conv_beams[i],
                   pa=0,
                   scale=0.0)

    # Regrid all images to Planck image
    print('\nRegridding images to Planck grid')

    images.append(im_co)

    for image in images:
        exists = check_file(image + '_regrid_planckres.mir', clobber=clobber)
        print('\t{:s}_smooth_planckres.mir\n'.format(image))

        if not exists:
            if image == im_co:
                regrid(image + '.mir',
                       out=image + '_regrid_planckres.mir',
                       tin=im_pl + '_5arcmin.mir',
                       axes=(1,2))
            else:
                regrid(image + '_smooth_planckres.mir',
                       out=image + '_regrid_planckres.mir',
                       tin=im_pl + '_5arcmin.mir',
                       axes=(1,2))

    # Write the images out to fits images
    print('\nWriting images to fits format')

    images = [
              im_pl + '_5arcmin',
              im_pl_err + '_5arcmin',
              im_pl2 + '_5arcmin',
              im_pl2_err + '_5arcmin',
              im_pl_co + '_5arcmin',
              im_pl_co_err + '_5arcmin',
              im_pl_co_10 + '_5arcmin',
              im_pl_co_10_err + '_5arcmin',
              im_pl_co_21+ '_5arcmin',
              im_pl_co_21_err+ '_5arcmin',
              im_pl_co_32+ '_5arcmin',
              im_pl_co_32_err+ '_5arcmin',
              im_k09 + '_regrid_planckres',
              im_p10 + '_regrid_planckres',
              im_co + '_regrid_planckres',
              im_hi + '_regrid_planckres']

    for image in images:
        check_file(image + '.fits', clobber=clobber)

        print('\t{:s}.mir\n'.format(image))

        fits(image + '.mir',
             out=image + '.fits',
             op='xyout')

if __name__ == '__main__':
    main()
