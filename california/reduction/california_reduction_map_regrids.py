#!/usr/bin/python

# Planck resolution is 5'. We will smooth and regrid all maps to the Planck map.

# Map                           Resolution      Convol size
# GALFA HI                      3.7'            3.363
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

    os.chdir('/d/bip3/ezbc/california/data')

    # If true, deletes files to be written
    clobber = 0

    in_images = (
              'hi/california_hi_galfa_cube',
              'av/california_av_planck_tau353',
              'av/california_av_error_planck_tau353',
              'dust_temp/california_dust_temp',
              'dust_temp/california_dust_temp_error',
              '/d/bip3/ezbc/multicloud/data/av/multicloud_av_k09_nan',
              'co/california_co_cfa_cube')

    im_pl = 'av/california_av_planck_tau353'
    im_pl_err = 'av/california_av_error_planck_tau353'
    im_Td = 'dust_temp/california_dust_temp'
    im_Td_err = 'dust_temp/california_dust_temp_error'
    im_k09 = 'av/california_av_k09'
    im_hi = 'hi/california_hi_galfa_cube'
    im_co = 'co/california_co_cfa_cube'

    # Load the images into miriad
    print('\nLoading images into miriad')
    out_images = (im_hi, im_pl, im_pl_err,
                  im_Td,
                  im_Td_err,
            im_k09, im_co)

    for i in xrange(len(in_images)):
        exists = check_file(out_images[i] + '.mir', clobber=clobber)
        print('\t{:s}.fits\n'.format(in_images[i]))
        if not exists:
            fits(in_images[i] + '.fits',
                    out=out_images[i] + '.mir',
                    op='xyin')

    # Regrid Planck images and HI image to have one beam/pixel
    print('\nRegridding Planck images')

    images = (im_pl, im_pl_err,
                  im_Td,
                  im_Td_err,
            im_hi)
    delta_ra = -0.083333333
    delta_dec = 0.083333333

    # Greater RA value comes first
    ref_pix, npix = calc_image_origin(x_limits=(81, 60),
                                      y_limits=(25, 39),
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

        print('\t{:s}_5arcmin.mir\n'.format(image))

        if not exists:
            regrid(image + '.mir',
                    out=image + '_5arcmin.mir',
                    desc=desc)

    # Determine convolution beam sizes to smooth to Planck res
    print('\nSmoothing images to Planck resolution')

    planck_beam = 5.0 # arcsec
    im_beams = np.array([3.7, 2.7]) # arcsec
    conv_beams = (planck_beam**2 - im_beams**2)**0.5

    images = [im_hi, im_k09]

    for i in xrange(len(images)):
        exists = check_file(images[i] + \
                 '_smooth_planckres.mir', clobber=clobber)

        print('\t{:s}.mir\n'.format(images[i]))

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

    images = [im_pl + '_5arcmin',
              im_pl_err + '_5arcmin',
              im_Td + '_5arcmin',
              im_Td_err + '_5arcmin',
              im_k09 + '_regrid_planckres',
              im_co + '_regrid_planckres',
              im_hi + '_regrid_planckres']

    for image in images:
        exists = check_file(image + '.fits', clobber=clobber)

        print('\t{:s}.mir\n'.format(image))

        fits(image + '.mir',
             out=image + '.fits',
             op='xyout')

if __name__ == '__main__':
    main()
