#!/usr/bin/python

# Planck resolution is 5'. We will smooth and regrid all maps to the Planck map.

# Map                           Resolution      Convol size
# Lee et al. (2012) Av          5'              5
# GALFA HI                      3.7'            3.363
# CfA CO Dame et al. (2001)     8'              No smoothing

import os
import mirpy
import numpy as np

def check_file(filename, clobber=False):

    exists = False

    if os.path.isfile(filename) or os.path.isdir(filename):
        exists = True
        print('\tImage {:s} exists'.format(filename))
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

    from mirpy import fits, regrid, smooth, imbin
    from mycoords import calc_image_origin
    from astropy.io import fits as pyfits
    import numpy as np

    os.chdir('/d/bip3/ezbc/multicloud/data')

    # If true, deletes files to be written
    clobber = False
    clobber_hi = False

    # First, change zeros in lee image to nans
    data, header = pyfits.getdata('av/multicloud_av_lee12_2mass.fits', header=True)
    data[data == 0] = np.nan
    pyfits.writeto('av/multicloud_av_lee12_2mass_nan.fits',
                   data,
                   header,
                   clobber=True)
    data, header = pyfits.getdata('av/multicloud_av_k09.fits',
                                  header=True)
    data[data == -1] = np.nan
    pyfits.writeto('av/multicloud_av_k09_nan.fits',
                   data,
                   header,
                   clobber=True)

    in_images = ('hi/multicloud_hi_galfa_cube',
              'av/multicloud_av_planck',
              'av/multicloud_av_error_planck',
              'dust_temp/multicloud_dust_temp',
              'dust_temp/multicloud_dust_temp_error',
              'dust_temp/multicloud_dust_beta',
              'dust_temp/multicloud_dust_beta_error',
              'av/multicloud_av_planck_radiance',
              'av/multicloud_av_error_planck_radiance',
              'co/multicloud_co_cfa_cube',
              'av/multicloud_av_k09_nan',
              'av/multicloud_av_lee12_2mass_nan',
              'av/multicloud_av_lee12_iris_masked')

    im_hi = 'hi/multicloud_hi_galfa_cube'
    im_pl = 'av/multicloud_av_planck'
    im_pl_err = 'av/multicloud_av_error_planck'
    im_Td = 'dust_temp/multicloud_dust_temp'
    im_Td_err = 'dust_temp/multicloud_dust_temp_error'
    im_beta = 'dust_temp/multicloud_dust_beta'
    im_beta_err = 'dust_temp/multicloud_dust_beta_error'
    im_pl2 = 'av/multicloud_av_planck_radiance'
    im_pl2_err = 'av/multicloud_av_error_planck_radiance'
    im_co = 'co/multicloud_co_cfa_cube'
    im_k09 = 'av/multicloud_av_k09_nan'
    im_lee_2mass = 'av/multicloud_av_lee12_2mass'
    im_lee_iris = 'av/multicloud_av_lee12_iris'

    # Load the images into miriad
    print('\nLoading images into miriad')
    out_images = (im_hi,
                  im_pl,
                  im_pl_err,
                  im_Td,
                  im_Td_err,
                  im_beta,
                  im_beta_err,
                  im_pl2,
                  im_pl2_err,
                  im_co,
                  im_k09,
                  im_lee_2mass,
                  im_lee_iris)

    for i in xrange(len(in_images)):
        if in_images[i] == im_hi:
            exists = check_file(out_images[i] + '.mir', clobber=clobber_hi)
        else:
            exists = check_file(out_images[i] + '.mir', clobber=clobber)

        print('\tLoading {:s}.fits\n'.format(in_images[i]))
        if not exists:
            fits(in_images[i] + '.fits',
                 out=out_images[i] + '.mir',
                 op='xyin')

    # Regrid Planck images and HI image to have one beam/pixel
    print('\nRegridding Planck images')

    images = (im_pl, im_pl_err, im_Td, im_Td_err, im_beta, im_beta_err, im_pl2, im_pl2_err, im_hi)

    desc = (59.75,0,-0.08333,180,26.05,0,0.08333,132)

    delta_ra = -0.083333333
    delta_dec = 0.083333333

    ref_pix, npix = calc_image_origin(x_limits=(85, 32),
                                      y_limits=(15, 48),
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
            exists = check_file(image + '_5arcmin.mir', clobber=clobber_hi)
        else:
            desc = desc_av
            exists = check_file(image + '_5arcmin.mir', clobber=clobber)

        if not exists:
            print('\tRegridding {:s}_5arcmin.mir\n'.format(image))
            regrid(image + '.mir',
                    out=image + '_5arcmin.mir',
                    desc=desc)

    # Determine convolution beam sizes to smooth to Planck res
    print('\nSmoothing images to Planck resolution')

    planck_beam = 5.0 # arcmin
    im_beams = np.array([3.7, 2.5, 2.5, 4.3]) # arcmin
    conv_beams = (planck_beam**2 - im_beams**2)**0.5

    images = [im_hi, im_k09, im_lee_2mass, im_lee_iris]

    for i, image in enumerate(images):
        if image == im_hi:
            desc = desc_hi
            exists = check_file(image + '_smooth_planckres.mir',
                                clobber=clobber_hi)
        else:
            desc = desc_av
            exists = check_file(image + '_smooth_planckres.mir',
                                clobber=clobber)
        print image

        if not exists:
            if images[i] == im_hi:
                image = images[i] + '_5arcmin'
            else:
                image = images[i]

            print('\tSmoothing {:s}.mir\n'.format(image))

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

        if not exists:
            print('\tRegridding {:s}_smooth_planckres.mir\n'.format(image))
            if image == im_co or image == im_lee_2mass or image == im_lee_iris:
                regrid(image + '.mir',
                       out=image + '_regrid_planckres.mir',
                       tin=im_pl + '_5arcmin.mir',
                       axes=(1,2))
            else:
                regrid(image + '_smooth_planckres.mir',
                       out=image + '_regrid_planckres.mir',
                       tin=im_pl + '_5arcmin.mir',
                       axes=(1,2))

    # Bin the images
    print('\nBin the images')

    images = [
              im_pl + '_5arcmin',
              im_pl_err + '_5arcmin',
              im_Td + '_5arcmin',
              im_Td_err + '_5arcmin',
              im_beta + '_5arcmin',
              im_beta_err + '_5arcmin',
              im_pl2 + '_5arcmin',
              im_pl2_err + '_5arcmin',
              im_k09 + '_regrid_planckres',
              im_lee_2mass + '_regrid_planckres',
              im_lee_iris + '_regrid_planckres',
              im_co + '_regrid_planckres',
              im_hi + '_regrid_planckres']

    for image in images:
        exists = check_file(image + '_bin.mir', clobber=clobber)

        if not exists:
            print('\tWriting {:s}.mir\n'.format(image))

            bin_size = 10
            if image == im_hi + '_regrid_planckres' or \
               image == im_co + '_regrid_planckres':
                bins = bin_size, bin_size, bin_size, bin_size,1,1
            else:
                bins = bin_size, bin_size, bin_size, bin_size,

            imbin(image + '.mir',
                  out=image + '_bin.mir',
                  bin=bins,
                  )

    # Write the images out to fits images
    print('\nWriting images to fits format')

    for image in images:
        if image == im_hi + '_regrid_planckres':
            exists = check_file(image + '.fits', clobber=clobber_hi)
        else:
            exists = check_file(image + '.fits', clobber=clobber)

        if not exists:
            print('\tWriting {:s}.fits\n'.format(image))

            fits(image + '.mir',
                 out=image + '.fits',
                 op='xyout')

        exists = check_file(image + '_bin.fits', clobber=clobber)

        if not exists:
            print('\tWriting {:s}_bin.fits\n'.format(image))

            fits(image + '_bin.mir',
                 out=image + '_bin.fits',
                 op='xyout')

if __name__ == '__main__':
    main()
