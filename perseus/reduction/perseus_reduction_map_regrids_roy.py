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

    from mirpy import fits, regrid, smooth
    from mycoords import calc_image_origin
    from astropy.io import fits as pyfits
    import numpy as np

    os.chdir('/d/bip3/ezbc/perseus/data')

    # If true, deletes files to be written
    clobber = 1
    clobber_hi = 1

    in_images = ('/d/bip3/ezbc/multicloud/data/hi/multicloud_hi_galfa_cube',
                 )

    im_hi = 'hi/perseus_hi_galfa_cube_roy'

    # Load the images into miriad
    print('\nLoading images into miriad...')
    out_images = (im_hi,
                  )

    for i in xrange(len(in_images)):
        if out_images[i] == im_hi:
            clobber_temp = clobber_hi
        else:
            clobber_temp = clobber
        exists = check_file(out_images[i] + '.mir', clobber=clobber_temp)
        if not exists:
            print('\tLoading {:s}.fits\n'.format(in_images[i]))
            fits(in_images[i] + '.fits',
                    out=out_images[i] + '.mir',
                    op='xyin')

    # Regrid Planck images and HI image to have one beam/pixel
    print('\nRegridding images')

    images = (im_hi,)

    # The center coordinate of the region is RA: 4h13m26s, Dec =33d22.
    # Will it be too much of asking if we want to get a coverage of around 20deg
    # X20 deg?

    # In terms of rectangular region it will be

    # DEC ---> 20deg to 40 deg
    # RA-------> 4h40m to 3h20m

    #desc = (59.75,0,-0.08333,180,26.05,0,0.08333,132)

    delta_ra = -0.083333333 / 5.0
    delta_dec = 0.083333333 / 5.0

    ref_pix, npix = calc_image_origin(x_limits=(15 * (4 + 40./60.),
                                                15 * (3 + 20./60.)),
                                      y_limits=(20, 38),
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
        if image in (im_hi,):
            desc = desc_hi
        else:
            desc = desc_av

        exists = check_file(image + '_regrid.mir', clobber=clobber)


        if not exists:
            print('\tRegridding {:s}_regrid.mir\n'.format(image))
            regrid(image + '.mir',
                    out=image + '_regrid.mir',
                    desc=desc)

    images = [im_hi,]

    # Write the images out to fits images
    print('\nWriting images to fits format')

    images = [
              im_hi + '_regrid',
              ]

    for image in images:
        exists = check_file(image + '.fits', clobber=clobber)

        if not exists:
            print('\tWriting {:s}.mir\n'.format(image))

            fits(image + '.mir',
                 out=image + '.fits',
                 op='xyout')

if __name__ == '__main__':
    main()
