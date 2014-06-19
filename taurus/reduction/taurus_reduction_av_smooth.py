#!/usr/bin/python

''' Script to smooth and regrid av images onto GALFA cube. Must start ipython
in c-shell for miriad wrapper to work.
'''

def __id_generator(size=6,):
    ''' Creates random string of letters and numbers. Useful for creating
    temporary files.'''

    # import external modules
    import string
    import random

    chars=string.ascii_uppercase + string.digits

    return ''.join(random.choice(chars) for x in range(size))

def smooth_image(input_file=None, output_file=None, final_res=None,
        clobber=True, verbose=True):
    ''' Smooths an image to a desired resolution in arcsec.
    '''

    # import external modules
    from scipy.ndimage.filters import gaussian_filter as smooth
    import pyfits as pf

    # Load fits file
    image, header = pf.getdata(input_file, header=True)

    # Smooth the image
    # convert from degrees/pix to "/pix
    fwhm = final_res / (header['CDELT2'] * 3600.)
    image_smooth = smooth(image,sigma=(fwhm/2.,fwhm/2.))

    # Write out the fits file
    pf.writeto(output_file, image_smooth, header=header, clobber=clobber)

def regrid_image(input_file=None, output_file=None, template_file=None,
        clobber=True, verbose=True):
    ''' Grids the input_file onto the template_file and writes output to
    output_file. Uses the miriad wrapper task regrid.
    '''

    # import external modules
    from mirexec import TaskRegrid as regrid
    from mirexec import TaskFits as fits
    from os import system,path

    # Check if output_file exists
    execute = True
    if path.lexists(input_file):
        if clobber:
            if verbose:
                print('Overwriting file ' + output_file)
            system('rm -rf ' + output_file)
        elif not clobber:
            if verbose:
                print('File exists, cannot proceed ' + output_file)
            execute = False

    # If the file doesn't exist
    if execute:
        # generate id
        input_file_temp1,input_file_temp2 = __id_generator(), __id_generator()
        template_file_temp = __id_generator()

        # load input file
        mir_task = fits(in_=input_file,
                out=input_file_temp1,
                op='xyin')
        mir_task.run()

        # load template file
        mir_task = fits(in_=template_file,
                out=template_file_temp,
                op='xyin')
        mir_task.run()

        # regrid file
        mir_task = regrid(in_=input_file_temp1,
                tin=template_file_temp,
                out=input_file_temp2,
                axes=(1,2))
        mir_task.run()

        # write regridded file
        mir_task = fits(in_=input_file_temp2,
                out=output_file,
                op='xyout')
        mir_task.run()
        system('rm -rf ' + input_file_temp1 + ' ' + \
                           input_file_temp2 + ' ' + \
                           template_file_temp)

def main():
    ''' Executes the script.
    '''

    # import external modules
    import numpy as np
    import pyfits as pf

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'
    co_dir = '/d/bip3/ezbc/taurus/data/cfa/'

    smooth_image(input_file = av_dir + 'taurus_av_kainulainen2009.fits',
                 output_file = av_dir + 'taurus_av_k09_smooth.fits',
                 final_res=240.) #240", 4'

    regrid_image(av_dir + 'taurus_av_k09_smooth.fits',
        output_file = av_dir + 'taurus_av_k09_regrid.fits',
        template_file = hi_dir + 'taurus.galfa.cube.bin.4arcmin.fits')

    # Now perform for Pineda image
    # Pineda et al. 2010, ApJ, 721, 686 
    # resolution of Av image is 200", see caption of figure 4
    # thus we will convolve the Av image to the Gaussian filer size

    # check size of pineda image, originally had an extra axis
    av_image,av_header = pf.getdata(av_dir + 'taurus_av_pineda2010.fits',
            header=True)
    if len(av_image.shape) == 3:
        av_image = av_image[0,:,:]
        av_header['NAXIS'] = 2
        del av_header['NAXIS3']
        pf.writeto(av_dir + 'taurus_av_pineda2010.fits',av_image,
                header=av_header,
                clobber=True)

    smooth_image(input_file = av_dir + 'taurus_av_pineda2010.fits',
                 output_file = av_dir + 'taurus_av_p10_smooth.fits',
                 final_res=240.,) #240", 4'

    regrid_image(input_file = av_dir + 'taurus_av_p10_smooth.fits',
        output_file = av_dir + 'taurus_av_p10_regrid.fits',
        template_file = hi_dir + 'taurus.galfa.cube.bin.4arcmin.fits')

if __name__ == '__main__':
    main()




