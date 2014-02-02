#!/usr/bin/python

''' Script to smooth and regrid av, cfa, and hi images to 36' resolution to
compare with Paradis et al. (2012). Must start ipython in c-shell for miriad
wrapper to work.  '''

def __id_generator(size=10,):
    ''' Creates random string of letters and numbers. Useful for creating
    temporary files.'''

    # import external modules
    import string
    import random

    chars=string.ascii_uppercase + string.digits

    return ''.join(random.choice(chars) for x in range(size))

def bin_image(input_file=None, output_file=None, bin_input=None,
        clobber=True, verbose=True):
    ''' Grids the input_file onto the template_file and writes output to
    output_file. Uses the miriad wrapper task regrid.
    '''

    # import external modules
    from mirexec import TaskRegrid as regrid
    from mymirexec import TaskImBin as imbin
    from mirexec import TaskFits as fits
    from os import system,path

    # Check if output_file exists
    execute = True
    if path.lexists(input_file):
        if clobber and input_file != output_file:
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

        # load input file
        mir_task = fits(in_=input_file,
                out=input_file_temp1,
                op='xyin')
        mir_task.run()

        # regrid file
        mir_task = imbin(in_=input_file_temp1,
                out=input_file_temp2,
                bin=bin_input)
        mir_task.run()

        # write regridded file
        mir_task = fits(in_=input_file_temp2,
                out=output_file,
                op='xyout')
        mir_task.run()

        system('rm -rf ' + input_file_temp1 + ' ' + \
                           input_file_temp2)

def smooth_image(input_file=None, output_file=None, final_res=None,
        clobber=True, verbose=True):
    ''' Smooths an image to a desired resolution in arcsec.
    '''

    # import external modules
    from scipy.ndimage.filters import gaussian_filter as smooth
    import pyfits as pf
    import numpy as np

    # Load fits file
    image, header = pf.getdata(input_file, header=True)

    print input_file,image.shape

    # Smooth the image
    # convert from degrees/pix to "/pix
    fwhm = final_res / (header['CDELT2'] * 3600.)
    if image.ndim > 2:
        image_smooth = np.zeros(image.shape)
        for i in range(image.shape[0]):
            image_smooth[i,:,:] = smooth(image[i,:,:],sigma=(fwhm/2.,fwhm/2.))
    else:
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

    # Define resolution used by Paradis et al. (2009)
    final_res = 36. * 60. # arcmin

    # Smooth Av image
    smooth_image(input_file = av_dir + 'taurus_av_kainulainen2009.fits',
                 output_file = av_dir + 'taurus_av_k09_paradis_smooth.fits',
                 final_res=final_res)
    # Smooth GALFA HI image
    smooth_image(input_file = hi_dir + 'taurus.galfa.cube.bin.4arcmin.fits',
                 output_file = hi_dir + 'taurus_galfa_paradis_smooth.fits',
                 final_res=final_res)
    # Smooth CfA CO image
    smooth_image(input_file = co_dir + 'taurus.cfa.galfaBin.4arcmin.fits',
                 output_file = co_dir + 'taurus_cfa_paradis_smooth.fits',
                 final_res=final_res)

    # 1 arcmin pixels bin 13 pixels to get 13.7' pixels
    # performed in command line
    #bin_image(input_file = av_dir + 'taurus_av_k09_paradis_smooth.fits',
    #         output_file = av_dir + 'taurus_av_k09_paradis_bin.fits',
    #         bin_input = (13,13,13,13))

    # Regrid GALFA HI data onto binned Av image
    regrid_image(hi_dir + 'taurus_galfa_paradis_smooth.fits',
        output_file = hi_dir + 'taurus_galfa_paradis_bin.fits',
        template_file = av_dir + 'taurus_av_k09_paradis_bin.fits')

    # Regrid CfA CO data onto binned Av image
    regrid_image(co_dir + 'taurus_cfa_paradis_smooth.fits',
        output_file = co_dir + 'taurus_cfa_paradis_bin.fits',
        template_file = av_dir + 'taurus_av_k09_paradis_bin.fits')

if __name__ == '__main__':
    main()




