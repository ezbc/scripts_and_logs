#!/usr/bin/python


import os

def make_model_dir(cubes):

    #model_dir = '/d/bip3/ezbc/shield/749237_lowres/modeling_cmode1_3inc/'
    model_dir = '/d/bip3/ezbc/shield/749237_lowres/modeling_highres'

    if 1:
        try:
            os.system('mkdir ' + model_dir)
        except:
           pass

    os.system('cp create_cola_script.pro ' + model_dir)
    #os.system('cp gipsy_script.col ' + model_dir)
    for cube in cubes:
        os.system('cp  ' + cube + ' ' + model_dir)
        os.system('cp  ' + cube.replace('cube','nhi') + ' ' + model_dir)

def main():

    os.chdir('/d/cosmos/ezbc/research/scripts/shield/')
    cube_dir = '/d/bip3/ezbc/shield/749237_lowres/'
    cubes = (
             cube_dir + '749237_rebin_cube_regrid.fits',
             cube_dir + '749237_rebin_cube_error_regrid.fits',
             cube_dir + '749237_rebin_cube.fits',
             )

    make_model_dir(cubes)

if __name__ == '__main__':
    main()



