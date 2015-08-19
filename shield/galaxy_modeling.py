#!/usr/bin/python


import os

def make_model_dir(cubes):

    model_dir = '/d/bip3/ezbc/shield/modeling/'

    if 0:
        try:
            os.system('mkdir ' + model_dir)
        except:
            pass

    os.system('cp create_cola_script.pro ' + model_dir)
    os.system('cp gipsy_script.col ' + model_dir)
    for cube in cubes:
        os.system('cp  ' + cube + ' ' + model_dir)
        os.system('cp  ' + cube.replace('cube','nhi') + ' ' + model_dir)

def main():

    cube_dir = '/d/bip3/ezbc/shield/'
    cubes = (
             cube_dir + '749237_cube_regrid.fits',
             cube_dir + '749237_cube_error_regrid.fits',
             cube_dir + '749237_cube.fits',
             )

    make_model_dir(cubes)

if __name__ == '__main__':
    main()



