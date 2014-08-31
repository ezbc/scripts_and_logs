#!/usr/bin/csh

cd /d/bip3/ezbc/california/data/galfa/

puthd in=california_galfa_cube.mir/restfreq value=1.420405752E+09 type='double'

# beamsize of cube: 4' x 3.5'
# geometric mean = (4*3.5)**0.5 = 3.7'
# degrees: 

regrid in=california_galfa_cube.mir \
    out=california_galfa_cube_bin_3.7arcmin.mir \
    desc=73,0,-0.0616,300,22,0,0.0616,265,0,100,1,200

fits in=california_galfa_cube_bin_3.7arcmin.mir \
    out=california_galfa_cube_bin_3.7arcmin.fits \
    op=xyout




