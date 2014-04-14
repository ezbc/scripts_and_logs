#!/bin/csh

# Regrids the planck Av image onto the 3.7" pixel GALFA cube.

cd /d/bip3/ezbc/california/data/galfa

set planck_dir=../av
set galfa_dir=../galfa

fits in=$planck_dir/california_planck_av.fits \
        out=$planck_dir/california_planck_av.mir \
        op=xyin

regrid in=$planck_dir/california_planck_av.mir \
    tin=$galfa_dir/california_galfa_cube_bin_3.7arcmin.mir \
    out=$planck_dir/california_planck_av_regrid.mir \
    axes=1,2

fits in=$planck_dir/california_planck_av_regrid.mir \
        out=$planck_dir/california_planck_av_regrid.fits \
        op=xyout


