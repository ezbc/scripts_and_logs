#!/usr/bin/csh

# Miriad script for regridding the cfa cube onto the GALFA cube

cd /d/bip3/ezbc/perseus/data/galfa

set galfa_dir=../galfa
set cfa_dir=../cfa

fits in=/d/bip3/ezbc/taurus/data/cfa/taurus_cfa_cube_eqcoord.fits\
    out=$cfa_dir/perseus_cfa_cube_eqcoord.mir \
    op=xyin

fits in=$galfa_dir/perseus_galfa_cube_bin_4arcmin.fits\
    out=$galfa_dir/perseus_galfa_cube_bin_4arcmin.mir \
    op=xyin

regrid in=$cfa_dir/perseus_cfa_cube_eqcoord.mir \
    out=$cfa_dir/perseus_cfa_cube_galfa_regrid.mir \
    tin=$galfa_dir/perseus_galfa_cube_bin_3.7arcmin.mir \
    axes=1,2,3

fits in=$cfa_dir/perseus_cfa_cube_galfa_regrid.mir \
    out=$cfa_dir/perseus_cfa_cube_galfa_regrid.fits \
    op=xyout

