#!/usr/bin/csh

# Miriad script for regridding the cfa cube onto the GALFA cube

cd /d/bip3/ezbc/taurus/data/galfa

set galfa_dir=../galfa
set cfa_dir=../cfa

fits in=$cfa_dir/taurus_cfa_cube_eqcoord.mir\
    out=$cfa_dir/taurus_cfa_cube_eqcoord.fits \
    op=xyout

fits in=$galfa_dir/taurus_galfa_cube_bin_3.7arcmin.fits\
    out=$galfa_dir/taurus_galfa_cube_bin_3.7arcmin.mir \
    op=xyin

regrid in=$cfa_dir/taurus_cfa_cube_eqcoord.mir \
    out=$cfa_dir/taurus_cfa_cube_galfa_regrid.mir \
    tin=$galfa_dir/taurus_galfa_cube_bin_3.7arcmin.mir \
    axes=1,2,3

fits in=$cfa_dir/taurus_cfa_cube_galfa_regrid.mir \
    out=$cfa_dir/taurus_cfa_cube_galfa_regrid.fits \
    op=xyout


