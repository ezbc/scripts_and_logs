#!/usr/bin/csh

cd /d/bip3/ezbc/perseus/data

# Load fits files
fits in=av/perseus_av_lee12_eq.fits \
    out=av/perseus_av_lee12_eq.mir \
    op=xyin

fits in=hi/perseus_hi_galfa_cube.fits \
    out=hi/perseus_hi_galfa_cube.mir \
    op=xyin

# regrid galfa hi image to have lower velocity resolution
regrid in=hi/perseus_hi_galfa_cube.mir \
    out=hi/perseus_hi_galfa_cube_lee12_eq_regrid.mir \
    tin=av/perseus_av_lee12_eq.mir \
    axes=1,2

fits in=hi/perseus_hi_galfa_cube_lee12_eq_regrid.mir \
    out=hi/perseus_hi_galfa_cube_lee12_eq_regrid.fits \
    op=xyout


