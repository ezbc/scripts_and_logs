#!/usr/bin/csh

# Planck resolution is 5'. We will smooth and regrid all maps to the Planck map.

# Convol size is determined by (5'^2 - convolsize'^2)^0.5

# Map                           Resolution      Convol size
# GALFA HI                      3.7'            3.363
# Lee et al. (2012) Av          5'              5
# CfA CO Dame et al. (2001)     8'              No smoothing

cd /d/bip3/ezbc/perseus/data

# Load fits files
fits in=av/perseus_av_lee12_masked.fits \
    out=av/perseus_av_lee12_masked.mir \
    op=xyin

fits in=av/perseus_av_planck.fits \
    out=av/perseus_av_planck.mir \
    op=xyin

fits in=hi/perseus_hi_galfa_cube.fits \
    out=hi/perseus_hi_galfa_cube.mir \
    op=xyin

fits in=../../taurus/data/co/taurus_cfa_cube_eqcoord.fits \
    out=co/perseus_co_cfa_cube.mir \
    op=xyin

# Regrid Planck image to have pixels = 1 resolution element = 5'
# Regrid task needs pixel size in degrees
# 5' / 60"/deg = 0.08333 deg
regrid in=av/perseus_av_planck.mir \
    out=av/perseus_av_planck_5arcmin.mir \
    desc=59.75,0,-0.08333,180,26.05,0,0.08333,132

# regrid galfa hi image to have lower velocity resolution
regrid in=hi/perseus_hi_galfa_cube.mir \
    out=hi/perseus_hi_galfa_cube_1kms.mir \
    desc=59.75,0,-0.08333,180,26.05,0,0.08333,132,0,100,1,200

# Smooth the files
smooth in=av/perseus_av_lee12_masked.mir fwhm=5 pa=0 \
    out=av/perseus_av_lee12_masked_smooth_planckres.mir \
    scale=0.0 # keep units, do not change to /beam

smooth in=hi/perseus_hi_galfa_cube_1kms.mir fwhm=3.363 pa=0 \
    out=hi/perseus_hi_galfa_cube_smooth_planckres.mir \
    scale=0.0 # keep units, do not change to /beam

# Regrid the images onto the planck grid
regrid in=av/perseus_av_lee12_masked_smooth_planckres.mir \
    out=av/perseus_av_lee12_masked_regrid_planckres.mir \
    tin=av/perseus_av_planck_5arcmin.mir

regrid in=hi/perseus_hi_galfa_cube_smooth_planckres.mir \
    out=hi/perseus_hi_galfa_cube_regrid_planckres.mir \
    tin=av/perseus_av_planck_5arcmin.mir \
    axes=1,2

regrid in=co/perseus_co_cfa_cube.mir \
    out=co/perseus_co_cfa_cube_regrid_planckres.mir \
    tin=av/perseus_av_planck_5arcmin.mir \
    axes=1,2

# write fits images out
fits in=av/perseus_av_planck_5arcmin.mir \
    out=av/perseus_av_planck_5arcmin.fits \
    op=xyout

fits in=av/perseus_av_lee12_masked_regrid_planckres.mir \
    out=av/perseus_av_lee12_maksed_regrid_planckres.fits \
    op=xyout

fits in=hi/perseus_hi_galfa_cube_regrid_planckres.mir \
    out=hi/perseus_hi_galfa_cube_regrid_planckres.fits \
    op=xyout

fits in=co/perseus_co_cfa_cube_regrid_planckres.mir \
    out=co/perseus_co_cfa_cube_regrid_planckres.fits \
    op=xyout


