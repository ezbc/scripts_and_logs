#!/usr/bin/csh

# Planck resolution is 5'. We will smooth and regrid all maps to the Planck map.

# Convol size is determined by (5'^2 - convolsize'^2)^0.5

# Map                           Resolution      Convol size
# GALFA HI                      3.7'            3.363
# CfA CO Dame et al. (2001)     8'              No smoothing

cd /d/bip3/ezbc/california/data

# Load fits files
fits in=av/california_av_planck.fits \
    out=av/california_av_planck.mir \
    op=xyin

fits in=av/california_av_error_planck.fits \
    out=av/california_av_error_planck.mir \
    op=xyin

fits in=hi/california_hi_galfa_cube.fits \
    out=hi/california_hi_galfa_cube.mir \
    op=xyin

fits in=../../california/data/co/california_cfa_cube_eqcoord.fits \
    out=co/california_co_cfa_cube.mir \
    op=xyin

# Regrid Planck image to have pixels = 1 resolution element = 5'
# Regrid task needs pixel size in degrees
# 5' / 60"/deg = 0.08333 deg
regrid in=av/california_av_planck.mir \
    out=av/california_av_planck_5arcmin.mir \
    desc=74.75,0,-0.08333,159,28.066,0,0.08333,144

regrid in=av/california_av_error_planck.mir \
    out=av/california_av_error_planck_5arcmin.mir \
    desc=77.5,0,-0.08333,240,18,0,0.08333,180

# regrid galfa hi image to have lower velocity resolution
regrid in=hi/california_hi_galfa_cube.mir \
    out=hi/california_hi_galfa_cube_1kms.mir \
    desc=74.75,0,-0.08333,159,28.066,0,0.08333,144,0,100,1,200

# Smooth the files
smooth in=hi/california_hi_galfa_cube_1kms.mir fwhm=3.363 pa=0 \
    out=hi/california_hi_galfa_cube_smooth_planckres.mir \
    scale=0.0 # keep units, do not change to /beam

# Regrid the images onto the planck grid
regrid in=hi/california_hi_galfa_cube_smooth_planckres.mir \
    out=hi/california_hi_galfa_cube_regrid_planckres.mir \
    tin=av/california_av_planck_5arcmin.mir \
    axes=1,2

regrid in=co/california_co_cfa_cube.mir \
    out=co/california_co_cfa_cube_regrid_planckres.mir \
    tin=av/california_av_planck_5arcmin.mir \
    axes=1,2

# write fits images out
fits in=av/california_av_planck_5arcmin.mir \
    out=av/california_av_planck_5arcmin.fits \
    op=xyout

fits in=av/california_av_error_planck_5arcmin.mir \
    out=av/california_av_error_planck_5arcmin.fits \
    op=xyout

fits in=hi/california_hi_galfa_cube_regrid_planckres.mir \
    out=hi/california_hi_galfa_cube_regrid_planckres.fits \
    op=xyout

fits in=co/california_co_cfa_cube_regrid_planckres.mir \
    out=co/california_co_cfa_cube_regrid_planckres.fits \
    op=xyout


