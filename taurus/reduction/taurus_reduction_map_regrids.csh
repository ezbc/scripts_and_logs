#!/usr/bin/csh

# Planck resolution is 5'. We will smooth and regrid all maps to the Planck map.

# Map                           Resolution      Convol size
# GALFA HI                      3.7'            3.363
# Kainulainen et al. (2009) Av  2.4'            4.386
# Pineda et al. (2010) Av       3.3'            3.75
# CfA CO Dame et al. (2001)     8'              No smoothing

cd /d/bip3/ezbc/taurus/data

# Load fits files
fits in=av/taurus_av_kainulainen2009_nan.fits \
    out=av/taurus_av_k09.mir \
    op=xyin

fits in=av/taurus_av_pineda2010.fits \
    out=av/taurus_av_p10.mir \
    op=xyin

fits in=hi/taurus_hi_galfa_cube.fits \
    out=hi/taurus_hi_galfa_cube.mir \
    op=xyin

fits in=av/taurus_av_planck.fits \
    out=av/taurus_av_planck.mir \
    op=xyin

fits in=av/taurus_av_error_planck.fits \
    out=av/taurus_av_error_planck.mir \
    op=xyin

fits in=co/taurus_cfa_cube_eqcoord.fits \
    out=co/taurus_co_cfa_cube.mir \
    op=xyin

# Regrid Planck image to have pixels = 1 resolution element = 5'
# Regrid task needs pixel size in degrees
# 5' / 60"/deg = 0.08333 deg
regrid in=av/taurus_av_planck.mir \
    out=av/taurus_av_planck_5arcmin.mir \
    desc=77.5,0,-0.08333,240,18,0,0.08333,180

regrid in=av/taurus_av_error_planck.mir \
    out=av/taurus_av_error_planck_5arcmin.mir \
    desc=77.5,0,-0.08333,240,18,0,0.08333,180

# regrid galfa hi image to have lower velocity resolution
regrid in=hi/taurus_hi_galfa_cube.mir \
    out=hi/taurus_hi_galfa_cube_1kms.mir \
    desc=77.5,0,-0.08333,240,18,0,0.08333,180,0,100,1,200

# Smooth the files
smooth in=av/taurus_av_k09.mir fwhm=4.386 pa=0 \
    out=av/taurus_av_k09_smooth_planckres.mir \
    scale=0.0 # keep units, do not change to /beam

smooth in=av/taurus_av_p10.mir fwhm=3.75 pa=0 \
    out=av/taurus_av_p10_smooth_planckres.mir \
    scale=0.0 # keep units, do not change to /beam

smooth in=hi/taurus_hi_galfa_cube_1kms.mir fwhm=3.363 pa=0 \
    out=hi/taurus_hi_galfa_cube_smooth_planckres.mir \
    scale=0.0 # keep units, do not change to /beam

# Regrid the images onto the planck grid
regrid in=av/taurus_av_k09_smooth_planckres.mir \
    out=av/taurus_av_k09_regrid_planckres.mir \
    tin=av/taurus_av_planck_5arcmin.mir

regrid in=av/taurus_av_p10_smooth_planckres.mir \
    out=av/taurus_av_p10_regrid_planckres.mir \
    tin=av/taurus_av_planck_5arcmin.mir

regrid in=hi/taurus_hi_galfa_cube_smooth_planckres.mir \
    out=hi/taurus_hi_galfa_cube_regrid_planckres.mir \
    tin=av/taurus_av_planck_5arcmin.mir \
    axes=1,2

regrid in=co/taurus_co_cfa_cube.mir \
    out=co/taurus_co_cfa_cube_regrid_planckres.mir \
    tin=av/taurus_av_planck_5arcmin.mir \
    axes=1,2

# write fits images out
fits in=av/taurus_av_planck_5arcmin.mir \
    out=av/taurus_av_planck_5arcmin.fits \
    op=xyout

fits in=av/taurus_av_error_planck_5arcmin.mir \
    out=av/taurus_av_error_planck_5arcmin.fits \
    op=xyout

fits in=av/taurus_av_k09_regrid_planckres.mir \
    out=av/taurus_av_k09_regrid_planckres.fits \
    op=xyout

fits in=av/taurus_av_p10_regrid_planckres.mir \
    out=av/taurus_av_p10_regrid_planckres.fits \
    op=xyout

fits in=hi/taurus_hi_galfa_cube_regrid_planckres.mir \
    out=hi/taurus_hi_galfa_cube_regrid_planckres.fits \
    op=xyout

fits in=co/taurus_co_cfa_cube_regrid_planckres.mir \
    out=co/taurus_co_cfa_cube_regrid_planckres.fits \
    op=xyout


