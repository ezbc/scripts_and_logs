#!/usr/bin/csh

cd /d/bip3/ezbc/perseus/data/galfa/

puthd in=perseus_galfa_cube.mir/restfreq value=1.420405752E+09 type='double'

# beamsize of cube: 4' x 3.5'
# geometric mean = (4*3.5)**0.5 = 3.7'
# degrees: 

regrid in=perseus_galfa_cube_sub.mir \
    out=perseus_galfa_cube_bin_3.7arcmin.mir \
    desc=65,0,-0.0616,400,21,0,0.0616,338,0,100,1,200

imbin in=perseus_galfa_cube_sub.mir \
    out=perseus_galfa_cube_bin_4arcmin.mir bin=2,2,2,2,2,2

fits in=perseus_galfa_cube_bin_4arcmin.mir \
    out=perseus_galfa_cube_bin_4arcmin.fits op=xyout

fits in=perseus_galfa_cube_bin_3.7arcmin.mir \
    out=perseus_galfa_cube_bin_3.7arcmin.fits op=xyout


fits in=perseus_galfa_cube.mir out=perseus_galfa_cube.fits op=xyout

set cfa=../cfa
set tauruscfa=/d/bip3/ezbc/taurus/data/cfa
set galfa=../galfa

# grid spatial axes of cfa onto galfa 4arcmin
regrid in=$tauruscfa/taurus.cfa.cube.mir \
    tin=$galfa/perseus.galfa.cube.bin.4arcmin.mir \
    out=$cfa/perseus.cfa.cube.galfaBin.4arcmin.mir \
    axes=1,2

# Regrid CfA map Min's HI image
fits in=$galfa/perseus.lee12.fits out=$galfa/perseus.lee12.mir op=xyin

puthd in=perseus.lee12.mir/restfreq value=1.420405752E+09 type='double'

# regrid Min's HI image to galfa 4arcmin
regrid in=$galfa/perseus.lee12.mir \
    tin=$galfa/perseus.galfa.cube.bin.4arcmin.mir \
    out=$galfa/perseus.lee12.galfaBin.4arcmin.mir axes=1,2

# regrid cfa to min's hi map
regrid in=../../../taurus/data/cfa/taurus.cfa.cube.eqcoord.mir tin=$galfa/perseus.lee12.mir out=$cfa/perseus.cfa.cube.eqcoord.mir axes=1,2

# write fits files
fits in=perseus.cfa.cube.galfaBin.4arcmin.mir out=perseus.cfa.cube.galfaBin.4arcmin.fits op=xyout 

fits in=$cfa/perseus.cfa.cube.galfaBin.4arcmin.mir \
    out=$cfa/perseus.cfa.cube.galfaBin.4arcmin.fits op=xyout 

fits in=$galfa/perseus.lee12.galfaBin.4arcmin.mir \
    out=$galfa/perseus.lee12.galfaBin.4arcmin.fits op=xyout 


# import and grid the Galfa cube used in Lee et al. (2012)
fits in=/d/leffe2/lee/perseus_cloud/galfa/galfa_perseus.sub2.fits \
    out=perseus.lee12.cube.mir op=xyin

regrid in=$galfa/perseus.lee12.cube.mir \
    tin=$galfa/perseus.galfa.cube.bin.4arcmin.mir \
    out=$galfa/perseus.lee12.cube.galfaBin.4arcmin.mir axes=1,2

fits in=$galfa/perseus.lee12.cube.galfaBin.4arcmin.mir \
    out=$galfa/perseus.lee12.cube.galfaBin.4arcmin.fits op=xyout 








