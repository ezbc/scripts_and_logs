#!/usr/bin/csh

# relevant cube: taurus.galfa.cube.mir
# We need to regrid to 4"/pix
# Cell size is 1"/pix, thus we bin 4 pixels together

puthd in=taurus.galfa.cube.mir/restfreq value=1.420405752E+09 type='double'

puthd in=taurus.galfa.cube.bin.mir/restfreq value=1.420405752E+09 type='double'

imbin in=/d/bip3/ezbc/taurus/data/galfa/taurus.galfa.cube.mir out=taurus.galfa.cube.bin.4arcmin.mir bin=4,4,4,4,4,4

#imbin in=taurus.galfa.cube.mir out=taurus.galfa.cube.bin.4arcmin.mir 'bin='

imbin in=taurus.galfa.cube.mir out=taurus.galfa.cube.bin.16arcmin.mir 'bin=8,8,8,8,4,4'

imbin in=taurus.galfa.cube.mir out=taurus.galfa.cube.bin.32arcmin.mir 'bin=16,16,16,16,4,4'

fits in=taurus.galfa.cube.mir out=taurus.galfa.cube.fits op=xyout

fits in=taurus.galfa.cube.bin.4arcmin.mir out=taurus.galfa.cube.bin.4arcmin.fits op=xyout

fits in=taurus.galfa.cube.bin.16arcmin.mir op=xyout out=taurus.galfa.cube.bin.16arcmin.fits

fits in=taurus.galfa.cube.bin.32arcmin.mir op=xyout out=taurus.galfa.cube.bin.32arcmin.fits

set cfa=../cfa
set galfa=../galfa

# Regrid CfA map to thi new cube
fits in=$cfa/DHT21_Taurus_interp.fits op=xyin out=$cfa/temp.cfacube.mir

reorder in=$cfa/temp.cfacube.mir out=$cfa/taurus.cfa.cube.mir mode=231

puthd in=$cfa/taurus.cfa.cube.mir/crval3 value=-0.02178E+03 type='double'

puthd in=$cfa/taurus.cfa.cube.mir/cdelt3 value=.6502 type='double'

puthd in=$cfa/taurus.cfa.cube.mir/restfreq value=115.2712018E+09 type='double'

regrid in=$cfa/taurus.cfa.cube.mir out=$cfa/taurus.cfa.cube.eqcoord.mir op=galeqsw

regrid in=$cfa/taurus.cfa.cube.eqcoord.mir tin=$galfa/taurus.galfa.cube.bin.4arcmin.mir out=$cfa/taurus.cfa.galfaBin.4arcmin.mir axes=1,2

regrid in=$cfa/taurus.cfa.cube.eqcoord.mir tin=$galfa/taurus.galfa.cube.bin.16arcmin.mir out=$cfa/taurus.cfa.galfaBin.16arcmin.mir axes=1,2

fits in=taurus.cfa.galfaBin.4arcmin.mir out=taurus.cfa.galfaBin.4arcmin.fits op=xyout 

fits in=taurus.cfa.galfaBin.16arcmin.mir out=taurus.cfa.galfaBin.16arcmin.fits op=xyout 



