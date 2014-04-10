#!/bin/csh

cd /d/bip3/sstanimi/Taurus/data/

foreach filename (GALFA*.fits) 
  set base=`basename ${filename} .fits` 
  fits in=${filename} op=xyin out=${base}.mir 
end

foreach mirdata (GALFA_HI*.mir)
  set base=`basename ${mirdata} .mir`
  imsub in=${mirdata} "region=images(460,1590)" out=${base}.sub.mir
end

foreach mirdata (GALFA_HI*sub.mir)
  set base=`basename ${mirdata} .mir`
  imbin in=${mirdata} 'bin=2,2,2,2,4,4' out=${base}.sub.bin.mir
end

imcomb "in=*.sub.bin.mir" out=taurus_galfa_cube_bin.mir

fits in=taurus_galfa_cube_sub_bin.mir \
out=taurus_galfa_cube_sub_bin.fits \
op=xyout

#imcomb "in=*.sub.mir" out=taurus.galfa.cube.mir













# Making column density maps using different integrated velocity ranges

# velocity range: -100 to 100 km/s
#---------------------------------
# create moment 0 map
moment in=taurus.galfa.cube.mir out=taurus.galfa.mom0.ve1.mir mom=0

# change to column density map
maths "exp=taurus.galfa.mom0.vel1.mir*1.823/100" out=taurus.galfa.colDens.vel1.mir
puthd in=taurus.galfa.colDens.mir/bunit value='10^20cm^-2'

# write out as fits file:
fits in=taurus.galfa.colDens.mir op=xyout out=taurus.galfa.colDens.fits

# velocity range: -35 to 14.7 km/s
#---------------------------------
# create moment 0 map

moment in=taurus.galfa.cube.mir out=taurus.galfa.mom0.vel2.mir mom=0 region=images"(375,644)"

maths "exp=taurus.galfa.mom0.vel2.mir*1.823/100" out=taurus.galfa.colDens.vel2.mir

puthd in=taurus.galfa.colDens.vel2.mir/bunit value='10^20cm^-2'

# write out as fits file:
fits in=taurus.galfa.colDens.vel2.mir op=xyout out=taurus.galfa.colDens.vel2.fits

# velocity range: -6.5 to 11.5 km/s
#----------------------------------
# create moment 0 map
moment in=taurus.galfa.cube.mir out=taurus.galfa.mom0.vel3.mir mom=0 region=images"(529,627)"

# change to column density map
maths "exp=taurus.galfa.mom0.vel3.mir*1.823/100" out=taurus.galfa.colDens.vel3.mir
puthd in=taurus.galfa.colDens.vel3.mir/bunit value='10^20cm^-2'

# write out as fits file:
fits in=taurus.galfa.colDens.vel3.mir op=xyout out=taurus.galfa.colDens.vel3.fits

# create moment 1 map

moment in=taurus.galfa.cube.mir out=taurus.galfa.mom1.mir mom=1

# make smaller cube to been loaded into kvis

# imsub in=taurus.galfa.cube.mir "region=images(0,200)" out=taurus.galfa.cube.sub.mir


# Convolve cube to be compared with 8.4' res CfA CO data
# Requires for bmaj and bmin information first
# GALFA res: 3.9' x 4.1' --> radians = 3.9' * 60"/' /206265 "/rad

puthd in=taurus.galfa.cube.mir/bmaj value=0.0011926405
puthd in=taurus.galfa.cube.mir/bmin value=0.0011344629
puthd in=taurus.galfa.cube.mir/bpa value=0.

convol map="taurus.galfa.cube.mir" options="final" fwhm=504 pa=0 out=taurus.galfa.cube.cfaRes.mir

# put galfa cube on same grid as cfa

fits in=taurus.cfa.cube.fits op=xyin out=taurus.cfa.cube.mir

regrid in=taurus.galfa.cube.cfaRes.mir options="" axes="1,2" out=taurus.galfa.cube.cfaRes.regrid.mir tin=taurus.cfa.cube.mir

# write out fits file

fits in=taurus.galfa.cube.cfaRes.regrid.mir op=xyout out=taurus.galfa.cube.cfaRes.regrid.fits


