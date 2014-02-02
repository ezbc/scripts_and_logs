#!/bin/csh

cd /d/bip3/sstanimi/Taurus/data/IRIS

foreach filename (IRIS*.fits)
  set base=`basename ${filename} .fits`
  fits in=${filename} op=xyin out=${base}.mir
end

imcomb "in=IRIS*.mir" out=IRIS.comboImage.mir

# imsub in=taurus.galfa.cube.mir "region=images(0,200)" out=taurus.galfa.cube.sub.mir

# put galfa cube on same grid as cfa
\op
fits in=taurus.cfa.cube.fits op=xyin out=taurus.cfa.cube.mir

regrid in=taurus.galfa.cube.cfaRes.mir options="" axes="1,2" out=taurus.galfa.cube.cfaRes.regrid.mir tin=taurus.cfa.cube.mir

# write out fits file

fits in=taurus.galfa.cube.cfaRes.regrid.mir op=xyout out=taurus.galfa.cube.cfaRes.regrid.fits


