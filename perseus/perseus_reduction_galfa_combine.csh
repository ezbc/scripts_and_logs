#!/bin/csh

cd /d/bip3/ezbc/perseus/data/galfa/

foreach filename (GALFA*.fits)
  set base=`basename ${filename} .fits`
  fits in=${filename} op=xyin out=${base}.mir
end

foreach mirdata (*.mir)
  set base=`basename ${mirdata} .mir`
  imsub in=${mirdata} "region=images(460,1590)" out=${base}.sub.mir
end

foreach filename (GALFA*.mir)
  set base=`basename ${filename} .fits`
  puthd in=${filename}/restfreq value=1.420405752E+09 type='double'
end

imcomb "in=*44*34*.sub.mir,*044*26*.sub.mir,*52*26*sub.mir,*52*34*sub.mir" out=perseus.galfa.cube.sub.mir

regrid in=perseus.galfa.cube.sub.mir out=perseus.galfa.cube.sub.regrid.mir op=galeqsw

#imcomb "in=*.mir" out=perseus.galfa.cube.mir

foreach filename (*galfa*.mir)
  set base=`basename ${filename} .fits`
  puthd in=${filename}/restfreq value=1.420405752E+09 type='double'
end

fits in=perseus.galfa.cube.sub.regrid.mir \
    out=perseus.galfa.cube.sub.regrid.fits op=xyout

fits in=perseus.galfa.cube.mir out=perseus.galfa.cube.fits \
    op=xyout






