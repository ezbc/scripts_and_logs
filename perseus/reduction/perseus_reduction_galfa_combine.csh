#!/bin/csh

cd /d/bip3/ezbc/perseus/data/hi/archive_files/

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

imcomb "in=*DEC_060*.sub.mir,*DEC_052*.sub.mir,*DEC_044*.sub.mir,*DEC_036*.sub.mir"\
out=../perseus_hi_galfa_cube.mir

cd ..

regrid in=perseus_hi_galfa_cube.mir \
out=perseus_hi_galfa_cube_galcoord.mir \
op=galeqsw

#imcomb "in=*.mir" out=perseus.galfa.cube.mir

foreach filename (*galfa*.mir)
  set base=`basename ${filename} .fits`
  puthd in=${filename}/restfreq value=1.420405752E+09 type='double'
end

fits in=perseus_hi_galfa_cube.mir \
    out=perseus_hi_galfa_cube.fits \
    op=xyout

fits in=perseus_hi_galfa_cube_galcoord.mir \
    out=perseus_hi_galfa_cube_galcoord.fits \
    op=xyout






