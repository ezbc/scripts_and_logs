#!/bin/csh

cd /d/bip3/ezbc/california/data/hi/archive_files/

foreach filename (GALFA*.fits)
  set base=`basename ${filename} .fits`
  fits in=${filename} op=xyin out=${base}.mir
end

foreach filename (GALFA*.fits)
  set base=`basename ${filename} .fits`
  puthd in=${base}.mir/restfreq value=1.420405752E+09 type='double'
end

foreach mirdata (GALFA_HI*.mir)
  set base=`basename ${mirdata} .mir`
  imsub in=${mirdata} "region=images(460,1590)" out=${base}_sub.mir
end

imcomb "in=*_sub.mir"\
out=../california_hi_galfa_cube.mir

fits in=../california_hi_galfa_cube.mir \
out=../california_hi_galfa_cube.fits \
op=xyout






