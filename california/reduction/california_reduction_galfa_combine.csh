#!/bin/csh

cd /d/bip3/ezbc/california/data/hi/archive_files/

set data = (GALFA_HI_RA+DEC_052.00+26.35_N.fits GALFA_HI_RA+DEC_052.00+34.35_N.fits GALFA_HI_RA+DEC_060.00+26.35_N.fits GALFA_HI_RA+DEC_060.00+34.35_N.fits GALFA_HI_RA+DEC_068.00+26.35_N.fits GALFA_HI_RA+DEC_068.00+34.35_N.fits GALFA_HI_RA+DEC_076.00+26.35_N.fits GALFA_HI_RA+DEC_076.00+34.35_N.fits)

set data_dir = /d/bip2/DR2W_v1/Narrow/

foreach filename ($data)
  set base=`basename $data_dir${filename} .fits`
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

imcomb "in=*sub.mir"\
out=../california_hi_galfa_cube.mir

fits in=../california_hi_galfa_cube.mir \
out=../california_hi_galfa_cube.fits \
op=xyout






