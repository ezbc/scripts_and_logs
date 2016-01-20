#!/bin/csh

cd /d/bip3/ezbc/taurus/data/hi/DR2_archive_files/

set data = (GALFA_HI_RA+DEC_060.00+18.35_N.fits GALFA_HI_RA+DEC_060.00+26.35_N.fits GALFA_HI_RA+DEC_060.00+34.35_N.fits GALFA_HI_RA+DEC_068.00+18.35_N.fits GALFA_HI_RA+DEC_068.00+26.35_N.fits GALFA_HI_RA+DEC_068.00+34.35_N.fits GALFA_HI_RA+DEC_076.00+18.35_N.fits GALFA_HI_RA+DEC_076.00+26.35_N.fits GALFA_HI_RA+DEC_076.00+34.35_N.fits)

# Create DR1 cube
# ------------------------------------------------------------------------------
set data_dir = /d/bip2/lee/perseus_cloud/galfa/
cd /d/bip3/ezbc/taurus/data/hi/DR1_archive_files/

foreach filename ($data)
    ls $data_dir$filename
end

foreach filename ($data)
  set base=`basename $data_dir${filename} .fits`
  fits in=$data_dir${filename} op=xyin out=${base}.mir
end

foreach filename (GALFA*.mir)
  set base=`basename ${filename} .mir`
  puthd in=${base}.mir/restfreq value=1.420405752E+09 type='double'
end


foreach mirdata (GALFA_HI*.mir)
  set base=`basename ${mirdata} .mir`
  imsub in=${mirdata} "region=images(460,1590)" out=${base}.sub.mir
end

ls *.sub.mir

imcomb "in=*.sub.mir"\
out=../taurus_hi_galfa_dr1_cube.mir

fits in=../taurus_hi_galfa_dr1_cube.mir \
out=../taurus_hi_galfa_dr1_cube.fits \
op=xyout

# Create DR2 cube
# ------------------------------------------------------------------------------

cd /d/bip3/ezbc/taurus/data/hi/DR2_archive_files/

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
  imsub in=${mirdata} "region=images(460,1590)" out=${base}.sub.mir
end

imcomb "in=*.sub.mir"\
out=../taurus_hi_galfa_cube.mir

fits in=../taurus_hi_galfa_cube.mir \
out=../taurus_hi_galfa_cube.fits \
op=xyout


