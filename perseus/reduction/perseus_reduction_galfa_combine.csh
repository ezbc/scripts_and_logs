#!/bin/csh

# Define the cubes needed to stitch together
set data = (GALFA_HI_RA+DEC_020.00+18.35_N.fits GALFA_HI_RA+DEC_020.00+26.35_N.fits GALFA_HI_RA+DEC_020.00+34.35_N.fits GALFA_HI_RA+DEC_028.00+18.35_N.fits GALFA_HI_RA+DEC_028.00+26.35_N.fits GALFA_HI_RA+DEC_028.00+34.35_N.fits GALFA_HI_RA+DEC_036.00+18.35_N.fits GALFA_HI_RA+DEC_036.00+26.35_N.fits GALFA_HI_RA+DEC_036.00+34.35_N.fits GALFA_HI_RA+DEC_044.00+18.35_N.fits GALFA_HI_RA+DEC_044.00+26.35_N.fits GALFA_HI_RA+DEC_044.00+34.35_N.fits GALFA_HI_RA+DEC_052.00+18.35_N.fits GALFA_HI_RA+DEC_052.00+26.35_N.fits GALFA_HI_RA+DEC_052.00+34.35_N.fits GALFA_HI_RA+DEC_060.00+18.35_N.fits GALFA_HI_RA+DEC_060.00+26.35_N.fits GALFA_HI_RA+DEC_060.00+34.35_N.fits)

# Create DR2 cube
# ------------------------------------------------------------------------------

cd /d/bip3/ezbc/perseus/data/hi/DR2_archive_files/

set data_dir = /d/bip2/DR2W_v1/Narrow/

foreach filename ($data)
  set base=`basename $data_dir${filename} .fits`
  echo loading: $data_dir${filename}
  echo to: ${base}.mir
  fits in=$data_dir${filename} op=xyin out=${base}.mir
end

foreach filename ($data)
  set base=`basename ${filename} .fits`
  echo loading cube: ${base}.mir
  puthd in=${base}.mir/restfreq value=1.420405752E+09 type='double'
end

foreach mirdata (GALFA_HI*.mir)
  set base=`basename ${mirdata} .mir`
  imsub in=${mirdata} "region=images(460,1590)" out=${base}.sub.mir
end

imcomb "in=*.sub.mir"\
out=../perseus_hi_galfa_cube.mir

fits in=../perseus_hi_galfa_cube.mir \
out=../perseus_hi_galfa_cube.fits \
op=xyout


# Create DR1 cube
# ------------------------------------------------------------------------------
#set data_dir = /d/bip2/lee/perseus_cloud/galfa/
set data_dir = /d/bip3/ezbc/galfa/DR1/
cd /d/bip3/ezbc/perseus/data/hi/DR1_archive_files/

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
out=../perseus_hi_galfa_dr1_cube.mir

fits in=../perseus_hi_galfa_dr1_cube.mir \
out=../perseus_hi_galfa_dr1_cube.fits \
op=xyout


