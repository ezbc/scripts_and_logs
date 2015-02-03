#!/bin/csh

cd /d/bip3/ezbc/multicloud/data/hi/archive_files/

foreach filename (GALFA*.mir)
  set base=`basename ${filename} .mir`
  puthd in=${filename}/restfreq value=1.420405752E+09 type='double'
end

imcomb "in=*.sub.mir"\
out=../multicloud_hi_galfa_cube.mir

cd ..

regrid in=multicloud_hi_galfa_cube.mir \
out=multicloud_hi_galfa_cube_galcoord.mir \
op=galeqsw

fits in=multicloud_hi_galfa_cube.mir \
    out=multicloud_hi_galfa_cube.fits \
    op=xyout

fits in=multicloud_hi_galfa_cube_galcoord.mir \
    out=multicloud_hi_galfa_cube_galcoord.fits \
    op=xyout






