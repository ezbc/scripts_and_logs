#!/bin/csh

cd /d/bip3/ezbc/taurus/data/2mass

foreach filename (aH*.fits)
  set base=`basename ${filename} .fits`
  fits in=${filename} op=xyin out=${base}.mir
end

imcomb "in=*.mir" out=taurus.2mass.h.mir

# Making column density maps using different integrated velocity ranges

# write out as fits file:
fits in=taurus.2mass.h.mir op=xyout out=taurus.2mass.h.fits

foreach filename (aJ*.fits)
  set base=`basename ${filename} .fits`
  fits in=${filename} op=xyin out=${base}.mir
end

imcomb "in=*.mir" out=taurus.2mass.j.mir

# Making column density maps using different integrated velocity ranges

# write out as fits file:
fits in=taurus.2mass.j.mir op=xyout out=taurus.2mass.j.fits

foreach filename (aK*.fits)
  set base=`basename ${filename} .fits`
  fits in=${filename} op=xyin out=${base}.mir
end

imcomb "in=*.mir" out=taurus.2mass.k.mir

# Making column density maps using different integrated velocity ranges

# write out as fits file:
fits in=taurus.2mass.k.mir op=xyout out=taurus.2mass.k.fits

