cd /d/bip3/ezbc/galactic_anticenter

# Load fits files
fits in=anti_center_celestial.fits \
    out=anti_center_celestial.mir \
    op=xyin

puthd in=anti_center_celestial.mir/restfreq \
    value=1.420405752E+09 \
    type='double'

regrid in=anti_center_celestial.mir \
    out=anti_center_galactic.mir \
    axes=1,2 \
    options=galeqsw

rm -rf anti_center_galactic_crop*

regrid in=anti_center_galactic.mir \
    out=anti_center_galactic_crop.mir \
    desc=123.31427,0,-0.01666,3120,58.330146,0,0.01666,1250,-.09901,0,0.00215,68

fits in=anti_center_galactic_crop.mir \
    out=anti_center_galactic_crop.fits \
    op=xyout


