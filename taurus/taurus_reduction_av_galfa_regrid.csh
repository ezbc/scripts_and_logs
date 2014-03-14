#!/usr/bin/csh
# Av map is from Kainulainen et al. (2009, A&A, 508, 35)
# Provided by Jouni
# Extinction derived from background stars
# resolution is 2.4' 
# final res = (3.7^2 - 2.4^2)^0.5 = 2.816'

cd /d/bip3/ezbc/taurus/data/av

set av_dir=/d/bip3/ezbc/taurus/data/av
set galfa_dir=../galfa

fits in=$av_dir/taurus_av_kainulainen2009.fits \
    out=$av_dir/taurus_av_k09.mir \
    op=xyin

fits in=$galfa_dir/taurus_galfa_cube_bin_3.7arcmin.fits\
    out=$galfa_dir/taurus_galfa_cube_bin_3.7arcmin.mir \
    op=xyin

smooth in=$av_dir/taurus_av_k09.mir fwhm=2.816 pa=0 \
    out=$av_dir/taurus_av_k09_smooth.mir \
    scale=0.0 # keep units, do not change to /beam

regrid in=$av_dir/taurus_av_k09_smooth.mir \
    out=$av_dir/taurus_av_k09_regrid.mir \
    tin=$galfa_dir/taurus_galfa_cube_bin_3.7arcmin.mir \
    axes=1,2

fits in=$av_dir/taurus_av_k09_regrid.mir \
    out=$av_dir/taurus_av_k09_regrid.fits \
    op=xyout

# Pineda, J.~L., Goldsmith, P.~F., Chapman, N., et al.\ 2010, \apj, 721, 686 
# resolution of Av image is 200", see caption of figure 4
# thus we will convolve the Av image to the Gaussian filer size

fits in=$av_dir/taurus_av_pineda2010.fits \
    out=$av_dir/taurus_av_p10.mir \
    op=xyin

fits in=$galfa_dir/taurus.galfa.cube.bin.4arcmin.fits \
    out=$galfa_dir/taurus.galfa.cube.bin.4arcmin.mir \
    op=xyin

smooth map=$av_dir/taurus_av_p10.mir fwhm=240 pa=0  \
    out=taurus_av_p10_smooth.mir \
    options=final

regrid in=$av_dir/taurus_av_p10_convol.mir \
    out=$av_dir/taurus_av_p10_regrid.mir \
    tin=$galfa_dir/taurus.galfa.cube.bin.4arcmin.mir \
    axes=1,2

fits in=$av_dir/taurus_av_p10_regrid.mir \
    out=$av_dir/taurus_av_pineda2010_regrid.fits \
    op=xyout







