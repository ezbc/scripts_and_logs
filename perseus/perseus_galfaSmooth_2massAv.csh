#!/usr/bin/csh

cd /d/bip3/ezbc/perseus/data/galfa

set twomass_dir=/d/bip3/ezbc/perseus/data/2mass
set galfa_dir=/d/bip3/ezbc/perseus/data/galfa

fits in=$twomass_dir/2mass_av_lee12.fits out=$twomass_dir/2mass_av_lee12.mir \
    op=xyin

fits in=$galfa_dir/perseus.galfa.cube.bin.4arcmin.fits \
    out=$galfa_dir/perseus.galfa.cube.bin.4arcmin.mir \
    op=xyin

smooth in=$galfa_dir/perseus.galfa.cube.bin.4arcmin.mir\
    out=$galfa_dir/perseus.galfa.2massSmooth.mir \
    fwhm=300\
    pa=0

regrid in=$twomass_dir/2mass_av_lee12.mir \
    out=$twomass_dir/2mass_av_lee12_regrid.mir \
    tin=$galfa_dir/perseus.galfa.2massSmooth.mir \
    axes=1,2

fits in=$twomass_dir/2mass_av_lee12_regrid.mir \
    out=$twomass_dir/2mass_av_lee12_regrid.fits \
    op=xyout

# Used incorrect Av image, previous Av image had already been 0 point calibrated
# Correct image needed for calculating the NHI/Av correlation is here:
# /d/leffe2/lee/perseus_cloud/Min/R_H2/isolate_HI_102311/PerA_Extn2MASS_F_Eq.congrid5.0.SNR.fits
# copied to /d/bip3/ezbc/perseus/data/2mass/2mass_av_lee12_nocal.fits
fits in=$twomass_dir/2mass_av_lee12_nocal.fits \
    out=$twomass_dir/2mass_av_lee12_nocal.mir \
    op=xyin

regrid in=$twomass_dir/2mass_av_lee12_nocal.mir \
    out=$twomass_dir/2mass_av_lee12_nocal_regrid.mir \
    tin=$galfa_dir/perseus.galfa.2massSmooth.mir \
    axes=1,2

fits in=$twomass_dir/2mass_av_lee12_nocal_regrid.mir \
    out=$twomass_dir/2mass_av_lee12_nocal_regrid.fits \
    op=xyout

fits in=$twomass_dir/2mass_av_lee12_nocal_SNR.fits \
    out=$twomass_dir/2mass_av_lee12_nocal_SNR.mir \
    op=xyin

regrid in=$twomass_dir/2mass_av_lee12_nocal_SNR.mir \
    out=$twomass_dir/2mass_av_lee12_nocal_SNR_regrid.mir \
    tin=$galfa_dir/perseus.galfa.2massSmooth.mir \
    axes=1,2

fits in=$twomass_dir/2mass_av_lee12_nocal_SNR_regrid.mir \
    out=$twomass_dir/2mass_av_lee12_nocal_SNR_regrid.fits \
    op=xyout







