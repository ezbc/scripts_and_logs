#!/usr/bin/csh

cd /d/bip3/ezbc/perseus/data/galfa/

set planck_dir=../av
set galfa_dir=../galfa

fits in=$planck_dir/perseus_planck_av.fits \
        out=$planck_dir/perseus_planck_av.mir \
        op=xyin

regrid in=$planck_dir/perseus_planck_av.mir \
    tin=$galfa_dir/perseus_galfa_cube_bin_3.7arcmin.mir \
    out=$planck_dir/perseus_planck_av_regrid.mir \
    axes=1,2

fits in=$planck_dir/perseus_planck_av_regrid.mir \
        out=$planck_dir/perseus_planck_av_regrid.fits \
        op=xyout

# Regrid of planck onto 5' pixels

regrid in=perseus_planck_av.mir \
    out=perseus_planck_av_5arcmin.mir \
    desc=65,0,-0.08324,300,21,0,0.08324,280

# 2mass Av image has same Planck 5' resolution

# regrid 2mass Av image from Lee+12 onto planck image
regrid in=../2mass/2mass_av_lee12.mir \
    out=../2mass/2mass_av_lee12_planck_regrid.mir \
    tin=../av/perseus_planck_av_5arcmin.mir

regrid in=../2mass/2mass_av_lee12_nocal.mir \
    out=../2mass/2mass_av_lee12_nocal_planck_regrid.mir \
    tin=../av/perseus_planck_av_5arcmin.mir

fits in=../2mass/2mass_av_lee12_planck_regrid.mir \
    out=../2mass/2mass_av_lee12_planck_regrid.fits \
    op=xyout

fits in=../2mass/2mass_av_lee12_nocal_planck_regrid.mir \
    out=../2mass/2mass_av_lee12_nocal_planck_regrid.fits \
    op=xyout

fits in=../av/perseus_planck_av_5arcmin.mir \
    out=../av/perseus_planck_av_5arcmin.fits \
    op=xyout

# regrid of planck test regions

set planck_dir=../av

fits in=$galfa_dir/perseus_nhi_lee12.fits \
    out=$galfa_dir/perseus_nhi_lee12.mir \
    op=xyin

imhead in=$galfa_dir/perseus_nhi_lee12.mir

# Write rest frequency to nhi
puthd in=$galfa_dir/perseus_nhi_lee12.mir/restfreq \
    value=1.420405752E+09 \
    type='double'

fits in=$planck_dir/perseus_planck_region.fits \
    out=$planck_dir/perseus_planck_region.mir \
    op=xyin

imhead in=$planck_dir/perseus_planck_region.mir

#regrid in=$planck_dir/perseus_planck_region.mir \
#    out=$planck_dir/perseus_planck_eqcoords.mir \
#    options=galeqsw 

regrid in=$planck_dir/perseus_planck_region.mir \
    tin=$galfa_dir/perseus_nhi_lee12.mir \
    out=$planck_dir/perseus_planck_lee12_regrid.mir \
    axes=1,2 

imhead in=$planck_dir/perseus_planck_region.mir

fits in=$planck_dir/perseus_planck_lee12_regrid.mir \
    out=$planck_dir/perseus_planck_lee12_regrid.fits \
    op=xyout

