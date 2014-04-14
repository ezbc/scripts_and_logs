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

