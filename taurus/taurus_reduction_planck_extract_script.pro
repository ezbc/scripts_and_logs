pro planckImageExtraction

.r rd_planckmap

rd_planckmap, '857', l, b, map, hdr

read_fits_map, 'HFI_SkyMap_857_2048_R1.10_nominal.fits', map, hdr
write_fits_map, 'HFI_SkyMap_857_2048_R1.10_nominal.fits', map, hdr, coordsys=C

rd_planckmap, '353', l, b, map, hdr
write_fits_map, 'planck.353ghz.fits', map, hdr

rd_planckmap, '545', l, b, map, hdr
write_fits_map, 'planck.545ghz.fits', map, hdr

end