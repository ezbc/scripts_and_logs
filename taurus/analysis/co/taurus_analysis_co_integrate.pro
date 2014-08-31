pro co_Int
Co = readfits('../data/DHT21_Taurus_interp.reorder.fits', hd)

Ico = total(Co,3)*0.65

sxaddpar, hd, 'NAXIS', 2
sxdelpar, hd, 'NAXIS3' 
sxdelpar, hd, 'CTYPE3'
sxdelpar, hd, 'CRVAL3'
sxdelpar, hd, 'CRPIX3'
sxdelpar, hd, 'CDELT3'
sxaddpar, hd, 'BUNIT', 'K km/s'

writefits, '../data/taurus.cfa.co_Int.fits', Ico, hd

stop
end
