# Leo P combination of GMRT and VLA data

importgmrt(fitsfile='/home/elijah/research/leoP/gmrt/casa/reductionFiles/leopgmrt.dbcon',vis='leopGMRT.ms')

# GMRT Observation Properties
# ---------------------------
# channel width: 8.13802051
# bandwidth: 1041.66663
# Chans: 128

# VLA Observation Properties
# --------------------------
# channel width: 3.90625
# 

# plotted time vs. amp
# flags needed

flagdata(vis='leopGMRT.ms',antenna='C10&C12',scan='1')

flagdata(vis='leopGMRT.ms',antenna='C05&C06')

flagdata(vis='leopGMRT.ms',antenna='C08&12,C12&C13',timerange='2012/12/28/18:10:51.9~2012/12/28/21:03:20.1')

flagdata(vis='leopGMRT.ms',antenna='')



listobs(vis='leopCombined.ms')

concat(vis=['leop.ms','leopGMRT.ms'],concatvis='leopCombined.ms')

clean(vis=['leopGMRT.ms','leop.ms'],
	imagename='../images/dirty.v1',
	imsize=512,
	spw='0',
	cell='5arcsec',
	interactive=True,
	niter=10,
        outframe='LSRK',
	mode='channel',
	weighting='briggs',
	imagermode='csclean',
	restfreq='1420.405751 MHz')
