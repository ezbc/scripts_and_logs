#!/bin/csh

################################################################
#  Log file for combination of VLA and GMRT datasets in miriad
################################################################

# the relevant continuum-subtracted datasets are in the following
# directory: /Users/research/data/leop/finalUVdataSets/
# leopGMRT.contsub.fits
# leopVLAd.contsub.fits
# leopVLAc.contsub.fits

# an important note: these have been written out from CASA using
# exportuvfits with the multisource parameter set to false 

miriad
cd /d/bip3/ezbc/leop/data/hi/miriad/reductionFiles/
uvDir=/d/bip3/ezbc/leop/data/hi/finalUVdataSets

# load the datasets in:
task fits
reset
in=$uvDir/leopGMRT.contsub.fits
op=uvin 
out=GMRT
inp
go

task fits
reset
unset options
in=$uvDir/leopGMRT.contsub.normwts.fits
op=uvin 
out=GMRTnormwts
inp
go

task fits 
reset
in=$uvDir/leopVLAd.contsub.fits
op=uvin 
out=VLAd
go

task fits 
reset
in=$uvDir/leopVLAc.contsub.fits
op=uvin 
out=VLAc
go

task fits 
reset
in=$uvDir/leopVLAb1.contsub.fits
op=uvin 
out=VLAb1
go

task fits 
reset
in=$uvDir/leopVLAb2.contsub.fits
op=uvin 
out=VLAb2
go





task invert
reset
vis=GMRT,VLAd,VLAc
map=leop.combine.2pt5kms.map
beam=leop.combine.2pt5kms.beam
imsize=512
cell=4
robust=2
line=velocity,70,185,2.5,2.5
inp
go

# Following notes frm Steve Warren aboyt MIRIAD invert of VLA 
task invert
reset
vis=GMRT,VLAd,VLAc
map=leop.combine.1pt8kms.map
beam=leop.combine.1pt8kms.beam
imsize=512
cell=4
robust=1
line=velocity,70,221,1.8,1.8
inp
go

# not good, VLAc+GMRT obs had:
#Visibilities accepted: 4917507
### Warning [invert]:  Visibilities rejected: 3664589

# and all three obs led to:
#Visibilities accepted: 2134718
### Warning [invert]:  Visibilities rejected: 1115710

# why do all three obs have fewer visibilities accepted and rejected?

# Compare against setting options=mosaic in invert
task invert
reset
vis=GMRT,VLAd,VLAc
map=leop.combine.mosaic.map
beam=leop.combine.mosaic.beam
imsize=512
cell=4
robust=1
line=velocity,70,221,1.8,1.8
options=mosaic
go

# mosaiced cube has higher noise level
# disregarding mosaic from now onwards

task clean
reset
map=leop.combine.1pt8kms.map
beam=leop.combine.1pt8kms.beam
model=
out=leop.combine.1pt8kms.clean 
cutoff=0.0023 # 2.5 sigma
niters=1000000
phat=0.5
go

# Now INVERT

task restor
reset
model=leop.combine.1pt8kms.clean
beam=leop.combine.1pt8kms.beam
map=leop.combine.1pt8kms.map
out=leop.combine.1pt8kms.restor
go




###################################################################################



# Beam size is about 35" x 35"...this is much larger than we got from CASA
task invert
reset
vis=GMRT,VLAd,VLAc
map=leop.combine.1pt8kms.r2.map
beam=leop.combine.1pt8kms.r2.beam
imsize=512
cell=4
robust=2
line=velocity,70,221,1.8,1.8
inp
go

# mosaiced cube has higher noise level
# disregarding mosaic from now onwards

task clean
reset
map=leop.combine.1pt8kms.r2.map
beam=leop.combine.1pt8kms.r2.beam
model=
out=leop.combine.1pt8kms.r2.clean 
cutoff=0.0025 # 2.5 sigma
niters=1000000
phat=0.5
go

# Now INVERT

task restor
reset
model=leop.combine.1pt8kms.r2.clean
beam=leop.combine.1pt8kms.r2.beam
map=leop.combine.1pt8kms.r2.map
out=leop.combine.1pt8kms.r2.restor
go

###################################################################################
# Try without phat parameter
task clean
reset
map=leop.combine.1pt8kms.r2.map
beam=leop.combine.1pt8kms.r2.beam
model=
out=leop.combine.1pt8kms.r2.nophat.clean 
cutoff=0.0025 # 2.5 sigma
niters=1000000
phat=0
go

# Now INVERT

task restor
reset
model=leop.combine.1pt8kms.r2.nophat.clean
beam=leop.combine.1pt8kms.r2.beam
map=leop.combine.1pt8kms.r2.map
out=leop.combine.1pt8kms.r2.nophat.restor
go

# unsuccessful...

###################################################################################
# Beam is still large.
# Try setting GMRT weights to unity in CASA, and reloading

task fits
reset
in=/Users/research/data/leop/finalUVdataSets/leopGMRT.contsub.nowts.fits
op=uvin 
out=GMRT.nowts
inp
go

task invert
reset
vis=GMRT.nowts,VLAd,VLAc
map=leop.combine.1pt8kms.nowts.map
beam=leop.combine.1pt8kms.nowts.beam
imsize=512
cell=4
robust=2
line=velocity,70,221,1.8,1.8
inp
go

# mosaiced cube has higher noise level
# disregarding mosaic from now onwards

task clean
reset
map=leop.combine.1pt8kms.nowts.map
beam=leop.combine.1pt8kms.nowts.beam
model=
out=leop.combine.1pt8kms.nowts.clean 
cutoff=0.0025 # 2.5 sigma
niters=1000000
phat=0.5
go

# Now INVERT

task restor
reset
model=leop.combine.1pt8kms.nowts.clean
beam=leop.combine.1pt8kms.nowts.beam
map=leop.combine.1pt8kms.nowts.map
out=leop.combine.1pt8kms.nowts.restor
go

# convolve to a circular beam
task convol
reset
map=leop.combine.1pt8kms.nowts.restor
options=final
fwhm=26
pa=0
out=leop.combine.1pt8kms.nowts.restor.convol
go

###################################################################################
################
# Moment 0 Map #
################

# integrate
task moment
reset
in=leop.combine.1pt8kms.nowts.restor.convol
out=leop.combine.1pt8kms.nowts.mom0
mom=0
go

# create column density map
#	Recall:  1 K = (7.354E-8)*[Bmaj(")*Bmin(")/lamda^2(m)] Jy/Bm

#	Here, units of images are Jy/Bm m/s; cellsize = 4"; 	
#	    lambda = 0.211061140507 m

#	Thus, for the 21 cm line of Hydrogen, we have:

#	    1 K = Bmaj(")*Bmin(")/(6.057493205E5) Jy/Bm
#			---- OR ----
#	    1 Jy/Bm = (6.057493205E5)/[Bmaj(")*Bmin(")]

#	Now, recall that:
#		N_HI = (1.8224E18 cm^-2)*[T_b (K)]*int(dv)
#	   -- For moment maps in K km/sec, just input the 
#	      values & multiply by coefficient.
#	   -- Assure that units are Jy/Bm km/sec (i.e., divide by 1000)
#	      Leave in units of 1E20 cm^-2 by dividing by 1E20:

#	   For a x beam:
#		N_HI (cm^-2) = (image) * [(6.057493205E5)/(*)] * (1/1000) * (1.8224E18 cm^-2) * (1/1E20)
#		N_HI (cm^-2) = (image)*(0.01633)

#	Then put units of images as "1E20/cm2" using "puthead".

# Use maths
task maths
reset
exp='leop.combine.1pt8kms.nowts.mom0*0.01633*1000' # units are in mJy/Bm km/s
out=leop.combine.1pt8kms.nowts.colDens
go

task puthd
reset
in=leop.combine.1pt8kms.nowts.colDens/bunit
value='1e20.cm2'
go

# Blanking is way too hard in miriad.
# export to CASA for rest of image analysis

task fits
reset 
in=leop.combine.1pt8kms.nowts.restor.convol
out=leop.combine.1pt8kms.nowts.restor.convol.fits
op=xyout
go

task fits
reset
in=leop.combine.1pt8kms.nowts.colDens
out=leop.combine.1pt8kms.nowts.colDens.fits
op=xyout
go

# moved fits files to:
# /Users/research/data/leop/obsCombine/casa/successfulCombinedFiles

###################################################################################
###################################################################################
###################################################################################

# John calibrated the B-config obs...we should first test to see if it
# increases the quality of the cubes or degrades the quality

task fits 
reset
in=/Users/research/data/leop/finalUVdataSets/leopVLAb1.contsub.fits
op=uvin 
out=VLAb1
go

task fits 
reset
in=/Users/research/data/leop/finalUVdataSets/leopVLAb2.contsub.fits
op=uvin 
out=VLAb2
go

task invert
reset
vis=GMRT.nowts,VLAd,VLAc,VLAb1,VLAb2
map=leop.combine.bconfigTest.map
beam=leop.combine.bconfigTest.beam
imsize=512
cell=4
robust=2
line=velocity,70,221,1.8,1.8
inp
go

task clean
reset
map=leop.combine.bconfigTest.map
beam=leop.combine.bconfigTest.beam
model=
out=leop.combine.bconfigTest.clean 
cutoff=0.0018 # 3 sigma
niters=1000
phat=0.5
go

# Now INVERT

task restor
reset
model=leop.combine.bconfigTest.clean
beam=leop.combine.bconfigTest.beam
map=leop.combine.bconfigTest.map
out=leop.combine.bconfigTest.restor
go

# convolve to a circular beam
task convol
reset
map=leop.combine.bconfigTest.restor
options=final
fwhm=25.97,23.37
pa=-48.83
out=leop.combine.bconfigTest.convol
go

# RMSs (Jy/Bm):
# with B-config: 4.855E-04
# without B-config: 5.269E-04

# B-config beam size = 7.160 by    5.444 arcsec with robust of 2

###################################################################################
###################################################################################
###################################################################################

# Now attempting to make four cubes with 4", 8", 16" and 32" resolution

###################################################################################
###################################################################################
###################################################################################

# We will determine the correct fwhm taper value by trial and error:
task invert
reset
vis=GMRT.nowts,VLAd,VLAc,VLAb1,VLAb2
map=fwhm.14.map
beam=fwhm.14.beam
imsize=340
cell=6
robust=0.5
line=velocity,1,280,1.8,1.8
fwhm=14
options=mosaic
inp
go

task restor
reset
model=fwhm.14.map
beam=fwhm.14.beam
map=fwhm.14.map
out=fwhm.14.restor
go


# iteration | fwhm | cell ("/px) | imsize (px) | beam size (")
#         0   0.3    1.1           1900          3.52          
#         1   4      2             1050          11.06
#         2   2      2             1050          5.74
#         3   3      2             1050          8.55 # Use fwhm=2.5 for 8" cube
#         4   7      2             1050          16.19
#         5   6.8    4             1050          19.353
#         6   5      4             1050          16.51
#         7   6      4             1050          18.127 
#         8   11     8             510           26.71
#         9   9.5    8             510           25.3
#        10   8.5    8             510           24.353 # Use fwhm=8 for 24" cube
#        11   5      4             510           16.196 # Use fwhm=4.8 for 16" cube
#        12   14     6             340           29.06
#        13   15     6             340           29.918
#        14   19     6             340           32.980 use 18 for 32"


#----------------------------------------------------------------------------------
# Now we begin the cube imaging
#----------------------------------------------------------------------------------

# We also need the lowest resolution cube possible to probe the
# highest column densities over the Halpha region

#----------------
# 4" cube
#----------------

# Try setting the mosaic option
task invert
reset
vis=GMRT.nowts,VLAd,VLAc,VLAb1,VLAb2
map=leop.mosaic.4arcsec.map
beam=leop.mosaic.4arcsec.beam
imsize=1900
cell=1.1
robust=0.5
line=velocity,70,221,1.8,1.8
fwhm=0.3
options=mosaic
inp
go

task clean
reset
map=leop.mosaic.4arcsec.map
beam=leop.mosaic.4arcsec.beam
out=leop..mosaic.4arcsec.clean 
cutoff=0.0024 # 3 sigma
niters=1000
phat=0.5
go

task restor
reset
model=leop.mosaic.4arcsec.clean
beam=leop.mosaic.4arcsec.beam
map=leop.mosaic.4arcsec.map
out=leop.mosaic.4arcsec.restor
go

task convol
reset
map=leop.mosaic.4arcsec.restor
options=final
fwhm=4
region='boxes(800,760,1130,1130)'
pa=0
out=leop.mosaic.4arcsec.convol
go

# this does not work when running miriad through the start_miriad script
linmos in=leop.mosaic.4arcsec.restor out=leop.mosaic.4arcsec.convol.pbcor
# received the error:
### Fatal Error [linmos]:  Bad size for mosaic table
# cannot proceed with 4arcsec pbcor

task fits
reset
in=leop.mosaic.4arcsec.convol
op=xyout
out=leop.mosaic.4arcsec.cube.fits
go

task moment
reset
in=leop.mosaic.4arcsec.convol
out=leop.mosaic.4arcsec.mom0
region='images()'
mom=0
go

# Use maths for column density
task maths
reset
exp='leop.mosaic.4arcsec.mom0*0.68995*1000' # units are in mJy/Bm km/s
out=leop.mosaic.4arcsec.colDens
go

task fits
reset
in=leop.mosaic.4arcsec.colDens
op=xyout
out=leop.mosaic.4arcsec.colDens.fits
go

#----------------
# 8" cube
#----------------

task invert
reset
vis=GMRT.nowts,VLAd,VLAc,VLAb1,VLAb2
map=leop.mosaic.8arcsec.map
beam=leop.mosaic.8arcsec.beam
imsize=1050
cell=2
robust=0.5
line=velocity,70,221,1.8,1.8
fwhm=2.5
inp
go

task clean
reset
map=leop.mosaic.8arcsec.map
beam=leop.mosaic.8arcsec.beam
model=
out=leop.mosaic.8arcsec.clean 
cutoff=0.0015 # 3 sigma
niters=1000
phat=0.5
go

task restor
reset
model=leop.mosaic.8arcsec.clean
beam=leop.mosaic.8arcsec.beam
map=leop.mosaic.8arcsec.map
out=leop.mosaic.8arcsec.restor
go

task convol
reset
map=leop.mosaic.8arcsec.restor
options=final
fwhm=8
region='boxes(450,450,600,600)'
pa=0
out=leop.mosaic.8arcsec.convol
go

# this does not work when running miriad through the start_miriad script
linmos in=leop.mosaic.8arcsec.convol out=leop.mosaic.8arcsec.convol.pbcor

task fits
reset
in=leop.mosaic.8arcsec.convol.pbcor
op=xyout
out=leop.mosaic.8arcsec.cube.fits
go

task moment
reset
in=leop.mosaic.8arcsec.convol.pbcor
out=leop.mosaic.8arcsec.mom0
region='images()'
mom=0
go

# Use maths for column density
task maths
reset
exp='leop.mosaic.8arcsec.mom0*0.172487*1000' # units are in mJy/Bm km/s
out=leop.mosaic.8arcsec.colDens
go

task fits
reset
in=leop.mosaic.8arcsec.colDens
op=xyout
out=leop.mosaic.8arcsec.colDens.fits
go

#----------------
# 16" cube
#----------------

task invert
reset
vis=GMRT.nowts,VLAd,VLAc,VLAb1,VLAb2
map=leop.mosaic.16arcsec.map
beam=leop.mosaic.16arcsec.beam
imsize=512
cell=4
robust=0.5
line=velocity,70,221,1.8,1.8
fwhm=4.8
inp
go

task clean
reset
map=leop.mosaic.16arcsec.map
beam=leop.mosaic.16arcsec.beam
model=
out=leop.mosaic.16arcsec.clean 
cutoff=0.0018 # 3 sigma
niters=1000
phat=0.5
go

task restor
reset
model=leop.mosaic.16arcsec.clean
beam=leop.mosaic.16arcsec.beam
map=leop.mosaic.16arcsec.map
out=leop.mosaic.16arcsec.restor
go

task convol
reset
map=leop.mosaic.16arcsec.restor
options=final
fwhm=16
region='boxes(200,200,325,325)'
pa=0
out=leop.mosaic.16arcsec.convol
go

# this does not work when running miriad through the start_miriad script
linmos in=leop.mosaic.16arcsec.convol out=leop.mosaic.16arcsec.convol.pbcor

task fits
reset
in=leop.mosaic.16arcsec.convol.pbcor
op=xyout
out=leop.mosaic.16arcsec.cube.fits
go

task moment
reset
in=leop.mosaic.16arcsec.convol.pbcor
out=leop.mosaic.16arcsec.mom0
region='images(15,40)'
mom=0
go

# Use maths for column density
task maths
reset
exp='leop.mosaic.16arcsec.mom0*0.04312177*1000' # units are in mJy/Bm km/s
out=leop.mosaic.16arcsec.colDens
go

task fits
reset
in=leop.mosaic.16arcsec.colDens
op=xyout
out=leop.mosaic.16arcsec.colDens.fits
go

#----------------
# 32" cube
#----------------

task invert
reset
vis=GMRT.nowts,VLAd,VLAc,VLAb1,VLAb2
map=leop.mosaic.32arcsec.map
beam=leop.mosaic.32arcsec.beam
imsize=340
cell=6
robust=0.5
line=velocity,70,221,1.8,1.8
fwhm=15
inp
go

task clean
reset
map=leop.mosaic.32arcsec.map
beam=leop.mosaic.32arcsec.beam
model=
out=leop.mosaic.32arcsec.clean 
cutoff=0.0024 # 3 sigma
niters=1000
phat=0.5
go

task restor
reset
model=leop.mosaic.32arcsec.clean
beam=leop.mosaic.32arcsec.beam
map=leop.mosaic.32arcsec.map
out=leop.mosaic.32arcsec.restor
go

task convol
reset
map=leop.mosaic.32arcsec.restor
options=final
fwhm=32
region='boxes(140,140,200,200)'
pa=0
out=leop.mosaic.32arcsec.convol
go

task convol
reset
map=leop.mosaic.32arcsec.restor
options=final
fwhm=32
pa=0
out=leop.mosaic.32arcsec.fullcube.convol
go

# this does not work when running miriad through the start_miriad script
linmos in=leop.mosaic.32arcsec.convol out=leop.mosaic.32arcsec.convol.pbcor

task fits
reset
in=leop.mosaic.32arcsec.convol.pbcor
op=xyout
out=leop.mosaic.32arcsec.cube.fits
go

task fits
reset
in=leop.mosaic.32arcsec.fullcube.convol
op=xyout
out=leop.mosaic.32arcsec.fullcube.fits
go

task moment
reset
in=leop.mosaic.32arcsec.convol.pbcor
out=leop.mosaic.32arcsec.mom0
region='images(15,40)'
mom=0
go

task maths
reset
exp='leop.mosaic.32arcsec.mom0*0.01078044*1000' # units are in mJy/Bm km/s
out=leop.mosaic.32arcsec.colDens
go

task fits
reset
in=leop.mosaic.32arcsec.colDens
op=xyout
out=leop.mosaic.32arcsec.colDens.fits
go


#----------------
# 24" cube   # needed for Steve Warren to perform HI temperature analysis
#----------------

task invert
reset
vis=GMRT.nowts,VLAd,VLAc,VLAb1,VLAb2
map=leop.mosaic.24arcsec.map
beam=leop.mosaic.24arcsec.beam
imsize=512
cell=4
robust=0.5
line=velocity,70,221,1.8,1.8
fwhm=10
options=mosaic
inp
go

task clean
reset
map=leop.mosaic.24arcsec.map
beam=leop.mosaic.24arcsec.beam
out=leop.mosaic.24arcsec.clean 
cutoff=0.00297 # 3 sigma
niters=1000
phat=0.5
go

task restor
reset
model=leop.mosaic.24arcsec.clean
beam=leop.mosaic.24arcsec.beam
map=leop.mosaic.24arcsec.map
out=leop.mosaic.24arcsec.restor
go

task convol
reset
map=leop.mosaic.24arcsec.restor
options=final
fwhm=23.6
region='boxes(200,200,325,325)'
pa=0
out=leop.mosaic.24arcsec.convol
go

task fits
reset
in=leop.mosaic.24arcsec.convol
op=xyout
out=leop.mosaic.24arcsec.cube.fits
go

task moment
reset
in=leop.mosaic.24arcsec.convol
out=leop.mosaic.24arcsec.mom0
mom=0
go

task fits
reset
in=leop.mosaic.24arcsec.mom0
op=xyout
out=leop.mosaic.24arcsec.mom0.fits
go





# regrid 32" cube onto 16" cube for making channael maps

task regrid
resest
in=leop.mosaic.32arcsec.convol
tin=leop.mosaic.16arcsec.convol
out=leop.mosaic.32arcsec.convol.regrid
go

task fits
reset
in=leop.mosaic.32arcsec.convol.regrid
op=xyout
out=leop.mosaic.32arcsec.cube.regrid.fits
go



################################################################################
# The previous imaging yielded increasing flux with resolution. This is bad.
#
# 


# VLA + GMRT
invert vis=VLAd,VLAc,VLAb1,VLAb2,GMRTnormwts map=leop.mosaic.vla+gmrt.map beam=leop.mosaic.vla+gmrt.beam imsize=600 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm= options=mosaic

clean map=leop.mosaic.vla+gmrt.map beam=leop.mosaic.vla+gmrt.beam out=leop.mosaic.vla+gmrt.clean cutoff= niters=500 phat=0.5

restor model=leop.mosaic.vla+gmrt.clean beam=leop.mosaic.vla+gmrt.beam map=leop.mosaic.vla+gmrt.map out=leop.mosaic.vla+gmrt.restor

set name=vla+gmrt

convol map=leop.mosaic.$name.restor options=final fwhm=11 pa=0 out=leop.mosaic.$name.convol

moment in=leop.mosaic.$name.restor out=leop.mosaic.$name.mom0 mom=0 region='images(35,85)'

fits in=leop.mosaic.$name.mom0 out=leop.mosaic.$name.mom0.fits op=xyout

fits in=leop.mosaic.$name.restor out=leop.mosaic.$name.restor.fits op=xyout

# VLA only
invert vis=VLAd,VLAc,VLAb1,VLAb2 map=leop.mosaic.vla.map beam=leop.mosaic.vla.beam imsize=600 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm= options=mosaic

clean map=leop.mosaic.vla.map beam=leop.mosaic.vla.beam out=leop.mosaic.vla.clean cutoff= niters=500 phat=0.5

restor model=leop.mosaic.vla.clean beam=leop.mosaic.vla.beam map=leop.mosaic.vla.map out=leop.mosaic.vla.restor

set name=vla

convol map=leop.mosaic.$name.restor options=final fwhm=11 pa=0 out=leop.mosaic.$name.convol

moment in=leop.mosaic.$name.restor out=leop.mosaic.$name.mom0 mom=0 region='images(35,85)'

fits in=leop.mosaic.$name.mom0 out=leop.mosaic.$name.mom0.fits op=xyout

fits in=leop.mosaic.$name.restor out=leop.mosaic.$name.restor.fits op=xyout




###############################################################################
# VLA/d only
###############################################################################

##############
# Resolution 1
##############
set name=res1
set telescope=vla-d

invert vis=VLAd map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=700 cell=5 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=40

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.474E-03
set threshold=0.003685

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=67 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 2
##############
set name=res2

invert vis=VLAd map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=20

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.445E-03
set threshold=0.0036125000000000003

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=57 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 3
##############
set name=res3

invert vis=VLAd map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=700 cell=4 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=1

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.451E-03
set threshold=0.0036275

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=54 pa=0 out=leop.mosaic.$telescope.$name.convol



###############################################################################
# GMRT only
###############################################################################

##############
# Resolution 1
##############
set name=res1
set telescope=gmrt

invert vis=GMRTnormwts map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=700 cell=5 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=40

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.24E-03
set threshold=0.0031

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=48 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 2
##############
set name=res2

invert vis=GMRTnormwts map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=20

imstat in=leop.mosaic.$telescope.$name.map

set rms=9.971E-04
set threshold=0.0025

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=31 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 3
##############
set name=res3

invert vis=GMRTnormwts map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=2000 cell=1 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=1

imstat in=leop.mosaic.$telescope.$name.map

set rms=5.402E-04
set threshold=0.00135

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=4 pa=0 out=leop.mosaic.$telescope.$name.convol

###############################################################################
# VLA D/C only
###############################################################################

##############
# Resolution 1
##############
set name=res1
set telescope=vla-c+d

invert vis=VLAd,VLAc map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=700 cell=5 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=40 options=mosaic

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.242E-03
set threshold=0.003105

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=58 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 2
##############
set name=res2

invert vis=VLAd,VLAc map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=20

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.056E-03
set threshold=0.00264

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=36 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 3
##############
set name=res3

invert vis=VLAd,VLAc map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=1500 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=1

imstat in=leop.mosaic.$telescope.$name.map

set rms=9.570E-04
set threshold=0.00239

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=22 pa=0 out=leop.mosaic.$telescope.$name.convol


###############################################################################
# VLA D/C/B 
###############################################################################

##############
# Resolution 1
##############
set name=res1
set telescope=vla-b+c+d

invert vis=VLAd,VLAc,VLAb1,VLAb2 map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=1000 cell=3 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=40 options=mosaic

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.209E-03
set threshold=0.0030225

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=52 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 2
##############
set name=res2

invert vis=VLAd,VLAc,VLAb1,VLAb2 map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=20

imstat in=leop.mosaic.$telescope.$name.map

set rms=9.962E-04
set threshold=0.0024905

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=30 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 3
##############
set name=res3

invert vis=VLAd,VLAc,VLAb1,VLAb2 map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=1

imstat in=leop.mosaic.$telescope.$name.map

set rms=7.897E-04
set threshold=0.00197425

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=11 pa=0 out=leop.mosaic.$telescope.$name.convol


###############################################################################
# VLA D/C + GMRT
###############################################################################

##############
# Resolution 1
##############
set name=res1
set telescope=vla-c+d.gmrt

invert vis=VLAd,VLAc,GMRTnormwts map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=700 cell=5 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=40 options=mosaic

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.233E-03
set threshold=0.00308

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=58 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 2
##############
set name=res2

invert vis=VLAd,VLAc,GMRTnormwts map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=20

imstat in=leop.mosaic.$telescope.$name.map

set rms=7.165E-04
set threshold=0.00179125

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=31 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 3
##############
set name=res3

invert vis=VLAd,VLAc,GMRTnormwts map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=1

imstat in=leop.mosaic.$telescope.$name.map

set rms=4.922E-04
set threshold=0.0012305

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=9 pa=0 out=leop.mosaic.$telescope.$name.convol

foreach filename (*.convol)
  set base=`basename ${filename}`
  fits in=${filename} op=xyout out=../fitsImages/${base}.fits
end

foreach filename (*.restor)
  set base=`basename ${filename}`
  fits in=${filename} op=xyout out=${base}.fits
end

fits in=leop.mosaic.vla-c+d.res1.restor op=xyout out=../fitsImages/leop.mosaic.vla-c+d.res1.restor.fits

fits in=leop.mosaic.vla-c+d.res1.convol op=xyout out=../fitsImages/leop.mosaic.vla-c+d.res1.convol.fits

################################################################################
# Write out fits files
################################################################################

foreach filename (*.restor)
  set base=`basename ${filename}`
  fits in=${filename} op=xyout out=../fitsImages/${base}.fits
end

foreach filename (*.convol)
  set base=`basename ${filename}`
  fits in=${filename} op=xyout out=../fitsImages/${base}.fits
end



################################################################################
################################################################################
################################################################################

# John requested the following:

# Make new cubes with SAME GRID for ALL FOUR resolutions.  Suggest 1 or 1.5
# arcseconds per pixel.  Then you will have 4 cubes, each on the same grid.
# The same RA/DEC corresponds to the same x/y.  You will over-sample the coarse
# beam cubes (e.g., the 32" cube will have hundreds of pixels per beam); this
# is ok, so long as you do not undersample the beam (e.g., the native beam size
# of the unconvolved 4" cube has < 10 pixels per beam; see above).

# Then, concolve the 32" image to a larger beam size.  Blank that cube.  Then
# use that cube as a mastr balnking mask against the other 4 cubes.  Blank both
# the original beam shape and the convolved beam shape cubes.  This way you
# have the same regions contributing at all resolutions (see Section 3.6 of
# Walter et al. 2008).  Then measure the fluxes in the resulting moment maps.  

################################################################################
################################################################################
################################################################################

# Make combined visibility dataset
uvaver vis=VLAc,VLAd,GMRT line=velocity,100,190,1.8,1.8 out=VLA+GMRT

##############
# Resolution 1
##############
set name=res1
set telescope=vla-c+d.gmrt.grid

invert vis=VLAc,VLAd,GMRT map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=2040 cell=1 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=40 options=mosaic

imstat in=leop.mosaic.$telescope.$name.map

set rms=2.137E-03
set threshold=0.0053425

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=46 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 2
##############
set name=res2
set telescope=vla-c+d.gmrt.grid

invert vis=VLAc,VLAd,GMRT map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=2040 cell=1 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=20 options=mosaic

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.433E-03
set threshold=0.0035825

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=26 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 3
##############
set name=res3
set telescope=vla-c+d.gmrt.grid

invert vis=VLAc,VLAd,GMRT map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=2040 cell=1 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=5 options=mosaic

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.182E-03
set threshold=0.002955

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=16 pa=0 out=leop.mosaic.$telescope.$name.convol

##############
# Resolution 4
##############
set name=res4
set telescope=vla-c+d.gmrt.grid

invert vis=VLAc,VLAd,GMRT map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam imsize=2040 cell=1 robust=0.5 line=velocity,100,190,1.8,1.8 fwhm=0.05 options=mosaic

imstat in=leop.mosaic.$telescope.$name.map

set rms=1.179E-03
set threshold=0.0029474999999999996

clean map=leop.mosaic.$telescope.$name.map beam=leop.mosaic.$telescope.$name.beam out=leop.mosaic.$telescope.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.mosaic.$telescope.$name.clean beam=leop.mosaic.$telescope.$name.beam map=leop.mosaic.$telescope.$name.map out=leop.mosaic.$telescope.$name.restor

convol map=leop.mosaic.$telescope.$name.restor options=final fwhm=16 pa=0 out=leop.mosaic.$telescope.$name.convol

foreach filename (*.restor)
  set base=`basename ${filename}`
  fits in=${filename} op=xyout out=../fitsImages/${base}.fits
end

foreach filename (*.convol)
  set base=`basename ${filename}`
  fits in=${filename} op=xyout out=../fitsImages/${base}.fits
end




################################################################################
################################################################################
################################################################################
# Perform UVaver on each dataset to see how the Miriad grids 
################################################################################
################################################################################
################################################################################

task uvaver
reset
vis=VLAc
line=velocity,70,185,2.5,2.5
out=VLAc_uvaver
go

task uvaver
reset
vis=VLAd
line=velocity,70,185,2.5,2.5
out=VLAd_uvaver
go

task uvaver
reset
vis=GMRT
line=velocity,70,185,2.5,2.5
out=GMRT_uvaver
go

set vis=VLAc
uvlist vis=$vis options=spectral
set vis=VLAc_uvaver
uvlist vis=$vis options=spectral

set vis=VLAd
uvlist vis=$vis options=spectral
set vis=VLAd_uvaver
uvlist vis=$vis options=spectral

set vis=GMRT
uvlist vis=$vis options=spectral
set vis=GMRT_uvaver
uvlist vis=$vis options=spectral



set data=VLAc_uvaver
set name=vlac

invert vis=$data map=leop.uvaver.$name.map beam=leop.uvaver.$name.beam imsize=1000 cell=2 robust=0.5 fwhm=20

imstat in=leop.uvaver.$name.map

set rms=1.355E-03
set threshold=0.0033875

clean map=leop.uvaver.$name.map beam=leop.uvaver.$name.beam out=leop.uvaver.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.uvaver.$name.clean beam=leop.uvaver.$name.beam map=leop.uvaver.$name.map out=leop.uvaver.$name.restor


set data=VLAd_uvaver
set name=vlad

invert vis=$data map=leop.uvaver.$name.map beam=leop.uvaver.$name.beam imsize=1000 cell=1 robust=0.5 fwhm=20

imstat in=leop.uvaver.$name.map

set rms=1.310E-03
set threshold=0.003275

clean map=leop.uvaver.$name.map beam=leop.uvaver.$name.beam out=leop.uvaver.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.uvaver.$name.clean beam=leop.uvaver.$name.beam map=leop.uvaver.$name.map out=leop.uvaver.$name.restor


set data=GMRT_uvaver
set name=gmrt

invert vis=$data map=leop.uvaver.$name.map beam=leop.uvaver.$name.beam imsize=1000 cell=1 robust=0.5 fwhm=20

imstat in=leop.uvaver.$name.map

set rms=9.117E-04
set threshold=0.00227925

clean map=leop.uvaver.$name.map beam=leop.uvaver.$name.beam out=leop.uvaver.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.uvaver.$name.clean beam=leop.uvaver.$name.beam map=leop.uvaver.$name.map out=leop.uvaver.$name.restor



################################################################################
################################################################################
################################################################################
# Load UV datasets keeping track of doppler shifting
################################################################################
################################################################################
################################################################################


cd /d/bip3/ezbc/leop/data/hi/miriad/images
set uvDir=/d/bip3/ezbc/leop/data/hi/finalUVdatasets

# load the datasets in:
# GMRT
set uvData=leopGMRT.contsub.normwts.fits
set uvOut=GMRTnormwts.doppler
fits in=$uvDir/$uvData op=uvin out=$uvOut velocity=lsr

# VLA/D
set uvData=leopVLAd.contsub.fits
set uvOut=VLAd.doppler
fits in=$uvDir/$uvData op=uvin out=$uvOut velocity=lsr

# VLA/C
set uvData=leopVLAc.contsub.fits
set uvOut=VLAc.doppler
fits in=$uvDir/$uvData op=uvin out=$uvOut velocity=lsr

# VLA/B-1
set uvData=leopVLAb1.contsub.fits
set uvOut=VLAb1.doppler
fits in=$uvDir/$uvData op=uvin out=$uvOut velocity=lsr

# VLA/B-2
set uvData=leopVLAb2.contsub.fits
set uvOut=VLAb2.doppler
fits in=$uvDir/$uvData op=uvin out=$uvOut velocity=lsr



################################################################################
# Now Image VLA/C, VLA/D, and GMRT separately

#######
# GMRT
#######
set data=GMRTnormwts.doppler
set name=gmrt

invert vis=$data map=leop.doppler.$name.map beam=leop.doppler.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,70,185,2.5,2.5

imstat in=leop.doppler.$name.map

set rms=5.292E-04
set threshold=0.0015875999999999998

clean map=leop.doppler.$name.map beam=leop.doppler.$name.beam out=leop.doppler.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.doppler.$name.clean beam=leop.doppler.$name.beam map=leop.doppler.$name.map out=leop.doppler.$name.restor

#######
# VLAd
#######
set data=VLAd.doppler
set name=vlad

invert vis=$data map=leop.doppler.$name.map beam=leop.doppler.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,70,185,2.5,2.5

imstat in=leop.doppler.$name.map

set rms=1.279E-03
set threshold=0.003837

clean map=leop.doppler.$name.map beam=leop.doppler.$name.beam out=leop.doppler.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.doppler.$name.clean beam=leop.doppler.$name.beam map=leop.doppler.$name.map out=leop.doppler.$name.restor

#######
# VLAc
#######
set data=VLAc.doppler
set name=vlac

invert vis=$data map=leop.doppler.$name.map beam=leop.doppler.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,70,185,2.5,2.5

imstat in=leop.doppler.$name.map

set rms=1.079E-03
set threshold=0.0032370000000000

clean map=leop.doppler.$name.map beam=leop.doppler.$name.beam out=leop.doppler.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.doppler.$name.clean beam=leop.doppler.$name.beam map=leop.doppler.$name.map out=leop.doppler.$name.restor


# The GMRT doppler-shifted data seems to align with the non-doppler-shifted VLA
# data. Rats. Image the VLA/D and VLA/C without doppler corrections.

# See the following link for a description of doppler tracking for the LVA
# https://science.nrao.edu/facilities/vla/docs/manuals/obsguide/modes/line/velocity
#
# The VLA does not support doppler tracking.
#
# Used task Cvel in CASA to put the three datasets in the LSRK frame

cd /d/bip3/ezbc/leop/data/hi/miriad/images
set uvDir=/d/bip3/ezbc/leop/data/hi/finalUVdatasets

# load the datasets in:
# GMRT
set uvData=leopGMRT.contsub.normwts.fits
set uvOut=GMRTnormwts.cvel
fits in=$uvDir/$uvData op=uvin out=$uvOut

# VLA/D
set uvData=leopVLAd.contsub.fits
set uvOut=VLAd.cvel
fits in=$uvDir/$uvData op=uvin out=$uvOut

# VLA/C
set uvData=leopVLAc.contsub.fits
set uvOut=VLAc.cvel
fits in=$uvDir/$uvData op=uvin out=$uvOut

################################################################################
# Now image

#######
# GMRT
#######

set data=GMRTnormwts.cvel
set name=gmrt

invert vis=$data map=leop.cvel.$name.map beam=leop.cvel.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,70,185,2.5,2.5 fwhm=40

imstat in=leop.cvel.$name.map

set rms=1.173E-03
set threshold=0.003519

clean map=leop.cvel.$name.map beam=leop.cvel.$name.beam out=leop.cvel.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.cvel.$name.clean beam=leop.cvel.$name.beam map=leop.cvel.$name.map out=leop.cvel.$name.restor

#######
# VLAd
#######
set data=VLAd.cvel
set name=vlad

invert vis=$data map=leop.cvel.$name.map beam=leop.cvel.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,70,185,2.5,2.5 

imstat in=leop.cvel.$name.map

set rms=1.280E-03
set threshold=0.0038400

clean map=leop.cvel.$name.map beam=leop.cvel.$name.beam out=leop.cvel.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.cvel.$name.clean beam=leop.cvel.$name.beam map=leop.cvel.$name.map out=leop.cvel.$name.restor

#######
# VLAc
#######
set data=VLAc.cvel
set name=vlac

invert vis=$data map=leop.cvel.$name.map beam=leop.cvel.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,70,185,2.5,2.5 fwhm=40

imstat in=leop.cvel.$name.map

set rms=2.040E-03
set threshold=0.00612

clean map=leop.cvel.$name.map beam=leop.cvel.$name.beam out=leop.cvel.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.cvel.$name.clean beam=leop.cvel.$name.beam map=leop.cvel.$name.map out=leop.cvel.$name.restor

# GMRT data is offset again, they switched peak flux velocities! Wtf?
#
# Try imaging the GMRT doppler-shifted data again.

#######
# GMRT
#######

set data=GMRTnormwts.doppler
set name=gmrt

invert vis=$data map=leop.cvel.doppler.$name.map beam=leop.cvel.doppler.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,70,185,2.5,2.5 fwhm=40

imstat in=leop.cvel.doppler.$name.map

set rms=1.173E-03
set threshold=0.003519

clean map=leop.cvel.doppler.$name.map beam=leop.cvel.doppler.$name.beam out=leop.cvel.doppler.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.cvel.doppler.$name.clean beam=leop.cvel.doppler.$name.beam map=leop.cvel.doppler.$name.map out=leop.cvel.doppler.$name.restor

# Works!

# Write out fits images
set name=vlac
fits in=leop.cvel.$name.restor op=xyout out=../fitsImages/leop.cvel.$name.restor.fits

set name=vlad
fits in=leop.cvel.$name.restor op=xyout out=../fitsImages/leop.cvel.$name.restor.fits

set name=doppler.gmrt
fits in=leop.cvel.$name.restor op=xyout out=../fitsImages/leop.cvel.$name.restor.fits




# Image original GMRT data

cd /d/bip3/ezbc/leop/data/hi/miriad/images
set uvDir=/d/bip3/ezbc/leop/data/hi/finalUVdatasets

# load the datasets in:
# GMRT
set uvData=leopGMRT.contsub.fits
set uvOut=GMRT.original
fits in=$uvDir/$uvData op=uvin out=$uvOut


set data=GMRT.original
set name=gmrt.original

invert vis=$data map=leop.$name.map beam=leop.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,70,185,2.5,2.5 

imstat in=leop.$name.map

set rms=5.293E-04
set threshold=0.00158790

clean map=leop.$name.map beam=leop.$name.beam out=leop.$name.clean cutoff=$threshold niters=10000 phat=0.5

restor model=leop.$name.clean beam=leop.$name.beam map=leop.$name.map out=leop.$name.restor

fits in=leop.$name.restor op=xyout out=../fitsImages/leop.$name.restor.fits

# Nope, giant flux
# Moment Image: leop.gmrt.original.4arcsec.mom0.blk.image
# Beamsize: 4.08997614772" X 3.6035479816"
# Flux: 2.16258864594 Jy km/s


################################################################################# Image VLA/B data to see if we see anything

set data=VLAb1,VLAb2
set name=VLAb

invert vis=$data map=leop.$name.map beam=leop.$name.beam imsize=1000 cell=2 robust=0.5 line=velocity,70,185,2.5,2.5 options=mosaic fwhm=20

fits in=leop.$name.map op=xyout out=../fitsImages/leop.$name.fits


################################################################################
################################################################################
################################################################################
# Imaging VLA/C + VLA/D + GMRT Once again, this time with a new method
################################################################################
################################################################################
################################################################################

# Steps for imaging heterogeneous interferometric observations in Miriad
# a) Image the three datasets together using invert. 
# b) Derive a combined primary beam PSF of the beam from a) using mospsf. 
# c) Create a model PSF from b) using imfit. 
# d) Use mossdi2 to deconvolve the inverted datasets using the output beam
#    from c).
# e) Restore the output model image from d) with restor.

cd /d/bip3/ezbc/leop/data/hi/miriad/images
set uvDir=/d/bip3/ezbc/leop/data/hi/finalUVdatasets

# load the datasets in:
# GMRT
set uvData=leopGMRT.contsub.normwts.fits
set uvOut=GMRTnormwts.doppler
fits in=$uvDir/$uvData op=uvin out=$uvOut velocity=lsr

# VLA/D
set uvData=leopVLAd.contsub.fits
set uvOut=VLAd.cvel
fits in=$uvDir/$uvData op=uvin out=$uvOut

# VLA/C
set uvData=leopVLAc.contsub.fits
set uvOut=VLAc.cvel
fits in=$uvDir/$uvData op=uvin out=$uvOut

set data=../images/GMRTnormwts.doppler,../images/VLAc.cvel,../images/VLAd.cvel
set name=mosaic.lsr
set res=res0

################################################################################
# Determine the needed FHWMs of Gaussian UV tapers
#

# See folder /d/bip3/ezbc/leop/data/hi/miriad/revisedImages/tests, and
# fwhmTests.csh for details

fwhm=0.1

invert vis=$data map=test.$fwhm.map beam=test.$fwhm.beam imsize=2000 cell=1 robust=0.5 line=velocity,1,270,1.8,1.8 options=mosaic fwhm=$fwhm 
restor map=test.$fwhm.map beam=test.$fwhm.beam model=test.$fwhm.map out=test.$fwhm.restor

# Fwhm | beamsize (")
# 0.8    3.989 x 3.335
# 4      7.995 x 7.489
# 9      15.999 x 14.067
# 22     31.940 x 30.643


##################################
# Without primary beam correction

set fwhm=0.8
set name=mosaic.lsr
set res=4arcsec
invert vis=$data map=leop.$name.$res.map beam=leop.$name.$res.beam imsize=2000 cell=1 robust=0.5 line=velocity,100,185,1.8,1.8 options=mosaic fwhm=$fwhm 

imstat in=leop.$name.$res.map

set cutoff=0.002367

mossdi2 map=leop.$name.$res.map beam=leop.$name.$res.beam out=leop.$name.$res.mossdi2 cutoff=$cutoff niters=5000

restor map=leop.$name.$res.map beam=leop.$name.$res.beam model=leop.$name.$res.mossdi2 out=leop.$name.$res.restor

restor map=leop.$name.$res.map beam=leop.$name.$res.beam model=leop.$name.$res.mossdi2 out=leop.$name.$res.restor.residual mode=residual


set fwhm=4
set name=mosaic.lsr
set res=8arcsec
invert vis=$data map=leop.$name.$res.map beam=leop.$name.$res.beam imsize=2000 cell=1 robust=0.5 line=velocity,100,185,1.8,1.8 options=mosaic fwhm=$fwhm 

imstat in=leop.$name.$res.map

set cutoff=0.002369

mossdi2 map=leop.$name.$res.map beam=leop.$name.$res.beam out=leop.$name.$res.mossdi2 cutoff=$cutoff niters=5000

restor map=leop.$name.$res.map beam=leop.$name.$res.beam model=leop.$name.$res.mossdi2 out=leop.$name.$res.restor

restor map=leop.$name.$res.map beam=leop.$name.$res.beam model=leop.$name.$res.mossdi2 out=leop.$name.$res.restor.residual mode=residual


set fwhm=9
set name=mosaic.lsr
set res=16arcsec
invert vis=$data map=leop.$name.$res.map beam=leop.$name.$res.beam imsize=2000 cell=1 robust=0.5 line=velocity,100,185,1.8,1.8 options=mosaic fwhm=$fwhm 

imstat in=leop.$name.$res.map

set cutoff=0.00240075

mossdi2 map=leop.$name.$res.map beam=leop.$name.$res.beam out=leop.$name.$res.mossdi2 cutoff=$cutoff niters=5000

restor map=leop.$name.$res.map beam=leop.$name.$res.beam model=leop.$name.$res.mossdi2 out=leop.$name.$res.restor

restor map=leop.$name.$res.map beam=leop.$name.$res.beam model=leop.$name.$res.mossdi2 out=leop.$name.$res.restor.residual mode=residual

set fwhm=22
set name=mosaic.lsr
set res=32arcsec
invert vis=$data map=leop.$name.$res.map beam=leop.$name.$res.beam imsize=2000 cell=1 robust=0.5 line=velocity,100,185,1.8,1.8 options=mosaic fwhm=$fwhm 

imstat in=leop.$name.$res.map

set cutoff=0.002669

mossdi2 map=leop.$name.$res.map beam=leop.$name.$res.beam out=leop.$name.$res.mossdi2 cutoff=$cutoff niters=5000

restor map=leop.$name.$res.map beam=leop.$name.$res.beam model=leop.$name.$res.mossdi2 out=leop.$name.$res.restor

restor map=leop.$name.$res.map beam=leop.$name.$res.beam model=leop.$name.$res.mossdi2 out=leop.$name.$res.restor.residual mode=residual

# to perform flux rescale we need to perform the following: G = D x C / (D - R)
# where G is true flux, D is dirty image, C is clean image, and R is residual
# image.

set name=mosaic.lsr
set res=32arcsec
fits in=leop.$name.$res.restor.residual op=xyout out=leop.$name.$res.residual.fits 
fits in=leop.$name.$res.restor op=xyout out=leop.$name.$res.clean.fits 
fits in=leop.$name.$res.map op=xyout out=leop.$name.$res.dirty.fits 
fits in=leop.$name.$res.beam op=xyout out=leop.$name.$res.psf.fits 


set name=mosaic.lsr
set res=16arcsec
fits in=leop.$name.$res.restor.residual op=xyout out=leop.$name.$res.residual.fits 
fits in=leop.$name.$res.restor op=xyout out=leop.$name.$res.clean.fits 
fits in=leop.$name.$res.map op=xyout out=leop.$name.$res.dirty.fits 
fits in=leop.$name.$res.beam op=xyout out=leop.$name.$res.psf.fits 


set name=mosaic.lsr
set res=8arcsec
fits in=leop.$name.$res.restor.residual op=xyout out=leop.$name.$res.residual.fits 
fits in=leop.$name.$res.restor op=xyout out=leop.$name.$res.clean.fits 
fits in=leop.$name.$res.map op=xyout out=leop.$name.$res.dirty.fits 
fits in=leop.$name.$res.beam op=xyout out=leop.$name.$res.psf.fits 


set name=mosaic.lsr
set res=4arcsec
fits in=leop.$name.$res.restor.residual op=xyout out=leop.$name.$res.residual.fits 
fits in=leop.$name.$res.restor op=xyout out=leop.$name.$res.clean.fits 
fits in=leop.$name.$res.map op=xyout out=leop.$name.$res.dirty.fits 
fits in=leop.$name.$res.beam op=xyout out=leop.$name.$res.psf.fits 




##################################
# With primary beam correction
set fwhm=0.8
set name=mosaic.lsr
set res=4arcsec
set data=../images/GMRTnormwts.doppler,../images/VLAc.cvel,../images/VLAd.cvel

invert vis=$data map=leop.$name.$res.map beam=leop.$name.$res.beam imsize=2000 cell=1 robust=0.5 line=velocity,100,185,1.8,1.8 options=mosaic fwhm=0.8 

mospsf beam=leop.$name.$res.beam out=leop.$name.$res.mospsf

imfit in=leop.$name.$res.mospsf object=beam > beamlog.int

set maj=`grep 'Major axis' beamlog.int | cut -c32-38`
set min=`grep 'Minor axis' beamlog.int | cut -c33-38`
echo $maj $min

imstat in=leop.$name.$res.map

set cutoff=0.002425999999

mossdi2 map=leop.$name.$res.map beam=leop.$name.$res.mospsf out=leop.$name.$res.mossdi2 cutoff=$cutoff niters=5000

restor map=leop.$name.$res.map beam=leop.$name.$res.beam model=leop.$name.$res.mossdi2 out=leop.$name.$res.restor











