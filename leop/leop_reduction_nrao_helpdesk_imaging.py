#!/usr/bin/python

################################################################################
# Image individual observations
################################################################################

#########
# VLA/C
#########
# image the flux calibrator
default(clean)
vis = vlaC
field = '1331+305=3C286'
mode = 'mfs'
imagename = 'leopVLAc.fluxCalibrator'
niter = 1000 # for light cleaning
threshold = '2mJy' # see note1
imsize = 512
cell = '2.8arcsec' # see note2
weighting = 'briggs'
robust = 0.5
clean()

# image the source continuum
default(clean)
vis = vlaC
imagename = 'leopVLAc.continuum'
mode = 'velocity'
restfreq = '1420.405752MHz'
niter = 0
imsize = 1024
cell = '2arcsec'
weighting = 'briggs'
robust = 0.5
clean()

#########
# VLA/D
#########

# image the flux calibrator
default(clean)
vis = vlaD
field = '1331+305=3C286'
mode = 'mfs'
imagename = 'leopVLAd.fluxCalibrator'
niter = 1000 # for light cleaning
threshold = '3mJy' # see note1
imsize = 512
cell = '9arcsec' # see note2
weighting = 'briggs'
robust = 0.5
clean()

# image the source continuum
default(clean)
vis = vlaD
imagename = 'leopVLAd.continuum'
mode = 'velocity'
restfreq = '1420.405752MHz'
niter = 0
imsize = 240
cell = '9arcsec'
weighting = 'briggs'
robust = 0.5
clean()

################################################################################
# Image observations together ################################################################################














