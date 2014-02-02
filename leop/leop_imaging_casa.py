#!/usr/bin/python
################################################################################
################################################################################
# Script to image GMRT + VLA/B/C/D observations in CASA
################################################################################
################################################################################

# Define uvdataset directory
uvDir = '/d/bip3/ezbc/leop/data/hi/finalUVdataSets/'

# First import uvdatasets
importgmrt(fitsfile=uvDir + 'leopGMRT.contsub.fits',
            vis='leopGMRT.contsub.ms')

importgmrt(fitsfile=uvDir + 'leopGMRT.contsub.nowts.fits',
            vis='leopGMRT.contsub.nowts.ms')

importuvfits(fitsfile=uvDir + 'leopVLAb1.contsub.fits',
             vis='leopVLAb1.contsub.ms')

importuvfits(fitsfile=uvDir + 'leopVLAb2.contsub.fits',
             vis='leopVLAb2.contsub.ms')

importuvfits(fitsfile=uvDir + 'leopVLAc.contsub.fits',
             vis='leopVLAc.contsub.ms')

importuvfits(fitsfile=uvDir + 'leopVLAd.contsub.fits',
             vis='leopVLAd.contsub.ms')


################################################################################
# Normalize GMRT weights
################################################################################

# First create GMRT copy


# Import numpy, which is handy for this sort of manipulation
import numpy as np

# Open the MS.  Note that this will have the ability to modify the values in the MS, so it might be a good idea to make a copy before proceeding in case something goes amiss.
ms.open('leopGMRT.contsub.ms', nomodify = False)

# Select your SPW / datadesc (see the DATA_DESCRIPTION table in the MS to relate the SPW ID and datadesc ID, which is the corresponding row index)
ms.selectinit(datadescid=0)

# Iterate in time order
ms.iterinit(columns='TIME')
ms.nrow()
ms.iterorigin()

# Note that you might want to play around with an ms.getdata() command at this point before going through and iterating / changing values

hasdata=True
while (hasdata):
# Gets the data and puts it in a dictionary called 'dat'
    dat = ms.getdata(['weight'])
# Do whatever you like with the manipulation here
    maxVal = max(dat['weight'][0])
    dat['weight'] = dat['weight'] / maxVal
# Puts the modified data back into the MS and moves to the next chunk
    ms.putdata(dat)
    hasdata=ms.iternext()

ms.close()




# Now try to clean the three uvdatasets together
# First we need to understand which uvtapers will return the desired beam sizes

################################################################################
# Now image VLA, GMRT, and VLA+GMRT separately
################################################################################

# define ms filenames
vlab1 = 'leopVLAb1.contsub.ms'
vlab2 = 'leopVLAb2.contsub.ms'
vlac = 'leopVLAc.contsub.ms'
vlad = 'leopVLAd.contsub.ms'
gmrt = 'leopGMRT.contsub.normwts.ms'

# VLA only
clean(vis=[vlab1,vlab2,vlac,vlad],
      imagename='images/leop.mosaic.vla.3',
      imagermode='mosaic',
      mode='velocity',
      outframe='lsrk',
      start='200km/s',
      width='1.8km/s',
      nchan=100,
      interpolation='nearest',
      imsize=512,
      cell='2arcsec',
      niter=1000,
      restfreq='1.420405752GHz',
      outertaper='2arcsec')

clean(vis=[gmrt],
      imagename='images/leop.mosaic.gmrt',
      mode='velocity',
      outframe='lsrk',
      start='240km/s',
      width='1.8km/s',
      nchan=20,
      interpolation='nearest',
      imsize=512,
      cell='2arcsec',
      niter=1000,
      restfreq='1.420405752GHz')

clean(vis=[vlab1,vlab2,vlac,vlad,gmrt],
      imagename='images/leop.mosaic.vla+gmrt',
      mode='velocity',
      outframe='lsrk',
      start='240km/s',
      width='1.8km/s',
      nchan=20,
      interpolation='nearest',
      imsize=512,
      cell='2arcsec',
      niter=1000,
      restfreq='1.420405752GHz')


################################################################################
# Now image VLA at different resolutions
################################################################################

import numpy as np

# concatenate VLA ms files together

concat(vis=[vlab1,vlab2,vlac,vlad],
        concatvis='leopVLAbcd.ms',
        dirtol='1arcsec')

uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/'

vla1 = uvdataDir + 'leopVLAbcd.ms'
vla2 = [uvdataDir + vlab1,uvdataDir + vlab2,uvdataDir + vlac,uvdataDir + vlad]


###################
# Res 1
###################

# Make dirty cube
name='test.res1'
clean(vis=uvdataDir + vla,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=4,
      imsize=1024,
      cell='1.5arcsec',
      niter=1000,
      threshold='1mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='50arcsec',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(imagesDir + 'leop.mosaic.vla.' + name + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value'] # arcsec
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 4.
imagesize = 2 * int(34*60. / cellsize) # images 2X entire primary beam (PB)

if imagesize%100 != 0: imagesize = imagesize + 100 - imagesize%100

# define multiscale 
multiscale = np.array([0.,1.,5.]) * beamsize / cellsize

# Clean full cube with concatenated visibilities
name='res1.noconcat'
clean(vis=vla1,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=2*imagesize,
      cell=cellsize,
      niter=1000,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      outertaper='50arcsec',
      multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# Clean full cube, testing without concatenated visibilities
name='res1.concat'
clean(vis=vla2,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=imagesize,
      cell=cellsize,
      niter=1000,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      outertaper='50arcsec',
      multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# this did not work!!!!

### NOTE: In order to avoid aliasing artifacts for ftmachine=’mosaic’ in the mosaic image, due to the discrete sampling of the mosaic pattern on the sky, you should make an image in which the desired unmasked part of the image (above minpb) lies within the inner quarter. In other words, make an image twice as big as necessary to encompass the mosaic.

### NOTE: Sault weighting implies a uniform noise mosaic  

### NOTE: that niter is set to large number so that stopping point is  
### controlled by threshold.

### NOTE: with pbcor=False, the final image is not "flux correct",  
### instead the image has constant noise despite roll off in power as  
### you move out from the phase center(s). Though this format makes it  
### "look nicest", for all flux density measurements, and to get an  
### accurate integrated intensity image, one needs to divide the  
### image.image/image.flux in order to correct for the mosaic  
### response pattern. One could also achieve this by setting pbcor=True  
### in clean.  


# Just perform FT to see if source is present, testing with concatenated
# visibilities
name='test1'
clean(vis=vla2,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=imagesize,
      cell=cellsize,
      niter=0,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      outertaper='50arcsec',
      multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################################
# Image only VLA D, no cleaning, multiscale=False
###################################

# Just perform FT to see if source is present, testing with only vlad
name='test2'
clean(vis=uvdataDir + vlad,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      #imagermode='mosaic',
      #ftmachine='mosaic',
      #mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=500,
      cell=10,
      niter=0,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='50arcsec',
      #multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################################
# Image only VLA D, no cleaning, multiscale=True
###################################

# Just perform FT to see if source is present, testing with only vlad
name='test3'
clean(vis=uvdataDir + vlad,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      #imagermode='mosaic',
      #ftmachine='mosaic',
      #mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=500,
      cell=10,
      niter=0,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='50arcsec',
      multiscale=[0,6,12],
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################################
# Image only VLA C/D, no cleaning
###################################

# Just perform FT to see if source is present, testing with only vlad
name='test4'
clean(vis=[uvdataDir + vlad,uvdataDir + vlac],
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='ft',
      mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=500,
      cell=10,
      niter=0,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='50arcsec',
      multiscale=[0,6,12],
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

######################
# Image only VLA C/D
######################

name='test5'
clean(vis=[uvdataDir + vlad,uvdataDir + vlac],
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=500,
      cell=10,
      niter=0,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='50arcsec',
      multiscale=[0,6,12],
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

#######################################
# Image VLA B/C/D with mosweight=True
#######################################

name='test6'
clean(vis=vla2,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=500,
      cell=10,
      niter=0,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='50arcsec',
      multiscale=[0,6,12],
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

########################################
# Image VLA B/C/D with mosweight=False
########################################

name='test7'
clean(vis=vla2,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=False,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=500,
      cell=10,
      niter=0,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='50arcsec',
      multiscale=[0,6,12],
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

####################
# Image only VLA/B
####################

name='test8'
clean(vis=[uvdataDir + vlab1,uvdataDir + vlab2],
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=False,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='270km/s',
      width='1.8km/s',
      nchan=10,
      imsize=500,
      cell=10,
      niter=0,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='50arcsec',
      multiscale=[0,6,12],
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# Clean full cube, testing without concatenated visibilities

# proceeding only with vla c/d observations

name='test.res1'
clean(vis=[uvdataDir + vlad,uvdataDir + vlac],
      imagename=imagesDir + 'leop.mosaic.vla-c+d.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=10,
      imsize=1000,
      cell=2,
      niter=0,
      #threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='50arcsec',
      #multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(imagesDir + 'leop.mosaic.vla-c+d.' + name + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value'] # arcsec
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 4.
imagesize = 2 * int(34*60. / cellsize) # images 2X entire primary beam (PB)

if imagesize%100 != 0: imagesize = imagesize + 100 - imagesize%100

# define multiscale 
multiscale = np.array([0.,1.,5.]) * beamsize / cellsize


name='res1'
clean(vis=[uvdataDir + vlad,uvdataDir + vlac],
      imagename=imagesDir + 'leop.mosaic.vla-c+d.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=100,
      imsize=imagesize,
      cell=cellsize,
      niter=10000,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='50arcsec',
      multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')





################################################################################
# Image gmrt at different resolution
################################################################################

################
# Resolution 1 #
################

name='test.res1'
clean(vis=[uvdataDir + gmrt],
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      #imagermode='mosaic',
      #ftmachine='mosaic',
      #mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=10,
      imsize=1000,
      cell=1,
      niter=500,
      #threshold=str(threshold) + 'Jy',
      robust=0.5,
      outertaper='40arcsec',
      #multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(imagesDir + 'leop.mosaic.vla-c+d.' + name + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value'] # arcsec
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 4.
imagesize = 2 * int(34*60. / cellsize) # images 2X entire primary beam (PB)

if imagesize%100 != 0: imagesize = imagesize + 100 - imagesize%100

# define multiscale 
multiscale = np.array([0.,1.,5.]) * beamsize / cellsize


name='res1'
clean(vis=[uvdataDir + gmrt],
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      #imagermode='mosaic',
      #ftmachine='mosaic',
      #mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=100,
      imsize=imagesize,
      cell=cellsize,
      niter=10000,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      outertaper='40arcsec',
      multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

################
# Resolution 2 #
################

name='test.res2'
clean(vis=[uvdataDir + gmrt],
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      #imagermode='mosaic',
      #ftmachine='mosaic',
      #mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=10,
      imsize=1000,
      cell=1,
      niter=500,
      #threshold=str(threshold) + 'Jy',
      robust=0.5,
      outertaper='10arcsec',
      #multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(imagesDir + 'leop.mosaic.gmrt.' + name + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value'] # arcsec
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 4.
imagesize = int(34*60. / cellsize) # images 2X entire primary beam (PB)

if imagesize%100 != 0: imagesize = imagesize + 100 - imagesize%100

# define multiscale 
multiscale = np.array([0.,1.,5.]) * beamsize / cellsize


name='res2'
clean(vis=[uvdataDir + gmrt],
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      #imagermode='mosaic',
      #ftmachine='mosaic',
      #mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=100,
      imsize=imagesize+4,
      cell=cellsize,
      niter=10000,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      outertaper='10arcsec',
      multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

################
# Resolution 3 #
################

name='test.res3'
clean(vis=[uvdataDir + gmrt],
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      #imagermode='mosaic',
      #ftmachine='mosaic',
      #mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=10,
      imsize=1000,
      cell=1,
      niter=500,
      #threshold=str(threshold) + 'Jy',
      robust=0.5,
      outertaper='25arcsec',
      #multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(imagesDir + 'leop.mosaic.gmrt.' + name + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value'] # arcsec
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 4.
imagesize = 2 * int(34*60. / cellsize) # images 2X entire primary beam (PB)

if imagesize%100 != 0: imagesize = imagesize + 100 - imagesize%100

# define multiscale 
multiscale = np.array([0.,1.,5.]) * beamsize / cellsize


name='res3'
clean(vis=[uvdataDir + gmrt],
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      #imagermode='mosaic',
      #ftmachine='mosaic',
      #mosweight=True,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=100,
      imsize=imagesize,
      cell=cellsize,
      niter=10000,
      threshold=str(threshold) + 'Jy',
      robust=0.5,
      outertaper='25arcsec',
      multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')






###################
# Res 2
###################

# Make dirty cube
name='test.res2'
clean(vis=uvdataDir + vla,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=2,
      imsize=1024,
      cell='1.5arcsec',
      niter=1000,
      threshold='1mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='7arcsec',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(imagesDir + 'leop.mosaic.vla.' + name)
threshold = ia.statistics()['sigma'] * 2.5
beamsize = ia.restoringbeam(channel=0)
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 4.
imagesize = 34*60. / cellsize # images entire PB

# define multiscale 
multiscale = [0,1,5] * beamsize / cellsize

# Clean full cube
name='res2'
clean(vis=uvdataDir + vla,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=100,
      imsize=1024,
      cell='1arcsec',
      niter=1000,
      threshold=threshold,
      robust=0.5,
      outertaper='7arcsec',
      multiscale=multiscale,
      cyclefactor=1.5,
      cyclespeedup=50,
      interpolation='nearest',
      restfreq='1.420405752GHz')


###################
# Res 3
###################

# Make dirty cube
name='test.res3'
clean(vis=uvdataDir + vla,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=4,
      imsize=1024,
      cell='4arcsec',
      niter=1000,
      threshold='1mJy',
      robust=0.5,
      outertaper='14arcsec',
      cyclefactor=1.5,
      cyclespeedup=50,
      interpolation='nearest',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(imagesDir + 'leop.mosaic.vla.' + name)
threshold = ia.statistics()['sigma'] * 2.5
beamsize = ia.restoringbeam(channel=0)
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 4.
imagesize = 34*60. / cellsize # images entire PB

# define multiscale 
multiscale = [0,1,5] * beamsize / cellsize

# Clean full cube
name='res3'
clean(vis=uvdataDir + vla,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      imagermode='mosaic',
      ftmachine='mosaic',
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=100,
      imsize=1024,
      cell='1arcsec',
      niter=1000,
      threshold=threshold,
      robust=0.5,
      outertaper='7arcsec',
      multiscale=multiscale,
      cyclefactor=1.5,
      cyclespeedup=50,
      interpolation='nearest',
      restfreq='1.420405752GHz')


# image GMRT without any uv weighting to find beamsize
name='test'
clean(vis=[uvdataDir + gmrt],
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      scaletype='SAULT',
      pbcor=False,
      minpb=0.3,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=10,
      imsize=1000,
      cell=1,
      niter=500,
      #threshold=str(threshold) + 'Jy',
      robust=0.5,
      #outertaper='40arcsec',
      #multiscale=multiscale.astype(int),
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.7,
      interpolation='nearest',
      restfreq='1.420405752GHz')
# beam = 3.6" x 2.8"


################################################################################
################################################################################
################################################################################
# begin measuring flux by reading in cubes imaged in miriad 
################################################################################
################################################################################
################################################################################

'''
Below is a list of the uvdata used to make cubes:
a) VLA D
b) VLA C/D
c) VLA B/C/D
d) GMRT
e) VLA C/D/GMRT

Each uvdataset is imaged at three different resolutions by using uvtaper
weights of 40, 20, and 1 arcsec. Now we will commence measuring the flux in
each cube.

'''

# import the relevant cubes to casa

imageDir = '/d/bip3/ezbc/leop/data/hi/casa/images/'
mirDir = '/d/bip3/ezbc/leop/data/hi/miriad/fitsImages/'

names = ['vla-d','vla-c+d','vla-b+c+d','vla-c+d.gmrt','gmrt']
fixednames = ['vla_d','vla_cd','vla_bcd','vla_cd.gmrt','gmrt']
filenames = []

for i, telescope in enumerate(names):
    for j, res in enumerate(['res1','res2','res3']):
        for k, image in enumerate(['convol','restor']):
            importfits(fitsimage=mirDir + 'leop.mosaic.' + \
                            telescope + '.' + res + '.' + image + '.fits',
                       imagename=imageDir + 'leop.mosaic.' + \
                            fixednames[i] + '.' + res + '.' + image + '.image')
            filenames.append('leop.mosaic.' + fixednames[i] + '.' + \
                    res + '.' + image)

for i, telescope in enumerate(names):
    for j, res in enumerate(['res1','res2','res3']):
        for k, image in enumerate(['convol','restor']):
            filenames.append('leop.mosaic.' + fixednames[i] + '.' + \
                    res + '.' + image)


dummyMs = uvdataDir + vla
dummyMs = '../../uvdata/leopVLAc.contsub.ms'

def makeImages(imageDir,image):
    # determine beamsize of cube
    ia.open(imageDir + image + '.image')
    beamsize = ia.restoringbeam(channel=0)['major']['value']
    beamsizeUnit = ia.restoringbeam(channel=0)['major']['unit']
    ia.close()

    # create cube for blanking
    immath(imagename=imageDir + image + '.image',
           outfile=imageDir + image + '.blk.image',
            mode='evalexpr',
            region='centerbox[[10h21m45s,18.05.14.9],[10arcmin,10arcmin]]',
            expr='IM0')

    # convolve cube to 2X beam for blanking 
    imsmooth(imagename=imageDir + image + '.image',
         outfile=imageDir + image + '.smooth.image',
         major=str(beamsize*2) + beamsizeUnit,
         minor=str(beamsize*2) + beamsizeUnit,
         pa=0,
         region='centerbox[[10h21m45s,18.05.14.9],[10arcmin,10arcmin]]',
         targetres=True)

    # determine threshold of cube
    ia.open(imageDir + image + '.smooth.image')
    threshold = ia.statistics()['sigma'][0] * 2.5
    ia.close()

    # blank the cube at 2.5sigma
    ia.open(imageDir + image + '.smooth.image')
    ia.calcmask(mask=image + '.smooth.image > ' + str(threshold),
             name='mask1')
    ia.close()

    # hand blank the cube
    im.open(dummyMs)
    im.drawmask(image=imageDir + image + '.smooth.image',
            mask=imageDir + image + '.mask')
    im.close

    # mask contains values of 0 and 1, change to a mask with only values of 1
    ia.open(image + '.mask')
    ia.calcmask(image + '.mask' + '>0.5')
    ia.close()

    # apply mask on smoothed image
    immath(imagename=imageDir + image + '.smooth.image',
       outfile=imageDir + image + '.smooth.blk.image',
       mode='evalexpr',
       mask='mask(' + image + '.mask)',
       expr='IM0')

    # apply mask on image
    ia.open(imageDir + image + '.blk.image')
    ia.maskhandler('copy',[image + '.smooth.blk.image:mask0', 'newmask'])
    ia.maskhandler('set','newmask')
    ia.done()

    # create moment 0 map
    immoments(imagename=imageDir + image + '.blk.image',
              moments=[0],
            axis='spectra',
              chans='',
          outfile=imageDir + image + '.mom0.blk.image')


for i, image in enumerate(filenames):
    makeImages(imageDir,image)


for i, image in enumerate(filenames):
    flux = imstat(imageDir + image + '.mom0.blk.image')['flux'][0]
    ia.open(imageDir + image + '.mom0.blk.image')
    beammaj = ia.restoringbeam(channel=0)['major']['value']
    beammin = ia.restoringbeam(channel=0)['minor']['value']
    beamsizeUnit = ia.restoringbeam(channel=0)['major']['unit']
    ia.close()
    print 'Image: ' + str(image)
    print 'Beamsize: ' + str(beammaj) + '" X ' + str(beammin) + '"'
    print 'Flux: ' + str(flux) + ' Jy km/s'

################################################################################
# VLA/C/D + GMRT data on same grid
# import the relevant cubes to casa

imageDir = '/d/bip3/ezbc/leop/data/hi/casa/images/'
mirDir = '/d/bip3/ezbc/leop/data/hi/miriad/fitsImages/'

names = ['vla-c+d.gmrt.grid']
fixednames = ['vla_cd.gmrt.grid']
filenames = []

for i, telescope in enumerate(names):
    for j, res in enumerate(['res1','res2','res3','res4']):
        for k, image in enumerate(['restor']):
            importfits(fitsimage=mirDir + 'leop.mosaic.' + \
                            telescope + '.' + res + '.' + image + '.fits',
                       imagename=imageDir + 'leop.mosaic.' + \
                            fixednames[i] + '.' + res + '.' + image + '.image')
            filenames.append('leop.mosaic.' + fixednames[i] + '.' + \
                    res + '.' + image)

filenames = []

for i, telescope in enumerate(names):
    for j, res in enumerate(['res1','res2','res3']):
        for k, image in enumerate(['convol','restor']):
            filenames.append('leop.mosaic.' + fixednames[i] + '.' + \
                    res + '.' + image)


for i, image in enumerate(filenames):
    makeImages(imageDir,image)



'''

Beamsizes and fluxes:
---------------------------------

Image: leop.mosaic.vla_d.res1.restor
Beamsize: 66.0609439016" X 60.2343812585"
Flux: 1.00384263385Jy km/s

Image: leop.mosaic.vla_d.res2.restor
Beamsize: 56.2864378093" X 49.7059296815"
Flux: 1.03179854026Jy km/s

Image: leop.mosaic.vla_d.res3.restor
Beamsize: 53.1238082796" X 46.4181751012"
Flux: 1.02284738326Jy km/s

Image: leop.mosaic.vla_cd.res1.restor
Beamsize: 57.6243974268" X 54.5605234801"
Flux: 0.625313997909Jy km/s

Image: leop.mosaic.vla_cd.res2.restor
Beamsize: 35.2064274251" X 33.6717043072"
Flux: 1.1801143653Jy km/s

Image: leop.mosaic.vla_cd.res3.restor
Beamsize: 21.4080378413" X 20.0760511681"
Flux: 1.33692634198Jy km/s

Image: leop.mosaic.vla_bcd.res1.restor
Beamsize: 51.5733305364" X 48.9740990101"
Flux: 1.06966806032Jy km/s

Image: leop.mosaic.vla_bcd.res2.restor
Beamsize: 29.9947176129" X 28.6223452538"
Flux: 1.16030827513Jy km/s

Image: leop.mosaic.vla_bcd.res3.restor
Beamsize: 10.6733738445" X 9.90427806972"
Flux: 1.04380022539Jy km/s

Image: leop.mosaic.vla_cd.gmrt.res1.restor
Beamsize: 57.0488288998" X 53.8925193248"
Flux: 0.854195508659Jy km/s

Image: leop.mosaic.vla_cd.gmrt.res2.restor
Beamsize: 30.6111596525" X 29.8303954303"
Flux: 1.09621401845Jy km/s

Image: leop.mosaic.vla_cd.gmrt.res3.restor
Beamsize: 8.6741818115" X 7.65103893354"
Flux: 1.22011245166Jy km/s

Image: leop.mosaic.gmrt.res1.restor
Beamsize: 47.0790915192" X 43.2303864509"
Flux: 0.837953120579Jy km/s

Image: leop.mosaic.gmrt.res2.restor
Beamsize: 30.1560558379" X 28.4224066883"
Flux: 0.866713690033Jy km/s

Image: leop.mosaic.gmrt.res3.restor
Beamsize: 3.73015017248" X 3.24872068595"
Flux: 1.14209258803Jy km/s
'''

################################################################################
# VLA/C/D + GMRT imaging without mosaicing
################################################################################

# The following if from a helpdesk ticket with Kristina Nyland

# 2) Just how extended is the emission you're trying to image with respect to
# the primary beam? For the VLA at L-band, the primary beam is 30 arcmin. I
# read something in your previous discussion with Miriam where you said you
# expect the real flux to be near the center of the image anyway. Do you even
# NEED to mosaic in the first place? If you don't, then it would be much
# simpler to just use imagermode = 'csclean', make a large image (imsize =
# 16000 is not unheard of - though I would start with something smaller - your
# computer might not have enough memory!), and use the w-projection algorithm
# to remove errors caused by the non-coplanar baselines of the VLA (I have no
# clue if this is a problem for GMRT or not). You could try feeding clean both
# datasets from the start using this method and see what happens. Mosaicing in
# CASA is not nearly as well-tested as standard imaging (and may have a few
# kinks that need to be worked out), so only use it if you *really* need to.




vlab1 = 'leopVLAb1.contsub.ms'
vlab2 = 'leopVLAb2.contsub.ms'
vlac = 'leopVLAc.contsub.ms'
vlad = 'leopVLAd.contsub.ms'
gmrt = 'leopGMRT.contsub.normwts.ms'

uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/'

uvdata = [uvdataDir + gmrt,
          uvdataDir + vlac,
          uvdataDir + vlad]


###################
# Res 1
###################

# Make dirty cube
name='res1'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.vla.gmrt.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      gridmode='widefield',
      outertaper='50arcsec',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################
# Res 2
###################

# Make dirty cube
name='res2'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.vla.gmrt.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='20arcsec',
      gridmode='widefield',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')



###################
# Res 3
###################

# Make dirty cube
name='res3'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.vla.gmrt.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='1arcsec',
      gridmode='widefield',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')


################################################################################
# Image only VLA data
################################################################################

vlab1 = 'leopVLAb1.contsub.ms'
vlab2 = 'leopVLAb2.contsub.ms'
vlac = 'leopVLAc.contsub.ms'
vlad = 'leopVLAd.contsub.ms'
gmrt = 'leopGMRT.contsub.normwts.ms'

uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/vla/'

uvdata = [uvdataDir + vlac,
          uvdataDir + vlad]


###################
# Res 1
###################

# Make dirty cube
name='res1'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      gridmode='widefield',
      outertaper='50arcsec',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################
# Res 2
###################

# Make dirty cube
name='res2'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='20arcsec',
      gridmode='widefield',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################
# Res 3
###################

# Make dirty cube
name='res3'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.vla.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='1arcsec',
      gridmode='widefield',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

################################################################################
# Image only VLA data, simply 
################################################################################

uvdir = '/d/bip3/ezbc/leop/data/hi/casa/reductionFiles/'
uvdata = [uvdir + 'vla.c/leopVLAc.lsr.ms.contsub',
          uvdir + 'vla.d/leopVLAd.lsr.ms.contsub']
imagename = 'leop.vlacd.mosaic'

clean(vis = uvdata,
    imagename = '.dirty',
    mode = 'velocity',
    outframe='lsrk',
    start='185km/s',
    width='1.8km/s',
    nchan=80,
    restfreq = '1420.405752MHz',
    niter = 0,
    imsize = 1024,
    cell = '2arcsec',
    weighting = 'briggs',
    robust = 0.5)

# determine threshold and beamsize of dirty cube
ia.open(imagename + '.dirty.image')
threshold = ia.statistics()['sigma'][0] * 2.5
ia.close()

# deconvolve
deconvolve(
    imagename = imagename + '.dirty.image',
    model = imagename + '.dirty.model',
    psf = imagename + '.dirty.psf',
    alg = 'clark',
    niter = 10000,
    threshold = str(threshold) + 'Jy')




################################################################################
# Image only VLA data with many pixels
################################################################################


###################
# Res 1
###################

# Make dirty cube
name='res1'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.vla.large.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=2500,
      cell='0.8arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      gridmode='widefield',
      outertaper='50arcsec',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################
# Res 2
###################

# Make dirty cube
name='res2'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.vla.large.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=2500,
      cell='0.8arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='20arcsec',
      gridmode='widefield',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################
# Res 3
###################

# Make dirty cube
name='res3'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.vla.large.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=2500,
      cell='0.8arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='1arcsec',
      gridmode='widefield',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')


################################################################################
# Image GMRT with no wts
################################################################################

vlab1 = 'leopVLAb1.contsub.ms'
vlab2 = 'leopVLAb2.contsub.ms'
vlac = 'leopVLAc.contsub.ms'
vlad = 'leopVLAd.contsub.ms'
gmrt = 'leopGMRT.contsub.nowts.ms'

uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/'

uvdata = [uvdataDir + gmrt
         ]


###################
# Res 1
###################

# Make dirty cube
name='nowts.res1'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      gridmode='widefield',
      outertaper='50arcsec',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################
# Res 2
###################

# Make dirty cube
name='nowts.res2'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='20arcsec',
      gridmode='widefield',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')



###################
# Res 3
###################

# Make dirty cube
name='nowts.res3'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.mosaic.gmrt.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='1arcsec',
      gridmode='widefield',
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')


################################################################################
# Image VLA with mosaic = True 
################################################################################

vlab1 = 'leopVLAb1.contsub.ms'
vlab2 = 'leopVLAb2.contsub.ms'
vlac = 'leopVLAc.contsub.ms'
vlad = 'leopVLAd.contsub.ms'
gmrt = 'leopGMRT.contsub.nowts.ms'

uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/'

uvdata = [uvdataDir + vlac,uvdataDir + vlad]


###################
# Res 1
###################

# Make dirty cube
name='res0'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.vlacd.mosaic.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='3mJy',
      robust=0.5,
      uvtaper=True,
      gridmode='widefield',
      outertaper='40arcsec',
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

###################
# Res 2
###################

# Make dirty cube
name='res1'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.vlacd.mosaic.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2.5mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='20arcsec',
      gridmode='widefield',
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')



###################
# Res 3
###################

# Make dirty cube
name='res2'
clean(vis=uvdata,
      imagename=imagesDir + 'leop.vlacd.mosaic.' + name,
      mode='velocity',
      outframe='lsrk',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=10000,
      threshold='2mJy',
      robust=0.5,
      uvtaper=True,
      outertaper='1arcsec',
      gridmode='widefield',
      imagermode='mosaic',
      ftmachine='mosaic',
      mosweight=True,
      cyclefactor=1.5,
      cyclespeedup=50,
      gain=0.5,
      interpolation='nearest',
      restfreq='1.420405752GHz')

################################################################################
################################################################################
################################################################################
# Grid the GMRT, VLA/C, and VLA/D observations onto the same frame.
################################################################################
################################################################################
################################################################################

vlab1 = 'leopVLAb1.contsub.ms'
vlab2 = 'leopVLAb2.contsub.ms'
vlac = 'leopVLAc.contsub.ms'
vlad = 'leopVLAd.contsub.ms'
gmrt = 'leopGMRT.contsub.normwts.ms'

uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/'

uvdata = [vlac,vlad,gmrt]
outputs = ['leopVLAc.contsub.cvel','leopVLAd.contsub.cvel',
           'leopGMRT.contsub.cvel']

for i in range(3):
    # change the source name of each dataset to be the same
    tb.open(uvdata[i] + '/FIELD',nomodify=False)
    st=tb.selectrows(0)
    st.putcol('NAME','LeoP')
    st.done()
    tb.close()
    # put the three datasets on the same grid
    cvel(vis=uvdataDir + uvdata[i],
         outputvis=uvdataDir + outputs[i] + '.ms',
         outframe='lsrk')
    # export the fits file
    exportuvfits(vis=uvdataDir + outputs[i] + '.ms',
                 fitsfile=uvdataDir + outputs[i] + '.fits',
                 multisource=False)


################################################################################################################################################################################################################################################

# import cube john made

importfits(fitsimage='/d/bip3/ezbc/leop/data/hi/casa/fitsImages/4.cube.fits',
           imagename='leop.gmrt.john.4arcsec.image')

# imaged cube in miriad with original GMRT observations

importfits(fitsimage='/d/bip3/ezbc/leop/data/hi/miriad/fitsImages/'+ \
                     'leop.gmrt.original.restor.fits',
           imagename='leop.gmrt.original.4arcsec.image')

imagedir = '/d/bip3/ezbc/leop/data/hi/casa/images/vla.gmrt/'

images = ['16','16.SR','32','32.SR','4','4.SR','8','8.SR']

for i, image in enumerate(images):
    importfits(fitsimage='/d/bip3/ezbc/leop/data/hi/casa/fitsImages/' + \
                         image + '.fits',
            imagename=imagedir + 'leop.aips.vla.gmrt.' + image + '.image')

from hianalysis import blankcube

dummyMS = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/leopVLAc.contsub.ms'

for i, image in enumerate(images):
    blankcube('./','leop.aips.vla.gmrt.' + image,dummyMS=dummyMS,pbcor=False,
            ruthless=True,smooth=True)


for i, image in enumerate(images):
    image = 'leop.aips.vla.gmrt.' + image
    flux = imstat(imagedir + image + '.mom0.blk.image')['flux'][0]
    ia.open(imagedir + image + '.mom0.blk.image')
    beammaj = ia.restoringbeam(channel=0)['major']['value']
    beammin = ia.restoringbeam(channel=0)['minor']['value']
    beamsizeUnit = ia.restoringbeam(channel=0)['major']['unit']
    ia.close()
    print 'Image: ' + str(image)
    print 'Beamsize: ' + str(beammaj) + '" X ' + str(beammin) + '"'
    print 'Flux: ' + str(flux) + ' Jy km/s'


'''
Regular cubes:
Image: leop.aips.vla.gmrt.32
Beamsize: 32.0000009251" X 31.9782040107"
Flux: 1.07569418939 Jy km/s

Image: leop.aips.vla.gmrt.16
Beamsize: 15.9999242929" X 15.8491495899"
Flux: 1.20993601001 Jy km/s

Image: leop.aips.vla.gmrt.8
Beamsize: 7.9999999931" X 7.99999895014"
Flux: 1.08020122198 Jy km/s

Image: leop.aips.vla.gmrt.4
Beamsize: 3.95968530752" X 3.7654933362"
Flux: 0.951149597574 Jy km/s

Flux rescaled cubes:
Image: leop.aips.vla.gmrt.32.SR
Beamsize: 32.0000009251" X 31.9782040107"
Flux: 1.03596761977 Jy km/s

Image: leop.aips.vla.gmrt.16.SR
Beamsize: 15.9999242929" X 15.8491495899"
Flux: 0.934338727197 Jy km/s

Image: leop.aips.vla.gmrt.8.SR
Beamsize: 7.9999999931" X 7.99999895014"
Flux: 0.646697581212 Jy km/s

Image: leop.aips.vla.gmrt.4.SR
Beamsize: 3.95968530752" X 3.7654933362"
Flux: 0.308148711513 Jy km/s

'''



################################################################################################################################################################################################################################################
# Perform residual 
################################################################################################################################################################################################################################################


imagedir = '/d/bip3/ezbc/leop/data/hi/casa/images/vla.gmrt/'

images = ['4','8','16','32']

imagenames = []

for i, image in enumerate(images):
    importfits(fitsimage='/d/bip3/ezbc/leop/data/hi/miriad/revisedImages/' + \
                         'leop.mosaic.lsr.' + image + 'arcsec.residual.fits',
            imagename=imagedir + 'leop.miriad.vla.gmrt.' + image + \
                      'arcsec.residual.image')
    importfits(fitsimage='/d/bip3/ezbc/leop/data/hi/miriad/revisedImages/' + \
                         'leop.mosaic.lsr.' + image + 'arcsec.clean.fits',
            imagename=imagedir + 'leop.miriad.vla.gmrt.' + image + \
                      'arcsec.clean.image')
    importfits(fitsimage='/d/bip3/ezbc/leop/data/hi/miriad/revisedImages/' + \
                         'leop.mosaic.lsr.' + image + 'arcsec.dirty.fits',
            imagename=imagedir + 'leop.miriad.vla.gmrt.' + image + \
                      'arcsec.dirty.image')
    importfits(fitsimage='/d/bip3/ezbc/leop/data/hi/miriad/revisedImages/' + \
                         'leop.mosaic.lsr.' + image + 'arcsec.psf.fits',
            imagename=imagedir + 'leop.miriad.vla.gmrt.' + image + \
                      'arcsec.psf.image')
    imagenames.append('leop.miriad.vla.gmrt.' + image + 'arcsec')

imagenames = ['leop.miriad.vla.gmrt.4arcsec',
             'leop.miriad.vla.gmrt.8arcsec',
             'leop.miriad.vla.gmrt.16arcsec',
             'leop.miriad.vla.gmrt.32arcsec']
# to perform flux rescale we need to perform the following: G = D x C / (D - R)
# where G is true flux, D is dirty image, C is clean image, and R is residual
# image.

def rescale_flux(imagename):
    from casa import immath
    from os import path
    if not path.isdir(imagename+'.rescaled.image.temp'):
        immath(imagename=[imagename+'.dirty.image',
                          imagename+'.clean.image',
                          imagename+'.residual.image'],
               outfile=imagename+'.rescaled.image.temp',
               expr='IM0*IM1 / (IM0 - IM2)')
    from casa import deconvolve
    from os import system
    deconvolve(imagename=imagename+'.rescaled.image.temp',
               model=imagename+'.rescaled',
               niter=0)
    system('rm -rf ' + imagename + '.rescaled.image.temp')
    system('mv ' + imagename + '.im.restored ' + imagename + '.rescaled.image')

for i, imagename in enumerate(imagenames):
    rescale_flux(imagename)

import hianalysis as im

dummyMS = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/leopGMRT.contsub.ms'

for i, imagename in enumerate(imagenames):
    im.blankcube('./',imagename + '.rescaled',
                 dummyMS=dummyMS,pbcor=False,
                 ruthless=True,smooth=True)






################################################################################
# VLA/C/D + GMRT imaging without mosaicing, with multiscale clean
################################################################################

vlac = 'leopVLAc.contsub.cvel.ms'
vlad = 'leopVLAd.contsub.cvel.ms'
gmrt = 'leopGMRT.contsub.normwts.ms'

uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/casaImages/'

uvdata = [uvdataDir + gmrt,
          uvdataDir + vlac,
          uvdataDir + vlad]

###################
# Res 1
###################

# Make dirty cube
name='res1'
image = 'leop.vla.gmrt.' + name

clean(vis=uvdata,
      imagename=image,
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1024,
      cell='2arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper='50arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,1,5]) * beamsize / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale)

# to perform flux rescale we need to perform the following: G = D x C / (D - R)
# where G is true flux, D is dirty image, C is clean image, and R is residual
# image.

import hianalysis as im

dummyMS = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/leopVLAc.contsub.cvel.ms'

im.blankcube('./',imagename + '.model.restored',
                 dummyMS=dummyMS,pbcor=False,
                 ruthless=True,smooth=True)

################################################################################
# performed uvtaper tests, imaged with various uvtapers using the following
# clean paramters:
#  clean(vis=uvdata,
#     imagename='uvtapertest.'+str(uvtaper),
#     mode='velocity',
#     start='270km/s',
#     width='1.8km/s',
#     nchan=1,
#     imsize=1600,
#     cell='1.1arcsec',
#     niter=0,
#     threshold='2mJy',
#     robust=0.5,
#     weighting='briggs',
#     uvtaper=True,
#     outertaper=str(uvtaper)+'arcsec',
#     restfreq='1.420405752GHz')

def makeimage(uvtaper='',iteration=0):
    vlac = 'leopVLAc.contsub.cvel.ms'
    vlad = 'leopVLAd.contsub.cvel.ms'
    gmrt = 'leopGMRT.contsub.normwts.ms'

    uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
    imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/casaImages/'

    uvdata = [uvdataDir + gmrt,
              uvdataDir + vlac,
              uvdataDir + vlad]

    from casa import clean
    clean(vis=uvdata,
      imagename='uvtapertest.'+str(uvtaper),
      mode='velocity',
      start='270km/s',
      width='1.8km/s',
      nchan=1,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      threshold='2mJy',
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

uvtapers = [0.3,0.4,0.5,0.6,0.8,1,1.5,2,3,4,6,8,10,12,15,18,20,22]
uvtapers = [0.2,24,26]
uvtapers = [0.02,0.05]


from os import system,chdir

system('mkdir uvtapertests')
chdir('uvtapertests')

for i, uvtaper in enumerate(uvtapers):
    makeimage(uvtaper=uvtaper,iteration=i)

for i, uvtaper in enumerate(uvtapers):
    ia.open('uvtapertest.'+str(uvtaper))
    beammaj = ia.restoringbeam(channel=0)['major']['value']
    beammin = ia.restoringbeam(channel=0)['minor']['value']
    beamsizeUnit = ia.restoringbeam(channel=0)['major']['unit']
    ia.close()
    print 'Beamsize: ' + str(beammaj) + '" X ' + str(beammin) + '"'



# uvtaper | beamsize (")
# 
# 1.5       7.58 x 7.32238
# 8         15.93 x 15.36
# 24        31.69 x 31.009 



################################################################################
# 4"

# Make dirty cube
name='4arcsec'
image = 'leop.vla.gmrt.' + name
uvtaper = 

clean(vis=uvdata,
      imagename=image
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,1,5]) * beamsize / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale)

################################################################################
# 8"

# Make dirty cube
name='8arcsec'
image = 'leop.vla.gmrt.' + name
uvtaper = 1.5

clean(vis=uvdata,
      imagename=image,
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,1,5]) * beamsize / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale)



################################################################################
# 16"

# Make dirty cube
name='16arcsec'
image = 'leop.vla.gmrt.' + name
uvtaper = 8

clean(vis=uvdata,
      imagename=image
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,1,5]) * beamsize / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale)



################################################################################
# 32"

# Make dirty cube
name='32arcsec'
image = 'leop.vla.gmrt.' + name
uvtaper = 24

clean(vis=uvdata,
      imagename=image
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,1,5]) * beamsize / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale)

################################################################################
################################################################################
################################################################################
# Now blank the cubes

import hianalysis as im

dummyMS = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/leopVLAc.contsub.cvel.ms'

reslist = ['4','8','16','32']
reslist = ['16']
for i, res in enumerate(reslist):
    imagename = 'leop.vla.gmrt.' + res + 'arcsec'
    im.blankcube('./',imagename,
                 dummyMS=dummyMS,extension='.model.restored',pbcor=False,
                 ruthless=True,smooth=True)

filenames=[]
for i, res in enumerate(reslist):
    filenames.append('leop.vla.gmrt.' + res + 'arcsec')

im.printProperties(filenames)


################################################################################
################################################################################
################################################################################
# 16" has too much flux, try using reference requency
################################################################################################################################################################################################################################################

vlac = 'leopVLAc.contsub.ms'
vlad = 'leopVLAd.contsub.ms'
gmrt = 'leopGMRT.contsub.normwts.ms'

uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/casaImages/'

uvdata = [uvdataDir + gmrt,
          uvdataDir + vlac,
          uvdataDir + vlad]



# 8"

# Make dirty cube
name='8arcsec'
image = 'leop.vla.gmrt.lsrk.' + name
uvtaper = 1.5

clean(vis=uvdata,
      imagename=image,
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      outframe='lsrk',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,1,5]) * beamsize / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale)



################################################################################
# 16"

# Make dirty cube
name='16arcsec'
image = 'leop.vla.gmrt.lsrk.' + name
uvtaper = 8

clean(vis=uvdata,
      imagename=image,
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      outframe='lsrk',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,1,5]) * beamsize / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale)



################################################################################
# 32"

# Make dirty cube
name='32arcsec'
image = 'leop.vla.gmrt.lsrk.' + name
uvtaper = 24

clean(vis=uvdata,
      imagename=image,
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      outframe='lsrk',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,1,5]) * beamsize / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale)

################################################################################
################################################################################
################################################################################
# Now blank the cubes

import hianalysis as im

dummyMS = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/leopVLAc.contsub.cvel.ms'

reslist = ['4','8','16','32']
reslist = ['8','16','32']
for i, res in enumerate(reslist):
    imagename = 'leop.vla.gmrt.lsrk.' + res + 'arcsec'
    im.blankcube('./',imagename,
                 dummyMS=dummyMS,extension='.model.restored',pbcor=False,
                 ruthless=True,smooth=True)

filenames=[]
for i, res in enumerate(reslist):
    filenames.append('leop.vla.gmrt.' + res + 'arcsec')

im.printProperties(filenames)




################################################################################
################################################################################
################################################################################
# 16" has too much flux, try implementing multiscale more similar to Rich et
# al. 2008
################################################################################################################################################################################################################################################

vlac = 'leopVLAc.contsub.ms'
vlad = 'leopVLAd.contsub.ms'
gmrt = 'leopGMRT.contsub.normwts.ms'

uvdataDir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'
imagesDir = '/d/bip3/ezbc/leop/data/hi/casa/images/casaImages/'

uvdata = [uvdataDir + gmrt,
          uvdataDir + vlac,
          uvdataDir + vlad]
# 8"

# Make dirty cube
name='8arcsec'
image = 'leop.vla.gmrt.lsrk.scale.' + name
uvtaper = 1.5

clean(vis=uvdata,
      imagename=image,
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      outframe='lsrk',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,10,30]) / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale,
           gain=0.7)



################################################################################
# 16"

# Make dirty cube
name='16arcsec'
image = 'leop.vla.gmrt.lsrk.scale.' + name
uvtaper = 8

clean(vis=uvdata,
      imagename=image,
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      outframe='lsrk',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,13.3,40]) / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale,
           gain=0.7)



################################################################################
# 32"

# Make dirty cube
name='32arcsec'
image = 'leop.vla.gmrt.lsrk.scale.' + name
uvtaper = 24

clean(vis=uvdata,
      imagename=image,
      mode='velocity',
      start='190km/s',
      width='1.8km/s',
      outframe='lsrk',
      nchan=80,
      imsize=1600,
      cell='1.1arcsec',
      niter=0,
      robust=0.5,
      weighting='briggs',
      uvtaper=True,
      outertaper=str(uvtaper)+'arcsec',
      restfreq='1.420405752GHz')

# determine threshold and beamsize of dirty cube
ia.open(image + '.image')
threshold = ia.statistics()['sigma'][0] * 2.5
beamsize = ia.restoringbeam(channel=0)['major']['value']
ia.close()

# define imagesize and cellsize for 
cellsize = beamsize / 5.
cellsize = 2
imagesize = 34*60. / cellsize # images entire PB

import numpy as np

# define multiscale 
multiscale = np.rint(np.array([0,20,60]) / cellsize).astype(int)

deconvolve(imagename=image + '.image',
           model=image + '.model',
           psf=image+'.psf',
           niter=10000,
           threshold=str(threshold) + 'Jy',
           alg='multiscale',
           scales=multiscale,
           gain=0.7)



################################################################################
# Now blank the cubes

import hianalysis as im

dummyMS = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/leopVLAc.contsub.cvel.ms'

reslist = ['8','16','32']
for i, res in enumerate(reslist):
    imagename = 'leop.vla.gmrt.lsrk.scale.' + res + 'arcsec'
    im.blankcube('./',imagename,
                 dummyMS=dummyMS,extension='.model.restored',pbcor=False,
                 ruthless=True,smooth=True)

filenames=[]
for i, res in enumerate(reslist):
    filenames.append('leop.vla.gmrt.lsrk.scale.' + res + 'arcsec')

im.printProperties(filenames)

# Huzzah!
#   Image: Multiscale clean resolution #0
#   Beamsize: 9.00000052972" X 9.00000043764"
#   Std: 0.000521610171802 Jy/Bm
#   Flux: 0.403214099764 Jy km/s
#
#   Image: Multiscale clean resolution #1
#   Beamsize: 15.999999265" X 15.9997815078"
#   Std: 0.000722466203225 Jy/Bm
#   Flux: 0.985104843301 Jy km/s
#
#   Image: Multiscale clean resolution #2
#   Beamsize: 31.0000028328" X 31.0000015405"
#   Std: 0.000855918069491 Jy/Bm
#   Flux: 1.03334659165 Jy km/s

# write out fits files:

images=[]
for i, res in enumerate(reslist):
    images.append('leop.vla.gmrt.lsrk.scale.' + res + 'arcsec')

for i, image in enumerate(images):
    exportfits(fitsimage='/d/bip3/ezbc/leop/data/hi/casa/fitsImages/' + \
                         image + '.cube.fits',
            imagename=image + '.image')
    exportfits(fitsimage='/d/bip3/ezbc/leop/data/hi/casa/fitsImages/' + \
                         image + '.mom0.fits',
            imagename=image + '.mom0.blk.image')








