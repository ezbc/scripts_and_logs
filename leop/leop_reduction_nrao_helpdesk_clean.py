#!/usr/bin/python

################################################################################
# Script for imaging HI observations of VLA/D + VLA/C configuration
# observations of Leo P
################################################################################

################################################################################
# Image VLA with mosaic = True 
################################################################################

vlac = 'leopVLAc.contsub.ms'
vlad = 'leopVLAd.contsub.ms'

uvdata = [vlac,vlad]

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


