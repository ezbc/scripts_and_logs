#!/usr/bin/python
################################################################################
################################################################################
# Script to analyze GMRT + VLA/C/D observations in CASA
# Cubes imaged in AIPS by John M. Cannon, jcannon@macalester.edu
################################################################################
################################################################################




################################################################################
################################################################################
################################################################################

from os import chdir
import hianalysis as im
import pyfits as pf

# import cubes John made
# the relevant cubes are:
# with residual rescaling:
#   leop.vla.gmrt.fluxRescale.16arcsec.fits
#   leop.vla.gmrt.fluxRescale.4arcsec.fits
#   leop.vla.gmrt.fluxRescale.32arcsec.fits
#   leop.vla.gmrt.fluxRescale.8arcsec.fits

# without residual rescaling:
#   leop.vla.gmrt.16arcsec.fits
#   leop.vla.gmrt.4arcsec.fits
#   leop.vla.gmrt.32arcsec.fits
#   leop.vla.gmrt.8arcsec.fits

# set the working computer
lappy = False
cosmos = True
# located on the laptop at:
if lappy:
    fitsDirOrig = '/home/elijah/research/leop/data/vla.gmrt/fitsImages/'
elif cosmos:
    fitsDirOrig = '/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/' + \
              'originalCubes/'
if lappy:
    fitsDir = '/home/elijah/research/leop/data/vla.gmrt/fitsImages/'
elif cosmos:
    fitsDir = '/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/'


if lappy:
    chdir('/home/elijah/research/leop/data/vla.gmrt/casaImages')
elif cosmos:
    chdir('/d/bip3/ezbc/leop/data/hi/casa/images/vla.gmrt/')

# without rescaling
# set the name variables
resolutions = [32,16,8,4]
extension = 'arcsec'
basenames = ['leop.','leop.fluxRescale.']
region='centerbox[[10h21m45s,18.05.14.9],[15arcmin,15arcmin]]'

# add reference frame information to all headers
# CASA task imhead changed the spectral values when the header info changed
# f-ing CASA
for j, basename in enumerate(basenames):
    images = []
    for i, resolution in enumerate(resolutions):
        images.append(basename + \
                      str(resolution) + extension + '.fits')
        f = pf.open(fitsDirOrig + 'noRefFrameCubes/' + images[i])
        h, d = f[0].header, f[0].data[:,:,:]
        h.update('specsys','LSRK')
        hdu = pyfits.PrimaryHDU(d,header=h)
        hdulist = pyfits.HDUList([hdu])
        hdulist.writeto(fitsDirOrig + images[i])


# import the fits images into CASA
for j, basename in enumerate(basenames):
    images = []
    for i, resolution in enumerate(resolutions):
        # import fits image
        importfits(fitsimage=fitsDirOrig + basename + \
                   str(resolution) + extension + '.fits',
                imagename=basename + str(resolution) + extension + \
                          '.large.image.tmp')
        # concatenate image names
        images.append(basename + str(resolution) + extension)
        # define rest frequency
        imhead(imagename=images[i] + '.large.image.tmp',
               mode='put',
               hdkey='restfreq',
               hdvalue='1.42040575177E+09')

        # define spectral axis as velocity
        ia.open(images[i] + '.large.image.tmp')
        ia.regrid(outfile=images[i] + '.large.image',
                  asvelocity=True)
        ia.close()

        # switch reference frame to galactic centric
        #imreframe(imagename=images[i] + '.large.image.tmp2',
        #          output=images[i] + '.large.image',
        #          outframe='GALACTO', # galactic centric
        #          restfreq='1.42040575177 GHz')

        # remove temporary image
        os.system('rm -rf ' + images[i] + '.large.image.tmp ')

        # create smaller, more manageable image
        immath(imagename=images[i] + '.large.image',
               outfile=images[i] + '.image',
               mode='evalexpr',
               region=region,
               expr='IM0')

        # write with velocity axis
        exportfits(fitsimage=fitsDir + images[i] + '.vel.fits',
                   imagename=images[i] + '.image',
                   velocity=True,
                   dropstokes=True)
        # Write without velocity axis
        exportfits(fitsimage=fitsDir + images[i] + '.fits',
                   imagename=images[i] + '.image',
                   velocity=False,
                   dropstokes=True)

# rename images
images = []
for j, basename in enumerate(basenames):
    for i, resolution in enumerate(resolutions):
        images.append(basename + str(resolution) + extension)

# create a small image for models
immath(imagename=images[5] + '.image',
               outfile=images[5] + '.crop.vel.image',
               mode='evalexpr',
               region='centerbox[[10h21m45s,18.05.14.9],[10arcmin,10arcmin]]',
               expr='IM0')
ia.open(images[5] + '.crop.vel.image')
ia.regrid(outfile=images[5] + '.crop.freq.image',
          asvelocity=False)
ia.close()
exportfits(fitsimage=fitsDir + images[5] + '.crop.freq.fits',
           imagename=images[5] + '.crop.freq.image',
           velocity=False,
           dropstokes=True)
exportfits(fitsimage=fitsDir + images[5] + '.crop.vel.fits',
           imagename=images[5] + '.crop.vel.image',
           velocity=True,
           dropstokes=True)

# make moment 0,1,2 images
if cosmos:
    dummyMS = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/leopVLAc.contsub.ms'
elif lappy:
    dummyMS = '/home/elijah/research/leop/data/vla/c_config/casa/' + \
              'reductionFiles/dataRR.ms'

for i, image in enumerate(images[4:]):
    im.blankcube(image,
              extension='.image',
              dummyMS=dummyMS,
              ruthless=True,
              smooth=True,
              blankThreshold=2.5,
              moments=[0])

ia.open(image + '.large.image')
threshold = ia.statistics()['sigma'][0] * 5
ia.close()

# create moment images from this threshold
for i in range(0,3):
    for j, image in enumerate(images[4:]):
        immoments(imagename=image + '.blk.image',
                  moments=[i],
                  mask='mask(' + image + '.blk.image)',
                  includepix=[threshold,1e10],
                  outfile=image + '.mom' + str(i) + '.2pt5sig.image')

# create moments blanked at higher threshold
for i, image in enumerate(images[4:6]):
    hianalysis.blankcube(image,
              extension='.image',
              dummyMS=dummyMS,
              ruthless=True,
              smooth=False,
              blankThreshold=5,
              moments=[1,2])

# print properties
for i, image in enumerate(images[4:]):
    flux = imstat(image + '.mom0.image')['flux'][0]
    ia.open(image + '.mom0.image')
    beammaj = ia.restoringbeam(channel=0)['major']['value']
    beammin = ia.restoringbeam(channel=0)['minor']['value']
    beamsizeUnit = ia.restoringbeam(channel=0)['major']['unit']
    ia.close()
    print 'Image: ' + str(image)
    print 'Beamsize: ' + str(beammaj) + '" X ' + str(beammin) + '"'
    print 'Flux: ' + str(flux) + ' Jy km/s'

# convert moment zero images to column density units.
#	Recall:  1 K = (7.354E-8)*[Bmaj(")*Bmin(")/lamda^2(m)] Jy/Bm

#	Here, units of images are Jy/Bm m/s; cellsize = 2"; 	
#	    lambda = 0.211061140507 m

#	Thus, for the 21 cm line of Hydrogen, we have:

#	    1 K = Bmaj(")*Bmin(")/(6.057493205E5) Jy/Bm
#			---- OR ----
#	    1 Jy/Bm = (6.057493205E5)/[Bmaj(")*Bmin(")]

#	Now, recall that: N_HI = (1.8224E18 cm^-2)*[T_b (K)]*int(dv)
#		-- For moment maps in K km/sec, just input the values
#		& multiply by coefficient.
#	   -- Assure that units are Jy/Bm km/sec (i.e., divide by 1000)
#	   Leave in units of 1E20 cm^-2 by dividing by 1E20:

#	   For a x beam: 
#               N_HI (cm^-2) = (image) *
#		[(6.057493205E5)/(*)] * (1/1000) * (1.8224E18 cm^-2) *
#		(1/1E20)
#		N_HI (cm^-2) = (image)*



# Calculate column density
for i, image in enumerate(images[4:]):
    n_hi = 6.057493205e5*1.8224e18/1e20/resolutions[i]**2
    immath(imagename=image + '.mom0.image',
           outfile=image + '.colDens.image',
           mode='evalexpr',
           expr='IM0*'+str(n_hi))
    imhead(imagename=image + '.colDens.image',
           mode='put',
           hdkey='bunit',
           hdvalue='10^20 cm^2')


# blank 16" mom0 map at 10^20 cm^2
# first calculate moment 0 map and noise
image = images[5]

immath(imagename=image + '.colDens.image',
   outfile=image + '.colDens.blk.image',
   mode='evalexpr',
   expr='IM0')

# blank the cube at threshold*sigma 
ia.open(image + '.colDens.blk.image')
ia.calcmask(mask=image + '.colDens.blk.image > ' + str(1),
         name='mask1')
wait = 'waits for calcmask to close'
ia.close()

# mask contains values of 0 and 1, change to a mask with only values of 1
ia.open(image + '.colDens.blk.image')
ia.calcmask(image + '.colDens.blk.image' + '>0.5')
ia.close()

for i in range(3):
    immath(imagename=image + '.mom' + str(i) + '.image',
       outfile=image + '.mom' + str(i) + '.blk.image',
       mode='evalexpr',
       mask='mask(' + image + '.colDens.blk.image)',
       expr='IM0')
    immath(imagename=image + '.mom' + str(i) + '.2pt5sig.image',
       outfile=image + '.mom' + str(i) + '.2pt5sig.blk.image',
       mode='evalexpr',
       mask='mask(' + image + '.colDens.blk.image)',
       expr='IM0')


# Convert to surface density in Msun/pc^2
# Calculate column density
for i, image in enumerate(images[4:]):
    mass = hianalysis.make_SBimage(image=image,
                            extension='.mom0.image',
                            cellsize=1.5,
                            beamsize=resolutions[i],
                            distance=1.76)



# export the fits images
for i, image in enumerate(images[4:]):
    for j in range(3):
        exportfits(fitsimage=fitsDir + image + '.mom' + str(j) + '.fits',
                   imagename=image + '.mom' + str(j) + '.image',
                   dropstokes=True)
        exportfits(fitsimage=fitsDir + image + '.mom' + str(j) + \
                             '.2pt5sig.fits',
                   imagename=image + '.mom' + str(j) + '.2pt5sig.image',
                   dropstokes=True)
        if j == 1:
            exportfits(fitsimage=fitsDir + image + '.mom' + str(j) + \
                                 '.blk.fits',
                       imagename=image + '.mom' + str(j) + '.blk.image',
                       dropstokes=True)
            exportfits(fitsimage=fitsDir + image + '.mom' + str(j) + \
                                 '.2pt5sig.blk.fits',
                       imagename=image + '.mom' + str(j) + '.2pt5sig.blk.image',
                       dropstokes=True)
    exportfits(fitsimage=fitsDir + image + '.colDens.fits',
               imagename=image + '.colDens.image',
               dropstokes=True)
    exportfits(fitsimage=fitsDir + image + '.sb.fits',
               imagename=image + '.sb.image',
               dropstokes=True)



# smooth 16" cube to 24" resolution for Steve's component analysis
if cosmos:
    dummyMS = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/leopVLAc.contsub.ms'
elif lappy:
    dummyMS = '/home/elijah/research/leop/data/vla/c_config/casa/' + \
              'reductionFiles/dataRR.ms'

imsmooth(imagename='leop.fluxRescale.16arcsec.image',
         outfile='leop.fluxRescale.24arcsec.image',
         targetres=True,
         major='24arcsec',
         minor='24arcsec',
         pa='0deg')

imsmooth(imagename='leop.16arcsec.image',
         outfile='leop.24arcsec.image',
         targetres=True,
         major='24arcsec',
         minor='24arcsec',
         pa='0deg')

imsmooth(imagename='leop.16arcsec.large.image',
         outfile='leop.24arcsec.large.image',
         targetres=True,
         major='24arcsec',
         minor='24arcsec',
         pa='0deg')


im.blankcube('leop.fluxRescale.24arcsec',
              extension='.image',
              dummyMS=dummyMS,
              ruthless=True,
              smooth=True,
              beamround=1.,
              blankThreshold=2.5,
              moments=[0])

image = 'leop.fluxRescale.24arcsec'

n_hi = 6.057493205e5*1.8224e18/1e20/24**2
immath(imagename=image + '.mom0.image',
       outfile=image + '.colDens.image',
       mode='evalexpr',
       expr='IM0*'+str(n_hi))
imhead(imagename=image + '.colDens.image',
       mode='put',
       hdkey='bunit',
       hdvalue='10^20 cm^2')


exportfits(fitsimage=fitsDir + 'leop.fluxRescale.24arcsec.fits',
           imagename='leop.fluxRescale.24arcsec.image',
           dropstokes=True,
           velocity=True)

exportfits(fitsimage=fitsDir + 'leop.24arcsec.fits',
           imagename='leop.24arcsec.image',
           dropstokes=True,
           velocity=True)

exportfits(fitsimage=fitsDir + 'leop.fluxRescale.24arcsec.mom0.fits',
           imagename='leop.fluxRescale.24arcsec.mom0.image',
           dropstokes=True,
           velocity=True,
           overwrite=True)

exportfits(fitsimage=fitsDir + 'leop.fluxRescale.24arcsec.colDens.fits',
           imagename='leop.fluxRescale.24arcsec.colDens.image',
           dropstokes=True,
           velocity=True,
           overwrite=True)











