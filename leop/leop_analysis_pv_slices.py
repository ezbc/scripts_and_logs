#!/usr/bin/python

################################################################################
################################################################################
# Log file for taking PV slices of 16" and 32" cubes along their
# semi-minor axes
################################################################################
################################################################################

import math
import numpy as np

# 16"

beamsize = 16./3600.
startRA = 155.438
startDec = 18.0883
PA = 57.

ras = np.zeros(21)
decs = np.zeros(21)

for i in xrange(0,10):
    ras[i+10] = startRA + beamsize * i * math.sin(math.radians(PA))
    decs[i+10] = startDec + beamsize * i * math.cos(math.radians(PA))

for i in xrange(0,11):
    ras[10 - i] = startRA - beamsize * (i) * math.sin(math.radians(PA))
    decs[10 - i] = startDec - beamsize * (i) * math.cos(math.radians(PA))

for i in xrange(0,20):     
    print(decs[i])

for i in xrange(0,20):     
    print(ras[i])

PA = 57. + 90.

for i in xrange(0,10):
    ras[i+10] = startRA + beamsize * i * math.sin(math.radians(PA))
    decs[i+10] = startDec + beamsize * i * math.cos(math.radians(PA))

for i in xrange(0,11):
    ras[10 - i] = startRA - beamsize * (i) * math.sin(math.radians(PA))
    decs[10 - i] = startDec - beamsize * (i) * math.cos(math.radians(PA))

for i in xrange(0,20):     
    print(decs[i])

for i in xrange(0,20):     
    print(ras[i])


#### 32"

beamsize = 32./3600.
startRA = 155.438
startDec = 18.0883
PA = 57.

for i in xrange(0,10):
    ras[i+10] = startRA + beamsize * i * math.sin(math.radians(PA))
    decs[i+10] = startDec + beamsize * i * math.cos(math.radians(PA))

for i in xrange(0,11):
    ras[10 - i] = startRA - beamsize * (i) * math.sin(math.radians(PA))
    decs[10 - i] = startDec - beamsize * (i) * math.cos(math.radians(PA))

for i in xrange(0,20):     
    print(decs[i])

for i in xrange(0,20):     
    print(ras[i])

ras = np.zeros(21)
decs = np.zeros(21)

PA = 57. + 90.


for i in xrange(0,10):
    ras[i+10] = startRA + beamsize * i * math.sin(math.radians(PA))
    decs[i+10] = startDec + beamsize * i * math.cos(math.radians(PA))

for i in xrange(0,11):
    ras[10 - i] = startRA - beamsize * (i) * math.sin(math.radians(PA))
    decs[10 - i] = startDec - beamsize * (i) * math.cos(math.radians(PA))

for i in xrange(0,20):     
    print(decs[i])

for i in xrange(0,20):     
    print(ras[i])



# Script stops here




# Using the bounding box: 184.74 X 272.9 arcsec for major axis slices
#                         154.95 X 110   arcsec for minor axis slices

##################
# 16" 
##################

# Using the bounding box: 184.74 X 272.9 arcsec for major axis slices
#                         154.95 X 110   arcsec for minor axis slices

# Using averaging width of 5 pixels --> 4"/pix

# major axis cuts:
# filename base: leop.mosaic.16arcsec.pvslice.major.
RA (deg)      Dec (deg)     File #
155.423090301 18.0786175283 01
155.426817726 18.0810381462 02
155.430545151 18.0834587641 03
155.434272575 18.0858793821 04
155.438       18.0883       05
155.441727425 18.0907206179 06
155.445454849 18.0931412359 07
155.449182274 18.0955618538 08

# minor axis cuts:
# filename base: leop.mosaic.16arcsec.pvslice.minor.
RA (deg)      Dec (deg)     File #
155.42589691  18.1069371237 01
155.428317528 18.103209699  02
155.430738146 18.0994822742 03
155.433158764 18.0957548495 04
155.435579382 18.0920274247 05
155.438       18.0883       06
155.440420618 18.0845725753 07
155.442841236 18.0808451505 08
155.445261854 18.0771177258 09
155.447682472 18.073390301  10
155.45010309  18.0696628763 11
155.452523708 18.0659354515 12

fileNums = np.array(['01','02','03','04','05',
                     '06','07','08','09','10',
                     '11','12','13','14','15',
                     '16','17','18','19','20'])

for i in xrange(0,7):
    exportfits(imagename='leop.mosaic.16arcsec.pvslice.major.' + fileNums[i],
               fitsimage='leop.mosaic.16arcsec.pvslice.major.' + \
                         fileNums[i] + '.fits',
               velocity=True)

for i in xrange(0,11):
    exportfits(imagename='leop.mosaic.16arcsec.pvslice.minor.' + fileNums[i],
               fitsimage='leop.mosaic.16arcsec.pvslice.minor.' + \
                         fileNums[i] + '.fits',
               velocity=True)

for i in xrange(0,7):
    exportfits(imagename='leop.mosaic.32arcsec.pvslice.major.' + fileNums[i],
               fitsimage='leop.mosaic.32arcsec.pvslice.major.' + \
                         fileNums[i] + '.fits',
               velocity=True)

for i in xrange(0,9):
    exportfits(imagename='leop.mosaic.32arcsec.pvslice.minor.' + fileNums[i],
               fitsimage='leop.mosaic.32arcsec.pvslice.minor.' + \
                         fileNums[i] + '.fits',
               velocity=True)



