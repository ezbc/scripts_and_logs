# Script for plotting spectra HI and CO spectra in regions with
# integrated CO greater than 3sigma of the Taurus cloud


import numpy as np # make sure vector and array arithmetic options are loaded
import matplotlib.pyplot as plt
import pyfits as pf

# read files into python
dirPath = '/d/bip3/ezbc/perseus/data/'
galfaData = pf.getdata(dirPath + 'perseus.galfa.cube.cfaRes.regrid.fits')
galfaHead = pf.getheader(dirPath + 'taurus.galfa.cube.cfaRes.regrid.fits')
cfaCubeData = pf.getdata(dirPath + 'taurus.cfa.cube.fits')
cfaCubeHead = pf.getheader(dirPath + 'taurus.cfa.cube.fits')

# We want to find the positions on the CO map with Tb > 3 K and
# integrate the galfa spectra coincident with these positions. Run a
# double for loop. First loop iterates through pixels in cfa cube
# finding 3sigma (3 K) emission. Second loop extracts the spectrum
# from the galfa cube. Append the spectrum to a 2D array list.

#-------------------------------------------------------------------------------
# Write out velocity array
galfa_vDelt = galfaHead['CDELT3']
galfa_vPix = int(round(galfaHead['CRPIX3']))
galfa_vVal = galfaHead['CRVAL3']

galfaVelocities = np.zeros(shape=(galfaData.shape[0]))
channel = galfa_vPix
while channel < galfaVelocities.shape[0]:
    galfaVelocities[channel] = galfa_vVal + (channel-galfa_vPix)*galfa_vDelt
    channel = channel + 1

channel = galfa_vPix - 1
while channel >= 0:
    galfaVelocities[channel] = galfa_vVal + (channel-galfa_vPix)*galfa_vDelt 
    channel = channel - 1

# Initialize list
spectrumList = [galfaVelocities]
spectrumListNoCo = [galfaVelocities]
spectrumListTot = [galfaVelocities]

# Find 3sigma points in cfa coldens map
for i in xrange(0,cfaIntData.shape[0]-1):
    for j in xrange(0,cfaIntData.shape[1]-1):
        if cfaIntData[i,j] >= 5:
            spectrumList.append(galfaData[:,i,j])
        else:
            spectrumListNoCo.append(galfaData[:,i,j])
        spectrumListTot.append(galfaData[:,i,j])

# average together the galfa spectra
avgSpectrum = np.zeros(shape=(2,galfaVelocities.shape[0]))
avgSpectrum[0] = galfaVelocities
avgSpectrumNoCo = np.zeros(shape=(2,galfaVelocities.shape[0]))
avgSpectrumNoCo[0] = galfaVelocities
avgSpectrumTot = np.zeros(shape=(2,galfaVelocities.shape[0]))
avgSpectrumTot[0] = galfaVelocities

# Hi Spectrum Not Coincident with CO
tempArray = np.zeros(shape=(galfaVelocities.shape[0]))
nanCount = 0

for i in xrange(0,spectrumListNoCo[0].shape[0]-1):
    for j in xrange(1,len(spectrumListNoCo)-1): # omit velocity axis

        if np.isnan(spectrumListNoCo[j][i]):
            nanCount = nanCount + 1
        else:
            tempArray[i] = tempArray[i] + spectrumListNoCo[j][i]

    tempArray[i] = tempArray[i] / (len(spectrumListNoCo) - nanCount)
    nanCount=0

avgSpectrumNoCo[1] = tempArray

# Total HI Spectrum
tempArray = np.zeros(shape=(galfaVelocities.shape[0]))
nanCount = 0

for i in xrange(0,spectrumListTot[0].shape[0]-1):
    for j in xrange(1,len(spectrumListTot)-1): # omit velocity axis

        if np.isnan(spectrumListTot[j][i]):
            nanCount = nanCount + 1
        else:
            tempArray[i] = tempArray[i] + spectrumListTot[j][i]

    tempArray[i] = tempArray[i] / (len(spectrumListTot) - nanCount)
    nanCount=0

avgSpectrumTot[1] = tempArray

# HI Spectrum Coincident with CO
tempArray = np.zeros(shape=(galfaVelocities.shape[0]))
nanCount = 0

for i in xrange(0,spectrumList[0].shape[0]-1):
    for j in xrange(1,len(spectrumList)-1): # omit velocity axis

        if np.isnan(spectrumList[j][i]):
            nanCount = nanCount + 1
        else:
            tempArray[i] = tempArray[i] + spectrumList[j][i]

    tempArray[i] = tempArray[i] / (len(spectrumList) - nanCount)
    nanCount=0

avgSpectrum[1] = tempArray

#-------------------------------------------------------------------------------
# Repeat for cfa data
# Write out velocity array for cfa data
cfa_vDelt = cfaCubeHead['CDELT3']
cfa_vPix = int(round(cfaCubeHead['CRPIX3']))
cfa_vVal = cfaCubeHead['CRVAL3']

cfaVelocities = np.zeros(shape=(cfaCubeData.shape[0]))
channel = cfa_vPix
while channel < cfaVelocities.shape[0]:
    cfaVelocities[channel] = cfa_vVal + (channel-cfa_vPix)*cfa_vDelt
    channel = channel + 1

channel = cfa_vPix - 1
while channel >= 0:
    cfaVelocities[channel] = cfa_vVal + (channel-cfa_vPix)*cfa_vDelt 
    channel = channel - 1

cfaSpectrum = np.zeros(shape=(2,cfaVelocities.shape[0]))
cfaSpectrum[0] = cfaVelocities


# Initialize list
spectrumList = [cfaVelocities]

# Find 3sigma points in cfa coldens map
for i in xrange(0,cfaIntData.shape[0]-1):
    for j in xrange(0,cfaIntData.shape[1]-1):
        if cfaIntData[i,j] >=  5:
            spectrumList.append(cfaCubeData[:,i,j])

# average together the galfa spectra
cfaSpectrum = np.zeros(shape=(2,cfaVelocities.shape[0]))
cfaSpectrum[0] = cfaVelocities

tempArray = np.zeros(shape=(cfaVelocities.shape[0]))
nanCount = 0

for i in xrange(0,spectrumList[0].shape[0]-1):
    for j in xrange(1,len(spectrumList)-1): # omit velocity axis

        if np.isnan(spectrumList[j][i]):
            nanCount = nanCount + 1
        else:
            tempArray[i] = tempArray[i] + spectrumList[j][i]

    tempArray[i] = tempArray[i] / (len(spectrumList) - nanCount)
    nanCount=0

# CfA antenna temperature Ta relates to brightness temperature as:
# Tb = Ta/n where n is the antenna efficiency, n = 0.84 +/- 0.04
cfaSpectrum[1] = tempArray/0.82



plt.clf() # clear the plot (figure)
plt.plot(avgSpectrum[0,:]/1000,avgSpectrum[1,:], 
         'r',linestyle='dashed',label="GALFA HI, CO Coincident")
plt.plot(avgSpectrumNoCo[0,:]/1000,avgSpectrumNoCo[1,:], 
         color='g',linestyle='dotted',label="GALFA HI, CO Incoincident",linewidth=3)
#plt.plot(avgSpectrumTot[0,:]/1000,avgSpectrumTot[1,:], 
#         color='k',linestyle='dashed',label="GALFA HI, Total")
plt.plot(cfaSpectrum[0,:],cfaSpectrum[1,:]*100, 
         'b', label="CfA CO x 100")
plt.xlabel("Velocity (km/s)")
plt.ylabel("Average Tb (K km/s)")
plt.xlim([-100,50])
plt.ylim([-10,250])
plt.legend(loc='upper left')
plt.savefig('spectra.5sigma.ps',dpi=500)
plt.show()
