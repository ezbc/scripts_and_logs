# Script for plotting spectra HI and CO spectra in regions with
# integrated CO greater than 3sigma of the Taurus cloud


import numpy as np # make sure vector and array arithmetic options are loaded
import matplotlib.pyplot as plt
import pyfits as pf

# read files into python
dirPath = '/d/bip3/ezbc/taurus/data/'
galfaData = pf.getdata(dirPath + 'galfa/taurus.galfa.cube.cfaRes.regrid.fits')
galfaHead = pf.getheader(dirPath + 'galfa/taurus.galfa.cube.cfaRes.regrid.fits')
cfaIntData = pf.getdata(dirPath + 'cfa/taurus.cfa.co_Int.fits')
cfaIntHead = pf.getheader(dirPath + 'cfa/taurus.cfa.co_Int.fits')
cfaCubeData = pf.getdata(dirPath + 'cfa/taurus.cfa.cube.fits')
cfaCubeHead = pf.getheader(dirPath + 'cfa/taurus.cfa.cube.fits')

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

################################################################################
# CO Spectrum
################################################################################

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
cfaSpectrum[1] = tempArray / 0.82

################################################################################
# HI Spectrum Coincident with CO
# (Includes gaussian fitting to components)
################################################################################

# Create spectrum
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



"""
# Fit spectrum
from agpy.gaussfitter import multigaussfit as gaussfit

# Initial guesses:
# params - Fit parameters: [amplitude, offset, width] * ngauss
# Visit the following url for Adam Ginsburg's gaussfitter:
# http://agpy.googlecode.com/svn/trunk/agpy/gaussfitter.py

# Note, ngauss parameter has a bug, need additional ngauss count for # of
# gaussians, e.g., for three gaussian components, ngauss=4

def modelData(x,amp,offset,width):
    profile = np.zeros(len(x))
    for i in xrange(len(x)):
        profile[i] = amp*np.exp(-(x[i]-offset)**2/(2*width**2))
    return profile

nGauss = 4
initialGuess = [5,-50,10,15,-10,20,45,5,20] 

fitData = multiGaussFit(avgSpectrum[0]/1000.,avgSpectrum[1],
                   ngauss=nGauss,params=initialGuess)

for i in xrange(nGauss):
comp = modelData(avgSpectrum[0]/1000,
                  fitData[0][0],fitData[0][1],fitData[0][2])

comp1 = modelData(avgSpectrum[0]/1000,
                  fitData[0][3],fitData[0][4],fitData[0][5])
comp2 = modelData(avgSpectrum[0]/1000,
                  fitData[0][6],fitData[0][7],fitData[0][8])

comp3 = modelData(avgSpectrum[0]/1000,
                  fitData[0][9],fitData[0][10],fitData[0][11])

test1 = gaussfit(avgSpectrum[0]/1000.,avgSpectrum[1],
                   ngauss=1,params=initialGuess)

################################################################################
# Now Plot
################################################################################

plt.clf() # clear the plot (figure)
plt.plot(avgSpectrum[0,:]/1000,avgSpectrum[1,:], 
         'b',label="GALFA HI, CO Coincident")

#plt.plot(cfaSpectrum[0,:],cfaSpectrum[1,:]*10, 
#         'b',label="CfA CO x 10")

plt.plot(avgSpectrum[0]/1000,fitData[1],
         color='k',
         linestyle='dotted',
         linewidth=3,
         label="Model Gaussian Fit, 3 components")

colors = ['g','r','c','m','y','k','b']
count = 0
compNum = 0

for i in xrange(nGauss):
    comp = modelData(avgSpectrum[0]/1000,
                     fitData[0][i+count],
                     fitData[0][i+count+1],
                     fitData[0][i+count+2])
    print fitData[0][i+count]
    if comp.max() != 0:
        print 'yes!'
        plt.plot(avgSpectrum[0]/1000,comp, 
                 color=colors[i],linestyle='dashed',
                 label="Model Gaussian Fit, comp #" + str(compNum))
        compNum = compNum + 1
    count = count + 2


plt.rc('legend',**{'fontsize':12})
plt.legend(loc='upper left')
plt.xlabel("Velocity (km/s)")
plt.ylabel("Average Tb (K km/s)")
plt.xlim([-100,50])
plt.ylim([-10,60])
plt.savefig('spectra.5sigma.ps',dpi=500)
plt.show()

"""
