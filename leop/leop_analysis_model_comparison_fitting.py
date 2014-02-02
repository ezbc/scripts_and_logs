#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Script for comparing model cubes with line profiles of 16" cube
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# INSTRUCTIONS (as given by Ed Elson):
#-------------------------------------
# 1. From your HI data cube pick out the line profiles that have a
#peak flux above 2 sigma.  So, from a line-free channel in your cube
#you'll first have to generate a good, robust estimate of the noise
#RMS.

#1.5 Alternatively, pick out line profiles in your cube that have at
#least 50% (or some other fraction) of their flux above 1 sigma or 2
#sigma.  So basically, you're looking for line profiles that contain a
#fair amount of signal, line profiles you can trust.

#2. Then calculate the (absolute?) residual between each of these line
#profiles in the data cube and the line profiles in the model.  Do you
#understand what I mean?  John, this is what I was saying we should do
#for the SHIELD galaxies.  Although I haven't tried all of this
#myself, I imagine it would make for fairly good data-model
#comparisons.

#3. You might consider parameterising each of your "trustable"
#profiles from the data cube, using a Gaussian.  Then compare those
#fitted Gaussians to the model profiles.  I say this because Gaussians
#are used to construct those models.  This method might make for
#cleaner comparisons.
#-------------------------------------

# Import libraries:
import numpy as np
import math
import pyfits as pf
from mycoords import make_velocityAxis
from agpy.gaussfitter import multigaussfit as gfit
import gaussFitting
import os

# Define paths:
fitsDir = '/d/bip3/ezbc/leop/data/hi/models/modelCreationFolder5/models/'
writeDir = '/d/bip3/ezbc/leop/data/hi/models/data/'
psDir =  '/d/bip3/ezbc/leop/figures/'

# Open data cube with flux scaling
f = pf.open('/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/' + \
            'leop.fluxRescale.16arcsec.vel.modelBin.fits')
h, d = f[0].header, f[0].data[:,:,:]

velArray = make_velocityAxis(h)

# Compute the final noise from data cube without flux scaling
def get_rms(x, axis=None):
    return np.sqrt(np.mean(x**2, axis=axis))

fnoise = pf.open('/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/' +\
                 'originalCubes/leop.16arcsec.fits')
hnoise,dnoise = fnoise[0].header, fnoise[0].data[:,:,:]
rms = get_rms(dnoise)
std = dnoise.std()
threshold = 3. * rms

# Define number of models being used:
models = [name for name in os.listdir(fitsDir) if name.endswith('.FITS')]
modelNum = len(models)
print 'Number of models: ' + str(modelNum)

'''
########### Gaussian Setup #############
# Define normalization parameters so that gaussian can be fit.

# Create cube for gaussian profiles
dGauss = np.zeros(d.shape)

########### Create convolved spectra #############
# Create array of integers from 1 to 70 to be the x-values for the
# gaussian fit. The velocities are too high 

# Establish which pixels are significant in the data cube:
# Create profile of each pixel in the data cube. First axis is frequency.
for j in xrange(0,d.shape[1]): # RA axis
    for k in xrange(0,d.shape[2]): # Dec axis

        # check if the sum of profile > threshold
        profile = d[:,j,k]
        fluxSum = sum(profile[profile == profile])
        if fluxSum > threshold and profile.max() > threshold:

            # Makes guesses
            peakVal = profile.max()
            peakPos = np.where(profile==peakVal)[0]
            peakVelocity = velArray[peakPos][0]
            width = gaussFitting.guess_width(velArray,profile,peakPos)
            if width > 20: width = 20
            guesses = [profile.max(),peakVelocity,width]
            #print guesses

            # perform fit
            fit = gfit(velArray,
                       profile,
                       ngauss=1,
                       err=std,
                       params=guesses,
                       limitedmin=[True,True,False],
                       limitedmax=[True,True,True],
                       minpars=[std,velArray.min(),0],
                       maxpars=[1.1*profile.max(),velArray.max(),20])

            # Write the gaussian profile to a new cube
            dGauss[:,j,k] = fit[1]

    if j%(d.shape[1]/10) == 0:
        print str(100./d.shape[1] * (j+1)) + '% done'
'''
# ==============================================================================
print 'Now calculating residuals...'
# ==============================================================================

# Define empty array for residual sum of each model
gofArray = np.zeros(modelNum)
residArray = np.zeros(modelNum)
gof_reduced_array = np.zeros(modelNum)

count = 0
fluxSum = 0.93 * 1.13*16**2/1.5**2 # for normalizing model fluxes
dataSum = d.sum()
dataError = d *0.1 + std # for calculating goodness of fit
for i in xrange(0,modelNum):
    # Load model cube
    if i < 32767:
        fmod = pf.open(fitsDir + 'model_' + str(i+1) + '.smooth.FITS')
        hmod, dmod = fmod[0].header, fmod[0].data
    elif i >= 32767:
        fmod = pf.open(fitsDir + 'model_-' + str(32768-count) + '.smooth.FITS')
        hmod, dmod = fmod[0].header, fmod[0].data
        count = count + 1
    # Create model for normalizing flux
    dmodNorm = dmod
    dmodNorm[dmodNorm != dmodNorm] = 0. # Set NaNs to 0, otherwise sum will fail
    dmodNorm = dmodNorm / dmodNorm.sum() * dataSum # normalize flux
    # calculate goodness of fit
    numerator = np.subtract(d, dmodNorm)**2
    denominator = dataError**2 # d.size = degrees of freedom
    gof = np.divide(numerator,denominator).sum()
    gof_reduced = gof / (d.size - 5.)

    if i%100 == 0:
        print 'Chi-squared calculated for the ' + str(i) + 'th image'
    gofArray[i] = gof
    gof_reduced_array[i] = gof_reduced

# Write residual array
np.savetxt(writeDir + 'models.gof.v5.txt',gofArray)
np.savetxt(writeDir + 'models.gof_reduced.v5.txt',gof_reduced_array)

gofArray = np.genfromtxt(writeDir + 'models.gof.v5.txt',delimiter='')
gof_reduced_array = np.genfromtxt(writeDir + 'models.gof_reduced.v5.txt',
                                  delimiter='')

# gof closest to 1 are best-fits
gof_reduced_array = abs(gof_reduced_array - 1)
sortIndices = gofArray.argsort()

# Read in parameter file for displaying model properties
paramDir = '/d/bip3/ezbc/leop/data/hi/models/modelCreationFolder5/' + \
           'parameters.16arcsec.txt'
paramData = np.genfromtxt(paramDir,delimiter='',usecols=(1,2,3,4,5,6))
modelNameData = np.genfromtxt(paramDir,delimiter='',usecols=(0),dtype=str)
# Data repeats each line, choose only unique lines
params = paramData[0::2]
modelNames = modelNameData[0::2]
# Sort Lines
modelNamesSorted = modelNames[sortIndices]
paramsSorted = params[sortIndices]
gofArraySorted = gofArray[sortIndices]

# Create a text file for writing out the best model data parameters
fileName = 'modelComparisonData.v5.txt'

file = open(writeDir + fileName,'wb')

# Write column headers
file.writelines(['model #  ','resid val \t','inc \t',
                 'pa \t','iflat \t','vflat \t','vdisp \t','z0 \n'])

# modelNum column data
for i in xrange(0,modelNum):
    file.writelines([str(modelNamesSorted[i]) + '\t',
                     str(gofArraySorted[i]) + '\t',
                     str(paramsSorted[i,0]) + '\t',
                     str(paramsSorted[i,1]) + '\t',
                     str(paramsSorted[i,2]) + '\t',
                     str(paramsSorted[i,3]) + '\t',
                     str(paramsSorted[i,4]) + '\t',
                     str(paramsSorted[i,5]) + '\n'])

file.close()

print 'File "' + fileName + '" has been written to ' + writeDir





