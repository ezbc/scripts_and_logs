#!/usr/bin/python

# import modules
from agpy.gaussfitter import multigaussfit as gfit
import pyfits as pf
import numpy as np
from mycoords import make_velocityAxis
import gaussFitting
import hianalysis
import pickle

# gridded cube in miriad using the following commands:
'''
fits in=leop.fluxRescale.16arcsec.vel.fits op=xyin \
    out=leop.cube.fluxRescale.mir

fits in=leop.16arcsec.vel.fits op=xyin out=leop.cube.mir

imbin in=leop.cube.fluxRescale.mir region='box(200,200,400,400)' \
    bin=4,4,4,4,1,1 out=leop.cube.fluxRescale.bin.mir

imbin in=leop.cube.mir \
    bin=4,4,4,4,1,1 out=leop.cube.bin.mir

fits in=leop.cube.fluxRescale.bin.mir op=xyout \
    out=leop.fluxRescale.16arcsec.vel.bin5arcsec.fits

fits in=leop.cube.bin.mir op=xyout \
    out=leop.16arcsec.vel.bin5arcsec.fits

rm -rf *.mir
'''

# Which computer?
lappy = False
cosmos = True

# Get fits cube without residual scaling to calculate noise
fitsfile = 'leop.16arcsec.vel.bin5arcsec.fits'
if lappy:
    fitsDir = '/home/elijah/research/leop/data/vla.gmrt/fitsImages/'
if cosmos:
    fitsDir = '/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/'
f = pf.open(fitsDir + fitsfile)
noiseCube = f[0].data

# Get fits cube
fitsfile = 'leop.fluxRescale.16arcsec.vel.bin5arcsec.fits'
if lappy:
    fitsDir = '/home/elijah/research/leop/data/vla.gmrt/fitsImages/'
if cosmos:
    fitsDir = '/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/'
f = pf.open(fitsDir + fitsfile)
header, cube = f[0].header, f[0].data

# Create subsection of cube
#subCube = cube[:,400:600,400:600]
#subCube = cube[:,450:550,450:550]
subCube = cube

# Define characteristics of cube
std = noiseCube.std() # standard deviation of cube without rescaling
velocities = make_velocityAxis(header) # writes velocities in km/s
velWidth = header['CDELT3']/1000.
beamsize = 16 # arcsec
cellsize = 6. # arcsec / pix
raSize = len(subCube[0,:])
decSize = len(subCube[0,0,:])
size = raSize * decSize
count = 0

# Set threshold
threshold = 5 # if max value of profile is above threshold*std then fit

# Perform monte carlo
simulationNum = 1000
simMeans = np.zeros(simulationNum)
simMedians = np.zeros(simulationNum)
simStds = np.zeros(simulationNum)
count = 0

for k in range(simulationNum):
    profileFits = []

    # Cube with new noise convolved with PSF
    newCube = hianalysis.addNoise2Cube(subCube, beamsize/cellsize, std=std)

    print 'Fitting cube ' + str(count) + ' of ' + str(simulationNum)
    # Fit each profile in cube with one Gaussian
    for i in range(raSize):
        for j in range(decSize):
            profile = newCube[:,i,j]
            if profile.max() > threshold*std:
                peakVal = profile.max()
                peakPos = np.where(profile==peakVal)[0]
                peakVelocity = velocities[peakPos][0]
                width = gaussFitting.guess_width(velocities,profile,peakPos)
                if peakVelocity + width < velocities.max() and \
                        peakVelocity - width > velocities.min():
                    if width < 2*velWidth:
                        width = 2*velWidth
                    elif width > 30:
                        width = 30
                    ncomp = 1
                    guesses = [profile.max(),peakVelocity,width]
                    fit = gfit(velocities,
                               profile,
                               ngauss=ncomp,
                               err=std,
                               params=guesses,
                               limitedmin=[True,True,False]*ncomp,
                               limitedmax=[True,True,True]*ncomp,
                               minpars=[std,velocities.min(),1*velWidth],
                               maxpars=[1.1*profile.max(),velocities.max(),30])
                    profileFits.append(fit)

    # calculate Gaussian width properties
    singleWidths = []
    for i in range(len(profileFits)):
        if profileFits[i][0][0] >= 5*std:
        #if True:
            singleWidths.append(profileFits[i][0][2])
        simMeans[k] = np.mean(singleWidths)
        simMedians[k] = np.median(singleWidths)
        simStds[k] = np.std(singleWidths)
    count += 1




#===============================================================================
# Histogram Plot: means
#===============================================================================

import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


figureDir = '/d/bip3/ezbc/leop/figures/'
figureName = 'leop.velDispersionFits.histogram'

# Now calculate histogram
hist, bin_edges = np.histogram(simMeans, density=True, bins=50)
hist /= len(hist)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
p0 = [0.1,8,0.2]

coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

# Get the fitted curve
fit_centres = np.arange(7,10,0.001)
hist_fit = gauss(fit_centres, *coeff)

# set up plot characteristics
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
size = 7
fig_size = [size*golden_mean*2,size]

fontScale = 12
params = {#'backend': 'eps',
          'axes.labelsize': fontScale,
          'text.fontsize': fontScale,
          'legend.fontsize': fontScale*3/4,
          'xtick.labelsize': fontScale,
          'ytick.labelsize': fontScale,
          'text.usetex': True, # NOTE: LaTeX prevents the use of tight_layout
          'figure.figsize': fig_size,
          }
plt.rcdefaults()
plt.rcParams.update(params)
fig, axes = plt.subplots(nrows=1, ncols=2, sharey=False)

axes[0].plot(bin_centres, hist, color='k',
        drawstyle='steps-pre')
axes[0].plot(fit_centres, hist_fit, color='k')
axes[0].annotate('(a)',
            (0.05,0.95),
             xycoords='axes fraction',
             textcoords='axes fraction')
axes[0].set_ylabel('Probability Density')
axes[0].set_xlabel('Mean Gaussian Width (km/s)')
axes[0].set_xlim(7.95,8.2)
#axes.legend(loc='upper right')

#===============================================================================
# Histogram Plot: stds
#===============================================================================

# Now calculate histogram
hist, bin_edges = np.histogram(simStds, density=True, bins=50)
hist /= len(hist)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2

# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
p0 = [0.1,1,0.2]

coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

# Get the fitted curve
fit_centres = np.arange(0,10,0.001)
hist_fit = gauss(fit_centres, *coeff)

# Now plot
axes[1].plot(bin_centres, hist, color='k',drawstyle='steps')
axes[1].plot(fit_centres, hist_fit, color='k')
axes[1].annotate('(b)',
               (0.05,0.95),
               xycoords='axes fraction',
               textcoords='axes fraction')
axes[1].set_ylabel('Probability Density')
axes[1].set_xlabel('Standard Deviation of Gaussian Widths (km/s)')
axes[1].set_xlim(1.3,1.42)

fig.savefig('/usr/users/ezbc/Dropbox/research/manuscripts/leop/paper5/' + \
            'post_circulation/finalImages/' + figureName + '.eps',dpi=600)
fig.savefig(figureDir + figureName + '.eps',dpi=600)
fig.show()
























