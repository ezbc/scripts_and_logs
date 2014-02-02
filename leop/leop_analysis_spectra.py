###################################################################################
###################################################################################
# Script for plotting VLA/C + VLA/D + GMRT spectra of Leo P with
# ALFALFA and LBW spectra
###################################################################################
###################################################################################

# Directions for exporting a spectra from a cube in CASA can be found here:
# http://casaguides.nrao.edu/index.php?title=MRK_6:_red-shifted_HI_absorption#Save_the_spectrum_to_a_text_file

# Beam Area = 102.261 pixels

# First lets make a cube with units of mJy instead of mJy / beam

import numpy as np # make sure vector and array arithmetic options are loaded
import matplotlib.pyplot as plt
import pyfits as pf
import math
from pylab import setp
import os

fudge = 0

#path = '/home/elijah/research/leop/data/spectra/'
path = '/d/bip3/ezbc/leop/data/hi/profiles/'
figureDir = '/d/bip3/ezbc/leop/figures/'
finalImagesDir = '/usr/users/ezbc/Dropbox/research/manuscripts/' + \
    'leop/paper5/post_circulation/finalImages/'


vla4 = np.genfromtxt(path + "leop.vla.gmrt.profile.4arcsec.txt",delimiter='')
vla8 = np.genfromtxt(path + "leop.vla.gmrt.profile.8arcsec.txt",delimiter='')
vla16 = np.genfromtxt(path + "leop.vla.gmrt.profile.16arcsec.txt",delimiter='')
vla32 = np.genfromtxt(path + "leop.vla.gmrt.profile.32arcsec.txt",delimiter='')
vla32noise = np.genfromtxt(path + \
        "leop.vla.gmrt.noiseprofile.32arcsec.txt",delimiter='')
lbwSpectrum = np.genfromtxt(path + "leop.lbw.profile.csv",delimiter='')
alfaSpectrum = pf.getdata(path + 'leop.alfalfa.profile.fits')

vla4 = vla4.T
vla4[1,:] = vla4[1,:] * 1000.

vla8 = vla8.T
vla8[1,:] = vla8[1,:] * 1000.

vla16 = vla16.T
vla16[1,:] = vla16[1,:] * 1000.

vla32 = vla32.T
vla32[1,:] = vla32[1,:] * 1000.

# Compute rms of 32"
def get_rms(x, axis=None):
    return np.sqrt(np.mean(x**2, axis=axis))
vla32rms = get_rms(vla32noise[:,1])*1000.

# Convert the ALFALFA observations to the Local standard of rest frame
l = math.radians(219.654) # galactic longitude
b = math.radians(54.430)  # galactic latitude
lsrCoeff = 9 * math.cos(l) * math.cos(b) + 12 * math.sin(l) \
    * math.cos(b) + 7 * math.sin(b)
lbwSpectrum = lbwSpectrum.T
lbwSpectrum[1,:] = lbwSpectrum[1,:] + lsrCoeff - fudge
lbwSpectrum[2,:] = lbwSpectrum[2,:] * 1000
alfaSpectrum[0,:] = alfaSpectrum[0,:] + lsrCoeff - fudge


# Plot

fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
#fig_size =  [fig_width,fig_height]
size = 7
fig_size = [size*golden_mean*2,size]
#fig_size = [3,4]

fontScale = 12
params = {#'backend': 'eps',
          'axes.labelsize': fontScale,
          'text.fontsize': fontScale,
          'legend.fontsize': fontScale*3/4,
          'xtick.labelsize': fontScale,
          'ytick.labelsize': fontScale,
          #'text.usetex': True, # NOTE: LaTeX prevents the use of tight_layout
          'figure.figsize': fig_size,
          }

plt.rcdefaults()
plt.rcParams.update(params)

fig, axes = plt.subplots(nrows=1, ncols=2, sharey=True)
fig.subplots_adjust(
    left = 0.2,      # the left side of the subplots of the figure
    right = 1,   # the right side of the subplots of the figure
    bottom = 0.0,  # the bottom of the subplots of the figure
    top = 1,    # the top of the subplots of the figure
    wspace = 0, # width reserved for blank space between subplots
    hspace = 0,    # height reserved for white space between subplots
    )


# ALFA, LBW, and VLA+GMRT spectra
axes[0].plot(vla32[0,:],vla32[1,:], 'r',label="VLA+GMRT, 32\'\' beam")
axes[0].plot(lbwSpectrum[1,:],lbwSpectrum[2,:],'g-.', label="LBW")
axes[0].plot(alfaSpectrum[0,:],alfaSpectrum[1,:], 'b--', label="ALFALFA")
axes[0].set_xlabel("LSR Velocity (km/s)")
axes[0].set_ylabel("Flux Density (mJy)")
axes[0].set_xlim([231,295])
axes[0].set_ylim([-10,65])
axes[0].legend(loc='upper left')
axes[0].annotate('(a)',xy=(0.05,0.05),xycoords='axes fraction',
    textcoords='axes fraction')

# 2nd: VLA+GMRT spectra at each resolution
setp(axes[1].get_yticklabels(), visible=False)
axes[1].plot(vla4[0,:],vla4[1,:], 'b--',label="4\'\' beam")
axes[1].plot(vla8[0,:],vla8[1,:], 'g-.',label="8\'\' beam")
axes[1].plot(vla16[0,:],vla16[1,:], 'k-..',label="16\'\' beam")
axes[1].plot(vla32[0,:],vla32[1,:], 'r',label="32\'\' beam")
axes[1].set_xlabel("LSR Velocity (km/s)")
axes[1].set_xlim([231,295])
axes[1].set_ylim([-10,65])
axes[1].legend(loc='upper left')
axes[1].annotate('(b)',xy=(0.05,0.05),xycoords='axes fraction',
    textcoords='axes fraction')

# Save figure
filename='leop.spectra.vla.gmrt'
fig.savefig(figureDir + filename + '.eps',dpi=600,bbox_inches='tight')
fig.savefig(finalImagesDir + filename + '.eps',dpi=600,bbox_inches='tight')

plt.show()
