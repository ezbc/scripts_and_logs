################################################################################
################################################################################
# Script for plotting HI and stellar mass profiles
################################################################################
################################################################################

import numpy as np # make sure vector and array arithmetic options are loaded
import matplotlib.pyplot as plt
import pyfits as pf
import math
from scipy.optimize import curve_fit


ellintDir = '/d/bip3/ezbc/leop/data/hi/gipsy/ellint/'
rotcurDir = '/d/bip3/ezbc/leop/data/hi/gipsy/rotcur/'

hiData32 = np.genfromtxt(
    ellintDir + 'leop.surfBrightProfile.32arcsec.inc65',
    delimiter='',
    comments='!').T

hiData16 = np.genfromtxt(
    ellintDir + 'leop.surfBrightProfile.16arcsec.inc65',
    delimiter='',
    comments='!').T

# Load the velocity dispersion profiles

vdispData32 = np.genfromtxt(
    ellintDir + 'leop.vdispProfile.32arcsec.inc65',
    delimiter='',
    comments='!').T

vdispData16 = np.genfromtxt(
    ellintDir + 'leop.vdispProfile.16arcsec.inc65',
    delimiter='',
    comments='!').T
'''
# Now load rotation curves

vrotData32 = np.genfromtxt(
    rotcurDir + 'leop.rotcur.32arcsec.inc65',
    delimiter='',
    skip_header=11).T

vrotData16 = np.genfromtxt(
    rotcurDir + 'leop.rotcur.16arcsec.inc65',
    delimiter='',
    skip_header=11).T
'''

################################################################################
# Begin fitting of profiles
################################################################################


# Header of hi profiles are as follows:
# radius radius    Msd-t   Mass-t Cum.Mass      Msd     Mass Cum.Mass      sum subpix     area  area-bl  segLO  segHI  width     pa 
# arcsec    Kpc  Mo/pc^2  10^9 M0  10^9 M0  Mo/pc^2  10^9 M0  10^9 M0 JY/BEAM.KM/S      #   pixels   pixels   deg.   deg. arcsec   d

# Header of vdisp profiles are as follows:
# radius   cum.sum       sum      mean       rms subpix unique      area   area-bl  segLO  segHI  width     pa   incl
# arcsec      KM/S      KM/S      KM/S      KM/S      # pixels    pixels    pixels   deg.   deg. arcsec   deg.   deg.

# Header of rotation curve:
#   radius    width systemic  error rotation  error expansion  error    pos.  error incli-  error x-pos.  error y-pos.  error npts   sigma
#                   velocity        velocity         velocity          angle        nation        centre        centre             velocity
# (arcsec) (arcsec)   (km/s) (km/s)   (km/s) (km/s)    (km/s) (km/s)  (deg.) (deg.) (deg.) (deg.) (grid) (grid) (grid) (grid)        (km/s)

# Get radii
# Same for each resolution
radii = hiData16[0]

# we want msd from the hi
hi16 = hiData16[5]
hiErradii = 0.1*hiData16[5]

vdisp16 = vdispData16[3]
vdispErradii = vdispData16[4]

# Define a decaying exponential function
def expo(x, a, b, c):
    return a * np.exp(-b * x) + c

# Define asymmetric drift function
def calc_pressure(r,hi,vdisp,hiErr=None,vdispErr=None):
    coeff, var_matrix = curve_fit(expo, r, hi, maxfev=int(1e8))
    hiFit = expo(r, coeff[0], coeff[1], coeff[2])

    if type(vdisp) is np.ndarray:
        coeff, var_matrix = curve_fit(expo, r, vdisp, maxfev=int(1e8))
        vdispFit = expo(r, coeff[0], coeff[1], coeff[2])
    else:
        vdispFit = np.ones(hi.shape)*vdisp

    hiFit = hi
    #vdispFit = vdisp

    # Initialize arrays for creating derivatives of the natural logs
    loghiFit = np.log(hiFit)
    logvdispFit = np.log(vdispFit**2)
    logr = np.log(r)

    # Find slope between each data point
    dloghiFit = np.gradient(loghiFit)
    dlogvdispFit = np.gradient(logvdispFit)
    dlogvdispFit = np.zeros(vdispFit.shape)
    dlogr = np.gradient(logr)

    # Initialize pressure component array
    pressureComp = np.zeros(hi.shape[0]-1)

    # Define pressure correction
    for i in xrange(0,hi.shape[0]-1):
       pressureComp[i] =  vdispFit[i] * \
          (dloghiFit[i] + dlogvdispFit[i] )

    dloghiFitErr = hiErr/hi
    dlogvdispFitErr = 2*vdispErr/vdisp
    print dlogvdispFitErr
    hiErrFit = np.gradient(dloghiFit) * dloghiFitErr

    if type(vdisp) is np.ndarray:
        vdispErrFit = np.gradient(dlogvdispFit) * dlogvdispFitErr
    elif type(vdisp) is float or type(vdisp) is int:
        vdispErrFit = np.ones(hi.shape)*vdispErr

    err = np.sqrt((dlogvdispFitErr/dlogvdispFit)**2 + \
            (dloghiFitErr/dloghiFit)**2) 

    err = np.sqrt((dlogvdispFitErr)**2 + \
            (dloghiFitErr/dloghiFit)**2) 
    return pressureComp, err

def rotcurve(vflat, iflat, radii, vflatErr, iflatErr):
    rotcurve = vflat * (1. - np.exp(-radii/iflat))
    rotcurveErr = ((1-np.exp(-radii/iflat)*vflatErr)**2 +\
        (vflat*radii*np.exp(-radii/iflat)*iflatErr/iflat**2)**2)**0.5
    return rotcurve,rotcurveErr

rotation_vel,rotation_vel_error = rotcurve(3,60,radii,3,32)

################################################################################
# Plot
################################################################################

from mpl_toolkits.axes_grid1 import Grid

#fig = plt.figure(figsize=(7,3))
fig = plt.figure()

size = 7
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_size = [7,3]

fontScale = 12
params = {'backend': 'eps',
          'axes.labelsize': fontScale,
          'text.fontsize': fontScale,
          'legend.fontsize': fontScale*3/4,
          'xtick.labelsize': fontScale,
          'ytick.labelsize': fontScale,
          #'text.usetex': True,
          'figure.figsize': fig_size
          }

plt.rcdefaults()
plt.rcParams.update(params)
plt.figure(1)
#plt.tight_layout(pad=0.0)
plt.clf()
plt.subplots_adjust(bottom=0.2)

grid = Grid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 direction='column',
                 axes_pad=0,
                 label_mode='L',
                 share_all=True)

rNum = 5

pressureComp,error = calc_pressure(radii[0:rNum],hi16[0:rNum],8.4,\
        hiErr=hiErradii[0:rNum],vdispErr=1.4)

markerSize = 5

# Pressure Component
if True:
    plt.errorbar(radii[0:rNum-1],pressureComp[0:rNum],
        label=r"$v_{\rm pressure}$",
        yerr=error[0:rNum-1],
        capsize=5,
        elinewidth=1.5,
        ecolor='k',
        marker='s',
        linestyle='.',
        color='b',
        markersize=markerSize)

# Circular Velocity
plt.errorbar(radii[0:rNum-1],rotation_vel[0:rNum-1] - pressureComp[0:rNum],
        color='k',
        yerr=np.sqrt(rotation_vel_error[0:rNum-1]**2 + error[0:rNum-1]**2),
        capsize=5,
        elinewidth=1.5,
        ecolor='k',
        label="$v_{c}$",
        linestyle='-.',
        marker='^',
        markersize=markerSize)

# Last point is 5.8 km/s...wow that's slow!

if True:
    # Rotational Velocity with error bars
    plt.errorbar(radii[0:rNum-1],rotation_vel[0:rNum-1],
            marker='o',
            linestyle='--',
            color='r',
            label='$v_{rot}$',
        markersize=markerSize,
        yerr=rotation_vel_error[0:rNum-1],
        capsize=5,
        elinewidth=1.5,
        ecolor='k',)

# Global plot properties
plt.xlim([0,70])
#plt.ylim([-7,16])

plt.xlabel('Radius (\'\')',family='serif')
plt.ylabel('Velocity (kms$^{-1}$)',family='serif')

plt.xticks((0,16,32,48,64))

#letters = ['(a)','(b)','(c)']
#plt.annotate(letters[i],
#            (0.05,0.9)
#            ,color='k'
#            ,textcoords='axes fraction',
#            xycoords='axes fraction')

plt.legend(loc='upper left')
plt.savefig('../figures/leop.rotationCurve.eps',dpi=600)
plt.show()



