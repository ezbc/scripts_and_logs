################################################################################
################################################################################
# Script for plotting mbary vs. vcirc of dwarfs, spirals + leop
################################################################################
################################################################################

import numpy as np # make sure vector and array arithmetic options are loaded
import matplotlib.pyplot as plt
import pyfits as pf
import math
from pylab import arange,pi,sin,cos,sqrt

#data = np.genfromtxt('/home/elijah/research//leop/data/obsCombine/data/'+
#                     'baryVsVcirc.mcGaugh.dat',delimiter='')

rawData = np.genfromtxt('/d/bip3/ezbc/leop/data/hi/mcgaugh/' +
                     'baryVsVcirc.mcGaugh.dat',delimiter='')
data = np.zeros(shape=(rawData.T.shape))

#gasdata = np.genfromtxt('/home/elijah/research/leop/data/obsCombine/data/'+
#                     'gasrichdatatable.txt',delimiter='')

gasdata = np.genfromtxt('/d/bip3/ezbc/leop/data/hi/mcgaugh/' +
                     'gasrichdatatable.mcgaugh.txt',delimiter='')

gasdata = gasdata.T

for j in xrange(data.shape[0]):
    for i in xrange(len(data[1])):
        data[j,i] = 10 ** rawData.T[j,i]
        #data[j,i] = rawData.T[j,i]


gasMb = 10**gasdata[4] + 10**gasdata[6]
gasVc = gasdata[2]

Vcerr = np.array([data[0]-data[2],-(data[0]-data[3])])
Mberr = np.array([data[1]-data[4],-(data[1]-data[5])])

gasVcerr = gasdata[3]
gasMberr = np.sqrt((np.log(10)*gasdata[5]*10**gasdata[4])**2 + \
        (np.log(10)*gasdata[7]*10**gasdata[6])**2)

# Leo P calculations
distErrP = 0.12
distErrM = 0.4
dist = 1.72
mhi = 9.15e5
S = 1.31
SErr = 0.1*1.31
mstar = 5.7e5
mstarErr = 0.05*5.9e5

# the eq for hi mass: 2.36e5 * D^2 * S
hiErrP = 2.36e5 * np.sqrt( 4* dist**2 * S**2 * distErrP**2 + \
        dist**2 * SErr**2)

hiErrM = 2.36e5 * np.sqrt( 4* dist**2 * S**2 * distErrM**2 + \
        dist**2 * SErr**2)

mtot = 1.35 * mhi + mstar

mtotErrP = 1.35 * np.sqrt((hiErrP/mhi)**2 + (mstarErr / mstar)**2) * mtot
mtotErrM = 1.35 * np.sqrt((hiErrM/mhi)**2 + (mstarErr / mstar)**2) * mtot


for i in range(3):
    # Now Plot
    fig = plt.figure()

    fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inch
    golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height = fig_width*golden_mean      # height in inches
    size = 8
    fig_size = [size*golden_mean,size]

    fontScale = 12
    params = {'backend': 'eps',
              'axes.labelsize': fontScale,
              'text.fontsize': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
              'text.usetex': False,
              'figure.figsize': fig_size
              }

    plt.rcParams.update(params)
    plt.clf()

    plt.subplot(1,1,1)
    plt.subplots_adjust(wspace=0.3)

    markerSize = 6
    # Spirals
    yMinPos=25
    yMaxPos=82

    plt.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
                 xerr=(Vcerr[0,yMinPos:yMaxPos],Vcerr[1,yMinPos:yMaxPos]),
                 yerr=(Mberr[0,yMinPos:yMaxPos],Mberr[1,yMinPos:yMaxPos]),
                       marker='o',color='r',ecolor='k',linestyle='none',
                       markersize=markerSize)
    if i == 0:
        # Dwarfs
        yMinPos=83
        yMaxPos=114
        plt.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
                     xerr=(Vcerr[0,yMinPos:yMaxPos],Vcerr[1,yMinPos:yMaxPos]),
                     yerr=(Mberr[0,yMinPos:yMaxPos],Mberr[1,yMinPos:yMaxPos]),
                           marker='s',color='b',ecolor='k',linestyle='none',
                           markersize=markerSize)

    # Gas-rich galaxies
    plt.errorbar(gasVc,gasMb,
                 xerr=gasVcerr,
                 yerr=gasMberr,
                 marker='^',color='g',ecolor='k',linestyle='none',
                 markersize=markerSize)

    # Plot line for fitting spirals: function = log(M_b) = x log(V_c) + A
    # Mb = V_c + (47 +/- 6)

    def line(array, slope):
        return slope * array**4

    vc_array = np.arange(0,1000,0.1)
    mb_array = line(vc_array, 47.)
    sigma = 3
    mb_array_high = line(vc_array,47 + sigma * 6.)
    mb_array_low = line(vc_array, 47 - sigma * 6.)

    plt.fill_between(vc_array, mb_array_low, mb_array_high, color='gray',
            alpha=0.7)

    plt.ylabel(r"$M_b$ (M$_{\odot}$)",family='serif')
    plt.xlabel(r"$v_c$ (km/s)",family='serif')
    plt.xscale('log')
    plt.yscale('log')
    if i < 4:
        plt.xlim([1,10**3])
        plt.ylim([50,10**13])
    else:
        plt.xlim([3,10**2.8])
        plt.ylim([10**5.6,10**12])
    plt.legend(loc='upper left')

    if i == 2:
        # Leo P
        plt.errorbar(16,mtot,
                     #xerr=([15,0]),yerr=([10**6,10**6]),
                     xerr=([5],[5]),yerr=([mtotErrM],[mtotErrP]),
                     xlolims=False,
                     marker='s',color='c',ecolor='k',elinewidth=2,capsize=3,
                     markersize=markerSize)

    if i == 0:
    	filename = 'leop_mcGaughPlot_dwarfs_samescale'
    elif i == 1:
    	filename = 'leop_mcGaughPlot_spirals_noLeoP_samescale'
    elif i == 2:
    	filename = 'leop_mcGaughPlot_spirals_LeoP_samescale'

    plt.savefig('/d/bip3/ezbc/leop/figures/' + filename + '.eps',dpi=600,
            bbox_inches='tight')
    plt.savefig('/d/bip3/ezbc/leop/figures/' + filename + '.png',dpi=600,
            bbox_inches='tight')

    plt.show()



