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

# Now Plot
fig = plt.figure(figsize=(8,9))

fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
#fig_size =  [fig_width,fig_height]
size = 7
fig_size = [size*golden_mean,size]

fontScale = 12
params = {'backend': 'eps',
          'axes.labelsize': fontScale,
          'text.fontsize': fontScale,
          'legend.fontsize': fontScale*3/4,
          'xtick.labelsize': fontScale,
          'ytick.labelsize': fontScale,
          'text.usetex': True,
          'figure.figsize': fig_size
          }

plt.rcParams.update(params)
plt.figure(1)
plt.clf()

plt.subplot(1,2,1)


markerSize = 6
# Spirals
yMinPos=25
yMaxPos=82

plt.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
             xerr=(Vcerr[0,yMinPos:yMaxPos],Vcerr[1,yMinPos:yMaxPos]),
             yerr=(Mberr[0,yMinPos:yMaxPos],Mberr[1,yMinPos:yMaxPos]),
                   marker='o',color='r',ecolor='k',linestyle='none',
                   markersize=markerSize)

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
lowPt = 10 ** ( 4.0 * math.log10(0.0001) + 1.65)
highPt = 10 ** ( 4.0 * math.log10(10**5) + 1.65)
plt.plot([0.0001,10**5],[lowPt,highPt],
         color='k',linestyle='dashed',linewidth=2)

plt.ylabel("$M_b$ (M$_{\odot}$)",family='serif')
plt.xlabel("$V_c$ (km s$^{-1}$)",family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1,10**3])
plt.ylim([50,10**13])
plt.legend(loc='upper left')
#plt.savefig('../figures/leop.mcGaughPlot.noLeoP.eps',dpi=600)

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


# Leo P
plt.errorbar(6,mtot,
             #xerr=([15,0]),yerr=([10**6,10**6]),
             xerr=([1],[1]),yerr=([mtotErrM],[mtotErrP]),
             xlolims=False,
             marker='s',color='c',ecolor='k',elinewidth=2,capsize=3,
             markersize=markerSize)

plt.annotate('(a)',
            (0.9,0.05)
            ,color='k'
            ,textcoords='axes fraction',
            xycoords='axes fraction')

#plt.savefig('../figures/leop.mcGaughPlot.dwarfs.eps',dpi=600)
#plt.show()




################################################################################
# Include only spirals and gas rich data
################################################################################

#plt.rcParams.update(params)
#plt.figure(1)
#plt.clf()

plt.subplot(1,2,2)

markerSize = 7
# Spirals
yMinPos=25
yMaxPos=82

plt.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
             xerr=(Vcerr[0,yMinPos:yMaxPos],Vcerr[1,yMinPos:yMaxPos]),
             yerr=(Mberr[0,yMinPos:yMaxPos],Mberr[1,yMinPos:yMaxPos]),
                   marker='o',color='r',ecolor='k',linestyle='none',
                   markersize=markerSize)

# Gas-rich galaxies
plt.errorbar(gasVc,gasMb,
             xerr=gasVcerr,
             yerr=gasMberr,
             marker='^',color='g',ecolor='k',linestyle='none',
             markersize=markerSize)

plt.ylabel("$M_b$ ($M_{\odot}$)",family='serif')
plt.xlabel("$v_c$ (km s$^{-1}$)",family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlim([3,10**2.8])
plt.ylim([10**5.6,10**12])
plt.legend(loc='upper left')
#plt.savefig('../figures/leop.mcGaughPlot.noLeoP.eps',dpi=600)

# Leo P
plt.errorbar(5.8,mtot,
             #xerr=([15,0]),yerr=([10**6,10**6]),
             xerr=([1],[1]),yerr=([mtotErrM],[mtotErrP]),
             xlolims=False,
             marker='s',color='c',ecolor='k',elinewidth=2,capsize=3,
             markersize=markerSize)

# Plot line for fitting spirals: function = log(M_b) = x log(V_c) + A
lowPt = 10 ** ( 4.0 * math.log10(0.0001) + 1.65)
highPt = 10 ** ( 4.0 * math.log10(10**5) + 1.65)
plt.plot([0.0001,10**5],[lowPt,highPt],
         color='k',linestyle='dashed',linewidth=2)
plt.annotate('(b)',
            (0.9,0.05)
            ,color='k'
            ,textcoords='axes fraction',
            xycoords='axes fraction')

plt.savefig('../figures/leop.mcGaughPlot.eps',dpi=600)
plt.show()













'''
plt.rcParams.update(params)
plt.figure(1)
plt.clf()

markerSize = 5
# Spirals
yMinPos=0
yMaxPos=81
plt.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
             xerr=(Vcerr[0,yMinPos:yMaxPos],Vcerr[1,yMinPos:yMaxPos]),
             yerr=(Mberr[0,yMinPos:yMaxPos],Mberr[1,yMinPos:yMaxPos]),
             marker='s',color='g',ecolor='k',linestyle='none',
             markersize=markerSize)

# Dwarfs
#yMinPos=82
#yMaxPos=114
#plt.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
#             xerr=(xerr[0,yMinPos:yMaxPos],xerr[1,yMinPos:yMaxPos]),
#yerr=(yerr[0,yMinPos:yMaxPos],yerr[1,yMinPos:yMaxPos]),
#             marker='o',color='r',ecolor='k',linestyle='none',
#             markersize=markerSize)

# Plot line for fitting spirals: function = log(M_b) = x log(V_c) + A
lowPt = 10 ** ( 4.0 * math.log10(0.0001) + 1.65)
highPt = 10 ** ( 4.0 * math.log10(10**5) + 1.65)
plt.plot([0.0001,10**5],[lowPt,highPt],
         color='b',linestyle='dashed',linewidth=2)

plt.ylabel("M$_b$ (M$_{\odot}$)",family='serif')
plt.xlabel("V$_c$ (km s$^{-1}$)",family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1,10**3])
plt.ylim([50,10**13])
plt.legend(loc='upper left')
#plt.savefig('../figures/leop.mcGaughPlot.noLeoP.nodSph.eps',dpi=600)

# Plot Leo P
plt.errorbar(15,1.82*10**6,
             #xerr=([15,0]),yerr=([10**6,10**6]),
             xerr=5,yerr=0.10**5,
             xlolims=True,
             marker='*',color='c',ecolor='k',elinewidth='2',
             markersize=markerSize*3.5)

plt.savefig('../figures/leop.mcGaughPlot.nodSph.eps',dpi=600)
plt.show()

############################################################################################
# Plot all four together
############################################################################################





plt.close('all')
plt.rcParams.update(params)
plt.figure(1)
plt.clf()



f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, sharex='col', sharey='row')

##########
# Plot 1 #
##########
# Spirals
ax1.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
             xerr=(xerr[0,yMinPos:yMaxPos],xerr[1,yMinPos:yMaxPos]),
             yerr=(yerr[0,yMinPos:yMaxPos],yerr[1,yMinPos:yMaxPos]),
             marker='s',color='g',ecolor='k',linestyle='none',
             markersize=markerSize)

# Dwarfs
yMinPos=82
yMaxPos=114
ax1.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
             xerr=(xerr[0,yMinPos:yMaxPos],xerr[1,yMinPos:yMaxPos]),
yerr=(yerr[0,yMinPos:yMaxPos],yerr[1,yMinPos:yMaxPos]),
             marker='o',color='r',ecolor='k',linestyle='none',
             markersize=markerSize)

# Plot line for fitting spirals: function = log(M_b) = x log(V_c) + A
lowPt = 10 ** ( 4.0 * math.log10(0.0001) + 1.65)
highPt = 10 ** ( 4.0 * math.log10(10**5) + 1.65)
ax1.plot([0.0001,10**5],[lowPt,highPt],
         color='b',linestyle='dashed',linewidth=2)

# Plot Leo P
ax1.errorbar(15,1.82*10**6,
             #xerr=([15,0]),yerr=([10**6,10**6]),
             xerr=5,yerr=0.10**5,
             xlolims=True,
             marker='*',color='c',ecolor='k',elinewidth='2',
             markersize=markerSize*3.5)

##########
# Plot 2 #
##########

# Spirals
ax2.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
             xerr=(xerr[0,yMinPos:yMaxPos],xerr[1,yMinPos:yMaxPos]),
             yerr=(yerr[0,yMinPos:yMaxPos],yerr[1,yMinPos:yMaxPos]),
             marker='s',color='g',ecolor='k',linestyle='none',
             markersize=markerSize)

# Dwarfs
yMinPos=82
yMaxPos=114
ax2.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
             xerr=(xerr[0,yMinPos:yMaxPos],xerr[1,yMinPos:yMaxPos]),
yerr=(yerr[0,yMinPos:yMaxPos],yerr[1,yMinPos:yMaxPos]),
             marker='o',color='r',ecolor='k',linestyle='none',
             markersize=markerSize)

# Plot line for fitting spirals: function = log(M_b) = x log(V_c) + A
lowPt = 10 ** ( 4.0 * math.log10(0.0001) + 1.65)
highPt = 10 ** ( 4.0 * math.log10(10**5) + 1.65)
ax2.plot([0.0001,10**5],[lowPt,highPt],
         color='b',linestyle='dashed',linewidth=2)


##########
# Plot 3 #
##########
# Spirals
ax3.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
             xerr=(xerr[0,yMinPos:yMaxPos],xerr[1,yMinPos:yMaxPos]),
             yerr=(yerr[0,yMinPos:yMaxPos],yerr[1,yMinPos:yMaxPos]),
             marker='s',color='g',ecolor='k',linestyle='none',
             markersize=markerSize)

# Plot line for fitting spirals: function = log(M_b) = x log(V_c) + A
lowPt = 10 ** ( 4.0 * math.log10(0.0001) + 1.65)
highPt = 10 ** ( 4.0 * math.log10(10**5) + 1.65)
ax3.plot([0.0001,10**5],[lowPt,highPt],
         color='b',linestyle='dashed',linewidth=2)

# Plot Leo P
ax3.errorbar(15,1.82*10**6,
             #xerr=([15,0]),yerr=([10**6,10**6]),
             xerr=5,yerr=0.10**5,
             xlolims=True,
             marker='*',color='c',ecolor='k',elinewidth='2',
             markersize=markerSize*3.5)

##########
# Plot 4 #
##########
# Spirals
ax4.errorbar(data[0,yMinPos:yMaxPos],data[1,yMinPos:yMaxPos],
             xerr=(xerr[0,yMinPos:yMaxPos],xerr[1,yMinPos:yMaxPos]),
             yerr=(yerr[0,yMinPos:yMaxPos],yerr[1,yMinPos:yMaxPos]),
             marker='s',color='g',ecolor='k',linestyle='none',
             markersize=markerSize)

# Plot line for fitting spirals: function = log(M_b) = x log(V_c) + A
lowPt = 10 ** ( 4.0 * math.log10(0.0001) + 1.65)
highPt = 10 ** ( 4.0 * math.log10(10**5) + 1.65)
ax4.plot([0.0001,10**5],[lowPt,highPt],
         color='b',linestyle='dashed',linewidth=2)




# Fine-tune figure; make subplots close to each other and hide x ticks for
# all but bottom plot.
f.subplots_adjust(hspace=0,vspace=0)

#plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
#plt.setp([a.get_yticklabels() for a in f.axes[:-1]], visible=False)

plt.ylabel("M$_b$ (M$_{\odot}$)",family='serif')
plt.xlabel("V$_c$ (km s$^{-1}$)",family='serif')
plt.xscale('log')
plt.yscale('log')
plt.xlim([1,10**3])
plt.ylim([50,10**13])
plt.legend(loc='upper left')

plt.savefig('../figures/leop.mcGaughPlot.allFour.eps',dpi=600)
plt.show()
'''
