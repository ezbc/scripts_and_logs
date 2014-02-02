# Script for plotting spectra HI and CO spectra in regions with
# integrated CO greater than 3sigma of the Taurus cloud


import numpy as np # make sure vector and array arithmetic options are loaded
import matplotlib.pyplot as plt
import pyfits as pf


# read files into python
plt.clf() # clear the plot (figure
colDens1 = pf.getdata('taurus.galfa.colDens.fits')
colDens2 = pf.getdata('taurus.galfa.colDens.vel2.fits')
colDens3 = pf.getdata('taurus.galfa.colDens.vel3.fits')

colDensList = [colDens1,colDens2,colDens3]
velocityStrings = ['-100 km s$^{-1}$ to 100 km s$^{-1}$',
                   '-35 km s$^{-1}$ to 15 km s$^{-1}$',
                   '-7 km s$^{-1}$ to 12 km s$^{-1}$']
# Log Plot
plt.clf() # clear the plot
fig = plt.figure()
ax = fig.add_subplot(111)    # The big subplot
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
axList = [ax1,ax2,ax3]

index = 0
for data in colDensList:
    whereAreNaNs = np.isnan(data);
    data[whereAreNaNs] = -1;
    hist, bins = np.histogram(data[0,:,:], bins=[i for i in xrange(1,70)])
    #width = 0.7*(bins[1]-bins[0])
    #center = (bins[:-1]+bins[1:])/2
    width = 1
    axList[index].bar(bins[:-1], hist, align = 'center', width=width, log=True)
    axList[index].set_xlim([0,70])
    #axList[index].set_ylim([0,2.5e5])
    axList[index].text(1,0.86,velocityStrings[index],
                       transform=axList[index].transAxes,
                       horizontalalignment='right')
    index += 1

# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

ax.set_title("Histogram of N(${HI})$ for 10 $\deg^2$ region of Taurus")
ax.set_ylabel("Frequency")
ax.set_xlabel("N(${HI}) (10^{20}$ cm$^{-2})$")
plt.savefig('colDensHistogram.logScale.ps',dpi=500)

# Linear Plot
plt.clf() # clear the plot
fig = plt.figure()
ax = fig.add_subplot(111)    # The big subplot
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
axList = [ax1,ax2,ax3]

index = 0
for data in colDensList:
    whereAreNaNs = np.isnan(data);
    data[whereAreNaNs] = -1;
    hist, bins = np.histogram(data[0,:,:], bins=[i for i in xrange(1,70)])
    #width = 0.7*(bins[1]-bins[0])
    #center = (bins[:-1]+bins[1:])/2
    width = 1
    axList[index].bar(bins[:-1], hist, align = 'center', width=width)
    axList[index].set_xlim([0,70])
    axList[index].set_ylim([0,2.5e5])
    axList[index].text(1,0.86,velocityStrings[index],
                       transform=axList[index].transAxes,
                       horizontalalignment='right')
    index += 1

# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

ax.set_title("Histogram of N(${HI})$ for 10 deg$^2$ region of Taurus")
ax.set_ylabel("Frequency",position=[0.5,0],transform=ax.transAxes)
ax.set_xlabel("N(${HI}) (10^{20}$ cm$^{-2})$")
ax2.set_ylabel("Frequency")
plt.savefig('colDensHistogram.linearScale.ps',dpi=500)


