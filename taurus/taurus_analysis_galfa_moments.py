# Script for plotting spectra HI and CO spectra in regions with
# integrated CO greater than 3sigma of the Taurus cloud


import numpy as np # make sure vector and array arithmetic options are loaded
import matplotlib.pyplot as plt
import pyfits as pf


# read files into python

mom0_1 = pf.getdata('taurus.galfa.mom0.fits')
mom0_2 = pf.getdata('taurus.galfa.mom0.vel2.fits')
mom0_3 = pf.getdata('taurus.galfa.mom0.vel3.fits')



plt.clf() # clear the plot (figure
plt.bar(center, hist, align = 'center', width = width)
plt.xlabel("Velocity (km/s)")
plt.ylabel("Integrated Tb (K km/s)")
plt.xlim([0,100])
plt.savefig('histogram.ps',dpi=500)
plt.show()
