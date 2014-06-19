#!/usr/bin/python

# Script for finding extent of HI associated with the Taurus molecular cloud
# Using methods from Andersson et al. 1991, ApJ, 366, 464

from pylab import *
import numpy as np # make sure vector and array arithmetic options are loaded
import matplotlib.pyplot as plt
import pyfits as pf

# Get the data
#----------------

# Read fits files
dirPath = '/d/bip3/ezbc/taurus/data/hi+co/'
hiData = pf.getdata(dirPath + 'taurus.galfa.cube.cfaRes.regrid.fits')
hiHdr = pf.getheader(dirPath + 'taurus.galfa.cube.cfaRes.regrid.fits')
coData = pf.getdata(dirPath + 'taurus.cfa.co_Int.fits')
coHdr = pf.getheader(dirPath + 'taurus.cfa.co_Int.fits')







# Begin plotting

golden_mean = (sqrt(5)-1.0)/2.0
fig_width  = 3.4 # column width in inches
fig_height = fig_width * golden_mean # Good enough for the greeks
fig_size   =  [fig_width,fig_height]

# display environment    
params = {'backend': 'eps',         # use eps output
          'axes.labelsize': 10,
          'font.family':'serif',    # font family
          'text.fontsize': 8,       # 8pt font
          'font.size':8,            
          'font.serif': 'Computer Modern Roman', # font type
          'legend.fontsize': 8,
          'legend.markersize': 2,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,   # will need TeX installtion, but will get LaTeX labels
          'figure.figsize': fig_size, # set above
          'figure.dpi': 600,          # resolution of figure in DPI (300 is ApJ for images, 600+ for line plots)
          'lines.markersize':3,
          'lines.linewidth': 1,
          'lines.dashes':(),
          }

#Set plot environment with the parameters defined above
rcParams.update(params)

# Create plotting space
figure()

ax1 = subplot(211)
ax2 = subplot(212)


# write out to EPS
savefig('plot_demo.eps', dpi=600)  # this takes the place of show() or draw()
