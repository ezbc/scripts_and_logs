#!/usr/bin/python

import numpy as np
import math
from mpl_toolkits.axes_grid1 import ImageGrid
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pyfits as pf
import matplotlib.pyplot as plt
import pywcsgrid2 as wcs
import pywcs

#fitsDir = '/home/elijah/research/leop/data/obsCombine/images/'
#psDir =  '/home/elijah/research/leop/figures/'

fitsDir = '/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/slices/'
figureDir = '/d/bip3/ezbc/leop/figures/'
finalImagesDir = '/usr/users/ezbc/Dropbox/research/manuscripts/' + \
    'leop/paper5/post_circulation/finalImages/'



fig = plt.figure(figsize=(7,4))

#16 beam 4"/pix
#32 beam 4"/pix

f0 = pf.open(fitsDir +  'leop.vla.gmrt.16arcsec.3beamMajorSlice.fits')

beamSize16 = 1.13 * 16**2 / 4.**2
n_hi16 = 6.057493205*10**5 / (16**2) *1/1000*1.8224*10**18*10**(-19)

h16, d16 = f0[0].header, f0[0].data[:,:]*1000

grid = ImageGrid(fig, (1,1,1),
                      nrows_ncols=(1,1),
                      ngrids=1,
                      direction='column',
                      cbar_mode="each",
                      cbar_location='top',
                      cbar_pad="3%",
                      cbar_size='6%',
                      axes_pad=0.0,
                      axes_class=(wcs.Axes, dict(header=h16)),
                      aspect=True,
                      label_mode='L')

im16 = grid[0].imshow(d16 * n_hi16, # mJy/bm
        origin='lower', # put x-axis below plot
        #cmap=plt.cm.gray_r,
        vmin=0.01,vmax=1.5)

ab = ['(a)','(b)']
im = [im16]

for i, ax in enumerate(grid):

    ax.set_xlabel('Offset (\'\')',
                  #size = 'small', 
                  family='serif')
    ax.set_ylabel('Velocity (km s$^{-1}$)',
                  #size = 'small',
                  family='serif')

    # Set tick and lines to white
    ax.axis[:].major_ticks.set_color("w")

    if i==1:
        ax.axis['left'].line.set_color("w")
    if i==0:
        ax.axis['right'].line.set_color("w")

    ax.set_ticklabel2_type('manual',
                           locs=[146,220,200],
                           labels=['a',220,'b'])

    grid[i].set_ylim(0,)
    #grid[i].set_xlim(

    grid[i].grid(color='w')

    #ax.annotate(ab[i],
    #            (10,10)
    #            ,color='w'
    #            ,textcoords='axes pixels')

    cb = grid[i].cax.colorbar(im[i])
    #cb = plt.colorbar(grid[i].images[0],orientation='horizontal')

    if i==0:
        cb.set_label_text('10$^{19}$ cm$^{-2}$',
                     size='small',
                     family='serif')
        grid[i].cax.set_xticks([0.2,0.6,1.0,1.4])

    if i==1:
        cb.set_label_text('10$^{19}$ cm$^{-2}$',
                     size='small',
                     family='serif')
        grid[i].cax.set_xticks([0.1,0.3,0.5])



# Save figure
filename='leop.16arcsec.majorSlice'
fig.savefig(figureDir + filename + '.eps',dpi=400,)#bbox_inches='tight',pad=0.1)
fig.savefig(finalImagesDir + filename + '.eps',dpi=400,)#bbox_inches='tight',

#plt.show()

