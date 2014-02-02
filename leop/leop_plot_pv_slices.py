#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid
import pyfits
import matplotlib.pyplot as plt
import pywcsgrid2 as wcs

#fitsDir = '/home/elijah/research/leop/data/obsCombine/images/sliceFiles/'
#psDir =  '/home/elijah/research/leop/obsCombine/figures/'

fitsDir = '/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/slices/'
psDir = '/d/bip3/ezbc/leop/figures/'

fig = plt.figure(figsize=(8,6))
#fig = plt.figure()
fontScale = 10

plt.rcdefaults()

params = {'backend': 'eps',
          'axes.labelsize': fontScale,
          'text.fontsize': fontScale,
          'legend.fontsize': fontScale*3/4,
          'xtick.labelsize': fontScale,
          'ytick.labelsize': fontScale,
          #'text.usetex': True
          }

plt.rcParams.update(params)
plt.figure(1)
plt.clf()

# Open a fits file to use the header for axes info
f = pyfits.open(fitsDir + 'leop.16arcsec.majorSlice.01.fits')
hmajor = f[0].header

f = pyfits.open(fitsDir + 'leop.16arcsec.minorSlice.01.fits')
hminor = f[0].header

gridMajor = ImageGrid(fig, (2,1,1),
                      nrows_ncols=(1,5),
                      ngrids=5,
                      direction='column',
                      cbar_mode="single",
                      cbar_location='right',
                      cbar_pad="3%",
                      cbar_size='12%',
                      axes_pad=0.0,
                      axes_class=(wcs.Axes, dict(header=hmajor)),
                      aspect=False,
                      share_all=True,
                      label_mode='L')

gridMinor = ImageGrid(fig, (2,1,2),
                      nrows_ncols=(1,7),
                      ngrids=7,
                      cbar_mode="single",
                      cbar_location='right',
                      cbar_pad="6%",
                      cbar_size='20%',
                      axes_pad=0.0,
                      axes_class=(wcs.Axes, dict(header=hminor)),
                      aspect=False,
                      label_mode='L',
                      share_all=True)

#gridMajor.set_aspect(1)
#gridMinor.set_aspect(25)

fileNums = np.array(['01','02','03','04','05',
                     '06','07','08','09','10',
                     '11','12','13','14','15'])

sliceNum = np.arange(0,20,1)

vmin = 0.0
vmax = 4
velRange = (6,25)

annoColor='w'

for i in xrange(5):
    f = pyfits.open(fitsDir +
                    'leop.16arcsec.majorSlice.' +
                    fileNums[i] + '.fits')

    h, d = f[0].header, f[0].data[:,:]

    ax = gridMajor[i]

    im = ax.imshow(d * 1000, # mJy/bm
                   origin='lower', # put x-axis below plot
                   #cmap=plt.cm.rainbow,
                   vmin=vmin,vmax=vmax,
                   aspect='auto'
                   )
    #ax.set_aspect(5)

    if i == 0:
        ax.annotate('(a)',
                    (0.03,0.05)
                    ,color=annoColor
                    ,size='x-large',
                    textcoords='axes fraction',
                    xycoords='axes fraction')
    ax.set_xlabel('Offset (\'\')',
                  size = 'small',
                  family='serif')
    ax.set_ylabel('Velocity \n(km s$^{-1}$)',
                  size = 'small',
                  family='serif')

    ax.set_ticklabel2_type('manual',
                           locs=[240e3,250e3,260e3,270e3,280e3],
                           labels=['240','250','260','270','280'])

    ax.set_ticklabel1_type('manual',
                           locs=[-48,-16,16,48],
                           labels=['-48','-16','16','48'])

    ax.set_ylim(velRange[0],velRange[1])
    ax.set_xlim(10,103)

    ax.annotate('(' + str(sliceNum[i]) + ')',
                (0.05,0.9)
                ,color=annoColor
                ,textcoords='axes fraction',
                xycoords='axes fraction')

    # Set tick and lines to white
    ax.axis[:].major_ticks.set_color("w")

    if i==0:
        ax.axis['right'].line.set_color("w")
    elif i==4:
        ax.axis['left'].line.set_color('w')
    else:
        ax.axis['right'].line.set_color("w")
        ax.axis['left'].line.set_color("w")

for i in xrange(1,8):
    f = pyfits.open(fitsDir + 
                    'leop.16arcsec.minorSlice.' + 
                    fileNums[i] + '.fits')

    h, d = f[0].header, f[0].data[:,:]

    ax = gridMinor[i-1]

    im = ax.imshow(d * 1000, # mJy/bm
                   origin='lower', # put x-axis below plot
                   #cmap=plt.cm.Blues,
                   vmin=vmin,vmax=vmax,
                   aspect='auto'
                   )

    if i == 1:
        ax.annotate('(b)',
                    (0.03,0.05)
                    ,color=annoColor
                    ,size='x-large',
                    textcoords='axes fraction',
                    xycoords='axes fraction')

    ax.set_ylim(velRange[0],velRange[1])

    ax.set_xlabel('Offset (\'\')',
                  size = 'small',
                  family='serif')
    ax.set_ylabel('Velocity \n(km s$^{-1}$)',
                  size = 'small',
                  family='serif')

    ax.set_ticklabel2_type('manual',
                           locs=[240e3,250e3,260e3,270e3,280e3],
                           labels=['240','250','260','270','280'])

    ax.set_ticklabel1_type('manual',
                           locs=[-32,0,32],
                           labels=['-32','0','32'])

    ax.annotate('(' + str(sliceNum[i-1]) + ')',
            (0.05,0.9)
            ,color=annoColor
            ,textcoords='axes fraction',
            xycoords='axes fraction')

    # Set tick and lines to white
    ax.axis[:].major_ticks.set_color("w")

    if i==0:
        ax.axis['right'].line.set_color("w")
    elif i==7:    
        ax.axis['left'].line.set_color('w')
    else:
        ax.axis['right'].line.set_color("w")
        ax.axis['left'].line.set_color("w")


    #ax.set_aspect(aspect=1/2.,
    #              adjustable='box-forced',
    #              anchor='C')


# Correspond color bar to data
cb1 = gridMajor[0].cax.colorbar(im)
# Write label to colorbar
cb1.set_label_text('mJy Bm$^{-1}$',
                   size='small',
                   family='serif')

cb2 = gridMinor[0].cax.colorbar(im)
# Write label to colorbar
cb2.set_label_text('mJy Bm$^{-1}$',
                   size='small',
                   family='serif')

plt.savefig('/d/bip3/ezbc/leop/figures/leop.16slices.color.eps',
            dpi=400)#,bbox_inches='tight')
plt.savefig('/usr/users/ezbc/Dropbox/research/manuscripts/leop/paper5/' + \
            'post_circulation/finalImages/' + \
            '/leop.16slices.color.eps',
            dpi=400)#,bbox_inches='tight')
plt.show()


