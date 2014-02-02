#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

from kapteyn import maputils
from kapteyn import wcs as kwcs
import numpy as np
import math
from mpl_toolkits.axes_grid1 import ImageGrid
import pyfits as pf
import matplotlib.pyplot as plt
import pywcsgrid2 as wcs
import pywcs

lappy = True
cosmos = False

if lappy:
    fitsDir = '/home/elijah/research/leop/data/vla.gmrt/fitsImages/'
    figureDir =  '/home/elijah/research/leop/figures/'
elif cosmos:
    fitsDir = '/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/'
    figureDir = '/d/bip3/ezbc/leop/figures/'

fig = plt.figure(figsize=(7,12))

imageName = 'leop.fluxRescale.16arcsec'

f0 = pf.open(fitsDir + imageName + '.mom0.blk.fits')
hmom0, dmom0 = f0[0].header, f0[0].data[0,:,:]

f1 = pf.open(fitsDir + imageName + '.mom1.2pt5sig.blk.fits')
hmom1, dmom1 = f1[0].header, f1[0].data[0,:,:]

f2 = pf.open(fitsDir + imageName + '.mom2.2pt5sig.blk.fits')
hmom2, dmom2 = f2[0].header, f2[0].data[0,:,:]

grid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(3,2),
                 ngrids=6,
                 direction='row',
                 cbar_mode="all",
                 cbar_location='right',
                 cbar_pad='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes, 
                             dict(header=hmom0)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

mom0a = grid[0]
mom0b = grid[1]
mom1a = grid[2]
mom1b = grid[3]
mom2a = grid[4]
mom2b = grid[5]

vminMom0=30
vmaxMom0=99
vminMom1=250
vmaxMom1=277
vminMom2=4
vmaxMom2=11

mom0aIm = mom0a.imshow(dmom0 * 1000, # mJy/bm
                       origin='lower', # put x-axis below plot
                       #cmap=plt.cm.gray_r,
                       vmin=vminMom0,vmax=vmaxMom0)

mom0aIm = mom0b.imshow(dmom0 * 1000, # mJy/bm
                       origin='lower', # put x-axis below plot
                       #cmap=plt.cm.gray_r,
                       vmin=vminMom0,vmax=vmaxMom0)

mom1aIm = mom1a.imshow(dmom1, # kms
                       origin='lower', # put x-axis below plot
                       #cmap=plt.cm.gray_r,
                       vmin=vminMom1,vmax=vmaxMom1)

mom1bIm = mom1b.imshow(dmom1, # kms
                       origin='lower', # put x-axis below plot
                       vmin=vminMom1,vmax=vmaxMom1)

mom2aIm = mom2a.imshow(dmom2, # kms
                       origin='lower', # put x-axis below plot
                       vmin=vminMom2,vmax=vmaxMom2)

mom2bIm = mom2b.imshow(dmom2, # kms
                       origin='lower', # put x-axis below plot
                       vmin=vminMom2,vmax=vmaxMom2)

momArray = np.array([mom0a,mom0b,mom1a,mom1b,mom2a,mom2b])

for i, ax in enumerate(momArray):

    xoffset = 55
    yoffset = 80
    ax.set_xlim(285 - xoffset, 285 + xoffset)
    ax.set_ylim(300 - yoffset, 300 + yoffset)
   
    ax.set_xlabel('Right Ascension (J2000)', 
                  size = 'x-small', 
                  family='serif')
    ax.set_ylabel('Declination (J2000)', 
                  size = 'x-small',
                  family='serif')

    # colorbar axes
    if i%2 == 0:
        from mpl_toolkits.axes_grid.inset_locator import inset_axes

        axins = inset_axes(grid[i],
                           width="95%",
                           height='4%',
                           loc=2,
                           )

        axins.axis[:].toggle(all=False)
        axins.axis["top"].toggle(all=True)
        axins.axis["top"].label.set_size(10)

        # colorbar
        cb = plt.colorbar(grid[i].images[0], cax = axins,
                          orientation='horizontal')

        #cb.makes_axes(orientation='horizontal')
        axins.axis[:].toggle(all=False)
        axins.axis["bottom"].toggle(all=True)

        if i==0:
            cb.set_label('mJy Bm$^{-1}$',
                   size='small',
                   family='serif')
            cb.set_ticks([40, 60, 80])
        elif i==2:
            cb.set_label('km s$^{-1}$',
                   size='small',
                   family='serif')
            cb.set_ticks([255,265,275])
        else:
            cb.set_label('kms s$^{-1}$',
                   size='small',
                   family='serif')
            cb.set_ticks([4,6,8,10])

# PA of major slice: 154.94 deg
# dimensions of major slice: 79.4862' x 160'
# seperation: 2.777'
# averaging width = 11 pix

# PA of minor slice: deg
# length of minor slice: '

import math
import numpy as np

beamsize = 16./3600.
startRA = 155.438
startDec = 18.0883
PA = 65.

# Minor axis cuts
ras = np.zeros(21)
decs = np.zeros(21)

for i in xrange(0,11):
    ras[i+10] = startRA + beamsize * i * math.sin(math.radians(PA))
    decs[i+10] = startDec + beamsize * i * math.cos(math.radians(PA))

for i in xrange(0,11):
    ras[10 - i] = startRA - beamsize * (i) * math.sin(math.radians(PA))
    decs[10 - i] = startDec - beamsize * (i) * math.cos(math.radians(PA))

majorCoordArray = np.array([ras,decs]).T

PA = PA + 90.

for i in xrange(0,11):
    ras[i+10] = startRA + beamsize * i * math.sin(math.radians(PA))
    decs[i+10] = startDec + beamsize * i * math.cos(math.radians(PA))

for i in xrange(0,11):
    ras[10 - i] = startRA - beamsize * (i) * math.sin(math.radians(PA))
    decs[10 - i] = startDec - beamsize * (i) * math.cos(math.radians(PA))


minorCoordArray = np.array([ras,decs]).T

'''


In [141]: majorCoordArray
Out[141]: 
array([[ 155.39771965,   18.06951697],
       [ 155.40174769,   18.07139527],
       [ 155.40577572,   18.07327357],
       [ 155.40980376,   18.07515188],
       [ 155.41383179,   18.07703018],
       [ 155.41785983,   18.07890848],
       [ 155.42188786,   18.08078679],
       [ 155.4259159 ,   18.08266509],

       [ 155.42994393,   18.08454339],
       [ 155.43397197,   18.0864217 ],
       [ 155.438     ,   18.0883    ],
       [ 155.44202803,   18.0901783 ],
       [ 155.44605607,   18.09205661],
       
       [ 155.4500841 ,   18.09393491],
       [ 155.45411214,   18.09581321],
       [ 155.45814017,   18.09769152],
       [ 155.46216821,   18.09956982],
       [ 155.46619624,   18.10144812],
       [ 155.47022428,   18.10332643],
       [ 155.47425231,   18.10520473],
       [ 155.47828035,   18.10708303]])

In [142]: minorCoordArray
Out[142]: 
array([[ 155.41921697,   18.12858035],
       [ 155.42109527,   18.12455231],
       [ 155.42297357,   18.12052428],
       [ 155.42485188,   18.11649624],
       [ 155.42673018,   18.11246821],
       [ 155.42860848,   18.10844017],

       [ 155.43048679,   18.10441214],
       [ 155.43236509,   18.1003841 ],
       [ 155.43424339,   18.09635607],
       [ 155.4361217 ,   18.09232803],
       [ 155.438     ,   18.0883    ],
       [ 155.4398783 ,   18.08427197],
       [ 155.44175661,   18.08024393],
       [ 155.44363491,   18.0762159 ],
       [ 155.44551321,   18.07218786],

       [ 155.44739152,   18.06815983],
       [ 155.44926982,   18.06413179],
       [ 155.45114812,   18.06010376],
       [ 155.45302643,   18.05607572],
       [ 155.45490473,   18.05204769],
       [ 155.45678303,   18.04801965]])


'''
# PA of major slice: 154.94 deg
# dimensions of major slice: 79.4862' x 160'
# seperation: 2.777'
# averaging width = 11 pix


# PA of major slice: 64.97 deg
# dimensions of major slice: 92.51 x 48'
# seperation: 1.536'
# averaging width = 11 pix

# Arrow is defined as origin and dx,dy
# A change in coordinates is required:

# Major slice calculation:
majorPosOffset = 2.77 / 60. / 2 # 1/2 of length in degrees
minorPosOffset = 1.536 / 60. / 2 # 1/2 of length in degrees
majorPA = 155 - 180
minorPA = 155 - 90

majorCoordOriginsArray = np.zeros(majorCoordArray.shape)
majorCoordTipsArray = np.zeros(majorCoordArray.shape)
minorCoordOriginsArray = np.zeros(minorCoordArray.shape)
minorCoordTipsArray = np.zeros(minorCoordArray.shape)

for i in xrange(0,majorCoordArray.shape[0]):
    majorCoordOriginsArray[i,0] = majorCoordArray[i,0] - \
        majorPosOffset * math.sin(math.radians(majorPA))

    majorCoordTipsArray[i,0] = majorCoordArray[i,0] + \
        majorPosOffset * math.sin(math.radians(majorPA))

for i in xrange(0,majorCoordArray.shape[0]):
    majorCoordOriginsArray[i,1] = majorCoordArray[i,1] - \
        majorPosOffset * math.cos(math.radians(majorPA))

    majorCoordTipsArray[i,1] = majorCoordArray[i,1] + \
        majorPosOffset * math.cos(math.radians(majorPA))

for i in xrange(0,minorCoordArray.shape[0]):
    minorCoordOriginsArray[i,0] = minorCoordArray[i,0] - \
        minorPosOffset * math.sin(math.radians(minorPA))

    minorCoordTipsArray[i,0] = minorCoordArray[i,0] + \
        minorPosOffset * math.sin(math.radians(minorPA))

for i in xrange(0,minorCoordArray.shape[0]):
    minorCoordOriginsArray[i,1] = minorCoordArray[i,1] - \
        minorPosOffset * math.cos(math.radians(minorPA))

    minorCoordTipsArray[i,1] = minorCoordArray[i,1] + \
        minorPosOffset * math.cos(math.radians(minorPA))

# Extract pixel coordinates from wcs coordinates
# Define header for calculations
proj1 = kwcs.Projection(hmom0)
# extract coordinate information of extra axes
extraCoords = proj1.toworld((1.,1.,1.))

def arrowHead(image,x,y,dx,dy,head='origin'):
    theta = math.atan2(-dy, -dx)
    # make body
    image.arrow(x, y, dx, dy)
    # Make head
    barb = 10
    phi = math.radians(10)
    if head=='origin':
        rho = theta + phi
    else:
        rho = math.atan2(dy,dx) + phi

    for i in xrange(0,2):
        dx2 = barb * math.cos(rho)
        dy2 = barb * math.sin(rho)
        if head=='origin':
            image.arrow(x + dx, y + dy, # start position
                        dx2   , dy2    # delta x and y
                        )
        else:
            image.arrow(x,y,dx2,dy2)
        if head=='origin':
            rho = theta - phi
        else:
            rho = math.atan2(dy,dx) - phi

# define number of major arrows
center = 10
majorArrowNum = 5
majorLowLim = center - majorArrowNum / 2
majorHighLim = center + majorArrowNum / 2 + 1
count = 0

for i in xrange(majorLowLim, majorHighLim):
    wcsCoords = np.array([[majorCoordOriginsArray[i,0],
                           majorCoordOriginsArray[i,1],
                           extraCoords[2],
                           ],
                          [majorCoordTipsArray[i,0],
                           majorCoordTipsArray[i,1],
                           extraCoords[2],
                           ]])

    pixCoords = np.zeros(wcsCoords.shape)
    pixCoords[0] = proj1.topixel(wcsCoords[0])
    pixCoords[1] = proj1.topixel(wcsCoords[1])

    arrowLengths = np.array([pixCoords[1,0] - pixCoords[0,0],
                         pixCoords[1,1] - pixCoords[0,1]])

    for j, ax in enumerate([mom0b,mom1b,mom2b]):
        arrowHead(ax,
                  pixCoords[0,0], pixCoords[0,1],
                  arrowLengths[0], arrowLengths[1],
                  head='')

        ax.annotate(str(count),xy=(pixCoords[0,0]+arrowLengths[0]+5,
                                   pixCoords[0,1]+arrowLengths[1]+5),
                    size='x-small')
    count += 1


# Minor arrows
center = 10
minorArrowNum = 7
minorLowLim = center - minorArrowNum / 2
minorHighLim = center + minorArrowNum / 2 + 1
count = 0

for i in xrange(minorLowLim, minorHighLim):
    wcsCoords = np.array([[minorCoordOriginsArray[i,0],
                           minorCoordOriginsArray[i,1],
                           extraCoords[2],
                           ],
                          [minorCoordTipsArray[i,0],
                           minorCoordTipsArray[i,1],
                           extraCoords[2],
                           ]])

    pixCoords = np.zeros(wcsCoords.shape)
    pixCoords[0] = proj1.topixel(wcsCoords[0])
    pixCoords[1] = proj1.topixel(wcsCoords[1])

    arrowLengths = np.array([pixCoords[1,0] - pixCoords[0,0],
                         pixCoords[1,1] - pixCoords[0,1]])


    for j, ax in enumerate([mom0b,mom1b,mom2b]):
        arrowHead(ax,
                  pixCoords[0,0], pixCoords[0,1],
                  arrowLengths[0], arrowLengths[1])

        ax.annotate(str(count),xy=(pixCoords[0,0]+10, pixCoords[0,1]-10),
                    size='x-small')
    count += 1
#print wcs1.wcs_pix2sky([20.,20.,0,0],1)

labelArray=['a','b','c','d','e','f']

for i, ax in enumerate(momArray):

    # Beam size
    ax.add_beam_size(14.54, 14.54, 0, loc=4)

    # Figure label
    ax.annotate('(' + labelArray[i] + ')',
                   xy=(10 - xoffset,50+yoffset),
                   size='large')

plt.savefig(figureDir + 'leop.16maps.mom012.eps',
            dpi=600)
plt.savefig('/home/elijah/Dropbox/research/manuscripts/leop/paper5/' + \
            'post_circulation/finalImages/' + \
            'leop.mom012.slices.eps',
            dpi=400)

plt.show()

