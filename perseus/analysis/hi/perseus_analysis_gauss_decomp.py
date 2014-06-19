#!/usr/bin/python

guesses = [50,0,10,20,-5,10,10,-10,10,30,-20,50,10,-10,10,30,-5,10,10,-20,10]

###############
# 16' cube
###############

reload(grid)

galfa4 = grid.SpectralGrid('perseus.galfa.cube.bin_4arcmin.sub.fits',
                        box=[100,100,180,180],
                        noiseScale=20.,
                        noiseRange=((-140,-90),(80,97)),
                        basesubtract=True)

guesses = [3,-30,15, 15,0,20, 3,0,20]

guesses = [50,5,15, 10,20,10, 4,-25,30]

galfa4.fit_profile(guesses,len(guesses)/3,166,158,coords='image')


galfa4.fit_profiles(
    growPos = (140,140),
    tileSaveFreq=10,
    threshold = 1,
    #filename='grid.4arcmin.noiseX1',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    coords='image',
    numberOfFits=4000)

grid.plot_spectrum(galfa4,166,158,coords='image')


grid.plot_fit(galfa4,166,158,coords='image')

grid.plot_fits(galfa4,(715,715,715,715),(414,415,416,417),coords='image')
grid.plot_fits(galfa4,(86,86,86,86),(83,84,85,86),coords='image')
grid.plot_fits(galfa4,(30,30,30,30),(29,30,31,32))
grid.plot_residualsImage(galfa4)
grid.plot_componentPV(galfa4,yslice=20,width=1)
grid.plot_ncompImage(galfa4)
grid.plot_velocityHistograms(galfa4)

grid.plot_NHI(galfa4,velocityrange=[-10,5])

reload(grid)

galfa4 = grid.SpectralGrid('perseus.galfa.cube.bin.16arcmin.fits',
                        box=[60,60,120,120],
                        noiseScale=50.,
                        noiseRange=((-110,-80),(80,110)),
                        basesubtract=True)

guesses = [  4.07575984, -26.97410452,  16.75406491,
 53.70105613,   4.56270824,   4.9002332 ]

galfa4.fit_profiles(
    growPos = (85,85),
    tileSaveFreq=10,
    threshold = 1,
    filename='grid.16arcsec.85_85.noiseX50',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    coords='image',
    numberOfFits=3600)

reload(grid)

galfa4 = grid.SpectralGrid('perseus.galfa.cube.bin.16arcmin.fits',
                        box=[60,60,120,120],
                        noiseScale=70.)

guesses = [  4.07575984, -26.97410452,  16.75406491,
 53.70105613,   4.56270824,   4.9002332 ]

galfa4.fit_profiles(
    growPos = (85,85),
    tileSaveFreq=10,
    threshold = 1,
    filename='grid.16arcsec.85_85.noiseX70',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    coords='image',
    numberOfFits=3600)


guesses = [  1.84033942, -26.40142723,   4.09068543,
 1.73964572, -42.0671819,    5.4270631 ,
  3.35117293, -12.63181254,  24.44756508,
 32.84281479,   0.5952481,    2.8300867 ,
44.84499815,   6.96796314,   3.34573227]

guesses = [  1.84033942, -26.40142723,   4.09068543,
 1.73964572, -42.0671819,    5.4270631 ,
  3.35117293, -12.63181254,  24.44756508,
 32.84281479,   0.5952481,    2.8300867 ,
44.84499815,   6.96796314,   3.34573227]

guesses = [50,0,10, 50,10,10, 5,-25,50, 5,-40,5, 10,-25,5]

guesses = [50,5,20, 5,-40,10, 5,-25,20]

galfa4.fit_profile(guesses,len(guesses)/3,25,25)


grid.plot_fit(galfa4,80,86,coords='image')
grid.plot_fits(galfa4,(84,85,85,85),(85,84,85,86),coords='image')
grid.plot_fits(galfa4,(85,85,85,85),(83,84,85,86),coords='image')
grid.plot_fits(galfa4,(80,80,80,80),(83,84,85,86),coords='image')

grid.plot_fits(galfa4,(5,5,5,5),(1,2,3,4),,coords='grid')

grid.plot_fit(galfa4,5,3)
grid.plot_fit(galfa4,85,83,coords='image')

grid.plot_fits(galfa4,(79,79,79,79),(74,75,76,77),coords='image')
grid.plot_fits(galfa4,(71,71,71,71),(74,75,76,77),coords='image')

grid.plot_residualsImage(galfa4)
grid.plot_componentPV(galfa4,yslice=25,width=3)
grid.plot_ncompImage(galfa4)
grid.plot_velocityHistograms(galfa4)

# for 75_75.5sigma.01alpha we see ncomp varying by 6 across (73,73:76)
grid.plot_fits(galfa4,(73,73,73,73),(73,74,75,76),coords='image')

# Try with initial guesses
galfa4 = grid.SpectralGrid('perseus.galfa.cube.bin.16arcmin.fits',box=[80,80,90,90])

# Save the grid
galfa4.write_grid('grid.16arcsec')

# Load the grid
# Allows changes to be made in SpectralGrid class without having to recreate the instance
galfa4 = grid.SpectralGrid('perseus.galfa.cube.bin.16arcmin.fits',gridfile='grid.16arcsec')



################################################################################
# Attempt decomposition on CfA 12CO J 2-->1 cube
###############################################################################

from mirexec import TaskRegrid as regrid

reload(grid)

cfa = grid.SpectralGrid('../cfa/perseus.cfa.cube.fits',
        noiseRange=((-0.0224,-0.01852915),(0.02438339,0.02958491)),
        basesubtract=False,
        box=[110,140,140,160])

cfa = grid.SpectralGrid('../cfa/perseus.cfa.cube.fits',
        noiseRange=((-0.0224,-0.01852915),(0.02438339,0.02958491)),
        basesubtract=False,
        box=[40,70,190,290])

guesses = [4,0.005,0.005]

cfa.fit_profile(guesses,len(guesses)/3,116,153,coords='image')

grid.plot_spectrum(cfa,116,153,coords='image')
grid.plot_spectrum(cfa,121,151,coords='image')


grid.plot_fit(cfa,116,153,coords='image')

cfa.fit_profiles(
    growPos = (121,151),
    tileSaveFreq=10,
    threshold = 1,
    filename='cfa.116_153',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=4000)

grid.plot_residualsImage(cfa)
grid.plot_componentPV(cfa,yslice=80,width=1)
grid.plot_ncompImage(cfa)
grid.plot_velocityHistograms(cfa)


grid.plot_fits(cfa,(116,116,116,116),(153,154,155,156),coords='image')


taurus = grid.SpectralGrid('taurus.cfa.cube.fits',
                        box=[60,60,90,90],
                        noiseScale=20.,
                        noiseRange=((-140,-90),(80,97)),
                        basesubtract=False)





################################################################################################################################################################################################################################################
# Performing decomposition on CfA map of Perseus
################################################################################
################################################################################
################################################################################


# redoing CfA fit 
reload(grid)

box=[0,36,180,180]
#box=[360,360,400,400]
cfa = grid.SpectralGrid('../cfa/perseus.cfa.cube.galfaBin.4arcmin.fits',
        noiseRange=((-21,-15),(21,30)),
        basesubtract=False,
        box=box,
        Tsys=400.)

# chose Tsys = 400 from Dame et al. (2001), ApJ, 547, 792

# Most of the data from the northern telescope were obtained with an extremely
# sensitive SIS heterodyne receiver which was installed in 1983 (Pan 1984). Its
# single-sideband noise temperature of  65 K yields total system temperatures
# referred to above the atmosphere of 400 - 800 K for the normal range of
# elevations observed, 30deg - 75deg

grid.plot_spectrum(cfa,121,84,coords='image')

guesses = [3.5,4,2]
cfa.fit_profile(guesses,len(guesses)/3,121,84,coords='image')

grid.plot_fit(cfa,121,84,coords='image')

guesses = [3.5,4,2]
cfa.fit_profiles(
    growPos = (121,84),
    tileSaveFreq=400,
    threshold = 5,
    filename='perseus.cfa.121_84.5sigma',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=1e6)

# failed on position 87,86 due to 0 initial components
grid.plot_spectrum(cfa,88,86,coords='image')

reload(grid)
box=[80,80,100,100]
cfa = grid.SpectralGrid('../cfa/perseus.cfa.cube.galfaBin.4arcmin.fits',
        noiseRange=((-21,-15),(21,30)),
        basesubtract=False,
        box=box,
        Tsys=400.)
guesses = [3.5,4,2]
cfa.fit_profiles(
    growPos = (90,90),
    tileSaveFreq=400,
    threshold = 3,
    #filename='perseus.cfa.121_84.5sigma',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.05,
    verbose=True,
    coords='image',
    numberOfFits=500)

# 88,86 has no emission, should not have been fit
grid.plot_fit(cfa,88,93,coords='image')

grid.plot_spectrum(cfa,94,91,coords='image')


grid.plot_fit(cfa,28,38,coords='image')



grid.plot_fit(cfa,89,82,coords='image')


grid.plot_fit(cfa,85,86,coords='image')
cfa.get_tile(85,86,coords='image').ncomps

grid.plot_fit(60,102,coords='image')
50, 80

'''
# redoing CfA fit 
reload(grid)

box=[0,36,180,180]
#box=[360,360,400,400]
cfa = grid.SpectralGrid('../cfa/perseus.cfa.cube.galfaBin.4arcmin.fits',
        noiseRange=((-21,-15),(21,30)),
        basesubtract=False,
        box=box,
        Tsys=400.)

guesses = [3.5,4,2]
cfa.fit_profiles(
    growPos = (121,84),
    tileSaveFreq=400,
    threshold = 5,
    filename='perseus.cfa.121_84.5sigma',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=1e6)

grid.plot_ncompImage(cfa)
grid.plot_velocityHistograms(cfa)

# redoing CfA fit 
reload(grid)

box=[0,60,10,80]
#box=[360,360,400,400]
cfa = grid.SpectralGrid('../cfa/perseus.cfa.cube.galfaBin.4arcmin.fits',
        noiseRange=((-21,-15),(21,30)),
        basesubtract=False,
        box=box,
        Tsys=400.)

guesses = [3.5,4,2]
cfa.fit_profiles(
    growPos = (2,65),
    tileSaveFreq=400,
    threshold = 5,
    filename='perseus.cfa.40_40.5sigma',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=1e6)
'''
# redoing CfA fit 
reload(grid)

box=[36,0,180,180]
cfa = grid.SpectralGrid('../cfa/perseus.cfa.cube.galfaBin.4arcmin.fits',
        noiseRange=((-21,-15),(21,30)),
        basesubtract=False,
        box=box,
        Tsys=400.)
grid.plot_spectrum(cfa,138,62,coords='image')

guesses = [10,8,2]
cfa.fit_profiles(
    growPos = (138,63),
    tileSaveFreq=20000,
    threshold = 5,
    filename='perseus.cfa.138_63.5sigma',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=1e6)

grid.plot_ncompImage(cfa)

################################################################################################################################################################
################################################################################
################################################################################
# Performed fits to GALFA HI cube using decomposed CfA CO cube to disclude
# velocities in fit
# script at:
# /d/bip3/ezbc/taurus/logFiles/perseus.gaussianDecomposition.fits1.py
################################################################################
################################################################################
################################################################################
################################################################################


grid.plot_spectrum(galfa0,40,40,coords='image')
guesses = [30,5,15]
grid.plot_spectrum(galfa0,100,100,coords='image')
guesses = [45,5,10, 4,-35,5]
grid.plot_spectrum(galfa0,140,155,coords='image')
guesses = [35,5,15]
grid.plot_spectrum(galfa0,10,120,coords='image')
guesses = [30,5,10]

grid.plot_ncompImage(galfa1)
grid.plot_NHI(galfa1,velocityrange=[-15,25],returnimage=True)
grid.plot_NHI(galfa1,velocityrange=[-15,25],returnimage=False)
grid.plot_velocityHistograms(galfa1)
grid.plot_fit(galfa1,90,90,coords='image')
grid.plot_spectrum(cfa,90,90,coords='image')

grid.plot_fit(galfa1,87,87,coords='image')
grid.plot_spectrum(cfa,87,87,coords='image')

grid.plot_fit(galfa1,127,62,coords='image')
grid.plot_spectrum(cfa,127,62,coords='image')

grid.plot_fit(cfa,40,40,coords='image')

reload(grid)
box=[0,36,180,180]
lee = grid.SpectralGrid('../galfa/perseus.lee12.cube.galfaBin.4arcmin.fits',
                        box=box,
                        noiseScale=10.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=False)

grid.display_image(lee)

grid.compare_grids(lee,galfa1)

box=[0,36,180,180]
reload(grid)
galfa0 = grid.SpectralGrid('../galfa/perseus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=10.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)
guesses = [48,5,10]
galfa0.fit_profiles(
        growPos = (138,62),
        tileSaveFreq=100,
        threshold = 1,
        filename='perseus.galfa0',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=1e6,
        COcube=cfa,
        COwidthScale=1.)

import pyfits as pf
from matplotlib import pyplot as plt
import numpy as np

# histogram of the differences between HI column densities calculated from
# integrating all emission in the GALFA cube betwen -5 and 15 km/s.
fitsfile='../galfa/perseus.galfa.cube.bin.4arcmin.fits'
f = pf.open(fitsfile)
h,galfacube = f[0].header, f[0].data
galfacubeNHI = galfacube[134:162,:,:].sum(axis=0) * \
        1.823 * h['cdelt3'] / 1000.* 1e-2

fitsfile='../galfa/perseus.lee12.galfaBin.4arcmin.fits'
f = pf.open(fitsfile)
h,leeNHI = f[0].header, f[0].data
diff = galfacubeNHI - leeNHI / 1e20
frac = galfacubeNHI / (leeNHI / 1e20)
plt.hist(frac[~np.isnan(diff)].ravel(),bins=50,normed=1)
plt.rc('text', usetex=True)
plt.xlabel(r'N(HI)$_{int}$ - N(HI)$_{Lee12}$ (1 $\times$ 10$^{20}$' + \
            'cm$^{-2}$)')
plt.ylabel('Fraction of pixels')
plt.title('Frequency of percent differences of N(HI) from \n' + \
          'GALFA DR1 cube and Lee et al. (2012) N(HI) map')
plt.savefig('/d/bip3/ezbc/perseus/figures/' + \
            'perseus.lee12.galfaInt.frac.histogram.png')
plt.show()

# histogram of the differences between HI column densities calculated from the
# decomposition and Lee et al. (2012).
galfaNHI = grid.plot_NHI(galfa0,velocityrange=[-5,15],returnimage=True)
fitsfile='../galfa/perseus.lee12.galfaBin.4arcmin.fits'
f = pf.open(fitsfile)
h,leeNHI = f[0].header, f[0].data
diff = galfaNHI - leeNHI / 1e20
plt.hist(diff[~np.isnan(diff)].ravel(),bins=100,normed=1)
plt.rc('text', usetex=True)
plt.xlabel(r'N(HI)$_{decomp}$ / N(HI)$_{Lee12}$ (1 $\times$ 10$^{20}$' + \
            'cm$^{-2}$')
plt.ylabel('Fraction of pixels')
plt.title('Frequency of differences of HI column density)')
plt.show()

# histogram of the differences between HI column densities calculated from
# integrating all emission in the GALFA cube used by Lee et al. (2012) between
# -5 and 15 km/s and the N(HI) from Lee et al. (2012)
from scipy.integrate import quad
fitsfile = 'perseus.lee12.cube.galfaBin.4arcmin.fits'
f = pf.open(fitsfile)
h,galfaLeeCube = f[0].header, f[0].data
dVel = h['cdelt3'] / 1000. # km/s
velRange = len(galfaLeeCube[:,0,0]) * dVel
galfaLeeNHI = galfaLeeCube.sum(axis=0) * 1.823 * dVel * 1e-2
galfaLeeNHI[~np.isnan(galfaLeeNHI)].max()
diff = galfaLeeNHI - (leeNHI / 1e20)
plt.hist(diff[~np.isnan(diff)].ravel(),bins=50,normed=1)
plt.rc('text', usetex=True)
plt.xlabel(r'N(HI)$_{int}$ / N(HI)$_{Lee12}$ (1 $\times$ 10$^{20}$' + \
            'cm$^{-2}$)')
plt.ylabel('Fraction of pixels')
plt.title('Frequency of percent differences of N(HI) from \n' + \
          'Lee et al. (2012) cube and Lee et al. (2012) N(HI) map')
plt.savefig('/d/bip3/ezbc/perseus/figures/' + \
            'perseus.lee12.leeGalfaInt.frac.histogram.png')
plt.show()

cfa = grid.load_grid('perseus.cfa.138_62.5sigma')

# trying with different noise profile scaling
box=[36,0,180,180]
reload(grid)
galfa0 = grid.SpectralGrid('../galfa/perseus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=30.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)
guesses = [48,5,10]
galfa0.fit_profiles(
        growPos = (138,62),
        tileSaveFreq=100,
        threshold = 1,
        filename='perseus.galfa.138_62.noiseX30',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=1e6,
        COcube=cfa,
        COwidthScale=1.)



################################################################################
# Loading cfa grid of persues
################################################################################

import grid
cfa = grid.load_grid('perseus.cfa.138_62.5sigma')

grid.plot_ncompImage(cfa)

box=[0,36,180,180]
reload(grid)
perseus_galfa = grid.SpectralGrid('../galfa/'+\
        'perseus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=10.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)
guesses = [48,5,10]
perseus_galfa.fit_profiles(
        growPos = (138,62),
        tileSaveFreq=100,
        threshold = 1,
        #filename='perseus.galfa.138_62.10',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=1e6,
        COcube=cfa,
        COwidthScale=1.)

perseus_galfa = grid.load_grid('perseus.galfa.138_62.10')

grid.plot_ncompImage(cfa,savedir='/d/bip3/ezbc/perseus/figures/',
                     filename='perseus.cfa.ncomps.png')

grid.plot_ncompImage(perseus_galfa,savedir='/d/bip3/ezbc/perseus/figures/',
                     filename='perseus.galfa.ncomps.png')
grid.plot_NHI(perseus_galfa,velocityrange=[-5,15],
              savedir='/d/bip3/ezbc/perseus/figures/',
              filename='perseus.galfa.NHI.velrange_neg5pos15.png',)
grid.plot_fit(perseus_galfa,93,138,coords='image',
              savedir='/d/bip3/ezbc/perseus/figures/',
              filename='perseus.galfa.spectrum.93_138.png')
grid.plot_fit(perseus_galfa,138,77,coords='image',
              savedir='/d/bip3/ezbc/perseus/figures/',
              filename='perseus.galfa.spectrum.138_77.png')


# histogram of the differences between HI column densities calculated from the
# decomposition and Lee et al. (2012).
galfaNHI = grid.plot_NHI(perseus_galfa,velocityrange=[-5,15],returnimage=True,
        show=False)
fitsfile='../galfa/perseus.lee12.galfaBin.4arcmin.fits'
f = pf.open(fitsfile)
h,leeNHI = f[0].header, f[0].data
diff = galfaNHI - leeNHI / 1e20
div = galfaNHI / (leeNHI / 1e20 *1.1)
plt.clf()
plt.hist(div[~np.isnan(div)].ravel(),bins=100,normed=False)
plt.rc('text', usetex=True)
plt.xlabel(r'N(HI)$_{decomp}$ / N(HI)$_{Lee12}$')
            #(1 $\times$ 10$^{20}$' + 'cm$^{-2}$')
plt.ylabel('Number of pixels')
plt.title('Frequency of differences of HI column density)')
plt.savefig('../../figures/perseus.histogram.lee12.galfaDecomp.fraction.png')
plt.show()

plot_fractionalDiff_NHI(perseus_galfa,leeNHI[0,:,:]/1e20*1.1,
        velocityrange=[-5,15],show=True,
        savedir='/d/bip3/ezbc/perseus/figures/',
        filename='perseus.image.lee12.galfaDecomp.fraction.png')


grid.plot_fits(perseus_galfa,(96,97,98,99),(96,96,96,96),coords='image',
                coRegion=False)


def plot_fractionalDiff_NHI(grid, image2, velocityrange=[],
        plotallpixels=False, savedir='./', filename=None, show=True,
        returnimage=False):
    ''' Shows and image of the fractional difference between a SpectralGrid
    calculated using components within the velocity range and image2. These
    must have the same shape.
    '''

    from kapteyn import wcs as kwcs
    import numpy as np
    import math
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs

    if len(velocityrange) != 2:
        sp = grid.spectralAxis
        velocityrange = [min(sp),max(sp)]

    fig = plt.figure(figsize=(8,8))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=grid.header)),
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    image = np.empty((grid.get_imagexsize(),grid.get_imageysize()))
    image[:,:] = np.NaN
    region = grid.region
    tilecount=0
    for i in xrange(region[0],region[2]):
        for j in xrange(region[1],region[3]):
            tile = grid.get_tile(i,j,coords='pixel')
            if tile.ncomps > 0:
                for k in range(tile.ncomps):                    
                    lowVel = tile.params[k][1] - tile.params[k][2]
                    highVel = tile.params[k][1] + tile.params[k][2]
                    if highVel < velocityrange[1] and lowVel > velocityrange[0]:
                        image[i,j] = tile.compAreaList[k]
            else:
                image[i,j] = np.NaN
            tilecount = tilecount + 1

    #print tilecount

    NHIimage = np.ma.array(image,mask=np.isnan(image)) * 1.823e-2

    #NHIimage[NHIimage == NHIimage.max()] = 1

    ax = imagegrid[0]
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")


    ax.set_xlabel('Right Ascension (J2000)', 
              size = 'small', 
              family='serif')
    ax.set_ylabel('Declination (J2000)', 
              size = 'small',
              family='serif')
    if not plotallpixels:
        ax.set_xlim(grid.region[0],grid.region[2])
        ax.set_ylim(grid.region[1],grid.region[3])


    #cmap = plt.cm.get_cmap('rainbow',NHIimage.max()-NHIimage.min()+1)
    #cmap = ax.cm.get_cmap('gray')
    #cmap.set_bad('w',1.)
    #ax.imshow(NHIimage, interpolation='nearest',cmap=cmap,origin='lower')    #cbar = ax.colorbar()
    #cbar.set_label(r'N(H\,\textsc{i}) $\times$ 10$^{18}$ cm$^{-2}$')

    image2 = np.ma.array(image2,mask=np.isnan(image2))

    print NHIimage.shape
    print image2.shape

    plt.rc('text', usetex=True)
    im = ax.imshow(NHIimage / image2, interpolation='nearest',origin='lower')
    cb = ax.cax.colorbar(im)
    # Write label to colorbar
    cb.set_label_text(r'N(HI)$_{decomp}$ / N(HI)$_{Lee12}$',
                   size='small',
                   family='serif')

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        plt.show()
    if returnimage:
        return NHIimage





