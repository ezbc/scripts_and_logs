#!/usr/bin/python

guesses = [50,0,10,20,-5,10,10,-10,10,30,-20,50,10,-10,10,30,-5,10,10,-20,10]

###############
# 16' cube
###############

reload(grid)

x16 = grid.SpectralGrid('taurus.galfa.cube.bin.16arcmin.fits',
                        box=[40,40,120,120],
                        noiseScale=40.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)

guesses = [  4.07575984, -26.97410452,  16.75406491,
 53.70105613,   4.56270824,   4.9002332 ]

x16.fit_profiles(
    growPos = (85,85),
    tileSaveFreq=10,
    threshold = 1,
    filename='grid.16arcsec.85_85.noiseX30',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    coords='image',
    numberOfFits=10000)

grid.plot_spectrum(x16,85,85,coords='image')
grid.plot_fit(x16,85,85,coords='image')
grid.plot_fit(x16,85,84,coords='image')
grid.plot_fits(x16,(85,85,85,85),(83,84,85,86),coords='image')
grid.plot_fits(x16,(86,86,86,86),(83,84,85,86),coords='image')
grid.plot_fits(x16,(30,30,30,30),(29,30,31,32))
grid.plot_residualsImage(x16)
grid.plot_componentPV(x16,yslice=40,width=1)
grid.plot_ncompImage(x16)
grid.plot_velocityHistograms(x16)


reload(grid)

x16 = grid.SpectralGrid('taurus.galfa.cube.bin.16arcmin.fits',
                        box=[60,60,120,120],
                        noiseScale=50.,
                        noiseRange=((-110,-80),(80,110)),
                        basesubtract=True)

guesses = [  4.07575984, -26.97410452,  16.75406491,
 53.70105613,   4.56270824,   4.9002332 ]

x16.fit_profiles(
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

x16 = grid.SpectralGrid('taurus.galfa.cube.bin.16arcmin.fits',
                        box=[60,60,120,120],
                        noiseScale=70.)

guesses = [  4.07575984, -26.97410452,  16.75406491,
 53.70105613,   4.56270824,   4.9002332 ]

x16.fit_profiles(
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

x16.fit_profile(guesses,len(guesses)/3,25,25)


grid.plot_fit(x16,80,86,coords='image')
grid.plot_fits(x16,(84,85,85,85),(85,84,85,86),coords='image')
grid.plot_fits(x16,(85,85,85,85),(83,84,85,86),coords='image')
grid.plot_fits(x16,(80,80,80,80),(83,84,85,86),coords='image')

grid.plot_fits(x16,(5,5,5,5),(1,2,3,4),,coords='grid')

grid.plot_fit(x16,5,3)
grid.plot_fit(x16,85,83,coords='image')

grid.plot_fits(x16,(79,79,79,79),(74,75,76,77),coords='image')
grid.plot_fits(x16,(71,71,71,71),(74,75,76,77),coords='image')

grid.plot_residualsImage(x16)
grid.plot_componentPV(x16,yslice=25,width=3)
grid.plot_ncompImage(x16)
grid.plot_velocityHistograms(x16)

# for 75_75.5sigma.01alpha we see ncomp varying by 6 across (73,73:76)
grid.plot_fits(x16,(73,73,73,73),(73,74,75,76),coords='image')

# Try with initial guesses
x16 = grid.SpectralGrid('taurus.galfa.cube.bin.16arcmin.fits',box=[80,80,90,90])

# Save the grid
x16.write_grid('grid.16arcsec')

# Load the grid
# Allows changes to be made in SpectralGrid class without having to recreate the instance
x16 = grid.SpectralGrid('taurus.galfa.cube.bin.16arcmin.fits',gridfile='grid.16arcsec')


###############
# 4' cube
###############

reload(grid)

galfa = grid.SpectralGrid('../galfa/taurus.galfa.cube.bin.4arcmin.fits',
                        box=[300,300,400,400],
                        noiseScale=40.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)

galfa2 = grid.SpectralGrid('../galfa/taurus.galfa.cube.bin_4arcmin.fits',
                        box=[360,360,390,390],
                        noiseScale=40.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)


guesses = [8,-20,30, 45,10,15]

galfa2.fit_profile(guesses,len(guesses)/3,390,390,coords='image',COcube=cfa,
                    COwidthScale=3000.)

guesses = [4,-20,30, 45,10,15]

galfa2.fit_profile(guesses,len(guesses)/3,373,383,coords='image',COcube=cfa,
                    COwidthScale=5000.)

grid.plot_spectrum(galfa,373,383,coords='image')
guesses = [35,8,7, 5,-20,15]
galfa.fit_profile(guesses,len(guesses)/3,373,383,coords='image')
grid.plot_fit(galfa,373,383,coords='image')


grid.plot_spectrum(galfa,325,325,coords='image')
guesses = [4,-30,60, 45,10,20, 4,-45,6]
galfa.fit_profile(guesses,len(guesses)/3,325,325,coords='image')
grid.plot_fit(galfa,325,325,coords='image')


grid.plot_spectrum(galfa,325,475,coords='image')
guesses = [20,10,20, 5,5,50]
galfa.fit_profile(guesses,len(guesses)/3,325,475,coords='image')
grid.plot_fit(galfa,325,475,coords='image')


grid.plot_spectrum(galfa,475,325,coords='image')
guesses = [50,5,10, 10,-15,30, 6,-50,20]
galfa.fit_profile(guesses,len(guesses)/3,475,325,coords='image')
grid.plot_fit(galfa,475,325,coords='image')



grid.plot_fit(galfa2,390,390,coords='image')
grid.plot_fit(galfa2,373,383,coords='image')


galfa2.fit_profiles(
        growPos = (373,383),
        tileSaveFreq=10,
        threshold = 1,
        #filename='galfa.390_390',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=10000,
        COcube=cfa,
        COwidthScale=3000.)


cfa = grid.load_grid('cfa.373_383')

reload(grid)
guesses = [4,-20,30, 45,10,15]
box=[360,370,380,390]


galfa = grid.SpectralGrid('../galfa/taurus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=40.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)

galfa.fit_profiles(
        growPos = (375,388),
        tileSaveFreq=10,
        threshold = 1,
        #filename='galfa.390_390',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=100,
        COcube=cfa,
        COwidthScale=3.)

grid.plot_ncompImage(galfa)


reload(grid)
guesses = [4,-20,30, 45,10,15]
box=[360,310,380,330]


galfa0 = grid.SpectralGrid('../galfa/taurus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=40.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)

guesses = [4,-30,60, 45,10,20, 4,-45,6]
galfa0.fit_profiles(
        growPos = (368,317),
        tileSaveFreq=20,
        threshold = 1,
        #filename='galfa0',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=1e6,
        COcube=cfa,
        COwidthScale=3.)



grid.plot_spectrum(galfa3,390,390,coords='image')
grid.plot_fit(galfa1,373,383,coords='image')
grid.plot_fit(galfa0,368,317,coords='image')

grid.plot_fit(galfa2,85,84,coords='image')

grid.plot_fits(galfa1,(373,373,373),(382,383,384),coords='image')

grid.plot_fits(galfa1,(367,367,367,367),(433,434,435,436),coords='image')

grid.plot_fits(galfa,(367,367,367,367),(433,434,435,436),coords='image')
grid.plot_fits(galfa,(373,373,373),(382,383,384),coords='image')
grid.plot_fit(galfa,373,382,coords='image')


grid.plot_fits(galfa,(86,86,86,86),(83,84,85,86),coords='image')

grid.plot_fits(galfa,(30,30,30,30),(29,30,31,32))

grid.plot_residualsImage(galfa)

grid.plot_componentPV(galfa1,yslice=73,width=3)

grid.plot_ncompImage(galfa3)
grid.plot_ncompImage(galfa1)
grid.plot_NHI(galfa1)


grid.plot_velocityHistograms(galfa)

reload(grid)
box=[300,300,500,500]

galfa1 = grid.SpectralGrid('../galfa/taurus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=1.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)

guesses = [30,10,20, 5,-5,50]
galfa1.fit_profiles(
        growPos = (325,400),
        tileSaveFreq=500,
        threshold = 1,
        #filename='galfa1',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=100,
        COcube=cfa,
        COwidthScale=1.)

grid.plot_fit(galfa1,325,401,coords='image')
grid.plot_fits(galfa1,(326,325,324),(401,401,401),coords='image')









################################################################################################################################################################
################################################################################
################################################################################
# Performed fits to GALFA HI cube using decomposed CfA CO cube to disclude
# velocities in fit
# script at:
# /d/bip3/ezbc/taurus/logFiles/taurus.gaussianDecomposition.fits0.py
################################################################################
################################################################################
################################################################################
################################################################################


################################################################################
# Analysis of grids from three separate starting positions
################################################################################
cfa = grid.load_grid('cfa.373_383')

# Grid    | Init Pos
# galfa0    325,325
# galfa1    325,400
# galfa2    380,390
# galfa3    373,383

galfa0 = grid.load_grid('galfa0')
galfa1 = grid.load_grid('galfa1')
galfa2 = grid.load_grid('galfa2')
galfa3 = grid.load_grid('galfa3')

grid.compare_grids(galfa0,galfa2,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.residualNHI.galfa0_2.png')
grid.compare_grids(galfa0,galfa1,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.residualNHI.galfa0_1.png')
grid.compare_grids(galfa2,galfa1,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.residualNHI.galfa2_1.png')

# area around 319,372 is troubled, are these complicated spectra?

grid.plot_fit(galfa1,319,372,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.galfaHI.spectrum.319_372.png')
grid.plot_fit(galfa0,319,372,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.galfaHI.spectrum2.319_372.png')

# saved images as 

# indeed these are much different, a large component is missing in galfa1

# redoing CfA fit 
reload(grid)

box=[300,300,500,500]
#box=[360,360,400,400]
cfa = grid.SpectralGrid('../cfa/taurus.cfa.galfaBin.4arcmin.fits',
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

grid.plot_spectrum(cfa,373,383,coords='image')

guesses = [4,8,2]
cfa.fit_profile(guesses,len(guesses)/3,373,383,coords='image')

grid.plot_fit(cfa,373,383,coords='image')

cfa.fit_profiles(
    growPos = (373,383),
    tileSaveFreq=400,
    threshold = 5,
    filename='cfa.373_383.5sigma',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=1e6)

grid.plot_residualsImage(cfa)
grid.plot_ncompImage(cfa)
grid.plot_velocityHistograms(cfa)
grid.plot_fit(cfa,375,391,coords='image')


# Decompose GALFA HI cube
reload(grid)
galfa0 = grid.SpectralGrid('../galfa/taurus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=5.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)

guesses = [35,8,7, 5,-20,15]
galfa0.fit_profiles(
        growPos = (373,383),
        tileSaveFreq=500,
        threshold = 1,
        filename='galfa.373_383',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=1e6,
        COcube=cfa,
        COwidthScale=1.)


cfa = grid.load_grid('cfa.373_383.5sigma')
galfa0 = grid.load_grid('galfa0') # used CO width=1
galfa1 = grid.load_grid('galfa1')
galfa2 = grid.load_grid('galfa2')
galfa3 = grid.load_grid('galfa3')

grid.plot_residualsImage(galfa0)
grid.plot_ncompImage(galfa0,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.ncomp.galfa0.noise5X.png')
grid.plot_NHI(galfa0,velocityrange=[-5,15],
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.nhi.galfa0.noise5X.png')
# position 352,392 has 18 x 10^20 cm^-2, plot fit:
grid.plot_fit(galfa0,352,392,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.352_392.galfa0.noise5X.png')
# CO region excluded has width of ~100 km/s. wtf?
grid.plot_fit(cfa,352,392,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.352_392.cfa.png')
grid.plot_fit(cfa,352,393,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.352_393.cfa.png')
grid.plot_fit(galfa0,352,393,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.352_393.galfa0.png')
# glitch in code, the following are the shifts and widths of the CO spectrum
# shift        width
# 8.6357773 ,  0.44782122
# 5.73345853,  1.1119888 

# compare to galfa1
grid.plot_fit(galfa1,352,392,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.352_392.galfa1.noise5X.png')
galfa1.get_tile(352,392,coords='image').coVelocities

# checking other high column density pixel
grid.plot_fit(galfa0,446,389,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.446_389.galfa0.noise5X.png')
# same story, CO extent is ~50 km/s, centered on -10 km/s
grid.plot_fit(galfa0,446,389,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.446_389.galfa0.noise5X.png')
# check CO pixel
grid.plot_fit(cfa,446,389,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.446_389.cfa.png')
# multiple components in both problemed pixels, is this the problem
grid.plot_fit(galfa0,446,387,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.446_387.galfa0.noise5X.png')
# nope, this is very reasonable, though consider using only one CO width, and
# discarding velocities for each component.
grid.plot_fit(cfa,446,387,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.446_387.cfa.png')

# attempt new decomposition, discarding velocities in HI fit based on widths of
# each CO component
# Decompose GALFA HI cube
reload(grid)
box=[435,375,450,395]
galfa4 = grid.SpectralGrid('../galfa/taurus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=10.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)

grid.plot_spectrum(galfa4,445,386,coords='image')

guesses = [3,-70,20, 5,-50,10, 5,-25,20, 4,-25,3, 3,-15,3, 40,5,6]

galfa4.fit_profile(guesses,len(guesses)/3,445,386,coords='image')
grid.plot_fit(galfa4,445,386,coords='image')

galfa4.fit_profiles(
        growPos = (445,386),
        tileSaveFreq=100,
        threshold = 1,
        filename='galfa4.445_386',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=10,
        COcube=cfa,
        COwidthScale=1.)

grid.plot_fit(galfa4,446,387,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.spectrum.446_387.galfa4.png')
# huzzah! this is looking very good, three CO regions are clearly defined.

grid.plot_ncompImage(galfa4,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.ncomp.galfa4.noise5X.png')

################################################################################################################################################################
################################################################################
################################################################################
# Performed fits to GALFA HI cube using decomposed CfA CO cube to disclude
# velocities in fit
# script at:
# /d/bip3/ezbc/taurus/logFiles/taurus.gaussianDecomposition.fits1.py
################################################################################
################################################################################
################################################################################
################################################################################




























grid.plot_NHI(galfa1,velocityrange=[-5,15],
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.nhi.galfa1.noise5X.png')
grid.plot_velocityHistograms(galfa0)
grid.plot_componentPV(galfa0,yslice=100,width=4)
grid.plot_fit(galfa0,373,383,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.galfaHI.spectrum1.noiseX5.373_383.png')
grid.plot_fit(galfa0,374,384,coords='image',
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.galfaHI.spectrum2.noiseX5.374_384.png')

grid.compare_grids(galfa0,galfa1,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.residualNHI.galfa0_1.noise5X.png')
grid.compare_grids(galfa1,galfa2,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.residualNHI.galfa1_2.noise5X.png')
grid.compare_grids(galfa2,galfa3,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.residualNHI.galfa2_3.noise5X.png')
grid.compare_grids(galfa1,galfa3,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.residualNHI.galfa1_3.noise5X.png')
grid.compare_grids(galfa0,galfa3,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.residualNHI.galfa0_3.noise5X.png')
grid.compare_grids(galfa0,galfa2,
              savedir='/d/bip3/ezbc/taurus/figures/',
              filename='taurus.residualNHI.galfa0_3.noise5X.png')














################################################################################
# Attempt decomposition on CfA 12CO J 1-->0 cube
###############################################################################

from mirexec import TaskRegrid as regrid

reload(grid)

cfa = grid.SpectralGrid('../cfa/taurus.cfa.galfaBin.4arcmin.fits',
        noiseRange=((-0.0224e3,-0.01852915e3),(0.02438339e3,0.02958491e3)),
        basesubtract=False,
        box=[360,360,390,390])

cfa = grid.SpectralGrid('../cfa/taurus.cfa.galfaBin.4arcmin.fits',
        noiseRange=((-0.0224,-0.01852915),(0.02438339,0.02958491)),
        basesubtract=False,
        noiseScale=50.,
        box=[530,150,550,170])


guesses = [2,3,3, 6,8,5]

cfa.fit_profile(guesses,len(guesses)/3,390,390,coords='image')

cfa.fit_profile(guesses,len(guesses)/3,373,383,coords='image')


grid.plot_spectrum(cfa,390,390,coords='image')
grid.plot_spectrum(cfa,121,151,coords='image')


grid.plot_fit(cfa,390,390,coords='image')
grid.plot_fit(cfa,373,383,coords='image')

grid.plot_fit(cfa,370,379,coords='image')


cfa.fit_profiles(
    growPos = (373,383),
    tileSaveFreq=25,
    threshold = 1,
    filename='test',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=100)

cfa.fit_profiles(
    growPos = (549,151),
    tileSaveFreq=25,
    threshold = 1,
    #filename='cfa.390_390',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=10000)

grid.plot_fit(cfa,373,383,coords='image')
grid.plot_fit(cfa,534,190,coords='image')

grid.plot_spectrum(cfa,549,151,coords='image')
grid.plot_fit(cfa,549,151,coords='image')





grid.plot_residualsImage(cfa)
grid.plot_componentPV(cfa,yslice=30,width=1)
grid.plot_ncompImage(cfa)
grid.plot_velocityHistograms(cfa)


grid.plot_fits(cfa,(116,116,116,116),(153,154,155,156),coords='image')


reload(grid)

box=[360,370,380,390]

cfa = grid.SpectralGrid('../cfa/taurus.cfa.galfaBin.4arcmin.fits',
        noiseRange=((-0.0224e3,-0.01852915e3),(0.02438339e3,0.02958491e3)),
        basesubtract=False,
        box=box)

guesses = [2,3,3,6,8,5]

cfa.fit_profiles(
    growPos = (373,383),
    tileSaveFreq=400,
    threshold = 1,
    #filename='cfa.373_383',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=100)


reload(grid)

cfa = grid.SpectralGrid('../cfa/taurus.cfa.galfaBin.4arcmin.fits',
        noiseRange=((-0.0224e3,-0.01852915e3),(0.02438339e3,0.02958491e3)),
        basesubtract=False,
        box=[360,360,390,390])







###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Fitting large region of Taurus
###############################################################################
###############################################################################
###############################################################################

reload(grid)

box=[0,180,730,700]

cfa = grid.SpectralGrid('/d/bip3/ezbc/taurus/data/cfa/' + \
                        'taurus.cfa.galfaBin.4arcmin.fits',
        noiseRange=((-0.0224e3,-0.01852915e3),(0.02438339e3,0.02958491e3)),
        basesubtract=False,
        box=box,
        Tsys=400.)

guesses = [2,3,3,6,8,5]

cfa.fit_profiles(
    growPos = (373,384),
    tileSaveFreq=10000,
    threshold = 1,
    filename='taurus.cfa.373_384',
    guesses=guesses,
    ncomp=len(guesses)/3,
    alpha=0.001,
    verbose=True,
    coords='image',
    numberOfFits=1e100)

cfa = grid.load_grid('taurus.cfa.373_384')

reload(grid)

galfa = grid.SpectralGrid('/d/bip3/ezbc/taurus/data/galfa/' + \
        'taurus.galfa.cube.bin.4arcmin.fits',
                        box=box,
                        noiseScale=20.,
                        noiseRange=((-110,-90),(90,110)),
                        basesubtract=True)

guesses = [3,-70,20, 5,-50,10, 5,-25,20, 4,-25,3, 3,-15,3, 40,5,6]

galfa.fit_profiles(
        growPos = (445,386),
        tileSaveFreq=1000,
        threshold = 1,
        filename='galfa.445_386.3CO_widths',
        guesses=guesses,
        ncomp=len(guesses)/3,
        alpha=0.001,
        coords='image',
        numberOfFits=1e6,
        COcube=cfa,
        COwidthScale=3.)

galfa = grid.load_grid('galfa.445_386.3CO_widths')

grid.plot_NHI(galfa,velocityrange=[-10,20])
