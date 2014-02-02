##########################################################
# Calibration script for project 13A-026
# Leo P D configuration observations

##########################################################
# 1: Prepare MS
##########################################################

# Important observation details:

# Field 0 = J1021+2159      Phase calibrator (disregard)
# Field 1 = J1021+2159	    Phase calibrator
# Field 2 = Leo P           Source
# Field 3 = 1331+305=3C286  Primary Calibrator

# We are dropping continnuum observations for now
# Our relevant spectral window:
# SpwID  #Chans Frame Ch1(MHz)    ChanWid(kHz)  TotBW(kHz)  Corrs          
# 4        1024 TOPO  1417.06793  3.90625       4000        RR  LL  

# We will now split out the spectral line spectral window

# Wait, why the hell do we have 4000 kHz of bandwidth?
# This is absurd. I'm cutting it down. Other observations have
# 1000 kHz of bandwidth, lets do that. UV - continuum subtraction
# may be less robust but dealing with that many channels will take
# too long
# Using middle 256 channels

from os import system,chdir

chdir('/d/bip3/ezbc/leop/data/hi/casa/reductionFiles/vla.d')

system('mv 13A-026.sb14766460.eb19475734.56374.04000275463.ms 13a-026.ms')
system('mv 13A-026.sb14766460.eb19475734.56374.04000275463.ms.flagversions \
        13a-026.ms.flagversions')

split(vis='13a-026.ms',
      outputvis='13a-026.spectralLine.ms',
      datacolumn='data',
      spw='4:384~640')

# Define visibility name as something easy

myVis='13a-026.spectralLine.ms'


##########################################################
# 2: Flag the data 
##########################################################

# # First examined data in plotms
# checking short baselines for interference from the sun
# Observations finished at 5:07 AM Mar 23rd, sun rose at 7:05 AM according to google

# 

flagdata(vis=myVis,
	 antenna='ea07&ea15;ea13&ea15',
	 correlation='',
	 scan='32~35',
	 timerange='')

flagdata(vis=myVis,
         antenna='ea13;ea05',
         correlation='',
         scan='',
         timerange='2013/03/23/03:33:57.5~2013/03/23/03:34:57.5')

# ea13 misbehaving

flagdata(vis=myVis,
         antenna='ea13',
         correlation='',
         scan='',
         timerange='')

flagdata(vis=myVis,
         timerange='2013/03/23/05:04:17.5~2013/03/23/05/04:27.5')


# Did I forget how to flag? I didn't see one thing...



###########################################################
# 3: Set the flux scale
###########################################################

setjy(vis=myVis,
      field='3',
      modimage='3C286_L.im')

# 3C286 flux at 14.76 Jy

###########################################################
# 4: Acquire phase and bandpass 
###########################################################

# First apply ant position corrections:

gencal(vis=myVis,
       caltable='antpos.gcal',
       caltype='antpos')

# no offsets found, be sure to pay attention to this later

# Derive phase solutions for bandpass

gaincal(vis=myVis,
        caltable='bpphase.gcal',
        field='3',
        refant='ea03',
	calmode='p')

bandpass(vis=myVis,
         caltable='bandpass.bcal',
	 field='3',
         refant='ea03',
	 solint='inf',
	 solnorm=T,
         gaintable=['bpphase.gcal'],
	 gaincurve=T)

# ea25 has some high gain amps, i.e. 1.6

gaincal(vis=myVis,
        caltable='intphase.gcal',
        field='1,3',
	spw='0',
        refant='ea03',
	calmode='p',
	solint='int',
        gaintable=['bandpass.bcal'],
        gaincurve=T)

gaincal(vis=myVis,
        caltable='scanphase.gcal',
        field='1,3',spw='0',
        refant='ea02',
	calmode='p',
	solint='inf',
        gaintable=['bandpass.bcal'],
        gaincurve=T)


#apply phase solutions to get amp solutions for calibrators

gaincal(vis=myVis,
        caltable='amp.gcal',
        field='1,3',
	spw='0',
        refant='ea03',
	calmode='ap',
	solint='inf',
        gaintable=['bandpass.bcal','scanphase.gcal'],
        gaincurve=T)

#now set the fluxscale
fluxscale(vis=myVis,
          caltable='amp.gcal',
          fluxtable='flux.cal',
	  reference='3')

# J1021_2159 flux at 1.688 +/- 0.012 Jy
# VLA calibrator manual gives 0.6 Jy at 7cm

###########################################################
# 5: Apply calibrations
###########################################################

# Apply to primary calibrator

applycal(vis=myVis,
         field='3',
         gaintable=['bandpass.bcal','scanphase.gcal','flux.cal'],
         gainfield=['3','3','3'],
         gaincurve=T,
	 calwt=F)

# Apply to phase calibrator

applycal(vis=myVis,
         field='1',
         gaintable=['bandpass.bcal','intphase.gcal','flux.cal'],
         gainfield=['3','1','1'],
         gaincurve=T,
	 calwt=F)

# Apply to target

applycal(vis=myVis,
	 field='2',
         gaintable=['bandpass.bcal','scanphase.gcal','flux.cal'],
         gainfield=['3','1','1'],
         gaincurve=T,
	 calwt=F)

# Put the visibilities back together by first splitting out source

# Flagging time

flagdata(vis=myVis,
	 timerange='2013/03/23/05:04:17.5~2013/03/23/05/04:27.5')


# Successful calibration!

flagdata(vis=myVis,
	 scan='36~40',
	 antenna='ea07&ea15')

split(vis=myVis,
      field='2',
      outputvis='leopVLAd.ms',
      datacolumn='corrected')

myVis = 'leopVLAd.ms'

uvdir = '/d/bip3/ezbc/leop/data/hi/finalUVdatasets/'

# change the source name of each dataset to be the same
tb.open(myVis + '/FIELD',nomodify=False)
st=tb.selectrows(0)
st.putcol('NAME','LeoP')
st.done()
tb.close()
# put the three datasets on the same grid
myVisLSR = 'leopVLAd.lsr.ms'
cvel(vis=myVis,
     outputvis=myVisLSR,
     outframe='lsrk')
# export the fits file
exportuvfits(vis=myVisLSR,
             fitsfile=uvdir + 'leopVLAd.lsr.continuum.fits',
             multisource=False)

###########################################################
# 6: Subtract the continuum
###########################################################

# Make a dirty cube to see which channels have emission

#clean(vis='leopVLA-d.ms',
#      imagename='../images/dirty.v1',
#      mode='velocity',
#      imsize=512,
#      cell='3arcsec',
#      niter=1000,
#      threshold='2mJy',
#      restfreq='1.420405752GHz')

# Channels 1~90 and 150~254 are source free

uvcontsub(vis=myVisLSR,
	  fitspw='0:1~90;150~254',
	  fitorder=1)

# split off relevant channels

#split(vis='leopVLAd.ms.contsub',
#      outputvis='leopVLAd.ms.contsub')

# export as a uvfits file:

exportuvfits(vis=myVisLSR + '.contsub',
             fitsfile=uvdir + 'leopVLAd.lsr.contsub.fits')

# image the phase calibrator
default(clean)
vis = '13a-026.spectralLine.ms'
field = '1331+305=3C286'
mode = 'mfs'
imagename = 'leopVLAd.fluxCablibrator'
niter = 1000 # for light cleaning
threshold = '3mJy' # see note1
imsize = 512
cell = '9arcsec' # see note2
weighting = 'briggs'
robust = 0.5
clean()

#no
# image the continuum
default(clean)
vis = myVisLSR
imagename = 'leopVLAd.continuum'
mode = 'velocity'
restfreq = '1420.405752MHz'
niter = 0
imsize = 240
cell = '9arcsec'
weighting = 'briggs'
robust = 0.5
clean()



# Image the galaxy

# clean(vis='leopVLA-d.60chan.ms.contsub',
#      imagename='../images/leop.v1',
#      mode='velocity',
#      imsize=512,
#      cell='2arcsec',
#      threshold='3mJy',
#      restfreq='1.420405752GHz',
#      interactive=True)

# interactive mode does not allow for deconvolution. Screw CASA
# leop.v1 had emission in first channel, try again. Deleted leop*.ms
# files, changed uv cont sub channels

#clean(vis='leopVLA-d.60chan.ms.contsub',
#      imagename='../images/leop.v2',
#      mode='velocity',
#      imsize=512,
#      cell='2arcsec',
#      threshold='3mJy',
#      restfreq='1.420405752GHz')
#
# leop.v2 had emission in first channel, try again

#clean(vis='leopVLA-d.60chan.ms.contsub',
#      imagename='../images/leop.v3',
#      mode='velocity',
#      imsize=512,
#      cell='2arcsec',
#      threshold='3mJy',
#      restfreq='1.420405752GHz')
#
# leop.v3 had emission in first channel, try again


'''
clean(vis='leopVLA-d.60chan.ms.contsub',
      imagename='../images/leop.v4',
      mode='velocity',
      imsize=512,
      cell='2arcsec',
      threshold='3mJy',
      restfreq='1.420405752GHz')

# Success! ~1.8 mJy noise in source-free channels
# images in /Users/research/data/leop/vla/d_config/images/partialspwImages/

# Now smooth to a circular beam
# Beam size: 61.2854 arcsec, 52.2832 arcsec
# Convolve to 62 arcsec

# define image filepath variable

imDir='../images/partialspwImages/'

imsmooth(imagename=imDir + 'leop.v4.image',
         outfile=imDir + 'leop.v4.image.smooth',
         major='62arcsec',
         minor='62arcsec',
         targetres=True)

# Now we export the image to GIPSY to create moment maps!

exportfits(imagename=imDir + 'leop.v4.image.smooth',
           velocity=True,
           fitsimage=imDir + 'leop.fits')

# Need more channels for a spectrum

clean(vis='leopVLA-d.ms.contsub',
      imagename=imdir + 'leop.v5',
      mode='velocity',
      imsize=512,
      cell='2arcsec',
      threshold='4mJy',
      restfreq='1.420405752GHz')


# create a moment map

immoments(imagename=imDir + 'leop.v4.image.smooth',
          moments=[0],
          axis='spectral',
          chans='3~50',
          outfile=imDir + 'leop.v1.mom0')


'''




