#!/usr/bin/python

from os import chdir

chdir('/d/bip3/ezbc/leop/data/hi/casa/reductionFiles/vla.c')

#1. Examine the data and flag

flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea23',correlation='LL')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea05&ea13')
##
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea01&ea05',correlation='RR',timerange='23:43:12.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea01&ea25',correlation='RR',timerange='23:43:17.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea01&ea06',correlation='LL',timerange='03:24:37.5~03:26:02.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea05&ea26',scan='29~49')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea13&ea14',timerange='03:24:37.5~03:25:12.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea15&ea18',correlation='RR',timerange='23:43:12.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea19&ea25',correlation='RR',timerange='23:43:12.5~23:43:22.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea04&ea14',spw='0:181')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea25&ea26',correlation='RR',timerange='23:43:12.5~23:43:22.5')

####
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea01&ea06',correlation='LL',timerange='03:24:27.5~03:26:02.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea05&ea14',scan='45')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea05&ea25',correlation='RR',timerange='23:43:12.5~23:43:17.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea05&ea26',correlation='RR',timerange='23:43:12.5~23:43:17.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea13&ea14',timerange='03:24:37.5~03:26:02.5,03:54:17.5~03:54:32.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea19&ea26',correlation='RR',timerange='23:43:12.5~23:43:17.5')

########
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea01&ea06',correlation='LL',timerange='03:10:27.5~03:12:47.5,03:24:27.5~03:26:02.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea05&ea14',timerange='03:10:42.5~03:12:57.5,03:24:12.5~03:24:17.5,03:35:47.5~03:36:32.5,03:44:57.5~03:45:32.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         antenna='ea13&ea14',timerange='03:10:22.5~03:12:57.5,03:24:02.5~03:25:57.5,03:35:47.5~03:36:32.5,03:44:57.5~03:45:32.5')

#############################

#3. Set the flux scale
setjy(vis='12A-456_sb9233552_1.56035.96345914352.ms',field='3',
      modimage='3C286_L.im')



flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea13&ea14',scan='45')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea13&ea26',scan='37,45')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea05&ea25',timerange='23:43:22.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea05&ea26',timerange='23:43:22.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea07&ea26',timerange='23:43:12.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea19&ea25',timerange='23:43:27.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea19&ea26',timerange='23:43:22.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea25&ea26',timerange='23:43:27.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='LL',antenna='ea13&ea14',scan='45')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='LL',antenna='ea01&ea06',timerange='03:25:02.5~03:25:57.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='LL',antenna='ea14&ea19',scan='50,51')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea13&ea19',scan='50,51')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='LL',antenna='ea05&ea14',timerange='03:11:22.5~03:12:17.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea05&ea14',timerange='03:11:22.5~03:12:17.5')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='LL',antenna='ea05&ea14',scan='43,47,48')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea05&ea14',scan='43,47,48')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='LL',antenna='ea13&ea26',scan='43,47,48')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea13&ea26',scan='43,47,48')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='LL',antenna='ea13&ea14',scan='43,47,48')
flagdata(vis='12A-456_sb9233552_1.56035.96345914352.ms',correlation='RR',antenna='ea13&ea14',scan='43,47,48')

gaincal(vis='12A-456_sb9233552_1.56035.96345914352.ms',
        caltable='bpphaseredo.gcal',
        field='3',spw='0:118~138',
        refant='ea02',calmode='p',solint='int',minsnr=2.0,
        opacity=0.,gaincurve=T)

bandpass(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         caltable='bandpassredo.bcal',field='3',
         refant='ea02',solint='inf',solnorm=T,
         gaintable=['bpphaseredo.gcal'],
         opacity=0.0,gaincurve=T)

gaincal(vis='12A-456_sb9233552_1.56035.96345914352.ms',
        caltable='intphaseredo.gcal',
        field='1,3',spw='0:30~225',
        refant='ea02',calmode='p',solint='int',minsnr=2.0,
        gaintable=['bandpassredo.bcal'],
        opacity=0.0,gaincurve=T)

gaincal(vis='12A-456_sb9233552_1.56035.96345914352.ms',
        caltable='scanphaseredo.gcal',
        field='1,3',spw='0:30~225',
        refant='ea02',calmode='p',solint='inf',minsnr=2.0,
        gaintable=['bandpassredo.bcal'],
        opacity=0.0,gaincurve=T)

#apply phase solutions to get amp solutions for calibrators
gaincal(vis='12A-456_sb9233552_1.56035.96345914352.ms',
        caltable='ampredo.gcal',
        field='1,3',spw='0:30~225',
        refant='ea02',calmode='ap',solint='inf',minsnr=2.0,
        gaintable=['bandpassredo.bcal','intphaseredo.gcal'],
        opacity=0.0,gaincurve=T)

#now set the fluxscale
fluxscale(vis='12A-456_sb9233552_1.56035.96345914352.ms',
          caltable='ampredo.gcal',
          fluxtable='fluxredo.cal',reference='3')

fluxscale(vis='12A-456_sb9233552_1.56035.96345914352.ms',
          caltable='ampredo.gcal',
          fluxtable='fluxredo.cal',reference='3')

#apply to flux cal
applycal(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         field='3',
         gaintable=['bandpassredo.bcal','intphaseredo.gcal','fluxredo.cal'],
         gainfield=['3','3','3'],
         opacity=0.0,gaincurve=T,calwt=F)

#apply to phase cal
applycal(vis='12A-456_sb9233552_1.56035.96345914352.ms',
         field='1',
         gaintable=['bandpassredo.bcal','intphaseredo.gcal','fluxredo.cal'],
         gainfield=['3','1','1'],
         opacity=0.0,gaincurve=T,calwt=F)

#apply to target
applycal(vis='12A-456_sb9233552_1.56035.96345914352.ms',field='2',
        gaintable=['bandpassredo.bcal','scanphaseredo.gcal','fluxredo.cal'],
        gainfield=['3','1','1'],
        opacity=0.0,gaincurve=T,calwt=F)

#put the visibilities back together by first splitting out source

split(vis='12A-456_sb9233552_1.56035.96345914352.ms',field='2',
      outputvis='leopVLAc.ms',datacolumn='corrected')

uvdir = '/d/bip3/ezbc/leop/data/hi/casa/uvdata/'

myVis = 'leopVLAc.ms'

# change the source name of each dataset to be the same
tb.open(myVis + '/FIELD',nomodify=False)
st=tb.selectrows(2)
st.putcol('NAME','LeoP')
st.done()
tb.close()
# put the three datasets on the same grid
myVisLSR = 'leopVLAc.lsr.ms'
cvel(vis=myVis,
     outputvis=myVisLSR,
     outframe='lsrk')

uvcontsub(vis=myVisLSR,
          fitspw='0:10~60;160~240',
          fitorder=1)

uvfitsdir = '/d/bip3/ezbc/leop/data/hi/finalUVdatasets/'

exportuvfits(vis=myVisLSR + '.contsub',
             padwithflags=False,
             multisource=False,
             fitsfile=uvfitsdir + 'leopVLAc.lsr.contsub.fits')

exportuvfits(vis=myVisLSR,
             padwithflags=False,
             multisource=False,
             fitsfile=uvfitsdir + 'leopVLAc.lsr.continuum.fits')


#exportuvfits(vis='leopVLAc.ms.contsub',
#             fitsfile='leopVLAc.contsub.fits')


# image the flux calibrator
default(clean)
vis = '12A-456_sb9233552_1.56035.96345914352.ms'
field = '1331+305=3C286'
mode = 'mfs'
imagename = 'leopVLAc.fluxCabibrator'
niter = 1000 # for light cleaning
threshold = '2mJy' # see note1
imsize = 512
cell = '2.8arcsec' # see note2
weighting = 'briggs'
robust = 0.5
clean()



# image the continuum
default(clean)
vis = myVisLSR
imagename = 'leopVLAc.continuum'
mode = 'velocity'
restfreq = '1420.405752MHz'
niter = 0
imsize = 1024
cell = '2arcsec'
weighting = 'briggs'
robust = 0.5
clean()







