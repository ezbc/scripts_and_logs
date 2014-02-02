
vis='12A-456.sb10554071.eb11466720.56157.681197384256.ms'

listobs(vis)

# flag bad antenna from observing log
tflagdata(vis,
	   field='',
	   spw='',
	   antenna='ea27',
	   timerange='2012/08/18/20:56:00~2012/08/18/21:20:01',
	   flagbackup=True)

# flag the quack scan
tflagdata(vis,
	  mode='quack',
	  quackinterval=10.0,
	  quackmode='beg')

# flag phase and primary calibrators, rfi found in viewer
tflagdata(vis,
	  field='2',
          spw='',
	  timerange='2012/08/18/21:07:00~2012/08/18/21:08:42',
	  antenna='',
	  flagbackup=True)

# autoflag
  testautoflag(vis,
	       datacolumn='data',
	       ntime=100,
	       timecutoff=4.0,
	       freqcutoff=5.0,
	       datadisplay=False,
	       flagbackup=True,
	       writeflags=True)

# SET THE FLUX DENSITY SCALE OF THE PRIMARY CALIBRATOR
setjy(vis,
	field='2',
	scalebychan=False,
	modimage='3C286_L.im')

# DERIVE PHASE SOLUTION FOR BANDPASS CALIBRATION OF THE PRIMARY CAL
gaincal(vis,
	  caltable='bandpass.pcal',
	  intent='*BANDPASS*',
	  field='2',
	  spw='',
	  solint='int',
	  refant='ea01',
	  minblperant=3,
	  minsnr=3.0,
	  calmode='p',
	  gaintable=[''],
	  gaincurve=T)

# DERIVE AMP SOLUTION FOR BANDPASS CALIBRATION OF THE PRIMARY CAL
bandpass(vis,
	    caltable='bandpass.bcal',
      	    intent='*BANDPASS*',
	    field='2',
	    spw='',
	    solint='inf',
	    combine='scan',
	    refant='ea01,ea05,ea12',
	    minblperant=3,
	    minsnr=10.0,
	    interp=['nearest'],
            gaintable='bandpass.pcal',
	    solnorm=True,
	    gaincurve=T)

# CALCULATE THE PER-ANTENNA GAIN SOLUTIONS
# first to apply to the calibrators
gaincal(vis,
	  caltable='intphase.pcal',
	  field='0,2',
	  spw='0:40~210',
	  solint='int',
	  combine='',
	  refant='ea01,ea05,ea12',
	  minblperant=3,
	  calmode='p',
	  minsnr=3.0,
	  gaintable='bandpass.bcal',
	  gaincurve=T)

# second to apply to the source
gaincal(vis,
	  caltable='scanphase.pcal',
	  field='0,2',
	  spw='0:40~210',
	  solint='inf',
	  combine='combine',
	  refant='ea01,ea05,ea12',
	  minblperant=3,
	  calmode='p',
	  minsnr=3.0,
	  gaintable='bandpass.bcal',
	  gaincurve=T)

# third, write amp gain solutions to calibrators, applying intphase.pcal
gaincal(vis,
	  caltable='amp.cal',
	  field='0,2',
	  spw='0:40~210',
	  solint='inf',
	  refant='ea01,ea05,ea12',
	  calmode='ap',
	  minblperant=3,
	  minsnr=3.0,
	  gaintable=['bandpass.bcal','intphase.pcal'],
	  gaincurve=T)

# SET THE FLUX SCALE FOR THE GAIN SOLUTIONS
fluxscale(vis,
	    caltable='amp.cal',
	    fluxtable='fluxScale.cal',
	    reference='2',
	    transfer='0')
# APPLY THE CALIBRATION
applycal(vis,
	   field='2',
	   gaintable=['bandpass.bcal','fluxScale.cal','intphase.pcal'],
	   calwt=False,
	   gainfield=['2','2','2'],
	   gaincurve=T,
	   interp=['','nearest','nearest'])

applycal(vis,
	   field='0',
	   gaintable=['bandpass.bcal','fluxScale.cal','intphase.pcal'],
	   gainfield=['2','0','0'],
	   calwt=False,
	   gaincurve=T,
	   interp=['','nearest','nearest'])

# poor calibration, flag further
tflagdata(vis,
	  field='',
          spw='',
	  timerange='2012/08/18/17:48:00~2012/08/18/19:24:20',
	  antenna='ea16&ea25',
	  flagbackup=True)

tflagdata(vis,
	  field='',
          spw='',
	  timerange='2012/08/18/17:48:00~2012/08/18/19:24:20',
	  antenna='ea16&ea25',
	  flagbackup=True)

tflagdata(vis,
	  field='',
          spw='',
	  timerange='2012/08/18/17:49:00~2012/08/18/18:35:20',
	  antenna='ea05&ea25',
	  flagbackup=True)

tflagdata(vis,
	  field='',
          spw='',
	  timerange='2012/08/18/20:15:00~2012/08/18/21:03:27.5',
	  antenna='ea05&ea25',
	  flagbackup=True)

tflagdata(vis,
	  field='',
          spw='',
	  timerange='2012/08/18/17:48:00~2012/08/18/18:35:00',
	  antenna='ea01&ea16',
	  flagbackup=True)

tflagdata(vis,
	  field='',
          spw='',
	  timerange='2012/08/18/17:48:00~2012/08/18/18:35:00',
	  antenna='ea01&ea16',
	  flagbackup=True)

tflagdata(vis,
	  field='',
          spw='',
	  timerange='',
	  antenna='ea07',
	  flagbackup=True)

applycal(vis,
	 field='1',
	 gaintable=['bandpass.bcal','fluxScale.cal','scanphase.pcal'],
	 gainfield=['2','0','0'],
	 calwt=False,
	 gaincurve=T,
	 interp=['','nearest','nearest'])

