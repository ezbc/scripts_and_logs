#!/usr/bin/python

# Observational set-up
# --------------------
# Three spectral windows
#   0 : source frequency
#   1 : source 2 MHz negative offset, BP calibrator 2 MHz negative offset
#   2 : source 2 MHz positive offset, BP calibrator 2 MHz positive offset

# Two sources, three fields
# J0318+1628 = source
#   field 0 --> on source frequency
#   field 1 --> off-source frequencies
# J0137+3309 = BP calibrator
#   field 2 --> off-source frequencies

my_vis = '10C-196.sb4997632.eb5088657.55813.373845856484.ms'

################################################################################
# 1: Set the flux scale
################################################################################

setjy(vis = my_vis,
      field = '2',
      modimage = '3C48_L.im')

################################################################################
# 2: Derive the bandpass solution for the calibrators
################################################################################

# Derive bandpass to apply for gain calibrations
bandpass(vis = my_vis,
        caltable = 'bandpass.bcal',
	    field = '2',
	    spw = '1, 2',
        refant = 'ea05',
	    solint = 'inf',
	    combine = 'scan',
        gaintable = [''])

################################################################################
# 3: Derive the bandpass solution for the source
################################################################################

gaincal(vis = my_vis,
        caltable = 'bpphase.gcal',
        field = '2',
        spw = '1, 2',
        refant = 'ea05',
        calmode = 'p',
        solint = 'int',
        gaintable = [''])

bandpass(vis = my_vis,
        caltable = 'bandpass_source.bcal',
	    field = '2',
	    spw = '1, 2',
	    solint = 'inf',
	    combine = 'spw, scan',
        refant = 'ea05',
	    solnorm=True,
	    interp='linear',
        gaintable = ['bpphase.gcal'])

bandpass(vis = my_vis,
        caltable = 'bandpass_source.bcal',
	    field = '2',
	    spw = '1, 2',
	    solint = 'inf',
	    combine = 'spw, scan',
        refant = 'ea05',
	    solnorm=True,
	    interp='linear',
        gaintable = ['bpphase.gcal'])

bandpass(vis = my_vis,
        caltable = 'bandpass_source.bcal',
	    field = '2',
        spw = '1, 2',
	    solint = 'inf',
	    combine = 'spw, scan',
        refant = 'ea05',
	    solnorm=True,
	    interp='linear',
        gaintable = ['amp_source.gcal'])

################################################################################
# 4: Derive phase and amplitude calibrations for each calibrator
################################################################################

# Phases
gaincal(vis = my_vis,
        caltable = 'scanphase.gcal',
        field = '1, 2',
        spw = '1, 2',
        refant = 'ea05',
        calmode = 'p',
        solint = 'inf',
        gaintable = ['bandpass.bcal'])

# Amplitudes + phases
gaincal(vis = my_vis,
        caltable = 'amp.gcal',
        field = '1, 2',
        refant = 'ea05',
        calmode = 'ap',
        solint = 'inf',
        gaintable = ['bandpass.bcal', 'scanphase.gcal'])

################################################################################
# 5: Derive phase and amplitude calibrations for the source
################################################################################

# Phases
gaincal(vis = my_vis,
        caltable = 'scanphase_source.gcal',
        field = '0, 1',
        spw = '0, 1',
        combine = 'spw, field',
        refant = 'ea05',
        calmode = 'p',
        solint = 'inf',
        spwmap=[[1, 1]],
        gaintable = ['bandpass_source.bcal'])

gaincal(vis=my_vis,
        caltable='amp_source.gcal',
        field='0, 1',
        spw='0, 1',
        combine='spw, field',
        refant='ea05',
        calmode='ap',
        solint='inf',
        spwmap=[[1, 1,], [1, 1]],
        gaintable=['bandpass_source.bcal', 'scanphase_source.gcal'])

gaincal(vis = my_vis,
        caltable = 'amp_source.gcal',
        append=True,
        field = '2',
        refant = 'ea05',
        calmode = 'ap',
        solint = 'inf',
        gaintable = ['bandpass.bcal', 'scanphase.gcal'])

gaincal(vis = my_vis,
        caltable = 'amp_source.gcal',
        field = '0, 1, 2',
        refant = 'ea05',
        calmode = 'ap',
        solint = 'inf',
        )

################################################################################
# 6: Apply the flux scaling
################################################################################

fluxscale(vis = my_vis,
          caltable = 'amp.gcal',
          fluxtable = 'flux.cal',
          reference = '2')

fluxscale(vis = my_vis,
          caltable = 'amp_source.gcal',
          fluxtable = 'flux_source.cal',
          reference = '2',
          refspwmap = [0, 1, 2])

    # VLA calibrator manual reports J0138+1629 20cm flux as 7.81 Jy
    # Calculated 7.57 Jy

################################################################################
# 7: Apply the calibrations
################################################################################

# Apply to primary calibrator
applycal(vis = my_vis,
         field = '2',
         gaintable = ['bandpass.bcal', 'scanphase.gcal', 'flux.cal'],
         gainfield = ['2','2','2'])

# Apply to phase calibrator
applycal(vis = my_vis,
         field = '1',
         gaintable = ['bandpass_source.bcal',],
         spwmap=['1'],
         gainfield = ['2',])

# Apply to source
applycal(vis = my_vis,
        field = '0',
        spw = '0',
        spwmap = [[1], [0], [0]],
        gainfield = ['2', '0', '1'],
        gaintable = ['bandpass_source.bcal', 'amp_source.gcal',
            'flux_source.cal'],
        flagbackup = False)

applycal(vis = my_vis,
        field = '0',
        spw = '0',
        spwmap = [[1],],
        gainfield = ['2',],
        gaintable = ['bandpass_source.bcal',],
        flagbackup = False)

applycal(vis = my_vis,
        field = '0',
        spw = '0',
        gaintable = ['amp_source.gcal', 'bandpass_source.bcal'],
        gainfield = ['1', '2'],
        interp='linear',
        spwmap = [[], [1]])

applycal(vis = my_vis,
        field = '0',
        spw = '0',
        gaintable = ['amp_source.gcal','bandpass_source.bcal'],
        gainfield = ['1', '2'],
        interp='linear',
        spwmap = [[1], [1]])

###########################
# Plotting

# define calibrators
source_primary = '2'
source_phase = '1'
source_target = '0'

plotms(vis=my_vis,
       xaxis='freq',
       yaxis='amp',
       ydatacolumn='corrected',
       field='0',
       scan='0~10',
       avgtime='1e8',
       avgbaseline=True)

plotcal(caltable = 'amp_source.gcal',
        xaxis = 'time',
        yaxis = 'phase',
        spw = '1',
        field = source_primary,
        iteration = 'antenna',
        subplot = 221)

plotcal(caltable = 'amp_source.gcal',
        xaxis = 'time',
        yaxis = 'amp',
        spw = '1',
        field = '1',
        iteration = 'antenna',
        subplot = 221)

plotcal(caltable = 'amp_source.gcal',
        xaxis = 'time',
        yaxis = 'amp',
        spw = '0',
        field = '0',
        iteration = 'antenna',
        subplot = 221)


plotcal(caltable = 'bandpass_source.bcal',
        xaxis = 'freq',
        yaxis = 'amp',
        spw = '1',
        field = '2',
        iteration = 'antenna',
        subplot = 221)

applycal(vis = my_vis,
        field = '0',
        spw = '0',
        spwmap = [[1], [0], [0]],
        gainfield = ['2', '0', '0'],
        gaintable = ['bandpass_source.bcal', 'amp_source.gcal',
            'flux_source.cal'])






