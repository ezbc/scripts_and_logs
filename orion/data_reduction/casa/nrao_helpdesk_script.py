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
        refant = 'ea28',
	    solint = 'inf',
	    combine = 'scan',
        gaintable = [''])

################################################################################
# 3: Derive the bandpass solution for the source
################################################################################

bandpass(vis = my_vis,
        caltable = 'bandpass_source.bcal',
	    field = '2',
	    spw = '1, 2',
        refant = 'ea28',
	    solint = 'inf',
	    combine = 'spw, scan',
        gaintable = [''])

################################################################################
# 4: Derive phase and amplitude calibrations for each calibrator
################################################################################

# Phases
gaincal(vis = my_vis,
        caltable = 'scanphase.gcal',
        field = '1, 2',
        spw = '1, 2',
        refant = 'ea28',
        calmode = 'p',
        solint = 'inf',
        gaintable = ['bandpass.bcal'])

# Amplitudes + phases
gaincal(vis = my_vis,
        caltable = 'amp.gcal',
        field = '1, 2',
        refant = 'ea28',
        calmode = 'ap',
        solint = 'inf',
        gaintable = ['bandpass.bcal', 'scanphase.gcal'])

################################################################################
# 5: Derive phase and amplitude calibrations for the source
################################################################################

gaincal(vis=my_vis,
        caltable='amp_source.gcal',
        field='0, 1',
        combine='spw, field',
        refant='ea28',
        calmode='ap',
        solint='inf',
        spwmap=[[1, 1, 2], [1, 1, 2]],
        gaintable=['bandpass_source.bcal', 'scanphase.gcal'])

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
          transfer = '0',
          refspwmap = [0])

    # VLA calibrator manual reports J0138+1629 20cm flux as 7.81 Jy
    # Calculated 7.59 Jy for field 1

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
         gaintable = ['bandpass.bcal', 'scanphase.gcal', 'flux.cal'],
         gainfield = ['2','1','1'])

# Apply to source
applycal(vis = my_vis,
        field = '0',
        spw = '0',
        spwmap = [[1], [0], [0]],
        gainfield = ['2', '0', '0'],
        gaintable = ['bandpass_source.bcal', 'amp_source.gcal',
            'flux_source.cal'])


