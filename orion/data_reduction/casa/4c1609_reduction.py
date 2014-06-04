#!/usr/bin/python

# ==============================================================================
# 0: Prepare MS
# ==============================================================================

import os

working_dir = '/d/bip3/ezbc/orion/data/absorption/casa/'

os.chdir(working_dir)

os.system('tar -xvf /d/bip4/NRAO_data/ezbc/10C-196.sb4997632.eb5088657.55813.373845856484.ms.tar')

os.system('mv lustre/aoc/ftp/e2earchive/10C-196.sb4997632.eb5088657.55813.373845856484.ms* .')

os.system('rm -rf lustre')

# Define visibility name as something easy
my_vis = '10C-196.sb4997632.eb5088657.55813.373845856484.ms'

# ==============================================================================
# 1: Examine observational setup
# ==============================================================================

# List the observation details
listobs(vis = my_vis, listfile = '4c1609_listobs.txt')

# Plot the antenna locations
plotants(vis = my_vis, figfile = '4c1609_antennas_locs.png')

# Three spectral windows
#   0 : source frequency
#   1 : source 2 MHz negative offset, BP calibrator 2 MHz negative offset
#   2 : source 2 MHz positive offset, BP calibrator 2 MHz positive offset

# J0318+1628 = source
#   field 0 --> on source frequency
#   field 1 --> off-source frequencies
# J0137+3309 = BP calibrator
#   field 2

# define calibrators
source_primary = '2'
source_phase = '1'
source_target = '0'

# ==============================================================================
# 2: Flag the data
# ==============================================================================

# Plot the visibilities one correlation at a time, examine for bad visibilities
plotms(vis = my_vis,
        field = source_target,
        xaxis = 'time',
        yaxis = 'amp',
        correlation = 'RR',
        coloraxis = 'antenna1')

plotms(vis = my_vis,
        xaxis = 'time',
        yaxis = 'amp',
        correlation = 'LL',
        coloraxis = 'antenna')

# Flag any bad data. It is in your best interest to script any flags.
flagdata(vis = my_vis,
	 antenna = 'ea28',
	 correlation = '',
	 scan = '20',
     spw = '',
     timerange = '')

flagdata(vis = my_vis,
	 antenna = 'ea06',
	 correlation = '',
	 scan = '6',
     spw = '',
     timerange = '')

flagdata(vis = my_vis,
	 antenna = 'ea06',
	 correlation = '',
	 scan = '',
     spw = '',
     timerange = '')


flagdata(vis = my_vis,
	 antenna = '',
	 correlation = '',
	 scan = '',
     spw = '',
     timerange = '')

# ==============================================================================
# 4: Apply antenna position corrections
# ==============================================================================

gencal(vis = my_vis,
       caltable = 'antenna_positions.gencal',
       caltype = 'antpos')

    # no offsets found, be sure to pay attention to this later

# ==============================================================================
# 5: Acquire phase solutions for the bandpass derivation
# ==============================================================================

# use plotms to plot phase vs time averaged across all channels colorized by
# antennas to determine a stable reference antenna

ref_ant = 'ea28'

# Derive phase solutions for bandpass
gaincal(vis = my_vis,
        caltable = 'bandpass.gcal',
        field = source_primary,
        refant = ref_ant,
	    calmode = 'ap')

# Plot the phase solution, should be just a few points
plotcal(caltable = 'bandpass.gcal',
        xaxis = 'freq',
        yaxis = 'phase',
        iteration = 'antenna',
        subplot = 221)

# ==============================================================================
# 6: Derive the bandpass solution for the calibrators
# ==============================================================================

# Derive bandpass to apply for gain calibrations
bandpass(vis = my_vis,
        caltable = 'bandpass.bcal',
	    field = source_primary,
	    spw = '1, 2',
        refant = ref_ant,
	    solint = 'inf',
	    solnorm = True,
        gaintable = ['bandpass.gcal'])

# ==============================================================================
# 7: Check the bandpass solution
# ==============================================================================

plotcal(caltable = 'bandpass.bcal',
        xaxis = 'freq',
        yaxis = 'amp',
        spw = '1,2',
        iteration = 'antenna',
        subplot = 221)

plotcal(caltable = 'bandpass.bcal',
        xaxis = 'freq',
        yaxis = 'amp',
        spw = '1',
        iteration = 'antenna',
        subplot = 221)

# ==============================================================================
# 8: Derive the bandpass solution for the source
# ==============================================================================

# The following link shows how to interpolate the bandpasses:
# http://casaguides.nrao.edu/index.php?title=Combining_Bandpasses

# Derive the bandpass, interpolate across spw 1 and 2 to acquire a bandpass
# solution for 0. Do this by setting combine = 'spw, scans'
bandpass(vis = my_vis,
        caltable = 'bandpass_source.bcal',
	    field = source_primary,
	    spw = '1, 2',
        refant = ref_ant,
	    solint = 'inf',
	    combine = 'spw, scan',
	    solnorm = True,
        gaintable = ['bandpass.gcal'])

# Note that the combined bandpass of spw 1 and 2 will, by convention, inherit
# the frequency label of the first spw. The solution will be under spw 1.


# ==============================================================================
# 9: Check the bandpass solution
# ==============================================================================

# Plot the bandpass solution
# Blue and green points correspond to the RR and LL polarizations

plotcal(caltable = 'bandpass_source.bcal',
        xaxis = 'freq',
        yaxis = 'amp',
        spw = '1',
        iteration = 'antenna',
        subplot = 221)

# ==============================================================================
# 10: Derive phase and amplitude calibrations for each calibrator
# ==============================================================================

# Phases
gaincal(vis = my_vis,
        caltable = 'scanphase.gcal',
        field = '1, 2',
        spw = '1, 2',
        refant = ref_ant,
        calmode = 'p',
        solint = 'inf',
        spwmap = 1,
        gaintable = ['bandpass.bcal'])

# Amplitudes
gaincal(vis = my_vis,
        caltable = 'amp.gcal',
        field = '1, 2',
        spw = '1, 2',
        refant = ref_ant,
        calmode = 'ap',
        solint = 'inf',
        gaintable = ['bandpass.bcal', 'scanphase.gcal'])

# ==============================================================================
# 11: Derive phase and amplitude calibrations for the source
# ==============================================================================

gaincal(vis = my_vis,
        caltable = 'amp_source.gcal',
        field = '1, 2',
        combine = 'spw',
        refant = ref_ant,
        calmode = 'ap',
        solint = 'inf',
        gaintable = ['bandpass.bcal', 'scanphase.gcal'])

# ==============================================================================
# 12: Check the gain solutions
# ==============================================================================

# Plot the bandpass solution
# Blue and green points correspond to the RR and LL polarizations

plotcal(caltable = 'amp.gcal',
        xaxis = 'time',
        yaxis = 'amp',
        spw = '1',
        field = source_phase,
        iteration = 'antenna',
        subplot = 221)

plotcal(caltable = 'amp.gcal',
        xaxis = 'time',
        yaxis = 'amp',
        spw = '2',
        field = source_phase,
        iteration = 'antenna',
        subplot = 221)

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
        field = source_primary,
        iteration = 'antenna',
        subplot = 221)

# ==============================================================================
# 13: Flag visibilities based on phase and bandpass solutions
# ==============================================================================

flagdata(vis = my_vis,
	 antenna = '',
	 correlation = '',
	 scan = '',
     spw = '',
     timerange = '')

# ==============================================================================
# 14: Set the flux scale
# ==============================================================================

setjy(vis = my_vis,
      field = source_primary,
      modimage = '3C48_L.im')

    # flux measured as 16.034 Jy, well-matched to 16.50 Jy from calibrator
    # manual

# ==============================================================================
# 15: Apply the flux scaling
# ==============================================================================

fluxscale(vis = my_vis,
          caltable = 'amp.gcal',
          fluxtable = 'flux.cal',
          reference = source_primary)

fluxscale(vis = my_vis,
          caltable = 'amp_source.gcal',
          fluxtable = 'flux_source.cal',
          reference = source_primary,
          refspwmap = [1, 1, 2])

    # VLA calibrator manual reports J0138+1629 20cm flux as 7.81 Jy
    # Calculated 7.57 Jy

# ==============================================================================
# 16: Apply the calibrations
# ==============================================================================

# Apply to primary calibrator
applycal(vis = my_vis,
         field = source_primary,
         gaintable = ['bandpass.bcal', 'scanphase.gcal', 'flux.cal'],
         gainfield = [source_primary,source_primary,source_primary])

# Apply to phase calibrator
applycal(vis = my_vis,
         field = source_phase,
         gaintable = ['bandpass.bcal', 'scanphase.gcal', 'flux.cal'],
         gainfield = [source_primary,source_phase,source_phase])

# Apply to source
applycal(vis = my_vis,
        field = source_target,
        spw = '0',
        spwmap = [[1], [1], [1]],
        gainfield = [source_primary, source_phase, source_phase],
        gaintable = ['bandpass_source.bcal', 'amp_source.gcal', 'flux.cal'])

# Apply to source
applycal(vis = my_vis,
        field = source_target,
        spwmap = [[1], [1], [1]],
        gainfield = [source_primary, source_phase, source_phase],
        gaintable = ['bandpass_source.bcal', 'amp_source.gcal',
            'flux_source.cal'])

# Apply to source without flux scaling
applycal(vis = my_vis,
        field = source_target,
        spwmap = [[1], [1]],
        gainfield = [[source_primary], [source_primary]],
        gaintable = ['bandpass.bcal', 'amp.gcal'],
        flagbackup = False)

# ==============================================================================
# 17: Split the source
# ==============================================================================

split_vis = '4c1609_source.ms'

if os.path.isdir(split_vis): rmtables(split_vis)

split(vis = my_vis,
      field = source_target,
      outputvis = split_vis,
      datacolumn = 'corrected')





