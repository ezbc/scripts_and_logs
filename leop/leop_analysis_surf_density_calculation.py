#!/usr/bin/python

import numpy as np

def flux2Mass(distance=None,arcsecPerPix=None,beamsize=None,
        flux=None,fluxErr=None):
        '''
        Distance in units of Mpc
        flux in units of Jy km/s
        beamsize in units of "
        Returns mass in units of Msun
        '''
        distanceErr = 0.40
        beamarea = 1.13*beamsize**2/arcsecPerPix**2
        mass = 2.36e5 * distance**2 * flux / beamarea
        massErr = np.sqrt((mass * 2*distanceErr/distance)**2 + (fluxErr/flux)**2)
        return mass, massErr


def mass2SB(mass=None,massErr=None,area=None,arcsecPerPix=None,
        distance=None):
        '''
        Distance in units of Mpc
        mass in units of Msun
        beamsize in units of "
        '''
        pcPerArcsec = distance * 1e6 / 206265.
        pcPerPixel = arcsecPerPix * pcPerArcsec
        pixSB = mass / area # units in Msun per pixel
        return pixSB / pcPerPixel**2, massErr/area / pcPerPixel**2


flux2Maxx(distance=1.74e6,arcsecPerPix=1.5,beamsize=4,flux=0.00823,fluxErr=0.1*0.00823)
-------------------------------------------------------------------
r                                 Traceback (most recent call last)
-input-9-50de6ea11929> in <module>()
flux2Maxx(distance=1.74e6,arcsecPerPix=1.5,beamsize=4,flux=0.00823,fluxErr=0.1*0.00823)

r: name 'flux2Maxx' is not defined

 flux2Mass(distance=1.74e6,arcsecPerPix=1.5,beamsize=4,flux=0.00823,fluxErr=0.1*0.00823)
 (731803406415929.4, 336461336.2831859)

 flux2Mass(distance=1.74,arcsecPerPix=1.5,beamsize=4,flux=0.00823,fluxErr=0.1*0.00823)
 (731.8034064159293, 336.46135114373408)

 mass,massErr = flux2Mass(distance=1.74,arcsecPerPix=1.5,beamsize=4,flux=0.00823,fluxErr=0.1*0.00823)

 mass2SB(mass=mass,massErr=massErr,area=1,arcsecPerPix=1.5,distance=1.74)
 (4.570508219414436, 2.1013831821981115)


