#!/usr/bin/python

import numpy as np

def transform_vel(v_hsr, l, b, vel_type='lsr'):
    from numpy import cos, sin

    apex_dict = {'gsr': {'l_apex': 87.8,
                         'b_apex': 1.7,
                         'v_apex': 232.3,
                         'reference': ''},
                 'lsr': {'l_apex': 30,
                         'b_apex': 23,
                         'v_apex': -19.8,
                         'reference': ''},
                 'lg':  {'l_apex': 93,
                         'b_apex': -4,
                         'v_apex': 316,
                         'reference': 'AJ 111, 794, 1996'}
                 }

    l_apex = np.radians(apex_dict[vel_type]['l_apex'])
    b_apex = np.radians(apex_dict[vel_type]['b_apex'])
    v_apex = apex_dict[vel_type]['v_apex']

    return v_hsr + v_apex*(sin(b)*sin(b_apex) + \
                           cos(b)*cos(b_apex)*cos(l-l_apex))

def main():

    '''

    Most of the velocity measurements are published in heliocentric (or
    barycentric) reference frame. You find the V=260.8+-2.5 in  the Local
    Standard of Rest Kinematic (LSRK) frame. I suppose the transformation from
    Vh to LSRK is
    -20.0*(sin(DEC.obj))*sin(DEC.apex)+cos(DEC.obj)*cos(DEC.apex)*cos(RA.obj-RA.apex)
    where RA.obj, DEC.obj are coordinates of an object on the sky in radians
    and RA.apex,DEC.apex are coordinates of the sun apex respect the LSRK. So,
    I use Sun motion respect to LSRK: V.apex=-20,
    RA.apex=18.063955555555555h=18 03 50.24 DEC.apex=+30.004666666666665deg=+30
    00 16.8 for J2000 equinox.

    In case Leo P I received dV=3.988 km/s. If I am right I have to remove that
    correction from V_LSRK to get heliocentric velocity. So, Vh=256.8. From
    other hand, Giovanelli et al. 2013 published Vh=264+-2 km/s. Obviously, the
    problem is in the sign of the correction.

    converting dmitry's ra and dec into galactic
    l = 56.31196923, b = 22.3745517


    NED converts velocities from one reference frame to another using the
    standard equation



    where l and b are the object's longitude and latitude, V is its unconverted
    velocity, and the apices (with Galactic coordinates) of the various motions
    are given as

    Conversion	            lapex	    bapex	    Vapex	        Source
    Helio to Galactic (GSR)	87.8 deg	+1.7 deg	232.3 km/sec	RC3
    Helio to Local Group	93 deg	    -4 deg	    316 km/sec	    AJ 111, 794, 1996
    Helio to 3K Background	264.14 deg	+48.26 deg	371.0 km/sec	ApJ 473, 576, 1996
    Helio to LSR            30          23          19.8 km/sec


    Courteau and van den Bergh AJ 118, 337, 1999 give reference frame
    definitions
    http://iopscience.iop.org/1538-3881/118/1/337/pdf/990091.web.pdf



    response to dmitry:

    From the following reference: Courteau and van den Bergh, 1999, AJ, 118,
    337 http://adsabs.harvard.edu/cgi-bin/bib_query?1999AJ....118..337C

    equation 4 gives:
    v_lsr = v_helio + 9*cos(l)*cos(b) + 12*sin(l)*cos(b) + 7*sin(b)

    where
    v_lsr = local standard of rest velocity
    v_helio = heliocentric velocity
    l = galactic longitude
    b = galactic latitude

    Using the position and heliocentric velocity of Leo P
    l = 219.654 degrees
    b = 54.430 degrees
    v_helio = 264

    and converting into radians I find a dV of -2.79 km/s, revealing v_lsr =
    261.21 km/s. This is within the published errors, however I suspect the
    difference of the calculation in the publication arises from the
    inhomogeneous data. The doppler tracking of the Earth's motion from the
    GMRT and the VLA may have been different, thus the radial velocities could
    differ between the VLA observations in Giovanelli et al. (2013) and my
    paper.
    '''


    from numpy import cos, sin

    l = 219.654 # degrees
    b = 54.430 # degrees

    l = np.radians(l)
    b = np.radians(b)

    v_lsr = v_helio + 9*cos(l)*cos(b) + 12*sin(l)*cos(b) + 7*sin(b)
    #v_lsr = transform_vel(v_helio, l, b, vel_type='lsr')

    dV = v_lsr - v_helio

    print 'v_lsr = %.2f' % v_lsr
    print 'v_helio = %.2f' % v_helio
    print 'dV = %.2f' % dV

    # an attempt at a function
    v_gsr = v_helio + 9*cos(l)*cos(b) + 232*sin(l)*cos(b) + 7*sin(b)
    v_gsr = transform_vel(v_helio, l, b, vel_type='gsr')

    #print 'v_gsr = %.2f' % v_gsr

if __name__ == '__main__':
    main()

