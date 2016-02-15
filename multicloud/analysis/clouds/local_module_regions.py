#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import myimage_analysis as myia
from multiprocessing.queues import Queue
import mygeometry as myg
import scipy
from local_module_plotting import *
from local_module_dataprep import *
from local_module_multiprocessing import *
from local_module_regions import *
from local_module_fitting import *
from local_module_bootstrapping import *

def read_ds9_region(filename):

    ''' Converts DS9 region file into format for plotting region.

    Need the following format:
        angle : degrees
        xy : pixels
        width : pixels
        height : pixels

    Region file provides following format:
        # Region file format: DS9 version 4.1
        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        fk5
        box(4:17:04.740,+29:20:31.32,5854.33",11972.7",130) # text={test}

    pyregion module reads DS9 regions:
    http://leejjoon.github.io/pyregion/users/overview.html


    '''

    # Import external modules
    import pyregion as pyr

    # Read region file
    try:
        region = pyr.open(filename)
    except IOError:
        return None

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    return region

def load_ds9_core_region(cores, filename='',
        header=None):

    from myimage_analysis import get_pix_coords

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    regions = read_ds9_region(filename)

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        core = tag[tag.find('{')+1:tag.find('}')]

        if core in cores:
            # Format vertices to be 2 x N array
            poly_verts = []
            for i in xrange(0, len(region.coord_list)/2):
                poly_verts.append((region.coord_list[2*i],
                                   region.coord_list[2*i+1]))

            poly_verts_pix = []
            for i in xrange(0, len(poly_verts)):
                poly_verts_pix.append(get_pix_coords(ra=poly_verts[i][0],
                                          dec=poly_verts[i][1],
                                          header=header)[:-1][::-1].tolist())

            cores[core]['poly_verts'] = {}
            cores[core]['poly_verts']['wcs'] = poly_verts
            cores[core]['poly_verts']['pixel'] = poly_verts_pix

    return cores

def load_wedge_core_region(cores, filename='', header=None):

    from myimage_analysis import get_pix_coords
    import pickle
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from astropy.io import fits
    from astropy.wcs import WCS

    # Create WCS object
    wcs_header = WCS(header)

    with open(filename, 'rb') as f:
        region_dict = pickle.load(f)

    for core_name in region_dict:
        if core_name in cores:
            core = region_dict[core_name]
            # Format vertices to be 2 x N array
            poly_verts = np.array((core['ra'], core['dec']))

            # Make a galactic coords object and convert to Ra/dec
            coords_fk5 = SkyCoord(core['ra'] * u.deg,
                                  core['dec'] * u.deg,
                                  frame='fk5',
                                  )
            # convert to pixel
            coords_pixel = np.array(coords_fk5.to_pixel(wcs_header))

            #print coords_pixel
            #print coords_pixel.shape

            # write data to dataframe
            poly_verts_pix = np.array((coords_pixel[1], coords_pixel[0])).T

            #print poly_verts_pix.shape

            #poly_verts_pix = np.array((core['ypix'], core['xpix'])).T

            cores[core_name]['poly_verts'] = {}
            cores[core_name]['poly_verts']['wcs'] = poly_verts
            cores[core_name]['poly_verts']['pixel'] = poly_verts_pix

    return cores

def get_cores_to_plot():

    '''

    '''

    # Which cores to include in analysis?
    cores_to_keep = [ # N(H2) cores
                     'G168.54-6.22',
                     'G168.12-6.42',
                     'G166.91-7.76',
                     'G165.36-7.51',
                     'G164.70-7.63',
                     'G164.65-8.12',
                     'G165.71-9.15',
                     'G164.99-8.60',
                     'G164.26-8.39',
                     'G164.18-8.84',
                     'G174.40-13.45',
                     'G174.70-15.47',
                     'G174.05-15.82',
                     'G171.49-14.91',
                     'G172.93-16.73',
                     'G171.00-15.80',
                     'G172.12-16.94',
                     'G171.14-17.57',
                     'G169.32-16.17',
                     'G168.10-16.38',
                     'G160.49-16.81',
                     'G160.46-17.99',
                     'G159.80-18.49',
                     'G160.14-19.08',
                     'G160.53-19.73',
                     'G159.19-20.11',
                     'G159.17-21.09',
                     'G158.39-20.72',
                     'G158.89-21.60',
                     'G158.26-21.81',
                     ]

    if 0: # Random cores
        cores_to_keep = [
                 'G166.83-8.68',
                 'G168.82-6.37',
                 'G168.05-7.01',
                 'G164.16-8.46',
                 'G165.23-8.78',
                 'G167.06-7.77',
                 'G168.12-6.42',
                 'G167.58-6.64',
                 'G164.70-7.63',
                 'G166.35-8.77',
                 'G166.73-15.06',
                 'G173.08-16.50',
                 'G172.74-14.53',
                 'G169.44-16.18',
                 'G173.86-17.65',
                 'G173.71-13.91',
                 'G171.75-14.18',
                 'G173.70-15.21',
                 'G170.28-19.48',
                 'G171.00-15.80',
                 'G158.23-20.15',
                 'G159.01-22.19',
                 'G159.19-20.11',
                 'G157.12-23.49',
                 'G160.10-19.90',
                 'G160.34-18.42',
                 'G158.40-21.86',
                 'G159.79-21.32',
                 'G158.89-21.60',
                 'G159.51-18.41',
                ]

    if 0:
        cores_to_keep = [# taur
                         'L1495',
                         'L1495A',
                         'B213',
                         'L1498',
                         'B215',
                         'B18',
                         'B217',
                         'B220-1',
                         'B220-2',
                         'L1521',
                         'L1524',
                         'L1527-1',
                         'L1527-2',
                         # Calif
                         'L1536',
                         'L1483-1',
                         'L1483-2',
                         'L1482-1',
                         'L1482-2',
                         'L1478-1',
                         'L1478-2',
                         'L1456',
                         'NGC1579',
                         #'L1545',
                         #'L1517',
                         #'L1512',
                         #'L1523',
                         #'L1512',
                         # Pers
                         'B5',
                         'IC348',
                         'B1E',
                         'B1',
                         'NGC1333',
                         'B4',
                         'B3',
                         'L1455',
                         'L1448',
                         ]

    return cores_to_keep

def get_core_properties(data_dict, cloud_name):

    from myimage_analysis import load_ds9_region, get_pix_coords
    import json

    box_method = 'ds9'
    core_dir = '/d/bip3/ezbc/multicloud/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'

    header = data_dict['av_header']
    if 0:
        # define core properties
        with open(core_dir + cloud_name + '_core_properties.txt', 'r') as f:
            cores = json.load(f)

        #cores = convert_core_coordinates(cores, header)

        # convert center WCS coordinate to pixel
        for core in cores:
            cores[core].update({'box_pixel': 0})
            cores[core].update({'center_pixel': 0})

            center_wcs = cores[core]['center_wcs']

            # convert centers to pixel coords
            center_pixel = get_pix_coords(ra=center_wcs[0],
                                          dec=center_wcs[1],
                                          header=header)[:2]
            cores[core]['center_pixel'] = center_pixel

    else:
        filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
            'planck11_coldclumps.pickle'

        import pickle
        with open(filename, 'rb') as f:
            core_sample = pickle.load(f)
        #core_sample = pd.load(filename)
        cores = {}
        df_cores = core_sample[cloud_name]
        for name in df_cores['Name'].values:
            cores[name] = {}
            ra = df_cores[df_cores['Name'] == name]['ra'].values[0]
            dec = df_cores[df_cores['Name'] == name]['dec'].values[0]
            cores[name]['center_wcs'] = (ra, dec)
            cores[name]['temp'] = \
                df_cores[df_cores['Name'] == name]['temp'].values[0]
            cores[name]['temp_error'] = \
                df_cores[df_cores['Name'] == name]['temp_error'].values[0]

            #print cores[name]['temp_error']

            # convert centers to pixel coords
            center_pixel = get_pix_coords(ra=ra,
                                          dec=dec,
                                          header=header)[:2]
            cores[name]['center_pixel'] = center_pixel[::-1]

    # load the bounding regions
    if 0:
        region_filename = region_dir + 'multicloud_coldclump_divisions.reg'
        cores = load_ds9_core_region(cores,
                                filename=region_filename,
                                header=header)
    else:
        filename = region_dir + 'multicloud_divisions_coldcore_wedges.pickle'
        cores = load_wedge_core_region(cores,
                                       filename=filename,
                                       header=header)

    # add indices of core in data
    add_core_mask(cores, data_dict['av_data'])

    return cores

def trim_cores_to_plot(cores, cores_to_plot):

    # Trim down cores to keep list to include only cores for cloud
    cores_to_keep_old = list(cores_to_plot)
    for core in cores_to_keep_old:
        if core not in cores or 'poly_verts' not in cores[core]:
            cores_to_plot.remove(core)

    return cores_to_plot

def add_core_mask(cores, data):

    for core in cores:
        try:
            vertices = cores[core]['poly_verts']['pixel']

            mask = np.logical_not(myg.get_polygon_mask(data,
                                                       vertices))

            cores[core]['indices_orig'] = np.where(mask == 0)

            cores[core]['mask'] = mask

            #print 'npix in core ' + core + ':'
            #print np.where(mask == 0)[0].size

        except KeyError:
            cores[core]['mask'] = None
            cores[core]['indices_orig'] = None

