#!/usr/bin/python

import numpy as np

def main():

    import json

    output_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
    filename = 'taurus_core_properties.txt'

    # define core properties
    cores = {'L1495':
                {'center_wcs': [(4,13,40), (28, 15, 0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(4,16,30), (27,44,30), (4,5,20), (28,28,33)]
                 },
             'L1495A':
                {'center_wcs': [(4,18,40), (28, 25., 0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(4,28,23),(28,12,50),(4,16,23),(29,46,5)],
                 },
             'B213':
                {'center_wcs': [(4, 19, 0), (27, 15, 0)],
                 'map': None,
                 'threshold': 4.75,
                 'box_wcs': [(4,22,27), (26,45,47),(4,5,25),(27,18,48)],
                },
             'B220':
                {'center_wcs': [(4, 41, 20.), (25,45,0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,47,49),(25,31,13),(4,40,37),(27,31,17)],
                 },
             'B220-1':
                {'center_wcs': [(4, 41, 20.), (26, 4, 0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,47,49),(25,31,13),(4,40,37),(27,31,17)],
                 },
             'B220-2':
                {'center_wcs': [(4, 41, 20.), (25,40,0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,47,49),(25,31,13),(4,40,37),(27,31,17)],
                 },
             'L1527':
                {'center_wcs': [(4, 39, 20.), (25, 35, 0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,40,13), (24,46,38), (4,34,35), (25,56,7)],
                 },
             'L1527-1':
                {'center_wcs': [(4, 49, 20.), (26, 12, 0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,40,13), (24,46,38), (4,34,35), (25,56,7)],
                 },
             'L1527-2':
                {'center_wcs': [(4, 39, 20.), (25, 49, 0)],
                 'map': None,
                 'threshold': 7,
                 'box_wcs': [(4,40,13), (24,46,38), (4,34,35), (25,56,7)],
                 },
             'B215':
                {'center_wcs': [(4, 23, 20), (25, 5, 0)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,24,51), (22,36,7), (4,20,54), (25,26,31)],
                 },
             'L1524':
                {'center_wcs': [(4,29,0.), (24,40.,0)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,31,0), (22,4,6), (4,25,33), (25,0,55)],
                 },
             'L1498':
                {'center_wcs': [(4,10,40.), (25,6,0)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,31,0), (22,4,6), (4,25,33), (25,0,55)],
                 },
             'B18':
                {'center_wcs': [(4,32,40.), (24,22,30)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,31,0), (22,4,6), (4,25,33), (25,0,55)],
                 },
             'L1521':
                {'center_wcs': [(4,33,25.), (26,10,0)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,31,0), (22,4,6), (4,25,33), (25,0,55)],
                 },
             'L1536':
                {'center_wcs': [(4,33,20.), (22,43,30)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,31,0), (22,4,6), (4,25,33), (25,0,55)],
                 },
             'B217':
                {'center_wcs': [(4,27,0.), (26,6,15)],
                 'map': None,
                 'threshold': 3,
                 'box_wcs': [(4,31,0), (22,4,6), (4,25,33), (25,0,55)],
                 },
                }

    with open(output_dir + filename, 'w') as f:
        json.dump(cores, f)

if __name__ == '__main__':
    main()


