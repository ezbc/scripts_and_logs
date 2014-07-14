#!/usr/bin/python

import numpy as np

def main():

    import json

    output_dir = '/d/bip3/ezbc/perseus/data/python_output/core_properties/'
    filename = 'perseus_core_properties.txt'

    # define core properties
    cores = {'IC348':
                {'center_wcs': [(3, 44, 0), (32, 8, 0)],
                 'map': None,
                 'threshold': None,
                 'box_wcs': [(3,46,13), (26,3,24), (3,43,4), (32,25,41)],
                 },
             'NGC1333':
                {'center_wcs': [(3, 29, 11), (31, 16, 53)],
                 'map': None,
                 'threshold': None,
                 'box_wcs': None,
                 },
             #'B4':
             #   {'center_wcs': [(3, 44, 18), (32, 05, 20)],
             #    'map': None,
             #    'threshold': None,
             #    'box_wcs': None,
             #    },
             'B5':
                {'center_wcs': [(3, 47, 34), (32, 48, 17)],
                 'map': None,
                 'threshold': None,
                 'box_wcs': None,
                 },
             'B1':
                {'center_wcs': [(3, 36, 27), (31, 07, 18)],
                 'map': None,
                 'threshold': None,
                 'box_wcs': None,
                 },
             'B1E':
                {'center_wcs': [(3, 33, 31), (31, 10, 31)],
                 'map': None,
                 'threshold': None,
                 'box_wcs': None,
                 },
             #'':
             #   {'center_wcs': [],
             #    'map': None,
             #    'threshold': None,
             #    'box_wcs': None,
             #    },
            }


    with open(output_dir + filename, 'w') as f:
        json.dump(cores, f)

if __name__ == '__main__':
    main()


