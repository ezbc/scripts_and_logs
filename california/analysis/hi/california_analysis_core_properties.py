#!/usr/bin/python

import numpy as np

def main():

    import json

    output_dir = '/d/bip3/ezbc/california/data/python_output/core_properties/'
    filename = 'california_core_properties.txt'

    # define core properties
    cores = {'L1482':
                {'center_wcs': [(4, 30, 41), (35, 14, 41)],
                 'map': None,
                 'threshold': None,
                 },
              'L1483':
                {'center_wcs': [(4, 35, 15), (36, 18, 12)],
                 'map': None,
                 'threshold': None,
                 },
              'L1478':
                {'center_wcs': [(4, 25, 7), (37, 9, 0)],
                 'map': None,
                 'threshold': None,
                 },
              'L1434':
                {'center_wcs': [(3, 55, 51), (37, 21, 0)],
                 'map': None,
                 'threshold': None,
                 },
              #'L1513':
              #  {'center_wcs': [(4, 51, 27), (30, 49, 12)],
              #   'map': None,
              #   'threshold': None,
              #   },
              #'L1517':
              #  {'center_wcs': [(4, 55, 51), (30, 41, 47)],
              #   'map': None,
              #   'threshold': None,
              #   },
                }


    with open(output_dir + filename, 'w') as f:
        json.dump(cores, f)

if __name__ == '__main__':
    main()


