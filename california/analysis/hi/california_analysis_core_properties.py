#!/usr/bin/python

import numpy as np

def main():

    import json

    output_dir = '/d/bip3/ezbc/california/data/python_output/core_properties/'
    filename = 'california_core_properties.txt'

    # define core properties
    cores = {'L1482':
                {'center_wcs': [(4, 30, 40), (35, 53, 57.60)],
                 'map': None,
                 'threshold': None,
                 },
              'L1482-1':
                {'center_wcs': [(4, 30, 40), (35, 15, 57.60)],
                 'map': None,
                 'threshold': None,
                 },
              'L1482-2':
                {'center_wcs': [(4, 30, 40), (35, 52, 57.60)],
                 'map': None,
                 'threshold': None,
                 },
              'L1483':
                {'center_wcs': [(4, 35, 15), (36, 18, 12)],
                 'map': None,
                 'threshold': None,
                 },
              'L1483-1':
                {'center_wcs': [(4, 35, 15), (36, 21, 12)],
                 'map': None,
                 'threshold': None,
                 },
              'L1483-2':
                {'center_wcs': [(4, 39, 15), (34, 58, 12)],
                 'map': None,
                 'threshold': None,
                 },
              'L1478':
                {'center_wcs': [(4, 25, 7), (37, 9, 0)],
                 'map': None,
                 'threshold': None,
                 },
              'L1478-1':
                {'center_wcs': [(4, 25, 7), (37, 10, 0)],
                 'map': None,
                 'threshold': None,
                 },
              'L1478-2':
                {'center_wcs': [(4, 21, 7), (37, 33, 0)],
                 'map': None,
                 'threshold': None,
                 },
              'L1473':
                {'center_wcs': [(4, 11, 0), (38, 8, 57.8)],
                 'map': None,
                 'threshold': None,
                 },
              'NGC 1579':
                {'center_wcs': [(4, 30, 0), (35, 13, 57.60)],
                 'map': None,
                 'threshold': None,
                 },
              'L1442':
                {'center_wcs': [(3, 51, 55.2), (39, 1, 12)],
                 'map': None,
                 'threshold': None,
                 },
                'L1449':
                {'center_wcs': [(3, 51, 45.6), (38, 16, 12)],
                 'map': None,
                 'threshold': None,
                 },
              'L1456':
                {'center_wcs': [(3, 54, 52.8), (37, 19, 48)],
                 'map': None,
                 'threshold': None,
                 },
              #'':
              # {'center_wcs': [(), ()],
              #  'map': None,
              #  'threshold': None,
              #  },
              #'L1434':
              #  {'center_wcs': [(3, 55, 51), (37, 21, 0)],
              #   'map': None,
              #   'threshold': None,
              #   },
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


