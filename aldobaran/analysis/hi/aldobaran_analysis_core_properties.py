#!/usr/bin/python

import numpy as np

def main():

    import json

    output_dir = '/d/bip3/ezbc/aldobaran/data/python_output/core_properties/'
    filename = 'aldobaran_core_properties.txt'

    # define core properties
    cores = {'L1517':
                {'center_wcs': [(4, 54, 54), (30, 36, 37)],
                 'map': None,
                 'threshold': None,
                 },
              'L1513':
                {'center_wcs': [(4, 52, 3), (30, 53, 14)],
                 'map': None,
                 'threshold': None,
                 },
              'L1512':
                {'center_wcs': [(5, 5, 3), (32, 47, 14)],
                 'map': None,
                 'threshold': None,
                 },
              'L1523':
                {'center_wcs': [(5, 5, 60), (31, 47, 46)],
                 'map': None,
                 'threshold': None,
                 },
              'L1545':
                {'center_wcs': [(5, 15, 60), (26, 31, 17)],
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


