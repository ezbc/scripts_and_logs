#!/usr/bin/python

import numpy as np

def main():

    import json

    output_dir = '/d/bip3/ezbc/california/data/python_output/'
    filename = 'california_global_properties.txt'

    # define core properties
    properties = {'dust2gas_ratio' : {'value' : 5.33e-2,
                                      'unit' : '10^-22 mag / 10^20 cm^-2'
                                      },
                  'dust2gas_ratio_error' : {'value' : 2.22e-2,
                                            'unit' : '10^-22 mag / 10^20 cm^-2'
                                      },
                  'metallicity' : {'value' : 1.00,
                                   'unit' : 'Z_solar',
                                   },
                  'region_limit' : {'wcs' : (((4, 50, 0), (32, 0, 0)),
                                             ((3, 40, 0), (44, 0, 0))),
                                    'pixel' : (),
                                    },
                  'co_noise_limits' : {'wcs' : [(((4, 31, 20), (32, 45, 0)),
                                                 ((4, 23, 20), (34, 40, 0))),
                                                (((4, 34, 20), (40, 15, 0)),
                                                 ((4, 25, 0), (42, 25, 0)))
                                                 ],
                                       'pixel' : (),
                                       }
                  }

    with open(output_dir + filename, 'w') as f:
        json.dump(properties, f)

if __name__ == '__main__':
    main()


