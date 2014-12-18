#!/usr/bin/python

import numpy as np

def main():

    import json

    output_dir = '/d/bip3/ezbc/perseus/data/python_output/'
    filename = 'perseus_global_properties.txt'

    # define core properties
    properties = {'dust2gas_ratio' : {'value' : 1.1e-1,
                                      'unit' : '10^-22 mag / 10^20 cm^-2'
                                      },
                  'dust2gas_ratio_error' : {'value' : 2.22e-2,
                                            'unit' : '10^-22 mag / 10^20 cm^-2'
                                      },
                  'metallicity' : {'value' : 1.00,
                                   'unit' : 'Z_solar',
                                   },
                  'region_limit' : {'wcs' : (((3, 58, 0), (26, 6, 0)),
                                             ((3, 0, 0), (35, 0, 0))),
                                    'pixel' : ()
                                    },
                  'plot_limit' : {'wcs' : (((3, 58, 0), (26, 6, 0)),
                                             ((3, 0, 0), (35, 0, 0))),
                                    'pixel' : ()
                                    },
                  'co_noise_limits' : {'wcs' : [(((3, 54, 0), (26, 13, 0)),
                                                  (3, 28, 20, 28, 48, 0)),
                                                  ],
                                       'pixel' : ()},
                  }

    with open(output_dir + filename, 'w') as f:
        json.dump(properties, f)

if __name__ == '__main__':
    main()


