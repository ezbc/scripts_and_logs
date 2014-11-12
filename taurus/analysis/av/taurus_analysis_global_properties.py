#!/usr/bin/python

import numpy as np

def main():

    import json

    output_dir = '/d/bip3/ezbc/taurus/data/python_output/'
    filename = 'taurus_global_properties.txt'

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
                  'region_limit' : {'wcs' : (((5, 10,  0), (19, 0, 0)),
                                             ((3, 50, 0), (32, 0, 0))),
                                    'pixel' : ()
                                    },
                  'plot_limit' : {'wcs' : (((5, 10,  0), (19, 0, 0)),
                                             ((3, 50, 0), (32, 0, 0))),
                                    'pixel' : ()
                                    },
                  'co_noise_limits' : {'wcs' : [(((4, 8, 0), (18, 15, 0)),
                                                ((3, 49, 0), (22, 10, 0))),
                                                ],
                                       'pixel' : ()
                                       }
                  }

    with open(output_dir + filename, 'w') as f:
        json.dump(properties, f)

if __name__ == '__main__':
    main()


