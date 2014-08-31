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
                                   'unit' : 'Z_solar'},
                  }

    with open(output_dir + filename, 'w') as f:
        json.dump(properties, f)

if __name__ == '__main__':
    main()


