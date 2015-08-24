#!/usr/bin/python

import os
import numpy as np

def rewrite_param_file(input_param_filename, output_param_filename):

    paramData = np.genfromtxt(input_param_filename,
                              delimiter='',
                              usecols=(1,2,3,4,5,6))
    modelNameData = np.genfromtxt(input_param_filename,
                                  delimiter='',
                                  usecols=(0),
                                  dtype=str)

    # Data repeats each line, choose only unique lines
    params = paramData[0::2]
    modelNames = modelNameData[0::2]

    param_names = ('modelnum', 'inc', 'pa','iflat','vflat','vdisp', 'z0')

    # Create a text file for writing out the best model data parameters
    file = open(output_param_filename,'wb')

    # Write column headers
    file.writelines(['model #  ', 'inc \t',
                     'pa \t','iflat \t','vflat \t','vdisp \t','z0 \n'])

    # modelNum column dat
    for i in xrange(0, params.shape[0]):
        file.writelines([str(modelNames[i]) + '\t',
                         str(params[i,0]) + '\t',
                         str(params[i,1]) + '\t',
                         str(params[i,2]) + '\t',
                         str(params[i,3]) + '\t',
                         str(params[i,4]) + '\t',
                         str(params[i,5]) + '\n'])

    file.close()

    return params, modelNames, param_names

def main():

    import mystats
    import myio
    import pickle

    os.chdir('/d/bip3/ezbc/shield/749237_lowres/modeling_fineres/')
    os.system('rm -rf *.image *.descr')

    input_param_filename = 'parameters.txt'
    output_param_filename = 'parameters.tsv'

    with open('statistics.pickle', 'rb') as f:
        stats = pickle.load(f)

    if 0:
        import matplotlib.pyplot as plt
        model_names_int = [int(model_name.replace('model_', '')) for model_name in model_names ]
        plt.plot(model_names_int)
        plt.plot(stats['logL'])
        plt.savefig('likelihood_test.png')


    stats['likelihoods'] = mystats.logL2L(stats['logL'])

    params, model_names, param_names = rewrite_param_file(input_param_filename,
                                                          output_param_filename)

    print np.nanmax(stats['logL'])
    print np.sort(stats['logL'])
    print np.argsort(stats['likelihoods'])
    print stats['likelihoods']
    print len(stats['likelihoods'])
    print param_names
    print model_names
    print stats['model_names'][np.argmax(stats['likelihoods'])]
    print params[np.argmax(stats['likelihoods'])]

    if 0:
        import matplotlib.pyplot as plt
        plt.plot(model_names_int)
        plt.plot(stats['likelihoods'])
        plt.savefig('likelihood_test.png')

if __name__ == '__main__':
    main()
