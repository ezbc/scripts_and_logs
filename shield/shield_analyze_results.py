#!/usr/bin/python

import os
import numpy as np


def plot_results(stats):

    import myplotting as myplt

    stats = format_distributions(stats)


    myplt.corner_plot(stats['likelihoods_reshaped'],
                      plot_grids=stats['param_grids'],
                      labels=stats['param_names'],
                      filename='/d/bip3/ezbc/shield/figures/749237_likelihoods.png',)

def format_distributions(stats):

    likes = stats['likelihoods']
    params = stats['params']
    param_names = stats['param_names']
    param_shapes = []
    param_grids = []
    keep_params = ['iflat', 'vflat', 'z0', 'inc']
    keep_params_sorted = []

    for i, param_name in enumerate(param_names[1:]):
        if param_name in keep_params:
            param_grid = np.sort(np.unique(params[:, i]))
            print param_name
            print param_grid
            param_grids.append(param_grid)
            param_shapes.append(len(param_grid))
            keep_params_sorted.append(param_name)


    stats['likelihoods_reshaped'] = likes.reshape(param_shapes)
    stats['param_names'] = keep_params_sorted
    stats['param_grids'] = param_grids

    return stats

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
    from astropy.io import fits

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

    stats['params'] = params
    stats['param_names'] = param_names

    #print np.sort(stats['logL'])
    model_best_name = stats['model_names'][np.argmax(stats['likelihoods'])]
    print params[np.argmax(stats['likelihoods'])]

    stats['model_best'] = fits.getdata('models/' + model_best_name)

    cube_name = '749237_rebin_cube_regrid.fits'
    cube_error_name = '749237_rebin_cube_error_regrid.fits'
    data = fits.getdata('/d/bip3/ezbc/shield/749237_lowres/' + \
                        'modeling_fineres/' + \
                        cube_name,)
    error = fits.getdata('/d/bip3/ezbc/shield/749237_lowres/' + \
                         'modeling_fineres/' + cube_error_name)

    stats['rescaled_std'] = np.nansum((stats['model_best'] - data)**2 / \
                            data.size)**0.5

    stats['chisq'] = np.nansum((stats['model_best'] - data)**2 / error**2) / \
                     data.size

    print('rescaled std = ', stats['rescaled_std'])
    print('chi^2 = ', stats['chisq'])
    print('best fit model = ' + model_best_name)

    if 0:
        import matplotlib.pyplot as plt
        plt.plot(model_names_int)
        plt.plot(stats['likelihoods'])
        plt.savefig('likelihood_test.png')

    plot_results(stats)

if __name__ == '__main__':
    main()
