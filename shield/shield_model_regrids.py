#!/usr/bin/python

# Planck resolution is 5'. We will smooth and regrid all maps to the Planck map.

# Map                           Resolution      Convol size
# Lee et al. (2012) Av          5'              5
# GALFA HI                      3.7'            3.363
# CfA CO Dame et al. (2001)     8'              No smoothing

import os
import numpy as np

def main():

    from myimage_analysis import bin_image, calculate_nhi
    from mycoords import make_velocity_axis
    from astropy.io import fits
    import numpy as np
    import mystats
    import myio
    import pickle

    os.chdir('/d/bip3/ezbc/shield/749237_lowres/modeling_fineres/models')

    # Ddelete gipsy files, leaving only fits files
    filename_delete_list = os.listdir('./')
    for filename in filename_delete_list:
        if '.image' in filename or '.descr' in filename:
            os.system('rm -rf ' + filename)

    # If true, deletes files to be written
    clobber = 0

    # get leftover fits files
    filename_init_list = os.listdir('./')

    filename_list = []
    for filename in filename_init_list:
        if 'model' in filename and 'regrid' not in filename:
            filename_list.append(filename)

    #filename_list = filename_list[:5]

    stats = {
            'logL': np.empty(len(filename_list)),
            'std': np.empty(len(filename_list)),
            'mean_abs_resid': np.empty(len(filename_list)),
            'sum_abs_resid': np.empty(len(filename_list)),
            }

    cube_name = '749237_rebin_cube_regrid.fits'
    unbinned_cube_name = '749237_rebin_cube.fits'
    cube_error_name = '749237_rebin_cube_error_regrid.fits'
    data, header = fits.getdata('/d/bip3/ezbc/shield/749237_lowres/modeling_fineres/' + \
                                cube_name,
                                header=True)
    error = fits.getdata('/d/bip3/ezbc/shield/749237_lowres/modeling_fineres/' + \
                         cube_error_name)

    data_sum = np.nansum(data)
    #binsize = 6

    _, unbinned_header = fits.getdata('/d/bip3/ezbc/shield/749237_lowres/' + \
                          unbinned_cube_name, header=True)
    beamsize = unbinned_header['BMAJ']
    cdelt = np.abs(unbinned_header['CDELT1'])
    binsize = int(beamsize / cdelt)

    mask = (np.isnan(data) | np.isnan(error))

    # Load the images into miriad
    for i, model_name in enumerate(filename_list):

        model_bin_name = model_name.replace('.FITS', '_regrid.FITS')

        exists = myio.check_file(model_bin_name, clobber=clobber)
        if exists:
            print('Loading cube:\n' + model_bin_name)
            model_bin = fits.getdata(model_bin_name)
        else:
            print('Binning cube:\n' + model_name)

            model = fits.getdata(model_name)

            print('\tBinsize = ' + str(binsize))

            # bin the model
            model_bin = bin_image(model,
                                  binsize=(1, binsize, binsize),
                                  statistic=np.nanmean,
                                  quick_bin=True
                                  )

            # normalize the model to have same sum as data
            model_bin = model_bin / np.nansum(model_bin) * data_sum
            #assert np.nansum(model_bin) == data_sum

            # write the model to a file
            fits.writeto(model_bin_name,
                         model_bin,
                         header,
                         clobber=clobber)

        residuals = model_bin[~mask] - data[~mask]
        stats['logL'][i] = mystats.calc_logL(model_bin[~mask],
                                             data[~mask],
                                             data_error=error[~mask])
        stats['std'][i] = np.nanstd(residuals)
        stats['mean_abs_resid'][i] = np.nanmean(np.abs(residuals))
        stats['sum_abs_resid'][i] = np.nansum(np.abs(residuals))

    with open('/d/bip3/ezbc/shield/749237_lowres/modeling_fineres/statistics.pickle', 'wb') as f:
        pickle.dump(stats, f)

    with open('/d/bip3/ezbc/shield/749237_lowres/modeling_fineres/statistics.pickle', 'rb') as f:
        stats = pickle.load(f)

if __name__ == '__main__':
    main()
