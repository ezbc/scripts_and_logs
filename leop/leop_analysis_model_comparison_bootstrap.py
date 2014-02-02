#!/usr/bin/python

import numpy as np
import math
import pyfits as pf
import mystats

# Define paths:
fitsDir = '/d/bip3/ezbc/leop/data/hi/models/modelCreationFolder5/models/'
paramsDir = '/d/bip3/ezbc/leop/data/hi/models/data/'
psDir =  '/d/bip3/ezbc/leop/figures/'

# Load in best-fit data
params = np.genfromtxt(paramsDir + 'modelComparisonData.v5.txt',delimiter='',
                       usecols=(1,2,3,4,5,6,7),skip_header=1,dtype=float)
model_names = np.genfromtxt(paramsDir + 'modelComparisonData.v5.txt',delimiter='',
                       usecols=(0),dtype=str,skip_header=1)
params_header = ['model #','resid val','inc',\
                 'pa','iflat','vflat','vdisp','z0']

# Load in best-fit model
f = pf.open(fitsDir + model_names[0] + '.smooth.FITS')
model_header, model = f[0].header, f[0].data

# Load in data cube
f = pf.open('/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/' + \
            'leop.fluxRescale.16arcsec.blk.modelBin.fits')
data_header, data = f[0].header, f[0].data[:,:,:]
data[data!=data] = 0.0

# load in cube without flux rescaling to calculate the STD
fnoise = pf.open('/d/bip3/ezbc/leop/data/hi/casa/fitsImages/fluxRescale/' +\
                 'originalCubes/leop.16arcsec.fits')
hnoise,dnoise = fnoise[0].header, fnoise[0].data[:,:,:]
std = dnoise.std()

# Initialize variables for bootstrapping
bootstrapNum = 10 # Define number of bootstrap runs
data_error = data * 0.1 + std # for goodness of fit calculation
model[model != model] = 0. # replace NaNs with 0
fluxSum = 0.93 * 1.13*16**2/1.5**2 # for normalizing model fluxes
model_normal = model / model.sum() * data.sum() #fluxSum # normalize model to data flux

# now boostrap the data with the model
sigma = 3
confid_int,gofArray = mystats.bootstrap(data,model_normal,num_samples=1000,
                                        data_error=data_error,sigma=sigma)

signif_indices = np.where(params[:,0] < params[0,0] + gofArray.std()*sigma)

params_signif = params[signif_indices]
model_names_signif = params[signif_indices]

print('The following are ' + str(sigma) + ' sigma confidence parameters:')
for i, head in enumerate(params_header):
    if i == 0:
        print(head + ' : ' + str(model_names_signif.min()) + ', ' + \
                str(model_names_signif.max()))
    elif i > 1:
        mean = params_signif[:,i-1].mean()
        low = mean - params_signif[:,i-1].min()
        high = params_signif[:,i-1].max() - mean
        print(head + ' : ' + str(mean) + ', + ' + str(high) + ', - ' + str(low))

