#!/usr/bin/python

def read_ds9_region(filename):

    ''' Converts DS9 region file into format for plotting region.

    Need the following format:
        angle : degrees
        xy : pixels
        width : pixels
        height : pixels

    Region file provides following format:
        # Region file format: DS9 version 4.1
        global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
        fk5
        box(4:17:04.740,+29:20:31.32,5854.33",11972.7",130) # text={test}

    pyregion module reads DS9 regions:
    http://leejjoon.github.io/pyregion/users/overview.html


    '''

    # Import external modules
    import pyregion as pyr

    # Read region file
    try:
        region = pyr.open(filename)
    except IOError:
        return None

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    return region

def load_ds9_core_region(cores, filename='',
        header=None):

    from myimage_analysis import get_pix_coords

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    regions = read_ds9_region(filename)

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        core = tag[tag.find('{')+1:tag.find('}')]

        if core in cores:

            # Format vertices to be 2 x N array
            poly_verts = []
            for i in xrange(0, len(region.coord_list)/2):
                poly_verts.append((region.coord_list[2*i],
                                   region.coord_list[2*i+1]))

            poly_verts_pix = []
            for i in xrange(0, len(poly_verts)):
                poly_verts_pix.append(get_pix_coords(ra=poly_verts[i][0],
                                          dec=poly_verts[i][1],
                                          header=header)[:-1][::-1].tolist())

            cores[core]['poly_verts'] = {}
            cores[core]['poly_verts']['wcs'] = poly_verts
            cores[core]['poly_verts']['pixel'] = poly_verts_pix

    return cores

def run_model_analysis(cloud_dict):



'''

'''

def load_results():

    import pickle
    from astropy.io import fits

    # define directory locations
    # --------------------------
    results_dir = '/d/bip3/ezbc/multicloud/data/python_output/' + \
                  'bootstrap_results/'

    cloud_names = ('perseus', 'california', 'taurus')
    global_results = {}

    for cloud_name in cloud_names:
        results_filename = results_dir + cloud_name + \
                           '_planck_noint_gaussrange_bootstrap_results.pickle'
        with open(results_filename, 'rb') as f:
            results = pickle.load(f)
        f.close()

        results['data']['av_data'], results['data']['av_header'] = \
                fits.getdata(results['filenames']['av_filename'],
                             header=True,
                             ignore_blank=False,
                             )
        results['data']['av_error_data'], results['data']['av_error_header'] = \
                fits.getdata(results['filenames']['av_error_filename'],
                             ignore_blank=False,
                             header=True)
        results['data']['av_data_ref'] = \
                fits.getdata(results['filenames']['av_ref_filename'])
        results['data']['hi_data'], results['data']['hi_header'] = \
                fits.getdata(results['filenames']['hi_filename'],
                             ignore_blank=False,
                             header=True)
        results['data']['co_data'], results['data']['co_header'] = \
                fits.getdata(results['filenames']['co_filename'],
                             ignore_blank=False,
                             header=True)

        global_results[cloud_name] = results

    return global_results

def add_cores_to_plot(results):

    '''

    '''

    for cloud in results:
        results[cloud]['plot_kwargs']['cores_to_plot'] = \
                {'taurus': {
                     'L1495',
                     'L1495A',
                     'B213',
                     'L1498',
                     'B215',
                     'B18',
                     'B217',
                     'B220-1',
                     'B220-2',
                     'L1521',
                     'L1524',
                     'L1527-1',
                     'L1527-2',},
                  'california': {
                     'L1536',
                     'L1483',
                     'L1478',
                     'L1456',
                     'NGC1579',
                     'L1545',
                     'L1517',
                     'L1512',
                     'L1523',
                     'L1512', },
                  'perseus': {
                     'B5',
                     'IC348',
                     'B1E',
                     'B1',
                     'NGC1333',
                     'L1482',},
                  }

    return None

def add_init_model_params(results):

    '''

    '''
    for cloud_name in results:
        N_monte_carlo_runs = 1000 # Number of monte carlo runs
        vary_alphaG = True # Vary alphaG in S+14 fit?
        vary_Z = False # Vary metallicity in S+14 fit?
        vary_phi_g = False # Vary phi_g in S+14 fit?
        # Error method:
        # options are 'edges', 'bootstrap'
        error_method = 'edges'
        alpha = 0.32 # 1 - alpha = confidence
        guesses=(1.0, 1.0, 1.0) # Guesses for (alphaG, Z, phi_g)
        h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

        # Monte carlo results file bases
        results_filename = '/d/bip3/ezbc/multicloud/' + \
                           '/data/python_output/' + \
                           'monte_carlo_results/' + \
                           cloud_name + '_mc_results_' + \
                           'planck' + '_'

        sternberg_params = {}
        sternberg_params['N_monte_carlo_runs'] = N_monte_carlo_runs
        sternberg_params['param_vary'] = [vary_alphaG, vary_Z, vary_phi_g]
        sternberg_params['error_method'] = error_method
        sternberg_params['alpha'] = alpha
        sternberg_params['guesses'] = guesses
        sternberg_params['h_sd_fit_range'] = h_sd_fit_range
        sternberg_params['results_filename'] = results_filename
        sternberg_params['parameters'] = ['alphaG', 'Z', 'phi_g']

        # Krumholz Parameters
        # --------------------
        vary_phi_cnm = True # Vary phi_cnm in K+09 fit?
        vary_Z = False # Vary metallicity in K+09 fit?
        vary_phi_mol = False # Vary phi_mol in K+09 fit?
        # Error method:
        # options are 'edges', 'bootstrap'
        error_method = 'edges'
        alpha = 0.32 # 1 - alpha = confidence
        guesses=(10.0, 1.0, 10.0) # Guesses for (phi_cnm, Z, phi_mol)
        h_sd_fit_range = [0.001, 1000] # range of fitted values for sternberg model

        krumholz_params = {}
        krumholz_params['N_monte_carlo_runs'] = N_monte_carlo_runs
        krumholz_params['param_vary'] = [vary_phi_cnm, vary_Z, vary_phi_mol]
        krumholz_params['error_method'] = error_method
        krumholz_params['alpha'] = alpha
        krumholz_params['guesses'] = guesses
        krumholz_params['h_sd_fit_range'] = h_sd_fit_range
        krumholz_params['results_filename'] = results_filename
        krumholz_params['parameters'] = ['phi_cnm', 'Z', 'phi_mol']

        results[cloud_name]['model_fitting'] = {}
        results[cloud_name]['model_fitting']['results_filename'] = \
                results_filename
        results[cloud_name]['model_fitting']['krumholz_params'] = \
                krumholz_params
        results[cloud_name]['model_fitting']['sternberg_params'] = \
                sternberg_params

def mask_data(results):

    from myimage_analysis import calc_region_mask
    import numpy as np

    for cloud_name in results:
        cloud_dict = results[cloud_name]
        region_mask = cloud_dict['data_products']['region_mask']

        cloud_dict['data']['av_data'][region_mask] = np.nan
        cloud_dict['data']['av_error_data'][region_mask] = np.nan
        cloud_dict['data']['hi_data'][:, region_mask] = np.nan
        cloud_dict['data']['co_data'][:, region_mask] = np.nan

def add_core_properties(results):

    from myimage_analysis import load_ds9_region, get_pix_coords
    import json

    box_method = 'ds9'
    core_dir = '/d/bip3/ezbc/multicloud/data/python_output/core_properties/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/core_properties/'

    for cloud_name in results:
        cloud_dict = results[cloud_name]

        # define core properties
        with open(core_dir + cloud_name + '_core_properties.txt', 'r') as f:
            cores = json.load(f)

        header = cloud_dict['data']['av_header']
        #cores = convert_core_coordinates(cores, header)

        # convert center WCS coordinate to pixel
        for core in cores:
            cores[core].update({'box_pixel': 0})
            cores[core].update({'center_pixel': 0})

            center_wcs = cores[core]['center_wcs']

            # convert centers to pixel coords
            center_pixel = get_pix_coords(ra=center_wcs[0],
                                          dec=center_wcs[1],
                                          header=header)[:2]
            cores[core]['center_pixel'] = center_pixel.tolist()

        # load the bounding regions
        cores = load_ds9_core_region(cores,
                                filename=region_dir + \
                                         cloud_name + '_av_poly_cores.reg',
                                header=header)

        # Trim down cores to keep list to include only cores for cloud
        cores_to_keep_old =\
                list(cloud_dict['plot_kwargs']['cores_to_plot'][cloud_name])
        for core in cores_to_keep_old:
            if core not in cores:
                clist = cloud_dict['plot_kwargs']['cores_to_plot'][cloud_name]
                clist.remove(core)

        cloud_dict['cores'] = cores

def get_core_results(results, clobber=False):

    import pickle
    import myio

    for cloud_name in results:
        cloud_dict = results[cloud_name]
        model_analysis_filename = \
                cloud_dict['model_fitting']['results_filename']
        exists = \
            myio.check_file(model_analysis_filename,
                            clobber=clobber)

        if not exists:
            run_model_analysis(cloud_dict)
        else:
            with open(model_analysis_filename, 'rb') as f:
                cloud_dict['model_fitting'] = pickle.load(f)
            f.close()

def main():

    clobber = True

    # Get the results for each cloud
    results = load_results()

    # select which cores to plot
    add_cores_to_plot(results)

    # Initialize guesses, and params necessary for model analysis
    add_init_model_params(results)

    # mask the data by region
    mask_data(results)

    # load core positions and regions
    add_core_properties(results)

    # get the core results
    get_core_results(results, clobber=clobber)


if __name__ == '__main__':
    main()

