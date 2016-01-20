#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import mygeometry as myg
import pickle
import warnings
warnings.filterwarnings('ignore')

def plot_params_map(header=None, av_image=None, df=None, limits=None,
        filename=None, vlimits=(None,None), contours=None, parameter='phi_cnm'):

    # Import external modules
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import AxesGrid
    import pywcsgrid2 as wcs
    from matplotlib.patches import Polygon
    import matplotlib.patheffects as PathEffects
    import myplotting as myplt

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Color map
    cmap = plt.cm.copper

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    font_scale = 9

    # Create figure instance
    fig = plt.figure(figsize=(3.6, 5))

    nrows_ncols=(3,1)
    ngrids=3
    axesgrid = AxesGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0.1,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=False)

    # ------------------
    # Av image
    # ------------------
    parameters = ['phi_cnm', 'alphaG', 'hi_transition']
    for i in xrange(ngrids):
        # create axes
        ax = axesgrid[i]
        # show the image
        X, Y = np.meshgrid(np.arange(av_image.shape[1]),
                           np.arange(av_image.shape[0]))

        im = ax.contour(X, Y, av_image,
                origin='lower',
                levels=contours,
                cmap=myplt.truncate_colormap(plt.cm.binary, minval=0.3,)
                #vmin=vlimits[0],
                #vmax=vlimits[1],
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension [J2000]',)
        ax.set_ylabel('Declination [J2000]',)

        ax.locator_params(nbins=4)

        # plot limits
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits, header)
            ax.set_xlim(limits_pix[0],limits_pix[1])
            ax.set_ylim(limits_pix[2],limits_pix[3])

        # Plot cores for each cloud
        # -------------------------
        from matplotlib.patches import Circle
        from matplotlib.collections import PatchCollection

        parameter = parameters[i]
        patches = get_patches(df, header)
        cmap = myplt.truncate_colormap(plt.cm.copper, minval=0.1, maxval=1.0)
        collection = PatchCollection(patches,
                                     cmap=cmap,
                                     edgecolors='none',
                                     zorder=1000,
                                     )

        collection.set_array(df[parameter])
        ax.add_collection(collection,
                          )

        # colorbar
        #cbar = axesgrid.cbar_axes[i].colorbar(collection)
        cbar = ax.cax.colorbar(collection)
        if parameter == 'phi_cnm':
            cbar.set_label_text(r'$\phi_{\rm CNM}$',)
        elif parameter == 'alphaG':
            cbar.set_label_text(r'$\alpha G$',)
        elif parameter == 'hi_transition':
            cbar.set_label_text(r'$\Sigma_{\rm HI,trans}$ ' + \
                                r'$[M_\odot\,{\rm pc}^{-2}]$',)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=100)

def plot_ISMparams_map(header=None, av_image=None, df=None, core_dict=None,
        limits=None, filename=None, vlimits=(None,None), contours=None,
        parameter='phi_cnm'):

    # Import external modules
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import AxesGrid
    import pywcsgrid2 as wcs
    from matplotlib.patches import Polygon
    import matplotlib.patheffects as PathEffects
    import myplotting as myplt

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Color map
    cmap = plt.cm.copper

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    font_scale = 9

    # Create figure instance
    fig = plt.figure(figsize=(3.6, 5))

    nrows_ncols=(3,1)
    ngrids=3
    axesgrid = AxesGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0.1,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=False)

    # ------------------
    # Av image
    # ------------------
    parameters = ['rad_field', 'n_H', 'T_H']
    for i in xrange(ngrids):
        # create axes
        ax = axesgrid[i]
        # show the image
        X, Y = np.meshgrid(np.arange(av_image.shape[1]),
                           np.arange(av_image.shape[0]))

        im = ax.contour(X, Y, av_image,
                origin='lower',
                levels=contours,
                cmap=myplt.truncate_colormap(plt.cm.binary, minval=0.3,)
                #vmin=vlimits[0],
                #vmax=vlimits[1],
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension [J2000]',)
        ax.set_ylabel('Declination [J2000]',)

        ax.locator_params(nbins=4)

        # plot limits
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits, header)
            ax.set_xlim(limits_pix[0],limits_pix[1])
            ax.set_ylim(limits_pix[2],limits_pix[3])

        # Plot cores for each cloud
        # -------------------------
        from matplotlib.patches import Circle
        from matplotlib.collections import PatchCollection

        parameter = parameters[i]
        patches = get_patches(df, header)
        cmap = myplt.truncate_colormap(plt.cm.copper, minval=0.1, maxval=1.0)
        collection = PatchCollection(patches,
                                     cmap=cmap,
                                     edgecolors='none',
                                     zorder=1000,
                                     )

        # set values in collection
        collection_values = []
        for core in df['core']:
            collection_values.append(core_dict[core][parameter])

        collection.set_array(np.array(collection_values))
        ax.add_collection(collection,
                          )

        # colorbar
        #cbar = axesgrid.cbar_axes[i].colorbar(collection)
        cbar = ax.cax.colorbar(collection)
        if parameter == 'rad_field':
            cbar.set_label_text(r'$I_{UV}$ [$I_{UV,D}$]',)
        elif parameter == 'n_H':
            cbar.set_label_text(r'$n_H$ [cm$^{-3}$]',)
        elif parameter == 'T_H':
            cbar.set_label_text(r'$T_H$ [1,000 K]',)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=100)

def plot_ISMparams_map(header=None, av_image=None, df=None, core_dict=None,
        limits=None, filename=None, vlimits=(None,None), contours=None,
        parameter='phi_cnm', models=['krumholz', 'sternberg']):

    # Import external modules
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import AxesGrid
    import pywcsgrid2 as wcs
    from matplotlib.patches import Polygon
    import matplotlib.patheffects as PathEffects
    import myplotting as myplt

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()

    # Color map
    cmap = plt.cm.copper

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    font_scale = 9

    # Create figure instance
    fig = plt.figure(figsize=(3.5, 10))


    if 0:
        parameters = []
        if 'krumholz' in models:
            parameters.append('phi_cnm', 'alphaG', 'hi_transition', 'n_H', 'T_H']

    parameters = ['phi_cnm', 'alphaG', 'hi_transition', 'n_H', 'T_H']


    nrows_ncols=(5,1)
    ngrids=5
    axesgrid = AxesGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0.1,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=False)

    # ------------------
    # Av image
    # ------------------
    for i in xrange(ngrids):
        # create axes
        ax = axesgrid[i]
        # show the image
        X, Y = np.meshgrid(np.arange(av_image.shape[1]),
                           np.arange(av_image.shape[0]))

        im = ax.contour(X, Y, av_image,
                origin='lower',
                levels=contours,
                cmap=myplt.truncate_colormap(plt.cm.binary, minval=0.3,)
                #vmin=vlimits[0],
                #vmax=vlimits[1],
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension [J2000]',)
        ax.set_ylabel('Declination [J2000]',)

        ax.locator_params(nbins=4)

        # plot limits
        if limits is not None:
            limits_pix = myplt.convert_wcs_limits(limits, header)
            ax.set_xlim(limits_pix[0],limits_pix[1])
            ax.set_ylim(limits_pix[2],limits_pix[3])

        # Plot cores for each cloud
        # -------------------------
        from matplotlib.patches import Circle
        from matplotlib.collections import PatchCollection

        parameter = parameters[i]
        patches = get_patches(df, header)
        cmap = myplt.truncate_colormap(plt.cm.copper, minval=0.1, maxval=1.0)
        collection = PatchCollection(patches,
                                     cmap=cmap,
                                     edgecolors='none',
                                     zorder=1000,
                                     )

        # set values in collection
        collection_values = []
        if parameter in ['phi_cnm', 'alphaG', 'hi_transition']:
            collection.set_array(df[parameter])
        else:
            for core in df['core']:
                collection_values.append(core_dict[core][parameter])


            collection.set_array(np.array(collection_values))

        ax.add_collection(collection,
                          )

        # colorbar
        #cbar = axesgrid.cbar_axes[i].colorbar(collection)
        cbar = ax.cax.colorbar(collection)
        if parameter == 'rad_field':
            cbar.set_label_text(r'$I_{UV}$ [$I_{UV,D}$]',)
        elif parameter == 'n_H':
            cbar.set_label_text(r'$n_H$ [cm$^{-3}$]',)
        elif parameter == 'T_H':
            cbar.set_label_text(r'$T_H$ [1,000 K]',)
        elif parameter == 'phi_cnm':
            cbar.set_label_text(r'$\phi_{\rm CNM}$',)
        elif parameter == 'alphaG':
            cbar.set_label_text(r'$\alpha G$',)
        elif parameter == 'hi_transition':
            cbar.set_label_text(r'$\Sigma_{\rm HI,trans}$ ' + \
                                r'$[M_\odot\,{\rm pc}^{-2}]$',)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=100)

def get_patches(df, header):

    df = add_pix_coords(df, header)

    from matplotlib.patches import Circle

    patches = []
    for i in xrange(len(df['ra_pix'])):
        patches.append(Circle((df['ra_pix'][i], df['dec_pix'][i]), radius=4))

    return patches

def calc_model_color_cycle(df):

    phi_cnm_max = np.max(df['phi_cnm'])
    phi_cnm_min = np.min(df['phi_cnm'])

    color_cycle = myplt.set_color_cycle(num_colors=len(df['phi_cnm']),
                                        cmap=plt.cm.gnuplot2,
                                        )

def add_pix_coords(df, header):

    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from astropy.io import fits
    from astropy.wcs import WCS

    wcs_header = WCS(header)

    df['ra_pix'] = np.empty(len(df['ra']))
    df['dec_pix'] = np.empty(len(df['dec']))

    coords_wcs = SkyCoord(df['ra'], df['dec'], unit='deg', frame='fk5')
    coords_pix = coords_wcs.to_pixel(wcs_header)
    df['ra_pix'], df['dec_pix'] = coords_pix

    return df

def load_table():

    # Load the table with the cores
    results_dir =  '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = results_dir + 'tables/multicloud_model_params.pickle'

    df = pd.load(filename)

    return df

def check_region(pos, regions):

    selected_region = None
    if pos[1] > 0:
        for region in regions:
            within_region = myg.point_in_polygon(pos,
                                                 regions[region]['vertices']['wcs'])
            if within_region:
                selected_region = region

            within_region = False



    result = selected_region

    return result

def load_regions():

    import pyregion as pyr

    data_dir = '/d/bip3/ezbc/multicloud/data/cold_clumps/'

    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = region_dir + 'multicloud_divisions_coldcore_selection.reg'

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    region_dict = {}

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        region_name = tag[tag.find('{')+1:tag.find('}')].lower()

        if region_name in ('taurus', 'california', 'perseus'):
            # Format vertices to be 2 x N array
            poly_verts = []
            for i in xrange(0, len(region.coord_list)/2):
                poly_verts.append((region.coord_list[2*i],
                                   region.coord_list[2*i+1]))

            if 0:
                poly_verts_pix = []
                for i in xrange(0, len(poly_verts)):
                    poly_verts_pix.append(get_pix_coords(ra=poly_verts[i][0],
                                        dec=poly_verts[i][1],
                                        header=header)[:-1][::-1].tolist())

            region_dict[region_name] = {}
            region_dict[region_name]['vertices'] = {}
            region_dict[region_name]['vertices']['wcs'] = np.array(poly_verts)
            #region_dict['regions'][region_name]['poly_verts']['pixel'] = \
            #    poly_verts_pix

    return region_dict

def load_cold_cores():

    # summary of cold clump data
    # http://wiki.cosmos.esa.int/planckpla2015/index.php/Catalogues#Individual_catalogues

    table_dir = '/d/bip3/ezbc/multicloud/data/cold_clumps/'
    df_dir = '/d/bip3/ezbc/multicloud/data/python_output/tables/'
    filename = table_dir + 'HFI_PCCS_GCC_R2.02.fits'

    if 0:
        print('\nAnalyzing table...')
        cc_hdu = fits.open(filename)

        cc_data = cc_hdu[1].data

        # get the region vertices
        regions = load_regions()

        df = dict()
        df['Glon'] = []
        df['Glat'] = []
        df['ra'] = []
        df['dec'] = []
        df['Region'] = []
        df['SNR'] = []
        for i in xrange(len(cc_data)):
            #if myg.point_in_polygon(ra, region_vertices):
            ra = cc_data[i][3]
            dec = cc_data[i][4]
            #if ra < 80 and ra > 40 and dec < 45 and dec > 15:
            region_check = check_region((ra, dec), regions)
            if region_check is not None:
                df['Glon'].append(cc_data.field('GLON')[i])
                df['Glat'].append(cc_data.field('GLAT')[i])
                df['ra'].append(cc_data.field('RA')[i])
                df['dec'].append(cc_data.field('DEC')[i])
                df['SNR'].append(cc_data.field('SNR')[i])
                df['Region'].append(region_check)

        df = pd.DataFrame(df)

        df.save(df_dir + 'multicloud_cold_clumps.pickle')
    else:
        df = pd.load(df_dir + 'multicloud_cold_clumps.pickle')

    print('\nFinished loading...')

    return df

def load_av_data():

    data_dir = '/d/bip3/ezbc/multicloud/data/av/'
    filename = data_dir + 'multicloud_av_planck_tau353_5arcmin.fits'

    av_data, av_header = fits.getdata(filename,
                                      header=True)

    return av_data, av_header

def calc_core_pixel_locations(df, header):

    # Make a galactic coords object and convert to Ra/dec
    coords_gal = SkyCoord(df['Glon'] * u.deg, df['Glat'] * u.deg,
                        frame='galactic',
                        )
    coords_fk5 = coords_gal.transform_to('fk5')

    # Create WCS object
    wcs_header = WCS(header)

    # convert to pixel
    coords_pixel = coords_fk5.to_pixel(wcs_header)

    # write data to dataframe
    df['ra'] = coords_fk5.ra.deg
    df['dec'] = coords_fk5.dec.deg
    df['xpix'] = coords_pixel[0]
    df['ypix'] = coords_pixel[1]

    return df

def crop_to_random_cores(df):

    core_sample = {}

    if 0:
        df.sort(['SNR'], ascending=[True])
        df = df._slice(slice(0, 40))
    #for region in ('TMC', 'PMC', 'CMC'):
    for region in ('taurus', 'california', 'perseus'):
        row_indices = np.where((df.Region == region))[0]
        if 0:
            row_indices = np.random.choice(row_indices,
                                           replace=False,
                                           size=15,
                                           #size=len(row_indices)
                                           )
            core_sample[region] = df._slice(row_indices)
        else:
            core_sample[region] = df._slice(row_indices)
            core_sample[region].sort(['SNR'], ascending=[True])
            core_sample[region] = core_sample[region]._slice(slice(0, 15))

    return core_sample

def load_cores():

    import pickle

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    filename = table_dir + 'tables/multicloud_model_summary.pickle'

    with open(filename, 'rb') as f:
        core_dict = pickle.load(f)

    if 'n_H' not in core_dict[core_dict.keys()[0]]:
        raise ValueError('Need interpretation, run' + \
                         'cloud_write_parameter_table.py first')

    return core_dict

def plot_cdfs(core_dict, df):

    import myplotting as myplt
    import matplotlib.pyplot as plt
    import mystats


    # load Stanimirovic temps
    filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
            'stanimirovic14_temps.npy'
    stani_temps = np.loadtxt(filename)
    filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
            'stanimirovic14_temp_errors.npy'
    stani_temp_errors = np.loadtxt(filename)
    stani_temp_errors = np.array((stani_temp_errors, stani_temp_errors)).T
    filename = \
            '/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
            'stanimirovic14_int_temps.npy'
    stani_tempKs = np.loadtxt(filename)

    # load the dict
    with open('/d/bip3/ezbc/multicloud/data/python_output/tables/' + \
              'stanimirovic14_temps.pickle', 'rb') as f:
        temp_dict = pickle.load(f)

    # concatenate the temps from each cloud
    Ts_list = []
    Ts_avg_list = []
    for cloud in temp_dict:
        for Ts in temp_dict[cloud]['Ts_list']:
            Ts_list.append(Ts)
        for Ts_avg in temp_dict[cloud]['Ts_avg_list']:
            Ts_avg_list.append(Ts_avg)

    print '\nNumber of LOS:', len(Ts_avg_list)

    # collect data
    T_cnms = np.empty(len(core_dict))
    T_Hs = np.empty(len(core_dict))
    T_cnm_errors = np.empty((len(core_dict), 2))
    T_H_errors = np.empty((len(core_dict), 2))
    n_cnms = np.empty(len(core_dict))
    n_Hs = np.empty(len(core_dict))
    n_cnm_errors = np.empty((len(core_dict), 2))
    n_H_errors = np.empty((len(core_dict), 2))


    for i, core_name in enumerate(core_dict):
        core = core_dict[core_name]
        T_cnms[i] = core['T_cnm']
        T_Hs[i] = core['T_H'] * 1000.0
        T_cnm_errors[i] = core['T_cnm_error']
        T_H_errors[i] = core['T_H_error'] * 1000.0
        n_cnms[i] = core['n_cnm']
        n_Hs[i] = core['n_H']
        n_cnm_errors[i] = core['n_cnm_error']
        n_H_errors[i] = core['n_H_error']

    # Make plot
    plt.close;plt.clf()

    # Create figure instance
    fig = plt.figure(figsize=(3.5, 2.5))
    ax = fig.add_subplot(111)
    c_cycle = myplt.set_color_cycle(4, cmap_limits=(0.0, 0.9))

    #data_list = [stani_temps, stani_tempKs, T_cnms, T_Hs]
    data_list = [Ts_list, Ts_avg_list, T_cnms, T_Hs]
    data_error_list = [stani_temp_errors, np.array(((5),)),
                       T_cnm_errors, T_H_errors]
    if 0:
        data_names = [r'Stanimirovic+14 $T_{\rm spin}$',
                      r'Stanimirovic+14 $T_K$',
                      r'$T_{\rm CNM}$',
                      r'$T_H$']
    else:
        data_names = [r'Stanimirovi$\acute{c}$+14' + '\n' + r'$T_{\rm s}$',
                      r'Stanimirovi$\acute{c}$+14' + '\n' + r'$<T_{\rm s}>$',
                      r'$T_{\rm CNM}$',
                      r'$T_H$']
    linestyles = ['--', '-.', '-', '-']

    # Plot MC sim of error?
    include_errors = True
    P_low = 1960
    P = 3000.
    P_hi = 4810

    for i, data in enumerate(data_list):
        data = np.array(data)
        if i < 2:
            include_errors = False
        else:
            include_errors = True
        if include_errors:
            if i < 2:
                n = 3000 / data
            elif i == 2:
                n = n_cnms
                n_error = n_cnm_errors[:, 0]
            elif i ==3:
                n = n_Hs
                n_error = n_H_errors[:, 0]

            n_sim = 100
            alpha = 1.0 / n_sim
            alpha = 0.1
            for j in xrange(n_sim):
                if 0:
                    error = data_error_list[i]
                    data_sim = data + np.random.normal(scale=error[:, 0])
                else:
                    #pressure = np.random.uniform(P_low, P_hi, size=data.size)
                    pressure = P + np.random.normal(scale=(P - P_low),
                                                    size=data.size)
                    data_sim = pressure / (n + np.random.normal(n_error))

                x = myplt.plot_cdf(data_sim,
                                   ax=ax,
                                   plot_kwargs={#'label': label,
                                                'color': c_cycle[i],
                                                'alpha': alpha})
        if i < 2:
            x = myplt.plot_cdf(data[data > 0],
                               ax=ax,
                               plot_kwargs={'label': data_names[i],
                                            'color': c_cycle[i],
                                            'linestyle': linestyles[i],
                                            'linewidth': 2,
                                            'zorder': 1000})
        else:
            x = myplt.plot_cdf((0,0),
                               ax=ax,
                               plot_kwargs={'label': data_names[i],
                                            'color': c_cycle[i],
                                            'linestyle': linestyles[i]})

        if 0:
            sort_indices = np.argsort(data)
            data_error = data_error_list[i][sort_indices]
            data_low = data - data_error[:, 0]
            data_hi = data + data_error[:, 1]

            cdf_low = mystats.calc_cdf(data_low)
            cdf_hi = mystats.calc_cdf(data_hi)

        if 0:
            ax.fill_between(x,
                            cdf_hi,
                            cdf_low,
                            #where=where,
                            facecolor=c_cycle[i],
                            edgecolor='none',
                            alpha=0.3,
                            interpolate=True,
                            zorder=0,
                            )

    ax.legend(loc='best')
    ax.set_xscale('log')
    ax.set_xlim([7*10**-1, 10**4])
    ax.set_ylim([0,1])

    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Cumulative Distribution')

    #plt.savefig('/d/bip3/ezbc/multicloud/figures/temps/temps_cdf.png')
    plt.savefig('/d/bip3/ezbc/multicloud/figures/temps/temps_cdf.pdf')

def main():

    # get core data
    df = load_table()

    # get core dict
    core_dict = load_cores()

    # get av data
    av_data, av_header = load_av_data()

    # compare cdfs of predicted vs. measured temperatures
    print('\nPlotting Temperature CDFs...')
    plot_cdfs(core_dict, df)

    # plot the cores
    print('\nPlotting maps...')
    filetypes = ['png', 'pdf']
    for filetype in filetypes:
        figure_dir = '/d/bip3/ezbc/multicloud/figures/'
        if 0:
            filename = figure_dir + 'maps/multicloud_av_modelparams_map.' + \
                       filetype
            plot_params_map(header=av_header,
                           av_image=av_data,
                           df=df,
                           #limits=[75, 50, 20, 37,],
                           limits=[76, 43.5, 19.5, 38,],
                           filename=filename,
                           contours=[2, 4, 8, 16],
                           #vlimits=[-0.1,18]
                           )

        elif 1:
            filename = figure_dir + 'maps/multicloud_av_ISMparams_map.' + \
                       filetype
            plot_ISMparams_map(header=av_header,
                               av_image=av_data,
                               df=df,
                               core_dict=core_dict,
                               #limits=[75, 50, 20, 37,],
                               limits=[76, 43.5, 19.5, 38,],
                               filename=filename,
                               contours=[2, 4, 8, 16],
                               #vlimits=[-0.1,18]
                               )
        else:
            filename = figure_dir + 'maps/multicloud_av_ISMparams_K+09_map.' + \
                       filetype
            plot_ISMparams_map(header=av_header,
                               av_image=av_data,
                               df=df,
                               core_dict=core_dict,
                               #limits=[75, 50, 20, 37,],
                               limits=[76, 43.5, 19.5, 38,],
                               filename=filename,
                               contours=[2, 4, 8, 16],
                               #vlimits=[-0.1,18]
                               models='krumholz',
                               )

            filename = figure_dir + 'maps/multicloud_av_ISMparams_S+14_map.' + \
                       filetype
            plot_ISMparams_map(header=av_header,
                               av_image=av_data,
                               df=df,
                               core_dict=core_dict,
                               #limits=[75, 50, 20, 37,],
                               limits=[76, 43.5, 19.5, 38,],
                               filename=filename,
                               contours=[2, 4, 8, 16],
                               #vlimits=[-0.1,18]
                               models='sternberg',
                               )

if __name__ == '__main__':
    main()

