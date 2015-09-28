#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import mygeometry as myg
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
    cmap = plt.cm.gnuplot2

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    font_scale = 9

    # Create figure instance
    fig = plt.figure(figsize=(7.5, 5))

    nrows_ncols=(1,1)
    ngrids=1

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
                 share_all=True)

    # ------------------
    # Av image
    # ------------------
    # create axes
    ax = axesgrid[0]
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

    ax.locator_params(nbins=6)


    # plot limits
    if limits is not None:
        limits = myplt.convert_wcs_limits(limits, header)
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Plot cores for each cloud
    # -------------------------
    if 0:
        for region in core_sample:
            for i, index in enumerate(core_sample[region].index):
                df = core_sample[region]
                xpix = df['xpix'][index]
                ypix = df['ypix'][index]

                anno_color = (0.3, 0.5, 1)

                if 1:
                    ax.scatter(xpix,ypix,
                            color='w',
                            s=40,
                            marker='+',
                            linewidths=0.75)
    if 1:
        from matplotlib.patches import Circle
        from matplotlib.collections import PatchCollection

        patches = get_patches(df, header)
        cmap = myplt.truncate_colormap(plt.cm.copper, minval=0.1)
        cmap = plt.cm.gnuplot
        collection = PatchCollection(patches,
                                     cmap=cmap,
                                     edgecolors='none',
                                     zorder=1000,
                                     )

        collection.set_array(df[parameter])
        ax.add_collection(collection,
                          )

        # colorbar
        cbar = ax.cax.colorbar(collection)
        if parameter == 'phi_cnm':
            cbar.set_label_text(r'$\phi_{\rm CNM}$',)
        else:
            cbar.set_label_text(r'$\alpha G$',)

        if 0:
            ax.scatter(df['ra_pix'], df['dec_pix'],
                    color=model_color_cycles,
                    s=10,
                    marker='^',
                    linewidths=0.75,
                    )

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=100)

def get_patches(df, header):

    df = add_pix_coords(df, header)

    from matplotlib.patches import Circle

    patches = []
    for i in xrange(len(df['ra_pix'])):
        patches.append(Circle((df['ra_pix'][i], df['dec_pix'][i]), radius=2))

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

def main():

    # get core data
    df = load_table()

    # get av data
    av_data, av_header = load_av_data()

    # plot the cores
    figure_dir = '/d/bip3/ezbc/multicloud/figures/'
    filename = figure_dir + 'maps/multicloud_av_modelparams_phicnm_map.png'
    plot_params_map(header=av_header,
                   av_image=av_data,
                   df=df,
                   limits=[75, 50, 20, 37,],
                   filename=filename,
                   contours=[2, 4, 8, 16],
                   parameter='phi_cnm',
                   #vlimits=[-0.1,18]
                   )
    filename = figure_dir + 'maps/multicloud_av_modelparams_alphaG_map.png'
    plot_params_map(header=av_header,
                   av_image=av_data,
                   df=df,
                   limits=[75, 50, 20, 37,],
                   filename=filename,
                   contours=[2, 4, 8, 16],
                   parameter='alphaG'
                   #vlimits=[-0.1,18]
                   )


if __name__ == '__main__':
    main()

