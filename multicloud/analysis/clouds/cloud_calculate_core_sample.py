#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import mygeometry as myg
import warnings
warnings.filterwarnings('ignore')

def plot_cores_map(header=None, av_image=None, core_sample=None, limits=None,
        filename=None, vlimits=(None,None)):

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
    im = ax.imshow(av_image,
            interpolation='nearest',
            origin='lower',
            cmap=cmap,
            vmin=vlimits[0],
            vmax=vlimits[1],
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("gal")
    #ax.set_ticklabel_type("hms", "dms")
    ax.set_ticklabel_type("dms", "dms")

    #ax.set_xlabel('Right Ascension [J2000]',)
    #ax.set_ylabel('Declination [J2000]',)
    ax.set_xlabel('Galactic Longitude [J2000]',)
    ax.set_ylabel('Galactic Latitude [J2000]',)
    ax.axis["top",].toggle(ticklabels=True)
    ax.grid(color='w', linewidth=0.5)
    ax.axis[:].major_ticks.set_color("w")

    ax.locator_params(nbins=6)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    cb.set_label_text(r'$A_V$ [mag]',)

    # plot limits
    if limits is not None:
        limits = myplt.convert_wcs_limits(limits, header, frame='fk5')
        ax.set_xlim(limits[0],limits[1])
        ax.set_ylim(limits[2],limits[3])

    # Plot cores for each cloud
    # -------------------------
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
                ax.annotate(df['Name'][index].replace('PGCC ', ''),
                            xy=(xpix, ypix),
                            xytext=(5, 5),
                            textcoords='offset points',
                            color='w',
                            fontsize=5,
                            zorder=10)

            print(df['Name'][index].replace('PGCC ','') + \
                  ': {0:.2f} deg'.format(df['ra'][index]) + \
                  ': {0:.2f} deg'.format(df['dec'][index]))


    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=300)

def load_table():

    # Load the table with the cores
    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/tables/'
    filename = table_dir + 'meng13_table1.txt'

    df = pd.read_csv(filename,
                     header=2,
                     skiprows=[4,],
                     delimiter='\t',
                     usecols=(1, 2, 7),
                     )
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

def load_cold_cores(load_raw_data=True):

    # summary of cold clump data
    # http://wiki.cosmos.esa.int/planckpla2015/index.php/Catalogues#Individual_catalogues

    table_dir = '/d/bip3/ezbc/multicloud/data/cold_clumps/'
    df_dir = '/d/bip3/ezbc/multicloud/data/python_output/tables/'
    filename = table_dir + 'HFI_PCCS_GCC_R2.02.fits'

    if not load_raw_data:
        print('\nAnalyzing table...')
        cc_hdu = fits.open(filename)

        cc_data = cc_hdu[1].data

        # get the region vertices
        regions = load_regions()

        df = dict()
        df['Name'] = []
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
                df['Name'].append(cc_data.field('NAME')[i].replace('PGCC ', ''))
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

def is_core_near_another(df, df_row, dist_thres=25./60.):

    ra_1 = df_row.get('ra').values
    dec_1 = df_row.get('dec').values
    core_near_another = False

    for row in df.index:
        if not np.isnan(df.loc[row].get('ra')):
            ra_0 = df.loc[row].get('ra')
            dec_0 = df.loc[row].get('dec')

            distance = ((ra_0 - ra_1)**2 + (dec_0 - dec_1)**2)**0.5

            if distance < dist_thres:
                core_near_another = True

    return core_near_another

def crop_to_random_cores(df, N_cores=15, load=False, previous_cores=None):

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/tables/'
    df_dir = '/d/bip3/ezbc/multicloud/data/python_output/tables/'
    filename = table_dir + 'planck11_coldclumps.txt'

    if not load:
        core_sample = {}

        if 0:
            df.sort(['SNR'], ascending=[True])
            df = df._slice(slice(0, 40))
        for region in ('taurus', 'california', 'perseus'):
            cloud_indices = np.where((df.Region == region))[0]
            row_indices = []
            df_new = pd.DataFrame(index=xrange(N_cores), columns=list(df.keys()))

            done = False
            row = 0
            while not done:
                if previous_cores is not None:
                    random_core_df = df.loc[df['Name'] == previous_cores.pop()]
                else:
                    # get a random core
                    random_core_index = np.random.choice(cloud_indices,
                                                   replace=False,
                                                   size=1,
                                                   #size=len(row_indices)
                                                   )
                    random_core_df = df.iloc[random_core_index]

                # check if core near another one, if not add it to the list
                if row > 0:
                    core_near_another = is_core_near_another(df_new,
                                                             random_core_df)
                    if not core_near_another:
                        df_new.loc[row] = random_core_df.values
                        row += 1
                else:
                    df_new.loc[row] = random_core_df.values
                    row += 1

                if row >= N_cores:
                    done = True

            # convert df to data frame
            core_sample[region] = df_new

            if 0:
                core_sample[region] = df._slice(cloud_indices)
                core_sample[region].sort(['SNR'], ascending=[True])
                core_sample[region] = core_sample[region]._slice(slice(0, 15))

        with open(filename, 'wb') as f:
            pickle.dump(core_sample, f)
        #df = pd.DataFrame(core_sample)
        #df.save(filename)
    else:
        with open(filename, 'rb') as f:
            core_sample = pickle.load(f)

        #df = pd.load(filename)

    return core_sample


def load_core_regions(core_sample):





    # convert regions to pixel coords
    core_sample = convert_region_wcs2pix(core_sample, av_header)

    return core_sample

def main():

    load_gcc_data = 1
    load_coresample_data = 1
    N_cores = 10

    cores = ['G166.83-8.68',
             'G168.82-6.37',
             'G168.05-7.01',
             'G164.16-8.46',
             'G165.23-8.78',
             'G167.06-7.77',
             'G168.12-6.42',
             'G167.58-6.64',
             'G164.70-7.63',
             'G166.35-8.77',
             'G166.73-15.06',
             'G173.08-16.50',
             'G172.74-14.53',
             'G169.44-16.18',
             'G173.86-17.65',
             'G173.71-13.91',
             'G171.75-14.18',
             'G173.70-15.21',
             'G170.28-19.48',
             'G171.00-15.80',
             'G158.23-20.15',
             'G159.01-22.19',
             'G159.19-20.11',
             'G157.12-23.49',
             'G160.10-19.90',
             'G160.34-18.42',
             'G158.40-21.86',
             'G159.79-21.32',
             'G158.89-21.60',
             'G159.51-18.41',
            ]

    # get core data
    df = load_table()

    # crop cores based on regions and data
    df = load_cold_cores(load_raw_data=load_gcc_data)

    # get av data
    av_data, av_header = load_av_data()

    # get core pixel locations
    df = calc_core_pixel_locations(df, av_header)

    # crop dataset to random cores
    core_sample = crop_to_random_cores(df,
                                       N_cores=N_cores,
                                       load=load_coresample_data,
                                       previous_cores=cores,
                                       )

    # load core regions
    #core_sample = load_core_regions(core_sample)

    # plot the cores
    figure_dir = '/d/bip3/ezbc/multicloud/figures/'
    filename = figure_dir + 'maps/multicloud_av_cores_meng13.png'
    plot_cores_map(header=av_header,
                   av_image=av_data,
                   core_sample=core_sample,
                   limits=[75, 45, 20, 38,],
                   filename=filename,
                   vlimits=[-0.1,18]
                   )


if __name__ == '__main__':
    main()

