#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
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
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension [J2000]',)
    ax.set_ylabel('Declination [J2000]',)

    ax.locator_params(nbins=6)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    cb.set_label_text(r'$A_V$ [mag]',)

    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

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
                        s=100,
                        marker='+',
                        linewidths=1.5)


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

def load_cold_cores():


    # summary of cold clump data
    # http://wiki.cosmos.esa.int/planckpla2015/index.php/Catalogues#Individual_catalogues

    data_dir = '/d/bip3/ezbc/multicloud/data/cold_clumps/'
    filename = data_dir + 'HFI_PCCS_GCC_R2.02.fits'

    cc_hdu = fits.open(filename)

    cc_data = cc_hdu[1].data

    df = pd.DataFrame()
    df['Glon'] = np.empty(len(cc_data))
    df['Glat'] = np.empty(len(cc_data))
    df['ra'] = np.empty(len(cc_data))
    df['dec'] = np.empty(len(cc_data))
    for i in xrange(len(cc_data)):
        #if myg.point_in_polygon(ra, region_vertices):
        df['Glon'][i] = cc_data[i][1]
        df['Glat'][i] = cc_data[i][2]
        df['ra'][i] = cc_data[i][3]
        df['dec'][i] = cc_data[i][4]

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

    for region in ('TMC', 'PMC', 'CMC'):
        row_indices = np.where((df.Region == region))[0]
        row_indices = np.random.choice(row_indices,
                                       replace=False,
                                       #size=10,
                                       size=len(row_indices)
                                       )

        core_sample[region] = df._slice(row_indices)

    return core_sample

def main():

    # get core data
    df = load_table()

    # get av data
    av_data, av_header = load_av_data()

    # get core pixel locations
    df = calc_core_pixel_locations(df, av_header)

    # crop cores based on regions and data
    load_cold_cores()

    # crop dataset to random cores
    core_sample = crop_to_random_cores(df)

    # plot the cores
    figure_dir = '/d/bip3/ezbc/multicloud/figures/'
    filename = figure_dir + 'maps/multicloud_av_cores_meng13.png'
    plot_cores_map(header=av_header,
                   av_image=av_data,
                   core_sample=core_sample,
                   #limits=[20, 36, 75, 45],
                   filename=filename,
                   vlimits=[0,20]
                   )


if __name__ == '__main__':
    main()

