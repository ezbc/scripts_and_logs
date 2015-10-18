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
        filename=None, vlimits=(None,None), region_dict=None,
        plot_regions=True, plot_names=True):

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
    cmap = plt.cm.gnuplot
    cmap = plt.cm.copper

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
    ax.set_xlabel('Galactic Longitude [deg]',)
    ax.set_ylabel('Galactic Latitude [deg]',)
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
    if plot_regions:
        plot_core_regions(ax, region_dict, core_sample)
    else:
        plot_core_locs(ax, region_dict, core_sample, plot_names=plot_names)

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=300)

def plot_core_regions(ax, region_dict, core_sample):

    import matplotlib.patheffects as PathEffects
    from matplotlib.patches import Polygon

    count = 0
    rects = []
    core_labels = []
    for cloud in core_sample:
        df = core_sample[cloud]
        df = df.sort(['ra'], ascending=False)
        for i, index in enumerate(df.index):
            region_name = df['Name'][index]

            if 1:
                xpix = df['xpix'][index]
                ypix = df['ypix'][index]
                print(df['Name'][index].replace('PGCC ',''))# + \
                      #': {0:.2f} deg'.format(df['ra'][index]) + \
                      #': {0:.2f} deg'.format(df['dec'][index]))

            region = region_dict[region_name]
            vertices = np.array((region['xpix'], region['ypix'])).T

            #[:, ::-1]
            rects.append(ax.add_patch(Polygon(
                            vertices,
                            facecolor='none',
                            edgecolor='w')
                            )
                        )
            core_labels.append(str(count) + ' - ' + region_name)

            n = float(vertices.shape[0])
            center_xy = (np.sum(vertices[:, 0]) / n,
                         np.sum(vertices[:, 1]) / n)

            ax.annotate(str(count),
                    #xy=[pix_coords[0], pix_coords[1]],
                    xy=center_xy,
                    xytext=(-4,-4),
                    label=region_name,
                    textcoords='offset points',
                    fontsize=9,
                    color='k',
                    path_effects=[PathEffects.withStroke(linewidth=2,
                                             foreground="w")])

            count += 1

    ax.legend(rects,
              core_labels,
              loc='lower right',
              ncol=3,
              columnspacing=0.01,
              handlelength=0.1)

def plot_core_locs(ax, region_dict, core_sample, plot_names=True):

    count = 0
    rects = []
    core_labels = []
    for cloud in core_sample:
        for i, index in enumerate(core_sample[cloud].index):
            df = core_sample[cloud]
            region_name = df['Name'][index]

            xpix = df['xpix'][index]
            ypix = df['ypix'][index]
            print(df['Name'][index].replace('PGCC ','') + \
                  ': {0:.2f} deg'.format(df['ra'][index]) + \
                  ': {0:.2f} deg'.format(df['dec'][index]))

            ax.scatter(xpix,ypix,
                    color='w',
                    s=40,
                    marker='+',
                    linewidths=0.75)
            if plot_names:
                ax.annotate(df['Name'][index].replace('PGCC ', ''),
                            xy=(xpix, ypix),
                            xytext=(5, 5),
                            textcoords='offset points',
                            color='w',
                            fontsize=5,
                            zorder=10)

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

    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'
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
        df['nh2'] = []
        df['temp'] = []
        df['temp_error'] = []
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
                df['nh2'].append(cc_data.field('NH2')[i])
                df['temp'].append(cc_data.field('TEMP_CLUMP')[i])
                temp_error = (cc_data.field('TEMP_CLUMP')[i] - \
                                cc_data.field('TEMP_CLUMP_LOW1')[i],
                              cc_data.field('TEMP_CLUMP_UP1')[i] - \
                                cc_data.field('TEMP_CLUMP')[i]
                              )
                df['temp_error'].append(temp_error)
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

def crop_to_core_sample(df, N_cores=10, load=False, previous_cores=None,
        sampling_type='cold_cores'):

    table_dir = '/d/bip3/ezbc/multicloud/data/python_output/tables/'
    df_dir = '/d/bip3/ezbc/multicloud/data/python_output/tables/'
    filename = table_dir + 'planck11_coldclumps.pickle'

    if not load:
        core_sample = {}

        for region in ('taurus', 'california', 'perseus'):
            cloud_indices = np.where((df.Region == region))[0]
            row_indices = []
            df_new = pd.DataFrame(index=xrange(N_cores),
                                  columns=list(df.keys()))

            # crop data to be within region
            df_region = df.iloc[cloud_indices]

            # organize previous core sample
            if previous_cores is not None:
                previous_cores_region = []
                for previous_core in previous_cores:
                    if previous_core in df_region['Name'].values:
                        previous_cores_region.append(previous_core)

            # Sort by temp
            if sampling_type == 'cold_cores':
                df_region = df_region.sort(['temp'], ascending=True)
                df_region = df_region[df_region['temp'] > 0]
            elif sampling_type == 'nh2_cores':
                df_region = df_region.sort(['nh2'], ascending=False)
                df_region = df_region[df_region['nh2'] > 0]


            print df_region['nh2']

            done = False
            row = 0
            i = 0
            while not done:
                if previous_cores is not None:
                    core_df = \
                        df.loc[df['Name'] == previous_cores_region.pop()]
                elif sampling_type == 'cold_cores' or \
                     sampling_type == 'nh2_cores':
                    # get next coldest core
                    core_df = df_region.iloc[[i]]
                else:
                    # get a random core
                    random_core_index = \
                        np.random.choice(cloud_indices,
                                         replace=False,
                                         size=1,
                                         )
                    core_df = df.iloc[random_core_index]

                    print core_df['nh2']

                # check if core near another one, if not add it to the list
                if row > 0:
                    core_near_another = is_core_near_another(df_new,
                                                             core_df)
                    if not core_near_another:
                        df_new.loc[row] = core_df.values
                        row += 1
                else:
                    df_new.loc[row] = core_df.values
                    row += 1

                if row >= N_cores:
                    done = True

                i += 1

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

def convert_region_wcs2pix(region_dict, header):

    for region in region_dict:
        df = region_dict[region]

        # Make a galactic coords object and convert to Ra/dec
        coords_fk5 = SkyCoord(df['ra'] * u.deg,
                              df['dec'] * u.deg,
                              frame='fk5',
                              )
        # Create WCS object
        wcs_header = WCS(header)

        # convert to pixel
        coords_pixel = coords_fk5.to_pixel(wcs_header)

        # write data to dataframe
        df['xpix'], df['ypix'] = coords_pixel[0], coords_pixel[1]

    return region_dict

def convert_region_pix2wcs(region_dict, header):

    for region in region_dict:
        df = region_dict[region]

        # Create WCS object
        wcs_header = WCS(header)

        coords_pix = np.array([df['xpix'], df['ypix']]).T

        # convert to pixel
        coords_wcs = wcs_header.all_pix2world(coords_pix, 0)

        # write data to dataframe
        df['ra'], df['dec'] = coords_wcs[:,0], coords_wcs[:,1]

    return region_dict

def load_core_regions(core_sample, header):

    import pyregion as pyr

    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'
    filename = region_dir + 'multicloud_coldclump_divisions.reg'

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    region_dict = {}
    for cloud in core_sample:
        df = core_sample[cloud]
        #df['region'] = []

        df = dict(df)

        for region in regions:
            # Cores defined in following format: 'tag={L1495A}'
            tag = region.comment
            region_name = tag[tag.find('text={')+6:tag.find('}')]

            #if region_name in df['Name'].values:
            if region_name in df['Name'].values:
                # Format vertices to be 2 x N array
                poly_verts = []
                for i in xrange(0, len(region.coord_list)/2):
                    poly_verts.append((region.coord_list[2*i],
                                       region.coord_list[2*i+1]))

                poly_verts = np.array(poly_verts)

                verts_dict = {'ra': poly_verts[:, 0], 'dec': poly_verts[:, 1]}

                region_dict[region_name] = verts_dict

    # convert regions to pixel coords
    region_dict = convert_region_wcs2pix(region_dict, header)

    return region_dict

def get_old_core_sample(use_old_sample=True):

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

    if use_old_sample:
        return cores
    else:
        return None

def derive_ideal_wedge(av_image, core_sample, wedge_angle=40, wedge_radius=10,
        av_image_error=None, core_rel_pos=0.1, angle_res=1.0, width=3):

    import mygeometry as myg
    import myimage_analysis as myim

    """
    Parameters
    ----------
    angle_res : float
        Resolution with which to rotate each new box in degrees. 1.0 degree
        gives 360 different box orientations.


    """

    angle_grid = np.arange(0, 360, angle_res)
    region_dict = {}

    for cloud_name in core_sample:
    #for cloud_name in ('perseus',):
        cloud_df = core_sample[cloud_name]
        gradient_sums_list = []
        for core_name in cloud_df['Name']:
        #for core_name in ('G158.26-21.81',):

            core = cloud_df[cloud_df['Name'] == core_name]

            print('Calculating optimal angle for core {:s}'.format(core_name))

            # Get center position in pixels
            core_pos = [core['xpix'].values[0], core['ypix'].values[0]][::-1]

            wedge_vertices = myg.create_wedge(core_pos,
                                              wedge_radius,
                                              wedge_angle,
                                              center_rel_pos=core_rel_pos,
                                              width=wedge_width,
                                              )

            gradient_sums = np.zeros((len(angle_grid)))

            for i, angle in enumerate(angle_grid):
                wedge_vertices_rotated = myg.rotate_wedge(wedge_vertices,
                                                          core_pos,
                                                          angle)

                try:
                    mask = \
                        myg.get_polygon_mask(av_image,
                                             wedge_vertices_rotated)
                    av_image_masked = np.copy(av_image)

                    # extract radial profile weighted by SNR
                    radii, profile = \
                        myim.get_radial_profile(av_image,
                                                binsize=1,
                                                center=core_pos, #[::-1],
                                                #weights=av_image_error,
                                                mask=mask
                                                )

                    if angle == 90:
                        av_image_masked = np.copy(av_image)
                        mask = myg.get_polygon_mask(av_image_masked,
                                                    wedge_vertices)
                        av_image_masked[mask==0] = np.NaN

                    indices = np.where((radii == radii) & \
                                       (profile == profile))
                    profile, radii = profile[indices], radii[indices]

                    # steeper gradients will have smaller sums
                    gradient_sum = np.sum(np.gradient(profile, radii))
                    gradient_sums[i] = gradient_sum
                except IndexError:
                    gradient_sums[i] = 0.

                gradient_sums_list.append(gradient_sums)

                #print wedge_vertices_rotated

            # find steepest profile and recreate the box mask
            angle_ideal = angle_grid[gradient_sums == np.min(gradient_sums)][0]

            wedge_vertices_rotated = myg.rotate_wedge(wedge_vertices,
                                                      core_pos,
                                                      angle_ideal)

            region_dict[core_name] = {}
            region_dict[core_name]['xpix'] = wedge_vertices_rotated[:,1]
            region_dict[core_name]['ypix'] = wedge_vertices_rotated[:,0]

    return region_dict

def load_wedges(av_image, core_sample, wedge_angle=40, wedge_radius=10,
        av_image_error=None, core_rel_pos=0.1, angle_res=1.0, wedge_width=3,):

    region_dict = {}
    core_wedge_angles = {}
    # angles are counter-clockwise from north
    core_wedge_angles = {
                         # Taurus
                         'G174.40-13.45': 45,
                         'G174.70-15.47': 100,
                         'G174.05-15.82': 20,
                         'G172.93-16.73': 160,
                         'G172.12-16.94': 190,
                         'G171.14-17.57': 200,
                         'G171.49-14.91': 0,
                         'G171.00-15.80': 195,
                         'G169.32-16.17': 210,
                         'G168.10-16.38': 230,
                         # Perseus
                         'G160.49-16.81': 20,
                         'G160.46-17.99': 90,
                         'G159.80-18.49': 0,
                         'G160.14-19.08': 120,
                         'G160.53-19.73': 170,
                         'G159.19-20.11': 10,
                         'G158.39-20.72': 0,
                         'G159.17-21.09': 100,
                         'G158.89-21.60': 170,
                         'G158.26-21.81': 270,
                         # California
                         'G168.54-6.22': 100,
                         'G168.12-6.42': 0,
                         'G166.91-7.76': 60,
                         'G165.71-9.15': 180,
                         'G165.36-7.51': 90,
                         'G164.99-8.60': 250,
                         'G164.70-7.63': 40,
                         'G164.18-8.84': 290,
                         'G164.26-8.39': 0,
                         'G164.65-8.12': 140,
                         }

    for cloud_name in core_sample:
    #for cloud_name in ('perseus',):
        cloud_df = core_sample[cloud_name]
        gradient_sums_list = []
        for core_name in cloud_df['Name']:
        #for core_name in ('G158.26-21.81',):

            core = cloud_df[cloud_df['Name'] == core_name]

            #print('Calculating optimal angle for core {:s}'.format(core_name))

            # Get center position in pixels
            core_pos = [core['xpix'].values[0], core['ypix'].values[0]][::-1]

            wedge_vertices = myg.create_wedge(core_pos,
                                              wedge_radius,
                                              wedge_angle,
                                              center_rel_pos=core_rel_pos,
                                              width=wedge_width,
                                              )

            # angle of wedge
            if core_name not in core_wedge_angles:
                core_wedge_angles[core_name] = 0
            angle_ideal = core_wedge_angles[core_name]

            wedge_vertices_rotated = myg.rotate_wedge(wedge_vertices,
                                                      core_pos,
                                                      angle_ideal)

            region_dict[core_name] = {}
            region_dict[core_name]['xpix'] = wedge_vertices_rotated[:,1]
            region_dict[core_name]['ypix'] = wedge_vertices_rotated[:,0]

    return region_dict

def calc_wedge_regions(core_sample, av_data, header):

    wedge_angle = 40.0 # degrees
    wedge_radius = 10.0 / 0.43 # pixels,
    wedge_width = 0.8 * wedge_radius # fraction of radius core is within wedge
    core_rel_pos = 0.30 # fraction of radius core is within wedge

    if 0:
        region_dict = derive_ideal_wedge(av_data,
                                         core_sample,
                                         wedge_angle=wedge_angle,
                                         wedge_radius=wedge_radius,
                                         wedge_width=wedge_width,
                                         core_rel_pos=core_rel_pos,
                                         angle_res=10.,
                                         )
    else:
        region_dict = load_wedges(av_data,
                                  core_sample,
                                  wedge_angle=wedge_angle,
                                  wedge_radius=wedge_radius,
                                  wedge_width=wedge_width,
                                  core_rel_pos=core_rel_pos,
                                  angle_res=10.,
                                  )

    region_dict = convert_region_pix2wcs(region_dict, header)


    return region_dict

def save_region_dict(region_dict):

    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'
    filename = region_dir + 'multicloud_divisions_coldcore_wedges.pickle'

    with open(filename, 'wb') as f:
        pickle.dump(region_dict, f)

def main():

    load_gcc_data = 0
    load_coresample_data = 0
    N_cores = 10
    use_old_sample = 0
    load_regions = 1
    region_type = 'wedges' # ds9 or wedges

    cores = get_old_core_sample(use_old_sample=use_old_sample)

    # get core data
    df = load_table()

    # crop cores based on regions and data
    df = load_cold_cores(load_raw_data=load_gcc_data)

    # get av data
    av_data, av_header = load_av_data()

    # get core pixel locations
    df = calc_core_pixel_locations(df, av_header)

    # crop dataset to random cores
    core_sample = crop_to_core_sample(df,
                                      N_cores=N_cores,
                                      load=load_coresample_data,
                                      previous_cores=cores,
                                      sampling_type='nh2_cores',
                                      )

    # load core regions
    if load_regions:
        if region_type == 'ds9':
            region_dict = load_core_regions(core_sample, av_header)
        else:
            region_dict = calc_wedge_regions(core_sample, av_data, av_header)

        save_region_dict(region_dict)
    else:
        region_dict = None

    # plot the cores
    print('\nPlotting...')
    figure_dir = '/d/bip3/ezbc/multicloud/figures/'
    filetypes = ['png', 'pdf']
    for filetype in filetypes:
        filename = figure_dir + 'maps/multicloud_av_cores_map.' + \
                   filetype
        plot_cores_map(header=av_header,
                       av_image=av_data,
                       core_sample=core_sample,
                       region_dict=region_dict,
                       limits=[76, 43.5, 19.5, 38,],
                       filename=filename,
                       plot_regions=load_regions,
                       plot_names=True,
                       vlimits=[0,15.5],
                       )


if __name__ == '__main__':
    main()

