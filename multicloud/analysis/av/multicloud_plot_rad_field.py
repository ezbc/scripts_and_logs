#!/usr/bin/python

''' Calculates the N(HI) map for multicloud

'''

import numpy as np

''' Plotting Functions
'''

def plot_rad_map(header=None, contour_image=None, rad_field=None, cores=None,
        cores_to_keep=None,
        props=None, cloud_dict=None, regions=None, title=None, limits=None,
        contours=None, boxes=False, savedir='./', filename=None, show=False,
        hi_vlimits=None, vlimits=None, vscale='linear'):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import AxesGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    import matplotlib.patheffects as PathEffects

    # Set up plot aesthetics
    # ----------------------
    #plt.close();plt.clf()

    # Color map
    cmap = plt.cm.copper

    # Color cycle, grabs colors from cmap

    # Create figure instance
    fig = plt.figure(figsize=(7.5, 5))

    if rad_field is not None:
        nrows_ncols=(1,1)
        ngrids=1

    imagegrid = AxesGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="single",
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
    ax = imagegrid[0]


    #ax = wcs.subplot(111, header=header)


    if vscale == 'log':
        norm = matplotlib.colors.LogNorm()
    else:
        norm = None

    # show the image
    im = ax.imshow(rad_field,
            interpolation='nearest',
            origin='lower',
            cmap=cmap,
            norm=norm,
            vmin=vlimits[0],
            vmax=vlimits[1],
            )

    # Asthetics
    #ax.set_display_coord_system("fk5")
    #ax.set_ticklabel_type("hms", "dms")

    #ax.set_xlabel('Right Ascension [J2000]',)
    #ax.set_ylabel('Declination [J2000]',)

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

    cmap.set_bad(color='w')

    #ax.locator_params(nbins=6)

    # colorbar
    from matplotlib.ticker import LogFormatter
    formatter = LogFormatter(10, labelOnlyBase=False)
    ticks = np.logspace(np.log10(vlimits[0]),
                                 np.log10(vlimits[1]),
                                 #np.log10(vlimits[1]) - np.log10(vlimits[0]),
                                 5,
                        )

    from matplotlib.ticker import LogFormatterMathtext

    if 1:
        cb = ax.cax.colorbar(im,
                             #ticks=ticks,
                             #format=LogFormatterMathtext(),
                             )
        #cb.set_label_text(r'log$_{10}$[$U_{0,M83}$]',)
        cb.set_label_text(r'$U_{0,M83}$',)
    else:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05)

        #plt.colorbar(im, cax=cax)

        #cb = fig.colorbar(im, ax=ax)
        #cb.ax.set_ylabel(r'$U_{0,M83}$',)

    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

        if props is not None:
            if 0:
                for region in props['region_name_pos']:
                    if region == 'taurus1':
                        region = 'taurus 1'
                    if region == 'taurus2':
                        region = 'taurus 2'
                    ax.annotate(region.capitalize(),
                                xy=props['region_name_pos'][region]['pixel'],
                                xytext=(0,0),
                                textcoords='offset points',
                                color='w',
                                fontsize=7,
                                zorder=10)

    if regions is not None:
        for region in props['region_name_pos']:
            vertices = np.copy(regions[region]['poly_verts']['pixel'])
            rect = ax.add_patch(Polygon(
                    vertices[:, ::-1],
                    facecolor='none',
                    edgecolor='w'))

   # ax.legend(rects, core_labels, bbox_to_anchor=(1.05, 1), loc=2,
   #         borderaxespad=0.)
    #ax.legend(rects, core_labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #       ncol=5, mode="expand", borderaxespad=0.)

    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        fig.show()

''' DS9 Region and Coordinate Functions
'''

def convert_limit_coordinates(prop_dict,
        coords=('region_limit', 'co_noise_limits', 'plot_limit'), header=None):

    # Initialize pixel keys
    for coord in coords:

        if coord == 'region_limit' or coord == 'plot_limit':
            prop_dict[coord].update({'pixel': []})
            limit_wcs = prop_dict[coord]['wcs']

            for limits in limit_wcs:
                # convert centers to pixel coords
                limit_pixels = get_pix_coords(ra=limits[0],
                                             dec=limits[1],
                                             header=header)[:2].tolist()

                prop_dict[coord]['pixel'].append(limit_pixels[0])
                prop_dict[coord]['pixel'].append(limit_pixels[1])
        elif coord == 'co_noise_limits':
            prop_dict[coord].update({'pixel': []})
            region_limits = prop_dict[coord]['wcs']

            # Cycle through each region, convert WCS limits to pixels
            for region in region_limits:
                region_pixels = []
                for limits in region:
                    # convert centers to pixel coords
                    limit_pixels = get_pix_coords(ra=limits[0],
                                                  dec=limits[1],
                                                  header=header)[:2].tolist()
                    region_pixels.append(limit_pixels)

                # Append individual regions back to CO noise
                prop_dict[coord]['pixel'].append(region_pixels)
        elif coord == 'region_name_pos':

            # convert centers to pixel coords
            for region in prop_dict[coord]:
                prop_dict[coord][region].update({'pixel': []})

                coord_wcs = prop_dict[coord][region]['wcs']

                coord_pixel = get_pix_coords(ra=coord_wcs[0],
                                             dec=coord_wcs[1],
                                             header=header)[:2].tolist()

                prop_dict[coord][region]['pixel'].append(coord_pixel[0])
                prop_dict[coord][region]['pixel'].append(coord_pixel[1])

    return prop_dict

def convert_core_coordinates(cores, header):

    for core in cores:
        cores[core].update({'box_pixel': 0})
        cores[core].update({'center_pixel': 0})

        #box_wcs = cores[core]['box_wcs']
        #box_pixel = len(box_wcs) * [0,]
        center_wcs = cores[core]['center_wcs']

        # convert centers to pixel coords
        center_pixel = get_pix_coords(ra=center_wcs[0],
                                      dec=center_wcs[1],
                                      header=header)[:2]
        cores[core]['center_pixel'] = center_pixel.tolist()

        # convert box corners to pixel coords
        #for i in range(len(box_wcs)/2):
        #    pixels = get_pix_coords(ra=box_wcs[2*i], dec=box_wcs[2*i + 1],
        #            header=header)
        #    box_pixel[2*i], box_pixel[2*i + 1] = int(pixels[0]), int(pixels[1])
        #cores[core]['box_pixel'] = box_pixel

    return cores

def load_fits(filename,return_header=False):
    ''' Loads a fits file.
    '''

    import pyfits as pf

    f = pf.open(filename)
    if return_header:
        return f[0].data,f[0].header
    else:
        return f[0].data

def get_sub_image(image, indices):

    return image[indices[1]:indices[3],
            indices[0]:indices[2]]

def get_pix_coords(ra=None, dec=None, header=None):

    ''' Ra and dec in (hrs,min,sec) and (deg,arcmin,arcsec), or Ra in degrees
    and dec in degrees.
    '''

    import pywcsgrid2 as wcs
    import pywcs

    # convert to degrees if ra and dec are array-like
    try:
        if len(ra) == 3 and len(dec) == 3:
            ra_deg, dec_deg = hrs2degs(ra=ra, dec=dec)
        else:
            raise ValueError('RA and Dec must be in (hrs,min,sec) and' + \
                    ' (deg,arcmin,arcsec) or in degrees.')
    except TypeError:
        ra_deg, dec_deg = ra, dec

    wcs_header = pywcs.WCS(header)
    pix_coords = wcs_header.wcs_sky2pix([[ra_deg, dec_deg, 0]], 0)[0]

    return pix_coords

def hrs2degs(ra=None, dec=None):
    ''' Ra and dec tuples in hrs min sec and deg arcmin arcsec.
    '''

    ra_deg = 15*(ra[0] + ra[1]/60. + ra[2]/3600.)
    dec_deg = dec[0] + dec[1]/60. + dec[2]/3600.

    return (ra_deg, dec_deg)

def load_ds9_core_region(cores, filename_base = 'taurus_av_boxes_',
        header=None):

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    regions = read_ds9_region(filename_base + '.reg')

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        core = tag[tag.find('{')+1:tag.find('}')]

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

        if core in cores:
            cores[core]['poly_verts'] = {}
            cores[core]['poly_verts']['wcs'] = poly_verts
            cores[core]['poly_verts']['pixel'] = poly_verts_pix
        else:
            pass

    return cores

def load_ds9_region(props, filename=None, header=None):

    import pyregion as pyr

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]

    regions = pyr.open(filename)

    props['regions'] = {}

    for region in regions:
        # Cores defined in following format: 'tag={L1495A}'
        tag = region.comment
        region_name = tag[tag.find('text={')+6:tag.find('}')].lower()

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

        props['regions'][region_name] = {}
        props['regions'][region_name]['poly_verts'] = {}
        props['regions'][region_name]['poly_verts']['wcs'] = poly_verts
        props['regions'][region_name]['poly_verts']['pixel'] = poly_verts_pix

    return props

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

'''
The main script
'''

def main(dgr=None, vel_range=(-5, 15), vel_range_type='single', region=None,
        av_data_type='planck'):
    ''' Executes script.

    Parameters
    ----------
    dgr : float
        If None, pulls best-fit value from properties.
    vel_range : tuple
        If None, pulls best-fit value from properties.
    '''

    # import external modules
    import pyfits as fits
    import numpy as np
    from mycoords import make_velocity_axis
    import mygeometry as myg
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    import json
    from os import system,path

    # Script parameters
    # -----------------
    # Name of noise cube
    noise_cube_filename = 'multicloud_hi_galfa_cube_regrid_planckres_noise.fits'

    # Use Planck dust Av map or Kainulainen 2009 optical extinction Av map?
    # options are 'planck' or 'lee12'
    #av_data_type = 'lee12'
    #av_data_type = 'planck'

    # Global parameter file
    prop_file = 'multicloud_global_properties'

    # Which cores to include in analysis?
    cores_to_keep = [# taur
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
                     'L1527-2',
                     # Calif
                     'L1536',
                     'L1483-1',
                     'L1483-2',
                     'L1482-1',
                     'L1482-2',
                     'L1478-1',
                     'L1478-2',
                     'L1456',
                     'NGC1579',
                     #'L1545',
                     #'L1517',
                     #'L1512',
                     #'L1523',
                     #'L1512',
                     # Pers
                     'B5',
                     'IC348',
                     'B1E',
                     'B1',
                     'NGC1333',
                     'B4',
                     'B3',
                     'L1455',
                     'L1448',
                     ]

    # Regions, regions to edit the global properties with
    if region == 1:
        region_limit = {'wcs' : (((5, 10, 0), (19, 0, 0)),
                                 ((4, 30, 0), (27, 0, 0))),
                          'pixel' : ()
                         }
    elif region == 2:
        region_limit = {'wcs' : (((4, 30, 0), (19, 0, 0)),
                                 ((3, 50, 0), (29, 0, 0))),
                          'pixel' : ()
                        }
    elif region == 3:
        region_limit = {'wcs' : (((4, 30, 0), (29, 0, 0)),
                                 ((3, 50, 0), (33, 0, 0))),
                          'pixel' : ()
                        }
    else:
        region_limit = None

    # define directory locations
    # --------------------------
    output_dir = '/d/bip3/ezbc/multicloud/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/multicloud/figures/'
    av_dir = '/d/bip3/ezbc/multicloud/data/av/'
    rad_dir = '/d/bip3/ezbc/multicloud/data/radiation_field/'
    hi_dir = '/d/bip3/ezbc/multicloud/data/hi/'
    co_dir = '/d/bip3/ezbc/multicloud/data/co/'
    core_dir = '/d/bip3/ezbc/multicloud/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/multicloud/data/python_output/'
    region_dir = '/d/bip3/ezbc/multicloud/data/python_output/regions/'

    rad_field, rad_field_header = load_fits(rad_dir + \
                'multicloud_rad_field_5arcmin.fits',
            return_header=True)
    rad_field_error, av_error_header = load_fits(rad_dir + \
            'multicloud_rad_field_error_5arcmin.fits',
            return_header=True)

    # Prepare data products
    # ---------------------
    # Load global properties of cloud
    # global properties written from script
    # 'av/multicloud_analysis_global_properties.txt'
    if region is not None:
        likelihood_filename += '_region{0:.0f}'.format(region)
        results_filename += '_region{0:.0f}'.format(region)

    print('\nLoading global property file {0:s}.txt'.format(prop_file))
    with open(property_dir + prop_file + '.txt', 'r') as f:
        props = json.load(f)

    # Change WCS coords to pixel coords of images
    props = convert_limit_coordinates(props,
                                      header=rad_field_header,
                                      coords=('region_limit',
                                              'plot_limit',
                                              'region_name_pos'))

    # Load cloud division regions from ds9
    #region_filename = region_dir + \
    #        'multicloud_divisions_coldcore_selection.reg'
    region_filename = region_dir + 'multicloud_divisions.reg'
    props = load_ds9_region(props,
                            filename=region_filename,
                            header=rad_field_header)

    # Write region name pos
    props['region_name_pos'] = {
             #'taurus 1' : {'wcs' : ((3, 50,  0),
             #                       (21.5, 0, 0)),
             #             },
             #'taurus 2' : {'wcs' : ((5, 10,  0),
             #                       (21.5, 0, 0)),
             #             },
             'taurus' : {'wcs' : ((4, 40,  0),
                                  (21, 0, 0)),
                          },
             'perseus' : {'wcs' : ((3, 30,  0),
                                   (26, 0, 0)),
                          },
             #'perseus 1' : {'wcs' : ((3, 0,  0),
             #                      (34, 0, 0)),
             #             },
             #'perseus 2' : {'wcs' : ((3, 10,  0),
             #                      (22.5, 0, 0)),
             #             },
             'california' : {'wcs' : ((4, 28,  0),
                                      (34, 0, 0)),
                             },
             }


    # Derive relevant region
    pix = props['region_limit']['pixel']
    region_vertices = ((pix[1], pix[0]),
                       (pix[1], pix[2]),
                       (pix[3], pix[2]),
                       (pix[3], pix[0])
                       )

    # block offregion
    region_mask = myg.get_polygon_mask(rad_field, region_vertices)

    print('\nRegion size = ' + \
          '{0:.0f} pix'.format(region_mask[region_mask == 1].size))

    cloud_dict = {'taurus' : {},
                  'perseus' : {},
                  'california' : {},
                  }

    # Plot
    figure_types = ['png', 'pdf']
    for figure_type in figure_types:
        filename = 'multicloud_rad_field_map' + \
                   '.{0:s}'.format(figure_type)

        print('\nSaving radiation field map to \n' + filename)

        rad_field[rad_field < 0] = np.nan

        plot_rad_map(header=rad_field_header,
                       #rad_field=np.log10(rad_field),
                       rad_field=rad_field,
                       limits=props['plot_limit']['pixel'],
                       regions=props['regions'],
                       cloud_dict=cloud_dict,
                       cores_to_keep=cores_to_keep,
                       props=props,
                       #vlimits=(0.5,10),
                       vlimits=(-0.5,9),
                       vscale='linear',
                       #vlimits=(0.1,30),
                       savedir=figure_dir + 'maps/',
                       filename=filename,
                       show=False)

if __name__ == '__main__':

    main(av_data_type='dust')



