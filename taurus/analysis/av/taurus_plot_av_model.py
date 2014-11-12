#!/usr/bin/python

''' Calculates the N(HI) map for Taurus

'''

import numpy as np

''' Plotting Functions
'''

def plot_nhi_image(nhi_image=None, header=None, contour_image=None,
        av_image=None,
        cores=None, title=None, limits=None,
        contours=None, boxes=False, savedir='./', filename=None, show=True):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 15
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (15, 7),
              'figure.titlesize': font_scale
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    if av_image is not None:
        nrows_ncols=(1,2)
        ngrids=2
    else:
        nrows_ncols=(1,1)
        ngrids=1

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode="each",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=1,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # ------------------
    # NHI image
    # ------------------
    # create axes
    ax = imagegrid[0]
    cmap = cm.Greys # colormap
    # show the image
    im = ax.imshow(nhi_image,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=0,
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',)
    ax.set_ylabel('Declination (J2000)',)

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Write label to colorbar
    cb.set_label_text(r'N(HI) $\times$ 10$^{20}$ cm$^{-2}$',)

    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    if type(cores) is dict:
        for core in cores:
            pix_coords = cores[core]['center_pixel']

            anno_color = (0.3, 0.5, 1)

            ax.scatter(pix_coords[0],pix_coords[1],
                    color=anno_color,
                    s=200,
                    marker='+',
                    linewidths=2)

            ax.annotate(core,
                    xy=[pix_coords[0], pix_coords[1]],
                    xytext=(5,5),
                    textcoords='offset points',
                    color=anno_color)

            if boxes:
                rect = ax.add_patch(Polygon(
                    cores[core]['box_vertices'][:, ::-1],
                        facecolor='none',
                        edgecolor=anno_color))

    # ------------------
    # Av image
    # ------------------
    if av_image is not None:
        # create axes
        ax = imagegrid[1]
        # show the image
        im = ax.imshow(av_image,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension (J2000)',)
        ax.set_ylabel('Declination (J2000)',)

        # colorbar
        cb = ax.cax.colorbar(im)
        cmap.set_bad(color='w')
        # plot limits
        if limits is not None:
            ax.set_xlim(limits[0],limits[2])
            ax.set_ylim(limits[1],limits[3])

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='r')

        # Write label to colorbar
        cb.set_label_text(r'$A_V$ (mag)',)

        # Convert sky to pix coordinates
        wcs_header = pywcs.WCS(header)
        if type(cores) is dict:
            for core in cores:
                pix_coords = cores[core]['center_pixel']

                anno_color = (0.3, 0.5, 1)

                if boxes:
                    rect = ax.add_patch(Polygon(
                        cores[core]['box_vertices'][:, ::-1],
                            facecolor='none',
                            edgecolor=anno_color))

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        fig.show()

def plot_av_model(av_image=None, header=None, contour_image=None,
        av_model=None, hi_velocity_axis=None, vel_range=None, hi_spectrum=None,
        cores=None, results=None, title=None, limits=None, contours=None,
        boxes=False, savedir='./', filename=None, show=True):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon

    # Set up plot aesthetics
    plt.clf()
    plt.close()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 12
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': False,
              'figure.figsize': (13, 10),
              'figure.titlesize': font_scale
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    # ==========================================================================
    # Av maps
    # ==========================================================================

    nrows_ncols=(1,2)
    ngrids=2

    imagegrid = ImageGrid(fig, (2,1,1),
                 nrows_ncols=nrows_ncols,
                 ngrids=ngrids,
                 cbar_mode='single',
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0.2,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=True,
                 label_mode='L',
                 share_all=True)

    # ------------------
    # Av model
    # ------------------
    # create axes
    ax = imagegrid[0]
    cmap = cm.Greys # colormap
    # show the image
    vmax = np.max((np.max(av_image[av_image == av_image]),
                   np.max(av_model[av_model == av_model])))
    vmax = 1.4
    im = ax.imshow(av_model,
            interpolation='nearest',origin='lower',
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            #norm=matplotlib.colors.LogNorm()
            )

    # Asthetics
    ax.set_display_coord_system("fk5")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',)
    ax.set_ylabel('Declination (J2000)',)

    ax.set_title(r'Model $A_V$')

    # colorbar
    cb = ax.cax.colorbar(im)
    cmap.set_bad(color='w')
    # plot limits
    if limits is not None:
        ax.set_xlim(limits[0],limits[2])
        ax.set_ylim(limits[1],limits[3])

    # Plot Av contours
    if contour_image is not None:
        ax.contour(contour_image, levels=contours, colors='r')

    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    if type(cores) is dict:
        for core in cores:
            pix_coords = cores[core]['center_pixel']

            anno_color = (0.3, 0.5, 1)

            ax.scatter(pix_coords[0],pix_coords[1],
                    color=anno_color,
                    s=200,
                    marker='+',
                    linewidths=2)

            ax.annotate(core,
                    xy=[pix_coords[0], pix_coords[1]],
                    xytext=(5,5),
                    textcoords='offset points',
                    color=anno_color)

            if boxes:
                rect = ax.add_patch(Polygon(
                    cores[core]['box_vertices'][:, ::-1],
                        facecolor='none',
                        edgecolor=anno_color))

    if results is not None:
        av_thres = results['av_threshold']['value']
        co_thres = results['co_threshold']['value']
        if av_thres > 15:
            av_thres = None
        text = ''
        text += r'N$_{\rm pix}$ = ' + \
                 '{0:.0f}'.format(results['npix'])
        text += '\n'
        if av_thres is not None:
            text += r'$A_V$ threshold = ' + \
                    '{0:.1f} mag'.format(av_thres)
        else:
            text += r'$A_V$ threshold = ' + \
                    '{0:s}'.format(av_thres)
        text += '\n'
        text += r'CO threshold = ' + \
                '{0:.1f} K km/s'.format(co_thres)
        text += '\n'
        text += r'DGR = {0:.2f} '.format(results['dust2gas_ratio']['value']) + \
                r'$\times$ 10$^{-20}$ (cm$^2$ mag$^1$)'
        text += '\n'
        vel_range = results['hi_velocity_range'][:2]
        if vel_range.ndim == 1:
            text += r'Velocity range = ' + \
                    '{0:.1f} to {1:.1f} km/s'.format(vel_range[0], vel_range[1])
            text += '\n'
        text += r'$\chi^2$ / $\nu$ = {0:.1f}'.format(results['chisq'])
        ax.annotate(text,
                xytext=(0.03, 0.95),
                xy=(0.03, 0.95),
                textcoords='axes fraction',
                xycoords='axes fraction',
                color='k',
                fontsize=font_scale*0.75,
                bbox=dict(boxstyle='round',
                          facecolor='w',
                          alpha=0.8),
                horizontalalignment='left',
                verticalalignment='top',
                )

    # ------------------
    # Av image
    # ------------------
    if av_image is not None:
        # create axes
        ax = imagegrid[1]
        # show the image
        im = ax.imshow(av_image,
                interpolation='nearest',origin='lower',
                cmap=cmap,
                vmin=0,
                vmax=vmax,
                #norm=matplotlib.colors.LogNorm()
                )

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        ax.set_xlabel('Right Ascension (J2000)',)
        ax.set_ylabel('Declination (J2000)',)

        ax.set_title(r'Observed $A_V$')

        # colorbar
        cb = ax.cax.colorbar(im)
        cmap.set_bad(color='c')

        # plot limits
        if limits is not None:
            ax.set_xlim(limits[0],limits[2])
            ax.set_ylim(limits[1],limits[3])

        # Plot Av contours
        if contour_image is not None:
            ax.contour(contour_image, levels=contours, colors='r')

        # Write label to colorbar
        cb.set_label_text(r'$A_V$ (mag)',)

    # ==========================================================================
    # HI Spectrum
    # ==========================================================================

    # create axes
    ax = fig.add_subplot(2,1,2)

    ax.plot(hi_velocity_axis,
            hi_spectrum,
            color='k',
            drawstyle = 'steps-mid'
            )

    # Plot velocity range
    if vel_range.ndim == 1:
    	ax.axvspan(vel_range[0], vel_range[1], color='k', alpha=0.3)
    elif vel_range.ndim == 2:
        for i in xrange(0, vel_range.shape[0]):
            ax.axvspan(vel_range[i, 0], vel_range[i, 1], color='k', alpha=0.3)

    ax.set_xlabel('Velocity (km/s)')
    ax.set_ylabel(r'T$_b$ (K)')

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        plt.show()

''' DS9 Region and Coordinate Functions
'''

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

    return region[0].coord_list

def load_ds9_region(cores, filename_base = 'taurus_av_boxes_', header=None):

    # region[0] in following format:
    # [64.26975, 29.342033333333333, 1.6262027777777777, 3.32575, 130.0]
    # [ra center, dec center, width, height, rotation angle]
    for core in cores:
        region = read_ds9_region(filename_base + core + '.reg')
        if region is not None:
            box_center_pixel = get_pix_coords(ra = region[0],
                                              dec = region[1],
                                              header = header)
            box_center_pixel = (int(box_center_pixel[1]),
                    int(box_center_pixel[0]))
            box_height = region[2] / header['CDELT1']
            box_width = region[3] / header['CDELT2']
            cores[core].update({'box_center_pix': box_center_pixel})
            cores[core].update({'box_width': box_width})
            cores[core].update({'box_height': box_height})
            cores[core].update({'box_angle': region[4]})

    return cores

'''
The main script
'''

def main(dgr=0.1, vel_range=(0,15), vel_range_type='single', region=None):
    ''' Executes script.
    '''

    # import external modules
    import pyfits as pf
    import numpy as np
    from mycoords import make_velocity_axis
    import mygeometry as myg
    from myimage_analysis import calculate_nhi, calculate_noise_cube, \
        calculate_sd, calculate_nh2, calculate_nh2_error
    import json

    # Script parameters
    # -----------------
    # Name of noise cube
    noise_cube_filename = 'taurus_hi_galfa_cube_regrid_planckres_noise.fits'

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
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/maps/av_models/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/hi/'
    co_dir = '/d/bip3/ezbc/taurus/data/co/'
    core_dir = '/d/bip3/ezbc/taurus/data/python_output/core_properties/'
    property_dir = '/d/bip3/ezbc/taurus/data/python_output/'
    region_dir = '/d/bip3/ezbc/taurus/data/python_output/ds9_regions/'

    # load Planck Av and GALFA HI images, on same grid
    av_image, av_header = load_fits(av_dir + \
                'taurus_av_planck_5arcmin.fits',
            return_header=True)

    av_error_image, av_error_header = load_fits(av_dir + \
                'taurus_av_error_planck_5arcmin.fits',
            return_header=True)

    hi_cube, hi_header = load_fits(hi_dir + \
                'taurus_hi_galfa_cube_regrid_planckres.fits',
            return_header=True)

    hi_noise_cube, noise_header = load_fits(hi_dir + noise_cube_filename,
            return_header=True)

    co_data, co_header = load_fits(co_dir + \
                'taurus_co_cfa_cube_regrid_planckres.fits',
            return_header=True)

    # Load global properties of cloud
    # global properties written from script
    # 'av/taurus_analysis_global_properties.txt'
    with open(property_dir + 'taurus_global_properties_region' + str(region) + \
              '.txt', 'r') as f:
        props = json.load(f)
        '''
        dgr = props['dust2gas_ratio_max']['value']
        vel_center = props['hi_velocity_center_max']['value']
        vel_width = props['hi_velocity_center_max']['value']
        vel_range = np.asarray(props['hi_velocity_range_max']['value'])
        dgr = dgr/3
        scale = 1.
        vel_range = (vel_center - scale*vel_width/2.0,
                     vel_center + scale*vel_width/2.0)
        vel_range = (0, 15)
        '''

    props['hi_velocity_range'] = vel_range
    props['dust2gas_ratio']['value'] = dgr

    # define core properties
    with open(core_dir + 'taurus_core_properties.txt', 'r') as f:
        cores = json.load(f)

    # make velocity axis for hi cube
    velocity_axis = make_velocity_axis(hi_header)

    # Write core coordinates in pixels
    cores = convert_core_coordinates(cores, hi_header)

    cores = load_ds9_region(cores,
            filename_base = region_dir + 'taurus_av_boxes_',
            header = hi_header)

    # create nhi image
    nhi_image, nhi_image_error = calculate_nhi(cube=hi_cube,
            velocity_axis=velocity_axis,
            velocity_range=vel_range,
            header=hi_header,
            noise_cube=hi_noise_cube)

    # create model av map
    av_model = nhi_image * dgr

    # Mask the images based on av trheshol
    co_data_nonans = np.copy(co_data)
    co_data_nonans[np.isnan(co_data_nonans)] = 0.0
    co_mom0 = np.sum(co_data_nonans, axis=0)
    mask = ((av_image > props['av_threshold']['value']) & \
            (co_mom0 > props['co_threshold']['value']))

    # Derive relevant region
    pix = props['region_limit']['pixel']
    region_vertices = ((pix[1], pix[0]),
                       (pix[1], pix[2]),
                       (pix[3], pix[2]),
                       (pix[3], pix[0])
                       )

    # block offregion
    region_mask = myg.get_polygon_mask(av_image, region_vertices)

    print('\nRegion size = ' + \
          '{0:.0f} pix'.format(region_mask[region_mask == 1].size))

    if vel_range_type == 'single':
        print('\nHI velocity integration range:')
        print('%.1f to %.1f km/s' % (vel_range[0],
                                     vel_range[1]))
    elif vel_range_type == 'multiple':
        print('\nHI velocity integration ranges:')
        for i in xrange(0, vel_range.shape[0]):
            print('%.1f to %.1f km/s' % (vel_range[i, 0],
                                         vel_range[i, 1]))

    print('\nDGR:')
    print('%.2f x 10^-20 cm^2 mag' % (dgr))

    # Get mask and mask images
    mask = np.asarray(props['mask'])

    av_image_masked = np.copy(av_image)
    #av_image_masked[(mask == 1) & (region_mask == 1)] = np.nan
    av_image_masked[mask == 1] = np.nan

    av_model_masked = np.copy(av_model)
    #av_model_masked[(mask == 1) & (region_mask == 1)] = np.nan
    av_model_masked[mask == 1] = np.nan

    indices = ((np.isnan(av_model_masked)) & \
               (np.isnan(av_image_masked)) & \
               (np.isnan(av_error_image)))

    print 'npix = ', mask[~mask].size
    print 'npix = ', props['npix']

    # import matplotlib.pyplot as plt
    # av_plot_data = np.copy(av_image)
    # av_plot_data[~indices] = np.nan
    # plt.imshow(av_plot_data, origin='lower')
    # plt.show()

    # Calc chi^2
    chisq = np.sum((av_image[~mask] - av_model[~mask])**2 / \
            av_error_image[~mask]**2) / props['npix']
    props['chisq'] = chisq

    # Create HI spectrum
    #hi_cube_copy = np.copy(hi_cube)
    hi_cube[:, mask==1] = 0
    #hi_spectrum = np.mean(hi_cube_copy, axis=(1,2))
    hi_spectrum = np.mean(hi_cube, axis=(1,2))

    # Plot
    figure_types = ['png',]
    for figure_type in figure_types:
        if vel_range_type == 'single':
            filename = 'single_vel_range/taurus_av_model_map_' + \
                    'dgr{0:.3f}_'.format(dgr) + \
                    '{0:.1f}to{1:.1f}kms'.format(vel_range[0], vel_range[1]) + \
                    '.%s' % figure_type
        elif vel_range_type == 'multiple':
            filename = 'multiple_vel_range/taurus_av_model_map_' + \
                       'dgr{0:.3f}'.format(dgr)
            for i in xrange(0, vel_range.shape[0]):
                filename += '_{0:.1f}to{1:.1f}kms'.format(vel_range[i, 0],
                                                          vel_range[i, 1])
            filename += '.%s' % figure_type
        plot_av_model(av_image=av_image_masked,
                      av_model=av_model_masked,
                      header=av_header,
                      results=props,
                      hi_velocity_axis=velocity_axis,
                      vel_range=vel_range,
                      hi_spectrum=hi_spectrum,
                      limits=props['region_limit']['pixel'],
                      savedir=figure_dir,
                      filename=filename,
                      show=False)

if __name__ == '__main__':
    dgrs = np.arange(0.05, 0.4, 0.025)
    vel_widths = (4, 6, 8, 10, 12, 18, 34, 50)
    vel_widths = (50,)
    vel_center = 6.15

    regions = (1,2,3)

    for region in regions:
        main(region=region)

    '''
    # Single velocity ranges
    for i in xrange(0, len(dgrs)):
        for j in xrange(0, len(vel_widths)):
            vel_range = (vel_center - vel_widths[j]/2.0,
                         vel_center + vel_widths[j]/2.0)
            vel_range = np.asarray(vel_range)
            main(dgr=dgrs[i], vel_range=vel_range,
                    vel_range_type='single')

    # Multiple velocity ranges
    for i in xrange(0, len(dgrs)):
        for j in xrange(0, len(vel_widths)):
            vel_range = ((-10, vel_center - vel_widths[j]/2.0),
                         (vel_center + vel_widths[j]/2.0, 20),)
            vel_range = np.asarray(vel_range)
            main(dgr=dgrs[i], vel_range=vel_range,
                    vel_range_type='multiple')


    '''



