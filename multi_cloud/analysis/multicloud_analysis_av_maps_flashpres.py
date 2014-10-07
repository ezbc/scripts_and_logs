#!/usr/bin/python

'''

Plots Av image for all three clouds.

'''

import pyfits as pf
import numpy as np

''' Plotting Functions
'''

def plot_av_images(cloud_dict, title=None, boxes=False, savedir='./',
        filename=None, show=True, rh2_limits=None, hi_sd_limits=None):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2
    import pywcs
    from pylab import cm # colormaps
    from matplotlib.patches import Polygon
    from matplotlib.gridspec import GridSpec

    # Set up plot aesthetics
    plt.clf()
    plt.rcdefaults()
    colormap = plt.cm.gist_ncar
    #color_cycle = [colormap(i) for i in np.linspace(0, 0.9, len(flux_list))]
    font_scale = 20
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'ytick.labelsize': font_scale,
              'font.weight': 500,
              'axes.labelweight': 500,
              'text.usetex': True,
              'figure.figsize': (20, 10),
              'figure.titlesize': font_scale
              #'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    fig = pywcsgrid2.plt.figure()

    fig.subplots_adjust(hspace=0.1)

    gs1 = GridSpec(1, 2)
    #gs1.update(left=0.02, right=0.6, bottom=0.02, top=0.98, wspace=0.00,
    #        hspace=0.12)

    #grid_helper = pywcsgrid2.GridHelper(wcs=cloud_dict[cloud]['av_header'])

    grid_pos = 0
    cloud_list = []

    cloud_dict = {'taurus':cloud_dict['taurus'],
                  'california':cloud_dict['california']}

    for i, cloud in enumerate(cloud_dict):
        av_image = cloud_dict[cloud]['av_data']
        av_header = cloud_dict[cloud]['av_header']

        ax = pywcsgrid2.subplot(gs1[grid_pos],
                header=av_header)

        # create axes
        cmap = cm.pink # colormap
        cmap.set_bad(color='w')
        # show the image
        vmin, vmax = cloud_dict[cloud]['color_scale_limits']
        im = ax.imshow(av_image,
                interpolation='nearest',
                origin='lower',
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                #norm=matplotlib.colors.LogNorm()
                )

        # Set limits
        limits = cloud_dict[cloud]['limit_pixels']
        ax.set_xlim([limits[0][0], limits[1][0]])
        ax.set_ylim([limits[0][1], limits[1][1]])

        # Asthetics
        ax.set_display_coord_system("fk5")
        ax.set_ticklabel_type("hms", "dms")

        if i != 0:
            ax.set_xlabel('')
        else:
            ax.set_xlabel('Right Ascension (J2000)',)
        ax.set_xlabel('Right Ascension (J2000)',)
        ax.set_ylabel('Declination (J2000)',)

        ax.annotate(cloud.capitalize(),
                xy=[0.02, 0.9],
                textcoords='axes fraction',
                xycoords='axes fraction',
                fontsize=font_scale * 1.5,
                color='w')

        # Convert sky to pix coordinates
        cores = cloud_dict[cloud]['cores']

        '''
        # Colorbar
        if cloud == 'perseus':
            cb = fig.colorbar(im, ax=ax)

            # Write label to colorbar
            cb.set_label(r'A$_V$ (Mag)',)
        else:
            cb = fig.colorbar(im, ax=ax)

            # Write label to colorbar
            cb.set_label(r'A$_V$ (Mag)',)
        '''

        for core in cores:
            if core in cloud_dict[cloud]['plot_cores']:
            	scale = 1.1
            	linewidth=3
            	linestyle='solid'
            else:
            	scale = 0.8
            	linewidth=1
            	linestyle='dashed'

            pix_coords = cores[core]['center_pixel']

            anno_color = (0.3, 0.5, 1)

            ax.scatter(pix_coords[0],pix_coords[1],
                    color=anno_color,
                    s=200,
                    marker='+',
                    linewidths=2)

            ax.annotate(core,
                    xy=[pix_coords[0], pix_coords[1]],
                    xytext=(5,10),
                    textcoords='offset points',
                    fontsize=font_scale*scale,
                    #arrowprops=dict(facecolor='w'),
                    color='w')

            if boxes:
                vertices = np.copy(cores[core]['box_vertices_rotated'])
                #[:, ::-1]
                rect = ax.add_patch(Polygon(
                        vertices[:, ::-1],
                        facecolor='none',
                        edgecolor=anno_color,
                        linewidth=linewidth,
                        linestyle=linestyle))

        cloud_list.append(cloud)
        grid_pos += 1

    #plt.tight_layout()

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(savedir + filename, bbox_inches='tight')
    if show:
        fig.show()

def plot_rh2_hi_vs_h(cloud_dict, fig, font_scale=10, rh2_limits=None,
        hi_sd_limits=None, cloud_list=None):

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib
    from mpl_toolkits.axes_grid1 import ImageGrid
    from myscience.krumholz09 import calc_T_cnm
    from matplotlib.gridspec import GridSpec

    nrows_ncols=(2,2)
    ngrids = 4

    grid_pos = 0

    for i, cloud in enumerate(cloud_list):
    	'''
        imagegrid = ImageGrid(fig, (3, 2, grid_pos),
                     nrows_ncols=nrows_ncols,
                     ngrids=ngrids,
                     cbar_location='right',
                     cbar_pad="2%",
                     cbar_size='3%',
                     axes_pad=0.0,
                     aspect='auto',
                     label_mode='L',
                     share_all=True)
        '''

        gs1 = GridSpec(2,2)
        gs1.update(left=0.63, right=0.98, bottom=0.69, top=1, wspace=0,
                hspace=0)
        gs2 = GridSpec(2,2)
        gs2.update(left=0.63, right=0.98, bottom=0.36, top=0.63, wspace=0,
                hspace=0)
        gs3 = GridSpec(2,2)
        gs3.update(left=0.63, right=0.98, bottom=0.03, top=0.30, wspace=0,
                hspace=0)

        # Plot rh2 and hsd for each cloud
        core_pos = 0
        for core in cloud_dict[cloud]['cores']:
            rh2 = cloud_dict[cloud]['cores'][core]['rh2']
            rh2_error = cloud_dict[cloud]['cores'][core]['rh2_error']
            h_sd = cloud_dict[cloud]['cores'][core]['h_sd']
            h_sd_error = cloud_dict[cloud]['cores'][core]['h_sd_error']
            hi_sd = cloud_dict[cloud]['cores'][core]['hi_sd']
            hi_sd_error = cloud_dict[cloud]['cores'][core]['hi_sd_error']
            rh2_fit = cloud_dict[cloud]['cores'][core]['rh2_fit']
            hi_sd_fit = cloud_dict[cloud]['cores'][core]['hi_sd_fit']
            h_sd_fit = cloud_dict[cloud]['cores'][core]['h_sd_fit']
            phi_cnm = cloud_dict[cloud]['cores'][core]['phi_cnm']
            phi_cnm_error = cloud_dict[cloud]['cores'][core]['phi_cnm_error']
            T_cnm = cloud_dict[cloud]['cores'][core]['T_cnm']
            T_cnm_error = cloud_dict[cloud]['cores'][core]['T_cnm_error']
            Z = cloud_dict[cloud]['cores'][core]['Z']
            Z_error = cloud_dict[cloud]['cores'][core]['Z_error']
            phi_mol = cloud_dict[cloud]['cores'][core]['phi_mol']
            phi_mol_error = cloud_dict[cloud]['cores'][core]['phi_mol_error']
            hi_sd_limits = cloud_dict[cloud]['hi_sd_limits']

            if core in cloud_dict[cloud]['plot_cores']:
                # Create plot
                if grid_pos == 0:
                    ax = fig.add_subplot(gs1[0:1, core_pos:core_pos+1])
                elif grid_pos == 1:
                    ax = fig.add_subplot(gs2[0:1, core_pos:core_pos+1])
                elif grid_pos == 2:
                    ax = fig.add_subplot(gs3[0:1, core_pos:core_pos+1])

                # Hide tick labels
                show_xylabel = [True, True]
                if core_pos == 1:
                    plt.setp(ax.get_yticklabels(), visible=False)
                    show_xylabel[1] = False
                plt.setp(ax.get_xticklabels(), visible=False)
                show_xylabel[0] = False

                plot_rh2_vs_h(ax,
                              rh2=np.asarray(rh2),
                              h_sd=np.asarray(h_sd),
                              rh2_error=np.asarray(rh2_error),
                              h_sd_error=np.asarray(h_sd_error),
                              rh2_fit=np.asarray(rh2_fit),
                              h_sd_fit=np.asarray(h_sd_fit),
                              core=core,
                              phi_cnm=phi_cnm,
                              phi_cnm_error=phi_cnm_error,
                              phi_mol=phi_mol,
                              phi_mol_error=phi_mol_error,
                              T_cnm=T_cnm,
                              T_cnm_error=T_cnm_error,
                              Z=Z,
                              Z_error=Z_error,
                              font_scale=font_scale,
                              rh2_limits=rh2_limits,
                              show_xylabel=show_xylabel,
                              )


                if grid_pos == 0:
                    ax = fig.add_subplot(gs1[1:2, core_pos:core_pos+1])
                elif grid_pos == 1:
                    ax = fig.add_subplot(gs2[1:2, core_pos:core_pos+1])
                elif grid_pos == 2:
                    ax = fig.add_subplot(gs3[1:2, core_pos:core_pos+1])

                # Hide tick labels
                show_xylabel = [True, True]
                if core_pos == 1:
                    plt.setp(ax.get_yticklabels(), visible=False)
                    show_xylabel[1] = False

                plot_hi_vs_h(ax,
                             hi_sd=np.asarray(hi_sd),
                             h_sd=np.asarray(h_sd),
                             hi_sd_error=np.asarray(hi_sd_error),
                             h_sd_error=np.asarray(h_sd_error),
                             hi_sd_fit=np.asarray(hi_sd_fit),
                             h_sd_fit=np.asarray(h_sd_fit),
                             core=core,
                             phi_cnm=phi_cnm,
                             phi_cnm_error=phi_cnm_error,
                             phi_mol=phi_mol,
                             phi_mol_error=phi_mol_error,
                             T_cnm=T_cnm,
                             T_cnm_error=T_cnm_error,
                             Z=Z,
                             Z_error=Z_error,
                             font_scale=font_scale,
                             hi_sd_limits=hi_sd_limits,
                             show_xylabel=show_xylabel,
                             )


                core_pos += 1

        grid_pos += 1

def plot_rh2_vs_h(ax, rh2_fit=None, h_sd_fit=None, rh2_error=None,
        h_sd_error=None, h_sd=None, rh2=None, core=None, phi_cnm=None,
        phi_cnm_error=None, phi_mol=None, phi_mol_error=None, T_cnm=None,
        T_cnm_error=None, Z=None, Z_error=None, font_scale=10,
        rh2_limits=None, show_xylabel=[1, 1]):

    # Drop the NaNs from the images
    if type(rh2_error) is float:
        indices = np.where((rh2 == rh2) &\
                           (h_sd == h_sd)&\
                           (h_sd > 0) &\
                           (rh2 > 0))

    if type(rh2_error) is np.ndarray or \
            type(rh2_error) is np.ma.core.MaskedArray or \
            type(h_sd_error) is np.ndarray or \
            type(h_sd_error) is np.ma.core.MaskedArray:
        indices = np.where((rh2 == rh2) &\
                           (h_sd == h_sd) &\
                           (h_sd_error == h_sd_error) &\
                           (rh2_error == rh2_error) &\
                           (h_sd > 0) &\
                           (rh2 > 0))

    rh2_nonans = rh2[indices]
    h_sd_nonans = h_sd[indices]

    if type(rh2_error) is np.ndarray:
        rh2_error_nonans = rh2_error[indices]
    else:
        rh2_error_nonans = np.array(rh2_error[indices])

    if type(h_sd_error) is np.ndarray or \
            type(h_sd_error) is np.ma.core.MaskedArray:
        h_sd_error_nonans = h_sd_error[indices]
    else:
        h_sd_error_nonans = h_sd_error * \
                np.ones(h_sd[indices].shape)

    image = ax.errorbar(h_sd_nonans.ravel(),
            rh2_nonans.ravel(),
            xerr=(h_sd_error_nonans.ravel()),
            yerr=(rh2_error_nonans.ravel()),
            alpha=0.75,
            color='k',
            marker='^',ecolor='k',linestyle='none',
            markersize=4
            )

    if rh2_fit is not None:
        ax.plot(h_sd_fit, rh2_fit,
                color = 'r', alpha=0.5)

    # Annotations
    anno_xpos = 0.95

    phi_cnm_text = r'\noindent$\phi_{\rm CNM}$ =' + \
                   r' %.2f' % (phi_cnm) + \
                   r'$^{+%.2f}_{-%.2f}$ \\' % (phi_cnm_error[0],
                                             phi_cnm_error[1])
    T_cnm_text = r'\noindent T$_{\rm CNM}$ =' + \
                 r' %.2f' % (T_cnm) + \
                 r'$^{+%.2f}_{-%.2f}$ \\' % (T_cnm_error[0],
                                             T_cnm_error[1])
    if Z_error == [0.0, 0.0] or Z_error == (0.0, 0.0):
        Z_text = r'Z = %.1f Z$_\odot$ \\' % (Z)
        Z_text = ''
    else:
        Z_text = r'Z = %.2f' % (Z) + \
        r'$^{+%.2f}_{-%.2f}$ Z$_\odot$ \\' % (Z_error[0],
                                              Z_error[1])
    if phi_mol_error == [0.0, 0.0] or phi_mol_error == (0.0, 0.0):
        phi_mol_text = r'\noindent $\phi_{\rm mol}$ = ' + \
                         '%.1f' % (phi_mol) + r''
        phi_mol_text = ''
    else:
        phi_mol_text = r'\noindent$\phi_{\rm mol}$ =' + \
                    r' %.2f' % (phi_mol) + \
                    r'$^{+%.2f}_{-%.2f}$ \\' % (phi_mol_error[0],
                                             phi_mol_error[1]) + \
                    r''

    ax.annotate(phi_cnm_text + T_cnm_text + Z_text + phi_mol_text,
            xytext=(anno_xpos, 0.05),
            xy=(anno_xpos, 0.05),
            textcoords='axes fraction',
            xycoords='axes fraction',
            size=font_scale*3/4.0,
            color='k',
            bbox=dict(boxstyle='round',
                      facecolor='w',
                      alpha=1),
            horizontalalignment='right',
            verticalalignment='bottom',
            )

    ax.set_xscale('linear', nonposx = 'clip')
    ax.set_yscale('log', nonposy = 'clip')

    if rh2_limits is not None:
        ax.set_xlim(rh2_limits[0],rh2_limits[1])
        ax.set_ylim(rh2_limits[2],rh2_limits[3])

    # Adjust asthetics
    if show_xylabel[0]:
        ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{H2}$ ' + \
                      '(M$_{\odot}$ / pc$^{2}$)',)
    if show_xylabel[1]:
        ax.set_ylabel(r'R$_{H2}$',)
    ax.set_title(core)
    ax.grid(False)

def plot_hi_vs_h(ax, hi_sd_fit=None, h_sd_fit=None, hi_sd_error=None,
        h_sd_error=None, h_sd=None, hi_sd=None, core=None, phi_cnm=None,
        phi_cnm_error=None, phi_mol=None, phi_mol_error=None, T_cnm=None,
        T_cnm_error=None, Z=None, Z_error=None, font_scale=10,
        hi_sd_limits=None, show_xylabel=[1, 1]):

    # Drop the NaNs from the images
    if type(hi_sd_error) is float:
        indices = np.where((hi_sd == hi_sd) &\
                           (h_sd == h_sd)&\
                           (h_sd > 0) &\
                           (hi_sd > 0))

    if type(hi_sd_error) is np.ndarray or \
            type(hi_sd_error) is np.ma.core.MaskedArray or \
            type(h_sd_error) is np.ndarray or \
            type(h_sd_error) is np.ma.core.MaskedArray:
        indices = np.where((hi_sd == hi_sd) &\
                           (h_sd == h_sd) &\
                           (h_sd_error == h_sd_error) &\
                           (hi_sd_error == hi_sd_error) &\
                           (h_sd > 0) &\
                           (hi_sd > 0))

    hi_sd_nonans = hi_sd[indices]
    h_sd_nonans = h_sd[indices]

    if type(hi_sd_error) is np.ndarray:
        hi_sd_error_nonans = hi_sd_error[indices]
    else:
        hi_sd_error_nonans = np.array(hi_sd_error[indices])

    if type(h_sd_error) is np.ndarray or \
            type(h_sd_error) is np.ma.core.MaskedArray:
        h_sd_error_nonans = h_sd_error[indices]
    else:
        h_sd_error_nonans = h_sd_error * \
                np.ones(h_sd[indices].shape)

    image = ax.errorbar(h_sd_nonans.ravel(),
            hi_sd_nonans.ravel(),
            xerr=(h_sd_error_nonans.ravel()),
            yerr=(hi_sd_error_nonans.ravel()),
            alpha=0.3,
            color='k',
            marker='^',ecolor='k',linestyle='none',
            markersize=3
            )

    if hi_sd_fit is not None:
        ax.plot(h_sd_fit, hi_sd_fit,
                color = 'r')

    ax.set_xscale('linear')
    ax.set_yscale('linear')

    if hi_sd_limits is not None:
        ax.set_xlim(hi_sd_limits[0], hi_sd_limits[1])
        ax.set_ylim(hi_sd_limits[2], hi_sd_limits[3])

    # Adjust asthetics
    if show_xylabel[0]:
        ax.set_xlabel('$\Sigma_{HI}$ + $\Sigma_{H2}$ (M$_\odot$ / pc$^2$)',)
    if show_xylabel[1]:
        ax.set_ylabel(r'$\Sigma_{HI}$ (M$_\odot$ / pc$^2$)',)

''' Calculations
'''

def create_box(core_pos, width, height, core_rel_pos=0.1):

    '''
    Parameters
    ----------
    core_pos : array-like
        x and y pixel coordinates of core
    width : int, float
        Width of box along x axis.
    height : int, float
        Height of box along y axis.
    core_rel_pos : float, optional
        Core position in box along y-axis as a fraction of the height with the
        origin in the south.

    Returns
    -------
    box_vertices : numpy.array
        4 x 2 array with box pixel vertices, starting from lower-left going
        clockwise.

    '''

    box_vertices = np.empty((4, 2))

    # x-coords
    box_vertices[:2, 0] = core_pos[0] - width / 2.
    box_vertices[2:4, 0] = core_pos[0] + width / 2.

    # y-coords
    offset = height * core_rel_pos
    box_vertices[1:3, 1] = core_pos[1] + height - offset
    box_vertices[(0, 3), 1] = core_pos[1] - offset

    return box_vertices

def rotate_box(box_vertices, anchor, angle):

    '''
    Parameters
    ----------
    box_vertices : numpy.array
        4 x 2 array with box pixel vertices, starting from lower-left going
        clockwise.
    anchor : tuple
        x and y coordinates of pivot point.
    angle : float
        Angle to rotate polygon clockwise from North.

    Returns
    -------
    box_vertices_rotated : numpy.array
        4 x 2 array with rotated box pixel vertices.
    '''

    from mygeometry import rotate_polygon

    box_vertices_rotated = rotate_polygon(box_vertices, anchor, angle)

    return box_vertices_rotated

def get_radial_profile(image, center=None, stddev=False, binsize=1,
        mask=None, weights=None):

    ''' Calculates radial profiles of an image at the center.

    '''

    # import external modules
    import numpy as np
    from scipy.optimize import curve_fit
    from agpy import azimuthalAverage as radial_average

    if stddev and weights is not None:
        weights=None

    result = radial_average(image, binsize=binsize, center=center,
            stddev=stddev, mask=mask, interpnan=False, returnradii=True,
            weights=weights)

    return result

def fit_profile(radii, profile, function, sigma=None):

    ''' Fits a radial profile with a power function A * radius**alpha where A
    and alpha are the fitted constants.
    '''

    # import external modules
    from scipy.optimize import curve_fit
    import numpy as np

    profile_fit = curve_fit(function, radii, profile, sigma=sigma,
            maxfev=1000000,)

    return profile_fit

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
        cores[core]['center_pixel'] = get_pix_coords(ra=center_wcs[0],
                                                     dec=center_wcs[1],
                                                     header=header)[:2].tolist()
        # convert box corners to pixel coords
        #for i in range(len(box_wcs)/2):
        #    pixels = get_pix_coords(ra=box_wcs[2*i], dec=box_wcs[2*i + 1],
        #            header=header)
        #    box_pixel[2*i], box_pixel[2*i + 1] = int(pixels[0]), int(pixels[1])
        #cores[core]['box_pixel'] = box_pixel

    return cores

def convert_limit_coordinates(clouds):

    for cloud in clouds:
        clouds[cloud].update({'limit_pixels': []})

        header = clouds[cloud]['av_header']

        limit_wcs = clouds[cloud]['limit_wcs']

        for limits in limit_wcs:
            # convert centers to pixel coords
            limit_pixels = get_pix_coords(ra=limits[0],
                                         dec=limits[1],
                                         header=header)[:2].tolist()

            clouds[cloud]['limit_pixels'].append(limit_pixels)

    return clouds

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

def main():

    import grid
    import numpy as np
    from os import system,path
    import mygeometry as myg
    from mycoords import make_velocity_axis
    import json
    #from astropy.io import fits
    import pyfits as fits

    # define directory locations
    figure_dir = '/d/bip3/ezbc/multicloud/figures/maps/'

    cloud_dict = {'taurus' : {},
                  'perseus' : {},
                  'california' : {},
                  }

    # load Planck Av and GALFA HI images, on same grid
    for cloud in cloud_dict:
        file_dir = '/d/bip3/ezbc/{0:s}/data/av/'.format(cloud)
        av_data, av_header = fits.getdata(file_dir + \
                    '{0:s}_av_planck_5arcmin.fits'.format(cloud),
                    header=True)
        av_error_data, av_error_header = fits.getdata(file_dir + \
                '{0:s}_av_error_planck_5arcmin.fits'.format(cloud),
                header=True)

        hi_dir = '/d/bip3/ezbc/{0:s}/data/hi/'.format(cloud)
        hi_data, hi_header = fits.getdata(hi_dir + \
                '{0:s}_hi_galfa_cube_regrid_planckres.fits'.format(cloud),
                header=True)

        cloud_dict[cloud]['av_data'] = av_data
        cloud_dict[cloud]['av_header'] = av_header
        cloud_dict[cloud]['av_error_data'] = av_error_data
        cloud_dict[cloud]['av_error_header'] = av_error_header
        cloud_dict[cloud]['hi_data'] = hi_data
        cloud_dict[cloud]['hi_header'] = hi_header

        # define core properties
        with open('/d/bip3/ezbc/{0:s}/data/python_output/'.format(cloud) + \
                  'core_properties/{0:s}_core_properties.txt'.format(cloud),
                  'r') as f:
             cores = json.load(f)

        cores = convert_core_coordinates(cores, av_header)

        cores = load_ds9_region(cores,
                filename_base = '/d/bip3/ezbc/{0:s}/data/'.format(cloud) + \
                                'python_output/ds9_regions/taurus_av_boxes_',
                header = av_header)

        cloud_dict[cloud]['cores'] = cores

    # Limits of cloud map in RA and Dec
    cloud_dict['taurus']['limit_wcs'] = (((4, 51, 0), (21, 54, 0)),
                                         ((4, 5, 0), (30, 25, 0)))
    cloud_dict['perseus']['limit_wcs'] = (((3, 58, 0), (27, 6, 0)),
                                          ((3, 20, 0), (35, 0, 0)))
    cloud_dict['california']['limit_wcs'] = (((4, 49, 0), (31, 54, 0)),
                                             ((3, 40, 0), (42, 24, 0)))

    # Which cores will be plotted for RH2 and HI surf dens?
    cloud_dict['taurus']['plot_cores'] = ('B213', 'L1495')
    cloud_dict['perseus']['plot_cores'] = ('B1', 'B5')
    cloud_dict['california']['plot_cores'] = ('L1482', 'L1456')

    # color scaling of Av images
    cloud_dict['taurus']['color_scale_limits'] = (1, 25)
    cloud_dict['perseus']['color_scale_limits'] = (0.1, 25)
    cloud_dict['california']['color_scale_limits'] = (0, 15)
    cloud_dict['taurus']['hi_sd_limits'] = [0, 79, 1, 5]
    cloud_dict['perseus']['hi_sd_limits'] = [0, 79, 3, 8.9]
    cloud_dict['california']['hi_sd_limits'] = [0, 79, 4, 8]

    cloud_dict = convert_limit_coordinates(cloud_dict)

    # Plot
    figure_types = ['png',]# 'pdf']
    for figure_type in figure_types:
        plot_av_images(cloud_dict=cloud_dict,
                       boxes=True,
                       #title=r'taurus: A$_V$ map with core boxed-regions.',
                       savedir=figure_dir,
                       filename='multicloud_av_cores_map_flashpres' + \
                                '.{0:s}'.format(figure_type),
                       show=0,
                       rh2_limits=[0, 79, 1.1*10**-3, 0.5*10**2],
                       #rh2_limits=[0.1, 99, 1.1*10**-3, 3*10**2],
                       #hi_sd_limits=[0.1, 99, 1, 12],
                       )

if __name__ == '__main__':
    main()



