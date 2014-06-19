#!/usr/bin/python

''' Calculates the N(HI) / Av correlation for the Taurus molecular cloud.

'''

def plot_core(core_image=None, header=None, savedir='./', filename=None,
        show=True):

    # Import external modules
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.axes_grid1 import ImageGrid
    import pyfits as pf
    import matplotlib.pyplot as plt
    import pywcsgrid2 as wcs
    import pywcs
    from pylab import cm # colormaps

    fig = plt.figure(figsize=(8,8))

    imagegrid = ImageGrid(fig, (1,1,1),
                 nrows_ncols=(1,1),
                 ngrids=1,
                 cbar_mode="single",
                 cbar_location='right',
                 cbar_pad="2%",
                 cbar_size='3%',
                 axes_pad=0,
                 axes_class=(wcs.Axes,
                             dict(header=header)),
                 aspect=False,
                 label_mode='L',
                 share_all=False)

    ax = imagegrid[0]
    im = ax.imshow(core_image,
            interpolation='nearest',origin='lower',
            cmap=cm.gray)

    # Asthetics
    ax.set_display_coord_system("fk4")
    ax.set_ticklabel_type("hms", "dms")

    ax.set_xlabel('Right Ascension (J2000)',
              size = 'small',
              family='serif')
    ax.set_ylabel('Declination (J2000)',
              size = 'small',
              family='serif')
    # colorbar
    cb = ax.cax.colorbar(im)
    # plot limits
    ax.set_xlim([112,296])
    ax.set_ylim([100,257])

    # Write label to colorbar
    cb.set_label_text(r'N(HI) $\times$ 10$^{20}$ cm$^{-2}$',
                   size='small',
                   family='serif')

    # Convert sky to pix coordinates
    wcs_header = pywcs.WCS(header)
    for key in cores:
        pix_coords = wcs_header.wcs_sky2pix([cores[key]['wcs_position']], 0)[0]

        ax.scatter(pix_coords[0],pix_coords[1], color='r', s=200, marker='+',
                linewidths=2)

        ax.annotate(key,
                xy=[pix_coords[0], pix_coords[1]],
                xytext=(5,5),
                textcoords='offset points',
                color='r')

    if filename is not None:
        plt.savefig(savedir + filename)
    if show:
        fig.show()

def main():

    import grid
    import numpy as np
    from os import system,path
    import myclumpfinder as clump_finder
    reload(clump_finder)

    # define directory locations
    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'

    # load 2mass Av and GALFA HI images, on same grid
    av_data, av_header = load_fits(av_dir + 'taurus_av_k09_regrid.fits',
            return_header=True)

    cores = {'L1495':
                {'wcs_position': [15*(4+14/60.), 28+11/60., 0],
                 'map': None,
                 'threshold': 4.75},
             'L1495A':
                {'wcs_position': [15*(4+18/60.), 28+23/60., 0],
                 'map': None,
                 'threshold': 4.75},
             'B213':
                {'wcs_position': [15*(4+19/60.), 27+15/60., 0],
                 'map': None,
                 'threshold': 4.75},
             'B220':
                {'wcs_position': [15*(4+41/60.), 26+7/60., 0],
                 'map': None,
                 'threshold': 7},
             'L1527':
                {'wcs_position': [15*(4+39/60.), 25+47/60., 0],
                 'map': None,
                 'threshold': 7},
             'B215':
                {'wcs_position': [15*(4+23/60.), 25+3/60., 0],
                 'map': None,
                 'threshold': 2},
             'L1524':
                {'wcs_position': [15*(4+29/60.), 24+31/60., 0],
                 'map': None,
                 'threshold': 3}}

    # calculate correlation for cores
    if not path.isfile(output_dir + 'core_arrays'):
    	system('mkdir ' + output_dir + 'core_arrays')
        for core in cores:
        	threshold = cores[core]['threshold']
        	av_clumps = clump_finder.find_clumps(av_data, threshold)
        	av_clumps_clean = clump_finder.find_clumps(av_clumps, 5)
        	print('Choose region for core ' + core)
        	region_mask = clump_finder.choose_region(av_clumps_clean,
        	        header=av_header, limits=[112,296,100,257])
        	cores[core]['map'] = np.copy(av_data)
        	cores[core]['map'][region_mask == 1] = np.NaN
        	np.save(output_dir + 'core_arrays/' + core, cores[core]['map'])

    for core in cores:
        plot_core(core_image=cores[core]['map'], header=av_header)

if __name__ == '__main__':
    main()



