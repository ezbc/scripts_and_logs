#!/usr/bin/python

''' Script for plotting histograms of Av magnitudes of two Av images before and
after binning.
'''

def plot_av_histogram(image1, image2, range=range, title='',
        image1_label=None, image2_label=None, savedir='./',
        filename=None, show=True, returnimage=False):
    ''' Plots histogram of pixel values.
    '''

    # Import external modules
    import numpy as np
    import math
    import pyfits as pf
    import matplotlib.pyplot as plt
    import matplotlib

    image1_nonans = image1[image1 == image1]
    image2_nonans = image2[image2 == image2]

    # Create figure
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.hist(image1_nonans,color='g',alpha=0.5,range=range,bins=50,
            label=image1_label,normed=True)
    ax.hist(image2_nonans,color='b',alpha=0.5,range=range,bins=50,
            label=image2_label,normed=True)

    # Adjust asthetics
    ax.set_ylabel('Probability Density',
              size = 'small',
              family='serif')
    ax.set_xlabel(r'A$_{\rm v}$ (mag)',
              size = 'small',
              family='serif')
    ax.set_title(title)
    if image1_label is not None or image2_label is not None:
        ax.legend(loc='upper right')
    ax.grid(True)

    if filename is not None:
        plt.savefig(savedir + filename,bbox_inches='tight')
    if show:
        fig.show()
    if returnimage:
        return None

def load_fits(filename,return_header=False):
    ''' Loads a fits file.
    '''

    import pyfits as pf

    f = pf.open(filename)
    if return_header:
        return f[0].data,f[0].header
    else:
        return f[0].data

def main(redo_cube_correlation_calculation = False,
        redo_grid_correlation_calculation = False):

    ''' Loads fits files and plots histograms.
    '''

    output_dir = '/d/bip3/ezbc/taurus/data/python_output/nhi_av/'
    figure_dir = '/d/bip3/ezbc/taurus/figures/'
    av_dir = '/d/bip3/ezbc/taurus/data/av/'
    hi_dir = '/d/bip3/ezbc/taurus/data/galfa/'

    # load 2mass Av from Kainulainen et al. (2009)
    av_bin_jouni = load_fits(av_dir + 'taurus_av_k09_regrid.fits')
    av_nobin_jouni = load_fits(av_dir + 'taurus_av_kainulainen2009.fits')
    # load Av image from goldsmith: Pineda et al. 2010, ApJ, 721, 686
    # also from Qian et al. (2012), ApJ, 760, 147
    av_bin_goldsmith = load_fits(av_dir + \
            'taurus_av_p10_regrid.fits')
    av_nobin_goldsmith = load_fits(av_dir + \
            'taurus_av_pineda2010.fits')

    plot_av_histogram(av_bin_jouni,av_nobin_jouni,
            range=[-1,20],
            title=r'A$_{\rm v}$ Kainulainen et al. (2009) Histogram',
            savedir=figure_dir,
            image1_label='After Binning',
            image2_label='Before Binning',
            filename='taurus.av_histogram_kainulainen2009.png',)

    plot_av_histogram(av_bin_goldsmith,av_nobin_goldsmith,
            range=[-1,20],
            title=r'A$_{\rm v}$ Pineda et al. (2012) Histogram',
            savedir=figure_dir,
            image1_label='After Binning',
            image2_label='Before Binning',
            filename='taurus.av_histogram_pineda2012.png',)

    plot_av_histogram(av_nobin_goldsmith,av_nobin_jouni,
            range=[-1,20],
            title=r'A$_{\rm v}$ Pineda+12 + Kainulainen+2009 Pre-bin ' + \
                    'Histogram',
            image1_label='Pineda',
            image2_label='Kainulainen',
            savedir=figure_dir,
            filename='taurus.av_histogram_prebin_pineda+kainulainen.png',)

    plot_av_histogram(av_bin_goldsmith,av_bin_jouni,
            range=[-1,20],
            title=r'A$_{\rm v}$ Pineda+12 + Kainulainen+2009 Post-bin ' + \
                    'Histogram',
            image1_label='Pineda',
            image2_label='Kainulainen',
            savedir=figure_dir,
            filename='taurus.av_histogram_postbin_pineda+kainulainen.png',)

if __name__ == '__main__':
    main()

