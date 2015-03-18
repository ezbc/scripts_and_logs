#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

import warnings
warnings.filterwarnings('ignore')

def plot_mask_residuals(residuals=None, x_fit=None, y_fit=None,
        residual_thres=None, filename=None, show=True, title=None):

    # Import external modules
    import matplotlib.pyplot as plt
    import matplotlib
    import numpy as np
    from scipy.integrate import simps as integrate

    # Set up plot aesthetics
    # ----------------------
    plt.close;plt.clf()
    plt.rcdefaults()

    # Color map
    cmap = plt.cm.gnuplot

    # Color cycle, grabs colors from cmap
    color_cycle = [cmap(i) for i in np.linspace(0, 0.8, 2)]
    font_scale = 9
    line_weight = 600
    font_weight = 600
    params = {#'backend': .pdf',
              'axes.labelsize': font_scale,
              'axes.titlesize': font_scale,
              'axes.weight': line_weight,
              'text.fontsize': font_scale,
              'legend.fontsize': font_scale*3/4,
              'xtick.labelsize': font_scale,
              'xtick.weight': line_weight,
              'ytick.labelsize': font_scale,
              'ytick.weight': line_weight,
              'font.weight': font_weight,
              'axes.labelweight': font_weight,
              'text.usetex': True,
              #'font.family': 'sans-serif',
              'figure.figsize': (7.3/2.0, 7.3/4.0),
              'figure.dpi': 600,
              'figure.titlesize': font_scale,
              'axes.color_cycle': color_cycle # colors of different plots
             }
    plt.rcParams.update(params)

    # Create figure instance
    fig = plt.figure()

    ax = fig.add_subplot(111)

    counts, bin_edges = \
        np.histogram(np.ravel(residuals[~np.isnan(residuals)]),
                                 bins=1000,
                                 )

    bin_edges_ext = np.zeros(len(counts) + 1)
    counts_ext = np.zeros(len(counts) + 1)

    bin_edges_ext[0] = bin_edges[0] - (bin_edges[1] - bin_edges[0])
    bin_edges_ext[1:] = bin_edges[:-1]
    counts_ext[0] = 0
    counts_ext[1:] = counts

    # Normalize so area = 1
    #counts_ext /= np.nansum(counts_ext) * (bin_edges_ext[2] - bin_edges_ext[1])
    counts_ext = counts_ext / integrate(counts_ext, )#x=bin_edges_ext)
    y_fit /= np.max(y_fit)
    y_fit *= np.max(counts_ext)

    ax.plot(bin_edges_ext, counts_ext, drawstyle='steps-mid',
            linewidth=1.5)
    ax.plot(x_fit, y_fit,
            linewidth=2,
            alpha=0.8)
    ax.set_xlim([np.nanmin(bin_edges_ext) - \
                 np.abs(0.8 * np.nanmin(bin_edges_ext)),4])
    #ax.set_ylim([-0.1, 1.1])
    ax.axvline(residual_thres,
               color='k',
               linestyle='--',
               linewidth=1.5)
    ax.set_xlabel(r'Residual $A_V$ [mag]')
    ax.set_ylabel('Normalized PDF')

    plt.locator_params(axis='y', nbins=5)

    if title is not None:
        fig.suptitle(title, fontsize=font_scale)
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight', dpi=600)
    if show:
        plt.show()


def main():

    import json
    import numpy as np

    with open('/d/bip3/ezbc/california/data/python_output/' + \
              'residual_parameter_results/' + \
              'california_residuals.json', 'r') as f:
        residuals_dict = json.load(f)

    residuals = np.array(residuals_dict['residuals'])
    residual_thres = 1.5 * residuals_dict['width'] + residuals_dict['x0']
    x_fit = np.array(residuals_dict['x_fit'])
    y_fit = np.array(residuals_dict['y_fit'])

    for filetype in ['.png', '.pdf']:
        plot_mask_residuals(residuals=residuals,
                            x_fit=x_fit,
                            y_fit=y_fit,
                            residual_thres=residual_thres,
                            filename='/d/bip3/ezbc/california/figures/' + \
                                     'likelihood/california_residual_pdf' + \
                                     filetype,
                            )

if __name__ == '__main__':
    main()
