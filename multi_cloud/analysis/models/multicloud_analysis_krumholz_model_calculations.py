#!/usr/bin/python

def plot_phi_vs_T(T_cnm, phi_cnm, scales=('linear', 'linear'), filename=None,
        G0s=None, limits=None,
        show=False):

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    ax.plot(T_cnm, phi_cnm, linestyle='-', marker='', color='k',)

    ax.set_xlabel(r'T$_{CNM}$ (K)')
    ax.set_ylabel(r'$\phi_{CNM}$')
    ax.legend()

    ax.set_xscale(scales[0])
    ax.set_yscale(scales[1])
    if limits is not None:
    	ax.set_xlim(limits[0], limits[1])
    	ax.set_ylim(limits[2], limits[3])

    if show:
        plt.show()
    if filename is not None:
    	plt.savefig(filename)

    plt.close()

def plot_phi_vs_n(T_cnm_list, phi_cnm_list, scales=('linear', 'linear'),
        filename=None, G0s=None, show=False, limits=None):

    import matplotlib.pyplot as plt
    from myscience.krumholz09 import calc_phi_cnm, calc_T_cnm, calc_n_cnm

    fig, ax = plt.subplots()

    n_cnm_list = []

    for i, G0 in enumerate(G0s):
    	n_cnm_list.append(calc_n_cnm(G_0=G0, T_cnm=T_cnm_list))

    for i, G0 in enumerate(G0s):
        ax.plot(n_cnm_list[i], phi_cnm_list, linestyle='-', marker='',
                label=r'$G_0$ = {0:.2f} $G_\odot$'.format(G0))

    ax.set_xlabel(r'n$_{CNM}$ (cm$^{-3}$)')
    ax.set_ylabel(r'$\phi_{CNM}$')
    ax.legend(loc='upper left')

    ax.set_xscale(scales[0])
    ax.set_yscale(scales[1])
    if limits is not None:
    	ax.set_xlim(limits[0], limits[1])
    	ax.set_ylim(limits[2], limits[3])

    if show:
        plt.show()
    if filename is not None:
    	plt.savefig(filename)

    plt.close()

def plot_phi_vs_P(T_cnm_list, phi_cnm_list, scales=('linear', 'linear'),
        filename=None, G0s=None, show=False, limits=None):

    import matplotlib.pyplot as plt
    from myscience.krumholz09 import calc_phi_cnm, calc_T_cnm, calc_n_cnm
    import numpy as np

    fig, ax = plt.subplots()

    n_cnm_list = []
    P_list = []

    for i, G0 in enumerate(G0s):
    	n_cnm_list.append(calc_n_cnm(G_0=G0, T_cnm=T_cnm_list))
    	P_list.append(np.asarray(n_cnm_list[i]) * np.asarray(T_cnm_list))

    for i, G0 in enumerate(G0s):
        ax.plot(P_list[i], phi_cnm_list, linestyle='-', marker='',
                label=r'$G_0$ = {0:.2f} $G_\odot$'.format(G0))

    ax.set_xlabel(r'P (K / cm$^{3}$)')
    ax.set_ylabel(r'$\phi_{CNM}$')
    ax.legend(loc='upper left')

    ax.set_xscale(scales[0])
    ax.set_yscale(scales[1])
    if limits is not None:
    	ax.set_xlim(limits[0], limits[1])
    	ax.set_ylim(limits[2], limits[3])

    if show:
        plt.show()
    if filename is not None:
    	plt.savefig(filename)

    plt.close()

def plot_spin_temps(spin_temps, bins=10, scales=('linear', 'linear'),
        filename=None, show=False):

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    ax.hist(spin_temps, color='k', bins=bins, alpha=0.2)

    ax.set_xlabel(r'T$_{s}$')
    ax.set_ylabel(r'Counts')

    ax.set_xscale(scales[0])
    ax.set_yscale(scales[1])

    if show:
        plt.show()
    if filename is not None:
    	plt.savefig(filename)

    plt.close

def main():

    from myscience.krumholz09 import calc_phi_cnm, calc_T_cnm, calc_n_cnm
    import numpy as np

    T_cnm = np.arange(30, 10000, 0.1)
    G0s = (0.01, 0.1, 1, 10)

    phi_cnm = calc_phi_cnm(T_cnm, Z=1.0)

    plot_phi_vs_T(T_cnm, phi_cnm,
            filename='/d/bip3/ezbc/multicloud/figures/phi_cnm_vs_T_cnm.png',
            scales=('linear','log'),
            limits=[10, 200, 0.1, 100],
            )

    plot_phi_vs_n(T_cnm, phi_cnm,
            G0s=G0s,
            scales=('log', 'linear'),
            limits=[10**-2, 10**4.5, 0.1, 160],
            filename='/d/bip3/ezbc/multicloud/figures/phi_cnm_vs_n_cnm.png',
            )

    plot_phi_vs_P(T_cnm, phi_cnm,
            G0s=G0s,
            scales=('log', 'log'),
            limits=[10**0, 10**6, 0.18, 160],
            filename='/d/bip3/ezbc/multicloud/figures/phi_cnm_vs_P.png',
            )

    '''
    T_cnm = calc_T_cnm(27.29, Z=1.0)
    n_cnm = calc_n_cnm(G_0=.5, T_cnm=T_cnm)
    print 'T_cnm =', T_cnm, 'phi_cnm =', 10.0, 'n_cnm =', n_cnm, ' cm^-3'
    print 'P/k =', T_cnm * n_cnm, ' K/cm^-3'
    T_cnm = calc_T_cnm(27.29 + 5.19, Z=1.0)
    print 'T_cnm =', T_cnm, 'phi_cnm =', 10.0
    T_cnm = calc_T_cnm(27.29 - 5.18, Z=1.0)
    print 'T_cnm =', T_cnm, 'phi_cnm =', 10.0
    phi_cnm = calc_phi_cnm(15.0, Z=1.0)

    perseus_phi_cnm = (9.61, 8.37, 8.32, 7.95, 6.26)
    perseus_T_cnm = []
    for phi_cnm in perseus_phi_cnm:
    	perseus_T_cnm.append(calc_T_cnm(phi_cnm, Z=1.0))

    print 'perseus T_cnm', perseus_T_cnm

    print 'T_cnm =', 15.0, 'phi_cnm =', phi_cnm

    # Plot the distribution of cold sources near Taurus from Heiles and Troland
    # 2003
    # --------------------------------------------------------------------------
    ht_dict = {'3c141.0' : {'l' : 174.5,
                            'b' : -1.3,
                            'cold Ts' : (9., 48., 23., 119.)},
               '3c133' : {'l' : 177.7,
                            'b' : -9.9,
                            'cold Ts' : ()},
               'a' : {'l' : 1,
                            'b' : 1,
                            'cold Ts' : (1)},
                                        }
    spin_temps = (9, 48, 23, 119, 18, 66, 44, 38, 80, 78, 45, 18, 22,
            41, 17, 35, 28, 21, 41, 17, 35, 28, 21, 41, 59, 26, 28, 33, 23)

    kinetic_temps = ()

    plot_spin_temps(spin_temps, bins=20,
            filename=\
                '/d/bip3/ezbc/multicloud/figures/heiles03_spin_temp_hist.png')
    '''

if __name__ == '__main__':
	main()

