#!/usr/bin/env python
"""
Plots the ccf in a 2D colour plot, as time bins vs either energy channel or keV
energy.

"""
import argparse
import numpy as np
import os.path
import subprocess
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator
import tools  # at https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2014-2017"


################################################################################
def make_plot(ccf, t_bins, mean_rate_ci, t_length, delay,
        plot_file, energies=[], prefix="--", tab_file=None):
    """
    Actually makes the plots.

    Parameters
    ----------
    ccf : np.array of floats
        2-D array of the cross-correlation function.

    t_bins : np.array of ints
        1-D array of integer time bins.

    mean_rate_ci : np.array of floats
        1-D array of the mean count rate of the channels of interest, in cts/s.
        Size = detchans.

    t_length : int
        The number of time bins before and after zero to plot (so, 2*t_length
        get plotted along the x-axis).

    delay : float
        The time delay bin size in milliseconds, or dt*1000.

    plot_file : str
        The name of the file to save the 2-D CCF plot to.

    energies : np.array of floats
        1-D array of the keV energy bounds of each detector energy channel.
        Size = (detchans+1)

    prefix : str
        The identifying prefix of the data (object nickname or data ID). [--]

    tab_file : str
        The file with the CCF table from ccf.py or multi_ccf.py, in FITS format,
        for writing to temp.dat file for making the FIMGCREATE image in the bash
        script.

    Returns
    -------
    nothing

    """
    ###################################################################
    ## Make a ratio of ccf to the mean count rate in the interest band
    ###################################################################

    # mean_ccf = np.mean(ccf, axis=0)
    # ccf_resid = ccf - mean_ccf
    a = np.array([mean_rate_ci, ] * (2 * t_length + 1)).T
    with np.errstate(all='ignore'):
        ratio = np.where(a != 0, ccf / a, 0)

    # print "\tMinimum value:", np.min(ratio)
    # 	print "\tMaximum value:", np.max(ratio)

    ######################################################
    ## Saving to a dat file so that we can use fimgcreate
    ######################################################

# #     ratio[27:, ] = 0
    if np.shape(ratio)[0] == 64:
        ratio[28:,] = 0
    elif np.shape(ratio)[0] == 32:
        ratio[26:,] = 0
    out_file = os.path.dirname(plot_file) + "/temp.dat"
    R = ratio.real.flatten('C')
    comment_str = "From %s" % tab_file
    np.savetxt(out_file, R, fmt="%.8f", comments=comment_str)

    #############
    ## Plotting!
    #############
    print("Plotting 2D CCF: %s" % plot_file)

    font_prop = font_manager.FontProperties(size=20)
    # fig, ax = plt.subplots(1, 1, figsize=(14, 8), dpi=300, tight_layout=True)
    fig, ax = plt.subplots(1, 1, figsize=(10, 7.5), dpi=300, tight_layout=True)

    if len(energies) > 0:  ## If energies exists as a variable
        # plt.pcolor(t_bins, energies, ratio, cmap='YlGnBu_r', vmin=-0.3, vmax=0.3)
        # plt.pcolor(t_bins, energies, ratio, cmap='YlGnBu_r', vmin=-0.05, vmax=0.05)

        plt.pcolor(t_bins, energies, ratio, cmap='spring')
        # plt.pcolor(t_bins, energies, ratio, cmap='hot', vmin=-0.26, vmax=0.42)
        # plt.pcolor(t_bins, energies, ratio, cmap='spring', vmin=-0.04,
        #         vmax=0.04)
    else:
# 		plt.pcolor(ratio, cmap='hot', vmin=-4.0, vmax=4.0)
        plt.pcolor(ratio, cmap='hot')
        plt.xlim(t_bins[0], t_bins[-1])
    # ax.vlines(0.0, 2, 31, linestyle='solid', color='black', lw=1.0)

    cbar = plt.colorbar()
    cbar.set_label('Ratio of CCF to mean count rate', \
            fontproperties=font_prop)
    cb_ax = cbar.ax
    cb_ax.tick_params(axis='y', labelsize=18)
    # cbar.set_ticks([-0.04, -0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03, 0.04])

    if len(energies) > 0:  ## If energies exists as a variable
        ax.set_ylabel('Energy (keV)', fontproperties=font_prop)
        if len(energies) == 65:
            ax.set_ylim(3, 20)
        # rect = patches.Rectangle((-t_length,energies[10]), 2*t_length, 0.41,
        #         facecolor="orange", ec="none")
            rect = patches.Rectangle((-t_length, energies[10]), 2 * t_length, 0.41,
                    facecolor="black", ec="none")
            ax.add_patch(rect)
    else:
        ax.set_ylabel('Energy channel', fontproperties=font_prop)
        ax.set_ylim(0, np.shape(ratio)[0])
        # rect = patches.Rectangle((-t_length,10), 2*t_length, 1, ec="none")

    zero_outline = patches.Rectangle((0, 2), 0.5, 26, edgecolor="black",
            facecolor="none")
    ax.add_patch(zero_outline)

    ax.set_xlim(-t_length, t_length)
    # ax.set_xlabel(r'Time-delay ($\times\,$8.15$\,$ms)', fontproperties=font_prop)
    # ax.set_xlabel('Time-delay bins ', fontproperties=font_prop)
    ax.set_xlabel(r'Time-delay ($\times\,$%.2f$\,$ms)' % delay,
                  fontproperties=font_prop)

    ## Setting the axes' minor ticks. It's complicated.
    x_maj_loc = ax.get_xticks()
    y_maj_loc = ax.get_yticks()
    # y_maj_loc = [5, 10, 15, 20]
    # ax.set_yticks(y_maj_loc)
    x_min_mult = 0.2 * (x_maj_loc[1] - x_maj_loc[0])
    y_min_mult = 0.2 * (y_maj_loc[1] - y_maj_loc[0])
    xLocator = MultipleLocator(x_min_mult)  ## loc of minor ticks on x-axis
    yLocator = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax.xaxis.set_minor_locator(xLocator)
    ax.yaxis.set_minor_locator(yLocator)

    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.tick_params(which='major', width=1.5, length=7)
    ax.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    # plt.show()
    plt.savefig(plot_file)
    plt.close()


################################################################################
def main(tab_file, plot_file, prefix, t_length, chan_to_en=None):
    """
    Makes a 2D plot of the cross-correlation functions per energy channel, plots
    it via keV energy if given, and saves plot to a specified file.

    Parameters
    ----------
    tab_file : str
        The file with the CCF table from ccf.py saved as an Astropy Table, with
        column 1 as the CCF and column 2 as the errors.

    plot_file : str
        The file name to save the plot to.

    prefix : str
        The identifying prefix of the data (object nickname or data ID).

    t_length : int
        Number of time bins to use along the x-axis.

    chan_to_en : str
        Table of actual energy boundaries for energy channels as made in
        channel_to_energy.py. If not given, plot will be vs detector mode energy
        channel.

    Returns
    -------
    nothing

    """

    ##########################################
    ## Load data from FITS table and txt file
    ##########################################

    if chan_to_en != None:  ## If chan_to_en exists as a variable
        energies = np.loadtxt(chan_to_en)
    else:
        energies = []

    try:
        in_table = Table.read(tab_file)
    except IOError:
        print("\tERROR: File does not exist: %s" % tab_file)
        exit()

    ## Need to transpose it here so that it plots with time on the x-axis
    ## and energy on the y-axis
    ccf = in_table['CCF'].T
    # error = in_table['ERROR']

    n_bins = in_table.meta['N_BINS']
    means_str = in_table.meta['RATE_CI']
    dt = in_table.meta['DT']
    print_dt = dt * 10**3  ## Putting it in milliseconds units.

    mean_rate_ci = [float(num.strip()) for num in means_str[1:-1].split(',')]

    pos_time_ccf = ccf[:,0:n_bins/2]
    neg_time_ccf = ccf[:,n_bins/2:]

    ccf = np.hstack((neg_time_ccf, pos_time_ccf))

    time_bins = np.arange(n_bins)
    pos_time_bins = time_bins[0:n_bins/2]
    neg_time_bins = time_bins[n_bins/2:] - n_bins
    time_bins = np.append(neg_time_bins, pos_time_bins)

    ccf = ccf[:, n_bins/2-t_length:n_bins/2+t_length+1]
    t_bins = time_bins[n_bins/2-t_length:n_bins/2+t_length+1]

    make_plot(ccf, t_bins, mean_rate_ci, t_length, print_dt, plot_file,
            energies=energies, prefix=prefix, tab_file=tab_file)


################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python plot_2d.py tab_file " \
            "[OPTIONAL ARGUMENTS]", description=__doc__,
            epilog="For optional arguments, default values are given in "\
            "brackets at end of description.")

    parser.add_argument('tab_file', help="The .fits file with the CCF table "\
            "from ccf.py saved as an Astropy Table, with column 1 as the CCF "\
            "and column 2 as the errors.")

    parser.add_argument('-o', '--outfile', dest='plot_file', default=""\
            "./ccf_2D.png", help="The file name to save the plot to. "\
             "[./ccf_2D.png]")

    parser.add_argument('-p', '--prefix', dest='prefix', default="--", \
            help="The identifying prefix of the data (object nickname or "\
             "data ID). [--]")

    parser.add_argument('-l', '--length', dest='t_length', default=100, \
            type=tools.type_positive_int, help="Number of time bins to use "\
            "along the x-axis. [100]")

    parser.add_argument('-e', dest='chan_to_en', help="Table of actual energy "\
            "boundaries for energy channels as made in channel_to_energy.py. "\
            "If not given, plot will be vs detector mode energy channel. "\
            "[no default]")

    args = parser.parse_args()

    ## Idiot checks
    assert args.tab_file[-4:].lower() == "fits", "ERROR: Data file must be in "\
            ".fits format."
    if args.chan_to_en != None:  ## If args.chan_to_en exists as a variable
        assert os.path.isfile(args.chan_to_en), "ERROR: Table of real energy "\
                "per energy channel does not exist."

    ## Calling main
    main(args.tab_file, args.plot_file, args.prefix, args.t_length, \
        chan_to_en=args.chan_to_en)

################################################################################