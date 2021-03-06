#!/usr/bin/env python
"""
Plots the cross-correlation per energy channel, in the time domain.

Enter   python plot_ccf.py -h   at the command line for help.
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator

import tools

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2014-2017"


################################################################################
def type_positive_int(num):
    """
    Check if an input is a positive integer, as an argparse type.

    Parameters
    ----------
    num : int, long, float, or double
        The number in question.

    Returns
    -------
    n : int
        The input number, if it's a positive integer

    Raises
    ------
    ArgumentTypeError if n isn't a real number or a positive integer.

    """
    try:
        n = int(num)
    except ValueError or TypeError:
        message = "Input is not a positive integer."
        raise argparse.ArgumentTypeError(message)

    if n >= 0:
        return n
    else:
        message = "%d is not a positive integer." % n
        raise argparse.ArgumentTypeError(message)


################################################################################
def make_plot(x_bins, ccf_amps, ccf_err, n_bins, prefix, plot_file, chan, \
        delay, t_length=30):
    """
    Actually makes the plot.

    Parameters
    ----------
    x_bins : np.array of ints
        1-D array of the time bin values to plot along the x-axis.

    ccf_amps, ccf_err : np.arrays of floats
        1-D arrays of the CCF in a specific energy channel, to plot as the y-
        values, and the error on those CCF values.

    n_bins : int
        The number of time bins in one segment of the light curve; the length of
         x_bins, ccf_amps and ccf_err.

    prefix : str
        The identifying prefix of the data (object nickname or proposal ID).

    plot_file : str
        The full path name of the file to save the plot to, including extension.

    chan : int
        The energy channel of the CCF plotted.

    delay : float
        The time delay bin size in milliseconds, or dt*1000.

    t_length : int
        The number of time bins to plot on either side of 0. Optional,
        default=30.

    Returns
    -------
    nothing
    """

    font_prop = font_manager.FontProperties(size=20)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6), dpi=300, tight_layout=True)
    # fig, ax = plt.subplots(1, 1, figsize=(8,5), tight_layout=True)

    ax.plot(x_bins, ccf_amps, lw=2, c='black')
    ax.vlines(0.0, -2., 4.0, linestyle='dotted', color='black', lw=1.5)
    ax.hlines(0.0, -t_length, t_length, linestyle='dashed', color='black',
              lw=1.5)

    # ax.errorbar(x_bins, ccf_amps, yerr=ccf_err, lw=2, c='black',
    #             drawstyle='steps-mid', elinewidth=1.5, capsize=1.5)
    # ax.plot([-10], [ccf_amps[n_bins/2-10]], "o", mfc='red', mew=1, mec='black',
    #         ms=25)
    # ax.plot([-5], [ccf_amps[n_bins/2-5]],"*", mfc='orange', mew=1, mec='black',
    #         ms=35)
    # ax.plot([1], [ccf_amps[n_bins/2+1]], "^", mfc='green', mew=1, mec='black',
    #         ms=25)
    # ax.plot([6], [ccf_amps[n_bins/2+6]], 's', mfc='blue', mew=1, mec='black',
    #         ms=25)
    # ax.set_xlabel(r'Time ($\times\,$8.15 ms)', fontproperties=font_prop)
    ax.set_xlabel(r'Time-delay ($\times\,$%.2f$\,$ms)' % delay,
                  fontproperties=font_prop)
    ax.set_ylabel(r'Deviation from mean (cts / s)', fontproperties=font_prop)
    ax.set_xlim(-t_length, t_length)
    # ax.set_xlim(0, t_length)
    # ax.set_ylim(-1.5, 3.0)
    # ax.set_ylim(-0.5, 0.5)

    ## Setting the axes' minor ticks. It's complicated.
    x_maj_loc = ax.get_xticks()
    # x_maj_loc = [0, 50, 100]
    # ax.set_xticks(x_maj_loc)
    y_maj_loc = ax.get_yticks()

    x_min_mult = 0.2 * (x_maj_loc[1] - x_maj_loc[0])
    y_min_mult = 0.2 * (y_maj_loc[1] - y_maj_loc[0])
    xLocator = MultipleLocator(x_min_mult)  ## loc of minor ticks on x-axis
    yLocator = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax.xaxis.set_minor_locator(xLocator)
    ax.yaxis.set_minor_locator(yLocator)

    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    ax.tick_params(which='major', width=1.5, length=7)
    ax.tick_params(which='minor', width=1.5, length=4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    title="%s, Energy channel %d" % (prefix, chan)
    # ax.set_title(title, fontproperties=font_prop)

    plt.savefig(plot_file)
    # plt.show()
    plt.close()


################################################################################
if __name__ == "__main__":

    ###########################
    ## Parsing input arguments
    ###########################

    parser = argparse.ArgumentParser(usage="python plot_ccf.py ccf_file "\
            "[OPTIONAL ARGUMENTS]", description=__doc__, epilog="For optional "\
            "arguments, default values are given in brackets at end of "\
            "description.")

    parser.add_argument('ccf_file', help="The CCF file, saved from an astropy "\
            "table, in .fits format.")

    parser.add_argument('-o', '--plotroot', dest='plot_root', default="./ccf",
            help="The root of the filename to save the 1-D ccf plot to. Energy"\
            " channel will be appended to name before saving. [./ccf]")

    parser.add_argument('-p', '--prefix', dest='prefix', default="--",
            help="The identifying prefix of the data (object nickname or data "\
            "ID). [--]")

    parser.add_argument('-e', '--ext', dest='plot_ext', default='eps',
            help="File extension for the plot. Do not include the '.' [eps]")

    parser.add_argument('-l', '--length', dest='t_length', default=30,
            type=type_positive_int,
            help="Number of time bins to plot on either side of 0. [30]")

    args = parser.parse_args()

    print("Plotting the ccf: %s_chan_xx.%s" % (args.plot_root, args.plot_ext))

    try:
        in_table = Table.read(args.ccf_file)
    except IOError:
        print("\tERROR: File does not exist: %s" % args.ccf_file)
        exit()

    ccf = in_table['CCF']
    error = in_table['ERROR']
    # print np.shape(ccf)
    # print np.shape(error)
    n_bins = in_table.meta['N_BINS']
    n_seconds = in_table.meta['SEC_SEG']
    dt = in_table.meta['DT']
    frac_time = int(1.0/dt)  ## each time bin represents 1/frac_time seconds

    for i in range(15, 16):
        ccf_i = ccf[:, i]
        ccf_err_i = error[:, i]
        # print np.shape(ccf_i)
        # print np.shape(ccf_err_i)
        time_bins = np.arange(n_bins)

        pos_time_bins = time_bins[0:n_bins/2]
        neg_time_bins = time_bins[n_bins/2:] - n_bins
        time_bins = np.append(neg_time_bins, pos_time_bins)

        pos_time_ccf = ccf_i[0:n_bins/2]
        neg_time_ccf = ccf_i[n_bins/2:]
        ccf_i = np.append(neg_time_ccf, pos_time_ccf)

        pos_time_ccf_err = ccf_err_i[0:n_bins/2]
        # neg_time_ccf_err = pos_time_ccf_err[::-1]
        neg_time_ccf_err = ccf_err_i[n_bins/2:]
        ccf_err_i = np.append(neg_time_ccf_err, pos_time_ccf_err)

        if i < 10:
            plot_file = args.plot_root + "_chan_" + str(0) + str(i) + "." + \
                        args.plot_ext
        else:
            plot_file = args.plot_root + "_chan_" + str(i) + "." + \
                        args.plot_ext


        make_plot(time_bins, ccf_i, ccf_err_i, n_bins, args.prefix, plot_file,
                i, dt*1000., args.t_length)

################################################################################
