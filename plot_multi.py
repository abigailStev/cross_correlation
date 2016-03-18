#!/usr/bin/env python
"""
Plots 1-D CCFs of multiple energy channels together on one plot.
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import subprocess
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator

__author__ = "Abigail Stevens <A.L.Stevens@uva.nl>"
__year__ = "2014-2016"

################################################################################
def main(ccf_file, plot_file, prefix):
    """
			main
	
	
    """
	
    try:
        in_table = Table.read(ccf_file)
    except IOError:
        print("\tERROR: File does not exist: %s" % ccf_file)
        exit()

    ccf = in_table['CCF']
    error = in_table['ERROR']
    n_bins = in_table.meta['N_BINS']
    # dt = in_table.meta['DT']

    time_bins = np.arange(n_bins)
    pos_time_bins = time_bins[0:n_bins/2]
    neg_time_bins = time_bins[n_bins/2:] - n_bins
    time_bins = np.append(neg_time_bins, pos_time_bins)

    ccf_3 = ccf[:, 3]
    err_3 = error[:, 3]
    pos_t_ccf_3 = ccf_3[0:n_bins/2]
    neg_t_ccf_3 = ccf_3[n_bins/2:]
    ccf_3 = np.append(neg_t_ccf_3, pos_t_ccf_3)
    pos_t_err_3 = err_3[0:n_bins/2]
    neg_t_err_3 = err_3[n_bins/2:]
    err_3 = np.append(neg_t_err_3, pos_t_err_3)

    ccf_9 = ccf[:, 9]
    err_9 = error[:, 9]
    pos_t_ccf_9 = ccf_9[0:n_bins/2]
    neg_t_ccf_9 = ccf_9[n_bins/2:]
    ccf_9 = np.append(neg_t_ccf_9, pos_t_ccf_9)
    pos_t_err_9 = err_9[0:n_bins/2]
    neg_t_err_9 = err_9[n_bins/2:]
    err_9 = np.append(neg_t_err_9, pos_t_err_9)

    ccf_15 = ccf[:, 15]
    err_15 = error[:, 15]
    pos_t_ccf_15 = ccf_15[0:n_bins/2]
    neg_t_ccf_15 = ccf_15[n_bins/2:]
    ccf_15 = np.append(neg_t_ccf_15, pos_t_ccf_15)
    pos_t_err_15 = err_15[0:n_bins/2]
    neg_t_err_15 = err_15[n_bins/2:]
    err_15 = np.append(neg_t_err_15, pos_t_err_15)

    ccf_24 = ccf[:, 24]
    err_24 = error[:, 24]
    pos_t_ccf_24 = ccf_24[0:n_bins/2]
    neg_t_ccf_24 = ccf_24[n_bins/2:]
    ccf_24 = np.append(neg_t_ccf_24, pos_t_ccf_24)
    pos_t_err_24 = err_24[0:n_bins/2]
    neg_t_err_24 = err_24[n_bins/2:]
    err_24 = np.append(neg_t_err_24, pos_t_err_24)
	
    #############
	## Plotting!
	#############
	
    font_prop = font_manager.FontProperties(size=20)
    fig, ax = plt.subplots(1, 1, figsize=(10, 6), tight_layout=True, dpi=300)

    ax.vlines(0.0, -20, 20.0, linestyle='dotted', color='black', lw=1.5)
    ax.hlines(0.0, -100, 100, linestyle='dashed', color='black', lw=1.5)
    ax.errorbar(time_bins, ccf_3, yerr=err_3, linewidth=2, color='orange',
                elinewidth=2, capsize=2, label="3.5 keV")
    ax.errorbar(time_bins, ccf_9, yerr=err_9, linewidth=2, elinewidth=2,
                capsize=2, label="6 keV")
    ax.errorbar(time_bins, ccf_15, yerr=err_15, linewidth=2, elinewidth=2,
                capsize=2, label="10.5 keV")
    ax.errorbar(time_bins, ccf_24, yerr=err_24, linewidth=2, elinewidth=2,
                capsize=2, label="18 keV")

    # ax.set_xlabel(r'Time ($\times\,$8.15$\,$ms)', fontproperties=font_prop)
    ax.set_xlabel('Time-delay bins', fontproperties=font_prop)
    ax.set_ylabel('Deviation from mean (cts$\;$s$^{-1}$)', \
    	fontproperties=font_prop)
    ax.set_xlim(-50, 50)
    ax.set_ylim(-4, 16)

    ## Setting the axes' minor ticks. It's complicated.
    x_maj_loc = ax.get_xticks()
    # x_maj_loc = [0, 50, 100]
    # ax.set_xticks(x_maj_loc)
    y_maj_loc = ax.get_yticks()

    x_min_mult = 0.25 * (x_maj_loc[1] - x_maj_loc[0])
    y_min_mult = 1.0
    xLocator = MultipleLocator(x_min_mult)  ## loc of minor ticks on x-axis
    yLocator = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax.xaxis.set_minor_locator(xLocator)
    ax.yaxis.set_minor_locator(yLocator)

    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
#     ax.set_title("%s, CCF per energy channel" % prefix, fontproperties=font_prop)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right', fontsize=18,
            borderpad=0.5, labelspacing=0.5, borderaxespad=0.5)
	
    plt.savefig(plot_file)
    # plt.show()
    plt.close()
    # subprocess.call(['cp', plot_file, \
    #         "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])


################################################################################
if __name__ == "__main__":
	
    parser = argparse.ArgumentParser(usage="python plot_multiccfs.py ccf_table"\
            " plot_file [-p prefix]", description="Plots CCFs of multiple "\
            "energy channels on one plot.")
	
    parser.add_argument('ccf_table', help="Name of file with CCF amplitudes in"\
            "a table.")
	
    parser.add_argument('plot_file', help="The output file name for the plot.")
	
    parser.add_argument('-p', '--prefix', dest='prefix', default="--",
            help="The identifying prefix of the data (object nickname or "\
            "proposal ID). [--]")

    args = parser.parse_args()

    # print("Input file: %s" % args.ccf_table)

    main(args.ccf_table, args.plot_file, args.prefix)
	
################################################################################
