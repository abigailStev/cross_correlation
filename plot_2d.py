import argparse
import numpy as np
from astropy.io import fits
import os.path
import subprocess
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator
import tools  # at https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2014-2015"
"""
Plots the ccf in a 2D colour plot, as time bins vs either energy channel or keV
energy.

"""

################################################################################
def make_plot(ccf, t_bins, mean_rate_ci, t_length, frac_time, plot_file,
        energies=None, prefix="", tab_file=None):

    ###################################################################
    ## Make a ratio of ccf to the mean count rate in the interest band
    ###################################################################

    mean_ccf = np.mean(ccf, axis=0)
    ccf_resid = ccf - mean_ccf
    print "CCF shape:", np.shape(ccf)
    a = np.array([mean_rate_ci, ] * 2*t_length).T
    print "A shape:", np.shape(a)
    with np.errstate(all='ignore'):
        ratio = np.where(a != 0, ccf / a, 0)

    # print "\tMinimum value:", np.min(ratio)
    # 	print "\tMaximum value:", np.max(ratio)

    ######################################################
    ## Saving to a dat file so that we can use fimgcreate
    ######################################################

# #     ratio[27:, ] = 0
    ratio[28:,] = 0
    out_file = os.path.dirname(plot_file) + "/temp.dat"
    print out_file
    R = ratio.flatten('C')
    comment_str = "From %s" % tab_file
    np.savetxt(out_file, R, comments=comment_str)

    #############
    ## Plotting!
    #############
    print("Plotting 2D CCF: %s" % plot_file)

    font_prop = font_manager.FontProperties(size=18)
    fig, ax = plt.subplots(1, 1, figsize=(10, 7.5), dpi=300, tight_layout=True)

    if energies[0] != None:  ## If energies exists as a variable
        plt.pcolor(t_bins, energies, ratio, cmap='hot')
        # plt.pcolor(t_bins, energies, ratio, cmap='hot', vmin=-0.26, vmax=0.42)
#         plt.pcolor(t_bins, energies, ratio, cmap='spring', vmin=-0.04, vmax=0.04)
    else:
# 		plt.pcolor(ratio, cmap='hot', vmin=-4.0, vmax=4.0)
        plt.pcolor(ratio, cmap='hot')
    # ax.vlines(0.0, 2, 31, linestyle='solid', color='black', lw=1.0)

    cbar = plt.colorbar()
    cbar.set_label('Ratio of CCF to mean count rate', \
            fontproperties=font_prop)
    # cbar.set_ticks([-0.04, -0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03, 0.04])
    # ax.set_xlabel(r'Time ( $\times\,\frac{1}{%d}\,$s)' % frac_time, \
    ax.set_xlabel(r'Time ( $\times\,1/%d\,$s)' % frac_time, \
            fontproperties=font_prop)

    if energies[0] != None:  ## If energies exists as a variable
        ax.set_ylabel('Energy (keV)', fontproperties=font_prop)
        # ax.set_ylim(3, 20)
        ax.set_ylim(3, 20)
        rect = patches.Rectangle((-t_length,energies[10]), 2*t_length, 0.41, \
                ec="none")
    else:
        ax.set_ylabel('Energy channel', fontproperties=font_prop)
        # ax.set_ylim(2, 31)
        ax.set_yscale('log')
        rect = patches.Rectangle((-t_length,10), 2*t_length, 1, ec="none")

    ax.add_patch(rect)
    zero_outline = patches.Rectangle((0, 2), 1, 26, edgecolor="black",
            facecolor="none")
    ax.add_patch(zero_outline)

    ax.set_xlim(-t_length, t_length)
    ## Setting the axes' minor ticks. It's complicated.
    x_maj_loc = ax.get_xticks()
    y_maj_loc = ax.get_yticks()
    # y_maj_loc = [5, 10, 15, 20]
    # ax.set_yticks(y_maj_loc)
    x_min_mult = 0.1 * (x_maj_loc[1] - x_maj_loc[0])
    y_min_mult = 0.2 * (y_maj_loc[1] - y_maj_loc[0])
    xLocator = MultipleLocator(x_min_mult)  ## loc of minor ticks on x-axis
    yLocator = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax.xaxis.set_minor_locator(xLocator)
    ax.yaxis.set_minor_locator(yLocator)

    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    # ax.set_title(prefix)

    # plt.show()
    plt.savefig(plot_file)
    plt.close()
    # subprocess.call(['cp', plot_file, "/Users/abigailstevens/Dropbox/Research/CCF_paper1/"])


################################################################################
def main(tab_file, plot_file, prefix, t_length, chan_to_en=None):
    """
    Makes a 2D plot of the cross-correlation functions per energy channel, plots
    it via keV energy if given, and saves plot to a specified file.

    Parameters
    ----------
    tab_file : string
        The file with the CCF table from ccf.py or multi_ccf.py, in FITS format.

    plot_file : string
        The file name to save the plot to.

    prefix : string
        The identifying prefix of the data (object nickname or proposal ID).

    t_length : int
        Number of time bins to use along the x-axis.

    chan_to_en : string
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

    try:
        file_hdu = fits.open(tab_file)
    except IOError:
        print("\tERROR: File does not exist: %s" % tab_file)
        exit()

    n_bins = int(file_hdu[0].header["N_BINS"])
    means_str = file_hdu[0].header["RATE_CI"]
    detchans = int(file_hdu[0].header["DETCHANS"])
    dt = float(file_hdu[0].header["DT"])
    frac_time = int(1.0 / dt)  ## Each time bin represents 1/frac_time sec

    table = file_hdu[1].data
    file_hdu.close()

    mean_rate_ci = [float(num.strip()) for num in means_str[1:-1].split(',')]
    ccf = np.reshape(table.field('CCF'), (n_bins, detchans), order='C')
    ccf = ccf.T  ## Transpose it to get the axes we want

    pos_time_ccf = ccf[:,0:n_bins/2]
    neg_time_ccf = ccf[:,n_bins/2:]

    ccf = np.hstack((neg_time_ccf, pos_time_ccf))

    time_bins = np.arange(n_bins)
    pos_time_bins = time_bins[0:n_bins/2]
    neg_time_bins = time_bins[n_bins/2:] - n_bins
    time_bins = np.append(neg_time_bins, pos_time_bins)

    ccf = ccf[:, n_bins/2-t_length:n_bins/2+t_length]
    t_bins = time_bins[n_bins/2-t_length:n_bins/2+t_length]

    make_plot(ccf, t_bins, mean_rate_ci, t_length, frac_time, plot_file,
            energies=energies, prefix=prefix, tab_file=tab_file)

################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python plot_2d.py tab_file [-o " \
            "plot_file] [-p prefix] [-e chan_to_en]", description="Plots the "\
            "ccf of multiple energy channels in a 2D colour plot.", epilog="" \
            "For optional arguments, default values are given in brackets at " \
            "end of description.")

    parser.add_argument('tab_file', help="The file with the CCF table from "\
            "ccf.py or multi_ccf.py, in .fits format.")

    parser.add_argument('-o', '--outfile', dest='plot_file', default="" \
            "./ccf_2D.png", help="The file name to save the plot to. " \
             "[./ccf_2D.png]")

    parser.add_argument('-p', '--prefix', dest='prefix', default="--", \
            help="The identifying prefix of the data (object nickname or "\
             "proposal ID). [--]")

    parser.add_argument('-l', '--length', dest='t_length', default=100, \
            type=tools.type_positive_int, help="Number of time bins to use " \
            "along the x-axis. [100]")

    parser.add_argument('-e', dest='chan_to_en', help="Table of actual energy " \
            "boundaries for energy channels as made in channel_to_energy.py. " \
            "If not given, plot will be vs detector mode energy channel. " \
            "[no default]")

    args = parser.parse_args()

    ## Idiot checks
    assert args.tab_file[-4:].lower() == "fits", "ERROR: Data file must be in "\
            ".fits format."
    if args.chan_to_en != None:  ## If args.chan_to_en exists as a variable
        assert os.path.isfile(args.chan_to_en), "ERROR: Table of real energy " \
                "per energy channel does not exist."

    ## Calling main
    main(args.tab_file, args.plot_file, args.prefix, args.t_length, \
        args.chan_to_en)

################################################################################