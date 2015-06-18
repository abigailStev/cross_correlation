import argparse
import numpy as np
from astropy.io import fits
import os.path
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator
import tools  # at https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens, A.L.Stevens at uva.nl"

"""
Plots the ccf of multiple energy channels in a 2D colour plot.

2014-2015

"""

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

    print "Plotting the 2D CCF: %s" % plot_file

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
    table = file_hdu[1].data
    file_hdu.close()

    mean_rate_ci = [float(num.strip()) for num in means_str[1:-1].split(',')]
    ccf = np.reshape(table.field('CCF'), (n_bins, detchans), order='C')

    ccf = ccf[0:t_length, ].T  ## Transpose it to get the axes we want
    time_bins = np.arange(n_bins)
    frac_time = int(1.0 / dt)  ## Each time bin represents 1/frac_time sec

    ###################################################################
    ## Make a ratio of ccf to the mean count rate in the interest band
    ###################################################################

    mean_ccf = np.mean(ccf, axis=0)
    ccf_resid = ccf - mean_ccf

    a = np.array([mean_rate_ci, ] * t_length).T
    with np.errstate(all='ignore'):
        ratio = np.where(a != 0, ccf / a, 0)

    # print "\tMinimum value:", np.min(ratio)
    # 	print "\tMaximum value:", np.max(ratio)

    ######################################################
    ## Saving to a dat file so that we can use fimgcreate
    ######################################################

# 	ratio[np.where(np.abs(ratio) > 4.0)] = 0
#     ratio[27:, ] = 0
    ratio[32:,] = 0
    out_file = "./temp.dat"
    R = ratio.flatten('C')
    comment_str = "From %s" % tab_file
    np.savetxt(out_file, R, comments=comment_str)

    #############
    ## Plotting!
    #############

    t_bins = np.arange(t_length + 1)
    font_prop = font_manager.FontProperties(size=16)
    fig, ax = plt.subplots(1, 1)

    if energies[0] != None:  ## If energies exists as a variable
# 		plt.pcolor(t_bins, energies, ratio, cmap='hot')
#         plt.pcolor(t_bins, energies, ratio, cmap='hot', vmin=-1, vmax=1.5)
        plt.pcolor(t_bins, energies, ratio, cmap='spring', vmin=-0.04, vmax=0.04)
    else:
# 		plt.pcolor(ratio, cmap='hot', vmin=-4.0, vmax=4.0)
        plt.pcolor(ratio, cmap='hot')

    cbar = plt.colorbar()
    cbar.set_label('Ratio of CCF to mean count rate', \
        fontproperties=font_prop)
    cbar.set_ticks([-0.04, -0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03, 0.04])
    ax.set_xlabel(r'Time ( $\times\,\frac{1}{%d}\,$s)' % frac_time, \
        fontproperties=font_prop)

    if energies[0] != None:  ## If energies exists as a variable
        ax.set_ylabel('Energy (keV)', fontproperties=font_prop)
        ax.set_ylim(3, 25)
        # ax.set_ylim(3, 20)
    else:
        ax.set_ylabel('Energy channel', fontproperties=font_prop)
        ax.set_ylim(2, 31)
        ax.set_yscale('log')

    ## Setting the axes' minor ticks. It's complicated.
    x_maj_loc = ax.get_xticks()
    y_maj_loc = ax.get_yticks()

    # print x_maj_loc
    # print y_maj_loc
    # y_maj_loc = [5, 10, 15, 20, 25]
    # ax.set_yticks(y_maj_loc)
    x_min_mult = 0.1 * (x_maj_loc[1] - x_maj_loc[0])
    y_min_mult = 0.2 * (y_maj_loc[1] - y_maj_loc[0])
    xLocator = MultipleLocator(x_min_mult)  ## loc of minor ticks on x-axis
    yLocator = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax.xaxis.set_minor_locator(xLocator)
    ax.yaxis.set_minor_locator(yLocator)

    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    # ax.set_title(prefix)

    fig.set_tight_layout(True)
    plt.savefig(plot_file, dpi=200)
    # 	plt.show()
    plt.close()


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