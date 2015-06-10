import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator

__author__ = "Abigail Stevens, A.L.Stevens at uva.nl"

"""
        plot_ccf.py

Plots the cross-correlation per energy channel, in the time domain.

Enter   python plot_ccf.py -h   at the command line for help.

Abigail Stevens, A.L.Stevens at uva.nl, 2014-2015

"""

################################################################################
def make_plot(x_bins, ccf_amps, ccf_err, prefix, plot_file, chan, frac_time):
    """
            make_plot

    Actually makes the plot.

    """

    font_prop = font_manager.FontProperties(size=20)

    fig, ax = plt.subplots(1,1, figsize=(12,6))
# 	ax.plot(x_bins, ccf_amps, lw=2, c='black')
    print( len(ccf_err) )
    print( len(ccf_amps) )
    ax.errorbar(x_bins, ccf_amps, yerr=ccf_err, lw=1.5, c='black', \
        drawstyle='steps-mid', elinewidth=1, capsize=1)
    ax.vlines(0.0, -1.5, 2.5, linestyle='dashed')
# 	ax.plot([6], [ccf_amps[6]], "o", mfc='red', mew=1, mec='black', ms=14)
# 	ax.plot([13], [ccf_amps[13]],"*",  mfc='orange', mew=1, mec='black', ms=20)
# 	ax.plot([19], [ccf_amps[19]], "^", mfc='green', mew=1, mec='black', ms=14)
# 	ax.plot([24], [ccf_amps[24]], 's', mfc='blue', mew=1, mec='black', ms=14)
    ax.set_xlabel(r'Time ($\times\,\frac{1}{%d}\,$s)' % frac_time, \
        fontproperties=font_prop)
    ax.set_ylabel(r'Deviation from mean rate (cts / s)', \
        fontproperties=font_prop)
    ax.set_xlim(-100, 100)
    ax.set_ylim(-1.5, 2.5)

    ## Setting the axes' minor ticks. It's complicated.
    x_maj_loc = ax.get_xticks()
    y_maj_loc = ax.get_yticks()
    x_min_mult = 0.1 * (x_maj_loc[1] - x_maj_loc[0])
    y_min_mult = 0.2 * (y_maj_loc[1] - y_maj_loc[0])
    xLocator = MultipleLocator(x_min_mult)  ## loc of minor ticks on x-axis
    yLocator = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
    ax.xaxis.set_minor_locator(xLocator)
    ax.yaxis.set_minor_locator(yLocator)

    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    # ax.set_title(prefix + ", Energy channel " + str(chan), fontproperties=font_prop)

    fig.set_tight_layout(True)
    plt.savefig(plot_file, dpi=200)
# 	plt.show()
    plt.close()

## End of function 'make_plot'


################################################################################
if __name__ == "__main__":

    ###########################
    ## Parsing input arguments
    ###########################

    parser = argparse.ArgumentParser(usage="python plot_ccf.py tab_file [-o \
out_root] [-p prefix]", description="Plots the cross-correlation per energy \
channel, in the time domain.", epilog="For optional arguments, default values \
are given in brackets at end of description.")

    parser.add_argument('tab_file', help="The table file, in .dat or .fits \
format.")

    parser.add_argument('-o', '--outroot', dest='out_root', default="./ccf", \
help="The root of the filename to save the plot to. Energy channel will be \
appended to name before saving. [./ccf]")

    parser.add_argument('-p', '--prefix', dest='prefix', default="--", \
help="The identifying prefix of the data (object nickname or proposal ID). \
[--]")

    parser.add_argument('-e', '--ext', dest='plot_ext', default='png', \
help="File extension for the plot. Do not include the '.' [png]")

    args = parser.parse_args()


    print( "Plotting the ccf: %s_chan_xx.%s" % (args.out_root, args.plot_ext) )

    # if args.tab_file[-4:].lower() == ".dat":
    #
    # 	table = np.loadtxt(args.tab_file, comments='#')  ## table of CCF
    # 	time_bins = table[:,0]
    # 	n_bins = len(time_bins)
    # 	ccf = np.zeros((n_bins, 64))
    # 	ccf_err = np.zeros((n_bins, 64))
    # 	for i in range(11, 12):
    # 		ccf[:,i] = table[:, i+1]
    # 		ccf_err[:,i] = table[:, i+65]
    # 		c = i
    # 		if c >= 64:
    # 			c -= 64
    # 		if c < 10:
    # 			plot_file = args.out_root + "_chan_" + str(0) + str(c) + ".png"
    # 		else:
    # 			plot_file = args.out_root + "_chan_" + str(c) + ".png"
    #
    #
    #
    # 		make_plot(time_bins, ccf[:, i], ccf_err[:, i], args.prefix, plot_file, c, 0)

    # elif args.tab_file[-5:].lower() == ".fits":

    try:
        file_hdu = fits.open(args.tab_file)
    except IOError:
        print( "\tERROR: File does not exist: %s" % args.tab_file )
        exit()

    table = file_hdu[1].data
    header = file_hdu[0].header
    file_hdu.close()

    dt = float(header["DT"])
    num_seconds = int(header['SEC_SEG'])
    n_bins = int(header['N_BINS'])
    frac_time = int(1.0/dt)  ## each time bin represents 1/frac_time seconds

    for i in range(15, 16):
#     for i in range(6, 7):
        channel_mask = table.field('CHANNEL') == i  ## make data mask for the energy channels i want
        table_i = table[channel_mask]
        ccf = table_i.field('CCF')
        ccf_err = table_i.field('ERROR')
        time_bins = table_i.field('TIME_BIN')

        pos_time_bins = time_bins[0:n_bins/2]
        neg_time_bins = time_bins[n_bins/2:] - n_bins
        time_bins = np.append(neg_time_bins, pos_time_bins)

        pos_time_ccf = ccf[0:n_bins/2]
        neg_time_ccf = ccf[n_bins/2:]
        ccf = np.append(neg_time_ccf, pos_time_ccf)

        pos_time_ccf_err = ccf_err[0:n_bins/2]
        neg_time_ccf_err = pos_time_ccf_err[::-1]
        # neg_time_ccf_err = ccf_err[n_bins/2:]

        print( neg_time_ccf_err )
        print( pos_time_ccf_err )
        ccf_err = np.append(neg_time_ccf_err, pos_time_ccf_err)

        if i < 10:
            plot_file = args.out_root + "_chan_" + str(0) + str(i) + "." + args.plot_ext
        else:
            plot_file = args.out_root + "_chan_" + str(i) + "." + args.plot_ext


        make_plot(time_bins, ccf, ccf_err, args.prefix, plot_file, i, frac_time)

    # else:

        # raise Exception('ERROR: File type not recognized. Must have extension .dat or .fits.')

################################################################################
