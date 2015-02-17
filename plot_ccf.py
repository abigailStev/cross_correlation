import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.font_manager as font_manager

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2014"
__description__ = "Plots the cross-correlation per energy channel, in the time domain."

"""
		plot_ccf.py

Enter   python plot_ccf.py -h   at the command line for help.

Written in Python 2.7.

All scientific modules imported above, as well as python 2.7, can be downloaded
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

"""

################################################################################
def make_plot(x_bins, ccf_amps, ccf_err, prefix, plot_file, chan):
	"""
			make_plot
		
	Actually makes the plot.
	
	"""
	
	font_prop = font_manager.FontProperties(size=16)

	fig, ax = plt.subplots(1,1)
	ax.plot(x_bins, ccf_amps, lw=2, c='black')
# 	ax.errorbar(x_bins, ccf_amps, yerr=ccf_err, lw=2, c='black', elinewidth=2, capsize=2)
	ax.set_xlabel('Arbitrary time bins', fontproperties=font_prop)
	ax.set_ylabel('Deviation from mean count rate [photons / s]', fontproperties=font_prop)
	ax.set_xlim(0, 2500)
# 	ax.set_ylim(-0.45, 0.45)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
# 	title = prefix + ", Energy channel " + str(chan)
# 	ax.set_title("Periodic signal at ~7 keV", fontproperties=font_prop)
	ax.set_title(prefix + ", Energy channel " + str(chan), fontproperties=font_prop)

	## The following legend code was found on stack overflow
	##  or a pyplot tutorial
# 	legend = ax.legend(loc='lower right')
# 	## Set the fontsize
# 	for label in legend.get_texts():
# 		label.set_fontsize('small')
# 	for label in legend.get_lines():
# 		label.set_linewidth(2)  # the legend line width

	fig.set_tight_layout(True)
	plt.savefig(plot_file, dpi=120)
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

	parser.add_argument('-p', '--prefix', dest='prefix', default="x", \
help="The identifying prefix of the data (object nickname or proposal ID). [x]")

	args = parser.parse_args()


	print "Plotting the ccf: %s_chan_xx.png" % args.out_root
	
	if args.tab_file[-4:].lower() == ".dat":
	
		table = np.loadtxt(args.tab_file, comments='#')  ## table of CCF
		time_bins = table[:,0] 
		n_bins = len(time_bins)
		ccf = np.zeros((n_bins, 64))
		ccf_err = np.zeros((n_bins, 64))
		for i in range(11, 12):
			ccf[:,i] = table[:, i+1] 
			ccf_err[:,i] = table[:, i+65]
			c = i
			if c >= 64: 
				c -= 64
			if c < 10:
				plot_file = args.out_root + "_chan_" + str(0) + str(c) + ".png"
			else:
				plot_file = args.out_root + "_chan_" + str(c) + ".png"
	
			make_plot(time_bins, ccf[:, i], ccf_err[:, i], args.prefix, plot_file, c)
			
	elif args.tab_file[-5:].lower() == ".fits":
	
		file_hdu = fits.open(args.tab_file)
		table = file_hdu[1].data
		file_hdu.close()
		
		for i in range(6, 7):
			channel_mask = table.field('CHANNEL') == i  ## make data mask for the energy channels i want
			table_i = table[channel_mask]
			ccf = table_i.field('CCF')
			ccf_err = table_i.field('ERROR')
			time_bins = table_i.field('TIME_BIN')

			if i < 10:
				plot_file = args.out_root + "_chan_" + str(0) + str(i) + ".png"
			else:
				plot_file = args.out_root + "_chan_" + str(i) + ".png"
	
			make_plot(time_bins, ccf, ccf_err, args.prefix, plot_file, i)
	
	else:
	
		raise Exception('ERROR: File type not recognized. Must have extension .dat or .fits.')
		
## End of function '__main__'
		
## End of program 'plot_ccf.py'
################################################################################

	