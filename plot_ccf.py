import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

"""
		plot_ccf.py

Plots the cross-correlation amplitudes in the time domain.

Enter   python plot_ccf.py -h   at the command line for help.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

"""
###############################################################################
def make_plot(x_bins, ccf_amps, y_err, propID, plot_file, chan):
		n_bins = len(x_bins)
		x_err = np.zeros(n_bins)

		fig, ax = plt.subplots()
		ax.errorbar(x_bins, ccf_amps, xerr=x_err, yerr=y_err, lw=2, c='black')
		plt.xlabel('Arbitrary time bins')
		plt.ylabel('Deviation from mean count rate [photons / s]')
		plt.xlim(15, 50)
		plt.ylim(-0.45, 0.45)
		#plt.xscale('symlog') # this works much better than 'log'
		#plt.yscale('symlog')
		title = propID + ", Energy channel " + str(chan)
		title = "Periodic signal at ~7 keV"
		plt.title(title)

		## The following legend code was found on stack overflow
		##  or a pyplot tutorial
# 		legend = ax.legend(loc='lower right')
# 		## Set the fontsize
# 		for label in legend.get_texts():
# 			label.set_fontsize('small')
# 		for label in legend.get_lines():
# 			label.set_linewidth(2)  # the legend line width

		plt.savefig(plot_file, dpi=120)
# 		plt.show()
		plt.close()
		

###############################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('tab_file', help="The table file, in .dat or .fits \
		format.")
	parser.add_argument('-o', '--out_root', dest='out_root', default="./ccf", \
		help="The root of the filename to save the plot to. Energy channel will\
		be appended to name before saving.")
	parser.add_argument('-p', '--propID', dest='propID', default="Pxxxxx", \
		help="The proposal ID of the data.")
	args = parser.parse_args()
	print args.tab_file
	print "\nPlotting the ccf: %s_chan_xx.png" % args.out_root
	
	if args.tab_file[-4:].lower() == ".dat":
		table = np.loadtxt(args.tab_file, comments='#') # table of CCF
		time_bins = table[:,0] 
		n_bins = len(time_bins)
		CCF = np.zeros((n_bins, 64))
		y_err = np.zeros((n_bins, 64))
		for i in range(11, 12):
			CCF[:,i] = table[:, i+1] 
			y_err[:,i] = table[:, i+65]
			c = i
			if c >= 64: 
				c -= 64
			if c < 10:
				plot_file = args.out_root + "_chan_" + str(0) + str(c) + ".png"
			else:
				plot_file = args.out_root + "_chan_" + str(c) + ".png"
	
			make_plot(time_bins, CCF[:, i], y_err[:, i], args.propID, plot_file, c)
			
	elif args.tab_file[-5:].lower() == ".fits":
		file_hdu = fits.open(args.tab_file)
		table = file_hdu[1].data
		file_hdu.close()
		
		for i in range(6, 7):
			channel_mask = table.field('CHANNEL') == i 		# make data mask for the energy channels i want
			table_i = table[channel_mask]
			CCF = table_i.field('CCF')
			y_err = table_i.field('ERROR')
			time_bins = table_i.field('TIME_BIN')

			if i < 10:
				plot_file = args.out_root + "_chan_" + str(0) + str(i) + ".png"
			else:
				plot_file = args.out_root + "_chan_" + str(i) + ".png"
	
			make_plot(time_bins, CCF, y_err, args.propID, plot_file, i)
	else:
		raise Exception('ERROR: File type not recognized. Must have extension .dat or .fits.')
		
		
## End of program 'plot_ccf.py'
	