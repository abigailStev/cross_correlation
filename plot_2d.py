import argparse
import numpy as np
from scipy import fftpack
from datetime import datetime
from astropy.io import fits
import os.path
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.ticker import MultipleLocator

import tools

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens at uva.nl"
__year__ = "2014-2015"
__description__ = "Plots the ccf of multiple energy channels in a 2D colour \
plot."

"""
        plot_2d.py

Written in Python 2.7.

"""

################################################################################
if __name__ == "__main__":
	
	###########################
	## Parsing input arguments
	###########################
	
	parser = argparse.ArgumentParser(usage="python plot_2d.py tab_file [-o \
plot_file] [-p prefix] [-e chan_to_en]", description="Plots the ccf of multiple\
 energy channels in a 2D colour plot.", epilog="For optional arguments, default\
 values are given in brackets at end of description.")

	parser.add_argument('tab_file', help="The table file, in .dat or .fits \
format.")

	parser.add_argument('-o', '--outfile', dest='plot_file', default="\
./ccf_2D.png", help="The file name to save the plot to. [./ccf_2D.png]")
	
	parser.add_argument('-p', '--prefix', dest='prefix', default="--", \
help="The identifying prefix of the data (object nickname or proposal ID). \
[--]")
	
	parser.add_argument('-l', '--length', dest='t_length', \
type=tools.type_positive_int, default=100, help="Number of time bins to use \
along the x-axis. [100]")
	
	parser.add_argument('-e', dest='chan_to_en', help='Table of actual energy \
boundaries for energy channels. If not given, plot will be vs energy channel. \
[no default]')
	
	args = parser.parse_args()


	print "Plotting the 2D CCF: %s" % args.plot_file
	
	## Idiot checks
	assert args.tab_file[-4:].lower() == "fits", "ERROR: Data file must be in \
.fits format."
	if args.chan_to_en != None:  ## If args.chan_to_en exists as a variable
		assert os.path.isfile(args.chan_to_en), "ERROR: Table of real energy \
per energy channel does not exist."
	
	
	##########################################
	## Load data from FITS table and txt file
	##########################################
	
	if args.chan_to_en != None:  ## If args.chan_to_en exists as a variable
		energies = np.loadtxt(args.chan_to_en)	

	try:
		file_hdu = fits.open(args.tab_file)
	except IOError:
		print "\tERROR: File does not exist: %s" % args.tab_file
		exit()
	
	n_bins = int(file_hdu[0].header["N_BINS"])
	means_str = file_hdu[0].header["RATE_CI"]
	detchans = int(file_hdu[0].header["DETCHANS"])	
	dt = float(file_hdu[0].header["DT"])
	table = file_hdu[1].data
	file_hdu.close()
	
	mean_rate_ci = [ float(num.strip()) for num in means_str[1:-1].split(',') ]
	ccf = np.reshape(table.field('CCF'), (n_bins, detchans), order='C')
	
	ccf = ccf[0:args.t_length,].T  ## Transpose it to get the axes we want
	time_bins = np.arange(n_bins)
	frac_time = int(1.0/dt)  ## each time bin represents 1/frac_time sec

	###################################################################
	## Make a ratio of ccf to the mean count rate in the interest band
	###################################################################
	
	mean_ccf = np.mean(ccf, axis=0)
	ccf_resid = ccf - mean_ccf

	a = np.array([mean_rate_ci,]*args.t_length).T
	with np.errstate(all='ignore'):
		ratio = np.where(a != 0, ccf / a, 0)
	
# 	print "\tMinimum value:", np.min(ratio)
# 	print "\tMaximum value:", np.max(ratio)
	
	######################################################
	## Saving to a dat file so that we can use fimgcreate
	######################################################
	
# 	ratio[np.where(np.abs(ratio) > 0.75)] = 0
	ratio[32:,] = 0
	out_file = "./temp.dat"
	R = ratio.flatten('C')
	comment_str = "From %s" % args.tab_file
	np.savetxt(out_file, R, comments=comment_str)
	
	#############
	## Plotting!
	#############

	t_bins = np.arange(args.t_length+1)
	font_prop = font_manager.FontProperties(size=16)
	fig, ax = plt.subplots(1,1)

	if energies[0] != None:  ## If energies exists as a variable
		plt.pcolor(t_bins, energies, ratio, cmap='hot')
	else:
# 		ax.pcolor(ratio, cmap='hot', vmin=-1.2, vmax=1.6)
		plt.pcolor(ratio, cmap='hot')
# 		plt.pcolor(ccf_resid, cmap='hot')
	
	cbar = plt.colorbar()
	cbar.set_label('Ratio of CCF to mean count rate', \
		fontproperties=font_prop)
	ax.set_xlabel(r'Time ( $\times\,\frac{1}{%d}\,$s)' % frac_time, \
		fontproperties=font_prop)

	if energies[0] != None:  ## If energies exists as a variable
		ax.set_ylabel('Energy (keV)', fontproperties=font_prop)
		ax.set_ylim(3,25)
	else:
		ax.set_ylabel('Energy channel', fontproperties=font_prop)
		ax.set_ylim(2, 31)
# 	ax.set_yscale('log')

	## Setting the axes' minor ticks. It's complicated.
	x_maj_loc = ax.get_xticks()
	y_maj_loc = ax.get_yticks()
	x_min_mult = 0.1 * (x_maj_loc[1] - x_maj_loc[0])
	y_min_mult = 0.2 * (y_maj_loc[1] - y_maj_loc[0])
	xLocator = MultipleLocator(x_min_mult)  ## loc of minor ticks on x-axis
	yLocator = MultipleLocator(y_min_mult)  ## loc of minor ticks on y-axis
	ax.xaxis.set_minor_locator(xLocator)
	ax.yaxis.set_minor_locator(yLocator)

	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
# 	ax.set_title(args.prefix)

	fig.set_tight_layout(True)
	plt.savefig(args.plot_file, dpi=600)
# 	plt.show()
	plt.close()		
	
## End of program 'plot_2d.py'

################################################################################