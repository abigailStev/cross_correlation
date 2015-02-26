import argparse
import numpy as np
import sys
from scipy import fftpack
from datetime import datetime
import os
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
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
plot_file] [-p prefix]", description="Plots the ccf of multiple energy channels\
 in a 2D colour plot.", epilog="For optional arguments, default values are \
given in brackets at end of description.")

	parser.add_argument('tab_file', help="The table file, in .dat or .fits \
format.")

	parser.add_argument('-o', '--outfile', dest='plot_file', default="\
./ccf_2D.png", help="The file name to save the plot to. [./ccf_2D.png]")
	
	parser.add_argument('-p', '--prefix', dest='prefix', default="--", \
help="The identifying prefix of the data (object nickname or proposal ID). \
[--]")
	
	args = parser.parse_args()


	print "Plotting the 2D CCF: %s" % args.plot_file
	assert args.tab_file[-4:].lower() == "fits", "ERROR: Data file must be in \
	.fits format."
	

	#############################
	## Load data from FITS table
	#############################
	try:
		file_hdu = fits.open(args.tab_file)
	except IOError:
		print "\tERROR: File does not exist: %s" % args.tab_file
		sys.exit()
	
	n_bins = int(file_hdu[0].header["N_BINS"])
	means_str = file_hdu[0].header["RATE_CI"]
	detchans = int(file_hdu[0].header["DETCHANS"])	
	table = file_hdu[1].data
	file_hdu.close()
	
	mean_rate_ci = [ float(num.strip()) for num in means_str[1:-1].split(',') ]
	ccf = np.reshape(table.field('CCF'), (n_bins, detchans), order='C')

	t_length = 70  ## Number of time bins to use
	ccf = ccf[0:t_length,].T  ## Transpose it to get the axes we want
	
	###################################################################
	## Make a ratio of ccf to the mean count rate in the interest band
	###################################################################
	
	a = np.array([mean_rate_ci,]*t_length).T
	with np.errstate(all='ignore'):
		ratio = np.where(a != 0, ccf / a, 0)
	
	print "\tMinimum value:", np.min(ratio)
	print "\tMaximum value:", np.max(ratio)
	
	######################################################
	## Saving to a dat file so that we can use fimgcreate
	######################################################
	
# 	ratio[np.where(np.abs(ratio) > 0.3)] = 0
	print np.shape(ratio)
	ratio[32:,] = 0
	out_file = "out_ccf/temp.dat"
	R = ratio.flatten('C')
	comment_str = "From %s" % args.tab_file
	np.savetxt(out_file, R, comments=comment_str)
	
	#############
	## Plotting!
	#############
	
	font_prop = font_manager.FontProperties(size=16)
	fig, ax = plt.subplots(1,1)
	plt.pcolor(ratio, cmap='hot')  ## in future, can use vmin and vmax to set limits
	plt.colorbar()
	ax.set_xlabel('Arbitrary time bins', fontproperties=font_prop)
	ax.set_ylabel('Energy channel', fontproperties=font_prop)
	ax.set_ylim(3, 31)
# 	ax.set_ylim(0,detchans)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)
	ax.set_title(args.prefix)

	fig.set_tight_layout(True)
	plt.savefig(args.plot_file, dpi=150)
# 	plt.show()
	plt.close()		
	
## End of program 'plot_2d.py'

################################################################################