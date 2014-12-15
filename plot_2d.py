import argparse
import numpy as np
import sys
from scipy import fftpack
from datetime import datetime
import os
from astropy.io import fits
import matplotlib.pyplot as plt
from tools import get_key_val
import matplotlib.font_manager as font_manager

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2014"
__description__ = ""

"""
        plot_2d.py

Written in Python 2.7.

"""

###############################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('tab_file', help="The table file, in .dat or .fits \
		format.")
	parser.add_argument('-o', '--outfile', dest='plot_file', default="./ccf", \
		help="The file name to save the plot to.")
	args = parser.parse_args()


	print "\nPlotting the 2D CCF: %s" % args.plot_file
	print args.tab_file
	assert args.tab_file[-4:].lower() == "fits", "ERROR: Data file must be in .fits\
	format."
	
	n_bins = get_key_val(args.tab_file, 0, "N_BINS")
# 	print n_bins
	
	file_hdu = fits.open(args.tab_file)
	table = file_hdu[1].data
	file_hdu.close()

# 	ccf = np.reshape(table.field('CCF'), (64,n_bins), order='C')
# 	print np.shape(ccf)
# 	ccf = ccf[:,15:50]

	ccf = np.reshape(table.field('CCF'), (n_bins, 64), order='C')
	print ccf[0,0:4]  # printing the first time bin, energies 0 - 3 inclusive

	print np.shape(ccf)
	ccf = ccf[0:50,]
	
	font_prop = font_manager.FontProperties(size=16)

	fig, ax = plt.subplots(1,1)
	plt.pcolor(ccf, cmap='hot')
	plt.colorbar()
	ax.set_xlabel('Energy channel', fontproperties=font_prop)
	ax.set_ylabel('Arbitrary time bins', fontproperties=font_prop)
	ax.set_xlim(2, 26)
# 	ax.set_ylim(-0.45, 0.45)
	ax.tick_params(axis='x', labelsize=14)
	ax.tick_params(axis='y', labelsize=14)

	plt.savefig(args.plot_file, dpi=120)
# 	plt.show()
	plt.close()
		
		
		
		
		
		
		
		
		
