import argparse
import numpy as np
import matplotlib.pyplot as plt

"""
		plot_ccf.py

Plots the cross-correlation function.

tab_file - Name of file with the output from ccf.py.
plot_root - Root of filename which the plots will be saved to in this script.
propID - The proposal ID of the data.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

"""

###############################################################################
def main(tab_file, plot_root, propID):
	pass
	
	print "\nPlotting the ccf: %s_chan_xx.png" % plot_root

	table = np.loadtxt(tab_file, comments='#') # table of CCF
	phase_bins = table[:,0] 
	n_bins = len(phase_bins)
# 	print "n_bins = ", n_bins
	CCF = [0 for x in range(0, 64)]
# 	exit()
	
# 	for i in range(0, 64):
	for i in range(6, 7):
		CCF[i] = table[:, i+1] 
		c = i
		if c >= 64: 
			c -= 64
		if c < 10:
			plot_file = plot_root + "_chan_" + str(0) + str(c) + ".png"
		else:
			plot_file = plot_root + "_chan_" + str(c) + ".png"
			
# 		print "Band", i
# 		print "Mean =", np.mean(CCF[i])
# 		print "Max =", np.max(CCF[i])

		fig, ax = plt.subplots()
		ax.plot(phase_bins, CCF[i], linewidth=1.5, label="Filtered CCF")
		plt.xlabel('Phase bins')
		plt.ylabel('Photon count rate [photons / s]')
		plt.xlim(0, n_bins)
# 		plt.ylim(-100,100)
		#plt.xscale('symlog') # this works much better than 'log'
		#plt.yscale('symlog')
		title = propID + ", Energy channel " + str(c)
		plt.title(title)

		## The following legend code was found on stack overflow
		##  or a pyplot tutorial
# 		legend = ax.legend(loc='lower right')
# 		## Set the fontsize
# 		for label in legend.get_texts():
# 			label.set_fontsize('small')
# 		for label in legend.get_lines():
# 			label.set_linewidth(2)  # the legend line width

		plt.savefig(plot_file, dpi=80)
# 		plt.show()
		plt.close()


###############################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--intable', required=True, dest='tab_file',
        help="The input table file, in ASCII/txt/dat format, with frequency in \
        column 1, CCF amplitudes in columns 2-65, error in columns 66-129.")
	parser.add_argument('-o', '--outroot', required=True, dest='plot_root',
        help="The root of the filename to save the plot to. Energy channel will\
         be appended to name before saving.")
	parser.add_argument('-p', '--propID', required=True, dest='propID',
        help="The proposal ID of the data.")
	args = parser.parse_args()

	main(args.tab_file, args.plot_root, args.propID)
