import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2014"
__description__ = "Plots multiple CCFs together on one plot."

"""
		plot_multi.py

Written in Python 2.7.

"""

##########################
def main(file, plot_file, numsec):
	# table_1 = np.loadtxt(file_1, comments='#')
	# table_2 = np.loadtxt(file_2, comments='#')
	# table_3 = np.loadtxt(file_3, comments='#')
	# table_4 = np.loadtxt(file_4, comments='#')
	# table_5 = np.loadtxt(file_5, comments='#')
	# table_6 = np.loadtxt(file_6, comments='#')
	# 
	# data_1 = table_1[:,0]
	# data_2 = table_2[:,0]
	# data_3 = table_3[:,0]
	# data_4 = table_4[:,0]
	# data_5 = table_5[:,0]
	# data_6 = table_6[:,0]
	if file[-3:].lower() == "dat":
# 	file="/Users/abigailstevens/Dropbox/Research/power_spectra/out_ccf/140610_t1_4sec.dat"
		table = np.loadtxt(file, comments='#')
	# 	data_1 = table[:,65] # chan 1, abs 5
	# 	data_2 = table[:,70] # chan 6, abs 11
	# 	data_3 = table[:,76] # chan 11, abs 16
	# 	data_4 = table[:,78] # chan 14, abs 23
	# 	data_5 = table[:,81] # chan 18, abs 31
	# 	data_6 = table[:,87] # chan 21, abs 48
		bins = table[:,0]
		data_0 = table[:,1]
		data_1 = table[:,2] # chan 1, abs 5
		data_3 = table[:,4]
		data_6 = table[:,7] # chan 6, abs 11
		data_9 = table[:,10]
		data_11 = table[:,12] # chan 11, abs 16
		data_14 = table[:,15] # chan 14, abs 23
		data_16 = table[:,17]
		data_18 = table[:,19] # chan 18, abs 31
		data_19 = table[:,20]
		data_21 = table[:,22] # chan 21, abs 48
		data_23 = table[:,24]
		data_25 = table[:,26]
		# chan is from 0 to 63 inclusive
		# abs is from 0 to 254 inclusive
	elif file[-4:].lower() == "fits":
		file_hdu = fits.open(file)
		table = file_hdu[1].data
		file_hdu.close()
		bins = table[table.field('CHANNEL') == 0].field('TIME_BIN')
		data_0 = table[table.field('CHANNEL') == 1].field('CCF')
		data_1 = table[table.field('CHANNEL') == 2].field('CCF') # chan 1, abs 5
		data_3 = table[table.field('CHANNEL') == 4].field('CCF')
		data_6 = table[table.field('CHANNEL') == 7].field('CCF') # chan 6, abs 11
		data_9 = table[table.field('CHANNEL') == 10].field('CCF') 
		data_11 = table[table.field('CHANNEL') == 12].field('CCF') # chan 11, abs 16
		data_14 = table[table.field('CHANNEL') == 15].field('CCF') # chan 14, abs 23
		data_16 = table[table.field('CHANNEL') == 17].field('CCF')
		data_18 = table[table.field('CHANNEL') == 19].field('CCF') # chan 18, abs 31
		data_19 = table[table.field('CHANNEL') == 20].field('CCF')
		data_21 = table[table.field('CHANNEL') == 22].field('CCF') # chan 21, abs 48
		data_23 = table[table.field('CHANNEL') == 24].field('CCF')
		data_25 = table[table.field('CHANNEL') == 26].field('CCF')
		
	else:
		raise Exception("ERROR: File type not recognized. Must have extension .dat or .fits.")
		
	fig, ax = plt.subplots()
	ax.plot(bins, data_0, linewidth=2, label="Chan 0")
	ax.plot(bins, data_1, linewidth=2, label="Chan 1") #label="2.5 keV")
	ax.plot(bins, data_3, linewidth=2, label="Chan 3")
	ax.plot(bins, data_6, linewidth=2, label="Chan 6") #label="5 keV")
	ax.plot(bins, data_9, linewidth=2, label="Chan 9")
	ax.plot(bins, data_11, linewidth=2, label="Chan 11") #label="7 keV")
	ax.plot(bins, data_14, linewidth=2, label="Chan 14") #label="10 keV")
	ax.plot(bins, data_16, linewidth=2, ls='-.', label="Chan 16")
	ax.plot(bins, data_18, linewidth=2, ls='-.', label="Chan 18") #label="13 keV")
	ax.plot(bins, data_19, linewidth=2, ls='-.', label="Chan 19")
	ax.plot(bins, data_21, linewidth=2, ls='-.', label="Chan 21") #label="20 keV")
	ax.plot(bins, data_23, linewidth=2, ls='-.', label="Chan 23")
	ax.plot(bins, data_25, linewidth=2, ls='-.', label="Chan 25")

	plt.xlabel('Arbitrary time bins')
	plt.ylabel('Deviation from mean count rate [photons / s]')
	# plt.xlim(0,20000)
	plt.xlim(15,60)
# 	plt.ylim(-0.0005,0.0005)
	# plt.ylim(-2.5,3)
	# plt.xscale('symlog') # this works much better than 'log'
	# plt.yscale('symlog')
	title_str = "CCF per energy channel"
	plt.title(title_str)

	## The following legend code was found on stack overflow I think, or a pyplot tutorial
	legend = ax.legend(loc='upper right')
	for label in legend.get_texts():
		label.set_fontsize('small')
	for label in legend.get_lines():
		label.set_linewidth(2)  # the legend line width

	plt.savefig(plot_file, dpi=140)
# 	plt.show()
	plt.close()


##########################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(description="Plots multiple cross-correlation functions together.")
	parser.add_argument('ccf_table', \
		help="Name of file with CCF amplitudes in a table.")
	parser.add_argument('plot_file', \
		help="The output file name for the plot.")
	parser.add_argument('numsec', type=int, help="")
	args = parser.parse_args()

	main(args.ccf_table, args.plot_file, args.numsec)