import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib.font_manager as font_manager

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2014-2015"
__description__ = "Plots CCFs of multiple energy channels together on one plot."

"""
		plot_multi.py

Written in Python 2.7.

"""

################################################################################
def main(file, plot_file, prefix):
	"""
			main
	
	
	"""
	
	# table_1 = np.loadtxt(file_1, comments='#')
	# table_2 = np.loadtxt(file_2, comments='#')
	# table_3 = np.loadtxt(file_3, comments='#')
	# table_4 = np.loadtxt(file_4, comments='#')
	# table_5 = np.loadtxt(file_5, comments='#')
	# table_6 = np.loadtxt(file_6, comments='#')
	
	######################################
	## Read in data from dat or fits file
	#######################################
	
	if file[-3:].lower() == "dat":

		table = np.loadtxt(file, comments='#')
		
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
	
		try:
			file_hdu = fits.open(file)
		except IOError:
			print "\tERROR: File does not exist: %s" % file
			exit()
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
		raise Exception("ERROR: File type not recognized. Must have extension \
.dat or .fits.")
	
	#############
	## Plotting!
	#############
	
	font_prop = font_manager.FontProperties(size=16)
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

	ax.set_xlabel('Arbitrary time bins', fontproperties=font_prop)
	ax.set_ylabel('Deviation from mean (photons / s)', \
		fontproperties=font_prop)
	ax.set_xlim(0, 4000)
# 	ax.set_xlim(0, 200)
# 	ax.set_ylim(-0.0005, 0.0005)
# 	ax.set_ylim(-2.5, 3)
# 	ax.set_xscale('symlog')  ## this works much better than 'log'
# 	ax.set_yscale('symlog')
	ax.set_title("%s, CCF per energy channel" % prefix, fontproperties=font_prop)

	## The following legend code was found on stack overflow I think
	legend = ax.legend(loc='upper right')
	for label in legend.get_texts():
		label.set_fontsize('small')
	for label in legend.get_lines():
		label.set_linewidth(2)  # the legend line width
	
	fig.set_tight_layout(True)
	plt.savefig(plot_file, dpi=140)
# 	plt.show()
	plt.close()

## End of function 'main'


################################################################################
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser(usage="python plot_multi.py ccf_table plot_file [-p prefix]", description="Plots CCFs of multiple energy channels on one \
plot.")
	
	parser.add_argument('ccf_table', help="Name of file with CCF amplitudes in \
a table.")
	
	parser.add_argument('plot_file', help="The output file name for the plot.")
	
	parser.add_argument('-p', '--prefix', dest='prefix', default="--", help="The identifying prefix of the data (object nickname or proposal ID). [--]")

	args = parser.parse_args()

	main(args.ccf_table, args.plot_file, args.prefix)
	
################################################################################
