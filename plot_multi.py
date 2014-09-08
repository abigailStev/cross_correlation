import argparse
import numpy as np
import matplotlib.pyplot as plt


##########################
def main(file, plot_file, numsec):
	pass
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
	data_15 = table[:,4]
	data_2 = table[:,7] # chan 6, abs 11
	data_25 = table[:,10]
	data_3 = table[:,12] # chan 11, abs 16
	data_4 = table[:,15] # chan 14, abs 23
	data_45 = table[:,17]
	data_5 = table[:,19] # chan 18, abs 31
	data_55 = table[:,20]
	data_6 = table[:,22] # chan 21, abs 48
	data_65 = table[:,24]
	data_7 = table[:,26]
	# chan is from 0 to 63 inclusive
	# abs is from 0 to 254 inclusive

	fig, ax = plt.subplots()
	ax.plot(bins, data_0, linewidth=2, label="Chan 0")
	ax.plot(bins, data_1, linewidth=2, label="Chan 1") #label="2.5 keV")
	ax.plot(bins, data_15, linewidth=2, label="Chan 3")
	ax.plot(bins, data_2, linewidth=2, label="Chan 6") #label="5 keV")
	ax.plot(bins, data_25, linewidth=2, label="Chan 9")
	ax.plot(bins, data_3, linewidth=2, label="Chan 11") #label="7 keV")
	ax.plot(bins, data_4, linewidth=2, label="Chan 14") #label="10 keV")
	ax.plot(bins, data_45, linewidth=2, ls='-.', label="Chan 16")
	ax.plot(bins, data_5, linewidth=2, ls='-.', label="Chan 18") #label="13 keV")
	ax.plot(bins, data_55, linewidth=2, ls='-.', label="Chan 19")
	ax.plot(bins, data_6, linewidth=2, ls='-.', label="Chan 21") #label="20 keV")
	ax.plot(bins, data_65, linewidth=2, ls='-.', label="Chan 23")
	ax.plot(bins, data_7, linewidth=2, ls='-.', label="Chan 25")

	plt.xlabel('Phase bins')
	plt.ylabel('Photon count rate [photons / s]')
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