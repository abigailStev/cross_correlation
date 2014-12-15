import argparse
import numpy as np
import os
# import sys
from scipy import fftpack
from datetime import datetime
from astropy.io import fits
import tools  # https://github.com/abigailStev/whizzy_scripts
from ccf import read_and_use_segments, filter_freq, cs_to_ccf_w_err, make_lags

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2014"
__description__ = "Computes the cross-correlation function of a band of \
interest with a reference band, over multiple RXTE event-mode data files."

"""
		multi_ccf.py

Written in Python 2.7.

"""

###############################################################################
def dat_out(out_file, in_file_list, dt, n_bins, total_exposure, 
	mean_rate_total_ci, mean_rate_total_ref, t, ccf, ccf_error):
	"""
	dat_out
	Writes the cross-correlation function to an output file.
	"""
	if out_file[-4:].lower() == "fits":
		out_file = out_file[:-4]+"dat"
	print "Output sent to: %s" % out_file
	with open(out_file, 'w') as out:
		out.write("#\t\tCross correlation function of multiple data files")
		out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
		out.write("\n# List of event lists: %s" % in_file_list)
		out.write("\n# Time bin size = %.21f seconds" % dt)
		out.write("\n# Number of bins per segment = %d" % n_bins)
		out.write("\n# Total exposure time = %d seconds" % total_exposure)
		out.write("\n# Mean count rate of ci = %s" % str(list(mean_rate_total_ci)))
		out.write("\n# Mean count rate of ref band = %.5f" % mean_rate_total_ref)
		out.write("\n# ")
		out.write("\n# Column 1: Time bins")
		out.write("\n# Column 2-65: CCF per energy channel [count rate]")
		out.write("\n# Column 66-129: Error on ccf per energy channel [count rate]")
		out.write("\n# ")
		for j in xrange(0, n_bins):
			out.write("\n%d" % t[j])
			for i in xrange(0, 64):
				out.write("\t%.6e" % ccf[j][i].real)
			for i in xrange(0, 64):
				out.write("\t%.6e" % ccf_error[i].real)

        ## End of for-loops
    ## End of with-block
## End of function 'dat_out'


###############################################################################
def fits_out(out_file, in_file_list, dt, n_bins, total_exposure, total_segments, 
	mean_rate_total_ci, mean_rate_total_ref, t, ccf, ccf_error):
    """
            fits_out

    Writes the cross-correlation function to a .fits output file.
    
    """
    print "Output sent to: %s" % out_file

    ## Making header
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross-correlation function of multiple data files")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file_list)
    prihdr.set('DT', dt, "seconds")
    prihdr.set('N_BINS', n_bins, "time bins per segment")
    prihdr.set('SEGMENTS', total_segments, "segments, of all data")
    prihdr.set('EXPOSURE', total_exposure, \
    	"seconds, of all data")
    prihdr.set('RATE_CI', str(mean_rate_total_ci.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_total_ref, "counts/second")
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    chan = np.arange(0,64)
    energy_channels = np.tile(chan, len(t))
    ccf_error = np.tile(ccf_error, len(t))
    time_bins = np.repeat(t, len(chan))
#     print len(energy_channels)
#     print len(time_bins)
    assert len(energy_channels) == len(time_bins)
    
    ## Making FITS table
    col1 = fits.Column(name='TIME_BIN', format='K', array=time_bins)
    col2 = fits.Column(name='CCF', unit='Counts/second', format='D', \
    	array=ccf.real.flatten('C'))
    col3 = fits.Column(name='ERROR', unit='', format='D', \
    	array=ccf_error)
    col4 = fits.Column(name='CHANNEL', unit='', format='I', \
    	array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3, col4])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    
    ## If the file already exists, remove it (still working on just updating it)
    assert out_file[-4:].lower() == "fits", \
    	'ERROR: Output file must have extension ".fits".'
    if os.path.isfile(out_file):
    	os.remove(out_file)
    ## Writing to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file)	
## End of function 'fits_out'


###############################################################################
def main(in_file_list, out_file, num_seconds, dt_mult, test):
    """
			main
		
	Reads in a FITS file, takes FFT of data, makes power spectrum, writes to a
	file.
	
	"""

    ## Idiot checks, to ensure that our assumptions hold
    assert num_seconds > 0  # num_seconds must be a positive integer
    assert tools.power_of_two(num_seconds)  # num_seconds must be a power of 2
    assert dt_mult >= 1

    input_files = [line.strip() for line in open(in_file_list)]

    t_res = 1.0 / 8192.0
    dt = dt_mult * t_res
    n_bins = num_seconds * int(1.0 / dt)
    print "dt = %.21f seconds" % dt
    print "n_bins = %d" % n_bins

    ## Initializations -- 'total' is over all data files,
    ## 'whole' is over one data file
    total_segments = 0
    sum_rate_total_ci = np.zeros(64, dtype=np.float64)
    sum_rate_total_ref = 0
    mean_rate_total_ci = np.zeros(64, dtype=np.float64)
    mean_rate_total_ref = 0
    ccf = np.zeros((n_bins, 64))
    ccf_filtered = np.zeros((n_bins, 64))
    total_cs_sum = np.zeros((n_bins, 64), dtype=np.complex128)
    cs_avg = np.zeros((n_bins, 64), dtype=np.complex128)
    total_sum_power_ref = 0
    total_sum_power_ci = 0
    sum_rate_ci = np.zeros(64)

    ## Looping through all data files
    for in_file in input_files:
        cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
            sum_power_ci, sum_power_ref, sum_rate_ci = read_and_use_segments(in_file, n_bins,
            dt, test)
            
        total_segments += num_segments
        total_cs_sum += cs_sum
        sum_rate_total_ci += sum_rate_whole_ci
        sum_rate_total_ref += sum_rate_whole_ref
#         print sum_rate_total_ci[19:24]
        total_sum_power_ci += sum_power_ci
        total_sum_power_ref += sum_power_ref
        # 		print "cs sum shape", np.shape(cs_sum)
        # 		print "total cs sum shape", np.shape(total_cs_sum)
        sum_rate_whole_ci = None
        sum_rate_whole_ref = None
        cs_sum = None
        num_segments = None
        sum_power_ci = None
        sum_power_ref = None
    ## End of for-loop
	
	mean_ci = sum_rate_ci / float(total_segments)
# 	print "Mean rate of ci: ", mean_ci
	
    ## Dividing these (currently just a sum of the segments) by the number of
    ## segments to get an arithmetic average
    mean_rate_total_ci = sum_rate_total_ci / float(total_segments)
    mean_rate_total_ref = sum_rate_total_ref / float(total_segments)
    mean_power_ci = total_sum_power_ci / float(total_segments)
    mean_power_ref = total_sum_power_ref / float(total_segments)
    cs_avg = total_cs_sum / float(total_segments)
    
    make_lags(out_file, in_file_list, dt, n_bins, num_seconds, total_segments, \
    	mean_rate_total_ci, mean_rate_total_ref, cs_avg, mean_power_ci, \
    	mean_power_ref)
	
    ccf_filtered, ccf_error = cs_to_ccf_w_err(cs_avg, dt, n_bins, num_seconds, \
    	total_segments, mean_rate_total_ci, mean_rate_total_ref, mean_power_ci,\
    	mean_power_ref, True)
    
    exposure = total_segments * num_seconds  # Exposure time of data used
    print "Exposure_time = %.3f seconds" % exposure
    print "Total number of segments:", total_segments
#     print "Total mean rate for ci:", mean_rate_total_ci
#     print "Mean rate for ci:", np.mean(mean_rate_total_ci) * 64
    print "Mean rate for all of ci:", np.sum(mean_rate_total_ci)
    print "Mean rate for ref:", mean_rate_total_ref

    t = np.arange(0, n_bins)  # gives the 'front of the bin'
    dat_out(out_file, in_file_list, dt, n_bins, exposure,
		mean_rate_total_ci, mean_rate_total_ref, t, ccf_filtered, ccf_error)
    fits_out(out_file, in_file_list, dt, n_bins, exposure, total_segments, 
		mean_rate_total_ci, mean_rate_total_ref, t, ccf_filtered, ccf_error)
## End of the function 'main'


###############################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Computes the cross-\
        correlation function of a channel of interest with a reference band, \
        over multiple RXTE event-mode data files.')
    parser.add_argument('-i', '--infile_list', required=True,
        dest='infile_list', help="The full path of the (ASCII/txt/dat) input \
        file listing the event lists to be used. One file per line. Assuming \
        that both PCU0 and PCU2 are in the event list.")
    parser.add_argument('-o', '--outfile', required=True, dest='outfile',
        help="The full path of the (ASCII/txt/dat) output file to write the \
        cross-correlation function to.")
    parser.add_argument('-n', '--num_seconds', type=int, default=1,
        dest='num_seconds', help="Number of seconds in each Fourier segment. \
        Must be an integer power of two.")
    parser.add_argument('-m', '--dt_mult', type=int, default=1, dest='dt_mult',
        help="Multiple of 1/8192 seconds for timestep between bins.")
    # parser.add_argument('--filter', action='store_true', help="If present, \
    # filter will be applied above and below 401 Hz.")
    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1}, 
    	dest='test', help="1 if computing 1 segment for testing, 0 if \
    	computing all segments.")
    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True

    main(args.infile_list, args.outfile, args.num_seconds, args.dt_mult, test)

## End of the program 'multi_ccf.py'
