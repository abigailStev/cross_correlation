import argparse
import numpy as np
import os
# import sys
from scipy import fftpack
from datetime import datetime
from astropy.io import fits
import tools  # https://github.com/abigailStev/whizzy_scripts
import ccf as xcor

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens@uva.nl"
__year__ = "2014-2015"
__description__ = "Computes the cross-correlation function of a band of \
interest with a reference band, over multiple RXTE event-mode data files."

"""
		multi_ccf.py

Written in Python 2.7.

"""

################################################################################
def dat_out(out_file, in_file_list, bkgd_file, dt, n_bins, detchans, exposure, \
	total_segments, mean_rate_total_ci, mean_rate_total_ref, t, ccf, ccf_error,\
	filtering):
	"""
			dat_out
	
	Writes the cross-correlation function to a .dat output file.
	
	"""
	if out_file[-4:].lower() == "fits":
		out_file = out_file[:-4]+"dat"
		
	print "\nOutput sent to: %s" % out_file
	
	with open(out_file, 'w') as out:
		out.write("#\t\tCross correlation function of multiple data files")
		out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
		out.write("\n# List of event lists: %s" % in_file_list)
		out.write("\n# Background spectrum: %s" % bkgd_file)
		out.write("\n# Time bin size = %.21f seconds" % dt)
		out.write("\n# Number of bins per segment = %d" % n_bins)
		out.write("\n# DETCHANS = %d" % detchans)
		out.write("\n# Total exposure time = %d seconds" % exposure)
		out.write("\n# Total number of segments = %d " % total_segments)
		out.write("\n# Mean count rate of ci = %s" % \
			str(list(mean_rate_total_ci)))
		out.write("\n# Mean count rate of ref band = %.5f" % \
			mean_rate_total_ref)
		out.write("\n# Filter applied in frequency domain? %s" % str(filtering))
		out.write("\n# ")
		out.write("\n# Column 1: Time bins")
		out.write("\n# Column 2-65: CCF per energy channel [count rate]")
		out.write("\n# Column 66-129: Error on ccf per energy channel [count rate]")
		out.write("\n# ")
		for j in xrange(0, n_bins):
			out.write("\n%d" % t[j])
			for i in xrange(0, 64):
				out.write("\t%.6e" % ccf[j][i].real)
			if filtering:
				for i in xrange(0, 64):
					out.write("\t%.6e" % ccf_error[i].real)
			else:
				for i in xrange(0, 64):
					out.write("\t%.6e" % ccf_error[j][i].real)

        ## End of for-loops
    ## End of with-block
## End of function 'dat_out'


################################################################################
def fits_out(out_file, in_file_list, bkgd_file, dt, n_bins, detchans, exposure,\
	total_segments, mean_rate_total_ci, mean_rate_total_ref, t, ccf, ccf_error,\
	filtering):
    """
            fits_out

    Writes the cross-correlation function to a .fits output file.
    
    """
    print "\nOutput sent to: %s" % out_file
	
    chan = np.arange(0, detchans)
    energy_channels = np.tile(chan, len(t))
    if filtering:
    	ccf_error = np.tile(ccf_error, len(t))
    else:
    	ccf_error = ccf_error.real.flatten('C')
    time_bins = np.repeat(t, len(chan))
    assert len(energy_channels) == len(time_bins)

    ## Making FITS header (extension 0)
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross-correlation function of multiple data files")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file_list)
    prihdr.set('BKGD', bkgd_file)
    prihdr.set('DT', dt, "seconds")
    prihdr.set('N_BINS', n_bins, "time bins per segment")
    prihdr.set('SEGMENTS', total_segments, "segments, of all data")
    prihdr.set('EXPOSURE', exposure, "seconds, of all data")
    prihdr.set('DETCHANS', detchans, "Number of detector energy channels")
    prihdr.set('RATE_CI', str(mean_rate_total_ci.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_total_ref, "counts/second")
    prihdr.set('FILTER', str(filtering))
    prihdu = fits.PrimaryHDU(header=prihdr)
        
    ## Making FITS table (extension 1)
    col1 = fits.Column(name='TIME_BIN', format='K', array=time_bins)
    col2 = fits.Column(name='CCF', unit='Counts/second', format='D', \
    	array=ccf.real.flatten('C'))
    col3 = fits.Column(name='ERROR', unit='', format='D', \
    	array=ccf_error)
    col4 = fits.Column(name='CHANNEL', unit='', format='I', \
    	array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3, col4])
    tbhdu = fits.BinTableHDU.from_columns(cols)
        
    ## If the file already exists, remove it
    assert out_file[-4:].lower() == "fits", \
    	'ERROR: Output file must have extension ".fits".'
    if os.path.isfile(out_file):
    	os.remove(out_file)
    	
    ## Writing to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file)	
    
## End of function 'fits_out'


################################################################################
def main(in_file_list, out_file, bkgd_file, num_seconds, dt_mult, test, \
	filtering):
    """
			main
		
	Reads in multiple event lists, splits into two light curves, makes segments 
	and populates them to give them length n_bins, computes the cross spectrum 
	of each segment per energy channel and then averaged cross spectrum of all 
	the segments per energy channel, and then computes the cross-correlation
    function (ccf) per energy channel.
	
	"""
	
	#####################################################
    ## Idiot checks, to ensure that our assumptions hold
    #####################################################
    
    assert num_seconds > 0, "ERROR: num_seconds must be a positive integer."
    assert dt_mult >= 1, "ERROR: dt_mult must be a positive integer."
    
    ##############################################################
	## Getting the list of input files and putting them in a list
	##############################################################
	
    input_files = [line.strip() for line in open(in_file_list)]
    if not input_files:  ## If data_files is an empty list
    	raise Exception("ERROR: No files in the eventlist list.")
    old_settings = np.seterr(divide='ignore')
	
	###########################################################
    ## Initializations
    ## 'total' is over all data files (i.e., in multi_ccf.py)
    ## 'whole' is over one data file (i.e., in ccf.py)
	###########################################################
	
    t_res = float(tools.get_key_val(input_files[0], 0, 'TIMEDEL'))
    dt = dt_mult * t_res
    n_bins = num_seconds * int(1.0 / dt)
    detchans = float(tools.get_key_val(input_files[0], 0, 'DETCHANS'))
    nyquist_freq = 1.0 / (2.0 * dt)
    
    total_segments = 0
    total_cs_sum = np.zeros((n_bins, detchans), dtype=np.complex128)
    sum_rate_total_ci = np.zeros(detchans)
    sum_rate_total_ref = 0
    total_sum_power_ci = np.zeros((n_bins, detchans), dtype=np.float64)
    total_sum_power_ref = np.zeros(n_bins, dtype=np.float64)
    
    print "\nDT = %.15f" % dt
    print "N_bins = %d" % n_bins
    print "Nyquist freq = %f" % nyquist_freq
    print "Filtering?", filtering
    
	###################################################################
	## Reading in the background count rate from a background spectrum
	####################################################################
   
    if bkgd_file:
   		print "Using background spectrum: %s" % bkgd_file
		bkgd_rate = xcor.get_background(bkgd_file)
    else:
		bkgd_rate = np.zeros(detchans)
    print " "
	
	##################################
    ## Looping through all data files
    ##################################
    
    for in_file in input_files:
    
        cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
            sum_power_ci, sum_power_ref, sum_rate_ci = \
            xcor.read_and_use_segments(in_file, n_bins, detchans, dt, test)
            
        total_segments += num_segments
        total_cs_sum += cs_sum
        sum_rate_total_ci += sum_rate_whole_ci
        sum_rate_total_ref += sum_rate_whole_ref
        total_sum_power_ci += sum_power_ci
        total_sum_power_ref += sum_power_ref
        sum_rate_whole_ci = None
        sum_rate_whole_ref = None
        cs_sum = None
        num_segments = None
        sum_power_ci = None
        sum_power_ref = None
        
    ## End of for-loop
    print " "
	
	#########################################
    ## Turning sums over segments into means
	#########################################
	
    mean_ci = sum_rate_ci / float(total_segments)
    mean_rate_total_ci = sum_rate_total_ci / float(total_segments)
    mean_rate_total_ref = sum_rate_total_ref / float(total_segments)
    mean_power_ci = total_sum_power_ci / float(total_segments)
    mean_power_ref = total_sum_power_ref / float(total_segments)
    cs_avg = total_cs_sum / float(total_segments)
    
    ################################################################
    ## Printing the cross spectrum to a file, for plotting/checking
    ################################################################
    
    cs_out = np.column_stack((fftpack.fftfreq(n_bins, d=dt), cs_avg))
    np.savetxt('cs_avg.dat', cs_out)
    
    ##################################################################
    ## Subtracting the background count rate from the mean count rate
    ##################################################################
    
    mean_rate_total_ci -= bkgd_rate
    
    ## Need to use a background from ref pcu for the reference band...
#     ref_bkgd_rate = np.mean(bkgd_rate[2:26])
#     mean_rate_whole_ref -= ref_bkgd_rate    
#     print np.shape(mean_rate_whole_ci)
#     print np.shape(mean_rate_whole_ref)
    
    ######################
    ## Making lag spectra
    ######################
    
    xcor.make_lags(out_file, in_file_list, dt, n_bins, detchans, num_seconds, \
    	total_segments, mean_rate_total_ci, mean_rate_total_ref, cs_avg, \
    	mean_power_ci, mean_power_ref)
	
	##############################################
	## Computing ccf from cs, and computing error
	##############################################
    
    if filtering:
    	ccf_end, ccf_error = xcor.FILT_cs_to_ccf_w_err(cs_avg, dt, n_bins, \
    		detchans, num_seconds, total_segments, mean_rate_total_ci, \
    		mean_rate_total_ref, mean_power_ci, mean_power_ref, True)
    else:
    	ccf_end, ccf_error = xcor.UNFILT_cs_to_ccf_w_err(cs_avg, dt, n_bins, \
    		detchans, num_seconds, total_segments, mean_rate_total_ci, \
    		mean_rate_total_ref, mean_power_ci, mean_power_ref, True)
    	
    exposure = total_segments * num_seconds  ## Exposure time of data used
    print "Exposure_time = %.3f seconds" % exposure
    print "Total number of segments:", total_segments
    print "Mean rate for all of ci:", np.sum(mean_rate_total_ci)
    print "Mean rate for ref:", mean_rate_total_ref

    t = np.arange(0, n_bins)  ## gives the 'front of the bin'
    
    ##########
    ## Output
    ##########
    
#     dat_out(out_file, in_file_list, bkgd_file, dt, n_bins, detchans, exposure, \
#     	total_segments, mean_rate_total_ci, mean_rate_total_ref, t, ccf_end, \
#     	ccf_error, filtering)

    fits_out(out_file, in_file_list, bkgd_file, dt, n_bins, detchans, exposure,\
    	total_segments, mean_rate_total_ci, mean_rate_total_ref, t, ccf_end, \
    	ccf_error, filtering)
		
## End of the function 'main'


################################################################################
if __name__ == "__main__":
	
	##############################################
	## Parsing input arguments and calling 'main'
	##############################################
	
    parser = argparse.ArgumentParser(usage="python multi_ccf.py infile outfile \
[-b BKGD_SPECTRUM] [-n NUM_SECONDS] [-m DT_MULT] [-t {0,1}] [-f {0,1}]", \
description='Computes the cross-correlation function of a channel of interest \
with a reference band, over multiple RXTE eventlists.', epilog='For optional \
arguments, default values are given in brackets at end of description.')
    	
    parser.add_argument('infile_list', help="The full path of the (ASCII/txt/\
dat) input file listing the event lists to be used. One file per line. \
Assuming that both PCU0 and PCU2 are in the event list.")
        
    parser.add_argument('outfile', help="The full path of the (ASCII/txt/dat) \
output file to write the cross-correlation function to.")
        
    parser.add_argument('-b', '--bkgd', required=False, dest='bkgd_file', \
help="Name of the (pha/fits) background spectrum. [none]")
    	
    parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two, \
default=1, dest='num_seconds', help="Number of seconds in each Fourier segment.\
 Must be a power of 2, positive, integer. [1]")
        
    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two, \
default=1, dest='dt_mult', help="Multiple of dt (dt is from data file) for \
timestep between bins. Must be a power of 2, positive, integer. [1]")
        
    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1}, 
dest='test', help="Int flag: 0 if computing all segments, 1 if computing only \
one segment for testing. [0]")
    	
    parser.add_argument('-f', '--filter', type=int, default=0, choices={0,1},
dest='filter', help='Int flag: 0 if NOT applying a filter in frequency-space, \
1 if applying filter (around a coherent pulsation). [0]')
        
    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True
	
    filtering = False
    if args.filter == 1:
    	filtering = True
    	
    main(args.infile_list, args.outfile, args.bkgd_file, args.num_seconds, \
    	args.dt_mult, test, filtering)
		
## End of the program 'multi_ccf.py'
################################################################################
