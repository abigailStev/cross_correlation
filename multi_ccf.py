import argparse
import numpy as np
from scipy import fftpack
from datetime import datetime

import tools
from ccf import make_crossspec, filter_freq

"""
		multi_ccf.py

Computes the cross-correlation function of a band of interest with a reference
band, over multiple data files.

Arguments:
data_file_list - The full path of the (ASCII/txt/dat) input file listing the
    event lists to be used. One file per line. Assuming that both PCU0 and PCU2
    are in the event list.
out_file - Name of output file for cross-correlation function.
num_seconds - Number of seconds in a segment (must be a power of 2).
dt_mult - Multiple of 1/8192 seconds for timestep between bins.
test - 1 if only computing one segment for testing, 0 if computing all segments.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

'tools.py' is available at https://github.com/abigailStev/whizzy_scripts

"""

###############################################################################
def multi_output(out_file, in_file_list, dt, n_bins,
        total_duration, mean_rate_total_ci, mean_rate_total_ref, t,
        ccf_filtered, ccf_error):
    """
			multi_output
			
	Writes the cross-correlation function to an output file.
		
	Passed: out_file - Name of output file for cross-correlation function.
			data_file_list - Name of file containing lists of files with count
			    rate input data.
			reference_band - Name of file with the reference band count rate,
			    n-bins long.
			dt - Size of time bin, in seconds (must be power of 2).
			n_bins - Number of time bins in a segment (must be power of 2).
			total_duration - Total duration of all light curves used for primary
			    light curve.
			mean_rate_total_ci -  Mean count rate of all curves used for
			    channels of interest.
			mean_rate_total_ref - Mean count rate of all curves used for
			    reference band.
			t - Integer 'time' bins for CCF.
			ccf_filtered - Frequency-filtered CCF
			ccf_error - Error on frequency-filtered CCF.
			
	Returns: nothing
		
	"""

    print "Output sent to %s" % out_file

    with open(out_file, 'w') as out:
        out.write("#\t\tCross correlation function of multiple data files")
        out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
        out.write("\n# List of event lists: %s" % in_file_list)
        out.write("\n# Time bin size = %.21f seconds" % dt)
        out.write("\n# Number of bins per segment = %d" % n_bins)
        out.write("\n# Total duration of light curve = %d seconds" % \
                  total_duration)
        out.write("\n# Mean filtered count rate, curve 1 (per chan)= %s" % \
                  str(list(mean_rate_total_ci)))
        out.write("\n# Sum of mean filtered count rates, curve 1 = %.4f" % \
                  np.sum(mean_rate_total_ci))
        out.write("\n# Mean count rate, curve 2 = %.3f" % \
                  np.mean(mean_rate_total_ref))
        # 		if filter:
        # 			out.write("\n# Filter applied in frequency domain to eliminate \
        #               excess noise.")
        out.write("\n# ")
        out.write("\n# Column 1: Integer time bins")
        out.write("\n# Column 2-65: Filtered CCF per energy channel, real part")
        out.write("\n# Column 66-129: Error on filtered ccf per energy \
            channel, real part")
        out.write("\n# ")
        for j in xrange(0, n_bins):
            out.write("\n%d" % t[j])
            for i in xrange(0, 64):
                out.write("\t%.5f" % ccf_filtered[j][i].real)
            for i in xrange(0, 64):
                out.write("\t%.5f" % ccf_error[j][i].real)

        ## End of for-loops
    ## End of with-block


## End of function 'output'


###############################################################################
def main(in_file_list, out_file, num_seconds, dt_mult, test):
    """
			main
		
	Reads in a FITS file, takes FFT of data, makes power spectrum, writes to a
	file.
	
	Passed: in_file_list - Name of the file with a list of data files. 
				One file per line.
			out_file - Name of output file for standard power spectrum.
			num_seconds - Number of seconds for each segment of the light curve.
				Must be a power of 2.
			dt_mult - Multiple of 1/8192 seconds for timestep between bins.
			test - True if computing one segment, False if computing all.
			
		future option:
			filter - True if filtering the cross spectrum above and below 401Hz,
				False if leaving other frequencies untouched.
	
	Returns: nothing
	
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
    sum_rate_total_ci = np.zeros(64)
    sum_rate_total_ref = 0
    mean_rate_total_ci = np.zeros(64)
    mean_rate_total_ref = 0
    ccf = np.zeros((n_bins, 64))
    ccf_filtered = np.zeros((n_bins, 64))
    total_cs_sum = np.zeros((n_bins, 64), dtype=np.complex128)
    cs_avg = np.zeros((n_bins, 64), dtype=np.complex128)
    total_sum_power_ref = 0
    total_sum_power_ci = 0

    ## Looping through all data files
    for in_file in input_files:
        cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
            sum_power_ci, sum_power_ref = make_crossspec(in_file, n_bins,
            dt, test)
        total_segments += num_segments
        total_cs_sum += cs_sum
        sum_rate_total_ci += sum_rate_whole_ci
        sum_rate_total_ref += sum_rate_whole_ref
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

    ## Dividing these (currently just a sum of the segments) by the number of
    ## segments to get an arithmetic average
    mean_rate_total_ci = sum_rate_total_ci / float(total_segments)
    mean_rate_total_ref = sum_rate_total_ref / float(total_segments)
    mean_power_ci = total_sum_power_ci / float(total_segments)
    mean_power_ref = total_sum_power_ref / float(total_segments)
    cs_avg = total_cs_sum / float(total_segments)

    filtered_cs_avg, j_min, j_max = filter_freq(cs_avg, dt, n_bins, 401.0)
    assert np.shape(filtered_cs_avg) == np.shape(cs_avg)

    ## Absolute rms norms of poisson noise
    noise_ci = 2.0 * mean_rate_total_ci
    noise_ref = 2.0 * mean_rate_total_ref
    # 	print np.shape(noise_ci)
    # 	print np.shape(noise_ref)
    noise_ref_array = np.repeat(noise_ref, 64)

    old_settings = np.seterr(divide='ignore')

    df = 1.0 / float(num_seconds)  # in Hz
#     print "df =", df

    signal_ci_pow = np.complex128(mean_power_ci[j_min:j_max, :])
    signal_ref_pow = np.complex128(mean_power_ref[j_min:j_max])

    signal_ref_pow_stacked = signal_ref_pow
    for i in xrange(63):
        signal_ref_pow_stacked = np.column_stack(
            (signal_ref_pow, signal_ref_pow_stacked))

    assert np.shape(signal_ref_pow_stacked) == np.shape(signal_ci_pow)

    ## Putting powers into absolute rms2 normalization
    signal_ci_pow = signal_ci_pow * (2.0 * dt / float(n_bins)) - noise_ci
    print "signal ci pow:", signal_ci_pow[:, 0:4]
    signal_ref_pow = signal_ref_pow_stacked * (
        2.0 * dt / float(n_bins)) - noise_ref_array
    print "signal ref pow:", signal_ref_pow[:, 0:4]

    ## Getting rms of reference band, to normalize the ccf
    signal_variance = np.sum(signal_ref_pow * df)
    rms_ref = np.sqrt(
        signal_variance)  # should be a few percent in fractional rms units
    print "RMS of reference band:", rms_ref

    temp = np.square(noise_ci * signal_ref_pow) + \
           np.square(noise_ref * signal_ci_pow) + \
           np.square(noise_ci * noise_ref)
    # 	print "Shape of temp:", np.shape(temp)
    cs_noise_amp = np.sqrt(np.sum(temp) / float(total_segments)) * df
    # Might be off by a factor of 2 here...

#     cs_signal_amp = np.sum(cs_avg[j_min:j_max, :], axis=0)
    cs_signal_amp = np.sqrt(np.sum(cs_avg[j_min:j_max, :], axis=0) / float(total_segments)) * df

    print "Sum of cs signal amp:", np.sum(cs_signal_amp)

    # 	temp2 = np.sqrt(np.square(signal_ref_pow * signal_ci_pow)) * df
    # 	print "shape of temp2:", np.shape(temp2)
    # 	cs_signal_amp = np.sqrt(np.sum(temp2, axis=0) / float(num_segments))
    print "CS signal amp:", cs_signal_amp[0:4]
    print "shape of cs signal amp:", np.shape(cs_signal_amp)
    print "CS noise amp:", cs_noise_amp
    cs_error_ratio = cs_signal_amp / cs_noise_amp
    print "Shape cs error ratio:", np.shape(cs_error_ratio)
    # 	print "cs error ratio:", cs_error_ratio

    error_ratio_sigtop = cs_signal_amp / cs_noise_amp
    error_ratio_noisetop = cs_noise_amp / cs_signal_amp
    error_ratio_noisetop[10] = np.complex128(0)  
    # because that's the bin with no signal

    print "error ratio, signal on top:", error_ratio_sigtop
    print "error ratio, noise on top:", error_ratio_noisetop

    ## Taking the IFFT of the cross spectrum to get the CCF
    ccf = fftpack.ifft(cs_avg, axis=0)
    ccf_filtered = fftpack.ifft(filtered_cs_avg, axis=0)
    assert np.shape(ccf) == np.shape(ccf_filtered)

    ccf_error = np.absolute(error_ratio_noisetop) * np.absolute(ccf_filtered)

    ## Dividing ccf by integrated rms power of signal in reference band
    # ccf *= (2.0 / float(n_bins) / rms_ref)
    # 	ccf_filtered *= (2.0 / float(n_bins) / rms_ref)
    # 	ccf_error *= (2.0 / float(n_bins) / rms_ref)
    ccf /= rms_ref
    ccf_filtered /= rms_ref
    ccf_error /= rms_ref

    # 	other_ccf_error = np.abs(cs_error_ratio) * np.abs(ccf_filtered)


    print "unfilt_ccf sum =", np.sum(ccf)
    print "filt_ccf sum=", np.sum(ccf_filtered)

    print "filt CCF:", ccf_filtered[0, 0:4]
    print "Shape of ccf error:", np.shape(ccf_error)
    print "CCF error:", ccf_error[0, 0:4]
    # 	print "Other ccf error:", other_ccf_error[0,:]

    exposure = total_segments * num_seconds  # Exposure time of data used
    print "Exposure_time = %.3f seconds" % exposure
    print "Total number of segments:", total_segments
    # 	print "Total mean rate for ci:", mean_rate_total_ci
    # 	print "Mean mean rate for ci:", np.mean(mean_rate_total_ci)
    print "Sum of mean rate for ci:", np.sum(mean_rate_total_ci)
    print "Mean rate for ref:", mean_rate_total_ref

    t = np.arange(0, n_bins)  # gives the 'front of the bin'
    # time = t * dt  # Converting to seconds

    multi_output(out_file, in_file_list, dt, n_bins, exposure,
        mean_rate_total_ci, mean_rate_total_ref, t, ccf_filtered, ccf_error)

## End of the function 'main'



###############################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Computes the cross-\
        correlation function of a channel of interest with a reference band.')
    parser.add_argument('-i', '--infile_list', required=True,
        dest='infile_list', help="The full path of the (ASCII/txt/dat) input \
        file listing the event lists to be used. One file per line. Assuming \
        that both PCU0 and PCU2 are in the event list.")
    parser.add_argument('-o', '--outfile', required=True, dest='outfile',
        help="The full path of the (ASCII/txt/dat) output file to write the \
        cross-correlation function to.")
    parser.add_argument('-n', '--num_seconds', type=int, default=1,
        dest='num_seconds', help="Duration of each segment the light curve is \
        broken up into, in seconds. Must be an integer power of two.")
    parser.add_argument('-m', '--dt_mult', type=int, default=1, dest='dt_mult',
        help="Multiple of 1/8192 seconds for timestep between bins.")
    # parser.add_argument('-f', '--filter', help="True if applying filter \
    #   above and below 401Hz, False if not filtering in frequency.")
    parser.add_argument('-t', '--test', type=int, default=0,
        choices=xrange(0, 2), dest='test', help="1 if computing 1 segment for \
        testing, 0 if computing all segments.")
    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True

    main(args.infile_list, args.outfile, args.num_seconds, args.dt_mult, test)

## End of the program 'multi_CCF.py'
