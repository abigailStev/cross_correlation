import argparse
import numpy as np
import sys
from scipy import fftpack
from datetime import datetime
import os
from astropy.io import fits

import tools
import warnings

"""
        ccf.py

Computes the cross-correlation function of two light curves, to do phase-
resolved spectroscopy. Able to read light curves from data with the program
'populate_lightcurve'.

Required arguments:
in_file - Name of (ASCII/txt/dat) input file event list containing both the
    reference band and the channels of interest. Assumes ref band = PCU 0,
    interest = PCU 2.
out_file - Name of (ASCII/txt/dat) output file which the table of
    cross-correlation function data will be written to.
    
Optional arguments:
num_seconds - Number of seconds in each Fourier segment. Must be a power of 2.
dt_mult - Multiple of 1/8192 seconds for timestep between bins.
test - 1 if computing one segment for testing, 0 if doing full run.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded
in the Anaconda package, https://store.continuum.io/cshop/anaconda/

'tools.py' is available in my public GitHub repo 'whizzy_scripts'

"""


###############################################################################
def dat_output(out_file, in_file, dt, n_bins, num_seconds, num_segments,
        mean_rate_whole_ci, mean_rate_whole_ref, t, ccf_filtered, ccf_error):
    """
            dat_output

    Writes the cross-correlation function to a .dat output file.

    Passed: out_file - Name of output file.
            in_file - Name of event list with both reference and interest bands.
            dt - Size of each time bin, in seconds.
            n_bins - Number of (time) bins per segment.
            num_seconds - Number of seconds in each Fourier segment.
            num_segments - Number of segments the light curve was split up into.
            mean_rate_whole_ci - Mean count rate of light curve 1, averaged over
                all segments.
            mean_rate_whole_ref - Mean count rate of light curve 2, averaged
                over all segments.
            t - Integer time bins to plot against the ccf.
            ccf_filtered - CCF amplitudes, filtered in frequency space.
            ccf_error - Error on the filtered CCF.

    Returns: nothing

    """
    print "Output sent to: %s" % out_file

    with open(out_file, 'w') as out:
        out.write("#\t\tCross-correlation function data")
        out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
        out.write("\n# Event list: %s" % in_file)
        out.write("\n# Time bin size = %.21f seconds" % dt)
        out.write("\n# Number of bins per segment = %d" % n_bins)
        out.write("\n# Number of seconds per segment = %d" % num_seconds)
        out.write("\n# Number of segments per light curve = %d" % num_segments)
        out.write("\n# Exposure time = %d seconds" \
                  % (num_segments * num_seconds))
        out.write("\n# Mean count rate of ci = %s" \
        	% str(list(mean_rate_whole_ci)))
        out.write("\n# Mean count rate of ref band = %.8f" \
                  % np.mean(mean_rate_whole_ref))
        out.write("\n# ")
        out.write("\n# Column 1: Time bins")
        out.write("\n# Columns 2-65: Filtered ccf per energy channel, real \
            part [count rate]")
        out.write("\n# Columns 66-129: Error on filtered ccf per energy \
            channel, real part [count rate]")
        out.write("\n# ")
        for j in xrange(0, n_bins):
            out.write("\n%d" % t[j])
            for i in xrange(0, 64):
                out.write("\t%.6e" % ccf_filtered[j][i].real)
            for i in xrange(0, 64):
                out.write("\t%.6e" % ccf_error[i].real)
        ## End of for-loops
    ## End of with-block
    
## End of function 'dat_output'


###############################################################################
def fits_output(out_file, in_file, dt, n_bins, num_seconds, num_segments,
        mean_rate_whole_ci, mean_rate_whole_ref, t, ccf_filtered, ccf_error):
    """
            fits_output

    Writes the cross-correlation function to a .fits output file.

    Passed: out_file - Name of output file.
            in_file - Name of event list with both reference and interest bands.
            dt - Size of each time bin, in seconds.
            n_bins - Number of (time) bins per segment.
            num_seconds - Number of seconds in each Fourier segment.
            num_segments - Number of segments the light curve was split up into.
            mean_rate_whole_ci - Mean count rate of light curve 1, averaged over
                all segments.
            mean_rate_whole_ref - Mean count rate of light curve 2, averaged
                over all segments.
            t - Integer time bins to plot against the ccf.
            ccf_filtered - CCF amplitudes, filtered in frequency space.
            ccf_error - Error on the filtered CCF.

    Returns: nothing
    """
    print "Output sent to: %s" % out_file

    ## Making header for standard power spectrum
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross-correlation data")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file)
    prihdr.set('DT', dt, "seconds")
    prihdr.set('N_BINS', n_bins, "time bins per segment")
    prihdr.set('SEGMENTS', num_segments, "segments in the whole light curve")
    prihdr.set('EXPOSURE', num_segments * n_bins * dt, \
    	"seconds, of light curve")
    prihdr.set('RATE_CI', str(mean_rate_whole_ci.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_whole_ref, "counts/second")
    prihdu = fits.PrimaryHDU(header=prihdr)
    
    chan = np.arange(0,65)
    energy_channels = np.tile(chan, len(t))
    time_bins = np.repeat(t, len(chan))
#     print len(energy_channels)
#     print len(time_bins)
    assert len(energy_channels) == len(time_bins)
    
    ## Making FITS table for standard power spectrum
    col1 = fits.Column(name='TIME_BIN', format='K', array=time_bins)
    col2 = fits.Column(name='CCF', unit='', format='D', \
    	array=ccf_filtered.real.flatten('C'))
    col3 = fits.Column(name='ERROR', unit='', format='D', \
    	array=ccf_error.flatten('C'))
    col4 = fits.Column(name='CHANNEL', unit='', format='I', \
    	array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3, col4])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    
    ## If the file already exists, remove it (still working on just updating it)
    assert out_file[-4:].lower() == "fits", 'ERROR: Output file must have extension ".fits".'
    if os.path.isfile(out_file):
    	os.remove(out_file)
    ## Writing the standard power spectrum to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file)	
    
	## End of function 'fits_output'
	

###############################################################################
def stack_reference_band(rate_ref_2d, obs_epoch):
    """
            stack_reference_band

    Stacks the photons in the reference band from 3-20 keV to make one 'band'.
    Assumes that we are in epoch 5 or 3.
    Should I be using the FTOOL rbnpha or grppha to assist with this?
    
    Epoch 1: abs channels 10 - 74
    Epoch 2: abs channels 9 - 62
    Epoch 3: abs channels 7 - 54
    Epoch 4: abs channels 6 - 46
    Epoch 5: abs channels 6 - 48

    Passed: rate_ref_2d - The populated 2-dimensional light curve for the
                reference band, split up by energy channel (0-63 inclusive).
            obs_epoch - The RXTE observational epoch of the data

    Returns: rate_ref - The reference band, from 3 - 20 keV.

    """
    if obs_epoch == 5:
        rate_ref = np.sum(rate_ref_2d[:, 2:26], axis=1)  # EPOCH 5
        # channel 2 to 25 inclusive
    elif obs_epoch == 3:
        rate_ref = np.sum(rate_ref_2d[:, 3:29], axis=1)  # EPOCH 3
        # channel 3 to 28 inclusive
    else:
    	rate_ref = np.sum(rate_ref_2d[:, 3:32], axis=1)  # TEMP, NOT FOR REALS
    
    return rate_ref

## End of function 'stack_reference_band'


###############################################################################
def filter_freq(freq_space_array, dt, n_bins, signal_freq):
    """
            filter_freq

    Applying a filter to the averaged cross-spectrum per energy channel (in
    frequency space). Any cross spectrum amplitudes above or below signal_freq
    get zeroed out.

    Passed: freq_space_array - The average cross spectrum over all segments (and
                all files, if applicable).
            dt - Time step between bins, in seconds.
            n_bins - Number of bins in one segment.
            signal_freq - The frequency, in Hz, of the signal we want to filter
            around.

    Returns: filt_freq_space_array - The cross-spectrum filtered in 
    		 	Fourier-frequency-space.
    		 j_min - The first frequency bin in which we keep the signal.
    		 j_max - 1+ the last frequency bin in which we keep the signal.
    		 
    The signal goes from j_min to j_max+1, so freq[j_min:j_max] gives you only 
    the bins with the kept frequencies.

    """
    freq = fftpack.fftfreq(n_bins, d=dt)
    min_freq_mask = freq < signal_freq  # we want the last 'True' element
    max_freq_mask = freq > signal_freq  # we want the first 'True' element
    j_min = list(min_freq_mask).index(False)
    j_max = list(max_freq_mask).index(True)
    # 	print j_min, j_max
    # 	print freq[j_min]
    # 	print freq[j_max-1]
    zero_front = np.zeros((j_min, 64), dtype=np.complex128)
    zero_end = np.zeros((len(freq_space_array) - j_max, 64), dtype=np.complex128)
    # 	print np.shape(zero_front)
    # 	print np.shape(freq_space_array[j_min:j_max,:])
    # 	print np.shape(zero_end)
    filt_freq_space_array = np.concatenate((zero_front,
    freq_space_array[j_min:j_max, :], zero_end), axis=0)
    assert np.shape(freq_space_array) == np.shape(filt_freq_space_array)
    
    return filt_freq_space_array, j_min, j_max

## End of function 'filter_freq'


###############################################################################
def make_cs(rate_ci, rate_ref, n_bins, dt):
    """
            make_cs

    Generating the cross spectrum for one segment of the light curve.

    Passed: rate_ci - The count rate for this segment of the channels of 
    			interest. Index orientation: (n_bins, energy_channels)
            rate_ref - The count rate for this segment of the reference band.
            n_bins - Number of bins in one segment.
            dt - Time step between bins, in seconds.

    Returns: cs_segment - The cross spectrum for this segment.
             mean_rate_segment_ci - The mean photon count rate of the channels 
             	of interest.
             mean_rate_segment_ref - The mean photon count rate of the reference 
             	band.

    """
    ## Computing the mean count rate of the segment
    mean_rate_segment_ci = np.mean(rate_ci, axis=0)
#     print "CI make:", np.mean(rate_ci) * 64
#     print mean_rate_segment_ci
#     print "CI make:", np.sum(mean_rate_segment_ci)
    mean_rate_segment_ref = np.mean(rate_ref)
#     print "Ref make:", mean_rate_segment_ref
    # 	print "Sum of mean rate segment 1:", np.sum(mean_rate_segment_ci)
    # 	print "Shape of mean rate segment 1:", np.shape(mean_rate_segment_ci)
    # 	print "Mean rate segment 2:", mean_rate_segment_ref
    ## Subtracting the mean off each value of 'rate'
    # 	print "Shape of rate 1", np.shape(rate_ci)
    rate_sub_mean_ci = np.subtract(rate_ci, mean_rate_segment_ci)
    rate_sub_mean_ref = np.subtract(rate_ref, mean_rate_segment_ref)
    # 	print "Shape of rate_sub_mean 2", np.shape(rate_sub_mean_ref)
    # 	print "Shape of rate sub mean 1", np.shape(rate_sub_mean_ci)


    ## Taking the FFT of the time-domain photon count rate
    ##  Using the SciPy FFT algorithm, as it's faster than NumPy for large lists
    fft_data_ci = fftpack.fft(rate_sub_mean_ci, axis=0)
    fft_data_ref = fftpack.fft(rate_sub_mean_ref)

    power_ci = np.absolute(fft_data_ci) ** 2
    power_ref = np.absolute(fft_data_ref) ** 2

    fft_ref = np.repeat(fft_data_ref[:, np.newaxis], 64, axis=1)
    assert np.shape(fft_ref) == np.shape(fft_data_ci)
	
    cs_segment = np.multiply(fft_data_ci, np.conj(fft_ref))
    # 	print "Shape of cs_segment", np.shape(cs_segment)
    
    return cs_segment, mean_rate_segment_ci, mean_rate_segment_ref, power_ci, \
        power_ref

## End of function 'make_cs'


###############################################################################
def each_segment(time_ci, time_ref, energy_ci, energy_ref, n_bins, dt, start_time, end_time, obs_epoch, sum_rate_whole_ci, sum_rate_whole_ref, sum_power_ci, sum_power_ref, cs_sum, num_segments, sum_rate_ci):
	"""
			each_segment
	
	Turns the event list into a populated histogram, stacks the reference band, 
	and makes the cross spectrum, per segment of light curve.
	
	Passed: time_ci - List of times of events in ci.
			time_ref - List of times of events in ref.
			energy_ci - List of energies of events in ci.
			energy_ref - List of energies of events in ref.
			n_bins - Number of bins in one segment of the light curve. Must be
                a power of 2.
			dt - Time step between bins, in seconds. Must be a power of 2.
			start_time - The start time of the segment.
			end_time - The end time of the segment.
			obs_epoch - The RXTE observational epoch.
			sum_rate_whole_ci - The sum of the count rates in ci.
			sum_rate_whole_ref - The sum of the count rates in ref.
			sum_power_ci - The sum of the power spectra in ci.
			sum_power_ref - The sum of the power spectra in ref.
			cs_sum - The sum of the cross spectra.
			num_segments - The number of segments computed.
			sum_rate_ci - The sum of the count rates in ci, not split per energy 		
				band.
	
	Returns: cs_sum - The sum of all the cross spectra.
			 sum_rate_whole_ci - The sum of the count rates in ci.
			 sum_rate_whole_ref - The sum of the count rates in ref.
			 num_segments - The number of segments computed.
			 sum_power_ci - The sum of the power spectra in ci.
			 sum_power_ref - The sum of the power spectra in ref.
			 sum_rate_ci - The sum of the count rates in ci, not split per 
			 	energy band.
	
	"""
	if len(time_ci) > 0 and len(time_ref) > 0:
	# Only take a cross spectrum if there's stuff in both lists

		assert len(time_ci) == len(energy_ci)
		assert len(time_ref) == len(energy_ref)
		mean_rate_segment_ci = np.zeros(64, dtype=np.float64)
		mean_rate_segment_ref = np.zeros(64, dtype=np.float64)
		cs_segment = np.zeros((n_bins, 64), dtype=np.complex128)

		rate_ci_2d, rate_ci_1d = tools.make_lightcurve( \
			np.asarray(time_ci), np.asarray(energy_ci), n_bins, dt, start_time)
		rate_ref_2d, rate_ref_1d = tools.make_lightcurve( \
			np.asarray(time_ref), np.asarray(energy_ref), n_bins, dt, \
			start_time)

		rate_ref = stack_reference_band(rate_ref_2d, obs_epoch)
		cs_segment, mean_rate_segment_ci, mean_rate_segment_ref, power_ci, \
			power_ref = make_cs(rate_ci_2d, rate_ref, n_bins, dt)
# 		print mean_rate_segment_ci
# 		print "CI each:", np.sum(mean_rate_segment_ci)
# 		print "Ref each:", mean_rate_segment_ref
# 		print "From rateci2d:", np.sum(np.mean(rate_ci_2d, axis=0))
		sum_rate_whole_ci += mean_rate_segment_ci
# 		print sum_rate_whole_ci
		sum_rate_whole_ref += mean_rate_segment_ref
		sum_power_ci += power_ci
		sum_power_ref += power_ref
		cs_sum += cs_segment  # This adds indices
		
		
		sum_rate_ci += np.mean(rate_ci_1d)
# 		print "CI real:", np.mean(rate_ci_1d)
		num_segments += 1
# 		print "1d", np.sum(rate_ci_1d)
# 		print "2d", np.sum(rate_ci_2d)

		
		
	## End of 'if there are events in the segment time'
	
# 	print np.sum(sum_rate_whole_ci)
# 	print sum_rate_ci
	
	## Clearing variables from memory
	time_ci = None
	energy_ci = None
	time_ref = None
	energy_ref = None
	mean_rate_segment_ci = None
	mean_rate_segment_ref = None
	cs_segment = None
	rate_ci_2d = None
	rate_ref_2d = None
	rate_ci_1d = None
	rate_ref_1d = None
	rate_ref = None
		
	return cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
		sum_power_ci, sum_power_ref, sum_rate_ci
		
## End of function 'each_segment'


###############################################################################
def fits_read_and_use_segments(in_file, n_bins, dt, test):


###############################################################################
def read_and_use_segments(in_file, n_bins, dt, test):
    """
            read_and_use_segments

    Reading in the eventlist to make the cross spectrum. Reads in a clock- 	
    corrected GTI'd event list, populates the light curves, computes cross 
    spectra per energy channel and keeps running average of cross spectra.

    Passed: in_file - The name of the event list with both reference and
                interest events.
            n_bins - Number of bins in one segment of the light curve. Must be
                a power of 2.
            dt - Time step between bins, in seconds. Must be a power of 2.
            test - True if computing one segment, False if computing all.

    Returns: cs_sum - The sum of the cross spectra over all segments of the
                light curves, one per energy channel.
             sum_rate_whole_ci - Sum of the mean count rates of all segments of
                light curve 1, per energy channel.
             sum_rate_whole_ref - Sum of the mean count rates of all segments
                of light curve 2, per energy channel.
             num_segments - Number of segments computed.

    """
    assert tools.power_of_two(n_bins)

    print "Input file: %s" % in_file

    if n_bins == 32768:
        print_iterator = int(10)
    elif n_bins < 32768:
        print_iterator = int(10)
    else:
        print_iterator = int(1)
    # 	print "Print iterator =", print_iterator

    fits_file = in_file[0:-4] + ".fits"
    # 	print fits_file
    obs_epoch = tools.obs_epoch_rxte(fits_file)
#     print "Observation epoch:", obs_epoch

    ## Initializations
    num_segments = 0
    sum_rate_whole_ci = np.zeros(64, dtype=np.float64)
    sum_rate_whole_ref = 0
    cs_sum = np.zeros((n_bins, 64), dtype=np.complex128)
    time_ci = []
    energy_ci = []
    time_ref = []
    energy_ref = []
    start_time = -99
    final_time = -99
    sum_power_ci = np.zeros((n_bins, 64), dtype=np.float64)
    sum_power_ref = np.zeros(n_bins, dtype=np.float64)
    sum_rate_ci = 0

    ## Reading only the first line of data to get the start time of the file
    with open(in_file, 'r') as fo:
        for line in fo:
            if line[0].strip() != "#":
                line = line.strip().split()
                start_time = np.float64(line[0])
                break
    end_time = start_time + (dt * n_bins)
#     print "Start time is %.21f" % start_time
# 	print "End time of first seg is %.21f" % end_time

    print "Segments computed:"
    with open(in_file, 'r') as f:
        for line, next_line in tools.pairwise(f):
            if line[0].strip() != "#":  # If the line is not a comment
                line = line.strip().split()
                next_line = next_line.strip().split()
                current_time = np.float64(line[0])
#                 print "%.21f" % current_time
                current_chan = np.int8(line[1])
                current_pcu = np.int8(line[2])
                next_time = np.float64(next_line[0])

                if current_pcu == 2:  ## If pcu = 2
                    time_ci.append(current_time)
                    energy_ci.append(current_chan)

                elif current_pcu == 0:  ## pcu = 0
                    time_ref.append(current_time)
                    energy_ref.append(current_chan)

                else:
                    print "\n\tERROR: PCU is not 0 or 2. Exiting."
                    sys.exit()

                if (next_time >= end_time):  # Triggered at end of a segment
				# Because we want start_time <= times < end_time
					cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, \
						num_segments, sum_power_ci, sum_power_ref, sum_rate_ci = \
						each_segment(time_ci, time_ref, energy_ci, energy_ref, \
						n_bins, dt, start_time, end_time, obs_epoch, \
						sum_rate_whole_ci, sum_rate_whole_ref, sum_power_ci, \
						sum_power_ref, cs_sum, num_segments, sum_rate_ci)
					
					if num_segments % print_iterator == 0:
						print "\t", num_segments
					if test is True and num_segments == 1:  # For testing
						break
					## Incrementing counters and loop control variables.
					start_time = end_time
					end_time += (n_bins * dt)
					time_ci = []
					energy_ci = []
					time_ref = []
					energy_ref = []
						
#                 elif next_time == final_time:
# 					print "Last segment!"
# # 					print "%.21f" % current_time
# 					next_chan = np.int8(next_line[1])
# 					next_pcu = np.int8(next_line[2])
#                 	
# 					if next_pcu == 2:  ## If pcu = 2
# 						time_ci.append(next_time)
# 						energy_ci.append(next_chan)
# 					elif next_pcu == 0:  ## pcu = 0
# 						time_ref.append(next_time)
# 						energy_ref.append(next_chan)
# 					else:
# 						print "\n\tERROR: PCU is not 0 or 2. Exiting."
# 						sys.exit()
# 						
# 					cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, \
# 						num_segments, sum_power_ci, sum_power_ref = \
# 						each_segment(time_ci, time_ref, energy_ci, energy_ref, \
# 						n_bins, dt, start_time, end_time, obs_epoch, \
# 						sum_rate_whole_ci, sum_rate_whole_ref, sum_power_ci, \
# 						sum_power_ref, cs_sum, num_segments)

				## End of 'if it`s at the end of a segment'
			## End of 'if the line is not a comment'
		## End of for-loop
	## End of with-block
		
    return cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
        sum_power_ci, sum_power_ref, sum_rate_ci
        
## End of function 'read_and_use_segments'


###############################################################################
def cs_to_ccf_w_err(cs_avg, dt, n_bins, num_seconds, total_segments, \
	mean_rate_total_ci, mean_rate_total_ref, mean_power_ci, mean_power_ref):
	"""
			cs_to_ccf_w_err
	
	Takes the iFFT of the filtered cross spectrum to get the cross-correlation 
	function, and computes the error on the cross-correlation function.
	
	Passed: cs_avg - Averaged cross spectrum of the data.
			dt - Time step between bins, in seconds.
			n_bins - Number of bins per Fourier segment. Must be a power of 2.
			num_seconds - Number of seconds per Fourier segment.
			total_segments - Total number of Fourier segments computed.
			mean_rate_total_ci - 
			mean_rate_total_ref - 
			mean_power_ci - 
			mean_power_ref - 
	
	Returns: ccf_filtered - The cross-correlation function of the frequency-
				filtered data.
			 ccf_error - The error on the frequency-filtered cross-correlation
			 	function amplitudes.
	
	"""
	filtered_cs_avg, j_min, j_max = filter_freq(cs_avg, dt, n_bins, 401.0)
	assert np.shape(filtered_cs_avg) == np.shape(cs_avg)
    
    ## Absolute rms norms of poisson noise
	noise_ci = 2.0 * mean_rate_total_ci
	noise_ref = 2.0 * mean_rate_total_ref
	noise_ref_array = np.repeat(noise_ref, 64)

	df = 1.0 / float(num_seconds)  # in Hz
#     print "df =", df
	
    ## Extracting only the signal frequencies of the mean powers
	signal_ci_pow = np.float64(mean_power_ci[j_min:j_max, :])
	signal_ref_pow = np.float64(mean_power_ref[j_min:j_max])
#     print j_min, j_max
	
    ## Putting powers into absolute rms2 normalization
	signal_ci_pow = signal_ci_pow * (2.0 * dt / float(n_bins)) - noise_ci
#     print "signal ci pow:", signal_ci_pow[:, 2:5]
	signal_ref_pow = signal_ref_pow * (2.0 * dt / float(n_bins)) - noise_ref
#     print "signal ref pow:", signal_ref_pow[2:5]
	
	
	## Getting rms of reference band, to normalize the ccf
	signal_variance = np.sum(signal_ref_pow * df)
	rms_ref = np.sqrt(signal_variance)  # should be a few percent in fractional rms units -- here it's in absolute rms units!!
	print "RMS of reference band:", rms_ref
    
    ## Putting signal_ref_pow in same shape as signal_ci_pow
	signal_ref_pow_stacked = signal_ref_pow
	for i in xrange(63):
		signal_ref_pow_stacked = np.column_stack(
			(signal_ref_pow, signal_ref_pow_stacked))
	assert np.shape(signal_ref_pow_stacked) == np.shape(signal_ci_pow)

	temp = (noise_ci * signal_ref_pow_stacked) + \
		(noise_ref * signal_ci_pow) + \
		(noise_ci * noise_ref)
# 	print "Shape of temp:", np.shape(temp)
	cs_noise_amp = np.sqrt(np.sum(temp, axis=0) / float(total_segments))
#     print "cs noise amp:", cs_noise_amp[2:5]

	temp1 = np.absolute(cs_avg[j_min:j_max, :]) * (2.0 * dt / float(n_bins))
	cs_signal_amp = np.sum(temp1, axis=0)
#     other_sig = np.sqrt(np.square(signal_ci_pow * signal_ref_pow_stacked) / \
#     	float(total_segments))
	
#     print "Shape of cs signal amp:", np.shape(cs_signal_amp)
#     print "Shape of other signal:", np.shape(other_sig)
#     print "CS signal amp:", cs_signal_amp
#     print "other signal amp:", other_sig[:,2:5]
#     print "shape of cs signal amp:", np.shape(cs_signal_amp)
#     print "CS noise amp:", cs_noise_amp

	## Assuming that cs_noise_amp and cs_signal_amp are float arrays, size 64
	error_ratio = np.zeros(64, dtype=np.float64)
	## A divide-by-zero RuntimeWarning in the following two lines means that there's no signal in that energy band, which is ok!
	error_ratio[:10] = cs_noise_amp[:10] / cs_signal_amp[:10]
	error_ratio[11:] = cs_noise_amp[11:] / cs_signal_amp[11:]
	
	error_ratio[np.where(np.isnan(error_ratio))] = 0.0
	
#     print "error ratio, noise on top:", error_ratio
#     print "Filtered cs, un-norm:", filtered_cs_avg[j_min:j_max,:]
#     print "Shape filt cs avg:", np.shape(filtered_cs_avg)
    
    ## Taking the IFFT of the cross spectrum to get the CCF
	ccf = fftpack.ifft(cs_avg, axis=0)
	ccf_filtered = fftpack.ifft(filtered_cs_avg, axis=0)
	assert np.shape(ccf) == np.shape(ccf_filtered)
    
    ## Dividing ccf by rms of signal in reference band
	ccf *= (2.0 / float(n_bins) / rms_ref)
	ccf_filtered *= (2.0 / float(n_bins) / rms_ref)
#     print "Unfilt norm CCF, 2-4:", ccf[0,2:5]
#     print "Filt norm ccf, 2-4:", ccf_filtered[0,2:5]
    
    ## Computing the error on the ccf
	ccf_rms_ci = np.sqrt(np.var(ccf_filtered, axis=0, ddof=1))
#     print "Shape of rms ci:", np.shape(ccf_rms_ci)
#     print "CCF rms ci:", ccf_rms_ci
#     print "Shape of error ratio:", np.shape(error_ratio)
	ccf_error = ccf_rms_ci * error_ratio
    
#     ccf_error *= (2.0 / float(n_bins) / rms_ref)

#     print "CCF:", ccf_filtered[0, 2:5]
#     print "CCF error:", ccf_error[2:5]
#     print "Shape of ccf error:", np.shape(ccf_error)
	
	return ccf_filtered, ccf_error
    
## End of function 'cs_to_ccf_w_err'


###############################################################################
def main(in_file, out_file, num_seconds, dt_mult, test):
    """
            main

    Reads in one event list, splits into two light curves, makes segments and
    populates them to give them length n_bins, computes the cross spectrum of
    each segment per energy channel and then averaged cross spectrum of all the
    segments per energy channel, and then computes the cross-correlation
    function (ccf) per energy channel.

    Passed: in_file - Name of (ASCII/txt/dat) input file event list containing
                both the reference band and the channels of interest. Assumes
                ref band = PCU 0, interest = PCU 2.
            out_file - Name of (ASCII/txt/dat) output file for ccf.
            num_seconds - Number of seconds in each Fourier segment. Must be a
                power of 2.
            dt_mult - Multiple of 1/8192 seconds for timestep between bins.
            test - True if computing one segment, False if computing all.

    Returns: nothing

    """
    t_res = 1.0 / 8192.0
    dt = dt_mult * t_res
    n_bins = num_seconds * int(1.0 / dt)

    ## Idiot checks, to ensure that our assumptions hold
    # 	assert num_seconds > 0 # num_seconds must be a positive integer
    assert n_bins > 0  # number of bins must be a positive integer
    assert dt > 0
    assert tools.power_of_two(n_bins)  # n_bins must be a power of 2 for the FFT

    # 	print "dt =", dt, "seconds"
    # 	print "bins =", n_bins

    mean_rate_whole_ci = np.zeros(64)
    mean_rate_whole_ref = 0
    ccf_filtered = np.zeros((n_bins, 64))
    cs_avg = np.zeros((n_bins, 64), dtype=np.complex128)

    cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, sum_power_ci, sum_power_ref, sum_rate_ci = read_and_use_segments(in_file, n_bins, dt, test)

    mean_rate_whole_ci = sum_rate_whole_ci / float(num_segments)
    mean_rate_whole_ref = sum_rate_whole_ref / float(num_segments)
    mean_power_ci = sum_power_ci / float(num_segments)
    mean_power_ref = sum_power_ref / float(num_segments)
    cs_avg = cs_sum / float(num_segments)

    ccf_filtered, ccf_error = cs_to_ccf_w_err(cs_avg, dt, n_bins, num_seconds, num_segments, mean_rate_whole_ci, mean_rate_whole_ref, mean_power_ci, mean_power_ref)
	
#     print "CCF:", ccf_filtered[0, 2:5]
#     print "Shape of ccf error:", np.shape(ccf_error)
#     print "CCF error:", ccf_error[2:5]
	
    print "Number of segments:", num_segments
    # 	print "Mean rate for ci:", mean_rate_whole_ci
    # 	print "Mean mean rate for ci:", np.mean(mean_rate_whole_ci)
    print "Sum of mean rate for ci:", np.sum(mean_rate_whole_ci)
    print "Mean rate for ref:", np.mean(mean_rate_whole_ref)
	
    t = np.arange(0, n_bins)
    # time = t * dt  # Converting to seconds

    ## Calling output function
    fits_output(out_file, in_file, dt, n_bins, num_seconds, num_segments,
        mean_rate_whole_ci, mean_rate_whole_ref, t, ccf_filtered, ccf_error)
	
## End of function 'main'


###############################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage='ccf.py infile outfile [-n \
    	NUM_SECONDS] [-m DT_MULT] [-t {0,1}]', description='Computes the cross-\
        correlation function of a channel of interest with a reference band.', \
        epilog='For optional arguments, default values are given in brackets \
        at end of description.')
    parser.add_argument('infile', help='Name of (ASCII/txt/dat) event list \
    	containing both the reference band and the channels of interest. \
    	Assumes ref band = PCU 0, interest = PCU 2.')
    parser.add_argument('outfile', help='The full path of the (ASCII/txt) file \
    	to write the frequency and power to.')
    parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two, \
    	default=1, dest='num_seconds', help='Number of seconds in each Fourier \
    	segment. Must be a power of 2. [1]')
    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two, \
    	default=1, dest='dt_mult', help='Multiple of 1/8192 seconds for \
    	timestep between bins. [1]')
    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
        dest='test', help='1 if only computing one segment for testing, 0 if \
        computing all segments. [0]')
    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True

    main(args.infile, args.outfile, args.num_seconds, args.dt_mult, test)


## End of program 'ccf.py'
