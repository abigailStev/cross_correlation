import argparse
import numpy as np
from scipy import fftpack
from astropy.io import fits  
from datetime import datetime

import populate_lightcurve as lc
import tools
import powerspec

"""
		ccf.py
		
Computes the cross-correlation function of two light curves to do phase-resolved 
spectroscopy. Able to read light curves from data with the program 'populate_lightcurve'.

Arguments:
in_file - Name of (ASCII/txt/dat) input file event list containing both the reference
	band and the channels of interest. Assumes ref band = PCU 0, interest = PCU 2.
out_file - Name of (ASCII/txt/dat) output file which the table of cross-correlation 
	function data will be written to.
num_seconds - Duration of each segment of the light curve, in seconds. Must be a power of 
	2.
dt_mult - Multiple of 1/8192 seconds for timestep between bins.
short_run - 1 if only computing one segment for testing, 0 if computing all segments.

Written in Python 2.7 by A.L. Stevens, A.L.Stevens@uva.nl, 2014

All scientific modules imported above, as well as python 2.7, can be downloaded in the 
Anaconda package, https://store.continuum.io/cshop/anaconda/

'tools.py' and 'populate_lightcurve.py' are available on GitHub at 
https://github.com/abigailStev/power_spectra

"""


#########################################################################################
def output(out_file, in_file, dt, n_bins, num_seconds, num_segments, \
	mean_rate_whole_ci, mean_rate_whole_ref, t, time, ccf, ccf_filtered, ccf_error):
	"""
			output
	
	Writes the cross-correlation function to an output file.
	
	Passed: out_file - Name of output file.
			in_file - Name of event list with both reference and interest bands.
			dt - Size of each time bin, in seconds.
			n_bins - Number of (time) bins per segment.
			num_seconds - Duration of each segment of the light curve, in seconds.
			num_segments - Number of segments the light curve was split up into.
			mean_rate_whole_ci - Mean count rate of light curve 1, averaged over all 
				segments.
			mean_rate_whole_ref - Mean count rate of light curve 2, averaged over all 
				segments.
			t - Integer time bins to plot against the ccf.
			time - Time (starting at arbitrary 0) to plot against the ccf. The time is 
				the length of one Fourier segment.
			ccf - CCF amplitudes, not filtered in frequency space.
			ccf_filtered - CCF amplitudes, filtered in frequency space.
			ccf_error - Error on the filtered CCF.
			
	Returns: nothing
			
	"""
	print "Output sent to %s" % out_file

	with open(out_file, 'w') as out:
		out.write("#\t\tCross-correlation function data")
		out.write("\n# Event list: %s" % in_file)
		out.write("\n# Time bin size = %.21f seconds" % dt)
		out.write("\n# Number of bins per segment = %d" % n_bins)
		out.write("\n# Number of seconds per segment = %d" % num_seconds)
		out.write("\n# Number of segments per light curve = %d" % num_segments)
		out.write("\n# Duration of light curve used = %d seconds" \
			% (num_segments * num_seconds))
		out.write("\n# Mean count rate = %.8f, over light curve 1" \
			% np.mean(mean_rate_whole_ci))
		out.write("\n# Mean count rate = %.8f, over light curve 2" \
			% np.mean(mean_rate_whole_ref))
		out.write("\n# ")
		out.write("\n# Column 1: Integer time bins")
		out.write("\n# Columns 2-65: Filtered ccf per energy channel, real part")
		out.write("\n# Columns 66-129: ccf per energy channel, real part")
		out.write("\n# ")
		for j in xrange(0, n_bins):
			out.write("\n%d" % t[j])
			for i in xrange(0,64):
				out.write("\t%.5f" % (ccf_filtered[j][i].real))
			for i in xrange(0,64):
				out.write("\t%.5f" % (ccf_error[j][i].real))
			## End of for-loops
		## End of with-block
	## End of function 'output'


#########################################################################################
def filter_freq_noise(cross_avg, dt, n_bins):
	"""
			filter_freq_noise
	
	Applying a filter to the averaged cross-spectrum per energy channel (in frequency 
	space). The current values are specific to SAX J1808, so any amplitudes above or 
	below 401 get zeroed out for now.
	
	Passed: cross_avg - The average cross spectrum over all segments (and all files, 
				if applicable).
			dt - Time step between bins, in seconds. Must be a power of 2.
	
	Returns: smoothed_cross_avg_i
	
	"""
	
# 	print "Length of cross spectrum before filter =", len(cross_avg_i)
	freq = fftpack.fftfreq(n_bins, d=dt)

	min_freq_mask = freq <= 400.75 # we want the last 'True' element
# 	print min_freq_mask
	max_freq_mask = freq >= 401.25 # we want the first 'True' element
	j_min = list(min_freq_mask).index(False)-1 # If this is the first False element, 
											   # then the last True one is -1.
	j_max = list(max_freq_mask).index(True)
# 	print j_min, j_max
# 	print "Shape of cross avg", np.shape(cross_avg)
	smooth_front = np.zeros(j_min)
# 	print "Length of cross avg", len(cross_avg)
	smooth_end = np.zeros((len(cross_avg) - j_max))
# 	print len(smooth_front)+len(cross_avg[j_min:j_max])+len(smooth_end)
# 	smoothed_cross_avg = np.zeros((n_bins,64), dtype=np.complex128)
	
	smooth_front_alt = np.reshape(np.tile(smooth_front, 64), (len(smooth_front), 64))
	smooth_end_alt = np.reshape(np.tile(smooth_end, 64), (len(smooth_end), 64))
	smoothed_cross_avg = np.concatenate((smooth_front_alt, cross_avg[j_min:j_max,:],
                                         smooth_end_alt))
# 	print "Length of cross spectrum after filter =", len(smoothed_cross_avg_i)
# 	print "Shape of smoothed cross avg", np.shape(smoothed_cross_avg)
	return smoothed_cross_avg
	## End of function 'filter_freq_noise'
	

#########################################################################################
def stack_reference_band(rate_ref_2d, obs_epoch):
	"""
			stack_reference_band
			
	Stacks the photons in the reference band from 3-20 keV to make one 'band'.
	Epoch 1: abs channels 10 - 74
	Epoch 2: abs channels 9 - 62
	Epoch 3: abs channels 7 - 54
	Epoch 4: abs channels 6 - 46
	Epoch 5: abs channels 6 - 48
	
	Passed: rate_ref_2d - The populated 2-dimensional light curve for the 
				second (reference) curve, split up by energy channel (0-63 inclusive).
			obs_epoch - The RXTE observational epoch of the data
	
	Returns: rate_ref - The reference band, from 3 - 20 keV.
	
	"""
	
	if obs_epoch == 5:
		rate_ref = np.sum(rate_ref_2d[:, 2:26], axis = 1) # EPOCH 5 -- channel 2 to 25, incl.
	else:
		rate_ref = np.sum(rate_ref_2d[:, 3:29], axis = 1) # EPOCH 3 -- channel 3 to 28, incl.
	
# 	print "SHAPE OF RATE 2", np.shape(rate_ref)
	return rate_ref
	## End of function 'stack_reference_band'
	

#########################################################################################
def each_segment(rate_ci, rate_ref, n_bins, dt):
	"""
			each_segment
	
	Generating the cross spectrum for each segment of the light curve.
	
	Passed: rate_ci - The count rate for this segment of light curve 1.
			rate_ref - The count rate for this segment of light curve 2.
	
	Returns: cross_segment - The cross spectrum for this segment.
			 mean_rate_segment_ci - The mean photon count rate of light curve 1.
			 mean_rate_segment_ref - The mean photon count rate of light curve 2.
	
	"""
# 	print "Each segment"
	## Computing the mean count rate of the segment
	mean_rate_segment_ci = np.mean(rate_ci, axis = 0)
	mean_rate_segment_ref = np.mean(rate_ref)
# 	print "Sum of mean rate segment 1:", np.sum(mean_rate_segment_ci)
# 	print "Shape of mean rate segment 1:", np.shape(mean_rate_segment_ci)
# 	print "Mean rate segment 2:", mean_rate_segment_ref
	## Subtracting the mean off each value of 'rate' to eliminate the spike at 0 Hz 
# 	print "Shape of rate 1", np.shape(rate_ci)
	rate_sub_mean_ci = np.subtract(rate_ci, mean_rate_segment_ci)
	rate_sub_mean_ref = np.subtract(rate_ref, mean_rate_segment_ref)
# 	print "Shape of rate_sub_mean 2", np.shape(rate_sub_mean_ref)
# 	print "Shape of rate sub mean 1", np.shape(rate_sub_mean_ci)

	
	## Taking the FFT of the time-domain photon count rate
	##  Using the SciPy FFT algorithm, as it is faster than NumPy for large lists
	fft_data_ci = fftpack.fft(rate_sub_mean_ci, axis = 0)
	fft_data_ref = fftpack.fft(rate_sub_mean_ref)
	
	power_ci = np.absolute(fft_data_ci) ** 2
	power_ref = np.absolute(fft_data_ref) ** 2
	
	fft_ref = np.repeat(fft_data_ref[:, np.newaxis], 64, axis = 1)
	assert np.shape(fft_ref) == np.shape(fft_data_ci)	
	
	cross_segment = np.multiply(fft_data_ci, np.conj(fft_ref))

# 	print "Shape of cross_segment", np.shape(cross_segment)
	return cross_segment, mean_rate_segment_ci, mean_rate_segment_ref, power_ci, power_ref
	## End of function 'each_segment'


#########################################################################################
def make_crossspec(in_file, n_bins, dt, short_run):
	"""
			make_crossspec
			
	Making the cross spectrum. Reads in a clock-corrected GTI'd event list, populates the 
	light curves, computes cross spectra per energy channel and keeps running average of 
	cross spectra.
	
	Passed: in_file - The name of the event list with both reference and interest events.
			n_bins - Number of bins in one segment of the light curve. Must be a power of 
				2.
			dt - Time step between bins, in seconds. Must be a power of 2.
			short_run - True if computing one segment, False if computing all.
	
	Returns: cross_sum - The sum of the cross spectra over all segments of the light
				curves, one per energy channel.
			 sum_rate_whole_ci - Sum of the mean count rates of all segments of light 
			 	curve 1, per energy channel.
			 sum_rate_whole_ref - Sum of the mean count rates of all segments of light 
			 	curve 2, per energy channel.
			 num_segments - Number of segments computed.
	
	"""
	assert tools.power_of_two(n_bins)
	assert tools.power_of_two(int(1.0 / dt))
	
	print "Input file: %s" % in_file

	if n_bins == 32768:
		print_iterator = int(10)
	elif n_bins < 32768:
		print_iterator = int(10)
	else:
		print_iterator = int(1)
# 	print "Print iterator =", print_iterator

	fits_file = in_file[0:-4]+".fits"
# 	print fits_file
	obs_epoch = tools.obs_epoch_rxte(fits_file)
	print "Observation epoch:", obs_epoch
	
	## Initializations
	num_segments = 0
	sum_rate_whole_ci = np.zeros(64)
	sum_rate_whole_ref = 0
	cross_sum = np.zeros((n_bins, 64), dtype = np.complex128)
	time_ci = []
	energy_ci = []
	time_ref = []
	energy_ref = []
	start_time = -99
	sum_power_ci = 0
	sum_power_ref = 0

	## Reading only the first line of data to get the start time of the file
	with open(in_file, 'r') as fo:
		for line in fo:
			if line[0].strip() != "#":
				line = line.strip().split()
				start_time = np.float64(line[0])
				break
	end_time = start_time + (dt * n_bins)
# 	print "Start time is %.21f" % start_time
# 	print "End time of first seg is %.21f" % end_time
	assert end_time > start_time	
	
	print "Segments computed:"
	with open(in_file, 'r') as f:
		for line, next_line in tools.pairwise(f):		
			if line[0].strip() != "#":  # If the line is not a comment
				line = line.strip().split()
				next_line = next_line.strip().split()
				current_time = np.float64(line[0])
				current_chan = np.int8(line[1])
				current_pcu = np.int8(line[2])
				next_time = np.float64(next_line[0])
# 				print "%.21f %d %d" % (current_time, current_chan, current_pcu)
	
				if current_pcu == 2:  ## If pcu = 2
					time_ci.append(current_time)
					energy_ci.append(current_chan)
					
				elif current_pcu == 0:  ## pcu = 0
					time_ref.append(current_time)
					energy_ref.append(current_chan)
				
				else:
					print "\n\tERROR: PCU is not 0 or 2. Exiting."
					exit()
			
				if (next_time > end_time):  # Triggered at end of a segment
					# Only take a cross spectrum if there's stuff in the list
					if len(time_ci) > 0 and len(time_ref) > 0:  
					
						assert len(time_ci) == len(energy_ci)
						assert len(time_ref) == len(energy_ref)
						mean_rate_segment_ci = np.zeros(64)
						mean_rate_segment_ref = np.zeros(64)
						cross_segment = np.zeros((n_bins, 64), dtype = np.complex128)

						rate_ci_2d, rate_ci_1d = lc.make_lightcurve(np.asarray(time_ci), \
							np.asarray(energy_ci), n_bins, dt, start_time)
						rate_ref_2d, rate_ref_1d = lc.make_lightcurve(np.asarray(time_ref), \
							np.asarray(energy_ref), n_bins, dt, start_time)
# 						print "Light curves populated."
						
						rate_ref = stack_reference_band(rate_ref_2d, obs_epoch)
						cross_segment, mean_rate_segment_ci, mean_rate_segment_ref, \
							power_ci, power_ref = each_segment(rate_ci_2d, rate_ref, \
							n_bins, dt)

						sum_rate_whole_ci += mean_rate_segment_ci
						sum_rate_whole_ref += mean_rate_segment_ref
						sum_power_ci += power_ci
						sum_power_ref += power_ref
						cross_sum += cross_segment  # This adds indices instead of concatenating arrays

						num_segments += 1

						if num_segments % print_iterator == 0:
							print "\t", num_segments
							
					## Clearing variables from memory
					time_ci = None
					energy_ci = None
					time_ref = None
					energy_ref = None
					mean_rate_segment_ci = None
					mean_rate_segment_ref = None
					cross_segment = None
					rate_ci_2d = None
					rate_ref_2d = None
					rate_ci_1d = None
					rate_ref_1d = None
					rate_ref = None
					## Incrementing counters and loop control variables.
					start_time += (n_bins * dt)
					end_time += (n_bins * dt)
					time_ci = []
					energy_ci = []
					time_ref = []
					energy_ref = []
					
					if short_run == True and num_segments == 20:  ## For testing purposes only
						break
						
				## End of 'if the line is not a comment'
			## End of for-loop 
		## End of with-block	
	return cross_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, sum_power_ci, \
		sum_power_ref
	## End of function 'make_crossspec'


#########################################################################################
def main(in_file, out_file, num_seconds, dt_mult, short_run): 
	"""	
			main
		
	Reads in one event list, splits into two light curves, makes segments and populates 
	them to give them length n_bins, computes the cross spectrum of each segment per 
	energy channel and then averaged cross spectrum of all the segments per energy 
	channel, and then computes the cross-correlation function (ccf) per energy channel.

	Passed: in_file - Name of (ASCII/txt/dat) input file event list containing both 
				the reference band and the channels of interest. Assumes ref band = PCU 0,
				interest = PCU 2.
			out_file - Name of (ASCII/txt/dat) output file for ccf.
			num_seconds - Number of seconds each segment should last. Must be a power of 
				2. 
			dt_mult - Multiple of 1/8192 seconds for timestep between bins. 
			short_run - True if computing one segment, False if computing all.
			
	Returns: nothing
	
	"""
	pass
	
	t_res = 1.0 / 8192.0
	dt = dt_mult * t_res
	n_bins = num_seconds * int(1.0 / dt)
	
	## Idiot checks, to ensure that our assumptions hold
# 	assert num_seconds > 0 # num_seconds must be a positive integer
	assert n_bins > 0 # number of bins must be a positive integer
	assert dt > 0
	assert tools.power_of_two(n_bins) # n_bins must be a power of 2 for the FFT
	
# 	print "dt =", dt, "seconds"
# 	print "bins =", n_bins
	
	mean_rate_whole_ci = np.zeros(64)
	mean_rate_whole_ref = 0
	ccf = np.zeros((n_bins, 64))
	ccf_filtered = np.zeros((n_bins, 64))
	cross_avg = np.zeros((n_bins, 64), dtype = np.complex128)

	cross_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, sum_power_ci, \
		sum_power_ref  = make_crossspec(in_file, n_bins, dt, short_run)
	
	mean_rate_whole_ci = sum_rate_whole_ci / float(num_segments)
	mean_rate_whole_ref = sum_rate_whole_ref / float(num_segments)
	mean_power_ci = sum_power_ci / float(num_segments)
	mean_power_ref = sum_power_ref / float(num_segments)
	cross_avg = cross_sum / float(num_segments)
	
	## Applying absolute rms normalization, and subtracting the noise
	cross_avg = cross_avg * (2.0 * dt / float(n_bins)) - (2.0 * mean_rate_whole_ref)
	
	old_settings = np.seterr(divide='ignore')
	
	## Normalizing ref band power to noise-subtracted fractional rms2, 
	## integrating signal power (i.e. computing the variance), sqrt'ing that to get rms
	freq = fftpack.fftfreq(n_bins, d=dt)
	min_freq_mask = freq < 401 # we want the last 'True' element
	max_freq_mask = freq > 401 # we want the first 'True' element
	j_min = list(min_freq_mask).index(False)-1
	j_max = list(max_freq_mask).index(True)
	df = freq[1] - freq[0]  # in Hz
	signal_ref_rms2 = 2.0 * mean_power_ref[j_min:j_max] * dt / float(n_bins) \
						/ (mean_rate_whole_ref ** 2)
	signal_ref_rms2 -= (2.0 / mean_rate_whole_ref)
	signal_variance = np.sum(signal_ref_rms2 * df) 
	rms_ref = np.sqrt(signal_variance)  # should be a few percent in fractional rms units
	
	## Filtering the cross spectrum in frequency
	
	np.savetxt('cs_avg.dat', cross_avg.real)
	
	
	zero_front = np.zeros((j_min, 64))
	zero_end = np.zeros((len(cross_avg) - j_max, 64))
	filtered_cross_avg = np.concatenate((zero_front, cross_avg[j_min:j_max,:], zero_end), \
						axis=0)
	
	print n_bins
	print freq[j_min]
	print freq[j_max-1]
	print np.shape(zero_front)
	print j_min
	print np.shape(cross_avg[j_min:j_max,:])
	print j_max
	print np.shape(zero_end)
	old_filtered_cs_avg = filter_freq_noise(cross_avg, dt, n_bins)
	print filtered_cross_avg == old_filtered_cs_avg
	assert np.shape(filtered_cross_avg) == np.shape(cross_avg)
		
	signal_ci_pow = np.complex128(mean_power_ci[j_min:j_max, :])
	signal_ref_pow = np.complex128(mean_power_ref[j_min:j_max])
# 	pow_noise_ci = np.concatenate((mean_power_ci[:j_min, :], mean_power_ci[j_max:, :]), axis=0)
# 	pow_noise_ref = np.concatenate((mean_power_ref[:j_min], mean_power_ref[j_max:]))
		
	signal_ref_pow_stacked = signal_ref_pow
# 	pow_noise_ref_stacked = pow_noise_ref
	for i in xrange(63):
		signal_ref_pow_stacked = np.column_stack((signal_ref_pow, signal_ref_pow_stacked))
# 		pow_noise_ref_stacked = np.column_stack((pow_noise_ref, pow_noise_ref_stacked))
	
	assert np.shape(signal_ref_pow_stacked) == np.shape(signal_ci_pow)
# 	assert np.shape(pow_noise_ref_stacked) == np.shape(pow_noise_ci)
	
	## Putting powers into absolute rms normalization
	signal_ci_pow = signal_ci_pow * (2.0 * dt / float(n_bins)) 
	signal_ci_pow -= (2.0 * mean_rate_whole_ci)
	print "signal ci pow:", signal_ci_pow[:, 5:8]
	signal_ref_pow = signal_ref_pow_stacked * (2.0 * dt / float(n_bins)) 
	signal_ref_pow -= (2.0 * mean_rate_whole_ref)
	print "signal ref pow:", signal_ref_pow[:, 5:8]
# 	pow_noise_ci = np.sqrt(pow_noise_ci * (2.0 * dt / float(n_bins)) - (2.0 * mean_rate_whole_ci))
# 	pow_noise_ref = np.sqrt(pow_noise_ref * (2.0 * dt / float(n_bins)) - (2.0 * mean_rate_whole_ref))
	noise_ci = 2.0 * mean_rate_whole_ci
	noise_ref = 2.0 * mean_rate_whole_ref
	
	
	temp = np.sqrt(np.square(noise_ci * signal_ref_pow) + 
			np.square(noise_ref * signal_ci_pow) + 
			np.square(noise_ci * noise_ref))
# 	print "Shape of temp:", np.shape(temp)
	cross_noise_amp = np.sqrt(np.sum(temp) / float(n_bins))
	# Might be off by a factor of 2 here...
	
	print "RMS of reference band:", rms_ref
	cross_signal_amp = np.sum(cross_avg[j_min:j_max, :], axis=0)
	print "Sum of cross signal amp:", np.sum(cross_signal_amp)

# 	temp2 = np.sqrt(np.square(signal_ref_pow * signal_ci_pow)) * df
# 	print "shape of temp2:", np.shape(temp2)
# 	cross_signal_amp = np.sqrt(np.sum(temp2, axis=0) / float(num_segments))
	print "cross signal amp:", cross_signal_amp[5:8]
	print "shape of cross signal amp:", np.shape(cross_signal_amp)
	print "cross noise amp:", cross_noise_amp
	cross_error_ratio = cross_signal_amp / cross_noise_amp
	print "Shape cross error ratio:", np.shape(cross_error_ratio)
	print "cross error ratio:", cross_error_ratio[5:8]
	
	## Taking the IFFT of the cross spectrum to get the ccf
	ccf = fftpack.ifft(cross_avg, axis = 0)
	ccf_filtered = fftpack.ifft(filtered_cross_avg, axis = 0)
	assert np.shape(ccf) == np.shape(ccf_filtered)

	ccf_error = cross_error_ratio * ccf_filtered

	## Dividing ccf by integrated rms power of signal in reference band
	ccf /= rms_ref
	ccf_filtered /= rms_ref
	ccf_error /= rms_ref
	
	print "filt CCF:", ccf_filtered[0, 5:8]
	print "Shape of ccf error:", np.shape(ccf_error)
	print "CCF error:", ccf_error[0, 5:8]

	print "Number of segments:", num_segments
# 	print "Mean rate for curve 1:", mean_rate_whole_ci
# 	print "Mean mean rate for curve 1:", np.mean(mean_rate_whole_ci)
	print "Sum of mean rate for curve 1:", np.sum(mean_rate_whole_ci)
# 	print "Sum of means for curve 1 chan 0-12:", np.sum(mean_rate_whole_ci[0:13])
	print "Mean rate for curve 2:", np.mean(mean_rate_whole_ref)
	
	t = np.arange(0, n_bins)
	## Converting to seconds
	time = t * dt

	## Calling 'output' function
	output(out_file, in_file, dt, n_bins, num_seconds, num_segments, mean_rate_whole_ci, \
		mean_rate_whole_ref, t, time, ccf, ccf_filtered, ccf_error)
	
	## End of function 'main'
	

#########################################################################################
if __name__ == "__main__":
	"""
	Parsing cmd-line arguments and calling 'main'
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument('in_file', \
		help = "Name of (ASCII/txt/dat) input file event list containing both the \
		reference band and the channels of interest. Assumes ref band = PCU 0, interest \
		= PCU 2.")
	parser.add_argument('out_file', help = \
		"The full path of the (ASCII/txt) file to write the frequency and power to.")
	parser.add_argument('num_seconds', type = int, help = \
		"Duration of segments the light curve is broken up into, in seconds. \
		Must be a power of 2.")
	parser.add_argument('dt_mult', type = int, \
		help = "Multiple of 1/8192 seconds for timestep between bins.")
	parser.add_argument('short_run', type = int, help = \
		"1 if only computing one segment for testing, 0 if computing all segments.")
	args = parser.parse_args()
	
	assert args.short_run == 1 or args.short_run == 0
	if args.short_run == 1: 
		short_run = True
	else:
		short_run = False

	main(args.in_file, args.out_file, args.num_seconds, args.dt_mult, short_run)


## End of program 'ccf.py'
