import argparse
import numpy as np
import sys
from scipy import fftpack
from datetime import datetime
import os
from astropy.io import fits
import warnings
import tools  # https://github.com/abigailStev/whizzy_scripts
import timeit

__author__ = "Abigail Stevens"
__author_email__ = "A.L.Stevens at uva.nl"
__year__ = "2014-2015"
__description__ = "Computes the cross-correlation function of an interest band \
with a reference band from RXTE event-mode data."

"""
        ccf.py

Written in Python 2.7.

"""

################################################################################
def dat_out(out_file, in_file, bkgd_file, dt, n_bins, detchans, num_seconds, \
	num_segments, mean_rate_whole_ci, mean_rate_whole_ref, t, ccf, ccf_error, \
	filter):
    """
            dat_out

    Writes the cross-correlation function to a .dat output file.

    """
    if out_file[-4:].lower() == "fits":
		out_file = out_file[:-4]+"dat"
		
    print "\nOutput sent to: %s" % out_file

    with open(out_file, 'w') as out:
        out.write("#\t\tCross-correlation function")
        out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
        out.write("\n# Event list: %s" % in_file)
        out.write("\n# Background spectrum: %s" % bkgd_file)
        out.write("\n# Time bin size = %.21f seconds" % dt)
        out.write("\n# Number of bins per segment = %d" % n_bins)
        out.write("\n# DETCHANS = %d" % detchans)
        out.write("\n# Number of seconds per segment = %d" % num_seconds)
        out.write("\n# Number of segments per light curve = %d" % num_segments)
        out.write("\n# Exposure time = %d seconds" \
                  % (num_segments * num_seconds))
        out.write("\n# Mean count rate of ci = %s" \
        	% str(list(mean_rate_whole_ci)))
        out.write("\n# Mean count rate of ref band = %.8f" \
                  % np.mean(mean_rate_whole_ref))
        out.write("\n# Filter applied in frequency domain? %s" % str(filter))
        out.write("\n# ")
        out.write("\n# Column 1: Time bins")
        out.write("\n# Columns 2-65: CCF per energy channel, real \
            part [count rate]")
        out.write("\n# Columns 66-129: Error on CCF per energy \
            channel, real part [count rate]")
        out.write("\n# ")
        for j in xrange(0, n_bins):
            out.write("\n%d" % t[j])
            for i in xrange(0, 64):
                out.write("\t%.6e" % ccf[j][i].real)
            if filter:
				for i in xrange(0, 64):
					out.write("\t%.6e" % ccf_error[i].real)
            else:
				for i in xrange(0, 64):
					out.write("\t%.6e" % ccf_error[j][i].real)
        ## End of for-loops
    ## End of with-block
## End of function 'dat_out'


################################################################################
def fits_out(out_file, in_file, bkgd_file, dt, n_bins, detchans, num_seconds, \
	num_segments, mean_rate_whole_ci, mean_rate_whole_ref, t, ccf, ccf_error, \
	filter):
    """
            fits_out

    Writes the cross-correlation function to a .fits output file.
    
    """
    
    ## Getting data into a good output structure
    chan = np.arange(0, detchans)
    energy_channels = np.tile(chan, len(t))
    if filter:
    	ccf_error = np.tile(ccf_error, len(t))
    else:
    	ccf_error = ccf_error.real.flatten('C')
    time_bins = np.repeat(t, detchans)
    assert len(energy_channels) == len(time_bins)
    
    print "\nOutput sent to: %s" % out_file
	
    ## Making FITS header (extension 0)
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross-correlation function")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file)
    prihdr.set('BKGD', bkgd_file)
    prihdr.set('DT', dt, "seconds")
    prihdr.set('N_BINS', n_bins, "time bins per segment")
    prihdr.set('SEGMENTS', num_segments, "segments in the whole light curve")
    prihdr.set('EXPOSURE', num_segments * n_bins * dt, \
    	"seconds, of light curve")
    prihdr.set('DETCHANS', detchans, "Number of detector energy channels")
    prihdr.set('RATE_CI', str(mean_rate_whole_ci.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_whole_ref, "counts/second")
    prihdr.set('FILTER', str(filter))
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
    
    ## If the file already exists, remove it (still working on just updating it)
    assert out_file[-4:].lower() == "fits", \
    	'ERROR: Output file must have extension ".fits".'
    if os.path.isfile(out_file):
    	os.remove(out_file)
    	
    ## Writing to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file)	
## End of function 'fits_out'


################################################################################
def get_phase_err(cs_avg, power_ci, power_ref, n, M):
	"""
			get_phase_err
	
	Computes the error on the complex phase (in radians). Power should not be 
	noise-subtracted.
	
	"""
	with np.errstate(all='ignore'):
		a = power_ci * power_ref
		coherence = np.where(a != 0, np.abs(cs_avg)**2 / a, 0)
		phase_err = np.sqrt(np.where(coherence != 0, (1 - coherence) / \
			(2 * coherence * n * M), 0))

	return phase_err
## End of function 'get_phase_err'


################################################################################
def phase_to_tlags(phase, freq, detchans):
	"""
			phase_to_tlags
	
	Converts a complex phase (in radians) to a time lag (in seconds).
	
	"""
	f = np.repeat(freq[:, np.newaxis], detchans, axis=1)
	with np.errstate(all='ignore'):
		tlags =  np.where(f != 0, phase / (2.0 * np.pi * f), 0)
		
	return tlags
## End of function 'phase_to_tlags'


################################################################################
def make_lags(out_file, in_file, dt, n_bins, detchans, num_seconds, \
	num_segments, mean_rate_whole_ci, mean_rate_whole_ref, cs_avg, power_ci, \
	power_ref):
	"""
			make_lags
			
	Computes the phase lag and time lag from the average cross spectrum, and 
	writes the lag information to a .fits output file.
	
	"""
	
	freq = fftpack.fftfreq(n_bins, d=dt)
	max_index = np.argmax(freq)+1  # because in python, the scipy fft makes the 
								   # nyquist frequency negative, and we want it 
								   # to be positive! (it is actually both pos 
								   # and neg)
	freq = np.abs(freq[1:max_index + 1])  # because it slices at end-1, and we 
										  # want to include 'max_index'; abs is
										  # because the nyquist freq is both pos
										  # and neg, and we want it pos here.
										  # but we don't want freq=0 because 
										  # that gives errors.
	cs_avg = cs_avg[1:max_index + 1]
	power_ci = power_ci[1:max_index + 1]
	power_ref = power_ref[1:max_index + 1]
	phase = np.arctan2(cs_avg.imag, cs_avg.real)
	err_phase = get_phase_err(cs_avg, power_ci, \
		np.repeat(power_ref[:, np.newaxis], detchans, axis=1), 1, num_segments)
	tlag = phase_to_tlags(phase, freq, detchans)
	err_tlag = phase_to_tlags(err_phase, freq, detchans)
    
	chan = np.arange(0,detchans)
	energy_channels = np.tile(chan, len(freq))
	bins = np.repeat(freq, len(chan))    
    
	out_file = out_file.replace("cross_correlation/out_ccf", "lags/out_lags")
	print "Output sent to: %s" % out_file
	
	## Making FITS header (extension 0)
	prihdr = fits.Header()
	prihdr.set('TYPE', "Time lag data")
	prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
	prihdr.set('EVTLIST', in_file)
	prihdr.set('DT', dt, "seconds")
	prihdr.set('N_BINS', n_bins, "time bins per segment")
	prihdr.set('SEGMENTS', num_segments, "segments in the whole light curve")
	prihdr.set('EXPOSURE', num_segments * n_bins * dt, \
		"seconds, of light curve")
	prihdr.set('DETCHANS', detchans, "Number of detector energy channels")
	prihdr.set('RATE_CI', str(mean_rate_whole_ci.tolist()), "counts/second")
	prihdr.set('RATE_REF', mean_rate_whole_ref, "counts/second")
	prihdr.set('FILTER', str(filter))
	prihdu = fits.PrimaryHDU(header=prihdr)

    ## Making FITS table (extension 1)
	col1 = fits.Column(name='FREQUENCY', format='D', array=bins)
	col2 = fits.Column(name='PHASE', unit='radians', format='D', \
		array=phase.flatten('C'))
	col3 = fits.Column(name='PHASE_ERR', unit='radians', format='D', \
		array=err_phase.flatten('C'))
	col4 = fits.Column(name='TIME_LAG', unit='s', format='D', \
		array=tlag.flatten('C'))
	col5 = fits.Column(name='TIME_LAG_ERR', unit='s', format='D', \
		array=err_tlag.flatten('C'))
	col6 = fits.Column(name='CHANNEL', unit='', format='I', \
		array=energy_channels)
	cols = fits.ColDefs([col1, col2, col3, col4, col5, col6])
	tbhdu = fits.BinTableHDU.from_columns(cols)
    
	## If the file already exists, remove it
	assert out_file[-4:].lower() == "fits", \
		'ERROR: Output file must have extension ".fits".'
	if os.path.isfile(out_file):
		os.remove(out_file)
	
	## Writing to a FITS file
	thdulist = fits.HDUList([prihdu, tbhdu])
	thdulist.writeto(out_file)	
	
## End of function 'make_lags'


################################################################################
def find_pulse_freq(freq, power_ref):
	"""
			find_pulse_freq
	
	Determines the frequency of a coherent pulse above 100 Hz (to not confuse
	with broadband noise).
	
	"""
	
	## Only searching above 100 Hz for a coherent pulse signal
	hf = np.where(freq > 100)
	hf_power = power_ref[hf]
	hf_freq = freq[hf]
	
	## Assuming that the pulse frequency will have the most power
	pulse_freq = hf_freq[np.argmax(hf_power)]

	return pulse_freq
## End of function 'find_pulse_freq'


################################################################################
def filter_freq(freq_space_array, dt, n_bins, detchans, power_ref):
    """
            filter_freq

    Applying a filter to the averaged cross-spectrum per energy channel (in
    frequency space). Any cross spectrum amplitudes above or below pulse_freq
    get zeroed out.

    """
    ## Compute the Fourier frequencies
    freq = fftpack.fftfreq(n_bins, d=dt)
    
    ## Determine pulse frequency
    pulse_freq = find_pulse_freq(freq, power_ref)
    
    print "Pulse frequency:", pulse_freq
    print "Index of pulse freq:", np.where(freq == pulse_freq)
    
    ## Get the indices of the beginning and end of the signal
    min_freq_mask = freq < pulse_freq  # we want the last 'True' element
    max_freq_mask = freq > pulse_freq  # we want the first 'True' element
    j_min = list(min_freq_mask).index(False)
    j_max = list(max_freq_mask).index(True)
	
    print "j min =", j_min
    print "j max =", j_max
	## Make zeroed arrays to replace with
    zero_front = np.zeros((j_min, detchans), dtype=np.complex128)
    zero_end = np.zeros((len(freq_space_array) - j_max, detchans), \
    	dtype=np.complex128)
	
	## Concatenate the arrays together
    filt_freq_space_array = np.concatenate((zero_front, \
    	freq_space_array[j_min:j_max, :], zero_end), axis=0)
    
    ## Check that the original array is the same shape as the filtered one
    assert np.shape(freq_space_array) == np.shape(filt_freq_space_array), \
    	"ERROR: Frequency-filtered cross spectrum does not have the same size \
as the original cross spectrum. Something went wrong."
    
    return filt_freq_space_array, j_min, j_max

## End of function 'filter_freq'
	

################################################################################
def FILT_cs_to_ccf_w_err(cs_avg, dt, n_bins, detchans, num_seconds, \
	total_segments, countrate_ci, countrate_ref, power_ci, power_ref, noisy):
	"""
			FILT_cs_to_ccf_w_err
	
	Filters the cross-spectrum in frequency space, takes the iFFT of the 
	filtered cross spectrum to get the cross-correlation function, and computes
	the error on the cross-correlation function. Note that error is NOT
	independent between time bins due to the filtering! But is still independent
	between energy bins.
	
	"""
	## Filter the cross spectrum in frequency
	filtered_cs_avg, j_min, j_max = filter_freq(cs_avg, dt, n_bins, detchans, \
		power_ref)
    
    ## Absolute rms norms of poisson noise
	noise_ci = 2.0 * countrate_ci
	noise_ref = 2.0 * countrate_ref
	
	## If there's no noise in a (simulated) power spectrum, noise level = 0
	if not noisy:
		noise_ci = np.zeros(detchans)
		noise_ref = 0
	
	noise_ref_array = np.repeat(noise_ref, detchans)
	
	df = 1.0 / float(num_seconds)  # in Hz
	
    ## Extracting only the signal frequencies of the mean powers
	signal_ci_pow = np.float64(power_ci[j_min:j_max, :])
	signal_ref_pow = np.float64(power_ref[j_min:j_max])
	
    ## Putting powers into absolute rms2 normalization, subtracting noise
	signal_ci_pow = signal_ci_pow * (2.0 * dt / float(n_bins)) - noise_ci
	signal_ref_pow = signal_ref_pow * (2.0 * dt / float(n_bins)) - noise_ref
	
	## Getting rms of reference band, to normalize the ccf
	ref_variance = np.sum(signal_ref_pow * df)
	print "Reference band variance:", ref_variance
	rms_ref = np.sqrt(ref_variance)  
	print "Frac RMS of reference band:", rms_ref / countrate_ref  
	## in frac rms units here -- should be few percent
    
    ## Broadcasting signal_ref_pow into same shape as signal_ci_pow
	signal_ref_pow = np.resize(np.repeat(signal_ref_pow, detchans), \
		np.shape(signal_ci_pow))
	assert np.shape(signal_ref_pow) == np.shape(signal_ci_pow)

	temp = (noise_ci * signal_ref_pow) + \
		(noise_ref * signal_ci_pow) + \
		(noise_ci * noise_ref)
	cs_noise_amp = np.sqrt(np.sum(temp, axis=0) / float(total_segments))

	temp1 = np.absolute(cs_avg[j_min:j_max, :]) * (2.0 * dt / float(n_bins))
	cs_signal_amp = np.sum(temp1, axis=0)

	## Assuming that cs_noise_amp and cs_signal_amp are float arrays, size 64
	error_ratio = np.zeros(detchans, dtype=np.float64)
	with np.errstate(all='ignore'):
		error_ratio = np.where(cs_signal_amp != 0, cs_noise_amp / cs_signal_amp, 0)
	
    ## Taking the IFFT of the cross spectrum to get the CCF
	ccf_end = fftpack.ifft(filtered_cs_avg, axis=0)
    
    ## Dividing ccf by rms of signal in reference band
	ccf_end *= (2.0 / float(n_bins) / rms_ref)
    
    ## Computing the error on the ccf
	ccf_rms_ci = np.sqrt(np.var(ccf_end, axis=0, ddof=1))
	ccf_error = ccf_rms_ci * error_ratio
	
	return ccf_end, ccf_error
    
## End of function 'FILT_cs_to_ccf_w_err'


################################################################################
def UNFILT_cs_to_ccf_w_err(cs_avg, dt, n_bins, detchans, num_seconds, \
	num_segments, countrate_ci, countrate_ref, power_ci, power_ref, noisy):
	"""
			UNFILT_cs_to_ccf_w_err
	Takes the iFFT of the cross spectrum to get the cross-correlation function, 
	and computes the error on the cross-correlation function. This error is 
	independent both between energy bins and between time bins!
	
	"""
	
	freq = fftpack.fftfreq(n_bins, d=dt)
	
	######################################################
	## Take the IFFT of the cross spectrum to get the CCF
	######################################################
	
	ccf_end = fftpack.ifft(cs_avg, axis=0)
		
	## Absolute rms norms of poisson noise	
	## Was over-estimating the noise from the count rate; using mean power from 
	## higher frequencies as the noise level
# 	noise_ref = 2.0 * countrate_ref
	print "Estimated noise:", 2.0 * countrate_ref
	noise_ci = 2.0 * countrate_ci - 0.016*(2.0 * countrate_ci)
	noise_ref = 2.0 * countrate_ref - 0.016*(2.0 * countrate_ref)
	temp = power_ref * (2.0 * dt / float(n_bins))
	if np.max(freq) > 100:
		noise_ref = np.mean(temp[np.where(freq >= 100)])
	
# 	print "NOISE CI", noise_ci
	print "NOISE REF:", noise_ref
# 	print "ABS RMS POWER CI", power_ci * (2.0 * dt / float(n_bins))
# 	print "ABS RMS POWER REF:", power_ref * (2.0 * dt / float(n_bins))
# 	print "SHAPE NOISE CI", np.shape(noise_ci)
# 	print "SHAPE NOISE REF", np.shape(noise_ref)
	
	df = 1.0 / float(num_seconds)  # in Hz
	
	## Putting powers into absolute rms2 normalization, subtracting noise
	absrms_power_ci = power_ci * (2.0 * dt / float(n_bins)) - noise_ci
	absrms_power_ref = power_ref * (2.0 * dt / float(n_bins)) - noise_ref
# 	absrms_power_ci[np.where(absrms_power_ci < 0)] = 0.0
# 	absrms_power_ref[np.where(absrms_power_ref < 0)] = 0.0
	print "Mean of whole ps:", np.mean(absrms_power_ref)
	if np.max(freq) > 100:
		print "Mean of hf ps:", np.mean(absrms_power_ref[np.where(freq >= 100)])
	
	
	## Computing the autocorrelation functions from the power spectra
	acf_ci = fftpack.ifft(absrms_power_ci)
	acf_ref = fftpack.ifft(absrms_power_ref)
	
	## Broadcasting acf_ref into same shape as acf_ci
	acf_ref = np.resize(acf_ref, np.shape(acf_ci))
	assert np.shape(acf_ref) == np.shape(acf_ci), "ERROR: Array broadcasting \
failed."
	
	## Bartlett formula for computing variance on unfilfered CCF
	var_ccf = (acf_ci * acf_ref) / num_segments
	
# 	print "Datatype of var_ccf:", type(var_ccf[0][0])
# 	print "Variance of ccf computed with Bartlett formula:", var_ccf[0:5,0]
# 	print "CCF:", ccf_end[0:5,0]
	
	## Getting rms of reference band, to normalize the ccf
	rms_variance = np.sum(absrms_power_ref * df)
	print "RMS variance:", rms_variance
	rms_ref = np.sqrt(rms_variance)  
	print "Frac RMS of reference band:", rms_ref / countrate_ref  
	# in frac rms units here -- should be few percent
	
	## Dividing ccf by rms of signal in reference band
	ccf_end *= (2.0 / float(n_bins) / rms_ref)
	var_ccf *= (2.0 / float(n_bins))
	
# 	
# 	print "Normed ccf:", ccf_end[0:5,0]
# 	print "Normed var_ccf:", var_ccf[0:5,0]	
	
	var_ccf = ccf_end * 0.0
	
	return ccf_end, var_ccf
	
## End of function 'UNFILT_cs_to_ccf_w_err'	


################################################################################
def stack_reference_band(rate_ref_2d, obs_epoch):
    """
            stack_reference_band
	
	WARNING: Only tested with RXTE event-mode detector channels, 0-63 incl.
	
    Stacks the photons in the reference band from 3-20 keV to make one 'band'.
    Assumes that we are in epoch 5 or 3.
    Should I be using chan.txt to assist with this?
    
    Epoch 1: abs channels 10 - 74
    Epoch 2: abs channels 9 - 62
    Epoch 3: abs channels 7 - 54
    Epoch 4: abs channels 6 - 46
    Epoch 5: abs channels 6 - 48

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


################################################################################
def make_cs(rate_ci, rate_ref, n_bins, detchans):
    """
            make_cs

    Generating the cross spectrum for one segment of the light curve.

    """
    ## Computing the mean count rate of the segment
    mean_rate_segment_ci = np.mean(rate_ci, axis=0)
    mean_rate_segment_ref = np.mean(rate_ref)
    
    ## Subtracting the mean off each value of 'rate'
    rate_sub_mean_ci = np.subtract(rate_ci, mean_rate_segment_ci)
    rate_sub_mean_ref = np.subtract(rate_ref, mean_rate_segment_ref)

    ## Taking the FFT of the time-domain photon count rate
    ## SciPy is faster than NumPy or pyFFTW for my array sizes
    fft_data_ci = fftpack.fft(rate_sub_mean_ci, axis=0)
    fft_data_ref = fftpack.fft(rate_sub_mean_ref)
	
	## Computing the power from the fourier transform
    power_ci = np.absolute(fft_data_ci) ** 2
    power_ref = np.absolute(fft_data_ref) ** 2
	
    ## Broadcasting fft of ref into same shape as fft of ci
    fft_data_ref = np.resize(np.repeat(fft_data_ref, detchans), (n_bins, \
    	detchans))
	
	## Computing the cross spectrum from the fourier transform
    cs_segment = np.multiply(fft_data_ci, np.conj(fft_data_ref))
    
    return cs_segment, mean_rate_segment_ci, mean_rate_segment_ref, power_ci, \
    	power_ref

## End of function 'make_cs'


################################################################################
def each_segment(time_ci, time_ref, energy_ci, energy_ref, n_bins, detchans, \
	dt, start_time, end_time, obs_epoch, sum_rate_whole_ci, sum_rate_whole_ref,\
	sum_power_ci, sum_power_ref, cs_sum, sum_rate_ci):
	"""
			each_segment
	
	Turns the event list into a populated histogram, stacks the reference band, 
	and makes the cross spectrum, per segment of light curve.
	
	"""
	assert len(time_ci) == len(energy_ci)
	assert len(time_ref) == len(energy_ref)
	
	## Initializations
	mean_rate_segment_ci = np.zeros(64, dtype=np.float64)
	mean_rate_segment_ref = np.zeros(64, dtype=np.float64)
	cs_segment = np.zeros((n_bins, 64), dtype=np.complex128)
	
	##############################################################
	## Populate the light curves for interest and reference bands
	##############################################################
	
	rate_ci_2d = tools.make_2Dlightcurve(np.asarray(time_ci), \
		np.asarray(energy_ci), n_bins, detchans, dt, start_time)
	rate_ref_2d = tools.make_2Dlightcurve( np.asarray(time_ref), \
		np.asarray(energy_ref), n_bins, detchans, dt, start_time)
	
	## Stack the reference band
	rate_ref = stack_reference_band(rate_ref_2d, obs_epoch)
	
	###########################
	## Make the cross spectrum
	###########################
	
	cs_segment, mean_rate_segment_ci, mean_rate_segment_ref, power_ci, \
		power_ref = make_cs(rate_ci_2d, rate_ref, n_bins, detchans)
		
	## Sums across segments -- arrays, so it adds by index
	sum_rate_whole_ci += mean_rate_segment_ci
	sum_rate_whole_ref += mean_rate_segment_ref
	sum_power_ci += power_ci
	sum_power_ref += power_ref
	cs_sum += cs_segment
	sum_rate_ci += np.mean(rate_ci_2d)
	
	## Clearing variables from memory for the next iteration
	cs_segment = None
	mean_rate_segment_ci = None
	mean_rate_segment_ref = None
	power_ci = None
	power_ref = None
	rate_ref = None
	rate_ci_2d = None
	rate_ref_2d = None
		
	return cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, sum_power_ci, \
		sum_power_ref, sum_rate_ci
	
## End of function 'each_segment'


################################################################################
def fits_in(in_file, n_bins, detchans, dt, print_iterator, test, obs_epoch):
	"""
			fits_in
	
	Reading in an eventlist in .fits format to make the cross spectrum. Reads 
	in a clock-corrected GTI'd event list, populates the light curves, computes 
	cross spectrum per energy channel and keeps running average of the cross 
	spectra.
	
	I take the approach: start time <= segment < end_time, to avoid double-
	counting and/or skipping events.
	
	"""
	
	#######################################################
	## Check if the FITS file exists; if so, load the data
	#######################################################
	
	try:
		fits_hdu = fits.open(in_file)
	except IOError:
		print "\tERROR: File does not exist: %s" % in_file
		sys.exit()
	
	header = fits_hdu[0].header	 ## Header info is in ext 0, data is in ext 1
	data = fits_hdu[1].data
	fits_hdu.close()
	
	###################
	## Initializations
	###################
	
	num_segments = 0
	sum_rate_whole_ci = np.zeros(detchans, dtype=np.float64)
	sum_rate_whole_ref = 0
	cs_sum = np.zeros((n_bins, detchans), dtype=np.complex128)
	sum_power_ci = np.zeros((n_bins, detchans), dtype=np.float64)
	sum_power_ref = np.zeros(n_bins, dtype=np.float64)
	sum_rate_ci = 0
	
	start_time = data.field('TIME')[0]
	final_time = data.field('TIME')[-1]
	seg_end_time = start_time + (dt * n_bins)
	
	###################################
	## Selecting PCU for interest band
	###################################
	
	PCU2_mask = data.field('PCUID') == 2
	data_pcu2 = data[PCU2_mask]
	all_time_ci = np.asarray(data_pcu2.field('TIME'), dtype=np.float64)
	all_energy_ci = np.asarray(data_pcu2.field('CHANNEL'), dtype=np.float64)
	
	## Determining the next-most-prevalent pcu
# 	pcus_on, occurrences = np.unique(data.field('PCUID'), return_counts=True)
# 	## Getting rid of the PCU 2 element (since we don't want that as the ref)
# 	pcu2 = np.where(pcus_on == 2)
# 	pcus_on = np.delete(pcus_on, pcu2)
# 	occurrences = np.delete(occurrences, pcu2)
# 	## Next-most pcu is:
# 	most_pcu = np.argmax(occurrences)
# 	ref_pcu = pcus_on[most_pcu]
# 	print "Ref PCU =", ref_pcu
	
	####################################
	## Selecting PCU for reference band
	####################################
	
# 	refpcu_mask = data.field('PCUID') == ref_pcu
	refpcu_mask = data.field('PCUID') != 2
	data_ref = data[refpcu_mask]
	all_time_ref = np.asarray(data_ref.field('TIME'), dtype=np.float64)
	all_energy_ref = np.asarray(data_ref.field('CHANNEL'), dtype=np.float64)
	
	############################
	## Looping through segments
	############################
	
	while seg_end_time < final_time:
		
		## Get events for channels of interest 
		time_ci = all_time_ci[np.where(all_time_ci < seg_end_time)]
		energy_ci = all_energy_ci[np.where(all_time_ci < seg_end_time)]
		
		## Chop current segment off the rest of the list
		for_next_iteration_ci = np.where(all_time_ci >= seg_end_time)
		all_time_ci = all_time_ci[for_next_iteration_ci]
		all_energy_ci = all_energy_ci[for_next_iteration_ci]
		
		## Get events for reference band
		time_ref = all_time_ref[np.where(all_time_ref < seg_end_time)]
		energy_ref = all_energy_ref[np.where(all_time_ref < seg_end_time)]
		
		## Chop current segment off the rest of the list
		for_next_iteration_ref = np.where(all_time_ref >= seg_end_time)
		all_time_ref = all_time_ref[for_next_iteration_ref]
		all_energy_ref = all_energy_ref[for_next_iteration_ref]
		
		###########################
		## At the end of a segment
		###########################
		
		if len(time_ci) > 0 and len(time_ref) > 0:
			
			num_segments += 1
			
			cs_sum, sum_rate_whole_ci, sum_rate_whole_ref,  sum_power_ci, \
				sum_power_ref, sum_rate_ci = each_segment(time_ci, time_ref, \
				energy_ci, energy_ref, n_bins, detchans, dt, start_time, \
				seg_end_time, obs_epoch, sum_rate_whole_ci, sum_rate_whole_ref,\
				sum_power_ci, sum_power_ref, cs_sum, sum_rate_ci)

			if num_segments % print_iterator == 0:
				print "\t", num_segments
			if test is True and num_segments == 1:  # For testing
				break
	
			start_time += (n_bins * dt)
			seg_end_time += (n_bins * dt)
		
		## This next bit deals with gappy data
		elif len(time_ci) == 0 and len(time_ref) == 0:
		
			start_time = min(all_time_ci[0], all_time_ref[0])
			seg_end_time = start_time + (n_bins * dt)
		
		else:
		 	start_time += (n_bins * dt)
			seg_end_time += (n_bins * dt)
			
		## End of 'if there are counts in this segment'

	## End of while-loop
	
	return cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
		sum_power_ci, sum_power_ref, sum_rate_ci
    
## End of function 'fits_in'


################################################################################
def dat_in(in_file, n_bins, detchans, dt, print_iterator, test, obs_epoch):
    """
            dat_in

    Reading in the eventlist to make the cross spectrum. Reads in a clock- 	
    corrected GTI'd event list, populates the light curves, computes cross 
    spectrum per energy channel and keeps running average of the cross spectrum,
    mean count rate for interest band and reference band, and power spectra for
    the interest and reference bands.
    
    """
    ###################
    ## Initializations
    ###################
    
    num_segments = 0
    sum_rate_whole_ci = np.zeros(detchans, dtype=np.float64)
    sum_rate_whole_ref = 0
    cs_sum = np.zeros((n_bins, detchans), dtype=np.complex128)
    time_ci = np.asarray([])
    energy_ci = np.asarray([])
    time_ref = np.asarray([])
    energy_ref = np.asarray([])
    start_time = -99
    sum_power_ci = np.zeros((n_bins, detchans), dtype=np.float64)
    sum_power_ref = np.zeros(n_bins, dtype=np.float64)
    sum_rate_ci = 0

    ## Reading only the first line of data to get the start time of the file
    try:
		with open(in_file, 'r') as fo:
			for line in fo:
				if line[0].strip() != "#":
					line = line.strip().split()
					start_time = np.float64(line[0])
					break
    except IOError:
    	print "\tERROR: File does not exist: %s" % in_file
		
    end_time = start_time + (dt * n_bins)
	
	############################
	## Looping through segments
	############################
	
    with open(in_file, 'r') as f:
        for line, next_line in tools.pairwise(f):
            if line[0].strip() != "#":  # If the line is not a comment
            
                line = line.strip().split()
                next_line = next_line.strip().split()
                current_time = np.float64(line[0])
                current_chan = np.uint16(line[1])
                current_pcu = np.uint8(line[2])
                next_time = np.float64(next_line[0])
                
				########################################
				## Selecting band for data based on PCU
				########################################
				
                if current_pcu == 2:  ## If pcu = 2
                    time_ci = np.append(time_ci, current_time)
                    energy_ci = np.append(energy_ci, current_chan)
                else:
                    time_ref = np.append(time_ref, current_time)
                    energy_ref = np.append(energy_ref, current_chan)

                next_time = float(next_line[0])
                next_end_time = end_time + (dt * n_bins)
				
				###########################
				## At the end of a segment
				###########################
		
                if next_time >= end_time:
                	if len(time_ci) > 0 and len(time_ref) > 0:
                	
                		num_segments += 1
                		
                		cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, \
                			sum_power_ci, sum_power_ref, sum_rate_ci = \
                			each_segment(time_ci, time_ref, energy_ci, \
                			energy_ref, n_bins, detchans, dt, start_time, \
                			end_time, obs_epoch, sum_rate_whole_ci, \
                			sum_rate_whole_ref, sum_power_ci, sum_power_ref, \
                			cs_sum, sum_rate_ci)
                			
                		if num_segments % print_iterator == 0:
                			print "\t", num_segments
                		if test is True and num_segments == 1:  # For testing
                			break
                		time_ci = []
                		energy_ci = []
                		time_ref = []
                		energy_ref = []
                		
                		start_time += (n_bins * dt)
                		end_time += (n_bins * dt)
                	
                	## This next bit helps it handle gappy data
                	elif len(time_ci) == 0 and len(time_ref) == 0:
                		start_time = next_time
                		end_time = start_time + (n_bins * dt)
                	
                	else:
                		start_time += (n_bins * dt)
                		end_time += (n_bins * dt)
                	## End of 'if there are counts in this segment'
				## End of 'if it`s at the end of a segment'
				
			## End of 'if the line is not a comment'
		## End of for-loop
	## End of with-block
		
    return cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
        sum_power_ci, sum_power_ref, sum_rate_ci
        
## End of function 'dat_in'
    
    
################################################################################
def read_and_use_segments(in_file, n_bins, detchans, dt, test):
    """
			read_and_use_segments
	
	Reads in segments of a light curve from a .dat or .fits file. Split into 
	'fits_in' and 'dat_in' for easier readability.
		
    """
	
    assert tools.power_of_two(n_bins), "ERROR: n_bins must be a power of 2."
  	
    print "Input file: %s" % in_file
  	
  	## Determining print iterator for segments
    if n_bins == 32768:
    	print_iterator = int(10)
    elif n_bins < 32768:
        print_iterator = int(10)
    else:
    	print_iterator = int(1)
  	
    print "Segments computed:"
	
	#######################################################
	## Data is read in differently from fits and dat files
	#######################################################
	
    if (in_file[-5:].lower() == ".fits"):
	
        obs_epoch = tools.obs_epoch_rxte(in_file)
		
        cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
			sum_power_ci, sum_power_ref, sum_rate_ci = fits_in(in_file, \
			n_bins, detchans, dt, print_iterator, test, obs_epoch)
		
    elif (in_file[-4:].lower() == ".dat"):
	
        fits_file = in_file[0:-4] + ".fits"  ## Still need a fits file to get the observation time
        obs_epoch = tools.obs_epoch_rxte(fits_file)
  		
        cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
			sum_power_ci, sum_power_ref, sum_rate_ci = dat_in(in_file, \
			n_bins, detchans, dt, print_iterator, test, obs_epoch)
			
    else:
	    raise Exception("ERROR: Input file type not recognized. Must be .dat or\
.fits.")
	
    return cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, \
		sum_power_ci, sum_power_ref, sum_rate_ci
        
## End of function 'read_and_use_segments'


################################################################################
def get_background(bkgd_file):
	"""
			get_background
	
	Get the background count rate from a background spectrum file.
	
	"""
	
	try:
		fits_hdu = fits.open(bkgd_file)
	except IOError:
		print "\tERROR: File does not exist: %s" % bkgd_file
		sys.exit()
		
	header = fits_hdu[1].header
	data = fits_hdu[1].data
	fits_hdu.close()
	
	exposure = float(header['EXPOSURE'])	
	counts = data.field('COUNTS')
	
	rate = counts / exposure
	
	return rate
	
## End of function 'get_background'


################################################################################
def main(in_file, out_file, bkgd_file, num_seconds, dt_mult, test, filter):
    """
            main

    Reads in one event list, splits into two light curves, makes segments and
    populates them to give them length n_bins, computes the cross spectrum of
    each segment per energy channel and then averaged cross spectrum of all the
    segments per energy channel, and then computes the cross-correlation
    function (ccf) per energy channel.

    """
    
    #####################################################
    ## Idiot checks, to ensure that our assumptions hold
    #####################################################
    
    assert num_seconds > 0, "ERROR: num_seconds must be a positive integer."
    assert dt_mult >= 1, "ERROR: dt_mult must be a positive integer."
	
	##################################################
	## Initializations; 'whole' is over one data file
	##################################################
	
    t_res = float(tools.get_key_val(in_file, 0, 'TIMEDEL'))
    dt = dt_mult * t_res
    n_bins = num_seconds * int(1.0 / dt)
    nyquist_freq = 1.0 / (2.0 * dt)
    detchans = float(tools.get_key_val(in_file, 0, 'DETCHANS'))
    
    print "\nDT = %f" % dt
    print "N_bins = %d" % n_bins
    print "Nyquist freq =", nyquist_freq
    print "Filtering?", filter
	
	###################################################################
    ## Reading in the background count rate from a background spectrum
    ###################################################################
    
    if bkgd_file:
		bkgd_rate = get_background(bkgd_file)
    else:
		bkgd_rate = np.zeros(detchans)	
    print " "
    
	#################################################
	## Reading in data, computing the cross spectrum
	#################################################
	
    cs_sum, sum_rate_whole_ci, sum_rate_whole_ref, num_segments, sum_power_ci, \
    	sum_power_ref, sum_rate_ci = read_and_use_segments(in_file, n_bins, \
    	detchans, dt, test)
        
	#########################################
	## Turning sums over segments into means
	#########################################
	
    mean_rate_whole_ci = sum_rate_whole_ci / float(num_segments)
    mean_rate_whole_ref = sum_rate_whole_ref / float(num_segments)
    mean_power_ci = sum_power_ci / float(num_segments)
    mean_power_ref = sum_power_ref / float(num_segments)
    cs_avg = cs_sum / float(num_segments)
        
    ################################################################
    ## Printing the cross spectrum to a file, for plotting/checking
    ################################################################
    
    cs_out = np.column_stack((fftpack.fftfreq(n_bins, d=dt), cs_avg))
    np.savetxt('cs_avg.dat', cs_out)
    
    ##################################################################
    ## Subtracting the background count rate from the mean count rate
    ##################################################################
    
    mean_rate_whole_ci -= bkgd_rate
    
    ## Need to use a background from ref. PCU for the reference band...
    ref_bkgd_rate = np.mean(bkgd_rate[2:26])
    mean_rate_whole_ref -= ref_bkgd_rate    
    
    ######################
    ## Making lag spectra
    ######################
    
    make_lags(out_file, in_file, dt, n_bins, detchans, num_seconds, \
    	num_segments, mean_rate_whole_ci, mean_rate_whole_ref, cs_avg, \
    	mean_power_ci, mean_power_ref)
	
	##############################################
	## Computing ccf from cs, and computing error
	##############################################
	
    if filter:
    	ccf_end, ccf_error = FILT_cs_to_ccf_w_err(cs_avg, dt, n_bins, \
    		detchans, num_seconds, num_segments, mean_rate_whole_ci, \
    		mean_rate_whole_ref, mean_power_ci, mean_power_ref, True)
    else:
    	ccf_end, ccf_error = UNFILT_cs_to_ccf_w_err(cs_avg, dt, n_bins, \
    		detchans, num_seconds, num_segments, mean_rate_whole_ci, \
    		mean_rate_whole_ref, mean_power_ci, mean_power_ref, True)
	
    print "Number of segments:", num_segments
    print "Sum of mean rate for ci:", np.sum(mean_rate_whole_ci)
    print "Mean rate for ref:", np.mean(mean_rate_whole_ref)
	
    t = np.arange(0, n_bins)
	
	##########
    ## Output
    ##########
    
    fits_out(out_file, in_file, bkgd_file, dt, n_bins, detchans, num_seconds, \
    	num_segments, mean_rate_whole_ci, mean_rate_whole_ref, t, ccf_end, \
    	ccf_error, filter)
	
## End of function 'main'


################################################################################
if __name__ == "__main__":
	
	##############################################
	## Parsing input arguments and calling 'main'
	##############################################
	
    parser = argparse.ArgumentParser(usage="python ccf.py infile outfile [-b \
BKGD_SPECTRUM] [-n NUM_SECONDS] [-m DT_MULT] [-t {0,1}] [-f {0,1}]", \
description="Computes the cross-correlation function of a channel of interest \
with a reference band.", epilog="For optional arguments, default values are \
given in brackets at end of description.")

    parser.add_argument('infile', help='Name of (ASCII/txt/dat) event list \
containing both the reference band and the channels of interest. Assumes ref \
band = PCU 0, interest = PCU 2.')

    parser.add_argument('outfile', help='The full path of the (ASCII/txt) file\
 to write the frequency and power to.')

    parser.add_argument('-b', '--bkgd', required=False, dest='bkgd_file', \
help='Name of the background spectrum (in pha/fits format).')
    	
    parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two, \
default=1, dest='num_seconds', help='Number of seconds in each Fourier segment.\
 Must be a power of 2, positive, integer. [1]')
    	
    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two, \
default=1, dest='dt_mult', help='Multiple of dt (dt is from data file) for \
timestep between bins. Must be a power of 2, positive, integer. [1]')
    	
    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
dest='test', help='Int flag: 0 if computing all segments, 1 if only computing \
one segment for testing. [0]')
        
    parser.add_argument('-f', '--filter', type=int, default=0, choices={0,1},
dest='filter', help='Int flag: 0 if NOT applying a filter in frequency-space, \
1 if applying frequency filter (around a pulsation). [0]')
        
    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True
        
    filter = False
    if args.filter == 1:
    	filter = True
    
    main(args.infile, args.outfile, args.bkgd_file, args.num_seconds, args.dt_mult, test, filter)


## End of program 'ccf.py'
################################################################################
