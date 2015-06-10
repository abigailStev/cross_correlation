import argparse
import numpy as np
import sys
from scipy import fftpack
from datetime import datetime
import os.path
import subprocess
from astropy.io import fits
import tools  # https://github.com/abigailStev/whizzy_scripts

__author__ = "Abigail Stevens, A.L.Stevens at uva.nl"

"""
Computes the cross-correlation function of narrow energy channels of interest
with a broad energy reference band from RXTE event-mode data.

Abigail Stevens, A.L.Stevens at uva.nl, 2014-2015

"""

class PSD(object):
    def __init__(self):
        self.power = np.asarray([])
        self.noise = 0
        self.variance = 0
        self.rms = 0

class Lightcurve(object):
    def __init__(self):
        self.mean_rate = 0
        self.fullpower = np.asarray([])
        self.raw = PSD()
        self.absrms = PSD()
        self.fracrms = PSD()


################################################################################
def fits_out(out_file, in_file, bkgd_file, param_dict, mean_rate_ci, \
    mean_rate_ref, t, ccf, ccf_error, filter):
    """
    Writes the cross-correlation function to a .fits output file.

    """

    ## Getting data into a good output structure
    chan = np.arange(0, param_dict['detchans'])
    energy_channels = np.tile(chan, len(t))
    if filter:
        ccf_error = np.tile(ccf_error, len(t))
    else:
        ccf_error = ccf_error.flatten('C')
    time_bins = np.repeat(t, param_dict['detchans'])
    assert len(energy_channels) == len(time_bins)

    print "\nOutput sent to: %s" % out_file

    ## Making FITS header (extension 0)
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross-correlation function")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file)
    prihdr.set('BKGD', bkgd_file)
    prihdr.set('DT', param_dict['dt'], "seconds")
    prihdr.set('N_BINS', param_dict['n_bins'], "time bins per segment")
    prihdr.set('SEGMENTS', param_dict['num_seg'], "segments in the whole light"\
        " curve")
    prihdr.set('SEC_SEG', param_dict['num_seconds'], "seconds per segment")
    prihdr.set('EXPOSURE', param_dict['num_seg'] * param_dict['num_seconds'],
        "seconds, of light curve")
    prihdr.set('DETCHANS', param_dict['detchans'], "Number of detector energy"\
        " channels")
    prihdr.set('RATE_CI', str(mean_rate_ci.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_ref, "counts/second")
    prihdr.set('FILTER', str(filter))
    prihdu = fits.PrimaryHDU(header=prihdr)

    ## Making FITS table (extension 1)
    col1 = fits.Column(name='TIME_BIN', format='K', array=time_bins)
    col2 = fits.Column(name='CCF', unit='Counts/second', format='D',
        array=ccf.flatten('C'))
    col3 = fits.Column(name='ERROR', unit='', format='D',
        array=ccf_error)
    col4 = fits.Column(name='CHANNEL', unit='', format='I',
        array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3, col4])
    tbhdu = fits.BinTableHDU.from_columns(cols)

    ## If the file already exists, remove it (still working on just updating it)
    assert out_file[-4:].lower() == "fits", \
        'ERROR: Output file must have extension ".fits".'
    if os.path.isfile(out_file):
        subprocess.call(['rm', out_file])

    ## Writing to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file)


################################################################################
def get_phase_err(cs_avg, power_ci, power_ref, n, M):
    """
    Computes the error on the complex phase (in radians). Power should not be
    noise-subtracted.

    """
    with np.errstate(all='ignore'):
        a = power_ci * power_ref
        coherence = np.where(a != 0, np.abs(cs_avg)**2 / a, 0)
        print np.shape(coherence)
        phase_err = np.sqrt(np.where(coherence != 0, (1 - coherence) / \
            (2 * coherence * n * M), 0))

    return phase_err


################################################################################
def phase_to_tlags(phase, f, detchans):
    """
    Converts a complex phase (in radians) to a time lag (in seconds).

    """
    assert np.shape(phase) == np.shape(f), "ERROR: Phase must have same \
dimensions as f."

    with np.errstate(all='ignore'):
        tlags =  np.where(f != 0, phase / (2.0 * np.pi * f), 0)

    return tlags


################################################################################
def make_lags(out_file, in_file, dt, n_bins, detchans, num_seconds,
    num_seg, mean_rate_ci_whole, mean_rate_ref_whole, cs_avg, power_ci,
    power_ref):
    """
    Computes the phase lag and time lag from the average cross spectrum, and
    writes the lag information to a .fits output file.

    """
    assert np.shape(power_ci) == (n_bins, detchans)
    assert np.shape(power_ref) == (n_bins, )
    assert np.shape(cs_avg) == (n_bins, detchans)

# 	low_freq = 4.47
# 	hi_freq = 6.35
    low_freq = 4.0
    hi_freq = 7.0

    ## Getting the Fourier frequencies for the cross spectrum
    freq = fftpack.fftfreq(n_bins, d=dt)

    ## Only keeping the parts associated with positive Fourier frequencies
    nyq_ind = np.argmax(freq)+1  ## because in python, the scipy fft makes the
            ## nyquist frequency negative, and we want it to be positive! (it is
            ## actually both pos and neg)
    freq = np.abs(freq[1:nyq_ind + 1])  ## because it slices at end-1, and we
            ## want to include 'nyq_ind'; abs is because the nyquist freq is
            ## both pos and neg, and we want it pos here. but we don't want
            ## freq=0 because that gives errors.
    cs_avg = cs_avg[1:nyq_ind + 1, ]
    power_ci = power_ci[1:nyq_ind + 1, ]
    power_ref = power_ref[1:nyq_ind + 1]

    ## Getting lag and error for lag-frequency plot
    phase = -np.arctan2(cs_avg.imag, cs_avg.real) ## Negative sign is so that a
            ## positive lag is a hard energy lag?
    err_phase = get_phase_err(cs_avg, power_ci, np.resize(np.repeat(power_ref, \
            detchans), np.shape(power_ci)), 1, num_seg)
    print np.shape(err_phase)
    f = np.resize(np.repeat(freq, detchans), (len(freq), detchans))
    tlag = phase_to_tlags(phase, f, detchans)
    err_tlag = phase_to_tlags(err_phase, f, detchans)

    ## Getting lag and error for lag-energy plot (averaging over frequencies)
    f_span_low = np.argmax(freq == low_freq)
    f_span_hi = np.argmax(freq == hi_freq)
    f_span = f_span_hi - f_span_low
    print "Fspan low:", f_span_low
    print "Fspan hi:", f_span_hi
    print "F span:", f_span
    frange_freq = freq[f_span_low:f_span_hi+1]
    frange_cs = np.mean(cs_avg[f_span_low:f_span_hi+1, ], axis=0)
    frange_pow_ci = np.mean(power_ci[f_span_low:f_span_hi+1, ], axis=0)
    frange_pow_ref = np.repeat(np.mean(power_ref[f_span_low:f_span_hi+1]), detchans)
    print "Shape cs:", np.shape(frange_cs)
    e_phase = -np.arctan2(frange_cs.imag, frange_cs.real) ## Negative sign is
            ## so that a positive lag is a hard energy lag?
    e_err_phase = get_phase_err(frange_cs, frange_pow_ci, frange_pow_ref, f_span, num_seg)
    print "Shape err phase:", np.shape(e_err_phase)
    f = np.repeat(np.mean(frange_freq), detchans)
    print "Shape phase:", np.shape(e_phase)
    print "Shape f:", np.shape(f)
    e_tlag = phase_to_tlags(e_phase, f, detchans)
    print "Shape tlag:", np.shape(e_tlag)
    e_err_tlag = phase_to_tlags(e_err_phase, f, detchans)
    print "Shape err tlag:", np.shape(e_err_tlag)

    chan = np.arange(0, detchans)
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
    prihdr.set('SEGMENTS', num_seg, "segments in the whole light curve")
    prihdr.set('EXPOSURE', num_seg * n_bins * dt,
        "seconds, of light curve")
    prihdr.set('DETCHANS', detchans, "Number of detector energy channels")
    prihdr.set('LOWFREQ', low_freq)
    prihdr.set('HIGHFREQ', hi_freq)
    prihdr.set('RATE_CI', str(mean_rate_ci_whole.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_ref_whole, "counts/second")
    prihdr.set('FILTER', str(filter))
    prihdu = fits.PrimaryHDU(header=prihdr)

    ## Making FITS table for lag-frequency plot (extension 1)
    col1 = fits.Column(name='FREQUENCY', format='D', array=bins)
    col2 = fits.Column(name='PHASE', unit='radians', format='D',
        array=phase.flatten('C'))
    col3 = fits.Column(name='PHASE_ERR', unit='radians', format='D',
        array=err_phase.flatten('C'))
    col4 = fits.Column(name='TIME_LAG', unit='s', format='D',
        array=tlag.flatten('C'))
    col5 = fits.Column(name='TIME_LAG_ERR', unit='s', format='D',
        array=err_tlag.flatten('C'))
    col6 = fits.Column(name='CHANNEL', unit='', format='I',
        array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6])
    tbhdu1 = fits.BinTableHDU.from_columns(cols)

    ## Making FITS table for lag-energy plot (extension 2)
    col1 = fits.Column(name='PHASE', unit='radians', format='D', array=e_phase)
    col2 = fits.Column(name='PHASE_ERR', unit='radians', format='D', \
        array=e_err_phase)
    col3 = fits.Column(name='TIME_LAG', unit='s', format='D', array=e_tlag)
    col4 = fits.Column(name='TIME_LAG_ERR', unit='s', format='D', \
        array=e_err_tlag)
    col5 = fits.Column(name='CHANNEL', unit='', format='I', \
        array=chan)
    cols = fits.ColDefs([col1, col2, col3, col4, col5])
    tbhdu2 = fits.BinTableHDU.from_columns(cols)

    ## If the file already exists, remove it
    assert out_file[-4:].lower() == "fits", \
        'ERROR: Output file must have extension ".fits".'
    if os.path.isfile(out_file):
        subprocess.call(["rm", out_file])

    ## Writing to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu1, tbhdu2])
    thdulist.writeto(out_file)


################################################################################
def save_for_lags(out_file, in_file, param_dict, ci, ref, cs_avg):
    """
    Saving header data, the cross spectrum, CoI power spectrum, and reference
    band power spectrum to a FITS file to use in the program make_lags.py to get
    cross-spectral lags.

    """
    ## Getting the Fourier frequencies for the cross spectrum
    freq = fftpack.fftfreq(param_dict['n_bins'], d=param_dict['dt'])
    nyq_index = param_dict['n_bins']/2
    assert np.abs(freq[nyq_index]) == param_dict['nyquist']

    ## Only keeping the parts associated with positive Fourier frequencies
    freq = np.abs(freq[0:nyq_index + 1])  ## because it slices at end-1, and we
            ## want to include 'nyq_index'; abs is because the nyquist freq is
            ## both pos and neg, and we want it pos here.
    cs_avg = cs_avg[0:nyq_index + 1, ]
    ci.raw.power = ci.fullpower[0:nyq_index + 1, ]
    ref.raw.power = ref.fullpower[0:nyq_index + 1]

    chan = np.arange(0, param_dict['detchans'])
    energy_channels = np.tile(chan, len(freq))

    out_file = out_file.replace("cross_correlation/out_ccf", "lags/out_lags")
    out_file = out_file.replace(".", "_cs.")
    print "Output sent to: %s" % out_file

    ## Making FITS header (extension 0)
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross spectrum, power spectrum ci, and power spectrum "\
            "ref.")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file)
    prihdr.set('DT', param_dict['dt'], "seconds")
    prihdr.set('N_BINS', param_dict['n_bins'], "time bins per segment")
    prihdr.set('SEGMENTS', param_dict['num_seg'], "segments in the whole light"\
            " curve")
    prihdr.set('EXPOSURE', param_dict['num_seg'] * param_dict['num_seconds'],
            "seconds, of light curve")
    prihdr.set('DETCHANS', param_dict['detchans'], "Number of detector energy "\
            "channels")
    prihdr.set('RATE_CI', str(ci.mean_rate.tolist()), "counts/second")
    prihdr.set('RATE_REF', ref.mean_rate, "counts/second")
    prihdu = fits.PrimaryHDU(header=prihdr)

    ## Making FITS table for cross spectrum
    col1 = fits.Column(name='FREQUENCY', format='D', array=freq)
    col2 = fits.Column(name='CROSS', unit='raw', format='C',
        array=cs_avg.flatten('C'))
    col3 = fits.Column(name='CHANNEL', unit='', format='I',
        array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu1 = fits.BinTableHDU.from_columns(cols)

    ## Making FITS table for power spectrum of channels of interest
    col1 = fits.Column(name='FREQUENCY', format='D', array=freq)
    col2 = fits.Column(name='POWER', unit='raw', format='D',
        array=ci.raw.power.flatten('C'))
    col3 = fits.Column(name='CHANNEL', unit='', format='I',
        array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu2 = fits.BinTableHDU.from_columns(cols)

    ## Making FITS table for power spectrum of reference band
    col1 = fits.Column(name='FREQUENCY', format='D', array=freq)
    col2 = fits.Column(name='POWER', unit='raw', format='D',
        array=ref.raw.power)
    cols = fits.ColDefs([col1, col2])
    tbhdu3 = fits.BinTableHDU.from_columns(cols)

    ## If the file already exists, remove it
    assert out_file[-4:].lower() == "fits", \
        'ERROR: Output file must have extension ".fits".'
    if os.path.isfile(out_file):
        subprocess.call(["rm", out_file])

    ## Writing to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu1, tbhdu2, tbhdu3])
    thdulist.writeto(out_file)


################################################################################
def find_pulse_freq(freq, power_ref):
    """
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


################################################################################
def filter_freq(freq_space_array, dt, n_bins, detchans, power_ref):
    """
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
    zero_end = np.zeros((len(freq_space_array) - j_max, detchans),
        dtype=np.complex128)

    ## Concatenate the arrays together
    filt_freq_space_array = np.concatenate((zero_front,
        freq_space_array[j_min:j_max, :], zero_end), axis=0)

    ## Check that the original array is the same shape as the filtered one
    assert np.shape(freq_space_array) == np.shape(filt_freq_space_array), \
        "ERROR: Frequency-filtered cross spectrum does not have the same size \
as the original cross spectrum. Something went wrong."

    return filt_freq_space_array, j_min, j_max


################################################################################
def FILT_cs_to_ccf_w_err(cs_avg, param_dict, ci, ref, noisy):
    """
    Filters the cross-spectrum in frequency space, takes the iFFT of the
    filtered cross spectrum to get the cross-correlation function, and computes
    the error on the cross-correlation function. Note that error is NOT
    independent between time bins due to the filtering! But is still independent
    between energy bins.

    """
    ## Filter the cross spectrum in frequency
    filtered_cs_avg, j_min, j_max = filter_freq(cs_avg, param_dict['dt'], \
        param_dict['n_bins'], param_dict['detchans'], ref.fullpower)

    ## Absolute rms norms of poisson noise
    ci.absrms.noise = 2.0 * ci.mean_rate
    ref.absrms.noise = 2.0 * ref.mean_rate

    ## If there's no noise in a (simulated) power spectrum, noise level = 0
    if not noisy:
        ci.absrms.noise = np.zeros(param_dict['detchans'])
        ref.absrms.noise = 0

    noise_ref_array = np.repeat(noise_ref, param_dict['detchans'])

    ## Extracting only the signal frequencies of the mean powers
    signal_ci_pow = np.float64(ci.fullpower[j_min:j_max, :])
    signal_ref_pow = np.float64(ref.fullpower[j_min:j_max])

    ## Putting powers into absolute rms2 normalization, subtracting noise
    signal_ci_pow = signal_ci_pow * (2.0 * param_dict['dt'] / \
            float(param_dict['n_bins'])) - ci.absrms.noise
    signal_ref_pow = signal_ref_pow * (2.0 * param_dict['dt'] / \
            float(param_dict['n_bins'])) - ref.absrms.noise

    ## Getting rms of reference band, to normalize the ccf
    ref.absrms.variance = np.sum(signal_ref_pow * param_dict['df'])
    print "Reference band variance:", ref.absrms.variance
    rms.absrms.ref = np.sqrt(ref.absrms.variance)
    print "Frac RMS of reference band:", rms.absrms.ref / ref.mean_rate
    ## in frac rms units here -- should be few percent

    ## Broadcasting signal_ref_pow into same shape as signal_ci_pow
    signal_ref_pow = np.resize(np.repeat(signal_ref_pow, param_dict['detchans']),
        np.shape(signal_ci_pow))
    assert np.shape(signal_ref_pow) == np.shape(signal_ci_pow)

    temp = (ci.absrms.noise * signal_ref_pow) + \
        (ref.absrms.noise * signal_ci_pow) + \
        (ci.absrms.noise * ref.absrms.noise)
    cs_noise_amp = np.sqrt(np.sum(temp, axis=0) / float(param_dict['num_seg']))

    temp1 = np.absolute(cs_avg[j_min:j_max, :]) * (2.0 * param_dict['dt'] / float(param_dict['n_bins']))
    cs_signal_amp = np.sum(temp1, axis=0)

    ## Assuming that cs_noise_amp and cs_signal_amp are float arrays, size 64
    with np.errstate(all='ignore'):
        error_ratio = np.where(cs_signal_amp != 0, cs_noise_amp / cs_signal_amp, 0)

    ## Taking the IFFT of the cross spectrum to get the CCF
    ccf_end = fftpack.ifft(filtered_cs_avg, axis=0)

    ## Dividing ccf by rms of signal in reference band
    ccf_end *= (2.0 / float(param_dict['n_bins']) / rms_ref)

    ## Computing the error on the ccf
    ccf_rms_ci = np.sqrt(np.var(ccf_end, axis=0, ddof=1))
    ccf_error = ccf_rms_ci * error_ratio

    return ccf_end, ccf_error


################################################################################
def standard_ccf_err(param_dict):
    """
    Computes the standard error on each ccf bin from the segment-to-segment
    variations. Use this for *UNFILTERED* CCFs.

    S. Vaughan 2013, "Scientific Inference", equations 2.3 and 2.4.

    """
    standard_err = np.zeros((200, param_dict['detchans']))
    for i in range(param_dict['detchans']):
        in_file = "./out_ccf/ccf_segs_" + str(i) + ".dat"
        table_i = np.loadtxt(in_file)
        mean_ccf_i = np.mean(table_i, axis=0)
        ccf_resid_i = table_i - mean_ccf_i
        sample_var_i = np.sum(ccf_resid_i**2, axis=0) / float(param_dict['num_seg']-1)  ## eq 2.3
        standard_err_i = np.sqrt(sample_var_i/float(param_dict['num_seg']))  ## eq 2.4
        standard_err[:,i] = standard_err_i
    print(standard_err)
    return standard_err


################################################################################
def UNFILT_cs_to_ccf_w_err(cs_avg, param_dict, ci, ref, noisy):
    """
    Takes the iFFT of the cross spectrum to get the cross-correlation function,
    and computes the error on the cross-correlation function. This error is
    not correlated between energy bins but may be correlated between time bins.

    """

    if len(ref.fullpower) == param_dict['n_bins']:
        freq = fftpack.fftfreq(param_dict['n_bins'], d=param_dict['dt'])
        nyq_index = param_dict['n_bins'] / 2
        assert np.abs(freq[nyq_index]) == param_dict['nyquist']
        ref.raw.power = ref.fullpower[0:nyq_index+1]
        ci.raw.power = ci.fullpower[0:nyq_index+1, ]

# 	print "CI count rate:", ci.mean_rate
# 	print "Ref count rate:", ref.mean_rate

    ######################################################
    ## Take the IFFT of the cross spectrum to get the CCF
    ######################################################

    ccf = fftpack.ifft(cs_avg, axis=0).real

    ## Absolute rms norms of Poisson noise
    ## Was over-estimating the noise from the count rate; using mean power from
    ## higher frequencies as the noise level

    if noisy:
        ci.absrms.noise = 2.0 * ci.mean_rate #* 0.985
        ref.absrms.noise = 2.0 * ref.mean_rate #* 0.985
    else:
        ci.absrms.noise = 0.0
        ref.absrms.noise = 0.0

    ## Putting powers into absolute rms2 normalization, subtracting noise
    ci.absrms.power = ci.raw.power * (2.0 * param_dict['dt'] / \
            float(param_dict['n_bins'])) - ci.absrms.noise
    ref.absrms.power = ref.raw.power * (2.0 * param_dict['dt'] / \
            float(param_dict['n_bins'])) - ref.absrms.noise

# 	## Getting rms of reference band, to normalize the ccf and acf
    ref.absrms.variance = np.sum(ref.absrms.power * param_dict['df'])
    ref.absrms.rms = np.sqrt(ref.absrms.variance)
    print "Ref band var:", ref.absrms.variance, "(abs rms)"
    print "Ref band rms:", ref.absrms.rms, "(abs rms)"
    print "Ref band var:", ref.absrms.variance / (ref.mean_rate ** 2), "(frac rms)"
    print "Ref band rms:", ref.absrms.rms / ref.mean_rate, "(frac rms)"
#
# 	## Computing the autocorrelation functions from the power spectra
# 	acf_ci = fftpack.ifft(power_ci).real
# 	acf_ref = fftpack.ifft(power_ref).real
# 	acf_ci *= (2.0 / float(n_bins) / absrms_rms_ref)
# 	acf_ref *= (2.0 / float(n_bins) / absrms_rms_ref)
#
# 	## Broadcasting acf_ref into same shape as acf_ci
# 	acf_ref = np.resize(np.repeat(acf_ref, detchans), np.shape(acf_ci))
# 	assert np.shape(acf_ref) == np.shape(acf_ci), "ERROR: Array broadcasting \
# failed."
#
# 	## Bartlett formula for computing variance on unfilfered CCF
# 	## Box & Jenkens 1976 eqn 11.1.9; Smith & Vaughan 2007 eqn 3
# 	var_ccf = np.abs((acf_ci * acf_ref) / n_bins)
# 	rms_ccf = np.sqrt(var_ccf)

    ## Dividing ccf by rms of signal in reference band
    ccf *= (2.0 / float(param_dict['n_bins']) / ref.absrms.rms)

    ccf_err = standard_ccf_err(param_dict)

    print "CCF:", ccf[2:7, 6]
    print "Err:", ccf_err[2:7, 6]
    print "Countrate ci chan 6:", ci.mean_rate[6]
    print "Countrate ref:", ref.mean_rate

    return ccf, ccf_err


################################################################################
def stack_reference_band(rate_ref_2d, obs_epoch):
    """
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


################################################################################
def make_cs(rate_ci, rate_ref, param_dict):
    """
    Generating the cross spectrum for one segment of the light curve.

    """
    assert np.shape(rate_ci) == (param_dict['n_bins'], param_dict['detchans']),\
        "ERROR: CoI light curve has wrong dimensions. Must have size (n_bins, "\
        "detector channels)."
    assert np.shape(rate_ref) == (param_dict['n_bins'], ), "ERROR: Reference "\
        "light curve has wrong dimensions. Must have size (n_bins, )."

    ci_seg = Lightcurve()
    ref_seg = Lightcurve()

    ## Computing the mean count rate of the segment
    ci_seg.mean_rate = np.mean(rate_ci, axis=0)
    ref_seg.mean_rate = np.mean(rate_ref)

    ## Subtracting the mean off each value of 'rate'
    rate_sub_mean_ci = np.subtract(rate_ci, ci_seg.mean_rate)
    rate_sub_mean_ref = np.subtract(rate_ref, ref_seg.mean_rate)

    ## Taking the FFT of the time-domain photon count rate
    ## SciPy is faster than NumPy or pyFFTW for my array sizes
    fft_data_ci = fftpack.fft(rate_sub_mean_ci, axis=0)
    fft_data_ref = fftpack.fft(rate_sub_mean_ref)

    ## Computing the power from the fourier transform
    ci_seg.fullpower = np.absolute(fft_data_ci) ** 2
    ref_seg.fullpower = np.absolute(fft_data_ref) ** 2

    ## Broadcasting fft of ref into same shape as fft of ci
    fft_data_ref = np.resize(np.repeat(fft_data_ref, param_dict['detchans']), \
        (param_dict['n_bins'], param_dict['detchans']))

    ## Computing the cross spectrum from the fourier transform
    cs_seg = np.multiply(fft_data_ci, np.conj(fft_data_ref))

    return cs_seg, ci_seg, ref_seg


################################################################################
def print_seg_ccf(param_dict, meanrate_ref, power_ref, cs_seg):
    """
    Printing the first 200 values of each segment of the ccf so that we can
    compute the standard error on the mean ccf.

    """
    absrms_noise_ref = 2.0 * meanrate_ref * 0.985
    absrms_power_ref = power_ref * (2.0 * param_dict['dt'] / \
        float(param_dict['n_bins'])) - absrms_noise_ref
    absrms_var_ref = np.sum(absrms_power_ref * param_dict['df'])
# 	print "Ref band var:", absrms_var_ref, "(abs rms)"

    ccf = fftpack.ifft(cs_seg, axis=0).real
    absrms_rms_ref = np.sqrt(absrms_var_ref)
    ccf *= (2.0 / float(param_dict['n_bins']) / absrms_rms_ref)

    if not np.isnan(absrms_rms_ref):
        for i in range(0, param_dict['detchans']):
            out_file = "./out_ccf/ccf_segs_" + str(i) + ".dat"

            with open(out_file, 'a') as out:
                for element in ccf[0:200,i]:
                    out.write("%.6e\t" % element)
                out.write("\n")


################################################################################
def each_segment(time_ci, time_ref, energy_ci, energy_ref, param_dict,\
    start_time, ci_file, ref_file, cs_sum, sum_rate_ci):
    """
    Turns the event list into a populated histogram, stacks the reference band,
    and makes the cross spectrum, per segment of light curve.

    """
    assert len(time_ci) == len(energy_ci)
    assert len(time_ref) == len(energy_ref)

    ##############################################################
    ## Populate the light curves for interest and reference bands
    ##############################################################

    rate_ci_2d = tools.make_2Dlightcurve(np.asarray(time_ci),
        np.asarray(energy_ci), param_dict['n_bins'], param_dict['detchans'],
        param_dict['dt'], start_time)
    rate_ref_2d = tools.make_2Dlightcurve( np.asarray(time_ref),
        np.asarray(energy_ref), param_dict['n_bins'], param_dict['detchans'],
        param_dict['dt'], start_time)

    ## Stack the reference band
    rate_ref = stack_reference_band(rate_ref_2d, param_dict['obs_epoch'])

    ## Save the reference band light curve to a text file
# 	out_file="./GX339-BQPO_ref_lc.dat"
# 	f_handle = file(out_file, 'a')
# 	np.savetxt(f_handle, rate_ref)
# 	f_handle.close()

    ###########################
    ## Make the cross spectrum
    ###########################

    cs_seg, ci_seg, ref_seg = make_cs(rate_ci_2d, rate_ref, param_dict)

    #####################################################
    ## Printing ccf to a file to later get error for ccf
    #####################################################
    print_seg_ccf(param_dict, ref_seg.mean_rate, ref_seg.fullpower, cs_seg)

    ## Sums across segments -- arrays, so it adds by index
    ci_file.mean_rate += ci_seg.mean_rate
    ref_file.mean_rate += ref_seg.mean_rate
    ci_file.fullpower += ci_seg.fullpower
    ref_file.fullpower += ref_seg.fullpower
    cs_sum += cs_seg
    sum_rate_ci += np.mean(rate_ci_2d)

    return cs_sum, ci_file, ref_file, sum_rate_ci


################################################################################
def fits_in(in_file, param_dict):
    """
    Reading in an eventlist in .fits format to make the cross spectrum. Reads
    in a clock-corrected GTI'd event list, populates the light curves, computes
    cross spectrum per energy channel and keeps running average of the cross
    spectra.

    I take the approach: start time <= segment < end_time, to avoid double-
    counting and/or skipping events.

    """

    assert tools.power_of_two(param_dict['n_bins']), "ERROR: n_bins must be a "\
            "power of 2."

    print("Input file: %s" % in_file)

    ## Determining print iterator for segments
    if param_dict['n_bins'] == 32768:
        print_iterator = int(10)
    elif param_dict['n_bins'] < 32768:
        print_iterator = int(10)
    else:
        print_iterator = int(1)

    #######################################################
    ## Check if the FITS file exists; if so, load the data
    #######################################################

    try:
        fits_hdu = fits.open(in_file)
    except IOError:
        print("\tERROR: File does not exist: %s" % in_file)
        sys.exit()

    header = fits_hdu[0].header	 ## Header info is in ext 0, data is in ext 1
    data = fits_hdu[1].data
    fits_hdu.close()

    ###################
    ## Initializations
    ###################

    num_seg = 0

    ci = Lightcurve()
    ref = Lightcurve()

    ci.mean_rate = np.zeros(param_dict['detchans'], dtype=np.float64)
    ref.mean_rate = 0
    cs_sum = np.zeros((param_dict['n_bins'], param_dict['detchans']), dtype=np.complex128)
    ci.fullpower = np.zeros((param_dict['n_bins'], param_dict['detchans']), dtype=np.float64)
    ref.fullpower = np.zeros(param_dict['n_bins'], dtype=np.float64)
    sum_rate_ci = 0

    start_time = data.field('TIME')[0]
    final_time = data.field('TIME')[-1]
    seg_end_time = start_time + param_dict['num_seconds']

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
    print("Segments computed:")

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

            num_seg += 1

            cs_sum, ci, ref, sum_rate_ci = each_segment(time_ci, time_ref,
                energy_ci, energy_ref, param_dict, start_time,
                ci, ref, cs_sum, sum_rate_ci)

            if num_seg % print_iterator == 0:
                print("\t %d" % num_seg)
            if param_dict['test'] is True and num_seg == 1:  # For testing
                break

            start_time += param_dict['num_seconds']
            seg_end_time += param_dict['num_seconds']

        ## This next bit deals with gappy data
        elif len(time_ci) == 0 and len(time_ref) == 0:

            start_time = min(all_time_ci[0], all_time_ref[0])
            seg_end_time = start_time + param_dict['num_seconds']

        else:
            start_time += param_dict['num_seconds']
            seg_end_time += param_dict['num_seconds']

        ## End of 'if there are counts in this segment'

    ## End of while-loop

    return cs_sum, ci, ref, num_seg, sum_rate_ci


################################################################################
def get_background(bkgd_file):
    """
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
    Reads in one event list, splits into reference band and channels of
    interest (CoI), makes segments and populates them to give them length
    n_bins, computes the cross spectrum of each segment per energy channel and
    then averaged cross spectrum of all the segments per energy channel, and
    then computes the cross-correlation function (ccf) per energy channel.

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
    detchans = int(tools.get_key_val(in_file, 0, 'DETCHANS'))
    df = 1.0 / float(num_seconds)
    param_dict = {'dt': dt, 't_res': t_res, 'num_seconds': num_seconds, \
            'df': df, 'nyquist': nyquist_freq, 'n_bins': n_bins, \
            'detchans': detchans, 'test': test, 'num_seg': 0, \
            'obs_epoch': tools.obs_epoch_rxte(in_file)}

    print "\nDT = %f" % param_dict['dt']
    print "N_bins = %d" % param_dict['n_bins']
    print "Nyquist freq =", param_dict['nyquist']
    print "Filtering?", filter

    ###################################################################
    ## Reading in the background count rate from a background spectrum
    ###################################################################

    if bkgd_file:
        bkgd_rate = get_background(bkgd_file)
    else:
        bkgd_rate = np.zeros(param_dict['detchans'])
    print " "

    #################################################
    ## Reading in data, computing the cross spectrum
    #################################################

    cs_sum, ci, ref, num_seg, sum_rate_ci = fits_in(in_file, param_dict)

    param_dict['num_seg'] = num_seg

    #########################################
    ## Turning sums over segments into means
    #########################################

    ci.mean_rate /= float(param_dict['num_seg'])
    ref.mean_rate /= float(param_dict['num_seg'])
    ci.fullpower /= float(param_dict['num_seg'])
    ref.fullpower /= float(param_dict['num_seg'])
    cs_avg = cs_sum / float(param_dict['num_seg'])

    ################################################################
    ## Printing the cross spectrum to a file, for plotting/checking
    ################################################################

    cs_out = np.column_stack((fftpack.fftfreq(param_dict['n_bins'], \
                                              d=param_dict['dt']), cs_avg))
    np.savetxt('cs_avg.dat', cs_out)

    ##################################################################
    ## Subtracting the background count rate from the mean count rate
    ##################################################################

    ci.mean_rate -= bkgd_rate

    ## Need to use a background from ref. PCU for the reference band...
    ref_bkgd_rate = np.mean(bkgd_rate[2:26])
    ref.mean_rate -= ref_bkgd_rate

    ######################
    ## Making lag spectra
    ######################

    save_for_lags(out_file, in_file, param_dict, ci, ref, cs_avg)

# 	make_lags(out_file, in_file, dt, n_bins, detchans, num_seconds, num_seg, \
# 		mean_rate_ci_whole, mean_rate_ref_whole, cs_avg, mean_power_ci, \
# 		mean_power_ref)

    ##############################################
    ## Computing ccf from cs, and computing error
    ##############################################

    if filter:
        ccf_end, ccf_error = FILT_cs_to_ccf_w_err(cs_avg, param_dict,
                ci, ref, True)
    else:
        ccf_end, ccf_error = UNFILT_cs_to_ccf_w_err(cs_avg, param_dict,
                ci, ref, True)

    print "Number of segments:", param_dict['num_seg']
    print "Sum of mean rate for ci:", np.sum(ci.mean_rate)
    print "Mean rate for ref:", np.mean(ref.mean_rate)

    t = np.arange(0, param_dict['n_bins'])

    ##########
    ## Output
    ##########

    fits_out(out_file, in_file, bkgd_file, param_dict, ci.mean_rate, \
            ref.mean_rate, t, ccf_end, ccf_error, filter)


################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python ccf.py infile outfile [-b "\
        "BKGD_SPECTRUM] [-n NUM_SECONDS] [-m DT_MULT] [-t {0,1}] [-f {0,1}]",
        description="Computes the cross-correlation function of a channel of "\
        "interest with a reference band.", epilog="For optional arguments, "\
        "default values are given in brackets at end of description.")

    parser.add_argument('infile', help="The name of the GTI'd FITS event "\
        "list containing both the reference band and the channels of interest. "\
        "Assumes channels of interest = PCU 2, ref band = all other PCUs.")

    parser.add_argument('outfile', help="The name the FITS file to write the "\
        "cross-correlation function to.")

    parser.add_argument('-b', '--bkgd', required=False, dest='bkgd_file',
        help="Name of the background spectrum (in pha/fits format).")

    parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two,
        default=1, dest='num_seconds', help="Number of seconds in each Fourier"\
        " segment. Must be a power of 2, positive, integer. [1]")

    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two,
        default=1, dest='dt_mult', help="Multiple of dt (dt is from data file)"\
        " for timestep between bins. Must be a power of 2, positive, integer. "\
        "[1]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
        dest='test', help="Int flag: 0 if computing all segments, 1 if only "\
        "computing one segment for testing. [0]")

    parser.add_argument('-f', '--filter', type=int, default=0, choices={0,1},
        dest='filter', help="Int flag: 0 if NOT applying a filter in frequency"\
        "-space, 1 if applying frequency filter (around a pulsation). [0]")

    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True

    filter = False
    if args.filter == 1:
        filter = True

    main(args.infile, args.outfile, args.bkgd_file, args.num_seconds,
        args.dt_mult, test, filter)

################################################################################
