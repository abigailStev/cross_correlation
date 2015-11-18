#!/usr/bin/env python
"""
Computes the cross-correlation function of narrow energy channels of interest
with a broad energy reference band, using an RXTE .fits event list.

Use run_ccf.sh for an example.
"""
import argparse
import numpy as np
import sys
from scipy import fftpack
from datetime import datetime
import os
import subprocess
from astropy.io import fits
import tools  ## in https://github.com/abigailStev/whizzy_scripts
import ccf_lightcurves as ccf_lc

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2014-2015"


################################################################################
def fits_out(out_file, in_file, bkgd_file, meta_dict, mean_rate_ci_whole, \
    mean_rate_ref_whole, t, ccf, ccf_error, filter, lo_freq, hi_freq):
    """
    Writes the cross-correlation function to a .fits output file.

    Parameters
    ----------
    out_file : str
        Name of the FITS file to write the cross correlation function to,
        including the ".fits".

    in_file : str
        Name of the input event list, to save in the FITS header.

    bkgd_file : str
        Name of the background energy spectrum of the channels of interest,
        to save in the FITS header.

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    mean_rate_ci_whole : np.array of floats
        1-D array of the mean count rate of the channels of interest, in cts/s.
        Size = (detchans).

    mean_rate_ref_whole : float
        The mean count rate of the reference band, in cts/s.

    t : np.array of ints
        1-D array of integer time bins, size = n_bins.

    ccf : np.array of floats
        2-D array of the cross-correlation function. Size = (n_bins, detchans).

    ccf_error : np.array of floats
        2-D array of the error on the cross-correlation function.
        If filtered, size = detchans. If normal, size = (n_bins, detchans).

    filter : bool
        If True, the average cross spectrum was filtered in frequency before
        being iFFT'd to the cross-correlation function.

    lo_freq : float
        Low frequency bound of cross spectrum filter, in Hz.

    hi_freq : float
        High frequency bound of cross spectrum filter, in Hz.

    Returns
    -------
    nothing, but writes to a file.
    """

    ## Getting data into a good output structure
    chan = np.arange(0, meta_dict['detchans'])
    energy_channels = np.tile(chan, len(t))
    if filter:
        ccf_error = np.tile(ccf_error, len(t))
    else:
        ccf_error = ccf_error.flatten('C')
    time_bins = np.repeat(t, meta_dict['detchans'])
    assert len(energy_channels) == len(time_bins)

    print "\nOutput sent to: %s" % out_file

    ## Making FITS header (extension 0)
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross-correlation function")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file)
    prihdr.set('BKGD', bkgd_file)
    prihdr.set('DT', np.mean(meta_dict['dt']), "seconds")
    prihdr.set('N_BINS', meta_dict['n_bins'], "time bins per segment")
    prihdr.set('SEGMENTS', meta_dict['n_seg'], "segments in the whole light"\
        " curve")
    prihdr.set('SEC_SEG', meta_dict['n_seconds'], "seconds per segment")
    prihdr.set('EXPOSURE', meta_dict['exposure'], "seconds of data used")
    prihdr.set('DETCHANS', meta_dict['detchans'], "Number of detector energy"\
        " channels")
    prihdr.set('RATE_CI', str(mean_rate_ci_whole.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_ref_whole, "counts/second")
    prihdr.set('FILTER', str(filter))
    prihdr.set('FILTFREQ', "%f:%f" % (lo_freq, hi_freq))
    prihdr.set('ADJUST', "%d" % meta_dict['adjust_seg'])

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
def raw_to_absrms(power, mean_rate, n_bins, dt, noisy=True):
    """
    Normalizes the power spectrum to absolute rms^2 normalization.

    Parameters
    ----------
    power : np.array of floats
        The raw power at each Fourier frequency, as a 1-D or 2-D array.
        Size = (n_bins) or (n_bins, detchans).

    mean_rate : float
        The mean count rate for the light curve, in cts/s.

    n_bins : int
        Number of bins per segment of light curve.

    dt : float
        Timestep between bins in n_bins, in seconds.

    noisy : boolean
        True if there is Poisson noise in the power spectrum (i.e., from real
        data), False if there is no noise in the power spectrum (i.e.,
        simulations without Poisson noise). Default is True.

    Returns
    -------
    np.array of floats
        The noise-subtracted power spectrum in absolute rms^2 units, in the
        same size array as the input power.

    """
    if noisy:
        noise = 2.0 * mean_rate
    else:
        noise = 0.0
    # print "Power shape:", np.shape(power)
    # print "DT shape:", np.shape(dt)
    # print "Noise shape:", np.shape(noise)
    return power * (2.0 * dt / np.float(n_bins)) - noise


################################################################################
def raw_to_fracrms(power, mean_rate, n_bins, dt, noisy=True):
    """
    Normalizes the power spectrum to fractional rms^2 normalization.

    Parameters
    ----------
    power : np.array of floats
        The raw power at each Fourier frequency, as a 1-D or 2-D array.
        Size = (n_bins) or (n_bins, detchans).

    mean_rate : float
        The mean count rate for the light curve, in cts/s.

    n_bins : int
        Number of bins per segment of light curve.

    dt : float
        Timestep between bins in n_bins, in seconds.

    noisy : boolean
        True if there is Poisson noise in the power spectrum (i.e., from real
        data), False if there is no noise in the power spectrum (i.e.,
        simulations without Poisson noise). Default is True.

    Returns
    -------
    np.array of floats
        The noise-subtracted power spectrum in fractional rms^2 units, in the
        same size array as the input power.

    """
    if noisy:
        noise = 2.0 / mean_rate
    else:
        noise = 0.0

    return power * (2.0 * dt / np.float(n_bins) / (mean_rate ** 2)) - noise


################################################################################
def raw_to_leahy(power, mean_rate, n_bins, dt, noisy=True):
    """
    Normalizes the power spectrum to Leahy normalization.

    Parameters
    ----------
    power : np.array of floats
        The raw power at each Fourier frequency, as a 1-D or 2-D array.
        Size = (n_bins) or (n_bins, detchans).

    mean_rate : float
        The mean count rate for the light curve, in cts/s.

    n_bins : int
        Number of bins per segment of light curve.

    dt : float
        Timestep between bins in n_bins, in seconds.

    noisy : boolean
        True if there is Poisson noise in the power spectrum (i.e., from real
        data), False if there is no noise in the power spectrum (i.e.,
        simulations without Poisson noise). Default is True.

    Returns
    -------
    np.array of floats
        The noise-subtracted power spectrum in Leahy units, in the same size
        array as the input power.

    """
    if noisy:
        noise = 2.0
    else:
        noise = 0.0

    return power * (2.0 * dt / np.float(n_bins) / mean_rate) - noise


################################################################################
def var_and_rms(power, df):
    """
    Computes the variance and rms (root mean square) of a power spectrum.
    Assumes the negative-frequency powers have been removed. DOES NOT WORK ON
    2-D POWER ARRAYS! Not sure why.

    Parameters
    ----------
    power : np.array of floats
        1-D array (size = n_bins/2+1) of the raw power at each of the *positive*
        Fourier frequencies.

    df : float
        The step size between Fourier frequencies.

    Returns
    -------
    variance : float
        The variance of the power spectrum.

    rms : float
        The rms of the power spectrum.

    """

    # print "Shape power:", np.shape(power)
    # print "Nonzero power:", power[np.where(power<=0.0)]
    variance = np.sum(power * df, axis=0)
    # print np.shape(variance)
    # print "Variance:", variance
    # if variance > 0:
    #     rms = np.sqrt(variance)
    # else:
    #     rms = np.nan
    rms = np.where(variance >= 0, np.sqrt(variance), np.nan)
    # print "rms:", rms
    return variance, rms


################################################################################
def get_phase_err(cs_avg, power_ci, power_ref, n, m):
    """
    Computes the error on the complex phase (in radians). Power should not be
    noise-subtracted.

    Parameters
    ----------
    cs_avg : np.array of complex numbers
        2-D array of the averaged cross spectrum.

    power_ci : np.array of floats
        2-D array of the averaged power spectrum in the channels of interest.

    power_ref : np.array of floats
        1-D array of the averaged power spectrum in the reference band.

    n : int
        The number of neighbouring frequency bins averaged together to compute
        the lag.

    m : int
        The number of segments averaged over to make the cross spectra and power
        spectra.

    Returns
    -------
    phase_err : np.array of floats
        2-D array of the error on the phase of the lags.
    """
    with np.errstate(all='ignore'):
        a = power_ci * power_ref
        coherence = np.where(a != 0, np.abs(cs_avg)**2 / a, 0)
        # print np.shape(coherence)
        phase_err = np.sqrt(np.where(coherence != 0, (1 - coherence) / \
                (2 * coherence * n * m), 0))

    return phase_err


################################################################################
def phase_to_tlags(phase, f, detchans):
    """
    Converts a complex phase (in radians) to a time lag (in seconds).

    Parameters
    ----------
    phase : np.array of floats
        The complex-plane cross-spectral phase lags.

    f : np.array of floats
        Array of the Fourier frequencies for the cross-spectrum.

    detchans : int
        Number of detector energy channels for the data mode.

    Returns
    -------
    tlags : np.array of floats
        The lags converted to time (in seconds).
    """
    assert np.shape(phase) == np.shape(f), "ERROR: Phase must have same "\
            "dimensions as f."

    with np.errstate(all='ignore'):
        tlags =  np.where(f != 0, phase / (2.0 * np.pi * f), 0)

    return tlags


################################################################################
def make_lags(out_file, in_file, dt, n_bins, detchans, n_seconds,
    n_seg, mean_rate_ci_whole, mean_rate_ref_whole, cs_avg, power_ci,
    power_ref):
    """
    Computes the phase lag and time lag from the average cross spectrum, and
    writes the lag information to a .fits output file.

    Parameters
    ----------
    out_file : str
        Full path name of the output file for the cross spectrum and power
        spectra.

    in_file : str
        The input file / event list, for writing into the header.

    dt : float
        The timestep between bins in the light curve.

    n_bins : int
        The number of time bins in one Fourier segment of the data.

    detchans : int
        The number of detector energy channels for the data mode used.

    n_seconds : int
        The number of seconds in one Fourier segment of the data.

    n_seg : int
        The number of Fourier segments in all the data used.

    mean_rate_ci_whole : np.array of floats
        1-D array (size=detchans) of the mean count rate of the channels of
        interest over all the data used.

    mean_rate_ref_whole : float
        Mean count rate of the reference band over all the data used.

    cs_avg : np.array of complex numbers
        2-D array (size = n_bins x detchans) of the cross-spectrum averaged over
        segments of the light curve.

    power_ci : np.array of floats
        2-D array (size = n_bins x detchans) of the power in each energy
        channel. Includes both positive and negative Fourier frequencies.

    power_ref : np.array of floats
        1-D array (size = n_bins) of the power in the reference band. Includes
        both positive and negative Fourier frequencies.

    Returns
    -------
    nothing, but writes to the output file

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
        ## want to include 'nyq_ind'; abs is because the nyquist freq is both
        ## pos and neg, and we want it pos here. but we don't want freq=0
        ## because that gives errors.
    cs_avg = cs_avg[1:nyq_ind + 1, ]
    power_ci = power_ci[1:nyq_ind + 1, ]
    power_ref = power_ref[1:nyq_ind + 1]

    ## Getting lag and error for lag-frequency plot
    phase = -np.arctan2(cs_avg.imag, cs_avg.real) ## Negative sign is so that a
            ## positive lag is a hard energy lag?
    err_phase = get_phase_err(cs_avg, power_ci, np.resize(np.repeat(power_ref, \
            detchans), np.shape(power_ci)), 1, n_seg)
    # print np.shape(err_phase)
    f = np.resize(np.repeat(freq, detchans), (len(freq), detchans))
    tlag = phase_to_tlags(phase, f, detchans)
    err_tlag = phase_to_tlags(err_phase, f, detchans)

    ## Getting lag and error for lag-energy plot (averaging over frequencies)
    f_span_low = np.argmax(freq == low_freq)
    f_span_hi = np.argmax(freq == hi_freq)
    f_span = f_span_hi - f_span_low
    # print "Fspan low:", f_span_low
    # print "Fspan hi:", f_span_hi
    # print "F span:", f_span
    frange_freq = freq[f_span_low:f_span_hi+1]
    frange_cs = np.mean(cs_avg[f_span_low:f_span_hi+1, ], axis=0)
    frange_pow_ci = np.mean(power_ci[f_span_low:f_span_hi+1, ], axis=0)
    frange_pow_ref = np.repeat(np.mean(power_ref[f_span_low:f_span_hi+1]), \
            detchans)
    # print "Shape cs:", np.shape(frange_cs)
    e_phase = -np.arctan2(frange_cs.imag, frange_cs.real)
    e_err_phase = get_phase_err(frange_cs, frange_pow_ci, frange_pow_ref, \
            f_span, n_seg)
    # print "Shape err phase:", np.shape(e_err_phase)
    f = np.repeat(np.mean(frange_freq), detchans)
    # print "Shape phase:", np.shape(e_phase)
    # print "Shape f:", np.shape(f)
# 	exit()
    e_tlag = phase_to_tlags(e_phase, f, detchans)
    # print "Shape tlag:", np.shape(e_tlag)
    e_err_tlag = phase_to_tlags(e_err_phase, f, detchans)
    # print "Shape err tlag:", np.shape(e_err_tlag)

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
    prihdr.set('SEGMENTS', n_seg, "segments in the whole light curve")
    prihdr.set('EXPOSURE', n_seg * n_bins * dt,
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
def save_for_lags(out_file, in_file, meta_dict, mean_rate_ci, mean_rate_ref,
    cs_avg, power_ci, power_ref):
    """
    Saving header data, the cross spectrum, CoI power spectrum, and reference
    band power spectrum to a FITS file to use in the program make_lags.py to get
    cross-spectral lags.

    Parameters
    ----------
    out_file : str
        The name the FITS file to write the cross spectrum and power spectra to,
        for computing the lags.

    in_file : str
        The name of the data file (or filename containing list of data files).

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    mean_rate_ci : np.array of floats
        1-D array of the mean count rate of each of the channels of interest.
        Size = (detchans).

    mean_rate_ref : float
        Mean count rate of the reference band.

    cs_avg : np.array of complex numbers
        2-D array of the averaged cross spectrum. Size = (n_bins, detchans).

    power_ci : np.array of floats
        2-D array of the power in the channels of interest.
        Size = (n_bins, detchans).

    power_ref : np.array of floats
        1-D array of the power in the reference band. Size = (n_bins)

    Returns
    -------
    nothing, but writes to a file.

    """
    ## Getting the Fourier frequencies for the cross spectrum
    freq = fftpack.fftfreq(meta_dict['n_bins'], d=np.mean(meta_dict['dt']))
    nyq_index = meta_dict['n_bins']/2
    # print "Nyquist frequency:", freq[nyq_index]
    # print "One before nyquist:", freq[nyq_index-1]
    # print "One after nyquist:", freq[nyq_index+1]
    # print "Should be the nyquist frequency:", meta_dict['nyquist']
    # assert np.abs(freq[nyq_index]) == meta_dict['nyquist']

    ## Only keeping the parts associated with positive Fourier frequencies
    freq = np.abs(freq[0:nyq_index + 1])  ## because it slices at end-1, and we
            ## want to include 'nyq_index'; abs is because the nyquist freq is
            ## both pos and neg, and we want it pos here.
    cs_avg = cs_avg[0:nyq_index + 1, ]
    power_ci = power_ci[0:nyq_index + 1, ]
    power_ref = power_ref[0:nyq_index + 1]

    chan = np.arange(0, meta_dict['detchans'])
    energy_channels = np.tile(chan, len(freq))

    out_file = out_file.replace("cross_correlation/out_ccf", "lags/out_lags")
    out_file = out_file.replace(".", "_cs.")
    out_dir = out_file[0:out_file.rfind("/")+1]
    subprocess.call(['mkdir', '-p', out_dir])
    print "Output sent to: %s" % out_file

    ## Making FITS header (extension 0)
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross spectrum, power spectrum ci, and power spectrum "\
            "ref.")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file)
    prihdr.set('DT', np.mean(meta_dict['dt']), "seconds")
    prihdr.set('N_BINS', meta_dict['n_bins'], "time bins per segment")
    prihdr.set('SEGMENTS', meta_dict['n_seg'], "segments in the whole light"\
            " curve")
    prihdr.set('EXPOSURE', meta_dict['exposure'], "seconds of data used")
    prihdr.set('DETCHANS', meta_dict['detchans'], "Number of detector energy"\
            " channels")
    prihdr.set('RATE_CI', str(mean_rate_ci.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_ref, "counts/second")
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
            array=power_ci.flatten('C'))
    col3 = fits.Column(name='CHANNEL', unit='', format='I',
            array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3])
    tbhdu2 = fits.BinTableHDU.from_columns(cols)

    ## Making FITS table for power spectrum of reference band
    col1 = fits.Column(name='FREQUENCY', format='D', array=freq)
    col2 = fits.Column(name='POWER', unit='raw', format='D',
            array=power_ref)
    cols = fits.ColDefs([col1, col2])
    tbhdu3 = fits.BinTableHDU.from_columns(cols)

    ## If the file already exists, remove it
    assert out_file[-4:].lower() == "fits", "ERROR: Output file must have "\
            "extension '.fits'."
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

    Parameters
    ----------
    freq : np.array of floats
        1-D array (size = n_bins) of Fourier frequencies.

    power_ref : np.array of floats
        1-D array (size = n_bins) of power in the reference band.

    Returns
    -------
    pulse_freq : float
        The frequency at which there's a periodic pulsation. Assumes that the
        pulse frequency has the maximum power above 100 Hz.

    """

    ## Only searching above 100 Hz for a coherent pulse signal
    hf = np.where(freq > 100)
    hf_power = power_ref[hf]
    hf_freq = freq[hf]

    ## Assuming that the pulse frequency will have the most power
    pulse_freq = hf_freq[np.argmax(hf_power)]

    return pulse_freq


################################################################################
def filter_freq(freq_space_array, dt, n_bins, detchans, lo_freq, hi_freq,
        power_ref):
    """
    Applying a filter to the averaged cross-spectrum per energy channel (in
    frequency space). Any cross spectrum amplitudes above or below pulse_freq
    get zeroed out.

    Parameters
    ----------
    freq_space_array : np.array of complex numbers
        2-D array of the cross spectrum, in frequency space, to be filtered.
        Size = (n_bins, detchans).

    dt : float
        Timestep between bins in the light curve, in seconds.

    n_bins : int
        The number of bins in one Fourier segment of the light curve.

    detchans : int
        Number of energy channels of the detector's data mode.

    lo_freq : float
        The lower frequency bound for filtering the cross spectrum, in Hz. [0.0]

    hi_freq : float
        The upper frequency bound for filtering the cross spectrum, in Hz. [0.0]

    power_ref : np.array of floats
        1-D array of the power in the reference band. Size = (n_bins).

    Returns
    -------
    filt_freq_space_array : np.array of complex numbers
        2-D array of the cross spectrum, zeroed out at non-filtered frequencies.

    j_min : int
        Index of the minimum frequency for filtering the averaged cross spectrum
        (out of 0 to n_bins).

    j_max : int
        Index of the maximum frequency for filtering the averaged cross spectrum
        (out of 0 to n_bins).

    """
    ## Compute the Fourier frequencies
    freq = fftpack.fftfreq(n_bins, d=dt)

    ## Determine pulse frequency
    # pulse_freq = find_pulse_freq(freq, power_ref)

    # print "Pulse frequency:", pulse_freq
    # print "Index of pulse freq:", np.where(freq == pulse_freq)

    ## Get the indices of the beginning and end of the signal
    min_freq_mask = freq < lo_freq  # we want the last 'True' element
    max_freq_mask = freq > hi_freq  # we want the first 'True' element
    j_min = list(min_freq_mask).index(False)
    j_max = list(max_freq_mask).index(True)

    print "j min =", j_min
    print "j max =", j_max

    print freq[j_min]
    print freq[j_max]
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
def FILT_cs_to_ccf_w_err(cs_avg, meta_dict, countrate_ci, countrate_ref,
    power_ci, power_ref, lo_freq=0.0, hi_freq=0.0, noisy=True):
    """
    Filters the cross-spectrum in frequency space, takes the iFFT of the
    filtered cross spectrum to get the cross-correlation function, and computes
    the error on the cross-correlation function. Note that error is NOT
    independent between time bins due to the filtering! But is still independent
    between energy bins.

    Parameters
    ----------
    cs_avg : np.array of complex numbers
        2-D array of the segment-averaged cross spectrum of the channels of
        interest with the reference band. Size = (n_bins, detchans).

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    countrate_ci : np.array of floats
        1-D array of the mean count rate of each energy channel in the channels
        of interest. Size = (detchans).

    countrate_ref : float
        The mean count rate of the reference band.

    power_ci : np.array of floats
        2-D array of the raw power of the channels of interest (only positive
        Fourier frequencies). Size = (n_bins/2+1, detchans).

    power_ref: np.array of floats
        1-D array of the raw power in the reference band (only positive Fourier
        frequencies). Size = (n_bins/2+1).

    lo_freq : float
        The lower frequency bound for filtering the cross spectrum, in Hz. [0.0]

    hi_freq : float
        The upper frequency bound for filtering the cross spectrum, in Hz. [0.0]

    noisy : bool
        If True, data has Poisson noise in it that must be subtracted away. If
        False, using simulated data without Poisson noise.

    Returns
    -------
    ccf_end : np.array of floats
        2-D array of the cross-correlation function. Size = (n_bins, detchans).

    ccf_error : np.array of floats
        The error on the cross-correlation function.

    """
    ## Filter the cross spectrum in frequency
    filtered_cs_avg, j_min, j_max = filter_freq(cs_avg, meta_dict['dt'], \
            meta_dict['n_bins'], meta_dict['detchans'], lo_freq, hi_freq, \
            power_ref)

    ## Absolute rms norms of poisson noise
    noise_ci = 2.0 * countrate_ci
    noise_ref = 2.0 * countrate_ref

    ## If there's no noise in a (simulated) power spectrum, noise level = 0
    if not noisy:
        noise_ci = np.zeros(meta_dict['detchans'])
        noise_ref = 0

    noise_ref_array = np.repeat(noise_ref, meta_dict['detchans'])

    ## Extracting only the signal frequencies of the mean powers
    signal_ci_pow = np.float64(power_ci[j_min:j_max, :])
    signal_ref_pow = np.float64(power_ref[j_min:j_max])

    ## Putting powers into absolute rms2 normalization, subtracting noise
    signal_ci_pow = signal_ci_pow * (2.0 * meta_dict['dt'] / np.float(meta_dict['n_bins'])) - noise_ci
    signal_ref_pow = signal_ref_pow * (2.0 * meta_dict['dt'] / np.float(meta_dict['n_bins'])) - noise_ref

    ## Getting rms of reference band, to normalize the ccf
    ref_variance = np.sum(signal_ref_pow * meta_dict['df'])
    print "Reference band variance:", ref_variance
    rms_ref = np.sqrt(ref_variance)
    print "Frac RMS of reference band:", rms_ref / countrate_ref
    ## in frac rms units here -- should be few percent

    ## Broadcasting signal_ref_pow into same shape as signal_ci_pow
    signal_ref_pow = np.resize(np.repeat(signal_ref_pow, meta_dict['detchans']),
        np.shape(signal_ci_pow))
    assert np.shape(signal_ref_pow) == np.shape(signal_ci_pow)

    temp = (noise_ci * signal_ref_pow) + (noise_ref * signal_ci_pow) + \
            (noise_ci * noise_ref)
    cs_noise_amp = np.sqrt(np.sum(temp, axis=0) / np.float(meta_dict['n_seg']))

    temp1 = np.absolute(cs_avg[j_min:j_max, :]) * (2.0 * meta_dict['dt'] / \
            np.float(meta_dict['n_bins']))
    cs_signal_amp = np.sum(temp1, axis=0)

    ## Assuming that cs_noise_amp and cs_signal_amp are float arrays, size 64
    with np.errstate(all='ignore'):
        error_ratio = np.where(cs_signal_amp != 0, cs_noise_amp / \
                cs_signal_amp, 0)

    ## Taking the IFFT of the cross spectrum to get the CCF
    ccf_end = fftpack.ifft(filtered_cs_avg, axis=0)

    ## Dividing ccf by rms of signal in reference band
    ccf_end *= (2.0 / np.float(meta_dict['n_bins']) / rms_ref)

    ## Computing the error on the ccf
    ccf_rms_ci = np.sqrt(np.var(ccf_end, axis=0, ddof=1))
    ccf_error = ccf_rms_ci * error_ratio
    print "CCF end:", np.shape(ccf_end)
    print "CCF error:", np.shape(ccf_error)

    return ccf_end, ccf_error


################################################################################
def standard_ccf_err(cs_array, meta_dict, ref, noisy=True, absrms_var=None, \
        absrms_rms=None):
    """
    Computes the standard error on each ccf bin from the segment-to-segment
    variations. Use this for *UNFILTERED* CCFs. This error is not correlated
    between energy bins but may be correlated between time bins.

    S. Vaughan 2013, "Scientific Inference", equations 2.3 and 2.4.

    Parameters
    ----------
    cs_array : np.array of complex numbers
        Description.

    meta_dict : dictionary
        Dictionary of necessary meta-parameters for data analysis.

    ref : ccf_lc.Lightcurves object
        Description.

    Returns
    -------
    np.array of floats
        The standard error on the CCF from the segment-to-segment variations.

    """
    # print "Shape mean rate array:", np.shape(ref.mean_rate_array)
    # print "Shape power array:", np.shape(ref.power_array)

    if absrms_var == None and absrms_rms == None:
        absrms_power = np.asarray([raw_to_absrms(ref.power_array[:,i], \
                ref.mean_rate_array[i], meta_dict['n_bins'], meta_dict['dt'][i], \
                noisy) for i in range(meta_dict['n_seg'])])
        ## Note that here, the axes are weird, so it's size (n_seg, n_bins)

        # print "Shape absrms power:", np.shape(absrms_power)
        # print "Num seg:", meta_dict['n_seg']

        absrms_var, absrms_rms = var_and_rms(absrms_power.T, meta_dict['df'])

        # print absrms_rms
        mask = np.isnan(absrms_rms)
        # print "Mask:", mask
        # print "Shape absrms var:", np.shape(absrms_var)
        # print "Shape absrms rms:", np.shape(absrms_rms)

    # print "Shape cs:", np.shape(cs_array)

    ccf_array = fftpack.ifft(cs_array, axis=0).real
    # print "Shape ccf array:", np.shape(ccf_array)
    # print "Shape CCF array before nan mask:", np.shape(ccf_array)
    # ccf_array = ccf_array[:,:,~mask]
    # print "Shape CCF array after nan mask:", np.shape(ccf_array)
    # absrms_rms = absrms_rms[~mask]
    # n_nonnan = meta_dict['n_seg'] - np.count_nonzero(mask)
    # print "Number non-nan:", n_nonnan
    # print "Number nan:", np.count_nonzero(mask)
    ccf_array *= (2.0 / np.float(meta_dict['n_bins']) / absrms_rms)

    mean_ccf = np.mean(ccf_array, axis=2)
    mean_ccf_b = np.resize(np.repeat(mean_ccf, meta_dict['n_seg']), \
            np.shape(ccf_array))
    ccf_resid = ccf_array - mean_ccf_b

    sample_var = np.sum(ccf_resid**2, axis=2) / \
            np.float(meta_dict['n_seg'] - 1)  ## eq 2.3
    standard_error = np.sqrt(sample_var / \
            np.float(meta_dict['n_seg']))  ## eq 2.4

    return standard_error


################################################################################
def UNFILT_cs_to_ccf(cs_avg, meta_dict, ref, noisy, rms=None):
    """
    Takes the iFFT of the cross spectrum to get the cross-correlation function,
    and computes the error on the cross-correlation function.

    Parameters
    ----------
    cs_avg : np.array of complex numbers
        2-D array (size = n_bins x detchans) of the averaged cross-spectrum.

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    ref : ccf_lc.Lightcurves object
        The reference band light curve and power.

    noisy : bool
        True if there is noise in the light curve (i.e., using real data) --
        will subtract the Poisson noise level from the power spectrum.

    Returns
    -------
    np.array of floats
        2-D array (size = n_bins x detchans) of the cross-correlation function
        of the channels of interest with the reference band.

    """

    if len(ref.power) == meta_dict['n_bins']:
        nyq_index = meta_dict['n_bins'] / 2
        ref.pos_power = ref.power[0:nyq_index+1]

    ######################################################
    ## Take the IFFT of the cross spectrum to get the CCF
    ######################################################

    ccf = fftpack.ifft(cs_avg, axis=0).real

    if rms == None:
        ## Get the variance and rms of the reference band
        absrms = ccf_lc.NormPSD()
        fracrms = ccf_lc.NormPSD()

        absrms.power = raw_to_absrms(ref.pos_power, ref.mean_rate, \
                meta_dict['n_bins'], np.mean(meta_dict['dt']), noisy)
        fracrms.power = raw_to_fracrms(ref.pos_power, ref.mean_rate, \
                meta_dict['n_bins'], np.mean(meta_dict['dt']), noisy)

    # 	## Getting rms of reference band, to normalize the ccf
        absrms.var, absrms.rms = var_and_rms(absrms.power, \
                np.mean(meta_dict['df']))
        fracrms.var, fracrms.rms = var_and_rms(fracrms.power, \
                np.mean(meta_dict['df']))

        # print "Ref band var:", absrms.var, "(abs rms)"
        # print "Ref band rms:", absrms.rms, "(abs rms)"
        # print "Ref band var:", fracrms.var, "(frac rms)"
        # print "Ref band rms:", fracrms.rms, "(frac rms)"
    else:
        absrms = ccf_lc.NormPSD()
        absrms.rms = rms

    ## Dividing ccf by rms of signal in reference band
    ccf *= (2.0 / np.float(meta_dict['n_bins']) / absrms.rms)

    return ccf


################################################################################
def stack_reference_band(rate_ref_2d, instrument="PCA", obs_epoch=5):
    """
    Stacks the photons in the reference band from 3-20 keV to make one broad
    reference band.
    WARNING: Only tested with RXTE PCA event-mode detector channels, 0-63 incl.

    RXTE PCA, 3-20 keV:
    Epoch 1: abs channels 10 - 74
    Epoch 2: abs channels 9 - 62
    Epoch 3: abs channels 7 - 54
    Epoch 4: abs channels 6 - 46
    Epoch 5: abs channels 6 - 48

    TODO: Let user input the keV energy range for the reference band, use a
    variety of instruments, use chan.txt to convert from keV bounds to energy
    channels for the actual column selecting.

    Parameters
    ----------
    rate_ref_2d : np.array of floats
        2-D array of the reference band light curve, with
        size=(n_bins, detchans).

    instrument : str
        Name of the instrument of the reference band light curve. Currently
        only supports "PCA", as in the RXTE PCA. [PCA]

    obs_epoch : int
        If the instrument requires an observation epoch, note that here. For the
        RXTE PCA, the energy channel to keV mapping changes over the epochs. [5]

    Returns
    -------
    rate_ref : np.array of floats
        1-D array of the reference band light curve, with size=(n_bins), summed
        along the energy channel axis!

    """
    if instrument.upper() == "PCA":
        if obs_epoch == 5:
            rate_ref = np.sum(rate_ref_2d[:, 2:26], axis=1)  # EPOCH 5
            # channel 2 to 25 inclusive
            #
        elif obs_epoch == 3:
            rate_ref = np.sum(rate_ref_2d[:, 3:29], axis=1)  # EPOCH 3
            # channel 3 to 28 inclusive
        else:
            rate_ref = np.sum(rate_ref_2d, axis=1) # Summing all of it.
            raise Warning("Reference band is not being properly stacked. Need "\
                          "to put in channel information for your specific "\
                          "RXTE PCA epoch.")

    elif instrument.lower() == "NICER":
        rate_ref = np.sum(rate_ref_2d[:, 5:100], axis=1)

    else:
        rate_ref = np.sum(rate_ref_2d, axis=1) # Summing all of it.
        raise Warning("Reference band is not being properly stacked. Need "\
                          "to put in channel information for your specific "\
                          "instrument.")

    return rate_ref


################################################################################
def make_cs(rate_ci, rate_ref, meta_dict):
    """
    Generating the cross spectrum for one segment of the light curve.

    Parameters
    ----------
    rate_ci : np.array of floats
        2-D array of the channel of interest light curve,
        Size = (n_bins, detchans).

    rate_ref : np.array of floats
        1-D array of the reference band lightcurve, Size = (n_bins).

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    Returns
    -------
    cs_seg : np.array of complex numbers
        2-D array of the cross spectrum of each channel of interest with the
        reference band.

    ci_seg : ccf_lc.Lightcurves object
        The channel of interest light curve.

    ref_seg : ccf_lc.Lightcurves object
        The reference band light curve.

    """
    assert np.shape(rate_ci) == (meta_dict['n_bins'], meta_dict['detchans']),\
        "ERROR: CoI light curve has wrong dimensions. Must have size (n_bins, "\
        "detector channels)."
    assert np.shape(rate_ref) == (meta_dict['n_bins'], ), "ERROR: Reference "\
        "light curve has wrong dimensions. Must have size (n_bins, )."

    ci_seg = ccf_lc.Lightcurves(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ci')
    ref_seg = ccf_lc.Lightcurves(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ref')

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
    ci_seg.power = np.absolute(fft_data_ci) ** 2
    ref_seg.power = np.absolute(fft_data_ref) ** 2

    ## Broadcasting fft of ref into same shape as fft of ci
    fft_data_ref = np.resize(np.repeat(fft_data_ref, meta_dict['detchans']), \
        (meta_dict['n_bins'], meta_dict['detchans']))

    ## Computing the cross spectrum from the fourier transform
    cs_seg = np.multiply(fft_data_ci, np.conj(fft_data_ref))

    return cs_seg, ci_seg, ref_seg


################################################################################
def each_segment(time_ci, energy_ci, rate_ref, meta_dict,\
    start_time, end_time):
    """
    Turns the event list into a populated histogram, stacks the reference band,
    and makes the cross spectrum, per segment of light curve.

    Parameters
    ----------
    time_ci : np.array of floats
        1-D array of the photon arrival times of events in this segment for the
        channel of interest.

    time_ref : np.array of floats
        1-D array of the photon arrival times of events in this segment for the
        reference band.

    energy_ci : np.array of ints
        1-D array of the energy channels of events in this segment for the
        channel of interest.

    energy_ref : np.array of ints
        1-D array of the energy channels of events in this segment for the
        reference band.

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    start_time : float
        Starting time of the segment (front of bin), in whatever units time_ci
        and time_ref are in.

    end_time : float
        End time of the segment (back of bin), in whatever units time_ci,
        time_ref, and start_time are in.

    Returns
    -------
    cs_seg : np.array of complex numbers
        The cross spectrum of the channels of interest with the reference band
        for this segment of data. Size=(n_bins, detchans).

    ci_seg : ccf_lc.Lightcurves object
        The channels of interest light curve.

    ref_seg : ccf_lc.Lightcurves object
        The reference band light curve.

    np.mean(rate_ci_2d) : float
        The total mean count rate of the channels of interest.

    """
    assert len(time_ci) == len(energy_ci)

    ##############################################################
    ## Populate the light curves for interest and reference bands
    ##############################################################

    rate_ci_2d = tools.make_2Dlightcurve(np.asarray(time_ci),
            np.asarray(energy_ci), meta_dict['n_bins'], meta_dict['detchans'],
            start_time, end_time)

    print "CI shape:", np.shape(rate_ci_2d)
    print "ref shape:", np.shape(rate_ref)

    ## Save the reference band light curve to a text file
# 	out_file="./GX339-BQPO_ref_lc.dat"
# 	f_handle = file(out_file, 'a')
# 	np.savetxt(f_handle, rate_ref)
# 	f_handle.close()

    ###########################
    ## Make the cross spectrum
    ###########################

    cs_seg, ci_seg, ref_seg = make_cs(rate_ci_2d, rate_ref, meta_dict)

    return cs_seg, ci_seg, ref_seg, np.mean(rate_ci_2d)


################################################################################
def fits_in(in_file, ref_band_file, meta_dict, test=False):
    """
    Reading in an eventlist in .fits format to make the cross spectrum. Reads
    in a clock-corrected GTI'd event list, populates the light curves, computes
    cross spectrum per energy channel and keeps running average of the cross
    spectra.

    Assumes that the reference band times are the center of the bins, and the
    CI times are at the front of the bins.

    I take the approach: start time <= segment < end_time, to avoid double-
    counting and/or skipping events.

    Parameters
    ----------
    in_file : str
        The full path of the FITS data file being analyzed.

    ref_band_file : str
        Name of FITS optical or IR data file for reference band. This one file
        has the reference band for the whole data set. Gaps are ok.

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    test : boolean
        True if only running one segment of data for testing, False if analyzing
        the whole data file. Default=False

    Returns
    -------
   cross_spec:  np.array of floats
        3-D array of the raw cross spectrum, per segment.
        Dimensions: [n_bins, detchans, n_seg]

    ci_whole : ccf_lc.Lightcurves object
        Channel of interest for this data file.

    ref_whole : ccf_lc.Lightcurves object
        Reference band for this data file.

    n_seg : int
        Number of segments in this data file.

    dt_whole : np.array of floats
        1-D array of timestep between light curve bins for each segment. These
        will be different if artificially adjusting the QPO frequency in Fourier
        space (currently only doing that per data file, not per segment).

    df_whole : np.array of floats
        1-D array of frequency step between Fourier bins for each segment. These
        will be different if artificially adjusting the QPO frequency in Fourier
        space (currently only doing that per data file, not per segment).

    exposure : float
        The total (used) exposure of the data file.

    """

    assert tools.power_of_two(meta_dict['n_bins']), "ERROR: n_bins must be a "\
            "power of 2."
    meta_dict['obs_epoch'] = tools.obs_epoch_rxte(in_file)

    print("Input file: %s" % in_file)
    print("Ref band file : %s" % ref_band_file)

    ## Determining print iterator for segments
    if meta_dict['n_bins'] == 32768:
        print_iterator = int(10)
    elif meta_dict['n_bins'] < 32768:
        print_iterator = int(20)
    else:
        print_iterator = int(1)

    #######################################################
    ## Check if the FITS file exists; if so, load the CI data
    #######################################################

    try:
        ci_fits_hdu = fits.open(in_file)
    except IOError:
        print("\tERROR: File does not exist: %s" % in_file)
        sys.exit()

    # ci_header = ci_fits_hdu[0].header  ## Header is in ext 0
    ci_data = ci_fits_hdu[1].data  ## Data is in ext 1
    ci_fits_hdu.close()

    ###################
    ## Initializations
    ###################

    n_seg = 0
    ci_whole = ccf_lc.Lightcurves(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ci')
    ref_whole = ccf_lc.Lightcurves(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ref')
    cross_spec = np.zeros((meta_dict['n_bins'], meta_dict['detchans'], 1), \
            dtype=np.complex128)
    dt_whole = np.array([])
    df_whole = np.array([])
    exposure = 0

    ci_start_time = ci_data.field('TIME')[0]
    ci_final_time = ci_data.field('TIME')[-1]
    print("CI start: %.15f" % ci_start_time)
    print("CI final: %.15f" % ci_final_time)

    ###################################
    ## Selecting PCU for interest band
    ###################################

    PCU2_mask = ci_data.field('PCUID') == 2
    data_pcu2 = ci_data[PCU2_mask]
    all_time_ci = np.asarray(data_pcu2.field('TIME'), dtype=np.float64)
    all_energy_ci = np.asarray(data_pcu2.field('CHANNEL'), dtype=np.float64)

    ####################################
    ## Selecting PCU for reference band
    ####################################

    try:
        ref_fits_hdu = fits.open(ref_band_file)
    except IOError:
        print("\tERROR: File does not exist: %s" % ref_band_file)
        sys.exit()

    # ref_header = ref_fits_hdu[0].header	 ## Header info is in ext 0
    ref_data = ref_fits_hdu[1].data  ## Data is in ext 1
    ref_fits_hdu.close()

    ## Correcting times to make them at the front of the bin.
    all_time_ref = np.asarray(ref_data.field('TIME'), dtype=np.float64) \
            - meta_dict['dt'] / 2.0
    all_rate_ref = np.asarray(ref_data.field('RATE'), dtype=np.float64)

    ref_start_time = all_time_ref[0]
    ref_final_time = all_time_ref[-1]
    print("Ref start: %.15f" % ref_start_time)
    print("Ref final: %.15f" % ref_final_time)

    if ci_start_time > ref_start_time:
        index = np.min(np.where(all_time_ref > ci_start_time))
        start_time = all_time_ref[index]
        all_time_ref = all_time_ref[index:]
    else:
        start_time = ref_start_time

    if ci_final_time < ref_final_time:
        # I don't think this has actually been tested...
        index = np.max(np.where(all_time_ref < ci_final_time))
        print "\t", index
        final_time = all_time_ref[index]
        all_time_ref = all_time_ref[0:index]
    else:
        final_time = ref_final_time

    seg_end_time = start_time + meta_dict['n_seconds']

    print("Start: %.15f" % start_time)
    print("Final: %.15f" % final_time)

    ############################
    ## Looping through segments
    ############################
    print "Segments computed:"

    while (seg_end_time + (meta_dict['adjust_seg'] * meta_dict['dt'])) <= final_time:
        ## Adjusting segment length to artificially line up the QPOs
        seg_end_time += (meta_dict['adjust_seg'] * meta_dict['dt'])
        print "Segment end time: %.15f" % seg_end_time

        ## Get events for channels of interest
        time_ci = all_time_ci[np.where(all_time_ci < seg_end_time)]
        energy_ci = all_energy_ci[np.where(all_time_ci < seg_end_time)]

        ## Chop current segment off the rest of the list
        for_next_iteration_ci = np.where(all_time_ci >= seg_end_time)
        all_time_ci = all_time_ci[for_next_iteration_ci]
        all_energy_ci = all_energy_ci[for_next_iteration_ci]

        ## Get events for reference band
        time_ref = all_time_ref[np.where(all_time_ref < seg_end_time)]
        rate_ref = all_rate_ref[np.where(all_time_ref < seg_end_time)]

        ## Chop current segment off the rest of the list
        for_next_iteration_ref = np.where(all_time_ref >= seg_end_time)
        all_time_ref = all_time_ref[for_next_iteration_ref]
        all_rate_ref = all_rate_ref[for_next_iteration_ref]

        ###########################
        ## At the end of a segment
        ###########################

        if len(time_ci) > 0 and len(time_ref) > 0:

            cs_seg, ci_seg, ref_seg, rate_ci = each_segment(time_ci, \
                    energy_ci, rate_ref, meta_dict, start_time, \
                    seg_end_time)

            dt_seg = (seg_end_time - start_time) / float(meta_dict['n_bins'])
            df_seg = 1.0 / (meta_dict['n_bins'] * dt_seg)
            # print dt_seg, " s"
            # print df_seg, " Hz"
            dt_whole = np.append(dt_whole, dt_seg)
            df_whole = np.append(df_whole, df_seg)

            cross_spec = np.dstack((cross_spec, cs_seg))

            ci_whole.power_array = np.dstack((ci_whole.power_array, \
                    ci_seg.power))
            ci_whole.mean_rate_array = np.hstack((ci_whole.mean_rate_array, \
                    np.reshape(ci_seg.mean_rate, (meta_dict['detchans'],1)) ))
            # print ci_seg.mean_rate[1:4]
            # print np.shape(ref_whole.power_array)
            # print np.shape(ref_seg.power)
            ref_whole.power_array = np.hstack((ref_whole.power_array, \
                    np.reshape(ref_seg.power, (meta_dict['n_bins'],1)) ))
            # print np.shape(ref_whole.power_array)
            ref_whole.mean_rate_array = np.append(ref_whole.mean_rate_array, \
                    ref_seg.mean_rate)

            # print cs_seg[0:3, 0:3]
            # print cross_spec[0:3, 0:3, -1]

            ## Sums across segments -- arrays, so it adds by index
            exposure += (seg_end_time - start_time)
            n_seg += 1
            ci_whole.mean_rate += ci_seg.mean_rate
            ref_whole.mean_rate += ref_seg.mean_rate
            ci_whole.power += ci_seg.power
            ref_whole.power += ref_seg.power

            if n_seg % print_iterator == 0:
                print "\t", n_seg
            if test is True and n_seg == 1:  # For testing
                break

            start_time = seg_end_time
            seg_end_time += meta_dict['n_seconds']

        ## This next bit deals with gappy data
        elif len(time_ci) == 0 and len(time_ref) == 0:
            start_time = max(all_time_ci[0], all_time_ref[0])
            seg_end_time = start_time + meta_dict['n_seconds']

        else:
            start_time = seg_end_time
            seg_end_time += meta_dict['n_seconds']

        ## End of 'if there are counts in this segment'

    ## End of while-loop

    cross_spec = cross_spec[:,:,1:]
    ci_whole.power_array = ci_whole.power_array[:,:,1:]
    ci_whole.mean_rate_array = ci_whole.mean_rate_array[:,1:]
    ref_whole.power_array = ref_whole.power_array[:,1:]
    ref_whole.mean_rate_array = ref_whole.mean_rate_array[1:]
    # print dt_whole
    # print df_whole

    return cross_spec, ci_whole, ref_whole, n_seg, dt_whole, df_whole, \
            exposure


################################################################################
def seg_average(array):
    """
    Takes the average of an array over all the segments (i.e., the last axis).

    Parameters
    ----------
    array : np.array of numbers
        Arbitrary-dimensioned array of numbers, where the last axis is the
        segments.

    Returns
    -------
    np.mean(array, axis=-1)
        The mean of the array down the last axis, i.e. averaged over the
        segments of data.

    """
    return np.mean(array, axis=-1)

################################################################################
def get_background(bkgd_file):
    """
    Get the background count rate from a background spectrum file.

    Parameters
    ----------
    bkgd_file : string
        The full path name of the .pha file containing the background energy
        spectrum for the channels of interest, with energy channels binned in
        the same way as the data file.

    Returns
    -------
    rate : np.array of floats
        1-D array of the count rate per energy channel of the background energy
        spectrum for the channels of interest.

    """

    try:
        fits_hdu = fits.open(bkgd_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % bkgd_file
        sys.exit()

    header = fits_hdu[1].header
    data = fits_hdu[1].data
    fits_hdu.close()

    exposure = np.float(header['EXPOSURE'])
    counts = data.field('COUNTS')

    rate = counts / exposure

    return rate

################################################################################
def alltogether_means(cross_spec, ci, ref, meta_dict):
    """
    Takes the means of all the data (cross spectrum, power spectra, mean count
    rates) across the segments.

    Parameters
    ----------
    cross_spec : np.array of complex numbers
        3-dimensional array of the cross spectrum of the channels of interest
        with the reference band. Size = (n_bins, detchans, n_seg)

    ci : ccf_lc.Lightcurves object
        The channels of interest.

    ref : ccf_lc.Lightcurves object
        The reference band.

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    Returns
    -------
    avg_cross_spec : np.array of floats
        The segment-averaged cross spectrum of the channels of interest with the
        reference band. Size = (n_bins, detchans).

    cross_spec : np.array of floats
        The cross spectrum. If boot=False, has the bad segments removed. If
        boot=True, it's the same as the cross_spec that was input.
        Size = (n_bins, detchans, n_seg) where n_seg has removed bad segments.

    ci : ccf_lc.Lightcurves object
        Channels of interest, with power and mean_rate assigned. If boot=False,
        power_array and mean_rate_array have bad segments removed.

    ref : ccf_lc.Lightcurves object
        Reference band, with power and mean_rate assigned. If boot=False,
        power_array and mean_rate_array have bad segments removed.

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis. If
        boot=False, exposure and n_seg have been updated.

    """

    #########################################
    ## Turning sums over segments into means
    #########################################

    # print ci.mean_rate_array[1:4]
    # print np.shape(ci.mean_rate_array)

    ci.mean_rate = seg_average(ci.mean_rate_array)
    ci.power = seg_average(ci.power_array)
    ref.power = seg_average(ref.power_array)
    ref.mean_rate_array = np.reshape(ref.mean_rate_array, (meta_dict['n_seg']))
    ref.mean_rate = seg_average(ref.mean_rate_array)
    avg_cross_spec = seg_average(cross_spec)
    # print ci.mean_rate[1:3]
    # print bkgd_rate[1:3]

    ## Printing the cross spectrum to a file, for plotting/checking
    # cs_out = np.column_stack((fftpack.fftfreq(meta_dict['n_bins'], d=dt), \
    #         avg_cross_spec.real))
    # np.savetxt('cs_avg.dat', cs_out)

    ## Need to use a background from ref pcu for the reference band...
    # ref_total.mean_rate -= np.mean(bkgd_rate[2:26])

    return avg_cross_spec, cross_spec, ci, ref, meta_dict


################################################################################
def main(in_file, out_file, ref_band_file, bkgd_file=None, n_seconds=64,
        dt=0.0625, test=False, filtering=False, lo_freq=0.0, hi_freq=0.0,
        adjust_seg=0):
    """
    Reads in one event list, splits into reference band and channels of
    interest (CoI), makes segments and populates them to give them length
    n_bins, computes the cross spectrum of each segment per energy channel and
    then averaged cross spectrum of all the segments per energy channel, and
    then computes the cross-correlation function (ccf) per energy channel.

    Parameters
    ----------
    in_file : str
        The name of the .fits event list containing both the reference band and
        the channels of interest. Assumes channels of interest = PCU 2, ref
        band = all other PCUs.

    out_file : str
        The name the FITS file to write the cross-correlation function to.

    ref_band_file : str
        Name of FITS optical or IR data file for reference band. This one file
        has the reference band for the whole data set. Gaps are ok.

    bkgd_file : str
        Name of the background spectrum (in .pha format), with the same energy
        channel binning as the event list. [None]

    n_seconds : int
        Number of seconds in each Fourier segment. Must be a power of 2,
        positive. [64]

    dt : float
        Timestep between bins of the light curve, in seconds. Should be equal
        to the time binning of the reference band. [0.0625]

    test : bool
        If true, only computes one segment of data. If false, runs like normal.
        [False]

    filtering : bool
        If true, filters the cross spectrum in frequency space using lo_freq and
        hi_freq as boundaries. [False]

    lo_freq : float
        The lower frequency bound for filtering the cross spectrum, in Hz. Only
        in effect if filtering=True. [0.0]

    hi_freq : float
        The upper frequency bound for filtering the cross spectrum, in Hz. Only
        in effect if filtering=True. Must be hi_freq >= lo_freq (checked with
        assert statement). [0.0]

    adjust_seg : int
        How much to adjust each n_bin by to artificially shift the QPO in
        frequency. [0]

    Returns
    -------
    nothing

    """

    #####################################################
    ## Idiot checks, to ensure that our assumptions hold
    #####################################################

    assert n_seconds > 0, "ERROR: n_seconds must be a positive integer."
    assert dt > 0, "ERROR: dt must be a positive float."
    assert tools.power_of_two(n_seconds / dt), "ERROR: n_bins must be a power "\
            "of 2. Must change n_seconds."
    assert hi_freq >= lo_freq, "ERROR: Upper bound of frequency filtering must"\
            " be equal to or greater than the lower bound."

    ##################################################
    ## Initializations; 'whole' is over one data file
    ##################################################

    meta_dict = {'dt': dt,
            't_res': float(tools.get_key_val(in_file, 0, 'TIMEDEL')),
            'n_seconds': n_seconds,
            'df': 1.0 / np.float(n_seconds),
            'nyquist': 1.0 / (2.0 * dt),
            'n_bins': int(n_seconds / dt),
            'detchans': int(tools.get_key_val(in_file, 0, 'DETCHANS')),
            'filter': filtering,
            'adjust_seg': adjust_seg,
            'exposure': 0}

    print "\nDT = %f" % meta_dict['dt']
    print "N_bins = %d" % meta_dict['n_bins']
    print "Nyquist freq =", meta_dict['nyquist']
    print "Testing?", test
    print "Filtering?", meta_dict['filter']
    print "Adjust seg by: %d" % meta_dict['adjust_seg']

    ###################################################################
    ## Reading in the background count rate from a background spectrum
    ###################################################################

    if bkgd_file:
        bkgd_rate = get_background(bkgd_file)
    else:
        bkgd_rate = np.zeros(meta_dict['detchans'])
    print " "

    #################################################
    ## Reading in data, computing the cross spectrum
    #################################################

    cross_spec, ci_whole, ref_whole, n_seg, dt_whole, df_whole, exposure = \
            fits_in(in_file, ref_band_file, meta_dict, test)

    meta_dict['n_seg'] = n_seg
    meta_dict['exposure'] = exposure
    meta_dict['dt'] = dt_whole
    meta_dict['df'] = df_whole

    avg_cross_spec, cross_spec, ci_whole, ref_whole, meta_dict = \
            alltogether_means(cross_spec, ci_whole, ref_whole, meta_dict, \
            bkgd_rate, False)

    #########################################
    ## Turning sums over segments into means
    #########################################

    # avg_cross_spec = np.mean(cross_spec, axis=2)
    # ci_whole.mean_rate /= np.float(meta_dict['n_seg'])
    # ref_whole.mean_rate /= np.float(meta_dict['n_seg'])
    # ci_whole.power /= np.float(meta_dict['n_seg'])
    # ref_whole.power /= np.float(meta_dict['n_seg'])

    ################################################################
    ## Printing the cross spectrum to a file, for plotting/checking
    ################################################################

    cs_out = np.column_stack((fftpack.fftfreq(meta_dict['n_bins'], \
            d=np.mean(meta_dict['dt'])), avg_cross_spec))
    np.savetxt('cs_avg.dat', cs_out)

    ##################################################################
    ## Subtracting the background count rate from the mean count rate
    ##################################################################

    ci_whole.mean_rate -= bkgd_rate

    # ## Need to use a background from ref. PCU for the reference band...
    # ref_bkgd_rate = np.mean(bkgd_rate[2:26])
    # ref_whole.mean_rate -= ref_bkgd_rate

    ######################
    ## Making lag spectra
    ######################

    save_for_lags(out_file, in_file, meta_dict, ci_whole.mean_rate,
        ref_whole.mean_rate, avg_cross_spec, ci_whole.power, ref_whole.power)

    ##############################################
    ## Computing ccf from cs, and computing error
    ##############################################

    if filtering:
        ccf_avg, ccf_error = FILT_cs_to_ccf_w_err(avg_cross_spec, meta_dict,
            ci_whole.mean_rate, ref_whole.mean_rate, ci_whole.power,
            ref_whole.power, True, lo_freq, hi_freq)
    else:
        ccf_avg = UNFILT_cs_to_ccf(avg_cross_spec, meta_dict, ref_whole, True)
        ccf_error = standard_ccf_err(cross_spec, meta_dict, ref_whole, True)

    # print "CCF:", ccf_avg[2:7, 6]
    # print "Err:", ccf_error[2:7, 6]

    # ccf_should_be = [6.42431753, 3.42944342, 4.89985092, 3.15374201, -6.34984769]
    # err_should_be = [10.25228203, 8.59091733, 4.35719107, 4.24096029, 6.52679086]
    #
    # for (e1, e2) in zip(ccf_avg[2:7, 6], ccf_should_be):
    #     print "\t", round(e1, 8) == e2
    #
    # for (e1, e2) in zip(ccf_error[2:7, 6], err_should_be):
    #     print "\t", round(e1, 8) == e2

    print "Number of segments:", meta_dict['n_seg']
    print "Sum of mean rate for ci:", np.sum(ci_whole.mean_rate)
    print "Mean rate for ci chan 6:", ci_whole.mean_rate[6]
    print "Mean rate for ci chan 15:", ci_whole.mean_rate[15]
    print "Mean rate for ref:", np.mean(ref_whole.mean_rate)

    t = np.arange(0, meta_dict['n_bins'])

    ##########
    ## Output
    ##########

    fits_out(out_file, in_file, bkgd_file, meta_dict, ci_whole.mean_rate,
            ref_whole.mean_rate, t, ccf_avg, ccf_error, filtering, lo_freq,
             hi_freq)


################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python ccf.py infile outfile "\
            "[OPTIONAL ARGUMENTS]", description=__doc__, epilog="For optional "\
            "arguments, default values are given in brackets at the end of the"\
            " description.")

    parser.add_argument('infile', help="The name of the .fits event "\
            "list containing both the reference band and the channels of "\
            "interest. Assumes channels of interest = PCU 2, ref band = all "\
            "other PCUs.")

    parser.add_argument('outfile', help="The name the FITS file to write the "\
            "cross-correlation function to.")

    parser.add_argument('--ref', dest='ref_band_file', default=None, \
            help="Name of FITS optical or IR data file for reference band.")

    parser.add_argument('-b', '--bkgd', dest='bkgd_file', default=None,
            help="Name of the background spectrum (in .pha "\
            "format), with the same energy channel binning as the event list.")

    parser.add_argument('-n', '--n_seconds', type=tools.type_power_of_two,
            default=64, dest='n_seconds', help="Number of seconds in each "\
            "Fourier segment. Must be a power of 2, positive, integer. [64]")

    parser.add_argument('--dt', type=tools.type_positive_float, default=0.0625,
            dest='dt', help="Timestep between bins of light curve (should be "\
            "set by the reference band), in seconds. [0.0625]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
            dest='test', help="Int flag: 0 if computing all segments, 1 if "\
            "only computing one segment for testing. [0]")

    parser.add_argument('-f', '--filter', default="no", dest='filter',
            help="Filtering the cross spectrum: 'no' for QPOs, or 'lofreq:"\
            "hifreq' in Hz for coherent pulsations. [no]")

    parser.add_argument('-a', '--adjust', default=0, type=int, dest='adjust',
            help="How much to adjust each n_bin by to artificially shift the "\
            "QPO in frequency. [0]")

    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True

    filtering = False
    lo_freq = -1
    hi_freq = -1
    if args.filter.lower() != "no":
        filtering = True
        if ':' in args.filter:
            temp = args.filter.split(':')
            lo_freq = float(temp[0])
            hi_freq = float(temp[1])
        else:
            Exception("Filter keyword used incorrectly. Acceptable inputs are "\
                      "'no' or 'low:high'.")

    main(args.infile, args.outfile, args.ref_band_file,
            bkgd_file=args.bkgd_file, n_seconds=args.n_seconds,
            dt=args.dt, test=test, filtering=filtering,
            lo_freq=lo_freq, hi_freq=hi_freq, adjust_seg=args.adjust)

################################################################################
