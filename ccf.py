#!/usr/bin/env python
"""
Compute the cross-correlation function of narrow energy channels of interest
with a broad energy reference band, using an RXTE .fits event list.

Use run_ccf.sh for an example.

Be sure that 'tools.py' (from https://github.com/abigailStev/whizzy_scripts) is
downloaded, and it's directory is in your PYTHONPATH bash environment variable.

November 2015, Federico M. Vincentelli: Minor adjustements for application to IR
data

Files created
-------------
*_cs.fits : lag_spectra directory
    Output file with Fourier frequency, cross spectrum, power spectra of the
    channels of interest, and power spectrum of the reference band saved as an
    astropy table in FITS extension 1. Header info is also in FITS extension 1.

*.fits : cross_correlation directory
    Output file with cross-correlation function and error saved as an astropy
    table in FITS extension 1. Header info is also in FITS extension 1.

"""
from __future__ import print_function
import argparse
import numpy as np
import sys
from scipy import fftpack
from datetime import datetime
import os
import subprocess
from astropy.io import fits
from astropy.table import Table, Column
from astropy.modeling import powerlaws, models, fitting
from astropy.modeling.models import custom_model
import matplotlib.pyplot as plt

## These are things I've written.
## Their directories are in my PYTHONPATH bash environment variable.
import tools  ## in whizzy_scripts
import ccf_lightcurves as ccf_lc  ## in cross_correlation

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2014-2017"

################################################################################
def type_positive_float(num):
    """
    Check if an input is a positive float, as an argparse type.

    Parameters
    ----------
    num : int, long, float, or double
        The number in question.

    Returns
    -------
    n : float
        The input number, if it's positive

    Raises
    ------
    ArgumentTypeError if num isn't a positive float.

    """
    try:
        n = float(num)
    except ValueError or TypeError:
        message = "Input is not a positive float."
        raise argparse.ArgumentTypeError(message)

    if n >= 0:
        return n
    else:
        message = "%d is not a positive number." % n
        raise argparse.ArgumentTypeError(message)

################################################################################
def power_of_two(num):
    """
    Check if a positive integer is a power of 2 (1 <= num < 2147483648).

    Parameters
    ----------
    num : int
        The number in question.

    Returns
    -------
    bool
        True if 'num' is a power of two, False if 'num' is not.

    Raises
    ------
    assert error if number is below zero.

    """
    n = int(num)
    x = 2
    assert n > 0, "ERROR: Number must be positive."

    if n == 1:
        return True
    else:
        while x < n and x < 2147483648:
            x *= 2
        return n == x

################################################################################
def type_positive_int(num):
    """
    Check if an input is a positive integer, as an argparse type.

    Parameters
    ----------
    num : int, long, float, or double
        The number in question.

    Returns
    -------
    n : int
        The input number, if it's a positive integer

    Raises
    ------
    ArgumentTypeError if n isn't a real number or a positive integer.

    """
    try:
        n = int(num)
    except ValueError or TypeError:
        message = "Input is not a positive integer."
        raise argparse.ArgumentTypeError(message)

    if n >= 0:
        return n
    else:
        message = "%d is not a positive integer." % n
        raise argparse.ArgumentTypeError(message)

################################################################################
def type_power_of_two(num):
    """
    Check if an input is a power of 2 (1 <= num < 2147483648), as an argparse
    type.

    Parameters
    ----------
    num : int
        The number in question.

    Returns
    -------
    n : int
        The number in question, if it's a power of two

    Raises
    ------
    ArgumentTypeError if n isn't a power of two.

    """
    n = int(num)
    x = 2
    assert n > 0

    if n == 1:
        return n
    else:
        while x <= n and x < 2147483648:
            if n == x:
                return n
            x *= 2

    message = "%d is not a power of two." % n
    raise argparse.ArgumentTypeError(message)


################################################################################
def flx2xsp_out(file_root, ref, freq, n_seg, n_bins):
    """
    Writes the reference band power spectrum to a text file, to be read in by
    the FTOOLS command FLX2SXP, to be made into a .pha spectrum for fitting in
    XSPEC.

    Parameters
    ----------
    file_root : str
        The basic file output name. Should be the ccf file output name.

    ref : ccf_lc.Lightcurve()
        The reference band. Only needs to have pos_power assigned, the power
        at the positive Fourier frequencies (including Nyquist).

    freq : np.array of floats
        The positive Fourier frequencies, including Nyquist. Size n_bins/2+1.

    n_seg : int
        The number of segments averaged together to make the power spectrum.

    n_bins : int
        The number of frequency bins in one segment of data.

    Returns
    -------
    nothing, but writes to a file.

    """
    dir = os.path.dirname(file_root)
    base = os.path.basename(file_root)[:-5]
    out_file = dir + "/" + base + "_powref_flx2xsp.txt"
    print("Table for FLX2XSP:", out_file)

    # if len(freq) == n_bins:
    #     freq = np.abs(freq[0:n_bins/2+1])

    err_pow = ref.pos_power / np.sqrt(float(n_seg))

    delta_freq = freq[1] - freq[0]

    out_table = Table()
    out_table.add_column(Column(data=freq[1:]-delta_freq/2., name='FREQ_MIN'))
    out_table.add_column(Column(data=freq[1:]+delta_freq/2., name='FREQ_MAX'))
    out_table.add_column(Column(data=ref.pos_power[1:]*freq[1:]*delta_freq,
                                name='POW*FREQ*DF'))
    out_table.add_column(Column(data=err_pow[1:]*freq[1:]*delta_freq,
                                name='ERR_POW*FREQ*DF'))

    out_table.write(out_file, format='ascii.no_header', delimiter="\t")


################################################################################
def ci_flx2xsp_out(file_root, ci, ref, freq, n_seg):
    """
    For writing all the power spectra out to a text file, to be read in by
    the FTOOLS command FLX2XSP, to be made into a .pha spectrum for fitting in
    XSPEC. This is done to make energy-channel-dependent cross spectrum filters.

    Parameters
    ----------
    file_root : str
        The basic file output name. Should be the ccf file output name.

    ci : ccf_lc.Lightcurve()
        The channels of interest. Only needs to have pos_power assigned, the
        power at the positive Fourier frequencies (including Nyquist).

    ref : ccf_lc.Lightcurve()
        The reference band. Only needs to have pos_power assigned, the power
        at the positive Fourier frequencies (including Nyquist).

    freq : np.array of floats
        The positive Fourier frequencies, including Nyquist. Size n_bins/2+1.

    n_seg : int
        The number of segments averaged together to make the power spectrum.

    Returns
    -------
    nothing, but writes to a file.

    """
    dir = os.path.dirname(file_root)
    base = os.path.basename(file_root)[:-5]
    out_file = dir + "/" + base + "_allps.txt"
    print("Table for FLX2XSP:", out_file)

    delta_freq = freq[1] - freq[0]

    out_table = Table()
    out_table.add_column(Column(data=freq[1:]-delta_freq/2., name='FREQ_MIN'))
    out_table.add_column(Column(data=freq[1:]+delta_freq/2., name='FREQ_MAX'))

    out_table.add_column(Column(data=ref.pos_power[1:]*freq[1:]*delta_freq,
                                name='POW REF'))
    err_pow = ref.pos_power / np.sqrt(float(n_seg))

    out_table.add_column(Column(data=err_pow[1:]*freq[1:]*delta_freq,
                                name='ERR REF'))
    for i in range(2, 10):
        err_pow = ci.pos_power[:,i] / np.sqrt(float(n_seg))
        out_table.add_column(Column(data=ci.pos_power[1:,i] * freq[1:] * delta_freq,
                                name='POW CHAN %d' % i))
        out_table.add_column(Column(data=err_pow[1:] * freq[1:] * delta_freq,
                                name='ERR CHAN %d' % i))
    for i in range(11, 27):
        err_pow = ci.pos_power[:, i] / np.sqrt(float(n_seg))
        out_table.add_column(
            Column(data=ci.pos_power[1:, i] * freq[1:] * delta_freq,
                   name='POW CHAN %d' % i))
        out_table.add_column(Column(data=err_pow[1:] * freq[1:] * delta_freq,
                                    name='ERR CHAN %d' % i))

    out_table.write(out_file, format='ascii.no_header', delimiter="\t")


################################################################################
def find_nearest(array, value):
    """
    Thanks StackOverflow!

    Parameters
    ----------
    array : np.array of ints or floats
        1-D array of numbers to search through. Should already be sorted from
        low values to high values.

    value : int or float
        The value you want to find the closest to in the array.

    Returns
    -------
    array[idx] : int or float
        The array value that is closest to the input value.

    idx : int
        The index of the array of the closest value.
    """
    idx = np.searchsorted(array, value, side="left")
    if idx == len(array) or np.fabs(value - array[idx-1]) < \
            np.fabs(value - array[idx]):
        return array[idx-1], idx-1
    else:
        return array[idx], idx


@custom_model
def bbn(x, amp=1.0, x_0=0.9, fwhm=3.0):
    numerator = fwhm / (np.pi * 2.0)
    denominator = (x - x_0) ** 2 + (1.0 / 2.0 * fwhm) ** 2
    L = (numerator / denominator) * amp * x
    return L


@custom_model
def qpo_and_harmonic(x, amp_f=1., x_0_f=3., fwhm_f=0.5, amp_h=1., fwhm_h=0.8):
    numerator1 = fwhm_f / (np.pi * 2.0)
    denominator1 = (x - x_0_f) ** 2 + (1.0 / 2.0 * fwhm_f) ** 2
    L1 = (numerator1 / denominator1) * amp_f * x
    numerator2 = fwhm_h / (np.pi * 2.0)
    denominator2 = (x - 2 * x_0_f) ** 2 + (1.0 / 2.0 * fwhm_h) ** 2
    L2 = (numerator2 / denominator2) * amp_h * x
    return L1 + L2


@custom_model
def qpo_fundamental_only(x, amp_f=1., x_0_f=3., fwhm_f=0.5):
    numerator = fwhm_f / (np.pi * 2.0)
    denominator = (x - x_0_f) ** 2 + (1.0 / 2.0 * fwhm_f) ** 2
    L = (numerator / denominator) * amp_f * x
    return L


################################################################################
def fits_out(out_file, in_file, bkgd_file, meta_dict, mean_rate_ci,
        mean_rate_ref, rms_ref, ccf, ccf_error, file_desc):
    """
    Write the cross-correlation function to a .fits output file.

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

    mean_rate_ci : np.array of floats
        1-D array of the mean count rate of the channels of interest with
        background count rate per CI subtracted (if applicable), in cts/s.
        Size = (detchans).

    mean_rate_ref : float
        The mean count rate of the reference band, in cts/s.

    ccf : np.array of floats
        2-D array of the cross-correlation function. Size = (n_bins, detchans).

    ccf_error : np.array of floats
        2-D array of the error on the cross-correlation function.
        If filtered, size = detchans. If normal, size = (n_bins, detchans).

    file_desc : str
        Short description of the output file.

    Returns
    -------
    Nothing, but writes to out_file.fits.
    """
    print(out_file)
    ## Check that the output file name has FITS file extension
    assert out_file[-4:].lower() == "fits", "ERROR: Output file must have "\
            "extension '.fits'."

    # print("\nOutput sent to: %s" % out_file)

    out_table = Table()
    out_table.add_column(Column(data=ccf, name='CCF'))
    out_table.add_column(Column(data=ccf_error, name='ERROR'))

    out_table.meta['TYPE'] = file_desc
    out_table.meta['DATE'] = str(datetime.now())
    out_table.meta['EVTLIST'] = in_file
    out_table.meta['BKGD'] = bkgd_file
    out_table.meta['DT'] = np.mean(meta_dict['dt'])
    out_table.meta['N_BINS'] = meta_dict['n_bins']
    out_table.meta['SEGMENTS'] = meta_dict['n_seg']
    out_table.meta['SEC_SEG'] = meta_dict['n_seconds']
    out_table.meta['EXPOSURE'] = meta_dict['exposure']
    out_table.meta['DETCHANS'] = meta_dict['detchans']
    out_table.meta['RATE_CI'] = str(mean_rate_ci.tolist())
    out_table.meta['RATE_REF'] = mean_rate_ref
    out_table.meta['RMS_REF'] = rms_ref
    out_table.meta['NYQUIST'] = meta_dict['nyquist']
    out_table.meta['DF'] = np.mean(meta_dict['df'])
    out_table.meta['FILTER'] = str(meta_dict['filter'])
    out_table.meta['ADJUST'] = "%s" % str(meta_dict['adjust_seg'])
    # if meta_dict['ref_file']:
    #     out_table.meta['REF_FILE'] = meta_dict['ref_file']

    out_table.write(out_file, overwrite=True)


################################################################################
def rebin(freq, power, err_power, rebin_const):
    """
    Re-bin the power spectrum in frequency space by some re-binning constant
    (rebin_const > 1).

    Parameters
    ----------
    freq : np.array of floats
        1-D array of the Fourier frequencies.

    power : np.array of floats
        1-D array of the power at each Fourier frequency, with any/arbitrary
        normalization.

    err_power : np.array of floats
        1-D array of the error on the power at each Fourier frequency, with the
        same normalization as the power.

    rebin_const : float
        The constant by which the data were geometrically re-binned.

    Returns
    -------
    rb_freq : np.array of floats
        1-D array of the re-binned Fourier frequencies.

    rb_power : np.array of floats
        1-D array of the power at the re-binned Fourier frequencies, with the
        same normalization as the input power array.

    rb_err : np.array of floats
        1-D array of the error on the power at the re-binned Fourier
        frequencies, with the same normalization as the input error on power.

    freq_min : np.array of floats
        1-D array of the lower bounds of each re-binned frequency bin.

    freq_max : np.array of floats
        1-D array of the upper bounds of each re-binned frequency bin.

    """
    assert rebin_const >= 1.0

    ## Initialize variables
    rb_power = np.asarray([])  # List of re-binned power
    rb_freq = np.asarray([])   # List of re-binned frequencies
    rb_err = np.asarray([])	   # List of error in re-binned power
    real_index = 1.0		   # The unrounded next index in power
    int_index = 1			   # The int of real_index, added to current_m every
                               #  iteration
    current_m = 1			   # Current index in power
    prev_m = 0				   # Previous index m
    bin_power = 0.0			   # The power of the current re-binned bin
    bin_freq = 0.0			   # The frequency of the current re-binned bin
    err_bin_power2 = 0.0	   # The error squared on 'bin_power'
    bin_range = 0.0			   # The range of un-binned bins covered by this
                               #  re-binned bin
    freq_min = np.asarray([])
    freq_max = np.asarray([])

    ## Loop through the length of the array power, new bin by new bin, to
    ## compute the average power and frequency of that new geometric bin.
    ## Equations for frequency, power, and error are from A. Ingram's PhD thesis
    while current_m < len(power):

        ## Determine the range of indices this specific geometric bin covers
        bin_range = np.absolute(current_m - prev_m)
        ## Want mean power of data points contained within one geometric bin
        bin_power = np.mean(power[prev_m:current_m])
        ## Compute error in bin -- equation from Adam Ingram's thesis
        err_bin_power2 = np.sqrt(np.sum(err_power[prev_m:current_m] ** 2)) / \
            float(bin_range)

        ## Compute the mean frequency of a geometric bin
        bin_freq = np.mean(freq[prev_m:current_m])

        ## Append values to arrays
        rb_power = np.append(rb_power, bin_power)
        rb_freq = np.append(rb_freq, bin_freq)
        rb_err = np.append(rb_err, err_bin_power2)
        freq_min = np.append(freq_min, freq[prev_m])
        freq_max = np.append(freq_max, freq[current_m])

        ## Increment for the next iteration of the loop
        ## Since the for-loop goes from prev_m to current_m-1 (since that's how
        ## the range function and array slicing works) it's ok that we set
        ## prev_m = current_m here for the next round. This will not cause any
        ## double-counting bins or skipping bins.
        prev_m = current_m
        real_index *= rebin_const
        int_index = int(round(real_index))
        current_m += int_index
        bin_range = None
        bin_freq = None
        bin_power = None
        err_bin_power2 = None

    return rb_freq, rb_power, rb_err, freq_min, freq_max


################################################################################
def raw_to_absrms(power, mean_rate, n_bins, dt, noisy=True):
    """
    Normalize the power spectrum to absolute rms^2 normalization.

    TODO: cite paper.

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
    # print("Power shape:", np.shape(power))
    # print("DT shape:", np.shape(dt))
    # print("Noise shape:", np.shape(noise))
    return power * (2.0 * dt / np.float(n_bins)) - noise


################################################################################
def raw_to_fracrms(power, mean_rate, n_bins, dt, noisy=True):
    """
    Normalize the power spectrum to fractional rms^2 normalization.

    TODO: cite paper

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
    Normalize the power spectrum to Leahy normalization.

    TODO: cite paper

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

    TODO: cite textbook or paper.

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

    # print("Shape power:", np.shape(power))
    # print("Nonzero power:", power[np.where(power<=0.0)])
    variance = np.sum(power * df, axis=0)
    # print(np.shape(variance))
    # print("Variance:", variance)
    # if variance > 0:
    #     rms = np.sqrt(variance)
    # else:
    #     rms = np.nan
    rms = np.where(variance >= 0, np.sqrt(variance), np.nan)
    # print("rms:", rms)
    return variance, rms


################################################################################
def stack_reference_band(rate_ref_2d, instrument="PCA", obs_epoch=5):
    """
    Stack the photons in the reference band from 3-20 keV to make one broad
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
        elif obs_epoch == 3:
            rate_ref = np.sum(rate_ref_2d[:, 3:29], axis=1)  # EPOCH 3
            # channel 3 to 28 inclusive
        elif obs_epoch == 1:
            rate_ref = np.sum(rate_ref_2d[:, 6:36], axis=1)  # EPOCH 1
            # channel 6 to 35 inclusive
        elif obs_epoch == 0:
            rate_ref = np.sum(rate_ref_2d, axis=1)  # SUMMING ALL
        else:
            # rate_ref = np.sum(rate_ref_2d, axis=1)  # Summing all of it.
            rate_ref = np.sum(rate_ref_2d[:, 3:30], axis=1)  # EPOCH 5
            print("Reference band is not being properly stacked. Need "\
                          "to put in channel information for your specific "\
                          "RXTE PCA epoch.")

    elif instrument.lower() == "NICER":
        rate_ref = np.sum(rate_ref_2d[:, 5:100], axis=1)

    else:
        rate_ref = np.sum(rate_ref_2d, axis=1) # Summing all of it.
        print("Reference band is not being properly stacked. Need "\
                          "to put in channel information for your specific "\
                          "instrument.")

    return rate_ref


################################################################################
def optimal_filt(meta_dict, dt):
    """
    Reads in from a mod .xcm file, for fitting multiple power spectra in XSPEC.
    Expects powerlaw + lorf + lorf + lorf + lorf + lorf
    .xcm file name is hardcoded. Saves parameters into a ccf_lc.Filter object
    and makes filters.

    :param meta_dict:
    :param dt:
    :return:
    """
    pow_filt_file = "/Users/abigail/Dropbox/Research/power_spectra/out_ps/GX339-4HzCQPO/allps_rb_mod_linetied_sigmatied_normfree.xcm"
    # pow_filt_file = "/Users/abigailstevens/Dropbox/Research/power_spectra/out_ps/GX339-4HzCQPO/allps_rb_mod_linetied_sigmatied_normfree.xcm"

    ## The first one is the reference band,
    ## the rest are channels 2-25 inclusive, without channel 10 (since it's
    ## zeros)
    f = open(pow_filt_file, 'r')
    f.seek(140)
    big_filter = []
    i = 1  # lines to read in for one spec;
            # a counter so we know when to stop and make the filter
    j = 0  # parameters in one spec; index in 'pars' array
    k = 0
    pars = np.zeros(14)
    for line in f:
        element0 = line.split()[0]
        element1 = line.split()[1]
        # print("i=%d, j=%d" %(i,j))

        if element0 != '=':
            # print(element0)
            pars[j] = element0
            j += 1
        else:
            # print(element1)
            if i == 10 or i == 15 or i == 16: # the three values computed by p12 and p13
                pass
            else:
                j += 1
        # print(pars)
        if i == 17: # There are 17, but we start at 0
            this_filt = ccf_lc.Filter(pars, k, n_bins=meta_dict['n_bins'], dt=dt)
            big_filter.append(this_filt)
            j = 0
            i = 1
            # print(len(pars))
            # print(pars[10:14])
            # print(this_filt)
            k += 1
        else:
            i += 1

        # For testing: to stop after 2 bands (ref and one ci) have been read in
        # if k == 2:
        #     exit()

    return big_filter

################################################################################
def TypeB_optimal_filt(meta_dict, dt):
    """
    Reads in from a mod .xcm file, for fitting multiple power spectra in XSPEC.
    Expects powerlaw + lorf + lorf + lorf
    .xcm file name is hardcoded. Saves parameters into a ccf_lc.TypeBFilter object
    and makes filters.

    :param meta_dict:
    :param dt:
    :return:
    """
    pow_filt_file = "/Users/abigail/Dropbox/Research/power_spectra/out_ps/GX339-BQPO/GX339-BQPO_170221_t64_64sec_adj_rb_allps.xcm"
    # pow_filt_file = "/Users/abigailstevens/Dropbox/Research/power_spectra/out_ps/GX339-4HzCQPO/allps_rb_mod_linetied_sigmatied_normfree.xcm"

    ## The first one is the reference band,
    ## the rest are channels 2-25 inclusive, without channel 10 (since it's
    ## zeros)
    f = open(pow_filt_file, 'r')
    f.seek(138)
    big_filter = []
    i = 1  # lines to read in for one spec;
            # a counter so we know when to stop and make the filter
    j = 0  # parameters in one spec; index in 'pars' array
    k = 0
    pars = np.zeros(9)
    for line in f:
        element0 = line.split()[0]
        element1 = line.split()[1]
        # print(element0, element1)
        # print("i=%d, j=%d" %(i,j))

        if element0 != '=':
            # print(element0)
            pars[j] = element0
            j += 1
        else:
            # print(element1)
            if i == 9 or i == 10:
                pass
            else:
                j += 1
        # print(pars)
        if i == 11: #
            this_filt = ccf_lc.TypeBFilter(pars, k, n_bins=meta_dict['n_bins'], dt=dt)
            big_filter.append(this_filt)
            j = 0
            i = 1
            # print(len(pars))
            # print(pars[10:14])
            # print(this_filt)
            k += 1
        else:
            i += 1

        ## For testing: to stop after 2 bands (ref and one ci) have been read in
        # if k == 2:
        #     exit()

    return big_filter


################################################################################
def make_cs(rate_ci, rate_ref, meta_dict, dt_seg):
    """
    Generate the cross spectrum for one segment of the light curve.

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

    ci_seg : ccf_lc.Lightcurve object
        The channel of interest light curve.

    ref_seg : ccf_lc.Lightcurve object
        The reference band light curve.

    """
    assert np.shape(rate_ci) == (meta_dict['n_bins'], meta_dict['detchans']),\
        "ERROR: CoI light curve has wrong dimensions. Must have size (n_bins, "\
        "detector channels)."
    assert np.shape(rate_ref) == (meta_dict['n_bins'], ), "ERROR: Reference "\
        "light curve has wrong dimensions. Must have size (n_bins, )."

    ci_seg = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ci')
    ref_seg = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ref')
    cs_seg = ccf_lc.Cross(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'])

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

    ## Getting the optimal filter
    optimal_filter = optimal_filt(meta_dict, dt_seg)
    optimal_filter = TypeB_optimal_filt(meta_dict, dt_seg)

    ## Applying the optimal filter to the reference band
    fft_data_ref_fund = fft_data_ref * optimal_filter[0].fund_filt
    fft_data_ref_harm = fft_data_ref * optimal_filter[0].harm_filt
    fft_data_ref_both = fft_data_ref * optimal_filter[0].both_filt
    ref_seg.filt_qpo_power = np.absolute(fft_data_ref_both) ** 2
    ref_seg.filt_qpo_power = ref_seg.filt_qpo_power[0:meta_dict['n_bins']/2+1]

    ## Applying the optimal filter to the channels of interest
    ## I have filters for channels 2-25 inclusive (except 11), which start
    ## at filter index 1 (since 0 is the reference band)
    ci_fund_filt = np.ones(np.shape(fft_data_ci))
    # print(ci_fund_filt[1,1])
    for i in range(1,9):
        # print(optimal_filter[i].continuum)
        ci_fund_filt[:,i+1] = optimal_filter[i].fund_filt
    for i in range(9,25):
        ci_fund_filt[:,i+2] = optimal_filter[i].fund_filt
    ci_harm_filt = np.ones(np.shape(fft_data_ci))
    for i in range(1, 9):
        ci_harm_filt[:, i + 1] = optimal_filter[i].harm_filt
    for i in range(9, 25):
        ci_harm_filt[:, i + 2] = optimal_filter[i].harm_filt

    # np.savetxt("ref_fund_filt.txt", optimal_filter[0].fund_filt)
    # np.savetxt("ref_harm_filt.txt", optimal_filter[0].harm_filt)
    # print(optimal_filter[0].continuum)
    # np.savetxt("ref_continuum.txt", optimal_filter[0].continuum)
    # np.savetxt("ref_fund.txt", optimal_filter[0].fund)
    # np.savetxt("ci_chan2_continuum.txt", optimal_filter[1].continuum)
    # np.savetxt("ci_chan2_fund.txt", optimal_filter[1].fund)
    # np.savetxt("ci_chan22_continuum.txt", optimal_filter[-3].continuum)
    # np.savetxt("ci_chan22_fund.txt", optimal_filter[-3].fund)
    # np.savetxt("ci_fund_filt.txt", ci_fund_filt)
    # np.savetxt("ci_harm_filt.txt", ci_harm_filt)
    # np.savetxt("pos_freqs.txt", optimal_filter[0].pos_freq)
    # print("Ref fund:", optimal_filter[0].fund[240])
    # print("Ref con:", optimal_filter[0].continuum[240])
    # print("Ref ratio:", optimal_filter[0].fund[240]/optimal_filter[0].continuum[240])
    # print("2 fund:", optimal_filter[2].fund[240])
    # print("2 con:", optimal_filter[2].continuum[240])
    # print("2 ratio:", optimal_filter[2].fund[240]/optimal_filter[2].continuum[240])
    # print("6 fund:", optimal_filter[6].fund[240])
    # print("6 con:", optimal_filter[6].continuum[240])
    # print("6 ratio:", optimal_filter[6].fund[240]/optimal_filter[6].continuum[240])

    # print("Before filter:")
    # print("CI2:", fft_data_ci[240,2], phase_angle(fft_data_ci[240,2]))
    # print("Ref:", fft_data_ref[240], phase_angle(fft_data_ref[240]))

    fft_data_ci_fund = fft_data_ci * ci_fund_filt
    fft_data_ci_harm = fft_data_ci * ci_harm_filt

    # print("After filter:")
    # print("CI2:", fft_data_ci_fund[240,2], phase_angle(fft_data_ci_fund[240,2]))
    # print("Ref:", fft_data_ref_fund[240], phase_angle(fft_data_ref_fund[240]))

    # lightcurve1 = fftpack.ifft(fft_data_ci_fund, axis=0)
    # print(np.shape(lightcurve1))
    # np.savetxt("ci_filtered_lightcurve.txt", lightcurve1)
    # exit()

    ## Broadcasting fft of ref into same shape as fft of ci
    fft_data_ref_fund = np.resize(np.repeat(fft_data_ref_fund, \
            meta_dict['detchans']), (meta_dict['n_bins'], \
                                     meta_dict['detchans']))
    fft_data_ref_harm = np.resize(np.repeat(fft_data_ref_harm, \
            meta_dict['detchans']), (meta_dict['n_bins'], \
                                     meta_dict['detchans']))

    ## Computing the cross spectrum from the fourier transform
    ## For both the fundamental filtered part and the harmonic filtered part
    cs_seg.fund = np.multiply(fft_data_ci_fund, np.conj(fft_data_ref_fund))
    cs_seg.harm = np.multiply(fft_data_ci_harm, np.conj(fft_data_ref_harm))

    ## Computing the cross spectrum and power spectra for frequency lags
    f_lag_seg = ccf_lc.LagFreq(n_bins=meta_dict['n_bins'])
    lc_1 = np.sum(rate_sub_mean_ci[:, 2:8], axis=-1)
    lc_2 = np.sum(rate_sub_mean_ci[:, 13:22], axis=-1)
    fft_1 = fftpack.fft(lc_1)
    fft_2 = fftpack.fft(lc_2)
    f_lag_seg.cross = np.multiply(fft_1, np.conj(fft_2))
    f_lag_seg.pow1 = np.absolute(fft_1) ** 2
    f_lag_seg.pow2 = np.absolute(fft_2) ** 2

    return cs_seg, ci_seg, ref_seg, f_lag_seg


################################################################################
def fits_in(in_file, meta_dict, test=False):
    """
    Read in an eventlist in .fits format to make the cross spectrum. Read
    in a clock-corrected GTI'd event list, populate the light curves, compute
    cross spectrum per energy channel and keep running sum (later to be an
    average) of the cross spectra.

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

    ci_whole : ccf_lc.Lightcurve object
        Channel of interest for this data file.

    ref_whole : ccf_lc.Lightcurve object
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

    assert power_of_two(meta_dict['n_bins']), "ERROR: n_bins must be a "\
            "power of 2."
    meta_dict['obs_epoch'] = tools.obs_epoch_rxte(in_file)

    print("Input file: %s" % in_file)

    ## Determining print iterator for segments
    if meta_dict['n_bins'] == 32768:
        print_iterator = int(10)
    elif meta_dict['n_bins'] < 32768:
        print_iterator = int(10)   #TODO: change back to 20
    else:
        print_iterator = int(1)

    #######################################################
    ## Check if the FITS file exists; if so, load the data
    #######################################################

    time = np.asarray([])
    channel = np.asarray([])
    pcuid = np.asarray([])

    ## Reading an event list from an astropy table FITS file
    # try:
    #     data_table = Table.read(in_file)
    #     time = data_table['TIME']
    #     channel = data_table['CHANNEL']
    #     pcuid = data_table['PCUID']
    # except IOError:
    #     print("\tERROR: File does not exist: %s" % in_file)
    #     exit()

    ## Reading an event list from a normal FITS table
    try:
        fits_hdu = fits.open(in_file)
        time = fits_hdu[1].data.field('TIME')  ## Data is in ext 1
        channel = fits_hdu[1].data.field('CHANNEL')
        pcuid = fits_hdu[1].data.field('PCUID')
        fits_hdu.close()
    except IOError:
        print("\tERROR: File does not exist: %s" % in_file)
        exit()

    # try:
    #     fits_hdu = fits.open(in_file)
    # except IOError:
    #     print("\tERROR: File does not exist: %s" % in_file)
    #     sys.exit()
    #
    # data = fits_hdu[1].data
    # fits_hdu.close()

    ###################
    ## Initializations
    ###################

    n_seg = 0
    ci_whole = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ci')
    ref_whole = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ref')
    cs_whole = ccf_lc.Cross(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'])
    f_lag_whole = ccf_lc.LagFreq(meta_dict['n_bins'])
    dt_whole = np.array([])
    df_whole = np.array([])
    exposure = 0
    start_time = time[0]
    final_time = time[-1]

    ###################################
    ## Selecting PCU for interest band
    ###################################

    PCU2_mask = pcuid == 2
    time_pcu2 = time[PCU2_mask]
    chan_pcu2 = channel[PCU2_mask]

    all_time_ci = np.asarray(time_pcu2, dtype=np.float64)
    all_energy_ci = np.asarray(chan_pcu2, dtype=np.float64)

    ######################################
    ## Getting reference band light curve
    ######################################

    if not meta_dict['ref_file']:
        ## Determining the next-most-prevalent pcu
        # pcus_on, occurrences = np.unique(data.field('PCUID'), return_counts=True)
        ## Getting rid of the PCU 2 element (since we don't want that as the ref)
        # pcu2 = np.where(pcus_on == 2)
        # pcus_on = np.delete(pcus_on, pcu2)
        # occurrences = np.delete(occurrences, pcu2)
        ## Next-most pcu is:
        # most_pcu = np.argmax(occurrences)
        # ref_pcu = pcus_on[most_pcu]
        # print("Ref PCU =", ref_pcu)
    	# refpcu_mask = data.field('PCUID') == ref_pcu

        refpcu_mask = pcuid != 2
        all_time_ref = np.asarray(time[refpcu_mask], dtype=np.float64)
        all_energy_ref = np.asarray(channel[refpcu_mask], dtype=np.float64)
        all_rate_ref = None
        all_err_ref = None

    ###########################################################
    ## If separate ref band file exists, read in ref band data
    ###########################################################

    else:
        print("Ref band file : %s" % meta_dict['ref_file'])

        try:
            ref_fits_hdu = fits.open(meta_dict['ref_file'])
        except IOError:
            print("\tERROR: File does not exist: %s" % meta_dict['ref_file'])
            sys.exit()

        # ref_header = ref_fits_hdu[0].header	 ## Header info is in ext 0
        ref_data = ref_fits_hdu[1].data  ## Data is in ext 1
        ref_fits_hdu.close()

        ## Correcting times to make them at the front of the bin.
        all_time_ref = np.asarray(ref_data.field('TIME'), dtype=np.float64) \
                - meta_dict['dt'] / 2.0
        all_rate_ref = np.asarray(ref_data.field('RATE'), dtype=np.float64)
        all_err_ref = np.asarray(ref_data.field('ERROR'), dtype=np.float64)
        all_energy_ref = None

        ci_final_time = final_time
        ci_start_time = start_time
        ref_start_time = all_time_ref[0]
        ref_final_time = all_time_ref[-1]
    #    print("Ref start: %.15f" % ref_start_time)
    #    print("Ref final: %.15f" % ref_final_time)

        if ci_start_time > ref_start_time:
            index = np.min(np.where(all_time_ref > ci_start_time))
            start_time = all_time_ref[index]
            all_time_ref = all_time_ref[index:]
        else:
            start_time = ref_start_time

        if ci_final_time < ref_final_time:
            ## I don't think this has actually been tested...
            print("WARNING: I don't think this segment time selection has "\
                    "been tested.")
            index = np.max(np.where(all_time_ref < ci_final_time))
            # print("\t", index)
            final_time = all_time_ref[index]
            all_time_ref = all_time_ref[0:index]
        else:
            final_time = ref_final_time

    seg_end_time = start_time + meta_dict['n_seconds']

#    print("Start: %.15f" % start_time)
#    print("Final: %.15f" % final_time)

    ############################
    ## Looping through segments
    ############################

    print("Segments computed:")

    while (seg_end_time + (meta_dict['adjust_seg'] * meta_dict['dt'])) <= \
            final_time:

        ## Adjusting segment length to artificially line up the QPOs
        seg_end_time += (meta_dict['adjust_seg'] * meta_dict['dt'])

        ## Get events for channels of interest
        time_ci = all_time_ci[np.where(all_time_ci < seg_end_time)]
        energy_ci = all_energy_ci[np.where(all_time_ci < seg_end_time)]

        ## Chop current segment off the rest of the list
        for_next_iteration_ci = np.where(all_time_ci >= seg_end_time)
        all_time_ci = all_time_ci[for_next_iteration_ci]
        all_energy_ci = all_energy_ci[for_next_iteration_ci]

        ## Get events for reference band
        time_ref = all_time_ref[np.where(all_time_ref < seg_end_time)]
        if not meta_dict['ref_file']:
            energy_ref = all_energy_ref[np.where(all_time_ref < seg_end_time)]
            rate_ref = [0]
        else:
            rate_ref = all_rate_ref[np.where(all_time_ref < seg_end_time)]
            err_ref = all_err_ref[np.where(all_time_ref < seg_end_time)]

        ## Chop current segment off the rest of the list
        for_next_iteration_ref = np.where(all_time_ref >= seg_end_time)
        all_time_ref = all_time_ref[for_next_iteration_ref]
        if not meta_dict['ref_file']:
            all_energy_ref = all_energy_ref[for_next_iteration_ref]
        else:
            all_rate_ref = all_rate_ref[for_next_iteration_ref]
            all_err_ref = all_err_ref[for_next_iteration_ref]

        ########################################################################
        ## At the end of a segment, populate light curve and make cross spectrum
        ########################################################################

        if len(time_ci) > 0 and \
                (len(time_ref) > 0 or
                (meta_dict['ref_file'] and \
                len(rate_ref) == meta_dict['n_bins'])):

            ##############################################################
            ## Populate the light curves for interest and reference bands
            ##############################################################

            rate_ci_2d = tools.make_2Dlightcurve(np.asarray(time_ci),
                    np.asarray(energy_ci), meta_dict['n_bins'],
                    meta_dict['detchans'], start_time, seg_end_time)

            if not meta_dict['ref_file']:
                rate_ref_2d = tools.make_2Dlightcurve( np.asarray(time_ref),
                        np.asarray(energy_ref), meta_dict['n_bins'],
                        meta_dict['detchans'], start_time, seg_end_time)
                ## Stack the reference band
                rate_ref = stack_reference_band(rate_ref_2d, instrument="PCA",
                        obs_epoch=meta_dict['obs_epoch'])

            ## Save the reference band light curve to a text file
        	# out_file="./GX339-BQPO_ref_lc.dat"
        	# f_handle = file(out_file, 'a')
        	# np.savetxt(f_handle, rate_ref)
        	# f_handle.close()

            dt_seg = (seg_end_time - start_time) / float(meta_dict['n_bins'])
            df_seg = 1.0 / (meta_dict['n_bins'] * dt_seg)

            ###########################
            ## Make the cross spectrum
            ###########################

            cs_seg, ci_seg, ref_seg, f_lag_seg = make_cs(rate_ci_2d, rate_ref,
                                                          meta_dict, dt_seg)

            ###################################################################
            ## Compute variance and rms of the positive-frequency power in the
            ## reference band.
            ###################################################################

            if meta_dict['ref_file']:
                IR_poisson = np.sum(err_ref ** 2) / float(len(err_ref))
                print(IR_poisson)
                absrms_pow = ref_seg.power[0:meta_dict['n_bins'] / 2 + 1] * \
                        (2.0 * dt_seg / np.float(meta_dict['n_bins']))

            else:
                absrms_pow = raw_to_absrms(ref_seg.power[0:meta_dict['n_bins'] \
                        / 2 + 1], ref_seg.mean_rate, meta_dict['n_bins'],
                        dt_seg, noisy=True)

            var, rms = var_and_rms(absrms_pow, df_seg)

            ######################################################
            ## Only keep and use segments where the variance > 0.
            ######################################################

            # if var >= 0.0:

            dt_whole = np.append(dt_whole, dt_seg)
            df_whole = np.append(df_whole, df_seg)
            # print("%.3f\t%.1f\t%.1f" % \
            #       (np.sum(ci_seg.mean_rate[15:27]) / \
            #        np.sum(ci_seg.mean_rate[2:7]),
            #       np.sum(ci_seg.mean_rate[15:27]),
            #       np.sum(ci_seg.mean_rate[2:7])))

            ## Append segment to arrays
            cs_whole.total = np.dstack((cs_whole.total, cs_seg.total))
            cs_whole.fund = np.dstack((cs_whole.fund, cs_seg.fund))
            cs_whole.harm = np.dstack((cs_whole.harm, cs_seg.harm))
            ci_whole.mean_rate_array = np.hstack((ci_whole.mean_rate_array,
                    np.reshape(ci_seg.mean_rate, (meta_dict['detchans'],
                    1))))
            ref_whole.power_array = np.hstack((ref_whole.power_array,
                    np.reshape(ref_seg.power, (meta_dict['n_bins'], 1))))
            ref_whole.mean_rate_array = np.append(ref_whole.mean_rate_array,
                    ref_seg.mean_rate)
            ref_whole.var_array = np.append(ref_whole.var_array, var)

            ## Sum across segments -- arrays, so it adds by index
            exposure += (seg_end_time - start_time)
            n_seg += 1
            ci_whole.mean_rate += ci_seg.mean_rate
            ref_whole.mean_rate += ref_seg.mean_rate
            ci_whole.power += ci_seg.power
            ref_whole.power += ref_seg.power
            ref_whole.filt_qpo_power += ref_seg.filt_qpo_power
            f_lag_whole.cross += f_lag_seg.cross
            f_lag_whole.pow1 += f_lag_seg.pow1
            f_lag_whole.pow2 += f_lag_seg.pow2

            if n_seg % print_iterator == 0:
                print("\t", n_seg)

            if test is True and n_seg == 1:  # For testing
                break
            # else:
            #     print("Neg var")
            #     print(var)
            #     print("Start: %.15f" % start_time)
            #     print("End: %.15f" % seg_end_time)
            #     print("dt seg: %.15f" % dt_seg)
            #     # print(" ! %.3f\t%.1f\t%.1f !" % \
            #     #       (np.sum(ci_seg.mean_rate[15:27]) / \
            #     #        np.sum(ci_seg.mean_rate[2:7]),
            #     #       np.sum(ci_seg.mean_rate[15:27]),
            #     #       np.sum(ci_seg.mean_rate[2:7])))

            start_time = seg_end_time
            seg_end_time += meta_dict['n_seconds']

        ## This next bit deals with gappy data
        else:
            start_time = max(all_time_ci[0], all_time_ref[0])
            seg_end_time = start_time + meta_dict['n_seconds']
        #
        # else:
        #     start_time = seg_end_time
        #     seg_end_time += meta_dict['n_seconds']

        ## End of 'if there are counts in this segment'

    ## End of while-loop

    cs_whole.total = cs_whole.total[:,:,1:]
    cs_whole.fund = cs_whole.fund[:,:,1:]
    cs_whole.harm = cs_whole.harm[:,:,1:]
    ref_whole.var_array = ref_whole.var_array[1:]
    ci_whole.power_array = ci_whole.power_array[:,:,1:]
    ci_whole.mean_rate_array = ci_whole.mean_rate_array[:,1:]
    ref_whole.power_array = ref_whole.power_array[:,1:]
    ref_whole.mean_rate_array = ref_whole.mean_rate_array[1:]
    # print(dt_whole)
    # print(df_whole)

    return cs_whole, ci_whole, ref_whole, f_lag_whole, n_seg, dt_whole, \
           df_whole, exposure


################################################################################
def get_background(bkgd_file="evt_bkgd_rebinned.pha"):
    """
    Get the background count rate from a background spectrum file.

    Parameters
    ----------
    bkgd_file : str
        The full path name of the .pha file containing the background energy
        spectrum for the channels of interest, with energy channels binned in
        the same way as the data file. [evt_bkgd_rebinned.pha]

    Returns
    -------
    rate : np.array of floats
        1-D array of the count rate per energy channel of the background energy
        spectrum for the channels of interest, in cts/s.

    """

    try:
        fits_hdu = fits.open(bkgd_file)
    except IOError:
        print("\tERROR: File does not exist: %s" % bkgd_file)
        sys.exit()

    header = fits_hdu[1].header
    data = fits_hdu[1].data
    fits_hdu.close()

    exposure = np.float(header['EXPOSURE'])
    counts = data.field('COUNTS')

    rate = counts / exposure

    return rate

#
# ################################################################################
# def compute_coherence(in_table):
#     """
#     Compute the raw coherence of the cross spectrum. Coherence equation from
#     Uttley et al 2014 eqn 11, bias term equation from footnote 4 on same page.
#
#     Parameters
#     ----------
#     in_table : astropy Table
#
#     Returns
#     -------
#     coherence : np.array of floats
#         The raw coherence of the cross spectrum. (Uttley et al 2014, eqn 11)
#         Size = n_bins/2+1 (if avg over energy) or detchans (if avg over freq).
#
#     """
#
#     ## Reshaping (broadcasting) the ref to have same size as ci
#     # if np.shape(power_ref) != np.shape(power_ci):
#     #     power_ref = np.resize(np.repeat(power_ref, np.shape(power_ci)[1]),
#     #             np.shape(power_ci))
#
#     powers = in_table['POWER_1'] * in_table['POWER_2']
#     temp_2 = in_table['CROSS'] * np.conj(in_table['CROSS'])
#     with np.errstate(all='ignore'):
#         coherence = np.where(powers != 0, temp_2 / powers, 0)
#     # print "Coherence shape:", np.shape(coherence)
#     # print coherence
#     return np.real(coherence)
#
#
# ################################################################################
# def get_phase_err(in_table, n_range=1):
#     """
#     Compute the error on the complex phase (in radians) via the coherence.
#     Power should NOT be Poisson-noise-subtracted or normalized.
#
#     Parameters
#     ----------
#     in_table : Astropy table
#
#
#     n_range : int or np.array of ints
#         Number of frequency bins averaged over per new frequency bin for lags.
#         For energy lags, this is the number of frequency bins averaged over. For
#         frequency lags not re-binned in frequency, this is 1. Same as K in
#         equations in Section 2 of Uttley et al. 2014.
#
#     Returns
#     -------
#     phase_err : np.array of floats
#         1-D array of the error on the phase of the lag.
#
#     """
#
#     # coherence = np.where(a != 0, cs_avg * np.conj(cs_avg) / a, 0)
#     coherence = compute_coherence(in_table)
#
#     with np.errstate(all='ignore'):
#         phase_err = np.sqrt(np.where(coherence != 0, (1 - coherence) /
#                 (2 * coherence * n_range * in_table.meta['SEGMENTS']), 0))
#
#     return phase_err
#
#
# ################################################################################
# def phase_to_tlags(phase, f):
#     """
#     Convert a complex-plane cross-spectrum phase (in radians) to a time lag
#     (in seconds).
#
#     Parameters
#     ----------
#     phase : float or np.array of floats
#         1-D array of the phase of the lag, in radians.
#
#     f : float or np.array of floats
#         1-D array of the Fourier frequency of the cross-spectrum, in Hz.
#
#     Returns
#     -------
#     tlags : float or np.array of floats
#         1-D array of the time of the lag, in seconds.
#
#     """
#
#     if np.shape(phase) != np.shape(f):
#         ## Reshaping (broadcasting) f to have same size as phase
#         f = np.resize(np.repeat(f, np.shape(phase)[1]), np.shape(phase))
#
#     assert np.shape(phase) == np.shape(f), "ERROR: Phase array must have same "\
#             "dimensions as frequency array."
#
#     with np.errstate(all='ignore'):
#         tlags = np.where(f != 0, phase / (2.0 * np.pi * f), 0)
#
#     return tlags
#
#
# ################################################################################
# def calc_freq_lags(in_table, ccf_out_file):
#     """
#     Calculating the frequency lags over an energy range.
#     See steps in Uttley et al 2014 section 2.2.1.
#
#     Parameters
#     ----------
#
#
#     """
#     print(np.shape(in_table['CROSS']))
#     fake_err = np.ones(len(in_table['CROSS']))
#     ## Re-bin the PSDs and cross-spectrum in frequency
#     ## Uttley et al 2014 section 2.2.1 step 4
#     rb_freq_cs, rb_cs, temp1, bin_min, bin_max = rebin(in_table['FREQUENCY'],
#                                       in_table['CROSS'], fake_err, 1.04)
#     rb_freq_pow1, rb_pow_1, temp1, bin_min, bin_max = rebin(in_table['FREQUENCY'],
#                                              in_table['POWER_1'], fake_err,
#                                              1.04)
#     rb_freq_pow2, rb_pow_2, temp1, bin_min, bin_max = rebin(in_table['FREQUENCY'],
#                                              in_table['POWER_2'],fake_err, 1.04)
#
#     assert rb_freq_cs.all() == rb_freq_pow1.all(), "Something went wrong in re-binning " \
#                                         "for frequency lags."
#     # assert rb_freq_cs == rb_freq_powref, "Something went wrong in re-binning " \
#     #                                      "for frequency lags."
#     # assert bin_cs == bin_powci, "Something went wrong in re-binning for " \
#     #                             "frequency lags."
#     # assert bin_cs == bin_powref, "Something went wrong in re-binning for " \
#     #                              "frequency lags."
#
#     binning = bin_max - bin_min
#     rb_table = Table()
#     rb_table.meta = in_table.meta
#     rb_table['FREQUENCY'] = Column(rb_freq_cs)
#     rb_table['CROSS'] = Column(rb_cs)
#     rb_table['POWER_1'] = Column(rb_pow_1)
#     rb_table['POWER_2'] = Column(rb_pow_2)
#
#     ############################################
#     ## Getting lag and error for frequency lags
#     ############################################
#
#     ## Negative sign is so that a positive lag is a hard energy lag
#
#     freq_lags = Table()
#     freq_lags.meta = in_table.meta
#     freq_lags['FREQUENCY'] = Column(rb_table['FREQUENCY'], unit='Hz',
#                                     description='Fourier frequency')
#     freq_lags['PHASE_LAG'] = Column(-np.arctan2(rb_table['CROSS'].imag,
#                                                 rb_table['CROSS'].real),
#                                     unit='rad', description='Phase lag')
#     freq_lags['PHASE_ERR'] = Column(get_phase_err(rb_table, binning),
#                                     unit='rad',description='Error on phase lag')
#     freq_lags['TIME_LAG'] = Column(phase_to_tlags(freq_lags['PHASE_LAG'],
#                                                   rb_table['FREQUENCY']),
#                                    unit='s', description='Time lag')
#     freq_lags['TIME_ERR'] = Column(phase_to_tlags(freq_lags['PHASE_ERR'],
#                                                   rb_table['FREQUENCY']),
#                                    unit='s', description='Error on time lag')
#
#     out_file = ccf_out_file.replace("cross_correlation/out_ccf",
#             "lag_spectra/out_lags")
#     out_file = out_file.replace(".", "_lag-freq.")
#     out_dir = out_file[0:out_file.rfind("/")+1]
#     subprocess.call(['mkdir', '-p', out_dir])
#     print("Lag-freq written to: %s" % out_file)
#
#     freq_lags.write(out_file, format='fits')


################################################################################
def save_for_lags(out_file, in_file, meta_dict, cs_avg, ci, ref, f_lag):
    """
    Saving header data, the cross spectrum, CoI power spectrum, and reference
    band power spectrum to a FITS file to use in the program make_lags.py to get
    cross-spectral lags. Cross spectra and power spectra are raw, as in un-
    normalized, and not Poisson-noise-subtracted.

    Parameters
    ----------
    out_file : str
        The name the FITS file to write the cross spectrum and power spectra to,
        for computing the lags.

    in_file : str
        The name of the data file (or filename containing list of data files).

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    cs_avg : np.array of complex numbers
        2-D array of the averaged cross spectrum. Size = (n_bins, detchans).

    ci : ccf_lc.Lightcurve object
        The channel of interest light curve. Must already have mean_rate and
        pos_power assigned.

    ref : ccf_lc.Lightcurve object
        The reference band light curve. Must already have mean_rate, rms, and
        pos_power assigned.

    f_lag : ccf_lc.LagFreq object
        The cross spectrum and power spectra of two bands for lag-frequency
        spectra. Sizes are n_bins/2+1 (only for positive Fourier frequencies).

    Returns
    -------
    nothing, but writes to a file.

    """
    ## Computing the Fourier frequencies
    freq = fftpack.fftfreq(meta_dict['n_bins'], d=np.mean(meta_dict['dt']))
    nyq_index = meta_dict['n_bins'] / 2

    ## Only keeping the parts associated with positive Fourier frequencies
    freq = np.abs(freq[0:nyq_index + 1])  ## because it slices at end-1, and we
            ## want to include 'nyq_index'; abs is because the nyquist freq is
            ## both pos and neg, and we want it pos here.
    cs_avg = cs_avg[0:nyq_index + 1, :]

    flx2xsp_out(out_file, ref, freq, meta_dict['n_seg'], meta_dict['n_bins'])

    out_file = out_file.replace("cross_correlation/out_ccf",
            "lag_spectra/out_lags")
    out_file = out_file.replace(".", "_cs.")
    out_dir = out_file[0:out_file.rfind("/")+1]
    subprocess.call(['mkdir', '-p', out_dir])
    print("Lag stuff written to: %s" % out_file)

    out_table = Table()
    out_table.add_column(Column(data=freq, name='FREQUENCY', unit='Hz'))
    out_table.add_column(Column(data=cs_avg, name='CROSS'))
    out_table.add_column(Column(data=ci.pos_power, name='POWER_CI'))
    out_table.add_column(Column(data=ref.pos_power, name='POWER_REF'))
    out_table.add_column(Column(data=f_lag.cross, name='LAGF_CS'))
    out_table.add_column(Column(data=f_lag.pow1, name='POWER_1'))
    out_table.add_column(Column(data=f_lag.pow2, name='POWER_2'))

    out_table.meta['TYPE'] = "Cross spectrum and power spectra, saved for lags."
    out_table.meta['DATE'] = str(datetime.now())
    out_table.meta['EVTLIST'] = in_file
    out_table.meta['DT'] = np.mean(meta_dict['dt'])
    out_table.meta['DF'] = np.mean(meta_dict['df'])
    out_table.meta['N_BINS'] = meta_dict['n_bins']
    out_table.meta['SEGMENTS'] = meta_dict['n_seg']
    out_table.meta['SEC_SEG'] = meta_dict['n_seconds']
    out_table.meta['EXPOSURE'] = meta_dict['exposure']
    out_table.meta['DETCHANS'] = meta_dict['detchans']
    out_table.meta['RATE_CI'] = str(ci.mean_rate.tolist())
    out_table.meta['RATE_REF'] = ref.mean_rate
    out_table.meta['RMS_REF'] = float(ref.rms)
    out_table.meta['NYQUIST'] = meta_dict['nyquist']
    out_table.meta['ADJUST'] = "%s" % str(meta_dict['adjust_seg'])
    out_table.meta['RATE_1'] = np.mean(ci.mean_rate[2:8])
    out_table.meta['RATE_2'] = np.mean(ci.mean_rate[13:22])
    out_table.write(out_file, overwrite=True, format='fits')


################################################################################
def tophat_filt(cs_avg, ref, ci, freq, meta_dict, lo_freq, hi_freq):
    """
    Apply a tophat filter to the averaged cross-spectrum per energy channel (in
    frequency space). Any cross spectrum amplitudes above or below pulse_freq
    get zeroed out.

    Parameters
    ----------
    cs_avg : np.array of complex numbers
        2-D array of the cross spectrum, in frequency space, to be filtered.
        Size = (n_bins, detchans).

    ref : ccf_lc.Lightcurve object
        The reference band light curve. Must already have power assigned.

    ci : ccf_lc.Lightcurve object
        The channel of interest light curve. Must already have power assigned.

    freq : np.array of floats
        1-D array of the positive Fourier frequencies, including Nyquist.
        Size = (n_bins/2).

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    lo_freq : float
        The lower frequency bound for filtering the cross spectrum, in Hz. [0.0]

    hi_freq : float
        The upper frequency bound for filtering the cross spectrum, in Hz. [0.0]

    Returns
    -------
    full_filt_cs_avg : np.array of complex numbers
        2-D array of the cross spectrum, zeroed out at non-filtered frequencies.
        Size = (n_bins, detchans)

    signal_ci_pow : np.array of floats
        2-D array of the filtered power spectrum in the channels of interest.
        Size = (n_bins, detchans)

    signal_ref_pow : np.array of floats
        1-D array of the filtered power spectrum in the reference band.
        Size = (n_bins)

    """

    ## Get the indices of the beginning and end of the signal
    min_freq_mask = freq < lo_freq  # we want the last 'True' element
    max_freq_mask = freq > hi_freq  # we want the first 'True' element
    j_min = list(min_freq_mask).index(False)
    j_max = list(max_freq_mask).index(True)

    ## Make zeroed arrays to replace with
    zero_front = np.zeros((j_min, meta_dict['detchans']), dtype=np.complex128)
    zero_end = np.zeros((meta_dict['n_bins']/2 - j_max, meta_dict['detchans']),
            dtype=np.complex128)

    ## Concatenate the arrays together
    filt_cs_avg = np.concatenate((zero_front,
            cs_avg[j_min:j_max, :], zero_end), axis=0)
    print(np.shape(filt_cs_avg))
    ## Want to also keep negative frequencies
    full_filt_cs_avg = np.vstack((filt_cs_avg, filt_cs_avg[::-1,:]))
    print(np.shape(full_filt_cs_avg))
    ## Check that the original array is the same shape as the filtered one
    assert np.shape(cs_avg) == np.shape(full_filt_cs_avg), \
            "ERROR: Frequency-filtered cross spectrum does not have the same "\
            "size as the original cross spectrum. Something went wrong."

    ## Extract only the signal frequencies of the mean powers
    signal_ci_pow = np.vstack((ci.power[j_min:j_max, :],
                               ci.power[meta_dict['n_bins']-\
                                        j_max:meta_dict['n_bins']-j_min, :]))
    signal_ref_pow = np.hstack((ref.power[j_min:j_max],
                                ref.power[meta_dict['n_bins']-\
                                          j_max:meta_dict['n_bins']-j_min]))

    return full_filt_cs_avg, signal_ci_pow, signal_ref_pow


################################################################################
def tie_harmonic_centroid(model):
    """
    Used in optimal filtering if filtering the harmonic as well.

    How to use:
    tied_parameters = {'x_0_1': tie_harmonic_centroid}
    then within model definition, use tied=tied_parameters

    Parameters
    ----------
    model : an astropy model
    Returns
    -------
    x_0_1 : the centroid frequency of the harmonic.
    """
    x_0_1 = 2. * model.x_0_0
    return x_0_1


################################################################################
def phase_angle(complex_number):
    return np.arctan2(complex_number.imag, complex_number.real)


################################################################################
def vecrotate(theta, complex_number):
    # print("Theta:", theta)
    # print("Before rotation, 4th segment")
    # print("Abs:", np.abs(complex_number[3]))
    # print("Angle:", phase_angle(complex_number[3]))
    x = complex_number.real
    y = complex_number.imag
    xrot = x * np.cos(theta) - y * np.sin(theta)
    yrot = x * np.sin(theta) + y * np.cos(theta)
    rotated_complex_number = xrot + yrot*1j
    # print type(rotated_complex_number)
    # print type(rotated_complex_number[0])
    # print np.shape(rotated_complex_number)
    # print("After rotation, 4th segment")
    # print("Abs:", np.abs(rotated_complex_number[3]))
    # print("Angle:", phase_angle(rotated_complex_number[3]))
    return rotated_complex_number


################################################################################
def cs_to_ccf_w_err(cs_array, meta_dict, ref):
    """
    Take the iFFT of the cross spectrum to get the cross-correlation function,
    and compute the error on the cross-correlation function. Compute the
    standard error on each ccf bin from the segment-to-segment variations.
    This error is not correlated between energy bins but may be correlated
    between time bins.

    Keep in mind that if adjusting segments to line up QPOs that shift in
    frequency between segments, the rms of the avg ref power spectrum and the
    avg of the rmses of each segment's ref power spectrum are not the same. This
    has been accounted for here.

    S. Vaughan 2013, "Scientific Inference", equations 2.3 and 2.4.

    Parameters
    ----------
    cs_avg : np.array of complex numbers
        3-D array of the cross-spectrum. Size = (n_bins, detchans, n_seg).

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    ref : ccf_lc.Lightcurve object
        The reference band light curve and power.

    Returns
    -------
    ccf_avg, standard_error : np.arrays of floats
        2-D arrays of the cross-correlation function of the channels of interest
        with the reference band, and the error on the ccf.
        Sizes = (n_bins, detchans).
    """

    assert np.shape(cs_array) == (meta_dict['n_bins'], meta_dict['detchans'],
            meta_dict['n_seg'])
    if len(ref.power) == meta_dict['n_bins']:
        nyq_index = meta_dict['n_bins'] / 2
        ref.pos_power = ref.power[0:nyq_index + 1]

    ## Take the IFFT of the cross spectrum to get the CCF
    ccf_array = fftpack.ifft(cs_array, axis=0).real

    ## Average across segments and normalize by rms of averaged reference band
    ## absolute-rms-normalized power spectrum
    ccf_avg = np.mean(ccf_array, axis=-1)
    ccf_avg *= (2.0 / np.float(meta_dict['n_bins']) / ref.rms)
    ccf_array *= (2.0 / np.float(meta_dict['n_bins']) / np.sqrt(ref.var_array))

    ## Compute the standard error on each ccf bin from the segment-to-segment
    ## variations.
    mean_ccf = np.mean(ccf_array, axis=2)
    ccf_resid = (ccf_array.T - mean_ccf.T).T

    ## Eqn 2.3 from S. Vaughan 2013, "Scientific Inference"
    sample_var = np.sum(ccf_resid**2, axis=2) / (meta_dict['n_seg'] - 1)

    ## Eqn 2.4 from S. Vaughan 2013, "Scientific Inference"
    standard_error = np.sqrt(sample_var / meta_dict['n_seg'])

    return ccf_avg, standard_error


################################################################################
def main(input_file, out_file, ref_band="", bkgd_file="./evt_bkgd_rebinned.pha",
        n_seconds=64, dt_mult=64, test=False, filtering=False, adjusting=False):
    """
    Read in one event list, split into reference band and channels of
    interest (CoI), make segments and populates them to give them length
    n_bins, compute the cross spectrum of each segment per energy channel and
    then averaged cross spectrum of all the segments per energy channel, and
    then compute the cross-correlation function (ccf) per energy channel.

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
        has the reference band for the whole data set. Gaps are ok. []

    n_seconds : float
        Number of seconds in each Fourier segment. Must be positive. [64]

    dt_mult : int
        Multiple of dt (dt is from data file) for timestep between bins. Must be
        positive. [64]

    bkgd_file : str
        Name of the background spectrum (in .pha format), with the same energy
        channel binning as the event list. [None]

    test : bool
        If true, only computes one segment of data. If false, runs like normal.
        [False]

    filtering : bool
        If true, filters the Fourier transforms with an optimal filter
        (hardcoded the file name). [False]

    adjusting : bool
        If true, the QPOs are lined up (adjusted) in frequency per obsID.
        [False]

    Returns
    -------
    nothing
    """

    #####################################################
    ## Idiot checks, to ensure that our assumptions hold
    #####################################################

    assert n_seconds > 0, "ERROR: Number of seconds per segment must be a "\
            "positive integer."
    assert dt_mult > 0, "ERROR: Multiple of dt must be a positive integer."

    ########################
    ## Read in data file(s)
    ########################

    if ".txt" in input_file or ".lst" in input_file or ".dat" in input_file:
        data_files = [line.strip() for line in open(input_file)]
        if not data_files:  ## If data_files is an empty list
            raise Exception("ERROR: No files in the list of event lists.")
    else:
        data_files = [input_file]

    if "stevens" in os.path.expanduser('~'):
        for i in range(len(data_files)):
            if "stevens" not in data_files[i]:
                data_files[i] = data_files[i].replace("abigail", "abigailstevens")

    if adjusting is True:
        # adjust_segments = [932, 216, 184, 570, 93, 346, 860, 533, -324]
        out_dir = os.path.dirname(out_file)
        basename = os.path.basename(out_file)
        prefix = basename.split('_')[0]
        if "test" in basename:
            prefix = basename.split('_')[1]
        adjust_file = out_dir+"/"+prefix+"_t"+str(dt_mult)+"_"+str(n_seconds)+\
						   "sec_adjustby.txt"
        if os.path.isfile(adjust_file):
            adjust_segments = [int(line.strip()) for line in open(adjust_file)]
        else:
            adjust_file = out_dir + "/" + prefix + "_t" + str(
                dt_mult) + "_" + str(int(n_seconds)) + \
                          "sec_adjustby.txt"
            if os.path.isfile(adjust_file):
                adjust_segments = [int(line.strip()) for line in
                                   open(adjust_file)]
            else:
                print("WARNING: adjustby.txt file does not exist. NOT "\
                      "adjusting segments to line up QPO.")
                adjusting = False
                adjust_segments = np.zeros(len(data_files))
                out_file = out_file.replace("_adj","")
    else:
        adjust_segments = np.zeros(len(data_files))

    #############################################
    ## Initialize; 'whole' is over one data file
    #############################################

    try:
        t_res = float(tools.get_key_val(data_files[0], 0, 'TIMEDEL'))
    except KeyError:
        t_res = float(tools.get_key_val(data_files[0], 1, 'TIMEDEL'))

    try:
        detchans = int(tools.get_key_val(data_files[0], 0, 'DETCHANS'))
    except KeyError:
        detchans = int(tools.get_key_val(data_files[0], 1, 'DETCHANS'))
    except IOError:
        detchans = 64

    meta_dict = {'dt': dt_mult * t_res,
                 't_res': t_res,
                 'n_seconds': n_seconds,
                 'df': 1.0 / np.float(n_seconds),
                 'nyquist': 1.0 / (2.0 * dt_mult * t_res),
                 'n_bins': int(n_seconds * 1.0 / (dt_mult * t_res)),
                 'detchans': detchans,
                 'filter': filtering,
                 'exposure': 0,
                 'ref_file': ref_band}

    print("\nDT = %f" % meta_dict['dt'])
    print("N_bins = %d" % meta_dict['n_bins'])
    print("Nyquist freq =", meta_dict['nyquist'])
    print("Testing?", test)
    print("Filtering?", meta_dict['filter'])
    print("Adjusting QPO?", adjusting)

    ci_total = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ci')
    ref_total = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ref')
    total_seg = 0
    cs_total = ccf_lc.Cross(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'])
    f_lag_total = ccf_lc.LagFreq(n_bins=meta_dict['n_bins'])
    dt_total = np.array([])
    df_total = np.array([])
    total_exposure = 0

    ###############################
    ## Loop through all data files
    ###############################

    for in_file, adj_seg in zip(data_files, adjust_segments):

        meta_dict['adjust_seg'] = adj_seg

        ############################################
        ## Read in data, compute the cross spectrum
        ############################################

        cs_whole, ci_whole, ref_whole, f_lag_whole, n_seg, \
                dt_whole, df_whole, exposure = fits_in(in_file, meta_dict, test)

        print("Segments for this file: %d\n" % n_seg)

        cs_total.total = np.dstack((cs_total.total, cs_whole.total))
        cs_total.fund = np.dstack((cs_total.fund, cs_whole.fund))
        cs_total.harm = np.dstack((cs_total.harm, cs_whole.harm))
        ref_total.var_array = np.append(ref_total.var_array,
                ref_whole.var_array)
        ref_total.power_array = np.hstack((ref_total.power_array,
                        ref_whole.power_array))
        ref_total.mean_rate_array = np.append(ref_total.mean_rate_array,
                        ref_whole.mean_rate_array)
        ci_total.mean_rate_array = np.hstack((ci_total.mean_rate_array,
                        ci_whole.mean_rate_array))
        dt_total = np.append(dt_total, dt_whole)
        df_total = np.append(df_total, df_whole)
        ci_total.mean_rate += ci_whole.mean_rate
        ref_total.mean_rate += ref_whole.mean_rate
        ci_total.power += ci_whole.power
        ref_total.power += ref_whole.power
        ref_total.filt_qpo_power += ref_whole.filt_qpo_power
        f_lag_total.cross += f_lag_whole.cross
        f_lag_total.pow1 += f_lag_whole.pow1
        f_lag_total.pow2 += f_lag_whole.pow2
        total_exposure += exposure
        total_seg += n_seg

    meta_dict['n_seg'] = total_seg
    meta_dict['exposure'] = total_exposure
    meta_dict['dt'] = dt_total
    meta_dict['df'] = df_total
    meta_dict['adjust_seg'] = adjust_segments
    meta_dict['nyquist'] = 1. / (2. * np.mean(dt_total))
    print("Mean dt: %.15f" % np.mean(dt_total))
    print("Mean df: %.10f\n" % np.mean(df_total))
    # print("Total segments: %d" % meta_dict['n_seg'])

    ## Remove the first zeros from stacked arrays
    cs_total.total = cs_total.total[:,:,1:]
    cs_total.fund = cs_total.fund[:,:,1:]
    cs_total.harm = cs_total.harm[:,:,1:]
    ci_total.mean_rate_array = ci_total.mean_rate_array[:,1:]
    ref_total.power_array = ref_total.power_array[:,1:]
    ref_total.mean_rate_array = ref_total.mean_rate_array[1:]
    ref_total.var_array = ref_total.var_array[1:]
    # print("Mean of var_array:", np.mean(ref_total.var_array))
    # print("Sqrt of var_array:", np.sqrt(ref_total.var_array))
    # print("Mean of sqrt of var_array:", np.mean(np.sqrt(ref_total.var_array)))

    ######################################
    ## Turn sums over segments into means
    ######################################

    ci_total.mean_rate /= np.float(meta_dict['n_seg'])
    ci_total.power /= np.float(meta_dict['n_seg'])
    ref_total.power /= np.float(meta_dict['n_seg'])
    ref_total.filt_qpo_power /= np.float(meta_dict['n_seg'])
    ref_total.mean_rate /= np.float(meta_dict['n_seg'])
    cs_avg_total = np.mean(cs_total.total, axis=-1)
    f_lag_total.cross /= np.float(meta_dict['n_seg'])
    f_lag_total.pow1 /= np.float(meta_dict['n_seg'])
    f_lag_total.pow2 /= np.float(meta_dict['n_seg'])

    assert ci_total.mean_rate.all() == \
           np.mean(ci_total.mean_rate_array, axis=-1).all()

    ci_total.pos_power = ci_total.power[0:meta_dict['n_bins']/2+1, :]
    ref_total.pos_power = ref_total.power[0:meta_dict['n_bins']/2+1]
    f_lag_total.cross = f_lag_total.cross[0:meta_dict['n_bins']/2+1]
    f_lag_total.pow1 = f_lag_total.pow1[0:meta_dict['n_bins']/2+1]
    f_lag_total.pow2 = f_lag_total.pow2[0:meta_dict['n_bins']/2+1]

    # freq = fftpack.fftfreq(meta_dict['n_bins'], d=np.mean(meta_dict['dt']))
    # freq = np.abs(freq[0:meta_dict['n_bins']/2+1])
    # flx2xsp_out(out_file, ref_total, freq, n_seg, n_bins)
    # ci_flx2xsp_out(out_file, ci_total, ref_total, freq, n_seg)

    ## Compute the variance and rms of the absolute-rms-normalized reference
    ## band power spectrum
    if meta_dict['filter']:
        absrms_ref_pow = raw_to_absrms(ref_total.filt_qpo_power,
                                       ref_total.mean_rate, meta_dict['n_bins'],
                                       np.mean(meta_dict['dt']), noisy=False)
    else:
        absrms_ref_pow = raw_to_absrms(ref_total.pos_power,
                                       ref_total.mean_rate, meta_dict['n_bins'],
                                       np.mean(meta_dict['dt']), noisy=True)


    ref_total.var, ref_total.rms = var_and_rms(absrms_ref_pow,
                np.mean(meta_dict['df']))
    print("Ref variance:", ref_total.var)
    print("Ref rms:", ref_total.rms)

    #############################################################
    ## Print the cross spectrum to a file, for plotting/checking
    #############################################################

    # cs_out = np.column_stack((fftpack.fftfreq(meta_dict['n_bins'],
    #         d=np.mean(meta_dict['dt'])), avg_cross_spec))
    # np.savetxt('cs_avg.dat', cs_out)

    #####################################################################
    ## Read in the background count rate from a background spectrum, and
    ## subtract from the mean count rate.
    #####################################################################

    if bkgd_file:
        bkgd_rate = get_background(bkgd_file)
    else:
        bkgd_rate = np.zeros(meta_dict['detchans'])

    ci_total.mean_rate -= bkgd_rate

    ## Need to use a background from ref PCU for the reference band...
    # ref_total.mean_rate -= np.sum(bkgd_rate[2:26])

    ####################################################################
    ## Save cross spectra and power spectra for computing lags later in
    ## lag_spectra/get_lags.py
    ####################################################################

    save_for_lags(out_file, input_file, meta_dict, cs_avg_total, ci_total,
                  ref_total, f_lag_total)


    ##############################################
    ## Compute ccf from cs, and compute error
    ##############################################

    if meta_dict['filter']:
        # cs_h_rotated = vecrotate(1.2883218366, cs_total.harm)
        # cs_combined = cs_total.fund + cs_h_rotated
        print("WARNING: just adding them to test on Type B!")
        cs_combined = cs_total.fund + cs_total.harm
        # print(np.shape(cs_combined))
        ccf_avg, ccf_error = cs_to_ccf_w_err(cs_combined, meta_dict,
                                             ref_total)
    else:
        ccf_avg, ccf_error = cs_to_ccf_w_err(cs_total.total, meta_dict,
                                             ref_total)

    print("Exposure_time = %.3f seconds" % meta_dict['exposure'])
    print("Total number of segments:", meta_dict['n_seg'])
    print("Mean rate for all of ci:", np.sum(ci_total.mean_rate))
    print("Mean rate for ci chan 6:", ci_total.mean_rate[6])
    print("Mean rate for ci chan 15:", ci_total.mean_rate[15])
    print("Mean rate for ref:", ref_total.mean_rate)

    ##########
    ## Output
    ##########

    # print(avg_cross_spec[529:532,3])
    print(ccf_avg[0:4,3])

    if len(data_files) == 1:
        file_description = "Cross-correlation function of one observation"
    else:
        file_description = "Cross-correlation function of multiple observations"

    # print(ccf_avg[0,0])
    # print(ccf_avg[0,2])
    # print(ccf_avg[0,15])

    # if not test and adjusting and len(data_files) == 9:
    #     assert round(ccf_avg[0,0], 12) == 0.117937948428
    #     assert round(ccf_avg[0,2], 11) == 9.22641398474
    #     assert round(ccf_avg[0,15], 11) == 1.76422640304
    #     print("Passed!")
    # elif test and adjusting and len(data_files) == 9:
    #     assert round(ccf_avg[0,0], 12) == 0.106747663439
    #     assert round(ccf_avg[0,2], 11) == 9.56560710672
    #     assert round(ccf_avg[0,15], 11) == 0.88144237181
    #     print("Passed!")
    # else:
    #     print("Do not have values to compare against.")

    fits_out(out_file, input_file, bkgd_file, meta_dict, ci_total.mean_rate,
            ref_total.mean_rate, float(ref_total.rms), ccf_avg, ccf_error,
            file_description)


################################################################################
if __name__ == "__main__":

    #########################################
    ## Parse input arguments and call 'main'
    #########################################

    parser = argparse.ArgumentParser(usage="python ccf.py infile outfile "\
            "[OPTIONAL ARGUMENTS]", description=__doc__, epilog="For optional "\
            "arguments, default values are given in brackets at the end of the"\
            " description.")

    parser.add_argument('infile', help="Could be either: 1) The full path of "\
            "the .fits input file containing an event list with time in "\
            "column 1, energy channel in column 2, detector ID in column 3; "\
            "or 2) The full path of the (ASCII/txt) file with a list of the "\
            "input files (as described in 1). One file per line. The code "\
            "assumes that if a separate reference band file is not given "\
            "with --ref flag, both reference band and channels of interest "\
            "are in the event list(s), and channels of interest = PCU 2, ref "\
            "band = all other PCUs.")

    parser.add_argument('outfile', help="The name the FITS file to write the "\
            "cross-correlation function to.")

    parser.add_argument('--ref', dest='ref_band_file', default=None,
            help="Name of FITS optical or IR data file for reference band. "\
            "For now, this data must already be populated as a light curve, "\
            "with time in column 1 and count rate in column 2. Times must be "\
            "in same units as the data in 'infile'.")

    parser.add_argument('-b', '--bkgd', dest='bkgd_file',
            default=None, help="Name of the background spectrum (in .pha "\
            "format), with the same energy channel binning as the event list."\
            " [None]")

    parser.add_argument('-n', '--n_seconds', type=type_positive_float,
            default=64, dest='n_seconds', help="Number of seconds in each "\
            "Fourier segment. Must be a positive number. [64]")

    parser.add_argument('-m', '--dt_mult', type=type_positive_int,
            default=64, dest='dt_mult', help="Multiple of dt (dt is from data "\
            "file) for timestep between bins. Must be a positive "\
            "integer. [64]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
            dest='test', help="Int flag: 0 if computing all segments, 1 if "\
            "only computing one segment for testing. [0]")

    parser.add_argument('-f', '--filter', default=0, type=int, choices={0,1},
            dest='filter', help="Int flag: 0 if not filtering the Fourier "\
            "transforms, 1 if applying a hardcoded optimal filter to the FTs."\
            " [0]")

    parser.add_argument('-a', '--adjust', default=1, type=int, choices={0,1},
            dest='adjust', help="Int flag: 0 if not adjusting segment length, "\
            "1 if adjusting segment length to line up peak frequency of QPOs."\
            " [1]")

    # parser.add_argument('-a', '--adjust', default=False, action='store_true',
    #         dest='adjust', help="If present, artificially adjusts the "\
    #         "frequency of the QPO by changing the segment length. [False]")

    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True
    adjust = False
    if args.adjust == 1:
        adjust = True
    filter = False
    if args.filter == 1:
        filter = True

    main(args.infile, args.outfile, ref_band=args.ref_band_file,
            bkgd_file=args.bkgd_file, n_seconds=args.n_seconds,
            dt_mult=args.dt_mult, test=test, filtering=filter, adjusting=adjust)

################################################################################
