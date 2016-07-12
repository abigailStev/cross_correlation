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
import matplotlib.pyplot as plt

## These are things I've written.
## Their directories are in my PYTHONPATH bash environment variable.
import tools  ## in whizzy_scripts
import ccf_lightcurves as ccf_lc  ## in cross_correlation

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2014-2016"


################################################################################
def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx == len(array) or np.fabs(value - array[idx-1]) < \
            np.fabs(value - array[idx]):
        return array[idx-1], idx-1
    else:
        return array[idx], idx


################################################################################
def fits_out(out_file, in_file, bkgd_file, meta_dict, mean_rate_ci,
        mean_rate_ref, rms_ref, ccf, ccf_error, lo_freq, hi_freq, file_desc):
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

    lo_freq : float
        Low frequency bound of cross spectrum filter, in Hz.

    hi_freq : float
        High frequency bound of cross spectrum filter, in Hz.

    file_desc : str
        Short description of the output file.

    Returns
    -------
    Nothing, but writes to out_file.fits.
    """

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
    out_table.meta['FILTFREQ'] = "%f:%f" % (lo_freq, hi_freq)
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
def make_cs(rate_ci, rate_ref, meta_dict):
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

    ## print(fft_data_ref[528,3])
    ## print(fft_data_ci[254,3])
    # fft_data_ref[528:561,:] = fft_data_ci[254:287,:]
    ## print(fft_data_ref[528,3])
    ## print("\n")
    ## print(fft_data_ref[-529,3])
    ## print(fft_data_ci[-255,3])
    # fft_data_ref[-561:-528,:] = fft_data_ci[-287:-254,:]
    ## print(fft_data_ref[-529,3])
    ## Computing the cross spectrum from the fourier transform
    cs_seg = np.multiply(fft_data_ci, np.conj(fft_data_ref))
    ## print(cs_seg[529:532,3])
    ## exit()
    return cs_seg, ci_seg, ref_seg

#
# ################################################################################
# def each_segment(time_ci, time_ref, energy_ci, energy_ref, meta_dict,\
#     start_time, end_time):
#     """
#     Turn the event list into a populated histogram, stack the reference band,
#     and make the cross spectrum, per segment of light curve.
#
#     Parameters
#     ----------
#     time_ci : np.array of floats
#         1-D array of the photon arrival times of events in this segment for the
#         channel of interest.
#
#     time_ref : np.array of floats
#         1-D array of the photon arrival times of events in this segment for the
#         reference band.
#
#     energy_ci : np.array of ints
#         1-D array of the energy channels of events in this segment for the
#         channel of interest.
#
#     energy_ref : np.array of ints
#         1-D array of the energy channels of events in this segment for the
#         reference band.
#
#     meta_dict : dict
#         Dictionary of necessary meta-parameters for data analysis.
#
#     start_time : float
#         Starting time of the segment (front of bin), in whatever units time_ci
#         and time_ref are in.
#
#     end_time : float
#         End time of the segment (back of bin), in whatever units time_ci,
#         time_ref, and start_time are in.
#
#     Returns
#     -------
#     cs_seg : np.array of complex numbers
#         The cross spectrum of the channels of interest with the reference band
#         for this segment of data. Size=(n_bins, detchans).
#
#     ci_seg : ccf_lc.Lightcurve object
#         The channels of interest light curve.
#
#     ref_seg : ccf_lc.Lightcurve object
#         The reference band light curve.
#
#     np.mean(rate_ci_2d) : float
#         The total mean count rate of the channels of interest.
#
#     """
#     assert len(time_ci) == len(energy_ci)
#     assert len(time_ref) == len(energy_ref)
#
#     ##############################################################
#     ## Populate the light curves for interest and reference bands
#     ##############################################################
#
#     rate_ci_2d = tools.make_2Dlightcurve(np.asarray(time_ci),
#         np.asarray(energy_ci), meta_dict['n_bins'], meta_dict['detchans'],
#         start_time, end_time)
#     rate_ref_2d = tools.make_2Dlightcurve( np.asarray(time_ref),
#         np.asarray(energy_ref), meta_dict['n_bins'], meta_dict['detchans'],
#         start_time, end_time)
#
#     ## Stack the reference band
#     rate_ref = stack_reference_band(rate_ref_2d, instrument="PCA",
#                                     obs_epoch=meta_dict['obs_epoch'])
#
#     ## Save the reference band light curve to a text file
# # 	out_file="./GX339-BQPO_ref_lc.dat"
# # 	f_handle = file(out_file, 'a')
# # 	np.savetxt(f_handle, rate_ref)
# # 	f_handle.close()
#
#     ###########################
#     ## Make the cross spectrum
#     ###########################
#
#     cs_seg, ci_seg, ref_seg = make_cs(rate_ci_2d, rate_ref, meta_dict)
#
#     return cs_seg, ci_seg, ref_seg, np.mean(rate_ci_2d)


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

    assert tools.power_of_two(meta_dict['n_bins']), "ERROR: n_bins must be a "\
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
    cs_whole = np.zeros((meta_dict['n_bins'], meta_dict['detchans'], 1),
            dtype=np.complex128)
    dt_whole = np.array([])
    df_whole = np.array([])
    exposure = 0
    # print(set(pcuid))
    # exit()
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

            ###########################
            ## Make the cross spectrum
            ###########################

            cs_seg, ci_seg, ref_seg = make_cs(rate_ci_2d, rate_ref, meta_dict)

            dt_seg = (seg_end_time - start_time) / float(meta_dict['n_bins'])
            df_seg = 1.0 / (meta_dict['n_bins'] * dt_seg)

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

            if var >= 0.0:

                dt_whole = np.append(dt_whole, dt_seg)
                df_whole = np.append(df_whole, df_seg)
                # print("%.3f\t%.1f\t%.1f" % \
                #       (np.sum(ci_seg.mean_rate[15:27]) / \
                #        np.sum(ci_seg.mean_rate[2:7]),
                #       np.sum(ci_seg.mean_rate[15:27]),
                #       np.sum(ci_seg.mean_rate[2:7])))

                ## Append segment to arrays
                cs_whole = np.dstack((cs_whole, cs_seg))
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

                if n_seg % print_iterator == 0:
                    print("\t", n_seg)
                    print(cs_seg[529:532,3])

                if test is True and n_seg == 1:  # For testing
                    break
            else:
                print("Neg var")
                print(var)
                print("Start: %.15f" % start_time)
                print("End: %.15f" % seg_end_time)
                print("dt seg: %.15f" % dt_seg)
                # print(" ! %.3f\t%.1f\t%.1f !" % \
                #       (np.sum(ci_seg.mean_rate[15:27]) / \
                #        np.sum(ci_seg.mean_rate[2:7]),
                #       np.sum(ci_seg.mean_rate[15:27]),
                #       np.sum(ci_seg.mean_rate[2:7])))

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

    cs_whole = cs_whole[:,:,1:]
    # ref_whole.rms_array = ref_whole.rms_array[1:]
    ref_whole.var_array = ref_whole.var_array[1:]
    # ci_whole.power_array = ci_whole.power_array[:,:,1:]
    ci_whole.mean_rate_array = ci_whole.mean_rate_array[:,1:]
    ref_whole.power_array = ref_whole.power_array[:,1:]
    ref_whole.mean_rate_array = ref_whole.mean_rate_array[1:]
    # print(dt_whole)
    # print(df_whole)
    # print(ref_whole.rms_array)
    # print(np.shape(ref_whole.rms_array))

    return cs_whole, ci_whole, ref_whole, n_seg, dt_whole, df_whole, \
            exposure


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


################################################################################
def save_for_lags(out_file, in_file, meta_dict, cs_avg, ci, ref, lo_freq=-1.0,
        hi_freq=-1.0):
    """
    Saving header data, the cross spectrum, CoI power spectrum, and reference
    band power spectrum to a FITS file to use in the program make_lags.py to get
    cross-spectral lags. Cross spectra and power spectra are raw, as in un-
    normalized.

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

    lo_freq, hi_freq : floats
        The lower and upper frequency bounds being filtered over, if applicable.
        Default values are -1.

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

    out_file = out_file.replace("cross_correlation/out_ccf",
            "lag_spectra/out_lags")
    out_file = out_file.replace(".", "_cs.")
    out_dir = out_file[0:out_file.rfind("/")+1]
    subprocess.call(['mkdir', '-p', out_dir])
    print("Output sent to: %s" % out_file)

    out_table = Table()
    out_table.add_column(Column(data=freq, name='FREQUENCY', unit='Hz'))
    out_table.add_column(Column(data=cs_avg, name='CROSS'))
    out_table.add_column(Column(data=ci.pos_power, name='POWER_CI'))
    out_table.add_column(Column(data=ref.pos_power, name='POWER_REF'))

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
    out_table.meta['FILTER'] = str(meta_dict['filter'])
    out_table.meta['FILTFREQ'] = "%f:%f" % (lo_freq, hi_freq)
    out_table.meta['ADJUST'] = "%s" % str(meta_dict['adjust_seg'])
    out_table.write(out_file, overwrite=True)


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
def optimal_filt(cs_avg, ref, ci, freq, meta_dict, lo_freq, hi_freq,
                 harmonic=False):

    # err_pow = ref.pos_power / np.sqrt(float(meta_dict['n_seg']))

    freq_mask = (freq > lo_freq) & (freq < hi_freq)
    power_ref_lim = ref.pos_power[freq_mask]
    freq_lim = freq[freq_mask]

    npn = power_ref_lim * freq_lim
    noise_init = powerlaws.PowerLaw1D(amplitude=1E8, x_0=1., alpha=-1.,
                                      bounds={'alpha': (-1.2, -0.8),
                                              'x_0': (0.8, 1.2)})
                                  #fixed={'x_0': True, 'alpha': True})
    qpo_init = models.Lorentz1D(amplitude=1E11, x_0=4.3240, fwhm=0.4863,
                                bounds={'fwhm': (0, 2.0),
                                        'x_0': (lo_freq, hi_freq)})
    #                             fixed={'fwhm': True})
    # print("WARNING: fwhm is frozen at 0.4863.")
    ## TODO: take FWHM from the power spectrum fit.

    if harmonic is False:
        qpo_model = noise_init + qpo_init
    else:
        tied_parameters = {'x_0_2': tie_harmonic_centroid}
        harmonic_init = models.Lorentz1D(amplitude=1E11, x_0=hi_freq-1.,
                                         fwhm=0.3, tied=tied_parameters,
                                         bounds={'fwhm': (0, 2.0)})
        qpo_model = noise_init + qpo_init + harmonic_init

    np.random.seed(0)
    fit_qpo = fitting.LevMarLSQFitter()
    qpo_and_noise = fit_qpo(qpo_model, freq_lim, npn)

    print(fit_qpo.fit_info['message'])

    if harmonic is False:
        qpo_filter_model = models.Lorentz1D(amplitude=qpo_and_noise.amplitude_1,
                                            x_0=qpo_and_noise.x_0_1,
                                            fwhm=qpo_and_noise.fwhm_1)
    else:
        qpo_filter_model = models.Lorentz1D(amplitude=qpo_and_noise.amplitude_1,
                                            x_0=qpo_and_noise.x_0_1,
                                            fwhm=qpo_and_noise.fwhm_1) + \
                           models.Lorentz1D(amplitude=qpo_and_noise.amplitude_2,
                                            x_0=qpo_and_noise.x_0_2,
                                            fwhm=qpo_and_noise.fwhm_2)

    temp1 = qpo_and_noise.x_0_1 - qpo_and_noise.fwhm_1 / 2
    temp2 = qpo_and_noise.x_0_1 + qpo_and_noise.fwhm_1 / 2
    if harmonic is True:
        temp3 = qpo_and_noise.x_0_2 - qpo_and_noise.fwhm_2 / 2
        temp4 = qpo_and_noise.x_0_2 + qpo_and_noise.fwhm_2 / 2

    plt.figure(figsize=(10,7.5))
    plt.plot(freq, ref.pos_power * freq, 'ko', label="Data")
    plt.plot(freq, qpo_and_noise(freq), label="Filter", lw=2)
    plt.vlines(temp1, ymin=0, ymax=1.1E11, color='magenta')
    plt.vlines(temp2, ymin=0, ymax=1.1E11, color='magenta')
    if harmonic is True:
        plt.vlines(temp3, ymin=0, ymax=1.1E11, color='purple')
        plt.vlines(temp4, ymin=0, ymax=1.1E11, color='purple')
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Raw power * frequency")
    plt.xlim(0, freq[-1])
    plt.legend(loc=2)
    plt.savefig("./optimal_filter.png")
    # plt.show()

    filter_ratio = np.where(freq != 0, qpo_filter_model(freq) / \
                            qpo_and_noise(freq), 0.)

    shortened_filter_ratio = filter_ratio[1:-1]
    full_filt = np.append(filter_ratio, shortened_filter_ratio[::-1])

    ## Filtering the reference band power spectrum
    ref_pow_filt = ref.power * full_filt

    ## Broadcasting the filter for detchans
    full_filt = np.resize(np.repeat(full_filt, meta_dict['detchans']),
            (meta_dict['n_bins'], meta_dict['detchans']))

    ## Filtering the cross spectrum and ci power spectrum
    cs_filt = (cs_avg.real * full_filt) + (1j * cs_avg.imag)
    ci_pow_filt = ci.power * full_filt

    return cs_filt, ci_pow_filt, ref_pow_filt


################################################################################
def filt_cs_to_ccf_w_err(cs_avg, meta_dict, ci, ref, lo_freq=0.0, hi_freq=0.0,
        filter_harmonic=False, noisy=True):
    """
    Filter the cross-spectrum in frequency space, take the iFFT of the
    filtered cross spectrum to get the cross-correlation function, and compute
    the error on the cross-correlation function. Note that error is definitely
    NOT independent between time bins due to the filtering! But is still
    independent between energy bins.

    Keep in mind that if adjusting segments to line up QPOs that shift in
    frequency between segments, the rms of the avg ref power spectrum and the
    avg of the rmses of each segment's ref power spectrum are not the same. This
    has been accounted for here.

    Parameters
    ----------
    cs_avg : np.array of complex numbers
        2-D array of the segment-averaged cross spectrum of the channels of
        interest with the reference band. Size = (n_bins, detchans).

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    ci : ccf_lc.Lightcurve object
        Channel of interest light curve, with (averaged) mean_rate and
        (averaged, raw) power assigned.

    ref : ccf_lc.Lightcurve object
        Reference band light curve, with (averaged) mean_rate, (averaged, raw)
        power, and rms (of averaged power in abs rms units) assigned.

    lo_freq : float
        The lower frequency bound for filtering the cross spectrum, in Hz. [0.0]

    hi_freq : float
        The upper frequency bound for filtering the cross spectrum, in Hz. [0.0]

    filter_harmonic : bool
        If true, the filter will include the harmonic. The lo_freq - hi_freq
        range must then include the harmonic! This isn't checked. If
        filtering=False, this does nothing. But then we wouldn't get to this
        method anyways. [False]

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
    nyq_index = meta_dict['n_bins'] / 2
    if len(ref.power) == meta_dict['n_bins']:
        ref.pos_power = ref.power[0:nyq_index + 1]

    ## Compute the Fourier frequencies
    freq_long = fftpack.fftfreq(meta_dict['n_bins'], d=np.mean(meta_dict['dt']))
    freq = np.abs(freq_long[0:nyq_index + 1])

    ## Apply optimal filter to cross-spectrum
    filtered_cs_avg, signal_ci_pow, signal_ref_pow = optimal_filt(cs_avg,
            ref, ci, freq, meta_dict, lo_freq, hi_freq,
            harmonic=filter_harmonic)

    ## Apply tophat filter to cross-spectrum
    # filtered_cs_avg, signal_ci_pow, signal_ref_pow = tophat_filt(cs_avg, ref,
    #         ci, freq, meta_dict, lo_freq, hi_freq)

    ## Poisson noise level in absolute rms^2 norm
    noise_ci = 2.0 * ci.mean_rate
    noise_ci[noise_ci <= 0] = 0.0  ## Can't have a negative noise
    noise_ref = 2.0 * ref.mean_rate

    ## If there's no noise in a (simulated) power spectrum, noise level = 0
    if not noisy:
        noise_ci = np.zeros(meta_dict['detchans'])
        noise_ref = 0

    ## Apply absolute rms2 normalization to power spectra, subtract noise
    signal_ci_pow = raw_to_absrms(signal_ci_pow, ci.mean_rate,
            meta_dict['n_bins'], np.mean(meta_dict['dt']), noisy=noisy)
    signal_ref_pow = raw_to_absrms(signal_ref_pow, ref.mean_rate,
            meta_dict['n_bins'], np.mean(meta_dict['dt']), noisy=noisy)

    ## If the power is negative, set it equal to zero.
    signal_ref_pow[signal_ref_pow < 0] = 0.
    signal_ci_pow[signal_ci_pow < 0] = 0.

    # print("Frac RMS of reference band:", ref.rms / ref.mean_rate)
    ## in frac rms units here -- should be few percent

    ## Broadcast signal_ref_pow into same shape as signal_ci_pow
    signal_ref_pow = np.resize(np.repeat(signal_ref_pow, meta_dict['detchans']),
            np.shape(signal_ci_pow))
    assert np.shape(signal_ref_pow) == np.shape(signal_ci_pow)

    ## Compute amplitude of noise in the cross spectrum
    temp = (noise_ci * signal_ref_pow) + (noise_ref * signal_ci_pow) + \
            (noise_ci * noise_ref)

    cs_noise_amp = np.sqrt(np.sum(temp, axis=0) / np.float(meta_dict['n_seg']))

    ## Compute amplitude of (filtered) signal in the cross spectrum
    temp1 = cs_avg * (2.0 * np.mean(meta_dict['dt']) /
            np.float(meta_dict['n_bins']))
    cs_signal_amp = np.sum(temp1, axis=0)

    ## Assume cs_noise_amp and cs_signal_amp are float arrays, size DETCHANS
    with np.errstate(all='ignore'):
        error_ratio = np.where(cs_signal_amp != 0, cs_noise_amp / \
                cs_signal_amp, 0)

    ## Take the IFFT of the cross spectrum to get the CCF
    ccf_end = fftpack.ifft(filtered_cs_avg, axis=0)

    ## Divide ccf by rms of signal in reference band
    ccf_end *= (2.0 / np.float(meta_dict['n_bins']) / ref.rms)

    ## Compute the error on the ccf
    ccf_rms_ci = np.sqrt(np.var(ccf_end, axis=0, ddof=1))
    ccf_error = ccf_rms_ci * error_ratio

    ## Re-sizing the error array to have the same shape as ccf_end
    ccf_error = np.resize(np.repeat(ccf_error, meta_dict['n_bins']),
                          (meta_dict['detchans'], meta_dict['n_bins']))
    ccf_error = ccf_error.T

    return ccf_end.real, ccf_error.real


################################################################################
def unfilt_cs_to_ccf_w_err(cs_array, meta_dict, ref):
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
        n_seconds=64, dt_mult=64, test=False, filtering=False, lo_freq=0.0,
        hi_freq=0.0, adjust=False, harmonic=False):
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

    n_seconds : int
        Number of seconds in each Fourier segment. Must be a power of 2,
        positive. [64]

    dt_mult : int
        Multiple of dt (dt is from data file) for timestep between bins. Must be
        a power of 2, positive. [64]

    bkgd_file : str
        Name of the background spectrum (in .pha format), with the same energy
        channel binning as the event list. [None]

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

    adjust : bool
        If true, the QPOs are lined up (adjusted) in frequency per obsID.
        [False]

    harmonic : bool
        If true, the filter will include the harmonic. The lo_freq - hi_freq
        range must then include the harmonic! This isn't checked. If
        filtering=False, this does nothing. [False]

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

    assert hi_freq >= lo_freq, "ERROR: Upper bound of frequency filtering must"\
            " be equal to or greater than the lower bound."

    ########################
    ## Read in data file(s)
    ########################

    if ".txt" in input_file or ".lst" in input_file or ".dat" in input_file:
        data_files = [line.strip() for line in open(input_file)]
        if not data_files:  ## If data_files is an empty list
            raise Exception("ERROR: No files in the list of event lists.")
    else:
        data_files = [input_file]

    if adjust is True:
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
            print("adjustby.txt file does not exist. NOT adjusting segments to"\
                  " line up QPO.")
            adjust_segments = np.zeros(len(data_files))
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
                 'n_bins': n_seconds * int(1.0 / (dt_mult * t_res)),
                 'detchans': detchans,
                 'filter': filtering,
                 'exposure': 0,
                 'ref_file': ref_band}

    print("\nDT = %f" % meta_dict['dt'])
    print("N_bins = %d" % meta_dict['n_bins'])
    print("Nyquist freq =", meta_dict['nyquist'])
    print("Testing?", test)
    print("Filtering?", meta_dict['filter'])
    print("Adjusting QPO?", adjust)

    ci_total = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ci')
    ref_total = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ref')
    total_seg = 0
    total_cross_spec = np.zeros((meta_dict['n_bins'], meta_dict['detchans'], 1),
            dtype=np.complex128)
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

        cross_spec_whole, ci_whole, ref_whole, n_seg, dt_whole, df_whole, \
                exposure = fits_in(in_file, meta_dict, test)

        print("Segments for this file: %d\n" % n_seg)

        total_cross_spec = np.dstack((total_cross_spec, cross_spec_whole))
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
    total_cross_spec = total_cross_spec[:,:,1:]
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
    ref_total.mean_rate /= np.float(meta_dict['n_seg'])
    avg_cross_spec = np.mean(total_cross_spec, axis=-1)

    assert ci_total.mean_rate.all() == \
           np.mean(ci_total.mean_rate_array, axis=-1).all()

    ci_total.pos_power = ci_total.power[0:meta_dict['n_bins']/2+1, :]
    ref_total.pos_power = ref_total.power[0:meta_dict['n_bins']/2+1]

    ## Compute the variance and rms of the absolute-rms-normalized reference
    ## band power spectrum
    absrms_ref_pow = raw_to_absrms(ref_total.pos_power,
            ref_total.mean_rate, meta_dict['n_bins'],
            np.mean(meta_dict['dt']), noisy=True)

    ref_total.var, ref_total.rms = var_and_rms(absrms_ref_pow,
            np.mean(meta_dict['df']))

    #############################################################
    ## Print the cross spectrum to a file, for plotting/checking
    #############################################################
    #
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

    save_for_lags(out_file, input_file, meta_dict, avg_cross_spec, ci_total,
            ref_total, lo_freq, hi_freq)

    ##############################################
    ## Compute ccf from cs, and compute error
    ##############################################

    if meta_dict['filter']:
        ccf_avg, ccf_error = filt_cs_to_ccf_w_err(avg_cross_spec, meta_dict,
                ci_total, ref_total, lo_freq, hi_freq, filter_harmonic=harmonic,
                noisy=True)
    else:
        ccf_avg, ccf_error = unfilt_cs_to_ccf_w_err(total_cross_spec,
                meta_dict, ref_total)

    print("Exposure_time = %.3f seconds" % meta_dict['exposure'])
    print("Total number of segments:", meta_dict['n_seg'])
    print("Mean rate for all of ci:", np.sum(ci_total.mean_rate))
    print("Mean rate for ci chan 6:", ci_total.mean_rate[6])
    print("Mean rate for ci chan 15:", ci_total.mean_rate[15])
    print("Mean rate for ref:", ref_total.mean_rate)

    ##########
    ## Output
    ##########

    print(avg_cross_spec[529:532,3])
    print(ccf_avg[0:4,3])

    if len(data_files) == 1:
        file_description = "Cross-correlation function of one observation"
    else:
        file_description = "Cross-correlation function of multiple observations"

    # print(ccf_avg[0,0])
    # print(ccf_avg[0,2])
    # print(ccf_avg[0,15])

    if not test and adjust and len(data_files) == 9:
        assert round(ccf_avg[0,0], 12) == 0.117937948428
        assert round(ccf_avg[0,2], 11) == 9.22641398474
        assert round(ccf_avg[0,15], 11) == 1.76422640304
        print("Passed!")
    elif test and adjust and len(data_files) == 9:
        assert round(ccf_avg[0,0], 12) == 0.106747663439
        assert round(ccf_avg[0,2], 11) == 9.56560710672
        assert round(ccf_avg[0,15], 11) == 0.88144237181
        print("Passed!")
    else:
        print("Do not have values to compare against.")

    fits_out(out_file, input_file, bkgd_file, meta_dict, ci_total.mean_rate,
            ref_total.mean_rate, float(ref_total.rms), ccf_avg, ccf_error,
            lo_freq, hi_freq, file_description)


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

    parser.add_argument('-n', '--n_seconds', type=tools.type_power_of_two,
            default=64, dest='n_seconds', help="Number of seconds in each "\
            "Fourier segment. Must be a power of 2, positive, integer. [64]")

    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two,
            default=64, dest='dt_mult', help="Multiple of dt (dt is from data "\
            "file) for timestep between bins. Must be a power of 2, positive, "\
            "integer. [64]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
            dest='test', help="Int flag: 0 if computing all segments, 1 if "\
            "only computing one segment for testing. [0]")

    parser.add_argument('-f', '--filter', default="no", dest='filter',
            help="Filtering the cross spectrum: 'no' for QPOs, or 'lofreq:"\
            "hifreq' in Hz for coherent pulsations. [no]")

    parser.add_argument('-a', '--adjust', default=1, type=int, choices={0,1},
            dest='adjust', help="Int flag: 0 if not adjusting segment length, "\
            "1 if adjusting segment length to line up peak frequency of QPOs."\
            " [1]")

    parser.add_argument('--harmonic', default=0, type=int, choices={0,1},
                        dest='filter_harmonic', help="Int flag: 0 if not "\
                        "including harmonic in filter, 1 if yes. If filter="\
                        "'no' then this makes no difference. [0]")
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
    filter_harmonic = False
    if args.filter_harmonic == 1:
        filter_harmonic = True

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

    main(args.infile, args.outfile, ref_band=args.ref_band_file,
            bkgd_file=args.bkgd_file, n_seconds=args.n_seconds,
            dt_mult=args.dt_mult, test=test, filtering=filtering,
            lo_freq=lo_freq, hi_freq=hi_freq, adjust=adjust,
            harmonic=filter_harmonic)

################################################################################
