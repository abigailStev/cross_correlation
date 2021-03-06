#!/usr/bin/env python
"""
Compute the cross-correlation function of narrow energy channels of interest
with a broad energy reference band, using an RXTE .fits event list.

Be sure that 'tools.py' (from https://github.com/abigailStev/whizzy_scripts) is
downloaded, and it's directory is in your PYTHONPATH bash environment variable.

November 2015, Federico M. Vincentelli: Minor adjustements for application to IR
data

February 2016: Loading in a cross spectrum, collab. with Philippe Peille

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

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2014-2016"

import argparse
import numpy as np
import sys
from scipy import fftpack
from datetime import datetime
import os
import subprocess
from astropy.io import fits
from astropy.table import Table, Column
import pickle

## These are things I've written.
## Their directories are in my PYTHONPATH bash environment variable.
import tools  ## in whizzy_scripts
import ccf_lightcurves as ccf_lc  ## in cross_correlation
from plot_ccf import make_plot

EXE_DIR = "/Users/abigailstevens/Reduced_data/kHz_Philippe"
CS_FILE = EXE_DIR + \
        "/cross_spectra_Tseg1024s_Tint_1s_M0.5_SE1_4153340-415408a_cd_700_1000.p"
OUT_FILE = EXE_DIR + "/out/ccf.fits"
PLOT_FILE = EXE_DIR + "/out/ccf.eps"
REF_POW_FILE = EXE_DIR + \
        "/ref_pds_Tseg1024s_Tint_1s_M0.5_SE1_4153340-415408a_cd_700_1000.p"
BKGD_SPEC = EXE_DIR + \
        "/bksp_Tseg1024s_Tint_1s_M0.5_SE1_4153340-415408a_cd_700_1000.pha"
CI_MEAN_SPEC = EXE_DIR + \
        "/total_spectrum_Tseg1024s_Tint_1s_M0.5_SE1_4153340-415408a_cd_700_1000.pha"
N_SEC = 1
DT_MULT = 2
T_RES = 1./8192.
DETCHANS = 34
EXPOSURE = 7158.0
N_SEG = 1024

TEST = False
FILTERING = False
LO_FREQ = -1
HI_FREQ = -1

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

    # print "\nOutput sent to: %s" % out_file

    # print ccf[0,0]
    # print ccf[0,2]
    # print ccf[0,15]

    out_table = Table()
    out_table.add_column(Column(data=ccf, name='CCF'))
    out_table.add_column(Column(data=ccf_error, name='ERROR'))
    print np.shape(out_table)
    print out_table

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
    # out_table.meta['RMS_REF'] = rms_ref
    out_table.meta['NYQUIST'] = meta_dict['nyquist']
    out_table.meta['DF'] = np.mean(meta_dict['df'])
    out_table.meta['FILTER'] = str(meta_dict['filter'])
    out_table.meta['FILTFREQ'] = "%f:%f" % (lo_freq, hi_freq)
    # if meta_dict['ref_file']:
    #     out_table.meta['REF_FILE'] = meta_dict['ref_file']

    if not os.path.exists(out_file):
        os.makedirs(os.path.dirname(out_file))

    out_table.write(out_file, overwrite=True)


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
    # print "Power shape:", np.shape(power)
    # print "DT shape:", np.shape(dt)
    # print "Noise shape:", np.shape(noise)
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

    ## Computing the cross spectrum from the fourier transform
    cs_seg = np.multiply(fft_data_ci, np.conj(fft_data_ref))

    return cs_seg, ci_seg, ref_seg


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
        print_iterator = int(20)
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

    start_time = data.field('TIME')[0]
    final_time = data.field('TIME')[-1]

    ###################################
    ## Selecting PCU for interest band
    ###################################

    PCU2_mask = data.field('PCUID') == 2
    data_pcu2 = data[PCU2_mask]
    all_time_ci = np.asarray(data_pcu2.field('TIME'), dtype=np.float64)
    all_energy_ci = np.asarray(data_pcu2.field('CHANNEL'), dtype=np.float64)

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
        # print "Ref PCU =", ref_pcu
    	# refpcu_mask = data.field('PCUID') == ref_pcu

        refpcu_mask = data.field('PCUID') != 2
        data_ref = data[refpcu_mask]
        all_time_ref = np.asarray(data_ref.field('TIME'), dtype=np.float64)
        all_energy_ref = np.asarray(data_ref.field('CHANNEL'), dtype=np.float64)
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
            # print "\t", index
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

    print "Segments computed:"

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
                print IR_poisson
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
                    print "\t", n_seg

                if test is True and n_seg == 1:  # For testing
                    break

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
    # print dt_whole
    # print df_whole
    # print ref_whole.rms_array
    # print np.shape(ref_whole.rms_array)

    return cs_whole, ci_whole, ref_whole, n_seg, dt_whole, df_whole, \
            exposure


################################################################################
def get_rate_from_spec(spec_file="evt_bkgd_rebinned.pha"):
    """
    Get a count rate from a .pha spectrum file.

    Parameters
    ----------
    spec_file : str
        The full path name of the .pha file containing the energy
        spectrum, with energy channels binned in
        the same way as the data file. [evt_bkgd_rebinned.pha]

    Returns
    -------
    rate : np.array of floats
        1-D array of the count rate per energy channel for the channels of
        interest, in cts/s.

    """

    try:
        fits_hdu = fits.open(spec_file)
    except IOError:
        print "\tERROR: File does not exist: %s" % spec_file
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
    ## Getting the Fourier frequencies for the cross spectrum
    freq = fftpack.fftfreq(meta_dict['n_bins'], d=np.mean(meta_dict['dt']))
    nyq_index = meta_dict['n_bins'] / 2

    ## Only keeping the parts associated with positive Fourier frequencies
    freq = np.abs(freq[0:nyq_index])  ## because it slices at end-1, and we
            ## want to include 'nyq_index'; abs is because the nyquist freq is
            ## both pos and neg, and we want it pos here.
    # cs_avg = cs_avg[0:nyq_index + 1, :]

    out_file = out_file.replace("cross_correlation/out_ccf",
            "lag_spectra/out_lags")
    out_file = out_file.replace(".", "_cs.")
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
    print "Output sent to: %s" % out_file


    print len(freq)
    print len(cs_avg)
    out_table = Table()
    out_table.add_column(Column(data=freq, name='FREQUENCY', unit='Hz'))
    out_table.add_column(Column(data=cs_avg, name='CROSS'))
    # out_table.add_column(Column(data=ci.pos_power, name='POWER_CI'))
    # out_table.add_column(Column(data=ref.pos_power, name='POWER_REF'))

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
    out_table.write(out_file, overwrite=True)


################################################################################
def filter_freq(freq_space_array, dt, n_bins, detchans, lo_freq, hi_freq):
    """
    Apply a filter to the averaged cross-spectrum per energy channel (in
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

    ## Get the indices of the beginning and end of the signal
    min_freq_mask = freq < lo_freq  # we want the last 'True' element
    max_freq_mask = freq > hi_freq  # we want the first 'True' element
    j_min = list(min_freq_mask).index(False)
    j_max = list(max_freq_mask).index(True)

    # print "j min =", j_min
    # print "j max =", j_max

    ## Make zeroed arrays to replace with
    zero_front = np.zeros((j_min, detchans), dtype=np.complex128)
    zero_end = np.zeros((len(freq_space_array) - j_max, detchans),
            dtype=np.complex128)

    ## Concatenate the arrays together
    filt_freq_space_array = np.concatenate((zero_front,
            freq_space_array[j_min:j_max, :], zero_end), axis=0)

    ## Check that the original array is the same shape as the filtered one
    assert np.shape(freq_space_array) == np.shape(filt_freq_space_array), \
            "ERROR: Frequency-filtered cross spectrum does not have the same "\
            "size as the original cross spectrum. Something went wrong."

    return filt_freq_space_array, j_min, j_max


################################################################################
def filt_cs_to_ccf_w_err(cs_avg, meta_dict, ci, ref, lo_freq=0.0, hi_freq=0.0,
        noisy=True):
    """
    WARNING: Has not been tested lately!!

    Filter the cross-spectrum in frequency space, take the iFFT of the
    filtered cross spectrum to get the cross-correlation function, and compute
    the error on the cross-correlation function. Note that error is definitely
    NOT independent between time bins due to the filtering! But is still
    independent between energy bins.

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
            meta_dict['n_bins'], meta_dict['detchans'], lo_freq, hi_freq)

    ## Absolute rms norms of poisson noise
    noise_ci = 2.0 * ci.mean_rate
    noise_ref = 2.0 * ref.mean_rate

    ## If there's no noise in a (simulated) power spectrum, noise level = 0
    if not noisy:
        noise_ci = np.zeros(meta_dict['detchans'])
        noise_ref = 0

    noise_ref_array = np.repeat(noise_ref, meta_dict['detchans'])

    ## Extract only the signal frequencies of the mean powers
    signal_ci_pow = np.float64(ci.power[j_min:j_max, :])
    signal_ref_pow = np.float64(ref.power[j_min:j_max])

    ## Apply absolute rms2 normalization to power spectra, subtract noise
    signal_ci_pow = raw_to_absrms(signal_ci_pow, ci.mean_rate, \
            meta_dict['n_bins'], meta_dict['dt'], noisy=noisy)
    signal_ref_pow = raw_to_absrms(signal_ref_pow, ref.mean_rate, \
            meta_dict['n_bins'], meta_dict['dt'], noisy=noisy)

    print "Frac RMS of reference band:", ref.rms / ref.mean_rate
    ## in frac rms units here -- should be few percent

    ## Broadcast signal_ref_pow into same shape as signal_ci_pow
    signal_ref_pow = np.resize(np.repeat(signal_ref_pow, meta_dict['detchans']),
        np.shape(signal_ci_pow))
    assert np.shape(signal_ref_pow) == np.shape(signal_ci_pow)

    ## Compute amplitude of noise in the cross spectrum
    temp = (noise_ci * signal_ref_pow) + (noise_ref * signal_ci_pow) + \
            (noise_ci * noise_ref)
    cs_noise_amp = np.sqrt(np.sum(temp, axis=0) / np.float(meta_dict['n_seg']))

    ## Compute amplitude of signal in the cross spectrum
    temp1 = np.absolute(cs_avg[j_min:j_max, :]) * (2.0 * meta_dict['dt'] / \
            np.float(meta_dict['n_bins']))
    cs_signal_amp = np.sum(temp1, axis=0)

    ## Assume that cs_noise_amp and cs_signal_amp are float arrays, size 64
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

    return ccf_end, ccf_error


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

    # assert np.shape(cs_array) == (meta_dict['n_bins'], meta_dict['detchans'],
    #         meta_dict['n_seg'])

    ## Take the IFFT of the cross spectrum to get the CCF
    ccf_avg = fftpack.ifft(cs_array, axis=0).real
    print "CS array:", np.shape(cs_array)
    print "CCF array:", np.shape(ccf_avg)
    print ccf_avg
    ## Average across segments and normalize by rms of averaged reference band
    ## absolute-rms-normalized power spectrum
    print "ccf avg:", np.shape(ccf_avg)
    ccf_avg *= (2.0 / np.float(meta_dict['n_bins']) / ref.rms)
    print "ccf avg_2:", np.shape(ccf_avg)
    print ccf_avg

    ## Compute the standard error on each ccf bin from the segment-to-segment
    ## variations.
    # mean_ccf = np.mean(ccf_array, axis=2)
    # ccf_resid = (ccf_array.T - mean_ccf.T).T

    ## Eqn 2.3 from S. Vaughan 2013, "Scientific Inference"
    # sample_var = np.sum(ccf_resid**2, axis=2) / (meta_dict['n_seg'] - 1)

    ## Eqn 2.4 from S. Vaughan 2013, "Scientific Inference"
    # standard_error = np.sqrt(sample_var / meta_dict['n_seg'])

    standard_error = np.zeros(np.shape(ccf_avg))
    return ccf_avg, standard_error


################################################################################
def main(cs_file, out_file, ref_pow_file, bkgd_spec, ci_mean_spec, n_seconds=16,
     dt_mult=2, test=False, filtering=False, lo_freq=0.0, hi_freq=0.0):

    """

    Parameters
    ----------

    n_seconds : int
        Number of seconds in each Fourier segment. Must be a power of 2,
        positive. [16]

    dt_mult : int
        Multiple of dt (dt is from data file) for timestep between bins. Must be
        a power of 2, positive. [2]

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

    """

    #####################################################
    ## Idiot checks, to ensure that our assumptions hold
    #####################################################

    assert n_seconds > 0, "ERROR: Number of seconds per segment must be a "\
            "positive integer."
    assert dt_mult > 0, "ERROR: Multiple of dt must be a positive integer."

    assert hi_freq >= lo_freq, "ERROR: Upper bound of frequency filtering must"\
            " be equal to or greater than the lower bound."


    meta_dict = {'dt': dt_mult * T_RES,
                 't_res': T_RES,
                 'n_seconds': n_seconds,
                 'df': 1.0 / np.float(n_seconds),
                 'nyquist': 1.0 / (2.0 * dt_mult * T_RES),
                 'n_bins': n_seconds * int(1.0 / (dt_mult * T_RES)),
                 'detchans': DETCHANS,
                 'filter': filtering,
                 'exposure': 0.,
                 'ref_file': "None",
                 'n_seg': N_SEG}

    print "\nDT = %f" % meta_dict['dt']
    print "N_bins = %d" % meta_dict['n_bins']
    print "Nyquist freq =", meta_dict['nyquist']
    print "Testing?", test
    print "Filtering?", meta_dict['filter']

    ci_total = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ci')
    ref_total = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ref')

    #####################################################################
    ## Read in the background count rate from a background spectrum, and
    ## subtract from the mean count rate.
    #####################################################################

    ## For count rate of reference band, need to get the ci spectrum then
    ## remove the channel of interest then sum.

    bkgd_rate = get_rate_from_spec(bkgd_spec)
    ci_total.mean_rate = get_rate_from_spec(ci_mean_spec)
    ci_total.mean_rate -= bkgd_rate

    freq_1, cs_list = pickle.load(open(cs_file, 'r'))
    freq_2, ref_pow_list = pickle.load(open(ref_pow_file, 'r'))
    assert freq_1.all() == freq_2.all()

    # print len(cs_list)
    # print len(cs_list[0])
    # print len(ref_pow_list)
    # print len(ref_pow_list[0])

    ## ONLY DOING THIS FOR THE FIRST REFERENCE BAND AND CROSS SPECTRUM
    ref_total.power_array = np.array(ref_pow_list[-2])
    cross_spec = np.array(cs_list[-2]) * meta_dict['n_bins'] / 2.0 / meta_dict['dt']

    # print np.shape(ref_total.power_array)
    # print np.shape(cross_spec)

    ref_total.power = ref_total.power_array
    ref_total.pos_power = ref_total.power
    ref_total.mean_rate = 2500

    ## Compute the variance and rms of the absolute-rms-normalized reference
    ## band power spectrum

    ref_total.var, ref_total.rms = var_and_rms(ref_total.pos_power,
            meta_dict['df'])

    #############################################################
    ## Print the cross spectrum to a file, for plotting/checking
    #############################################################

    cs_out = np.column_stack((freq_1, cross_spec))
    np.savetxt(EXE_DIR + '/out/cs_avg.dat', cs_out)

    ####################################################################
    ## Save cross spectra and power spectra for computing lags later in
    ## lag_spectra/get_lags.py
    ####################################################################

    save_for_lags(out_file, cs_file, meta_dict, cross_spec, ci_total,
                  ref_total, lo_freq, hi_freq)

    ##############################################
    ## Compute ccf from cs, and compute error
    ##############################################

    if meta_dict['filter']:
        ccf_avg, ccf_error = filt_cs_to_ccf_w_err(cross_spec, meta_dict,
                ci_total, ref_total, lo_freq, hi_freq, noisy=True)
    else:
        ccf_avg, ccf_error = unfilt_cs_to_ccf_w_err(cross_spec,
                meta_dict, ref_total)

    print np.shape(ccf_avg)
    print np.shape(ccf_error)

    ##########
    ## Output
    ##########
    print "Exposure_time = %.3f seconds" % meta_dict['exposure']
    print "Total number of segments:", meta_dict['n_seg']
    print "Mean rate for all of ci:", np.sum(ci_total.mean_rate)
    print "Mean rate for ci chan 6:", ci_total.mean_rate[6]
    print "Mean rate for ci chan 15:", ci_total.mean_rate[15]
    print "Mean rate for ref:", ref_total.mean_rate

    file_description = "Cross-correlation function"

    fits_out(out_file, cs_file, bkgd_spec, meta_dict, ci_total.mean_rate,
            ref_total.mean_rate, float(ref_total.rms), ccf_avg, ccf_error,
            lo_freq, hi_freq, file_description)

    x_bins = np.arange(meta_dict['n_bins']/2)

    make_plot(x_bins, ccf_avg, ccf_error, meta_dict['n_bins'], " - ", PLOT_FILE, 0, 1.,
              t_length=20)
    subprocess.call(['open', PLOT_FILE])

################################################################################

## TODO: load avg cs, avg ref psd, count rate per ci: with pickle.
## don't need to check rms or variance. will need to filter the cross spectrum
## probably, or select short segments.
## mean count rate per ci will be background-subtracted

main(CS_FILE, OUT_FILE, REF_POW_FILE, BKGD_SPEC, CI_MEAN_SPEC, n_seconds=N_SEC,
     dt_mult=DT_MULT, test=TEST, filtering=FILTERING, lo_freq=LO_FREQ,
     hi_freq=HI_FREQ)