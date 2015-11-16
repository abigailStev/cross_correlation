#!/usr/bin/env
"""
Computes the cross-correlation function of narrow energy channels of interest
with a broad energy reference band, over multiple RXTE event lists.

Use run_multi_ccf.sh for an example.
"""
import argparse
import numpy as np
import os.path
import subprocess
from datetime import datetime
from astropy.io import fits
import tools  # at https://github.com/abigailStev/whizzy_scripts
import ccf as xcor
import ccf_lightcurves as ccf_lc

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2014-2015"

################################################################################
def dat_out(out_file, in_file_list, bkgd_file, meta_dict, mean_rate_ci_total, \
    mean_rate_ref_total, t, ccf, ccf_error, filtering):
    """
    Writes the cross-correlation function to a .dat output file.

    Parameters
    ----------
    out_file : str
        Description.

    in_file_list : str
        Description.

    out_file : str
        Description.

    bkgd_file : str
        Description.

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    mean_rate_ci_total : np.array of floats
        Description.

    mean_rate_ref_total : float
        Description.

    t : np.array of ints
        Description.

    ccf : np.array of floats
        Description.

    ccf_error : np.array of floats
        Description.

    filtering : boolean
        Description.


    Returns
    -------
    nothing
    """
    if out_file[-4:].lower() == "fits":
        out_file = out_file[:-4]+"dat"

    print "\nOutput sent to: %s" % out_file

    with open(out_file, 'w') as out:
        out.write("#\t\tCross correlation function of multiple data files")
        out.write("\n# Date(YYYY-MM-DD localtime): %s" % str(datetime.now()))
        out.write("\n# List of event lists: %s" % in_file_list)
        out.write("\n# Background spectrum: %s" % bkgd_file)
        out.write("\n# Time bin size = %.21f seconds" % meta_dict['dt'])
        out.write("\n# Number of bins per segment = %d" % meta_dict['n_bins'])
        out.write("\n# DETCHANS = %d" % meta_dict['detchans'])
        out.write("\n# Total exposure time = %d seconds" % \
            (meta_dict['n_seg'] * meta_dict['n_seconds']))
        out.write("\n# Total number of segments = %d " % meta_dict['n_seg'])
        out.write("\n# Mean count rate of ci = %s" % \
            str(list(mean_rate_ci_total)))
        out.write("\n# Mean count rate of ref band = %.5f" % \
            mean_rate_ref_total)
        out.write("\n# Filter applied in frequency domain? %s" % str(filtering))
        out.write("\n# ")
        out.write("\n# Column 1: Time bins")
        out.write("\n# Column 2-65: CCF per energy channel [count rate]")
        out.write("\n# Column 66-129: Error on ccf per energy channel [count "\
            "rate]")
        out.write("\n# ")
        for j in xrange(0, meta_dict['n_bins']):
            out.write("\n%d" % t[j])
            for i in xrange(0, meta_dict['detchans']):
                out.write("\t%.6e" % ccf[j][i].real)
            if filtering:
                for i in xrange(0, meta_dict['detchans']):
                    out.write("\t%.6e" % ccf_error[i].real)
            else:
                for i in xrange(0, meta_dict['detchans']):
                    out.write("\t%.6e" % ccf_error[j][i].real)


################################################################################
def fits_out(out_file, in_file_list, bkgd_file, meta_dict, mean_rate_ci_total,\
    mean_rate_ref_total, t, ccf, ccf_error, filtering, lo_freq, hi_freq):

    """
    Writes time bins and the cross-correlation function to a .fits output file.

    Parameters
    ----------
    out_file : string
        Description.

    in_file_list : string
        Description.

    out_file : string
        Description.

    bkgd_file : string
        Name of the background spectrum (in .pha format), with the same energy
        channel binning as the event list.

    meta_dict : dict
        Dictionary of necessary meta-parameters for data analysis.

    mean_rate_ci_total : np.array of floats
        1-D array of the mean count rate in each energy channel, across all
        data files.

    mean_rate_ref_total : float
        The mean count rate in the reference band, across all data files.

    t : np.array of ints
        1-D array of integer time bins, size = n_bins.

    ccf : np.array of floats
        2-D array of the cross correlation function. Size=(n_bins, detchans).

    ccf_error : np.array of floats
        2-D array of the error on the cross correlation function.
        Size=(n_bins, detchans).

    filtering : boolean
        If true, filters the cross spectrum in frequency space using lo_freq and
        hi_freq as boundaries. [False]

    lo_freq : float
        Low frequency bound of cross spectrum filter, in Hz.

    hi_freq : float
        High frequency bound of cross spectrum filter, in Hz.


    Returns
    -------
    nothing

    """
    print("\nOutput sent to: %s" % out_file)

    chan = np.arange(0, meta_dict['detchans'])
    energy_channels = np.tile(chan, len(t))
    if filtering:
        ccf_error = np.tile(ccf_error, len(t))
    else:
        ccf_error = ccf_error.flatten('C')
    time_bins = np.repeat(t, len(chan))
    assert len(energy_channels) == len(time_bins)

    ## Making FITS header (extension 0)
    prihdr = fits.Header()
    prihdr.set('TYPE', "Cross-correlation function of multiple data files")
    prihdr.set('DATE', str(datetime.now()), "YYYY-MM-DD localtime")
    prihdr.set('EVTLIST', in_file_list)
    prihdr.set('BKGD', bkgd_file)
    prihdr.set('DT', np.mean(meta_dict['dt']), "seconds")
    prihdr.set('DF', np.mean(meta_dict['df']), "Hz")
    prihdr.set('N_BINS', meta_dict['n_bins'], "time bins per segment")
    prihdr.set('SEC_SEG', meta_dict['n_seconds'], "seconds per segment")
    prihdr.set('NYQUIST', meta_dict['nyquist'], "Hz")
    prihdr.set('SEGMENTS', meta_dict['n_seg'], "segments, of all data")
    prihdr.set('EXPOSURE', meta_dict['exposure'], "seconds of data used")
    prihdr.set('DETCHANS', meta_dict['detchans'], "Number of detector energy "\
        "channels")
    prihdr.set('RATE_CI', str(mean_rate_ci_total.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_ref_total, "counts/second")
    prihdr.set('FILTER', str(filtering))
    prihdr.set('FILTFREQ', "%f:%f" % (lo_freq, hi_freq))
    prihdr.set('ADJUST', str(meta_dict['adjust_seg']))
    prihdu = fits.PrimaryHDU(header=prihdr)

    ## Making FITS table (extension 1)
    col1 = fits.Column(name='TIME_BIN', format='K', array=time_bins)
    col2 = fits.Column(name='CCF', unit='Counts/second', format='D',
        array=ccf.real.flatten('C'))
    col3 = fits.Column(name='ERROR', unit='', format='D',
        array=ccf_error)
    col4 = fits.Column(name='CHANNEL', unit='', format='I',
        array=energy_channels)
    cols = fits.ColDefs([col1, col2, col3, col4])
    tbhdu = fits.BinTableHDU.from_columns(cols)

    ## If the file already exists, remove it
    assert out_file[-4:].lower() == "fits", \
        'ERROR: Output file must have extension ".fits".'
    if os.path.isfile(out_file):
        subprocess.call(["rm", out_file])

    ## Writing to a FITS file
    thdulist = fits.HDUList([prihdu, tbhdu])
    thdulist.writeto(out_file)


################################################################################
def main(in_file_list, out_file, bkgd_file=None, n_seconds=64, dt_mult=64,
        test=False, filtering=False, lo_freq=0.0, hi_freq=0.0, adjust=False):
    """
    Reads in multiple event lists, splits into two light curves, makes segments
    and populates them to give them length n_bins, computes the cross spectrum
    of each segment per energy channel and then averaged cross spectrum of all
    the segments per energy channel, and then computes the cross-correlation
    function (ccf) per energy channel.

    Parameters
    ----------
    in_file_list : str
        Name of text file that contains a list of the input .fits event lists
        for analysis. Must be full path names. One file per line.

    out_file : str
        Name the FITS file to write the cross-correlation function to.

    bkgd_file : str
        Name of the background spectrum (in .pha format), with the same energy
        channel binning as the event list. [None]

    n_seconds : int
        Number of seconds in each Fourier segment. Must be a power of 2,
        positive. [64]

    dt_mult : int
        Multiple of dt (dt is from data file) for timestep between bins of light
        curve. Must be a power of 2, positive. [64]

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
        If True, adjust cross spectra to line up QPOs. Amounts to adjust each
        data file by are built-in as of 150828, for the 2010 outburst GX339-4
        Type B QPO data. [False]

    Returns
    -------
    nothing

    """

    #####################################################
    ## Idiot checks, to ensure that our assumptions hold
    #####################################################

    assert n_seconds > 0, "ERROR: n_seconds must be a positive integer."
    assert dt_mult >= 1, "ERROR: dt_mult must be a positive integer."

    ##############################################################
    ## Getting the list of input files and putting them in a list
    ##############################################################

    input_files = [line.strip() for line in open(in_file_list)]
    if not input_files:  ## If data_files is an empty list
        raise Exception("ERROR: No files in the eventlist list.")

    ###########################################################
    ## Initializations
    ## 'total' is over all data files (i.e., in multi_ccf.py)
    ## 'whole' is over one data file (i.e., in ccf.py)
    ###########################################################

    adjust_segments = [932, 216, 184, 570, 93, 346, 860, 533, -324]
    if not adjust:
        adjust_segments = [0,0,0,0,0,0,0,0,0]

    t_res = np.float64(tools.get_key_val(input_files[0], 0, 'TIMEDEL'))
    dt = dt_mult * t_res
    n_bins = n_seconds * np.int(1.0 / dt)
    try:
        detchans = np.int(tools.get_key_val(input_files[0], 0, 'DETCHANS'))
    except KeyError:
        detchans = 64

    nyq_freq = 1.0 / (2.0 * dt)
    df = 1.0 / np.float64(n_seconds)
    meta_dict = {'dt': dt, 't_res': t_res, 'n_seconds': n_seconds, \
                'df': df, 'nyquist': nyq_freq, 'n_bins': n_bins, \
                'detchans': detchans, 'filter': filtering, \
                'adjust_seg': 0, 'exposure': 0}

    ci_total = ccf_lc.Lightcurves(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ci')
    ref_total = ccf_lc.Lightcurves(n_bins=meta_dict['n_bins'],
            detchans=meta_dict['detchans'], type='ref')
    total_seg = 0
    total_cross_spec = np.zeros((meta_dict['n_bins'], meta_dict['detchans'], \
            1), dtype=np.complex128)
    # ci_total.power = np.zeros((meta_dict['n_bins'], meta_dict['detchans']), \
    #         dtype=np.float64)
    # ci_total.power_array = np.zeros((meta_dict['n_bins'], \
    #         meta_dict['detchans'], 1), dtype=np.float64)
    # ci_total.mean_rate = np.zeros(meta_dict['detchans'])
    # ci_total.mean_rate_array = np.zeros((meta_dict['detchans'], 1), \
    #         dtype=np.float64)
    # ref_total.power = np.zeros(meta_dict['n_bins'], dtype=np.float64)
    # ref_total.power_array = np.zeros((meta_dict['n_bins'], 1), \
    #         dtype=np.float64)
    # ref_total.mean_rate = 0
    # ref_total.mean_rate_array = 0
    dt_total = np.array([])
    df_total = np.array([])
    total_exposure = 0


    print "\nDT = %.15f" % meta_dict['dt']
    print "N_bins = %d" % meta_dict['n_bins']
    print "Nyquist freq = %f" % meta_dict['nyquist']
    print "Testing?", test
    print "Filtering?", filtering
    print "Adjusting QPO?", adjust

    ###################################################################
    ## Reading in the background count rate from a background spectrum
    ####################################################################

    if bkgd_file:
        print "Background spectrum: %s" % bkgd_file
        bkgd_rate = xcor.get_background(bkgd_file)
    else:
        print "No background spectrum."
        bkgd_rate = np.zeros(meta_dict['detchans'])
    print " "

    ##################################
    ## Looping through all data files
    ##################################
    i=0
    for in_file in input_files:

        if adjust:
            meta_dict['adjust_seg'] = adjust_segments[i]
        else:
            meta_dict['adjust_seg'] = 0

        cross_spec, ci_whole, ref_whole, n_seg, dt_whole, df_whole, exposure = \
                xcor.fits_in(in_file, meta_dict, test)

        print "Segments for this file: %d\n" % n_seg
        total_cross_spec = np.dstack((total_cross_spec, cross_spec))
        ci_total.power_array = np.dstack((ci_total.power_array, \
                ci_whole.power_array))
        ci_total.mean_rate_array = np.hstack((ci_total.mean_rate_array, \
                ci_whole.mean_rate_array))
        ref_total.power_array = np.hstack((ref_total.power_array, \
                ref_whole.power_array))
        ref_total.mean_rate_array = np.append(ref_total.mean_rate_array, \
                ref_whole.mean_rate_array)
        dt_total = np.append(dt_total, dt_whole)
        df_total = np.append(df_total, df_whole)
        total_exposure += exposure
        total_seg += n_seg
        i += 1
    ## End of for-loop
    print " "

    meta_dict['n_seg'] = total_seg
    meta_dict['exposure'] = total_exposure
    # print "DT array:", dt_total
    # print "df array:", df_total
    meta_dict['dt'] = dt_total
    meta_dict['df'] = df_total
    meta_dict['adjust_seg'] = adjust_segments

    print "Mean dt:", np.mean(dt_total)
    print "Mean df:", np.mean(df_total)

    ## Removing the first zeros from stacked arrays
    total_cross_spec = total_cross_spec[:,:,1:]
    ci_total.power_array = ci_total.power_array[:,:,1:]
    ci_total.mean_rate_array = ci_total.mean_rate_array[:,1:]
    ref_total.power_array = ref_total.power_array[:,1:]
    ref_total.mean_rate_array = ref_total.mean_rate_array[1:]

    ###################################################
    ## Array processing and making means from segments
    ###################################################

    avg_cross_spec, cross_spec, ci_total, ref_total, meta_dict = \
            xcor.alltogether_means(total_cross_spec, ci_total, ref_total, \
            meta_dict, bkgd_rate, False)

    # print ci_total.mean_rate[1:3]

    ######################
    ## Making lag spectra
    ######################

    xcor.save_for_lags(out_file, in_file_list, meta_dict, ci_total.mean_rate,
        ref_total.mean_rate, avg_cross_spec, ci_total.power, ref_total.power)

    ##############################################
    ## Computing ccf from cs, and computing error
    ##############################################

    if filtering:
        ccf_end, ccf_error = xcor.FILT_cs_to_ccf_w_err(avg_cross_spec, meta_dict,
            ci_total.mean_rate, ref_total.mean_rate, ci_total.power,
            ref_total.power, True, lo_freq, hi_freq)
    else:
        ccf_end = xcor.UNFILT_cs_to_ccf(avg_cross_spec, meta_dict, ref_total, \
                True)
        # print "Done with cs to ccf, going to compute error."
        ccf_error = xcor.standard_ccf_err(cross_spec, meta_dict, \
                ref_total, True)

    # print "ccf end:", ccf_end[1:3,1:3]

    print "e = ", meta_dict['n_seconds'] * meta_dict['n_seg']
    print "Exposure_time = %.3f seconds" % meta_dict['exposure']
    print "Total number of segments:", meta_dict['n_seg']
    print "Mean rate for all of ci:", np.sum(ci_total.mean_rate)
    print "Mean rate for ci chan 6:", ci_total.mean_rate[6]
    print "Mean rate for ci chan 15:", ci_total.mean_rate[15]
    print "Mean rate for ref:", ref_total.mean_rate

    t = np.arange(0, meta_dict['n_bins'])  ## gives the 'front of the bin'

    ##########
    ## Output
    ##########

    fits_out(out_file, in_file_list, bkgd_file, meta_dict, ci_total.mean_rate,\
            ref_total.mean_rate, t, ccf_end, ccf_error, filtering, lo_freq, \
            hi_freq)


################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python multi_ccf.py infile outfile"\
            " [OPTIONAL ARGUMENTS]", description=__doc__, epilog="For optional"\
            " arguments, default values are given in brackets at end of "\
            "description.")

    parser.add_argument('infile_list', help="The name of the txt "\
            "file listing the event lists to be used. One file per line. Each "\
            "event list must be .fits format containing both the reference"\
            " band and the channels of interest. Assumes channels of interest"\
            " = PCU 2, ref band = all other PCUs.")

    parser.add_argument('outfile', help="The name of the FITS file to write "\
            "the cross-correlation function to.")

    parser.add_argument('-b', '--bkgd', required=False, default=None,
            dest='bkgd_file', help="Name of the background spectrum (in .pha "\
            "format), with the same energy channel binning as the event list. "\
            "[None]")

    parser.add_argument('-n', '--n_seconds', type=tools.type_power_of_two,
            default=64, dest='n_seconds', help="Number of seconds in each "\
            "Fourier segment. Must be a power of 2, positive, integer. [64]")

    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two,
            default=64, dest='dt_mult', help="Multiple of dt (dt is from data "\
            "file) for timestep between bins. Must be a power of 2, positive, "\
            "integer. [64]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
            dest='test', help="Int flag: 0 if computing all segments, 1 if "\
            "computing only one segment for testing. [0]")

    parser.add_argument('-f', '--filter', default="no", dest='filter',
            help="Filtering the cross spectrum: 'no' for QPOs, or 'lofreq:"\
            "hifreq' in Hz for coherent pulsations. [no]")

    parser.add_argument('-a', '--adjust', action='store_true', default=False,
            dest='adjust', help="If present, adjust cross spectra to line up "\
            "QPOs. Values are built-in as of 150828. [False]")

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

    main(args.infile_list, args.outfile, bkgd_file=args.bkgd_file,
            n_seconds=args.n_seconds, dt_mult=args.dt_mult, test=test,
            filtering=filtering, lo_freq=lo_freq, hi_freq=hi_freq,
            adjust=args.adjust)

################################################################################
