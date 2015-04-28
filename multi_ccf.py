import argparse
import numpy as np
import os.path
import subprocess
from scipy import fftpack
from datetime import datetime
from astropy.io import fits
import tools  # https://github.com/abigailStev/whizzy_scripts
import ccf as xcor

__author__ = "Abigail Stevens"

"""
multi_ccf.py

Computes the cross-correlation function of narrow energy channels of interest
with a broad energy reference band, over multiple RXTE event-mode data files.

Abigail Stevens, A.L.Stevens at uva.nl, 2014-2015

"""

################################################################################
def dat_out(out_file, in_file_list, bkgd_file, param_dict, mean_rate_ci_total, \
    mean_rate_ref_total, t, ccf, ccf_error, filtering):
    """
    Writes the cross-correlation function to a .dat output file.

    Parameters
    ----------
    out_file : string
        Description.

    in_file_list : string
        Description.

    out_file : string
        Description.

    bkgd_file : string
        Description.

    param_dict : dict
        Description.

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
        out.write("\n# Time bin size = %.21f seconds" % param_dict['dt'])
        out.write("\n# Number of bins per segment = %d" % param_dict['n_bins'])
        out.write("\n# DETCHANS = %d" % param_dict['detchans'])
        out.write("\n# Total exposure time = %d seconds" % \
            (param_dict['num_seg'] * param_dict['num_seconds']))
        out.write("\n# Total number of segments = %d " % param_dict['num_seg'])
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
        for j in xrange(0, param_dict['n_bins']):
            out.write("\n%d" % t[j])
            for i in xrange(0, param_dict['detchans']):
                out.write("\t%.6e" % ccf[j][i].real)
            if filtering:
                for i in xrange(0, param_dict['detchans']):
                    out.write("\t%.6e" % ccf_error[i].real)
            else:
                for i in xrange(0, param_dict['detchans']):
                    out.write("\t%.6e" % ccf_error[j][i].real)


################################################################################
def fits_out(out_file, in_file_list, bkgd_file, param_dict, mean_rate_ci_total,\
    mean_rate_ref_total, t, ccf, ccf_error, filtering):
    """
    Writes the cross-correlation function to a .fits output file.

    Parameters
    ----------
    out_file : string
        Description.

    in_file_list : string
        Description.

    out_file : string
        Description.

    bkgd_file : string
        Description.

    param_dict : dict
        Description.

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
    print "\nOutput sent to: %s" % out_file

    chan = np.arange(0, param_dict['detchans'])
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
    prihdr.set('DT', param_dict['dt'], "seconds")
    prihdr.set('N_BINS', param_dict['n_bins'], "time bins per segment")
    prihdr.set('SEC_SEG', param_dict['num_seconds'], "seconds per segment")
    prihdr.set('SEGMENTS', param_dict['num_seg'], "segments, of all data")
    prihdr.set('EXPOSURE', param_dict['num_seg'] * param_dict['num_seconds'], \
        "seconds, of all data")
    prihdr.set('DETCHANS', param_dict['detchans'], "Number of detector energy "\
        "channels")
    prihdr.set('RATE_CI', str(mean_rate_ci_total.tolist()), "counts/second")
    prihdr.set('RATE_REF', mean_rate_ref_total, "counts/second")
    prihdr.set('FILTER', str(filtering))
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
def main(in_file_list, out_file, bkgd_file, num_seconds, dt_mult, test,
    filtering):
    """
    Reads in multiple event lists, splits into two light curves, makes segments
    and populates them to give them length n_bins, computes the cross spectrum
    of each segment per energy channel and then averaged cross spectrum of all
    the segments per energy channel, and then computes the cross-correlation
    function (ccf) per energy channel.

    Parameters
    ----------
    in_file_list : string
        Description.

    out_file : string
        Description.

    bkgd_file : string
        Description.

    num_seconds : int
        Description.

    dt_mult : int
        Description.

    test : boolean
        Description.

    filtering : boolean
        Description.


    Returns
    -------
    nothing

    """

    #####################################################
    ## Idiot checks, to ensure that our assumptions hold
    #####################################################

    assert num_seconds > 0, "ERROR: num_seconds must be a positive integer."
    assert dt_mult >= 1, "ERROR: dt_mult must be a positive integer."

    ##############################################################
    ## Getting the list of input files and putting them in a list
    ##############################################################

    input_files = [line.strip() for line in open(in_file_list)]
    if not input_files:  ## If data_files is an empty list
        raise Exception("ERROR: No files in the eventlist list.")
    # old_settings = np.seterr(divide='ignore')

    ###########################################################
    ## Initializations
    ## 'total' is over all data files (i.e., in multi_ccf.py)
    ## 'whole' is over one data file (i.e., in ccf.py)
    ###########################################################

    t_res = float(tools.get_key_val(input_files[0], 0, 'TIMEDEL'))
# 	print t_res
# 	print dt_mult
    dt = dt_mult * t_res
    n_bins = num_seconds * int(1.0 / dt)
    detchans = int(tools.get_key_val(input_files[0], 0, 'DETCHANS'))
    nyq_freq = 1.0 / (2.0 * dt)
    df = 1.0 / float(num_seconds)
    param_dict = {'dt': dt, 't_res': t_res, 'num_seconds': num_seconds, \
                 'df': df, 'nyquist': nyq_freq, 'n_bins': n_bins, \
                 'detchans': detchans}

    total_seg = 0
    cs_sum_total = np.zeros((param_dict['n_bins'], param_dict['detchans']), \
        dtype=np.complex128)
    sum_rate_ci_total = np.zeros(param_dict['detchans'])
    sum_rate_ref_total = 0
    sum_power_ci_total = np.zeros((param_dict['n_bins'], \
        param_dict['detchans']), dtype=np.float64)
    sum_power_ref_total = np.zeros(param_dict['n_bins'], dtype=np.float64)

    print "\nDT = %.15f" % param_dict['dt']
    print "N_bins = %d" % param_dict['n_bins']
    print "Nyquist freq = %f" % param_dict['nyquist']
    print "Filtering?", filtering

    ###################################################################
    ## Reading in the background count rate from a background spectrum
    ####################################################################
   
    if bkgd_file:
        print "Using background spectrum: %s" % bkgd_file
        bkgd_rate = xcor.get_background(bkgd_file)
    else:
        bkgd_rate = np.zeros(param_dict['detchans'])
    print " "

    ##################################
    ## Looping through all data files
    ##################################

    for in_file in input_files:

        cs_sum, sum_rate_ci_whole, sum_rate_ref_whole, num_seg, \
            sum_power_ci, sum_power_ref, sum_rate_ci = \
            xcor.read_and_use_segments(in_file, param_dict, test)

        print "Segments for this file: %d\n" % num_seg

        total_seg += num_seg
        cs_sum_total += cs_sum
        sum_rate_ci_total += sum_rate_ci_whole
        sum_rate_ref_total += sum_rate_ref_whole
        sum_power_ci_total += sum_power_ci
        sum_power_ref_total += sum_power_ref
        sum_rate_ci_whole = None
        sum_rate_ref_whole = None
        cs_sum = None
        num_seg = None
        sum_power_ci = None
        sum_power_ref = None

    ## End of for-loop
    print " "

    param_dict['num_seg'] = total_seg

    #########################################
    ## Turning sums over segments into means
    #########################################

    mean_ci = sum_rate_ci / float(param_dict['num_seg'])
    mean_rate_ci_total = sum_rate_ci_total / float(param_dict['num_seg'])
    mean_rate_ref_total = sum_rate_ref_total / float(param_dict['num_seg'])
    mean_power_ci = sum_power_ci_total / float(param_dict['num_seg'])
    mean_power_ref = sum_power_ref_total / float(param_dict['num_seg'])
    cs_avg = cs_sum_total / float(param_dict['num_seg'])

# 	print "CS avg:", cs_avg[1:5, 6]

    ## Printing the cross spectrum to a file, for plotting/checking
# 	cs_out = np.column_stack((fftpack.fftfreq(n_bins, d=dt), cs_avg))
# 	np.savetxt('cs_avg.dat', cs_out)

    ##################################################################
    ## Subtracting the background count rate from the mean count rate
    ##################################################################

    mean_rate_ci_total -= bkgd_rate

    ## Need to use a background from ref pcu for the reference band...
#     ref_bkgd_rate = np.mean(bkgd_rate[2:26])
#     mean_rate_ref_total -= ref_bkgd_rate    
#     print np.shape(mean_rate_ci_total)
#     print np.shape(mean_rate_ref_total)

    ######################
    ## Making lag spectra
    ######################

    xcor.save_for_lags(out_file, in_file_list, param_dict, mean_rate_ci_total,
        mean_rate_ref_total, cs_avg, mean_power_ci, mean_power_ref)

    ##############################################
    ## Computing ccf from cs, and computing error
    ##############################################

    if filtering:
        ccf_end, ccf_error = xcor.FILT_cs_to_ccf_w_err(cs_avg, param_dict,
            mean_rate_ci_total, mean_rate_ref_total, mean_power_ci,
            mean_power_ref, True)
    else:
        ccf_end, ccf_error = xcor.UNFILT_cs_to_ccf_w_err(cs_avg, param_dict,
            mean_rate_ci_total, mean_rate_ref_total, mean_power_ci,
            mean_power_ref, True)

    exposure = param_dict['num_seg'] * param_dict['num_seconds']  ## Exposure time of data used
    print "Exposure_time = %.3f seconds" % exposure
    print "Total number of segments:", param_dict['num_seg']
    print "Mean rate for all of ci:", np.sum(mean_rate_ci_total)
    print "Mean rate for ref:", mean_rate_ref_total

    t = np.arange(0, param_dict['n_bins'])  ## gives the 'front of the bin'

    ##########
    ## Output
    ##########

    fits_out(out_file, in_file_list, bkgd_file, param_dict, mean_rate_ci_total,\
        mean_rate_ref_total, t, ccf_end, ccf_error, filtering)


################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python multi_ccf.py infile outfile"\
        " [-b BKGD_SPECTRUM] [-n NUM_SECONDS] [-m DT_MULT] [-t {0,1}] [-f "\
        "{0,1}]", description="Computes the cross-correlation function of a "\
        "channel of interest with a reference band, over multiple RXTE "\
        "eventlists.", epilog="For optional arguments, default values are "\
        "given in brackets at end of description.")

    parser.add_argument('infile_list', help="The name of the ASCII "\
        "file listing the event lists to be used. One file per line. Each "\
        "event list must be .fits or .dat format containing both the reference"\
        " band and the channels of interest. Assumes channels of interest = "\
        "PCU 2, ref band = all other PCUs.")

    parser.add_argument('outfile', help="The name of the FITS file to write "\
        "the cross-correlation function to.")

    parser.add_argument('-b', '--bkgd', required=False, dest='bkgd_file',
        help="Name of the (pha/fits) background spectrum. [none]")

    parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two,
        default=1, dest='num_seconds', help="Number of seconds in each Fourier"\
        " segment. Must be a power of 2, positive, integer. [1]")

    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two,
        default=1, dest='dt_mult', help="Multiple of dt (dt is from data file)"\
        " for timestep between bins. Must be a power of 2, positive, integer. "\
        "[1]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
        dest='test', help="Int flag: 0 if computing all segments, 1 if "\
        "computing only one segment for testing. [0]")

    parser.add_argument('-f', '--filter', type=int, default=0, choices={0,1},
        dest='filter', help="Int flag: 0 if NOT applying a filter in frequency"\
        "-space, 1 if applying filter (around a coherent pulsation). [0]")

    args = parser.parse_args()

    test = False
    if args.test == 1:
        test = True

    filtering = False
    if args.filter == 1:
        filtering = True

    main(args.infile_list, args.outfile, args.bkgd_file, args.num_seconds,
        args.dt_mult, test, filtering)

################################################################################
