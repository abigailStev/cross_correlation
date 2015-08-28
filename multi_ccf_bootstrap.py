import argparse
import numpy as np
import os.path
import subprocess
from scipy import fftpack
from datetime import datetime
from astropy.io import fits
import tools  # at https://github.com/abigailStev/whizzy_scripts
import ccf as xcor
import multi_ccf as mxcor

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"

"""
Bootstraps the ccf to get errors for phase-resolved energy spectra.
Use with run_multi_ccf_bootstrap.sh and ccf_bootstrap.sh.

2015

"""


################################################################################
def main(in_file_list, out_file, bkgd_file, num_seconds, dt_mult, test,
    filtering, lo_freq, hi_freq, boot_num, adjust):
    """
    Reads in multiple event lists, splits into two light curves, makes segments
    and populates them to give them length n_bins, computes the cross spectrum
    of each segment per energy channel and then averaged cross spectrum of all
    the segments per energy channel, and then computes the cross-correlation
    function (ccf) per energy channel.

    Parameters
    ----------
    in_file_list : string
        Name of text file that contains a list of the input data files for
        analysis. Must be full path names. One file per line.

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

    boot_num : int
        Number of realizations to iterate for bootstrapping.

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

    ###########################################################
    ## Initializations
    ## 'total' is over all data files (i.e., in multi_ccf.py)
    ## 'whole' is over one data file (i.e., in ccf.py)
    ###########################################################

    adjust_segs = [932, 216, 184, 570, 93, 346, 860, 533, -324]

    t_res = np.float64(tools.get_key_val(input_files[0], 0, 'TIMEDEL'))
    dt = dt_mult * t_res
    n_bins = num_seconds * np.int(1.0 / dt)
    try:
        detchans = np.int(tools.get_key_val(input_files[0], 0, 'DETCHANS'))
    except KeyError:
        detchans = 64

    nyq_freq = 1.0 / (2.0 * dt)
    df = 1.0 / np.float64(num_seconds)
    param_dict = {'dt': dt, 't_res': t_res, 'num_seconds': num_seconds, \
                'df': df, 'nyquist': nyq_freq, 'n_bins': n_bins, \
                'detchans': detchans, 'filter': filtering, 'err_bin': 200, \
                'adjust_seg': 0}

    ci_total = xcor.Lightcurve()
    ref_total = xcor.Lightcurve()
    total_seg = 0
    total_cross_spec = np.zeros((param_dict['n_bins'], param_dict['detchans'], \
            1), dtype=np.complex128)
    ci_total.power = np.zeros((param_dict['n_bins'], param_dict['detchans']), \
            dtype=np.float64)
    ci_total.power_array = np.zeros((param_dict['n_bins'], \
            param_dict['detchans'], 1), dtype=np.float64)
    ci_total.mean_rate = np.zeros(param_dict['detchans'])
    ci_total.mean_rate_array = np.zeros((param_dict['detchans'], 1), \
            dtype=np.float64)
    ref_total.power = np.zeros(param_dict['n_bins'], dtype=np.float64)
    ref_total.power_array = np.zeros((param_dict['n_bins'], 1), \
            dtype=np.float64)
    ref_total.mean_rate = 0
    ref_total.mean_rate_array = 0

    print "\nDT = %.15f" % param_dict['dt']
    print "N_bins = %d" % param_dict['n_bins']
    print "Nyquist freq = %f" % param_dict['nyquist']
    print "Filtering?", filtering
    print "Adjusting QPO?", adjust

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
    i=0
    for in_file in input_files:

        if adjust:
            param_dict['adjust_seg'] = adjust_segs[i]
        else:
            param_dict['adjust_seg'] = 0

        cross_spec, ci_whole, ref_whole, num_seg  = xcor.fits_in(in_file, \
                param_dict, test)

        print "Segments for this file: %d\n" % num_seg
        total_cross_spec = np.dstack((total_cross_spec, cross_spec))
        ci_total.power_array = np.dstack((ci_total.power_array, \
                ci_whole.power_array))
        ci_total.mean_rate_array = np.hstack((ci_total.mean_rate_array, \
                ci_whole.mean_rate_array))
        ref_total.power_array = np.hstack((ref_total.power_array, \
                ref_whole.power_array))
        ref_total.mean_rate_array = np.append(ref_total.mean_rate_array, \
                ref_whole.mean_rate_array)
        total_seg += num_seg
        ci_total.mean_rate += ci_whole.mean_rate
        ref_total.mean_rate += ref_whole.mean_rate
        ci_total.power += ci_whole.power
        ref_total.power += ref_whole.power
        i += 1
    ## End of for-loop
    print " "

    param_dict['num_seg'] = total_seg

    #########################################
    ## Turning sums over segments into means
    #########################################

    ## Removing the first zeros from stacked arrays
    total_cross_spec = total_cross_spec[:,:,1:]
    ci_power_array = ci_total.power_array[:,:,1:]
    ci_mean_rate_array = ci_total.mean_rate_array[:,1:]
    ref_power_array = ref_total.power_array[:,1:]
    ref_mean_rate_array = ref_total.mean_rate_array[1:]

    print np.shape(ci_power_array)

    # # print "Shape mean rate array:", np.shape(ref_mean_rate_array)
    # # print "Shape power array:", np.shape(ref_power_array)
    # absrms_power = np.asarray([xcor.raw_to_absrms(ref_power_array[:,i], \
    #         ref_mean_rate_array[i], param_dict['n_bins'], param_dict['dt'], \
    #         True) for i in range(param_dict['num_seg'])])
    # ## Note that here, the axes are weird, so it's size (num_seg, n_bins)
    #
    # # print "Shape of absrms_power:", np.shape(absrms_power.T)
    # absrms_var, absrms_rms = xcor.var_and_rms(absrms_power.T, param_dict['df'])
    #
    # mask = np.isnan(absrms_rms)
    #
    # total_cross_spec = total_cross_spec[:,:,~mask]
    # ci_power_array = ci_power_array[:,:,~mask]
    # ci_mean_rate_array = ci_mean_rate_array[:,~mask]
    # ref_power_array = ref_power_array[:,~mask]
    # ref_mean_rate_array = ref_mean_rate_array[~mask]
    # print np.shape(absrms_rms)
    # ref_rms = absrms_rms[~mask]
    #
    # print np.shape(ci_power_array)
    # print np.shape(total_cross_spec)
    # print np.shape(ci_mean_rate_array)
    # print np.shape(ref_power_array)
    # print len(ref_mean_rate_array)
    #
    # param_dict['num_seg'] = param_dict['num_seg'] - np.count_nonzero(mask)

    ######################################################################
    ## Bootstrapping the data to get errors changes which segments we use
    ## Doing boot_num realizations of this
    ######################################################################
    if boot_num >= 1:
        for b in range(1, boot_num+1):

            # random_segs = np.random.randint(0, total_cross_spec.shape[2], \
            #         param_dict['num_seg'])  ## Draw with replacement
            # if test:
            #     print random_segs
            #
            # random_cross_spec = total_cross_spec[:,:,random_segs]
            # ci_total.power_array = ci_power_array[:,:,random_segs]
            # ci_total.mean_rate_array = ci_mean_rate_array[:,random_segs]
            # ref_total.power_array = ref_power_array[:,random_segs]
            # ref_total.mean_rate_array = ref_mean_rate_array[random_segs]



            random_cross_spec = total_cross_spec
            ci_total.power_array = ci_power_array
            ci_total.mean_rate_array = ci_mean_rate_array
            ref_total.power_array = ref_power_array
            ref_total.mean_rate_array = ref_mean_rate_array

            ## Printing for tests and checks
            # print "Total cs:", total_cross_spec[1,4,:]
            # print "Random cs:", random_cross_spec[1,4,:]
            # print "First slice of random:", random_cross_spec[0:3,0:3,1]
            # print "Corresponding slice of total:", total_cross_spec[0:3,0:3,random_segs[1]]

            ## Making sure it worked correctly
            # assert random_cross_spec[0:3,0:3,1].all() == total_cross_spec[0:3,0:3,random_segs[1]].all(), \
            #         "ERROR: Random draw-with-replacement of segments was not "\
            #         "successful."
            assert np.shape(random_cross_spec) == np.shape(total_cross_spec), \
                    "ERROR: Random draw-with-replacement of segments was not "\
                    "successful. Arrays are not the same size."

            ## Making means
            avg_cross_spec = np.mean(random_cross_spec, axis=2)
            # avg_ccf = np.mean(random_ccf, axis=2)
            # ref_rms = np.mean(ref_rms)
            ci_total.power = np.mean(ci_total.power_array, axis=2)
            ci_total.mean_rate = np.mean(ci_total.mean_rate_array, axis=1)
            ref_total.power = np.mean(ref_total.power_array, axis=1)
            ref_total.mean_rate = np.mean(ref_total.mean_rate_array)

            print np.shape(ci_total.mean_rate_array)
            print ci_total.mean_rate_array[1:4,:]
            print np.shape(ci_total.mean_rate)
            print ci_total.mean_rate[1:4]

            # ## Printing the cross spectrum to a file, for plotting/checking
            # cs_out = np.column_stack((fftpack.fftfreq(n_bins, d=dt), \
            #       avg_cross_spec.real))
            # np.savetxt('cs_avg_adj.dat', cs_out)


            ##################################################################
            ## Subtracting the background count rate from the mean count rate
            ##################################################################

            ci_total.mean_rate -= bkgd_rate

            ## Need to use a background from ref pcu for the reference band...
            # ref_total.mean_rate -= np.mean(bkgd_rate[2:26])

            ######################
            ## Making lag spectra
            ######################

            # xcor.save_for_lags(out_file, in_file_list, param_dict, \
            #         ci_total.mean_rate, ref_total.mean_rate, avg_cross_spec, \
            #         ci_total.power, ref_total.power)

            ##############################################
            ## Computing ccf from cs, and computing error
            ##############################################

            # avg_ccf /= ref_rms

            if filtering:
                ccf_end, ccf_error = xcor.FILT_cs_to_ccf_w_err(avg_cross_spec, \
                        param_dict, ci_total.mean_rate, ref_total.mean_rate, \
                        ci_total.power, ref_total.power, True, lo_freq, hi_freq)
            else:
                ccf_end = xcor.UNFILT_cs_to_ccf(avg_cross_spec, param_dict, \
                        ref_total, True)

                ccf_error = xcor.standard_ccf_err(total_cross_spec, param_dict, \
                        ref_total, True)

            # print "Avg ccf:", avg_ccf[1:3,1:3]
            print "ccf end:", ccf_end[1:3,1:3]

            exposure = param_dict['num_seg'] * param_dict['num_seconds']  ## Exposure time of data used
            print "Exposure_time = %.3f seconds" % exposure
            print "Total number of segments:", param_dict['num_seg']
            # print "Mean rate for all of ci:", np.sum(ci_total.mean_rate)
            # print "Mean rate for ci chan 6:", ci_total.mean_rate[6]
            # print "Mean rate for ci chan 15:", ci_total.mean_rate[15]
            print "Mean rate for ref:", ref_total.mean_rate

            t = np.arange(0, param_dict['n_bins'])  ## gives the 'front of the bin'

            ##########
            ## Output
            ##########

            out_file = out_file.replace("boot", "boot-%d" % b)
            print out_file

            mxcor.fits_out(out_file, in_file_list, bkgd_file, param_dict, \
                    ci_total.mean_rate, ref_total.mean_rate, t, ccf_end, \
                    ccf_error, filtering, lo_freq, hi_freq, adjust)
    else:

        pass

################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python multi_ccf_bootstrap.py "\
            "infile outfile [-b BKGD_SPECTRUM] [-n NUM_SECONDS] [-m DT_MULT] "\
            "[-t {0,1}] [-f {0,1}] [--bootstrap BOOT_NUM]",
            description="Bootstraps the cross-correlation function for a set "\
            "of RXTE eventlists to compute the error on the phase-resolved "\
            "energy spectra. See multi_ccf.py and ccf.py for more details.",
            epilog="For optional arguments, default values are given in "\
            "brackets at end of description.")

    parser.add_argument('infile_list', help="The name of the ASCII "\
            "file listing the event lists to be used. One file per line. Each "\
            "event list must be .fits or .dat format containing both the "\
            "reference band and the channels of interest. Assumes channels of "\
            "interest = PCU 2, ref band = all other PCUs.")

    parser.add_argument('outfile', help="The name of the FITS file to write "\
            "the cross-correlation function to.")

    parser.add_argument('-b', '--bkgd', required=False, dest='bkgd_file',
            help="Name of the (pha/fits) background spectrum. [none]")

    parser.add_argument('-n', '--num_seconds', type=tools.type_power_of_two,
            default=1, dest='num_seconds', help="Number of seconds in each "\
            "Fourier segment. Must be a power of 2, positive, integer. [1]")

    parser.add_argument('-m', '--dt_mult', type=tools.type_power_of_two,
            default=1, dest='dt_mult', help="Multiple of dt (dt is from data "\
            "file) for timestep between bins. Must be a power of 2, positive,"\
            " integer. [1]")

    parser.add_argument('-t', '--test', type=int, default=0, choices={0,1},
            dest='test', help="Int flag: 0 if computing all segments, 1 if "\
            "computing only one segment for testing. [0]")

    parser.add_argument('-f', '--filter', default="no", dest='filter',
            help="Filtering the cross spectrum: 'no' for QPOs, or 'lofreq:"\
            "hifreq' in Hz for coherent pulsations. [no]")

    parser.add_argument('--bootstrap', dest='boot_num', default=1,
            type=tools.type_positive_int, help="Number of realizations for "\
            "bootstrapping. Must be a positive integer. [1]")

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

    main(args.infile_list, args.outfile, args.bkgd_file, args.num_seconds,
        args.dt_mult, test, filtering, lo_freq, hi_freq, args.boot_num, args.adjust)

################################################################################
