#!/usr/bin/env python

import argparse
import numpy as np
import tools  # at https://github.com/abigailStev/whizzy_scripts
import ccf as xcor
import multi_ccf as mxcor
import datetime

__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2015"

"""
Bootstraps the ccf to get errors for phase-resolved energy spectra.
Use with run_multi_ccf_bootstrap.sh and ccf_bootstrap.sh.

"""


class Lightcurves(object):
    def __init__(self, n_bins=8192, detchans=64, type='ci'):
        self.type = type

        if type.lower() == "ci":
            self.power = np.zeros((n_bins, detchans), dtype=np.float64)
            self.power_array = np.zeros((n_bins, detchans, 1), dtype=np.float64)
            self.mean_rate = np.zeros(detchans)
            self.mean_rate_array = np.zeros((detchans, 1), dtype=np.float64)
        elif type.lower() == "ref":
            self.power = np.zeros(n_bins, dtype=np.float64)
            self.power_array = np.zeros((n_bins, 1), dtype=np.float64)
            self.mean_rate = 0
            self.mean_rate_array = 0
        else:
            self.mean_rate = 0
            self.mean_rate_array = 0
            self.power = 0
            self.power_array = 0
            self.pos_power = 0



################################################################################
def main(in_file_list, out_root, bkgd_file, n_seconds, dt_mult, test,
    filtering, lo_freq, hi_freq, boot_num, adjust):
    """
    Reads in multiple event lists, splits into two light curves, makes segments
    and populates them to give them length n_bins, computes the cross spectrum
    of each segment per energy channel and then averaged cross spectrum of all
    the segments per energy channel, and then computes the cross-correlation
    function (ccf) per energy channel.

    Parameters
    ----------
    in_file_list : str
        Name of text file that contains a list of the input data files for
        analysis. Must be full path names. One file per line.

    out_root : str
        Description.

    bkgd_file : str
        Description.

    n_seconds : int
        Description.

    dt_mult : int
        Description.

    test : boolean
        Description.

    filtering : boolean
        Description.

    boot_num : int
        Number of realizations to iterate for bootstrapping.

    adjust : bool
        True if artificially adjusting QPO centroid frequency.

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

    ###################
    ## Initializations
    ###################

    if adjust:
        adjust_segments = [932, 216, 184, 570, 93, 346, 860, 533, -324]
    else:
        adjust_segments = [0, 0, 0, 0, 0, 0, 0, 0, 0]

    t_res = np.float64(tools.get_key_val(input_files[0], 0, 'TIMEDEL'))
    dt = dt_mult * t_res
    n_bins = n_seconds * np.int(1.0 / dt)
    try:
        detchans = np.int(tools.get_key_val(input_files[0], 0, 'DETCHANS'))
    except KeyError:
        detchans = 64

    nyq_freq = 1.0 / (2.0 * dt)
    df = 1.0 / np.float64(n_seconds)
    param_dict = {'dt': dt, 't_res': t_res, 'n_seconds': n_seconds, \
                'df': df, 'nyquist': nyq_freq, 'n_bins': n_bins, \
                'detchans': detchans, 'filter': filtering, 'err_bin': 200, \
                'adjust_seg': 0}

    ## 'total' is over all data files (i.e., in multi_ccf.py)
    ## 'whole' is over one data file (i.e., in ccf.py)

    ci_total = Lightcurves(n_bins=param_dict['n_bins'], \
            detchans=param_dict['detchans'], type='ci')
    ref_total = Lightcurves(n_bins=param_dict['n_bins'], \
            detchans=param_dict['detchans'], type='ref')
    total_seg = 0
    total_cross_spec = np.zeros((param_dict['n_bins'], param_dict['detchans'], \
            1), dtype=np.complex128)
    dt_total = np.array([])
    df_total = np.array([])
    total_exposure = 0

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
        print "No background spectrum."
        bkgd_rate = np.zeros(param_dict['detchans'])
    print " "

    ##################################
    ## Looping through all data files
    ##################################
    i=0
    for in_file in input_files:

        if adjust:
            param_dict['adjust_seg'] = adjust_segments[i]
        else:
            param_dict['adjust_seg'] = 0

        cross_spec, ci_whole, ref_whole, n_seg, dt_whole, df_whole, exposure \
                = xcor.fits_in(in_file, param_dict, test)

        print "Segments for this file: %d\n" % (n_seg)
        total_cross_spec = np.dstack((total_cross_spec, cross_spec))
        ci_total.power_array = np.dstack((ci_total.power_array, \
                ci_whole.power_array))
        ci_total.mean_rate_array = np.hstack((ci_total.mean_rate_array, \
                ci_whole.mean_rate_array))
        ref_total.power_array = np.hstack((ref_total.power_array, \
                ref_whole.power_array))
        ref_total.mean_rate_array = np.append(ref_total.mean_rate_array, \
                ref_whole.mean_rate_array)
        total_seg += n_seg
        dt_total = np.append(dt_total, dt_whole)
        df_total = np.append(df_total, df_whole)
        total_exposure += exposure
        ci_total.mean_rate += ci_whole.mean_rate
        ref_total.mean_rate += ref_whole.mean_rate
        ci_total.power += ci_whole.power
        ref_total.power += ref_whole.power
        i += 1
    ## End of for-loop
    # print "\tFinished cross-spec:", datetime.datetime.now()

    param_dict['n_seg'] = total_seg
    param_dict['exposure'] = total_exposure
    param_dict['adjust_seg'] = adjust_segments
    # print "DT array:", dt_total
    # print "df array:", df_total
    param_dict['dt'] = dt_total
    param_dict['df'] = df_total

    # print "Mean dt:", np.mean(dt_total)
    # print "Mean df:", np.mean(df_total)

    ## Removing the first zeros from stacked arrays, and selecting the positive
    ## Fourier frequencies in the powers.
    total_cross_spec = total_cross_spec[:,:,1:]
    ci_total.power_array = ci_total.power_array[:,:,1:]
    ci_total.mean_rate_array = ci_total.mean_rate_array[:,1:]
    ref_total.power_array = ref_total.power_array[:,1:]
    ref_total.mean_rate_array = ref_total.mean_rate_array[1:]

    ##################################################
    ## Removing the segments with a negative variance
    ##################################################

    absrms_power = xcor.raw_to_absrms(ref_total.power_array[0:param_dict['n_bins']/2+1,:], \
            ref_total.mean_rate_array, param_dict['n_bins'], param_dict['dt'], \
            True)
    print np.shape(absrms_power)

    absrms_var, absrms_rms = xcor.var_and_rms(absrms_power, param_dict['df'])
    # print var
    print "\n\nPre:", absrms_rms
    mask = np.isnan(absrms_rms)
    absrms_rms = absrms_rms[~mask]
    print "RMSs:", absrms_rms
    # print absrms_rms
    # print mask
    # print mask.shape
    # print total_cross_spec.shape
    total_cross_spec = total_cross_spec[:,:,~mask]
    # print total_cross_spec[1,1,0]
    # print np.shape(total_cross_spec)
    ci_total.power_array = ci_total.power_array[:,:,~mask]
    ci_total.mean_rate_array = ci_total.mean_rate_array[:,~mask]
    ref_total.power_array = ref_total.power_array[:,~mask]
    ref_total.mean_rate_array = ref_total.mean_rate_array[~mask]

    temp = absrms_power[:,~mask]
    print "Temp:", temp[4,:]

    if param_dict['dt'] is not np.array([]):
        for element in param_dict['dt'][mask]:
                param_dict['exposure'] -= element * param_dict['n_bins']
    param_dict['dt'] = param_dict['dt'][~mask]
    param_dict['df'] = param_dict['df'][~mask]
    param_dict['n_seg'] = param_dict['n_seg'] - np.count_nonzero(mask)

    ######################################################################
    ## Bootstrapping the data to get errors changes which segments we use
    ## Doing boot_num realizations of this
    ######################################################################

    print "\nBootstrap realization:"
    if boot_num >= 1:
        for b in range(1, boot_num+1):

            if b % 20 == 0:
                print "\t%d" % b

            out_file = out_root.replace("_b-", "_b-%d" % b)

            # print "\n\tGetting random segs:", datetime.datetime.now()
            random_segs = np.random.randint(0, total_cross_spec.shape[2], \
                    param_dict['n_seg'])  ## Draw with replacement
            if test:
                print "Segments:", random_segs

            ci_boot = Lightcurves(n_bins=param_dict['n_bins'], \
                    detchans=param_dict['detchans'], type='ci')
            ref_boot = Lightcurves(n_bins=param_dict['n_bins'], \
                    detchans=param_dict['detchans'], type='ref')
            boot_cross_spec = np.zeros((param_dict['n_bins'], param_dict['detchans'], \
                    1), dtype=np.complex128)

            # boot_cross_spec = total_cross_spec[:,:,random_segs]
            # ci_boot.power_array = ci_total.power_array[:,:,random_segs]
            # ci_boot.mean_rate_array = ci_total.mean_rate_array[:,random_segs]
            # ref_boot.power_array = ref_total.power_array[:,random_segs]
            # ## I can confirm that it's selecting random_segs correctly.
            # ref_boot.mean_rate_array = ref_total.mean_rate_array[random_segs]

            boot_var = np.array([])
            boot_rms = np.array([])

            for i in random_segs:
                boot_cross_spec = np.dstack((boot_cross_spec, \
                        total_cross_spec[:,:,i]))
                ci_boot.power_array = np.dstack((ci_boot.power_array, \
                        ci_total.power_array[:,:,i]))
                ci_boot.mean_rate_array = np.hstack((ci_boot.mean_rate_array, \
                        np.reshape(ci_total.mean_rate_array[:,i], \
                        (param_dict['detchans'],1))))
                ref_boot.power_array = np.hstack((ref_boot.power_array, \
                        np.reshape(ref_total.power_array[:,i], \
                        (param_dict['n_bins'],1))))
                ref_boot.mean_rate_array = np.vstack((ref_boot.mean_rate_array,\
                        ref_total.mean_rate_array[i]))
                absrms_pow = xcor.raw_to_absrms(ref_total.power_array[0:param_dict['n_bins']/2+1, i], \
                        ref_total.mean_rate_array[i], param_dict['n_bins'], param_dict['dt'][i], \
                        True)
                # print "\tNew Pow:", absrms_pow[4], temp[4,i], absrms_pow[4] == temp[4,i]

                var, rms = xcor.var_and_rms(absrms_pow, param_dict['df'][i])
                boot_var = np.append(boot_var, var)
                boot_rms = np.append(boot_rms, rms)

                # print "\tNew RMSs:", rms, absrms_rms[i], rms == absrms_rms[i]

            boot_cross_spec = boot_cross_spec[:,:,1:]
            ci_boot.power_array = ci_boot.power_array[:,:,1:]
            ci_boot.mean_rate_array = ci_boot.mean_rate_array[:,1:]
            ref_boot.power_array = ref_boot.power_array[:,1:]
            ref_boot.mean_rate_array = ref_boot.mean_rate_array[1:]

            # absrms_poweer = xcor.raw_to_absrms(ref_boot.power_array[0:param_dict['n_bins']/2+1,:], \
            #         np.reshape(ref_boot.mean_rate_array, (param_dict['n_seg'])),
            #         param_dict['n_bins'],
            #         np.reshape(param_dict['dt'], (param_dict['n_seg'])))
            # print "Pow:", absrms_poweer[4,:]

            # variance, aremes = xcor.var_and_rms(absrms_poweer, param_dict['df'][i])

            # print aremes




            # for i in range(param_dict['n_seg']):
            #     absrms_pow = xcor.raw_to_absrms(ref_boot.power_array[0:param_dict['n_bins']/2+1, i], \
            #             ref_boot.mean_rate_array[i], param_dict['n_bins'], param_dict['dt'][i], \
            #             True)
            #     print "\tNew Pow:", absrms_pow[4], temp[4,i], absrms_pow[4] == temp[4,i]
            #
            #     var, rms = xcor.var_and_rms(absrms_pow, param_dict['df'][i])
            #
            #
            #     print "\tNew RMSs:", rms, absrms_rms[i], rms == absrms_rms[i]



            ## Making sure it worked correctly
            assert boot_cross_spec[0:3,0:3,1].all() == total_cross_spec[0:3,0:3,random_segs[1]].all(), \
                    "ERROR: Random draw-with-replacement of segments was not "\
                    "successful. Values of CS broke."
            assert np.shape(boot_cross_spec) == np.shape(total_cross_spec), \
                    "ERROR: Random draw-with-replacement of segments was not "\
                    "successful. Arrays are not the same size. CS shape broke."
            assert ref_boot.mean_rate_array.all() == ref_total.mean_rate_array[random_segs].all(), \
                    "ERROR: Random draw-with-replacement of segments was not "\
                    "successful. Ref mean rate array values broke."

            # exit()

            ##############################
            ## Making means from segments
            ##############################

            avg_boot_cross_spec, cross_spec, ci_boot, ref_boot, param_dict = \
                    xcor.alltogether_means(boot_cross_spec, ci_boot, \
                    ref_boot, param_dict, bkgd_rate, True)
            # print "\tMade means from segments:", datetime.datetime.now()

            ######################
            ## Making lag spectra
            ######################

            xcor.save_for_lags(out_file, in_file_list, param_dict, \
                    ci_boot.mean_rate, ref_boot.mean_rate, boot_cross_spec, \
                    ci_boot.power, ref_boot.power)

            ##############################################
            ## Computing ccf from cs, and computing error
            ##############################################

            if filtering:
                ccf_end, ccf_error = xcor.FILT_cs_to_ccf_w_err(avg_boot_cross_spec, \
                        param_dict, ci_boot.mean_rate, ref_boot.mean_rate, \
                        ci_boot.power, ref_boot.power, True, lo_freq, hi_freq)
            else:
                ccf_end = xcor.UNFILT_cs_to_ccf(avg_boot_cross_spec, param_dict, \
                        ref_boot, True)

                ccf_error = xcor.standard_ccf_err(boot_cross_spec, param_dict, \
                        ref_boot, noisy=True, absrms_var=boot_var, \
                        absrms_rms=boot_rms)

            # print ccf_error
            # print "ccf end:", ccf_end[1:3,1:3]
            # print "\tGoing to output:", datetime.datetime.now()
            # print "Exposure_time = %.3f seconds" % exposure
            # print "Total number of segments:", param_dict['n_seg']
            # print "Mean rate for all of ci:", np.sum(ci_total.mean_rate)
            # print "Mean rate for ci chan 6:", ci_total.mean_rate[6]
            # print "Mean rate for ci chan 15:", ci_total.mean_rate[15]
            # print "Mean rate for ref:", ref_total.mean_rate

            ##########
            ## Output
            ##########

            t = np.arange(0, param_dict['n_bins'])  ## gives the 'front of the bin'

            mxcor.fits_out(out_file, in_file_list, bkgd_file, param_dict, \
                        ci_boot.mean_rate, ref_boot.mean_rate, t, ccf_end, \
                        ccf_error, filtering, lo_freq, hi_freq)

    else:

        ##############################
        ## Making means from segments
        ##############################

        avg_cross_spec, cross_spec, ci_total, ref_total, param_dict = \
                xcor.alltogether_means(total_cross_spec, ci_total, \
                ref_total, param_dict, bkgd_rate, True)

        ##############################################
        ## Computing ccf from cs, and computing error
        ##############################################

        if filtering:
            ccf_end, ccf_error = xcor.FILT_cs_to_ccf_w_err(avg_cross_spec, \
                    param_dict, ci_total.mean_rate, ref_total.mean_rate, \
                    ci_total.power, ref_total.power, True, lo_freq, hi_freq)
        else:
            ccf_end = xcor.UNFILT_cs_to_ccf(avg_cross_spec, param_dict, \
                    ref_total, True)

            ccf_error = xcor.standard_ccf_err(total_cross_spec, param_dict, \
                    ref_total, True)

        # print "ccf end:", ccf_end[1:3,1:3]
        # print "Exposure_time = %.3f seconds" % param_dict['exposure']
        # print "Total number of segments:", param_dict['n_seg']
        # print "Mean rate for all of ci:", np.sum(ci_total.mean_rate)
        # print "Mean rate for ci chan 6:", ci_total.mean_rate[6]
        # print "Mean rate for ci chan 15:", ci_total.mean_rate[15]
        # print "Mean rate for ref:", ref_total.mean_rate

        ##########
        ## Output
        ##########

        out_file = out_root.replace("_b-", "")
        t = np.arange(0, param_dict['n_bins'])  ## gives the 'front of the bin'

        mxcor.fits_out(out_file, in_file_list, bkgd_file, param_dict, \
                    ci_total.mean_rate, ref_total.mean_rate, t, ccf_end, \
                    ccf_error, filtering, lo_freq, hi_freq)




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

    parser.add_argument('-n', '--n_seconds', type=tools.type_power_of_two,
            default=1, dest='n_seconds', help="Number of seconds in each "\
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
            "bootstrapping. Must be a positive integer. If boot_num=0, doesn't"\
            " bootstrap and just does regular CCF. [1]")

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

    main(args.infile_list, args.outfile, args.bkgd_file, args.n_seconds,
            args.dt_mult, test, filtering, lo_freq, hi_freq, args.boot_num, \
            args.adjust)

################################################################################