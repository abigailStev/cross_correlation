#!/usr/bin/env python
"""
Bootstraps the cross-correlation function for a set of RXTE event lists to
compute the errors on the phase-resolved energy spectra.
Use with run_multi_ccf_bootstrap.sh and ccf_bootstrap.sh. See multi_ccf.py and
ccf.py for more details.
"""
import argparse
import numpy as np
import tools  # at https://github.com/abigailStev/whizzy_scripts
import ccf as xcor
import ccf_lightcurves as ccf_lc
import datetime


__author__ = "Abigail Stevens <A.L.Stevens at uva.nl>"
__year__ = "2015"

################################################################################
def main(input_file, out_root, bkgd_file, n_seconds, dt_mult, test,
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
    print "\tStart of program:", datetime.datetime.now()

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

    if adjust is True and len(data_files) == 9:
        adjust_segments = [932, 216, 184, 570, 93, 346, 860, 533, -324]
    else:
        adjust_segments = np.zeros(len(data_files))

    #############################################
    ## Initialize; 'whole' is over one data file
    #############################################

    t_res = np.float(tools.get_key_val(data_files[0], 0, 'TIMEDEL'))
    try:
        detchans = int(tools.get_key_val(data_files[0], 0, 'DETCHANS'))
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
                 'ref_file': ""}

    print "\nDT = %f" % meta_dict['dt']
    print "N_bins = %d" % meta_dict['n_bins']
    print "Nyquist freq =", meta_dict['nyquist']
    print "Testing?", test
    print "Filtering?", meta_dict['filter']
    print "Adjusting QPO?", adjust

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

        cross_spec, ci_whole, ref_whole, n_seg, dt_whole, df_whole, exposure = \
                xcor.fits_in(in_file, meta_dict, test)

        # print "Segments for this file: %d\n" % n_seg

        total_cross_spec = np.dstack((total_cross_spec, cross_spec))
        ci_total.mean_rate_array = np.hstack((ci_total.mean_rate_array,
                ci_whole.mean_rate_array))
        ref_total.power_array = np.hstack((ref_total.power_array,
                ref_whole.power_array))
        ref_total.mean_rate_array = np.append(ref_total.mean_rate_array,
                ref_whole.mean_rate_array)
        ref_total.var_array = np.append(ref_total.var_array,
                ref_whole.var_array)
        dt_total = np.append(dt_total, dt_whole)
        df_total = np.append(df_total, df_whole)
        total_exposure += exposure
        total_seg += n_seg
    ## End of for-loop

    print ""
    meta_dict['n_seg'] = total_seg
    meta_dict['exposure'] = total_exposure
    meta_dict['dt'] = dt_total
    meta_dict['df'] = df_total
    meta_dict['adjust_seg'] = adjust_segments

    print "Mean dt:", np.mean(dt_total)
    print "Mean df:", np.mean(df_total)

    print np.shape(ci_total.mean_rate_array)

    ## Remove the first zeros from stacked arrays
    total_cross_spec = total_cross_spec[:,:,1:]
    ci_total.mean_rate_array = ci_total.mean_rate_array[:,1:]
    ref_total.power_array = ref_total.power_array[:,1:]
    ref_total.mean_rate_array = ref_total.mean_rate_array[1:]
    ref_total.var_array = ref_total.var_array[1:]
    print np.shape(ci_total.mean_rate_array)

    #####################################################################
    ## Read in the background count rate from a background spectrum, and
    ## subtract from the mean count rate.
    #####################################################################

    if bkgd_file:
        print "Using background spectrum: %s" % bkgd_file
        bkgd_rate = xcor.get_background(bkgd_file)
    else:
        print "No background spectrum."
        bkgd_rate = np.zeros(meta_dict['detchans'])

    ci_total.mean_rate -= bkgd_rate

    ## Need to use a background from ref_total. PCU for the reference band...
    # ref_total.mean_rate -= np.sum(bkgd_rate[2:26])

    ##################################################################
    ## Bootstrap the data to get errors changes which segments we use
    ## Do boot_num realizations of this
    ##################################################################

    print "\tStart bootstrapping:", datetime.datetime.now()

    print "\nBootstrap realization:"
    if boot_num >= 1:
        for b in range(1, boot_num + 1):

            if b % 5 == 0:
                print "\t%d" % b

            out_file = out_root.replace("_b-", "_b-%d" % b)

            random_segs = np.random.randint(0, total_cross_spec.shape[2], \
                    meta_dict['n_seg'])  ## Draw with replacement
            if test:
                print "Segments:", random_segs

            ref_boot = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'], \
                    detchans=meta_dict['detchans'], type='ref')
            ci_boot = ccf_lc.Lightcurve(n_bins=meta_dict['n_bins'], \
                    detchans=meta_dict['detchans'], type='ci')
            # print "\tSelect random_segs:", datetime.datetime.now()

            boot_cross_spec = total_cross_spec[:,:,random_segs]
            ci_boot.mean_rate_array = ci_total.mean_rate_array[:,random_segs]
            ref_boot.power_array = ref_total.power_array[:,random_segs]
            ref_boot.mean_rate_array = ref_total.mean_rate_array[random_segs]
            ref_boot.var_array = ref_total.var_array[random_segs]

            # ## Make sure it worked correctly
            # assert boot_cross_spec[0:3,0:3,1].all() == total_cross_spec[0:3,0:3,random_segs[1]].all(), \
            #         "ERROR: Random draw-with-replacement of segments was not "\
            #         "successful. Values of CS broke."
            # assert np.shape(boot_cross_spec) == np.shape(total_cross_spec), \
            #         "ERROR: Random draw-with-replacement of segments was not "\
            #         "successful. Arrays are not the same size. CS shape broke."
            # assert ref_boot.mean_rate_array.all() == ref_total.mean_rate_array[random_segs].all(), \
            #         "ERROR: Random draw-with-replacement of segments was not "\
            #         "successful. Ref mean rate array values broke."

            ############################
            ## Make means from segments
            ############################

            # print "\tMake means:", datetime.datetime.now()

            ref_boot.power = np.mean(ref_boot.power_array, axis=-1)
            ci_boot.mean_rate = np.mean(ci_boot.mean_rate_array, axis=-1)
            ref_boot.mean_rate = np.mean(ref_boot.mean_rate_array)

            ## Compute the variance and rms of the absolute-rms-normalized
            ## reference band power spectrum
            absrms_ref_pow = xcor.raw_to_absrms(ref_boot.power[0: \
                    meta_dict['n_bins']/2+1], ref_boot.mean_rate,
                    meta_dict['n_bins'], np.mean(meta_dict['dt']), noisy=True)

            ref_boot.var, ref_boot.rms = xcor.var_and_rms(absrms_ref_pow,
                    np.mean(meta_dict['df']))

            #####################################
            ## Compute ccf from cs and ccf error
            #####################################

            # print "\tCompute ccf and error:", datetime.datetime.now()

            ccf_avg, ccf_error = xcor.unfilt_cs_to_ccf_w_err(boot_cross_spec,
                    meta_dict, ref_boot)

            ##########
            ## Output
            ##########

            # print "\tOutput:", datetime.datetime.now()

            xcor.fits_out(out_file, input_file, bkgd_file, meta_dict,
                    ci_boot.mean_rate, ref_boot.mean_rate, float(ref_boot.rms),
                    ccf_avg, ccf_error, lo_freq, hi_freq,
                    "Bootstrapped ccf (from multiple observations)")

    else:

        ############################
        ## Make means from segments
        ############################

        ref_total.power = np.mean(ref_total.power_array, axis=-1)
        ref_total.mean_rate = np.mean(ref_total.mean_rate_array)
        ci_total.mean_rate = np.mean(ci_total.mean_rate_array, axis=-1)


        ## Compute the variance and rms of the absolute-rms-normalized reference
        ## band power spectrum
        absrms_ref_pow = xcor.raw_to_absrms(ref_total.power[0: \
                meta_dict['n_bins']/2+1], ref_total.mean_rate,
                meta_dict['n_bins'], np.mean(meta_dict['dt']), noisy=True)

        ref_total.var, ref_total.rms = xcor.var_and_rms(absrms_ref_pow,
                np.mean(meta_dict['df']))

        #####################################
        ## Compute ccf from cs and ccf error
        #####################################

        ccf_avg, ccf_error = xcor.unfilt_cs_to_ccf_w_err(total_cross_spec,
                meta_dict, ref_total)

        ##########
        ## Output
        ##########

        print ccf_avg[0,0]
        print ccf_avg[0,2]
        print ccf_avg[0,15]

        if not test and adjust:
            assert round(ccf_avg[0,0], 12) == 0.117937948428
            assert round(ccf_avg[0,2], 11) == 9.22641398474
            assert round(ccf_avg[0,15], 11) == 1.76422640304
            print "Passed!"
        elif test and adjust:
            assert round(ccf_avg[0,0], 12) == 0.106747663439
            assert round(ccf_avg[0,2], 11) == 9.56560710672
            assert round(ccf_avg[0,15], 11) == 0.88144237181
            print "Passed!"
        else:
            print "Do not have values to compare against."

        xcor.fits_out(out_root, input_file, bkgd_file, meta_dict,
                ci_total.mean_rate, ref_total.mean_rate, float(ref_total.rms),
                ccf_avg, ccf_error, lo_freq, hi_freq,
                "CCF (from multiple observations)")

    print "\tDone:", datetime.datetime.now()



################################################################################
if __name__ == "__main__":

    ##############################################
    ## Parsing input arguments and calling 'main'
    ##############################################

    parser = argparse.ArgumentParser(usage="python multi_ccf_bootstrap.py "\
            "infile outfile [OPTIONAL ARGUMENTS]", description=__doc__,
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
