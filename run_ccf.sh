#!/bin/bash

################################################################################
##
## Bash script to run ccf.py, plot_ccf.py, plot_multi.py, plot_2d.py,
## get_lags.py, and covariance_spectrum.py.
##
## Don't give command line arguments. Change things in this script below.
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: HEASOFT 6.11.*, bash 3.*, and Python 2.7.* (with supporting libraries) 
## 		  must be installed in order to run this script. 
##
## Author: Abigail Stevens <A.L.Stevens at uva.nl> 2014-2015
##
################################################################################

## Checking the number of input arguments; shouldn't have any!
if (( $# != 0 )); then
    echo -e "\tDo not give command line arguments. Usage: ./run_ccf.sh\n"
    exit
fi

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	. ${HEADAS}/headas-init.sh
fi

################################################################################

## Identifying prefix (object nickname or data ID)
prefix="GX339-BQPO"
## ObsID of the data, if using just one data file (otherwise, should be ignored)
obsID="95335-01-01-01"

## Multiple of time resolution of the data for binning the light curve
dt=512
## Number of seconds per Fourier segment
numsec=64
## RXTE observation epoch of your data
obs_epoch=5
## Whether or not to run 'testing'; 0 for no, 1 for yes
testing=0
## Number of time bins to plot for the two-dimensional ccf
tlen=20
## Whether or not you want to filter the cross spectrum in frequency.
## Expected format: "no" or "low:high" (e.g. "400:402")
filtfreq="no"
## Lower frequency bound for lag-energy spectra, in Hz
lag_lf=4.0
## Upper frequency bound for lag-energy spectra, in Hz
lag_uf=7.0
## Lower energy bound for lag-frequency spectra, in detector channels
lag_le=3
## Upper energy bound for lag-frequency spectra, in detector channels
lag_ue=15
## Desired plot file extension, without the dot
p_ext="eps"
## Today's date (gets automatically), for writing in file names
day=$(date +%y%m%d)
#day="151123"
## Local file name for output files
filename="${prefix}_${day}_t${dt}_${numsec}sec_adj"
#filename="${obsID}_${day}_t${dt}_${numsec}sec"

## Your computer's home directory (gets automatically)
home_dir=$(ls -d ~)
## Directory with the ccf code. Probably $home_dir/[something]/cross_correlation
ccf_exe_dir="$home_dir/Dropbox/Research/cross_correlation"
## Directory with the ccf code. Probably $home_dir/[something]/rxte_reduce
xte_exe_dir="$home_dir/Dropbox/Research/rxte_reduce"
## Directory with the ccf code. Probably $home_dir/[something]/lag_spectra
lag_exe_dir="$home_dir/Dropbox/Research/lag_spectra"
## Directory with the reduced data
red_dir="$home_dir/Reduced_data/${prefix}/$obsID"
## A folder of lists; tells which files to use for averaging together multiple
## observations
list_dir="$home_dir/Dropbox/Lists"
## File name of the GTI'd event list in fits format (from rxte_reduce/
## good_event.sh), or
#in_file="$red_dir/GTId_eventlist.fits"
in_file="$list_dir/${prefix}_eventlists_9.lst"
## File name of the background energy spectrum for the chans of interest,
## binned to the same energy resolution as the chans of interest.
## Created in rxte_reduce/rxte_reduce_data.sh.
bkgd_spec="$red_dir/evt_bkgd_rebinned.pha"
## Energy to channel conversion table from the heasoft website.
ec_table_file="${xte_exe_dir}/e-c_table.txt"
## Table telling how the RXTE Std2 energy channels are binned together for your
## particular data mode.
chan_bin_file="$home_dir/Reduced_data/${prefix}/chan.txt"
## Energy bounds for detector channels, in keV.
## Ex: if 4 channels, should have 5 energies in this file.
## If this doesn't exist, it's created in rxte_reduce/channel_to_energy.py
energies_file="$home_dir/Reduced_data/${prefix}/energies.txt"
## Response matrix for the chans of interest data.
## Gets copied to the local lag directory later.
rsp_matrix="$home_dir/Reduced_data/${prefix}/PCU2.rsp"

################################################################################
################################################################################

## Output directory for cross_correlation and lag_spectra products.
## I recommend leaving this, hence why it's below the double hash lines.
ccf_out_dir="${ccf_exe_dir}/out_ccf/${prefix}"
lag_out_dir="${lag_exe_dir}/out_lags/${prefix}"

if [ ! -d "${ccf_out_dir}" ]; then mkdir -p "${ccf_out_dir}"; fi
if [ ! -d "${lag_out_dir}" ]; then mkdir -p "${lag_out_dir}"; fi

if (( $testing == 0 )); then
	ccf_out_file="${ccf_out_dir}/$filename.fits"
	ccf_plot_root="${ccf_out_dir}/$filename"
	lag_out_file="$lag_out_dir/$filename"
	lag_plot_root="$lag_out_dir/$filename"
elif (( $testing == 1 )); then
	ccf_out_file="${ccf_out_dir}/test_$filename.fits"
	ccf_plot_root="${ccf_out_dir}/test_$filename"
    lag_out_file="$lag_out_dir/test_$filename"
	lag_plot_root="$lag_out_dir/test_$filename"
fi

ccf_multi_plot="${ccf_out_file}_multiccfs.${p_ext}"
ccf2d_plot="${ccf_plot_root}_2Dccf.${p_ext}"
ccf2d_fits_plot="${ccf_plot_root}_2Dccf.fits"

##################
## Running ccf.py
##################

cd "${ccf_exe_dir}"

#if [ -e "$in_file" ] && [ -e "$bkgd_spec" ]; then
#	time python "${ccf_exe_dir}"/ccf.py "${in_file}" "${ccf_out_file}" \
#		    -b "$bkgd_spec" -n "$numsec" -m "$dt" -t "$testing" -f "$filtfreq"
#
#elif [ -e "${in_file}" ]; then
#	time python "${ccf_exe_dir}"/ccf.py "${in_file}" "${ccf_out_file}" \
#		    -n "$numsec" -m "$dt" -t "$testing" -f "$filtfreq"
#else
#	echo -e "\tERROR: ccf.py was not run. Eventlist and/or background energy "\
#            "spectrum doesn't exist."
#fi


#################
## Plotting ccfs
#################

if [ -e "${ccf_out_file}" ]; then

	python "${ccf_exe_dir}"/plot_ccf.py "${ccf_out_file}" -o "${ccf_plot_root}" \
		    --prefix "${prefix}/${obsID}" --ext "${p_ext}"

#	if [ -e "${ccf_out_file}_chan_15.${p_ext}" ]; then
#	    open "${ccf_out_file}_chan_15.${p_ext}"; fi

	python "${ccf_exe_dir}"/plot_multi.py "${ccf_out_file}" "$ccf_multi_plot" \
		    --prefix "${prefix}/${obsID}"

# 	if [ -e "$ccf_multi_plot" ]; then open "$ccf_multi_plot"; fi

fi

###############################################
## Getting the energy list from a channel list
###############################################

if [ ! -e "$energies_file" ]; then
	if [ -e "$ec_table_file" ] && [ -e "$chan_bin_file" ]; then

		python "${xte_exe_dir}"/channel_to_energy.py "$ec_table_file" \
			    "$chan_bin_file" "$energies_file" "$obs_epoch"

	else
		echo -e "\tERROR: channel_to_energy.py not run. ec_table_file and/or "\
                "chan_bin_file do not exist."
	fi
fi

####################################
## Plotting the 2D ccf with colours
####################################

if [ -e "${ccf_out_file}" ]; then

	python "${ccf_exe_dir}"/plot_2d.py "${ccf_out_file}" -o "${ccf2d_plot}" \
		    -p "${prefix}/${obsID}" -l "$tlen" -e "$energies_file"

#	if [ -e "${ccf2d_plot}" ]; then open "${ccf2d_plot}"; fi
fi

fits_cmd="print int(tools.get_key_val('${ccf_out_file}', 0, 'DETCHANS'))"
detchans=$(python -c "import tools; ${fits_cmd}")
tlen2=$(( 2*tlen ))

if [ -e "${ccf_out_dir}/temp.dat" ]; then
	fimgcreate bitpix=-32 \
		naxes="${tlen2},${detchans}" \
		datafile="${ccf_out_dir}/temp.dat" \
		outfile="${ccf2d_fits_plot}" \
		nskip=1 \
		history=true \
		clobber=yes
else
	echo -e "\tERROR: FIMGCREATE did not run. 2Dccf temp file does not exist."
fi

if [ -e "${ccf2d_fits_plot}" ]; then
	echo "FITS 2D ccf image: ${ccf2d_fits_plot}"
else
	echo -e "\tERROR: FIMGCREATE was not successful."
fi

#####################
## Plotting the lags
#####################

cd "${lag_exe_dir}"
local_rsp_matrix="${lag_out_dir}/${prefix}.rsp"
cp "${rsp_matrix}" "${local_rsp_matrix}"

if [ -e "${lag_out_file}_cs.fits" ]; then

	python "${lag_exe_dir}"/get_lags.py "${lag_out_file}_cs.fits" \
			"${lag_out_file}_lag.fits" "${energies_file}" \
			-o "${lag_plot_root}" --prefix "$prefix" --ext "${p_ext}" \
			--lf "${lag_lf}" --uf "${lag_uf}" --le "${lag_le}" --ue "${lag_ue}"

	if [ -e "${lag_plot_root}"_lag-energy."${p_ext}" ]; then
	    open "${lag_plot_root}"_lag-energy."${p_ext}"
	fi

	python "${lag_exe_dir}"/covariance_spectrum.py "${lag_out_file}_cs.fits" \
			"${lag_out_file}_cov.fits" --prefix "$prefix" --ext "${p_ext}" \
			--rsp "${local_rsp_matrix}" --lf "${lag_lf}" --uf "${lag_uf}"

else
	echo -e "\tERROR: get_lags.py was not run. Cross spectrum output file does"\
	        " not exist."
fi

################################################################################
## All done!
################################################################################
