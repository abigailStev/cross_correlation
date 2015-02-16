#!/bin/bash

################################################################################
##
## Simple script to run ccf.py, plot_ccf.py, plot_multi.py, plot_2d.py, and 
## plot_lags.py
##
## Abigail Stevens, A.L.Stevens@uva.nl, 2014-2015
##
################################################################################

home_dir=$(ls -d ~)  # the -d flag is extremely important here
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'
exe_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exe_dir/out_ccf"
# prefix="P70080"
# obsID="70080-01-01-02"
prefix="GX339-BQPO"
# obsID="95335-01-01-00"
obsID="95409-01-18-00"

red_dir="$home_dir/Reduced_data/${prefix}/$obsID"
# red_dir="$home_dur/Dropbox/Research/sample_data"
# in_file="$red_dir/GTId_eventlist_1.dat"
in_file="$red_dir/GTId_eventlist_1.fits"
# in_file="$red_dir/GTId_eventlist_1.dat"
# bkgd_spec="$home_dir/Reduced_data/$prefix/evt_bkgd_rebinned.pha"

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi

dt=1
numsec=4
testing=0  # 0 for no, 1 for yes
filtering=0 # 0 for no, 1 for yes

if (( $testing == 0 )); then
	out_file="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
elif (( $testing == 1 )); then
	out_file="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
fi

tab_ext="fits"


##################
## Running ccf.py
##################

if [ -e "$in_file" ] && [ -e "$bkgd_spec" ]; then
	time python "$exe_dir"/ccf.py "${in_file}" "${out_file}.${tab_ext}" -b "$bkgd_spec" -n "$numsec" -m "$dt" -t "$testing" -f "$filtering"
elif [ -e "$in_file" ]; then
	time python "$exe_dir"/ccf.py "${in_file}" "${out_file}.${tab_ext}" -n "$numsec" -m "$dt" -t "$testing" -f "$filtering"
else
	echo "$in_file and/or $bkgd_spec does not exist. CCF.py was not run."
fi


#################
## Plotting ccfs
#################

if [ -e "${out_file}.${tab_ext}" ]; then
	python "$exe_dir"/plot_ccf.py "${out_file}.${tab_ext}" -o "${plot_root}" -p "${prefix}/${obsID}"
# 	if [ -e "${plot_root}_chan_06.png" ]; then open -a ImageJ "${plot_root}_chan_06.png"; fi

	ccfs_plot="$exe_dir/${day}_ccfs.png"
	python "$exe_dir"/plot_multi.py "${out_file}.${tab_ext}" "$ccfs_plot" "${numsec}"
# 	if [ -e "$ccfs_plot" ]; then open -a ImageJ "$ccfs_plot"; fi
fi


####################################
## Plotting the 2D ccf with colours
####################################

# plot_file="${plot_root}_2Dccf.png"
# if [ -e "${out_file}.${tab_ext}" ]; then
# 	python "$exe_dir"/plot_2d.py "${out_file}.${tab_ext}" -o "${plot_file}"
# fi
# if [ -e "${plot_file}" ]; then open -a ImageJ "${plot_file}"; fi


#####################
## Plotting the lags
#####################

# lag_exe_dir="$home_dir/Dropbox/Research/lags"
# lag_out_dir="$lag_exe_dir/out_lags"
# cd "$lag_exe_dir"
# 
# if (( $testing == 0 )); then
# 	out_file="$lag_out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
# 	plot_root="$lag_out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
# elif (( $testing == 1 )); then
# 	out_file="$lag_out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
# 	plot_root="$lag_out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
# fi
# 
# if [ -e "${out_file}.${tab_ext}" ]; then
# # 	echo python "$lag_exe_dir"/plot_lags.py "${out_file}.${tab_ext}" -o "${plot_root}"	
# 	python "$lag_exe_dir"/plot_lags.py "${out_file}.${tab_ext}" -o "${plot_root}"	
# else
# 	echo "${out_file}.${tab_ext} does not exist. plot_lags.py was not run."
# fi


################################################################################
## All done!
################################################################################
