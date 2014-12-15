#!/bin/bash

## Script to run multi_CCF.py and plot_CCF.py

home_dir=$(ls -d ~)  ## the -d flag is extremely important here

exe_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exe_dir/out_ccf" # for multiple input files with different obsIDs
day=$(date +%y%m%d)  ## make the date a string and assign it to 'day'
# day="141029"


if [ ! -d "$out_dir" ]; then 
	mkdir -p "$out_dir"
fi

file_list=$1
propID=$2
dt=$3
numsec=$4
testing=$5

tab_ext="fits"
plot_ext="png"

if (( $testing == 0 )); then
	out_file="$out_dir/${propID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/${propID}_${day}_t${dt}_${numsec}sec"
	saved_file_list="$out_dir/${propID}_${day}_t${dt}_${numsec}sec_filelist"
elif (( $testing == 1 )); then
	out_file="$out_dir/test_${propID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/test_${propID}_${day}_t${dt}_${numsec}sec"
	saved_file_list="$out_dir/test_${propID}_${day}_t${dt}_${numsec}sec_filelist"
fi

cp "$file_list" "$saved_file_list"

if [ -e "$saved_file_list" ]; then
	python "$exe_dir"/multi_CCF.py -i "$saved_file_list" -o "${out_file}.${tab_ext}" -n "$numsec" -m "$dt" -t "$testing"
fi

if [ -e "${out_file}.${tab_ext}" ]; then
	python "$exe_dir"/plot_CCF.py "${out_file}.${tab_ext}" -o "${plot_root}" -p "$propID"
	if [ -e "${plot_root}_chan_11.${plot_ext}" ]; then
		open -a ImageJ "${plot_root}_chan_11.${plot_ext}"
	fi
	
	ccfs_plot="$exe_dir/ccf_plot.${plot_ext}"
	python "$exe_dir"/plot_multi.py "${out_file}.${tab_ext}" "$ccfs_plot" "$numsec"
	if [ -e "$ccfs_plot" ]; then
		open -a ImageJ "$ccfs_plot"
	fi
fi

## Plotting the 2D ccf with colours

plot_file="${plot_root}_2Dccf.${plot_ext}"
if [ -e "${out_file}.${tab_ext}" ]; then
	python "$exe_dir"/plot_2d.py "${out_file}.${tab_ext}" -o "${plot_file}"
fi
if [ -e "${plot_file}" ]; then
	open -a ImageJ "${plot_file}"
fi


## Plotting the lags

lag_exe_dir="$home_dir/Dropbox/Research/lags"
lag_out_dir="$lag_exe_dir/out_lags"
cd "$lag_exe_dir"

if (( $testing == 0 )); then
	out_file="$lag_out_dir/${propID}_${day}_t${dt}_${numsec}sec"
	plot_root="$lag_out_dir/${propID}_${day}_t${dt}_${numsec}sec"
elif (( $testing == 1 )); then
	out_file="$lag_out_dir/test_${propID}_${day}_t${dt}_${numsec}sec"
	plot_root="$lag_out_dir/test_${propID}_${day}_t${dt}_${numsec}sec"
fi

if [ -e "${out_file}.${tab_ext}" ]; then
# 	echo python "$lag_exe_dir"/plot_lags.py "${out_file}.${tab_ext}" -o "${plot_root}"	
	python "$lag_exe_dir"/plot_lags.py "${out_file}.${tab_ext}" -o "${plot_root}"
else
	echo "${out_file}.${tab_ext} does not exist. plot_lags.py was not run."
fi
if [ -e "${plot_root}_lag-energy.png" ]; then
	open -a ImageJ "${plot_root}_lag-energy.png"
fi
