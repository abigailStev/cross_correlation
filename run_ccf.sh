#!/bin/bash

## Simple script to run ccf.py, plot_ccf.py, and plot_multi.py

home_dir=$(ls -d ~)  # the -d flag is extremely important here
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'
exe_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exe_dir/out_ccf"
propID="P70080"
obsID="70080-01-01-02"

in_file="$home_dir/Reduced_data/$propID/$obsID/eventlist_1.dat"
# in_file="$home_dir/Reduced_data/$propID/$obsID/eventlist_1.fits"

# in_file="$home_dir/Dropbox/Research/sample_data/eventlist_1.dat"

if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
fi

dt=1
numsec=4
testing=1  # 0 for no, 1 for yes

if (( $testing == 0 )); then
	out_file="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
elif (( $testing == 1 )); then
	out_file="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
fi

tab_ext="fits"

if [ -e "$in_file" ]; then
	time python "$exe_dir"/ccf.py "${in_file}" "${out_file}.${tab_ext}" -n "$numsec" -m "$dt" -t "$testing"
fi

ccfs_plot="$exe_dir/${day}_ccfs.png"

if [ -e "${out_file}.${tab_ext}" ]; then
	python "$exe_dir"/plot_ccf.py "${out_file}.${tab_ext}" -o "${plot_root}" -p "${propID}/${obsID}"
	open -a ImageJ "${plot_root}_chan_06.png"
	python "$exe_dir"/plot_ccf.py "${out_file}.dat" -o "${plot_root}" -p "${propID}/${obsID}"
	open -a ImageJ "${plot_root}_chan_06.png"
	
	python "$exe_dir"/plot_multi.py "${out_file}.${tab_ext}" "$ccfs_plot" "${numsec}"
	open -a ImageJ "$ccfs_plot"
fi
