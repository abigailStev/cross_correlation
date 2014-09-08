#!/bin/bash

## Script to run multi_CCF.py and plot_CCF.py

home_dir=$(ls -d ~)  ## the -d flag is extremely important here

exe_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exe_dir/out_ccf" # for multiple input files with different obsIDs
day=$(date +%y%m%d)  ## make the date a string and assign it to 'day'
# day="140718"

if [ ! -d "$out_dir" ]; then 
	mkdir -p "$out_dir"
fi

file_list=$1
propID=$2
dt=$3
numsec=$4
testing=$5

tab_ext="dat"
plot_ext="png"

if (( $testing == 0 )); then
	out_file="$out_dir/${propID}_${day}_t${dt}_${numsec}sec"
	plot_file="$out_dir/${propID}_${day}_t${dt}_${numsec}sec"
	saved_file_list="$out_dir/${propID}_${day}_t${dt}_${numsec}sec_filelist"
elif (( $testing == 1 )); then
	out_file="$out_dir/test_${propID}_${day}_t${dt}_${numsec}sec"
	plot_file="$out_dir/test_${propID}_${day}_t${dt}_${numsec}sec"
	saved_file_list="$out_dir/test_${propID}_${day}_t${dt}_${numsec}sec_filelist"
fi

cp "$file_list" "$saved_file_list"
ccfs_plot="$exe_dir/ccf_plot.png"

python "$exe_dir"/multi_CCF.py -i "$saved_file_list" -o "${out_file}.${tab_ext}" -n "$numsec" -m "$dt" -t "$testing"

# if [ -e "${out_file}.${tab_ext}" ]; then
# 	python "$exe_dir"/plot_CCF.py -i "${out_file}.${tab_ext}" -o "${plot_file}" -p "$propID"
# 	open -a ImageJ "${plot_file}_chan_06.${plot_ext}"
# 	
# 	python "$exe_dir"/plot_multi.py "${out_file}.${tab_ext}" "$ccfs_plot" "$numsec"
# 	open -a ImageJ "$ccfs_plot"
# fi
