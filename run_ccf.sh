# !/bin/bash

## Simple script to run ccf.py

home_dir=$(ls -d ~)  # the -d flag is extremely important here
day=$(date +%y%m%d)  # make the date a string and assign it to 'day', for filename
exec_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exec_dir/out_ccf"
propID="P70080"
obsID="70080-01-01-02"
in_file="$home_dir/Reduced_data/$propID/$obsID/eventlist_1.dat"
# in_file="$home_dir/Dropbox/Research/sample_data/eventlist_1.dat"

if test ! -d "$out_dir"
	then mkdir -p "$out_dir"
fi

dt=1
numsec=4
testing=0

if (( $testing == 0 )); then
	out_file="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
elif (( $testing == 1 )); then
	out_file="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
fi

tab_ext="dat"

if [ -e "$in_file" ]; then
	time python "$exec_dir"/ccf.py "${in_file}" "${out_file}.${tab_ext}" "$numsec" "$dt" "$testing"
fi

if [ -e "${out_file}.${tab_ext}" ]; then
	python "$exec_dir"/plot_ccf.py "${out_file}.${tab_ext}" "${plot_root}" "${propID}/${obsID}"
	open -a ImageJ "${plot_root}_chan_06.png"
fi

echo ""
