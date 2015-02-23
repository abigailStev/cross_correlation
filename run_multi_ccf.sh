#!/bin/bash

################################################################################
## 
## Bash script to run multi_CCF.py, plot_CCF.py, plot_multi.py, plot_2d.py, and 
## plot_lags.py
##
## Example call: ./run_multi_ccf.sh ./eventlists.lst J1808 4 16 0 150131
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: HEASOFT 6.11.*, bash 3.*, and Python 2.7.* (with supporting libraries) 
## 		  must be installed in order to run this script. 
##
## Written by Abigail Stevens, A.L.Stevens@uva.nl, 2014-2015
##
################################################################################

## Checking the number of input arguments
if (( $# != 6 )); then
    echo -e "\tUsage: ./run_multi_ccf.sh <file list> <prefix> <dt multiple> <num seconds> <testing> <date>\n"
    exit
fi

file_list=$1
prefix=$2
dt=$3
numsec=$4
testing=$5
day=$6

################################################################################

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	. $HEADAS/headas-init.sh
fi

home_dir=$(ls -d ~)

exe_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exe_dir/out_ccf" # for multiple input files with different obsIDs
lag_exe_dir="$home_dir/Dropbox/Research/lags"
lag_out_dir="$lag_exe_dir/out_lags"

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi
if [ ! -d "$lag_out_dir" ]; then mkdir -p "$lag_out_dir"; fi

filtering=0  ## 0 = no, 1 = yes; 0 is for QPOs, 1 is for coherent pulses
bkgd_spec="$home_dir/Reduced_data/$prefix/evt_bkgd_rebinned.pha"

lag_lf=4  ## Lower frequency bound for lag spectra, in Hz
lag_uf=7  ## Upper frequency bound for lag spectra, in Hz

tab_ext="fits"
plot_ext="png"

################################################################################
################################################################################

if (( $testing == 0 )); then
	out_file="$out_dir/${prefix}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/${prefix}_${day}_t${dt}_${numsec}sec"
	saved_file_list="$out_dir/${prefix}_${day}_t${dt}_${numsec}sec_filelist"
elif (( $testing == 1 )); then
	out_file="$out_dir/test_${prefix}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/test_${prefix}_${day}_t${dt}_${numsec}sec"
	saved_file_list="$out_dir/test_${prefix}_${day}_t${dt}_${numsec}sec_filelist"
fi

cp "$file_list" "$saved_file_list"

########################
## Running multi_ccf.py
########################

if [ -e "$saved_file_list" ] && [ -e "$bkgd_spec" ]; then

	python "$exe_dir"/multi_CCF.py "$saved_file_list" "${out_file}.${tab_ext}" \
		-b "$bkgd_spec" -n "$numsec" -m "$dt" -t "$testing" -f "$filtering"
		
elif [ -e "$saved_file_list" ]; then

	python "$exe_dir"/multi_CCF.py "$saved_file_list" "${out_file}.${tab_ext}" \
		-n "$numsec" -m "$dt" -t "$testing" -f "$filtering"

else 
	echo -e "\tERROR: multi_ccf.py was not run. List of eventlists and/or background energy spectrum doesn't exist."
fi

####################
## Plotting the ccf
####################

if [ -e "${out_file}.${tab_ext}" ]; then
	python "$exe_dir"/plot_CCF.py "${out_file}.${tab_ext}" -o "${plot_root}" \
		-p "$prefix"
# 	if [ -e "${plot_root}_chan_06.${plot_ext}" ]; then open -a ImageJ "${plot_root}_chan_06.${plot_ext}"; fi
	
	ccfs_plot="$exe_dir/ccf_plot.${plot_ext}"
	python "$exe_dir"/plot_multi.py "${out_file}.${tab_ext}" "$ccfs_plot" \
		"$numsec"
# 	if [ -e "$ccfs_plot" ]; then open -a ImageJ "$ccfs_plot"; fi
fi

#######################
## Plotting the 2D ccf
#######################

plot_file="${plot_root}_2Dccf.${plot_ext}"
if [ -e "${out_file}.${tab_ext}" ]; then
	python "$exe_dir"/plot_2d.py "${out_file}.${tab_ext}" -o "${plot_file}"
# 	if [ -e "${plot_file}" ]; then open -a ImageJ "${plot_file}"; fi
fi
	
plot_file="${plot_root}_2Dccf.fits"
detchans=$(python -c "from tools import get_key_val; print get_key_val('${out_file}.fits', 0, 'DETCHANS')")
echo "$detchans"

# if [ -e "$out_dir/temp.dat" ]; then
# 	fimgcreate bitpix=-32 \
# 		naxes="70,${detchans}" \
# 		datafile="$out_dir/temp.dat" \
# 		outfile="${plot_root}_2Dccf.fits" \
# 		nskip=1 \
# 		history=true \
# 		clobber=yes
# else
# 	echo -e "\tERROR: FIMGCREATE did not run. 2Dccf temp file does not exist."
# fi

if [ -e "${plot_root}_2Dccf.fits" ]; then
	echo "FITS 2D ccf ratio image: ${plot_root}_2Dccf.fits"
else
	echo -e "\tERROR: FIMGCREATE was not successful."
fi

#####################
## Plotting the lags
#####################

cd "$lag_exe_dir"

if (( $testing == 0 )); then
	out_file="$lag_out_dir/${prefix}_${day}_t${dt}_${numsec}sec"
	plot_root="$lag_out_dir/${prefix}_${day}_t${dt}_${numsec}sec"
elif (( $testing == 1 )); then
	out_file="$lag_out_dir/test_${prefix}_${day}_t${dt}_${numsec}sec"
	plot_root="$lag_out_dir/test_${prefix}_${day}_t${dt}_${numsec}sec"
fi

# if [ -e "${out_file}.${tab_ext}" ]; then
# 	python "$lag_exe_dir"/plot_lags.py "${out_file}.${tab_ext}" \
# 		-o "${plot_root}" -l "$lag_lf" -u "$lag_uf"
# # 	if [ -e "${plot_root}_lag-energy.png" ]; then open -a ImageJ "${plot_root}_lag-energy.png"; fi
# else
# 	echo -e "\tERROR: plot_lags.py was not run. Lag output file does not exist."
# fi

################################################################################
## All done!
################################################################################
