#!/bin/bash

################################################################################
## 
## Use with ccf_bootstrap.sh
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: HEASOFT 6.14 (or higher), bash 3.*, and Python 2.7.* (with supporting 
## 		  libraries) must be installed in order to run this script. 
##
## Written by Abigail Stevens, A.L.Stevens at uva.nl, 2015
##
################################################################################

## Checking the number of input arguments
if (( $# != 8 )); then
    echo -e "\tUsage: ./run_multi_ccf.sh <file list> <prefix> <dt multiple> "\
            "<num seconds> <testing> <date> <filtering> <boot_num>\n"
    exit
fi

file_list=$1
prefix=$2
dt=$3
numsec=$4
testing=$5
day=$6
filtering=$7 ## "no" for QPOs, or "lofreq-hifreq" in Hz for coherent pulsations
boot_num=$8  ## Number of bootstrap iterations to run

################################################################################

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	. ${HEADAS}/headas-init.sh
fi

home_dir=$(ls -d ~)

exe_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exe_dir/out_ccf/${prefix}/bootstrapped"
lag_exe_dir="$home_dir/Dropbox/Research/lags"
lag_out_dir="$lag_exe_dir/out_lags"
xte_exe_dir="$home_dir/Dropbox/Research/rxte_reduce"
red_dir="$home_dir/Reduced_data/${prefix}"

ec_table_file="$xte_exe_dir/e-c_table.txt"
chan_bin_file="$red_dir/chan.txt"
energies_file="$red_dir/energies.txt"

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi
if [ ! -d "$lag_out_dir" ]; then mkdir -p "$lag_out_dir"; fi

bkgd_spec="$home_dir/Reduced_data/$prefix/evt_bkgd_rebinned.pha"
filename="${prefix}_${day}_t${dt}_${numsec}sec_adj_b-"

lag_lf=4.0  ## Lower frequency bound for lag spectra, in Hz
lag_uf=7.0  ## Upper frequency bound for lag spectra, in Hz

tlen=30  ## Number of time bins to plot along the 2D CCF x-axis
obs_epoch=5

t_ext="fits"
p_ext="eps"

################################################################################
################################################################################

if (( $testing == 0 )); then
	out_file="$out_dir/$filename"
	plot_root="$out_dir/$filename"
	saved_file_list="$out_dir/${filename}_filelist"
elif (( $testing == 1 )); then
	out_file="$out_dir/test_$filename"
	plot_root="$out_dir/test_$filename"
	saved_file_list="$out_dir/test_${filename}_filelist"
fi

cp "$file_list" "$saved_file_list"

########################
## Running multi_ccf.py
########################

if [ -e "$saved_file_list" ] && [ -e "$bkgd_spec" ]; then
	python "$exe_dir"/multi_ccf_bootstrap.py "$saved_file_list" \
	        "${out_file}.${t_ext}" -b "$bkgd_spec" -n "$numsec" -m "$dt" \
	        -t "$testing" -f "$filtering" --bootstrap "$boot_num" -a
elif [ -e "$saved_file_list" ]; then
	python "$exe_dir"/multi_ccf_bootstrap.py "$saved_file_list" \
	        "${out_file}.${t_ext}" -n "$numsec" -m "$dt" -t "$testing" \
	        -f "$filtering" --bootstrap "$boot_num" -a
else
	echo -e "\tERROR: multi_ccf.py was not run. List of eventlists and/or "\
            "background energy spectrum doesn't exist."
fi
#
#################################
### Plotting the individual ccfs
#################################
#
#if [ -e "${out_file}1.${t_ext}" ]; then
#	python "$exe_dir"/plot_CCF.py "${out_file}1.${t_ext}" -o "${plot_root}1" \
#		-p "$prefix" --ext "$p_ext"
##	if [ -e "${plot_root}1_chan_15.${p_ext}" ]; then open "${plot_root}1_chan_15.${p_ext}"; fi
#
#	multi_plot="${plot_root}1_multiccfs.${p_ext}"
#	python "$exe_dir"/plot_multi.py "${out_file}1.${t_ext}" "$multi_plot" \
#		-p "$prefix"
##	if [ -e "$multi_plot" ]; then open "$multi_plot"; fi
#fi
#
################################################
### Getting the energy list from a channel list
################################################
#
#if [ ! -e "$energies_file" ]; then
#	if [ -e "$ec_table_file" ] && [ -e "$chan_bin_file" ]; then
#		python "$xte_exe_dir"/channel_to_energy.py "$ec_table_file" \
#			"$chan_bin_file" "$energies_file" "$obs_epoch"
#	else
#		echo -e "\tERROR: channel_to_energy.py not run. ec_table_file and/or "\
#                "chan_bin_file do not exist."
#	fi
#fi
#
########################
### Plotting the 2D ccf
########################
#
#plot_file="${plot_root}1_2Dccf.${p_ext}"
#if [ -e "${out_file}1.${t_ext}" ]; then
#	python "$exe_dir"/plot_2d.py "${out_file}1.${t_ext}" -o "${plot_file}" \
#		-p "$prefix" -l "$tlen" -e "$energies_file"
##	if [ -e "${plot_file}" ]; then
##		open "${plot_file}"
##        cp "$plot_file" "$home_dir/Dropbox/Research/CCF_paper1/"
##    fi
#fi
#
#plot_file="${plot_root}1_2Dccf.fits"
#detchans=$(python -c "import tools; print int(tools.get_key_val('${out_file}1.fits', 0, 'DETCHANS'))")
#tlen2=$(( 2*tlen ))
#curr_dir=$(pwd)
#cd "$out_dir"
#if [ -e "$out_dir/temp.dat" ]; then
#	fimgcreate bitpix=-32 \
#		naxes="${tlen2},${detchans}" \
#		datafile="temp.dat" \
#		outfile="${plot_root}1_2Dccf.fits" \
#		nskip=1 \
#		history=true \
#		clobber=yes
#else
#	echo -e "\tERROR: FIMGCREATE did not run. 2Dccf temp file does not exist."
#fi
#
#if [ -e "${plot_root}1_2Dccf.fits" ]; then
#	echo "FITS 2D ccf ratio image: ${plot_root}1_2Dccf.fits"
#else
#	echo -e "\tERROR: FIMGCREATE was not successful."
#fi
#cd "$curr_dir"

################################################################################
## All done!
################################################################################
