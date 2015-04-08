#!/bin/bash

################################################################################
##
## Simple script to run ccf.py, plot_ccf.py, plot_multi.py, plot_2d.py, and 
## plot_lags.py
##
## Don't give command line arguments. Change things in this script below.
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: HEASOFT 6.11.*, bash 3.*, and Python 2.7.* (with supporting libraries) 
## 		  must be installed in order to run this script. 
##
## Abigail Stevens, A.L.Stevens@uva.nl, 2014-2015
##
################################################################################

## Checking the number of input arguments
if (( $# != 0 )); then
    echo -e "\tDo not give command line arguments. Usage: ./run_ccf.sh\n"
    exit
fi

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	. ${HEADAS}/headas-init.sh
fi

################################################################################

home_dir=$(ls -d ~) 
day=$(date +%y%m%d)
exe_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exe_dir/out_ccf"
xte_exe_dir="$home_dir/Dropbox/Research/rxte_reduce"
lag_exe_dir="$home_dir/Dropbox/Research/lags"
lag_out_dir="$lag_exe_dir/out_lags"

# prefix="P70080"
# obsID="70080-01-01-02"
prefix="GX339-BQPO"
obsID="95335-01-01-06"

red_dir="$home_dir/Reduced_data/${prefix}/$obsID"
# red_dir="$home_dur/Dropbox/Research/sample_data"
in_file="$red_dir/GTId_eventlist.fits"
bkgd_spec="$home_dir/Reduced_data/$prefix/evt_bkgd_rebinned.pha"
ec_table_file="$xte_exe_dir/e-c_table.txt"
chan_bin_file="$red_dir/chan.txt"
energies_file="$red_dir/energies.txt"

dt=64
numsec=64
testing=0  # 0 for no, 1 for yes
filtering=0 # 0 for no, 1 for yes
t_len=80
obs_epoch=5

t_ext="fits"

################################################################################
################################################################################

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi
if [ ! -d "$lag_out_dir" ]; then mkdir -p "$lag_out_dir"; fi

if (( $testing == 0 )); then
	out_file="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/${obsID}_${day}_t${dt}_${numsec}sec"
elif (( $testing == 1 )); then
	out_file="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
	plot_root="$out_dir/test_${obsID}_${day}_t${dt}_${numsec}sec"
fi

##################
## Running ccf.py
##################

# for (( i=0; i<64; i++ )); do
# 	tmp_file="$out_dir/${obsID}_ccf_segs_${i}.dat"
# 	if [ -e "$tmp_file" ]; then rm "$tmp_file"; fi; touch "$tmp_file"
# done

if [ -e "$in_file" ] && [ -e "$bkgd_spec" ]; then
	time python "$exe_dir"/ccf.py "${in_file}" "${out_file}.${t_ext}" \
		-b "$bkgd_spec" -n "$numsec" -m "$dt" -t "$testing" -f "$filtering"
elif [ -e "$in_file" ]; then
	time python "$exe_dir"/ccf.py "${in_file}" "${out_file}.${t_ext}" \
		-n "$numsec" -m "$dt" -t "$testing" -f "$filtering"
else
	echo -e "\tERROR: ccf.py was not run. Eventlist and/or background energy \
spectrum doesn't exist."
fi


#################
## Plotting ccfs
#################

if [ -e "${out_file}.${t_ext}" ]; then

	python "$exe_dir"/plot_ccf.py "${out_file}.${t_ext}" -o "${plot_root}" \
		-p "${prefix}/${obsID}"
	
	if [ -e "${plot_root}_chan_06.${p_ext}" ]; then open "${plot_root}_chan_06.${p_ext}"; fi

	ccfs_plot="${plot_root}_multiccfs.${p_ext}"
	python "$exe_dir"/plot_multi.py "${out_file}.${t_ext}" "$ccfs_plot" \
		-p "${prefix}/${obsID}"

	if [ -e "$ccfs_plot" ]; then open "$ccfs_plot"; fi

fi

###############################################
## Getting the energy list from a channel list
###############################################

if [ ! -e "$energies_file" ]; then
	if [ -e "$ec_table_file" ] && [ -e "$chan_bin_file" ]; then
		python "$xte_exe_dir"/channel_to_energy.py "$ec_table_file" \
			"$chan_bin_file" "$energies_file" "$obs_epoch"
	else
		echo -e "\tERROR: channel_to_energy.py not run. ec_table_file and/or \
chan_bin_file do not exist."
	fi
fi	

####################################
## Plotting the 2D ccf with colours
####################################

plot_file="${plot_root}_2Dccf.${p_ext}"
if [ -e "${out_file}.${t_ext}" ]; then
	python "$exe_dir"/plot_2d.py "${out_file}.${t_ext}" -o "${plot_file}" \
		-p "${prefix}/${obsID}" -l "$tlen" -e "$energies_file"
# 	if [ -e "${plot_file}" ]; then open "${plot_file}"; fi
fi

plot_file="${plot_root}_2Dccf.fits"
detchans=$(python -c "import tools; print int(tools.get_key_val('${out_file}.fits', 0, 'DETCHANS'))")

if [ -e "$out_dir/temp.dat" ]; then
	fimgcreate bitpix=-32 \
		naxes="${tlen},${detchans}" \
		datafile="$exe_dir/temp.dat" \
		outfile="${plot_root}_2Dccf.fits" \
		nskip=1 \
		history=true \
		clobber=yes
else
	echo -e "\tERROR: FIMGCREATE did not run. 2Dccf temp file does not exist."
fi

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

if [ -e "${out_file}.${t_ext}" ]; then
	python "$lag_exe_dir"/plot_lags.py "${out_file}.${t_ext}" -o "${plot_root}" \
		-p "${prefix}/${obsID}"
	if [ -e "$plot_root"_lag-energy.png ]; then open "$plot_root"_lag-energy.png; fi
	if [ -e "$plot_root"_lag-freq_15.png ]; then open "$plot_root"_lag-freq_15.png; fi
else
	echo -e "\tERROR: plot_lags.py was not run. Lag output file does not exist."
fi

################################################################################
## All done!
################################################################################
