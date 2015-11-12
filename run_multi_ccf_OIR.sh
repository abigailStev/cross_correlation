#!/bin/bash

################################################################################
## 
## Bash script to run multi_ccf_OIR.py, plot_ccf.py, plot_multi.py, plot_2d.py,
## and plot_lags.py
##
## Example call:  ./run_multi_ccf_OIR.sh ./eventlists.lst ./IR_ref_band.fits \
##                      cygx1 4 16 0 150131 no
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: HEASOFT 6.14 (or higher), bash 3.*, and Python 2.7.* (with supporting 
## 		  libraries) must be installed in order to run this script. 
##
## Written by Abigail Stevens, A.L.Stevens at uva.nl, 2014-2015
##
################################################################################

## Checking the number of input arguments
if (( $# != 8 )); then
    echo -e "\tUsage: ./run_multi_ccf.sh <file list> <ref band file> <prefix>"\
            " <dt multiple> <num seconds> <testing> <date> <filtering>\n"
    exit
fi

file_list=$1
ref_band_file=$2
prefix=$3
dt=$4
numsec=$5
testing=$6
day=$7
filtering=$8 ## "no" for QPOs, or "lofreq-hifreq" in Hz for coherent pulsations

################################################################################

## If heainit isn't running, start it
if (( $(echo $DYLD_LIBRARY_PATH | grep heasoft | wc -l) < 1 )); then
	. ${HEADAS}/headas-init.sh
fi

home_dir=$(ls -d ~)

exe_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exe_dir/out_ccf/${prefix}"
lag_exe_dir="$home_dir/Dropbox/Research/lags"
lag_out_dir="$lag_exe_dir/out_lags/${prefix}"
xte_exe_dir="$home_dir/Dropbox/Research/rxte_reduce"
red_dir="$home_dir/Reduced_data/${prefix}"

ec_table_file="$xte_exe_dir/e-c_table.txt"
chan_bin_file="$red_dir/chan.txt"
energies_file="$red_dir/energies.txt"

bkgd_spec="$home_dir/Reduced_data/$prefix/evt_bkgd_rebinned.pha"

filename="${prefix}_${day}_t${dt}_${numsec}sec"

lag_lf=4.0  ## Lower frequency bound for lag spectra, in Hz
lag_uf=7.0  ## Upper frequency bound for lag spectra, in Hz
lag_le=3
lag_ue=9

tlen=30 ## Number of time bins to plot along the 2D CCF x-axis
obs_epoch=5

p_ext="eps"
#p_ext="pdf"

################################################################################
################################################################################

if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi
if [ ! -d "$lag_out_dir" ]; then mkdir -p "$lag_out_dir"; fi

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

if [ -e "$saved_file_list" ] && [ -e "$bkgd_spec" ] && [ -e "${ref_band_file}" ]
	then
	python "$exe_dir"/multi_ccf_OIR.py "$saved_file_list" \
			"${out_file}.fits" \
			"${ref_band_file}" \
			-b "$bkgd_spec" \
			-n "$numsec" \
			-m "$dt" \
			-t "$testing" \
			-f "$filtering"

elif [ -e "$saved_file_list" ] && [ -e "${ref_band_file}" ]; then
	python "$exe_dir"/multi_ccf_OIR.py "$saved_file_list" \
			"${out_file}.fits" \
			"${ref_band_file}" \
			-n "$numsec" \
			-m "$dt" \
			-t "$testing" \
			-f "$filtering"
else
	echo -e "\tERROR: multi_ccf.py was not run. List of eventlists and/or "\
			"reference band file and/or background energy spectrum doesn't "\
			"exist."
fi

################################
## Plotting the individual ccfs
################################

if [ -e "${out_file}.fits" ]; then
	python "$exe_dir"/plot_CCF.py "${out_file}.fits" \
	        -o "${plot_root}" -p "$prefix" --ext "$p_ext"

	if [ -e "${plot_root}_chan_15.${p_ext}" ]; then
	    open "${plot_root}_chan_15.${p_ext}"; fi

	multi_plot="${plot_root}_multiccfs.${p_ext}"
	python "$exe_dir"/plot_multi.py "${out_file}.fits" "$multi_plot" \
		    -p "$prefix"
#	if [ -e "$multi_plot" ]; then open "$multi_plot"; fi
fi

###############################################
## Getting the energy list from a channel list
###############################################

if [ ! -e "$energies_file" ]; then
	if [ -e "$ec_table_file" ] && [ -e "$chan_bin_file" ]; then
		python "$xte_exe_dir"/channel_to_energy.py "$ec_table_file" \
			"$chan_bin_file" "$energies_file" "$obs_epoch"
	else
		echo -e "\tERROR: channel_to_energy.py not run. ec_table_file and/or "\
                "chan_bin_file do not exist."
	fi
fi


#######################
## Plotting the 2D ccf
#######################

plot_file="${plot_root}_2Dccf.${p_ext}"
if [ -e "${out_file}.fits" ]; then
	python "$exe_dir"/plot_2d.py "${out_file}.fits" -o "${plot_file}" \
		    -p "$prefix" -l "$tlen" -e "$energies_file"

	if [ -e "${plot_file}" ]; then
		open "${plot_file}"
    fi
fi

plot_file="${plot_root}_2Dccf.fits"
detchans=$(python -c "import tools; print int(tools.get_key_val('${out_file}.fits', 0, 'DETCHANS'))")
tlen2=$(( 2*tlen ))

if [ -e "$exe_dir/temp.dat" ]; then
	fimgcreate bitpix=-32 \
		naxes="${tlen2},${detchans}" \
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
	out_file="$lag_out_dir/$filename"
	plot_root="$lag_out_dir/$filename"
elif (( $testing == 1 )); then
	out_file="$lag_out_dir/test_$filename"
	plot_root="$lag_out_dir/test_$filename"
fi

if [ -e "${out_file}_cs.fits" ]; then

	python "$lag_exe_dir"/get_lags.py "${out_file}_cs.fits" \
			"${out_file}_lag.fits" -o "${plot_root}" -p "$prefix" \
			-e "${p_ext}" --lf "${lag_lf}" --uf "${lag_uf}" --le "${lag_le}" \
			--ue "${lag_ue}"

	if [ -e "$plot_root"_lag-energy."${p_ext}" ]; then
	    open "$plot_root"_lag-energy."${p_ext}"
	fi
#	if [ -e "$plot_root"_lag-freq."${p_ext}" ]; then
#       open "$plot_root"_lag-freq."${p_ext}"
#   fi

else
	echo -e "\tERROR: plot_lags.py was not run. Lag output file does not exist."
fi

################################################################################
## All done!
################################################################################
