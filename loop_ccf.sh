#!/bin/bash

################################################################################
##
## Bash script to run ccf.py, plot_ccf.py, plot_multi.py, and plot_2Dccf.py.
##
## Runs ccf.py for many obsIDs.
## 
## Don't give command line arguments. Change things in this script below.
##
## Change the directory names and specifiers before the double '#' row to best
## suit your setup.
##
## Notes: bash 3.* and Python 2.7.* (with supporting libraries) must be 
##		  installed in order to run this script. For the gif-making to work, 
##		  gifsicle must be installed (open source, available on e.g. MacPorts 
##		  and HomeBrew)
## 
## Written by Abigail Stevens, A.L.Stevens at uva.nl, 2015
## 
################################################################################

##########################################
## Checking the number of input arguments
##########################################

if (( $# != 0 )); then
    echo -e "\tDo not give command line arguments. Usage: ./loop_ccf.sh\n"
    exit
fi

################################################################################

home_dir=$(ls -d ~)
day=$(date +%y%m%d)  # make the date a string and assign it to 'day'

exe_dir="$home_dir/Dropbox/Research/cross_correlation"
out_dir="$exe_dir/out_ccf"

# prefix="j1808-2002"
prefix="j1808-1HzQPO"
# prefix="GX339-BQPO"
# obsID="95335-01-01-06"
# prefix="4u1636superburst"

obsID_list="$home_dir/Dropbox/Lists/${prefix}_obsIDs_goodSN.lst"
# bkgd_spec="$home_dir/Reduced_data/$prefix/evt_bkgd_rebinned.pha"
ec_table_file="$xte_exe_dir/e-c_table.txt"
chan_bin_file="$home_dir/Reduced_data/${prefix}/chan.txt"
energies_file="$home_dir/Reduced_data/${prefix}/energies.txt"

dt=128
numsec=128
testing=0  # 0 for no, 1 for yes
filtering=0 # 0 for no, 1 for yes
tlen=100
obs_epoch=5

t_ext="fits"
p_ext="png"

plots_1d="$out_dir/${prefix}_gif_1d_goodSN.txt"
plots_2d="$out_dir/${prefix}_gif_2d_goodSN.txt"

gif_name_1d="$out_dir/${day}_t${dt}_${numsec}sec_1d_goodSN.gif"
gif_name_2d="$out_dir/${day}_t${dt}_${numsec}sec_2d_goodSN.gif"


################################################################################
################################################################################

if [ -e "$plots_1d" ]; then rm "$plots_1d"; fi; touch "$plots_1d"
if [ -e "$plots_2d" ]; then rm "$plots_2d"; fi; touch "$plots_2d"

if [ ! -e "$energies_file" ]; then
	if [ -e "$ec_table_file" ] && [ -e "$chan_bin_file" ]; then
		python "$xte_exe_dir"/channel_to_energy.py "$ec_table_file" \
			"$chan_bin_file" "$energies_file" "$obs_epoch"
	else
		echo -e "\tERROR: channel_to_energy.py not run. ec_table_file and/or \
chan_bin_file do not exist."
	fi
fi	

for obsID in $( cat $obsID_list ); do

	red_dir="$home_dir/Reduced_data/${prefix}/$obsID"
# 	red_dir="$home_dur/Dropbox/Research/sample_data"
	in_file="$red_dir/GTId_eventlist.fits"
		
	if [ ! -d "$out_dir" ]; then mkdir -p "$out_dir"; fi

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
	
	for (( i=0; i<64; i++ )); do
		tmp_file="$out_dir/ccf_segs_${i}.dat"
		if [ -e "$tmp_file" ]; then rm "$tmp_file"; fi; touch "$tmp_file"
	done
	
	if [ -e "$in_file" ] && [ -e "$bkgd_spec" ]; then
		time python "$exe_dir"/ccf.py "${in_file}" "${out_file}.${t_ext}" \
			-b "$bkgd_spec" -n "$numsec" -m "$dt" -t "$testing" -f "$filtering"
	elif [ -e "$in_file" ]; then
		time python "$exe_dir"/ccf.py "${in_file}" "${out_file}.${t_ext}" \
			-n "$numsec" -m "$dt" -t "$testing" -f "$filtering"
	else
		echo -e "\tERROR: ccf.py was not run. Eventlist and/or background \
energy spectrum doesn't exist."
	fi

	#############
	## Plotting 
	############

	if [ -e "${out_file}.${t_ext}" ]; then
		
		multi_plot="${plot_root}_multiccfs.${p_ext}"
		plot_file_2d="${plot_root}_2Dccf.${p_ext}"
		plot_fits="${plot_root}_2Dccf.fits"


		########################################
		## Plotting 1D single and multiple CCFs
		########################################
		
		python "$exe_dir"/plot_ccf.py "${out_file}.${t_ext}" \
			-o "${plot_root}" -p "${prefix}/${obsID}"
# 		if [ -e "${plot_root}_chan_06.${p_ext}" ]; then open "${plot_root}_chan_06.${p_ext}"; fi
		
		echo "${plot_root}_chan_06.${p_ext}" >> $plots_1d  
		## Could also use stars here instead of the chan num
		
		python "$exe_dir"/plot_multi.py "${out_file}.${t_ext}" "$multi_plot" \
			-p "${prefix}"
# 		if [ -e "$multi_plot" ]; then open "$multi_plot"; fi
		
		###################
		## Plotting 2D CCF
		###################
		
		if [ -e "${out_file}.${t_ext}" ]; then
		
			python "$exe_dir"/plot_2d.py "${out_file}.${t_ext}" \
				-o "${plot_file_2d}" -p "${prefix}" -l "$tlen" -e "$energies_file"
# 			if [ -e "${plot_file_2d}" ]; then open "${plot_file_2d}"; fi
			
			echo "$plot_file_2d" >> $plots_2d
			
		fi
		
		detchans=$(python -c "import tools; print int(tools.get_key_val('${out_file}.fits', 0, 'DETCHANS'))")
		
# 		if [ -e "$out_dir/temp.dat" ]; then
# 			fimgcreate bitpix=-32 \
# 				naxes="${tlen},${detchans}" \
# 				datafile="$out_dir/temp.dat" \
# 				outfile="$plot_fits" \
# 				nskip=1 \
# 				history=true \
# 				clobber=yes
# 		else
# 			echo -e "\tERROR: FIMGCREATE did not run. 2Dccf temp file does not exist."
# 		fi
#  
# 		if [ -e "$plot_fits" ]; then
# 			echo "FITS 2D ccf ratio image: $plot_fits"
# 		else
# 			echo -e "\tERROR: FIMGCREATE was not successful."
# 		fi
		
		
	else
		echo -e "\tERROR: Plots were not made. CCF output file does not exist."
	fi

done

###############################
## Making the plots into a gif
###############################

convert @"$plots_1d" "$gif_name_1d"
if [ -e "$gif_name_1d" ]; then
	echo "GIF made! $gif_name_1d"
	open "$gif_name_1d"
fi

convert @"$plots_2d" "$gif_name_2d"
if [ -e "$gif_name_2d" ]; then
	echo "GIF made! $gif_name_2d"
	open "$gif_name_2d"
fi

################################################################################
## All done!
################################################################################
