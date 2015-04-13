# cross_correlation
Computes the cross-correlation function of a broad reference energy band with 
narrow energy channels of interest.

## Contents

### ccf.py
Computes the cross-correlation function of a broad reference energy band with 
narrow energy channels of interest (detector energy channels) to do 
phase-resolved spectroscopy. Able to make light curves from an eventlist with 
helper methods in 'tools', in the whizzy_scripts repo.

### loop_ccf.sh
Bash script to compute the ccf of individual observations and stitch the plots 
together in a gif!

### multi_ccf.py
Same as ccf.py but over multiple data files.
Computes the cross-correlation function of channels of interest with a reference
band over multiple data files.

### plot_2d.py
Plots the two-dimensional cross correlation function.

### plot_ccf.py
Plots the cross-correlation function in one dimension (i.e. for one energy 
channel).

### plot_cs.py
Plots the cross spectrum. Only intended for testing purposes.

### plot_multi.py
Plots multiple 1-D cross-correlation functions from different energy channels on
one plot.

### run_ccf.sh
Bash script to run ccf.py and plotting scripts.

### run_multi_ccf.sh
Bash script to run multi_ccf.py and plotting scripts.



##### Disclaimer: This code comes with no legal guarantees.