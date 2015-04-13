# cross_correlation
Computes the cross-correlation function of a broad reference energy band with narrow energy channels of interest.

## Contents

### ccf.py
Computes the cross-correlation function of two light curves to do phase-resolved spectroscopy. 
Able to read light curves from data with the program 'populate_lightcurve'.

### loop_ccf.sh
Bash script to compute the ccf of individual observations and stitch the plots together in a gif!

### multi_ccf.py
Computes the cross-correlation function of a band of interest with a reference band over multiple data files.

### plot_2d.py
Plots the two-dimensional cross correlation function.

### plot_ccf.py
Plots the cross-correlation function in one dimension (i.e. for one energy channel).

### plot_multi.py
Plots multiple 1-D cross-correlation functions from different energy channels on one plot.

### run_ccf.sh
Bash script to run ccf.py and plotting scripts.

### run_multi_ccf.sh
Bash script to run multi_ccf.py and plotting scripts.



THIS CODE COMES WITH NO LEGAL GUARANTEES.