# cross_correlation
Computes the cross-correlation function of a broad reference energy band with 
narrow energy channels of interest.

## Contents

### ccf.py
Computes the cross-correlation function of a broad reference energy band with 
narrow energy channels of interest (detector energy channels) to do 
phase-resolved spectroscopy. Able to make light curves from an eventlist with 
helper methods in 'tools', in the whizzy_scripts repo.

### ccf_bootstrap.py
Like ccf.py, but to be used for bootstrapping the data, where it selects a 
certain subset (with replacement) of the segments and goes on from there. To be
used in run_multi_ccf_bootstrap.sh.

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

### README.md
This document.

### run_ccf.sh
Bash script to run ccf.py and plotting scripts.

### run_multi_ccf.sh
Bash script to run multi_ccf.py and plotting scripts.

### run_multi_ccf_bootstrap.sh
Bash script to run ccf_bootstrap.py and plotting scripts.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) 

## Authors and License
* Abigail Stevens (UvA)

Pull requests are welcome!

All code is Copyright 2014-2015 The Authors, and is distributed under the MIT 
Licence. See LICENSE for details. If you are interested in the further 
development of cross_correlation, please [get in touch via the issues](https://github.com/abigailstev/cross_correlation/issues)!

