# cross_correlation
Computes the cross-correlation function of a broad reference energy band with 
narrow energy channels of interest. Please see Stevens et al. (in prep) for 
reference. If you use this software, please cite that paper!!

## Contents

### ccf.py
Computes the cross-correlation function of a broad reference energy band with 
narrow energy channels of interest (detector energy channels) to do 
phase-resolved spectroscopy. Requires a FITS-format event list or light curve,
or list of FITS-format event lists or light curves, as input (as created in 
rxte_reduce/good_events.sh). Able to accept a separate file for the reference 
band, for a multi-wavelength light curve.
Requires whizzy_scripts/tools.py.

### bootstrap_ccf.py
Like ccf.py, but to be used for bootstrapping the data, where it selects a 
certain subset (with replacement) of the segments and goes on from there. To be
used in run_bootstrap_ccf.sh.

### loop_ccf.sh
Bash script to compute the ccf of individual observations and stitch the plots 
together in a gif! Requires ImageMagick to be installed, so that a gif can be 
created at the command line.

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

### run_bootstrap_ccf.sh
Bash script to run bootstrap_ccf.py.

### ccf_OIR.py, run_ccf_OIR.py
These are offshoots of their non-_OIR counterparts, created for collaboration
with Federico Vincentelli and Piergiorgio Casella. TODO: incorporate these in
the regular ccf.py and run_ccf.py.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)]
(http://www.astropy.org/) 

## Authors
* Abigail Stevens (UvA API)

Please cite Stevens et al. (in prep) if you use this paper for research!!

## Collaborators
* Phil Uttley (UvA API) -- Math and stats theory
* Federico Vincentelli (INAF Roma, INAF Brera) -- Adding option for IR and/or 
optical reference band

Pull requests are welcome!


## License
All code is Copyright 2014-2015 The Authors, and is distributed under the MIT 
Licence. See LICENSE for details. If you are interested in the further 
development of cross_correlation, please [get in touch via the issues]
(https://github.com/abigailstev/cross_correlation/issues)!



