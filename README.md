# cross_correlation
Computes the cross-correlation function of a broad reference energy band with 
narrow energy channels of interest. Please see [Stevens & Uttley 2016](https://ui.adsabs.harvard.edu/#abs/2016MNRAS.460.2796S/abstract)
for reference.

The functionality of this software will be folded into [Stingray](http://stingraysoftware.github.io/),
so please get involved there.

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

### run_ccf.sh
Bash script to run ccf.py and plotting scripts.

### run_bootstrap_ccf.sh
Bash script to run bootstrap_ccf.py.

### ccf_OIR.py, run_ccf_OIR.py
These are offshoots of their non-_OIR counterparts, created for collaboration
with Federico Vincentelli and Piergiorgio Casella. TODO: incorporate these in
the regular ccf.py and run_ccf.py.

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

## Authors
* Abigail Stevens (UvA API)


## Collaborators
* Phil Uttley (UvA API) -- Math and stats theory
* Federico Vincentelli (INAF Roma, INAF Brera) -- Adding option for IR and/or 
optical reference band

## Copyright
All content © 2014-2017 the Authors, and is distributed under the MIT
Licence. See LICENSE.md for details.

If you use this code, please cite [Stevens & Uttley 2016](https://ui.adsabs.harvard.edu/#abs/2016MNRAS.460.2796S/abstract).

The functionality of this software will be folded into [Stingray](http://stingraysoftware.github.io/),
so please get involved there.


