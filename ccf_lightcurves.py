"""
Class for CCF Lightcurves object and normalized power spectra.
Typically imported as ccf_lc. At top of program, along with other import
statements, write: "import ccf_lightcurves.py as ccf_lc", then you can call
ccf_lc.Lightcurve and ccf_lc.NormPSD!
"""
import numpy as np

__author__ = 'Abigail Stevens <A.L.Stevens at uva.nl>'
__year__ = "2015"


class Lightcurve(object):
    def __init__(self, n_bins=8192, detchans=64, type='ref'):
        self.type = type

        if type.lower() == "ci":
            self.power_array = np.zeros((n_bins, detchans, 1), dtype=np.float64)
            self.power = np.zeros((n_bins, detchans), dtype=np.float64)
            self.pos_power = np.zeros((n_bins / 2 + 1, detchans), \
                    dtype=np.float64)
            self.mean_rate_array = np.zeros((detchans, 1), dtype=np.float64)
            self.mean_rate = np.zeros(detchans)
        else:
            self.power_array = np.zeros((n_bins, 1), dtype=np.float64)
            self.power = np.zeros(n_bins, dtype=np.float64)
            self.pos_power = np.zeros(n_bins/2+1, dtype=np.float64)
            self.mean_rate_array = 0.0
            self.mean_rate = 0.0
            self.var = 0.0   # variance of the absolute-rms-normalized power
                             # spectrum of the ref band
            self.rms = 0.0   # rms of the absolute-rms-normalized power spectrum
                             #  of the ref band
            self.var_array = 0.0


class NormPSD(object):
    def __init__(self, n_bins=8192, detchans=64, type='ref'):
        self.type = type
        self.noise = 0.0

        if type.lower() == "ci":
            self.power = np.zeros((n_bins, detchans), dtype=np.float64)
            self.variance = np.zeros(detchans, dtype=np.float64)
            self.rms = np.zeros(detchans, dtype=np.float64)
        else:
            self.power = np.zeros(n_bins, dtype=np.float64)
            self.variance = 0.0
            self.rms = 0.0