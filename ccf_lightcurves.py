"""
Class for CCF Lightcurves object and normalized power spectra.
Typically imported as ccf_lc. At top of program, along with other import
statements, write: "import ccf_lightcurves.py as ccf_lc", then you can call
ccf_lc.Lightcurve and ccf_lc.NormPSD!
"""
import numpy as np

__author__ = 'Abigail Stevens <A.L.Stevens at uva.nl>'
__year__ = "2015-2016"


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


class LagFreq(object):
    def __init__(self, n_bins=8192):
        self.cross = np.zeros(n_bins, dtype=np.complex128)
        self.pow1 = np.zeros(n_bins, dtype=np.float64)
        self.pow2 = np.zeros(n_bins, dtype=np.float64)


class Cross(object):
    def __init__(self, n_bins=8192, detchans=64):
        self.freq = np.zeros(n_bins, dtype=np.float64)
        self.ci = np.zeros((n_bins, detchans), dtype=np.float64)
        self.ci_noise = np.zeros((n_bins, detchans), dtype=np.float64)
        self.ref = np.zeros(n_bins, dtype=np.float64)
        self.ref_noise = np.zeros(n_bins, dtype=np.float64)


class Filter(object):
    def __init__(self, pars, n_bins=8192, dt=0.0078125):
        """

        Parameters
        ----------
        pars : 1-D np.array of floats
            Parameters from fitting power spectra, for one power spectrum.
            pars[0] = power law photon index
            pars[1] = power law norm
            pars[2] = bbn1 sigma
            pars[3] = bbn1 centroid frequency
            pars[4] = bbn1 norm
            pars[5] = bbn2 sigma
            pars[6] = bbn2 centroid frequency
            pars[7] = bbn2 norm
            pars[8] = subharmonic sigma
            pars[9] = subharmonic norm
            pars[10] = qpo fundamental sigma
            pars[11] = qpo fundamental centroid frequency
            pars[12] = qpo fundamental norm
            pars[13] = qpo harmonic norm

            (subharmonic centroid and qpo harmonic centroid are computed from qpo
            fundamental centroid. qpo harmonic sigma is computed from qpo
            fundamental sigma.)

        n_bins : int
            Number of bins in one Fourier transform segment

        dt : float
            Sampling timestep for the light curve

        Attributes
        ----------
        pos_freq :
        powerlaw :
        bbn1 :
        bbn2 :
        subharm :
        fund :
        harm :
        continuum :
        whole_continuum :
        whole_fund :
        whole_harm :
        fund_filt :
        harm_filt :

        """
        self.pos_freq = np.abs(fftpack.fftfreq(n_bins, d=dt)[0:n_bins/2+1])

        self.powerlaw = __xspec_powerlaw(pars[0], pars[1])
        self.bbn1 = __xspec_lorf(pars[2], pars[3], pars[4])
        self.bbn2 = __xspec_lorf(pars[5], pars[6], pars[7])
        self.subharm = __xspec_lorf(pars[8], pars[11]/2.0, pars[9])
        self.fund = __xspec_lorf(pars[10], pars[11], pars[12])
        self.harm = __xspec_lorf(2.0*pars[10], 2.0*pars[11], pars[13])
        self.continuum = powerlaw + bbn1 + bbn2 + subharm + fund + harm

        nf_continuum = self.continuum[1:-1]
        self.whole_continuum = np.append(self.continuum, nf_continuum[::-1])
        nf_fund = self.fund[1:-1]
        self.whole_fund = np.append(self.fund, nf_fund[::-1])
        nf_harm = self.harm[1:-1]
        self.whole_harm = np.append(self.harm, nf_harm[::-1])

        self.fund_filt = self.whole_fund ** 2 / self.whole_continuum ** 2
        self.harm_filt = self.whole_harm ** 2 / self.whole_continuum ** 2

    def __xspec_lorf(self, sigma, lineE, norm):
        """
        The lorentz function as defined by xspec, times frequency.
        Note that sigma here is the full width half max, and lineE is the centroid
        frequency.

        :param freq:
        :param sigma:
        :param lineE:
        :param norm:

        Returns
        -------
        The lorentzian*f function evaluated at every input frequency.
        """
        return norm * self.pos_freq * sigma / (2 * np.pi) / \
               ((self.pos_freq - lineE) ** 2 + ((sigma / 2) ** 2))

    def __xspec_powerlaw(self, phoindex, norm):
        """
        The powerlaw function as defined by xspec.
        Note that phoindex is automatically made negative in here, so a negative
        phoindex input returns a positive slope!
        :param freq:
        :param phoindex:
        :param norm:

        Returns
        -------
        The powerlaw function evaluated at every input frequency.
        norm*freq**(-phoindex)
        """
        return norm * self.pos_freq ** (-phoindex)