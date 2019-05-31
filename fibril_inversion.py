# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield


Creates a class for a given fibril data.

Defines functions to fit trends and fit sunnysoids to detrended data.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import matplotlib.pyplot as plt
import alfven_speed_inversion as asi
import warnings


class Fibril:
    def __init__(self, yb, yt, t_vals, pixel_size=50., cadence=7.68,
                 trend_range=None):
        """
        Inputs:
            yb = numpy array of position values for the bottom boundary,
            yt = numpy array of position values for the top boundary,
            t_vals = numpy array of corresponding time values,
            pixel_size = size of pixels used in the observation (km),
            cadence = time between observation time frames (s),
            trend_range = slice(start, end).
        """
        if len(t_vals) != len(yb):
            raise ValueError("Each time value must have an associated yb and "
                             "yt value")
        self.yb = yb
        self.yt = yt
        self.t_vals = t_vals
        self.t_vals_cont = np.linspace(self.t_vals[0], self.t_vals[-1], 1000)
        self.pixel_size = pixel_size
        self.cadence = cadence
        self.trend_range = trend_range

    @property
    def yb_sub(self):
        if self.trend_range is None:
            return self.yb
        elif type(self.trend_range) is not slice:
            raise ValueError("trend_range must be a slice")
        else:
            return self.yb[self.trend_range]

    @property
    def yt_sub(self):
        if self.trend_range is None:
            return self.yt
        elif type(self.trend_range) is not slice:
            raise ValueError("trend_range must be a slice")
        else:
            return self.yt[self.trend_range]

    @property
    def t_vals_sub(self):
        if self.trend_range is None:
            return self.t_vals
        elif type(self.trend_range) is not slice:
            raise ValueError("trend_range must be a slice")
        else:
            return self.t_vals[self.trend_range]

    @property
    def t_vals_cont_sub(self):
        return np.linspace(self.t_vals_sub[0], self.t_vals_sub[-1], 1000)

    def pix_to_km(self):
        """
        Convert boundary data from units of pixels to km
        """
        self.yb = self.yb * self.pixel_size
        self.yt = self.yt * self.pixel_size

    def smooth(self, indices_of_anomalies):
        """
        Linearly smooth over any anomolies.

        Inputs:
            indicies_of_anomolies = list of indices of isolated anomolies or of
            lists of non-isolated anomolies.
        """
        for i in indices_of_anomalies:
            if type(i) is int:
                self.yb[i] = 0.5*(self.yb[i+1] + self.yb[i-1])
                self.yt[i] = 0.5*(self.yt[i+1] + self.yt[i-1])
            elif type(i) is list:
                length = len(i)
                prev = i[0] - 1
                after = i[-1] + 1
                grad_b = (self.yb[after] - self.yb[prev]) / (length + 1)
                grad_t = (self.yt[after] - self.yt[prev]) / (length + 1)
                for n, j in enumerate(i):
                    self.yb[j] = self.yb[prev] + n * grad_b
                    self.yt[j] = self.yt[prev] + n * grad_t

    def trend(self, data, N=3):
        """
        Return a trend function (for each boundary) via least-squares
        regression for a polynomial of degree N.

        Inputs:
            data = input id data as numpy array,
            N = degree of trend polynomial.
        """
        x_vals = np.arange(len(data))
        polyb = np.poly1d(np.polyfit(x_vals, data, N))
        polyt = np.poly1d(np.polyfit(x_vals, data, N))
        return [polyb(x_vals), polyt(x_vals)]

    def trend_plot(self, N=3, savefig=None):
        """
        Plot the top and bottom boundaries with overlayed trends

        Inputs:
            N = degree of trend polynomial,
            savefig = None (not saved) or saved name string.
        """
        plt.figure()
        plt.errorbar(self.t_vals, self.yb, yerr=self.pixel_size, color='black')
        plt.errorbar(self.t_vals, self.yt, yerr=self.pixel_size, color='black')
        plt.plot(self.t_vals_sub, self.trend(self.yb_sub, N=N)[0], color='red')
        plt.plot(self.t_vals_sub, self.trend(self.yt_sub, N=N)[1], color='red')
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (km)")
        if savefig is not None:
            plt.savefig(savefig)

    def width(self, N=3):
        """
        Calculate average with of the fibril through the given time values.

        Inputs:
            N = degree of trend polynomial,
        """
        return np.average(self.trend(self.yt_sub, N)[1] -
                          self.trend(self.yb_sub, N)[0])

    def detrend(self, N=3, trend_type="poly", plot=False, savefig=None):
        """
        Return data minus the trend.

        Inputs:
            N = degree of trend polynomial,
            trend_type = 'poly' uses polynomial trend of degree N,
                       = 'diff' uses difference trend,
            savefig = None (not saved) or saved name string.
        """
        if trend_type == "poly":
            detrend_b = self.yb_sub - self.trend(self.yb_sub, N)[0]
            detrend_t = self.yt_sub - self.trend(self.yt_sub, N)[1]
        elif trend_type == "diff":
            detrend_b = self.yb_sub[1:] - self.yb_sub[:-1]
            detrend_t = self.yt_sub[1:] - self.yt_sub[:-1]
        else:
            raise ValueError('trend_type must be either "poly" or "diff"')

        if plot is True:
            plt.figure()
            plt.errorbar(self.t_vals, detrend_b, yerr=self.pixel_size,
                         color='black')
            plt.errorbar(self.t_vals, detrend_t, yerr=self.pixel_size,
                         color='black')
            plt.xlabel("Time (s)")
            plt.ylabel("Distance (km)")
            if savefig is not None:
                plt.savefig(savefig)

        return [detrend_b, detrend_t]

    def sin_fitting(self, p0, N=3, trend_type="poly", plot=False,
                    savefig=None):
        """
        Fit sinusoid to each boundary oscillation

        Inputs:
            p0 = initial [amplitude, frequency, shift],
            N = degree of trend polynomial,
            trend_type = 'poly' uses polynomial trend of degree N,
                       = 'diff' uses difference trend,
            savefig = None (not saved) or list of saved name strings.
        """
        if trend_type == "diff":
            t_vals_sub = self.t_vals_sub[1:]
        else:
            t_vals_sub = self.t_vals_sub
        
        def sin_func(x, a, b, c):  # a=amplitude, b=angular frequency
            return a * np.sin(b * x + c)

        # different p0 for top and bottom
        if type(p0[0]) is not list:
            p0_b = p0
            p0_t = p0
        elif type(p0[0]) is list:
            p0_b = p0[0]
            p0_t = p0[1]

        yb_params, yb_params_covar = curve_fit(sin_func, t_vals_sub,
                                               self.detrend(N, trend_type)[0],
                                               p0=p0_b)
        yt_params, yt_params_covar = curve_fit(sin_func, t_vals_sub,
                                               self.detrend(N, trend_type)[1],
                                               p0=p0_t)
#        print("\nTOP: Amp = " "%.4g" % yt_params[0] + " km, Freq = "
#              "%.4g" % yt_params[1] + " s-1 \n\nBOTTOM: Amp = "
#              "%.4g" % yb_params[0] + " km, Freq = "
#              "%.4g" % yb_params[1] + " s-1\n")

        sin_fit = [sin_func(self.t_vals_cont_sub, yb_params[0], yb_params[1],
                            yb_params[2]), yb_params,
                   sin_func(self.t_vals_cont_sub, yt_params[0], yt_params[1],
                            yt_params[2]), yt_params]
        if plot is True:
            y_lim = max(max(self.detrend(N, trend_type)[0]),
                        max(self.detrend(N, trend_type)[1])) + self.pixel_size
            
            # Bottom boundary oscillation
            plt.figure()
            plt.errorbar(t_vals_sub, self.detrend(N, trend_type)[0],
                         yerr=self.pixel_size, color='black')
            plt.plot(self.t_vals_cont_sub, sin_fit[0], color='red')
            plt.ylim([-y_lim, y_lim])
            plt.xlabel("Time (s)")
            plt.ylabel("Distance (km)")
            if savefig is not None:
                plt.savefig(savefig[0])

            # Top boundary oscillation
            plt.figure()
            plt.errorbar(t_vals_sub, self.detrend(N, trend_type)[1],
                         yerr=self.pixel_size, color='black')
            plt.plot(self.t_vals_cont_sub, sin_fit[2], color='red')
            plt.ylim([-y_lim, y_lim])
            plt.xlabel("Time (s)")
            plt.ylabel("Distance (km)")

            if savefig is not None:
                plt.savefig(savefig[1])

        return sin_fit

    def AR_inversion(self, p0, N, vA_guess, c_phase, c0, R1, R2, mode):
        """
        Use amplitude ratio technique to acheive an inversion for the Alfven
        speed inside the fibril.

        Inputs:
            p0 = initial [amplitude, frequency, shift],
            N = degree of trend polynomial,
            vA_guess = initial Alfven speed guess,
            Prescribed parameters:
                c_phase = phase speed,
                c0 = sound speed,
                R1 = rho1 / rho0,
                R2 = rho2 / rho0,
                mode = identified mode, either 'saus' or 'kink'.
        """
        sin_fit = self.sin_fitting(p0, N)
        w = (sin_fit[1][1] + sin_fit[3][1]) / 2.  # freq
        k = w / c_phase  # wavenumber
        x0 = self.width(N) / 2  # half-width of structure
        # amplitude ratio
        RA = sin_fit[3][0] / sin_fit[1][0]
        if mode == "saus":
            RA = -abs(RA)
        elif mode == "kink":
            RA = abs(RA)

        if mode == "saus" and RA >= 0:
            warnings.warn("Warning: You have identified this mode as 'saus' "
                          "but it has a ...........Message")

        vA_sol = asi.alfven_AR_inversion(w, k, vA_guess, c0, R1, R2, x0, RA,
                                         mode)
        return vA_sol

    def AR_inversion_multiple_init(self, p0, N, vA_inits, c_phase, c0, R1,
                                   R2, mode):
        """
        Use amplitude ratio technique to acheive an inversion for the Alfven
        speed inside the fibril using multiple initial vA values to reduce the
        chance of multiple roots.

        Inputs:
            p0 = initial [amplitude, frequency, shift],
            N = degree of trend polynomial,
            vA_guess = numpy array of initial Alfven speed guesses,
            Prescribed parameters:
                c_phase = phase speed,
                c0 = sound speed,
                R1 = rho1 / rho0,
                R2 = rho2 / rho0,
                mode = identified mode, either 'saus' or 'kink'.

        Output:
            [Alven speed estimate,
             proportion of initial guesses that inverted to this value]
        """

        # Initialise the vector for inverted values
        vA_inversion_vec = np.zeros_like(vA_inits)

        for i, vA_guess in enumerate(vA_inits):
            vA_inversion_exact = self.AR_inversion(p0, N, vA_guess, c_phase,
                                                   c0, R1, R2, mode)[0]
            vA_inversion_rounded = round(vA_inversion_exact, 2)
            vA_inversion_vec[i] = abs(vA_inversion_rounded)

        # Calculate the modal value of the vector of inversions
        vA_inversion = stats.mode(vA_inversion_vec)

        # Set the threshold for necessary proportion of consistent inversions
        threshold_inversion_proportion = 0.75
        sample_inversion_proportion = vA_inversion[1][0] / len(vA_inits)

        # Indicators for multiple roots
        if len(vA_inversion[0]) > 1:
            print('WARNING: There was more than one mode.')
        elif sample_inversion_proportion < threshold_inversion_proportion:
            print('WARNING: The most common inversion accounted for less ' +
                  'than ' + str(threshold_inversion_proportion) + ' of the ' +
                  'inversions.')

        return np.array([vA_inversion[0][0], sample_inversion_proportion])
