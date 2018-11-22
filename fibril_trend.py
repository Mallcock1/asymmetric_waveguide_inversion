# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield


Creates a class for a given fibril data.

Defines functions to fit trends and fit sunisoids to detrended data.
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import alfven_speed_inversion as asi
import warnings


class Fibril:
    def __init__(self, yb, yt, t_vals, trend_range=None):
        """
        trend_range = [start, end]
        """
        self.yb = yb
        self.yt = yt
        self.t_vals = t_vals
        if len(t_vals) != len(yb):
            raise ValueError("Each time value must have an associated yb and "
                             "yt value")
        self.t_vals_cont = np.linspace(self.t_vals[0], self.t_vals[-1], 1000)
        self.trend_range = trend_range

    @property
    def yb_sub(self):
        if self.trend_range is None:
            return self.yb
        elif type(self.trend_range) != list or len(self.trend_range) != 2:
            raise ValueError("trend_range must be a list of length 2")
        else:
            return self.yb[self.trend_range[0]:self.trend_range[1]]

    @property
    def yt_sub(self):
        if self.trend_range is None:
            return self.yt
        elif type(self.trend_range) != list or len(self.trend_range) != 2:
            raise ValueError("trend_range must be a list of length 2")
        else:
            return self.yt[self.trend_range[0]:self.trend_range[1]]

    @property
    def t_vals_sub(self):
        if self.trend_range is None:
            return self.t_vals
        elif type(self.trend_range) != list or len(self.trend_range) != 2:
            raise ValueError("trend_range must be a list of length 2")
        else:
            return self.t_vals[self.trend_range[0]:self.trend_range[1]]

    @property
    def t_vals_cont_sub(self):
        return np.linspace(self.t_vals_sub[0], self.t_vals_sub[-1], 1000)

    # convert units from pixels to km (1px = 50km with ROSA)
    def pix_to_km(self):
        """
        Convert from units of pixels to km, assuming 50km pixel size
        """
        self.yb = self.yb * 50
        self.yt = self.yt * 50

    # smooth over anomalous data points
    def smooth(self, indices_of_anomalies):
        """
        Linearly smooth over any known isolated anomolies.
        """
        for i in indices_of_anomalies:
            self.yb[i] = 0.5*(self.yb[i+1] + self.yb[i-1])
            self.yt[i] = 0.5*(self.yt[i+1] + self.yt[i-1])

    # fit data with an Nth degree polynomial with least squares regression.
    def trend(self, data, N=3):
        """
        Return a trend function (for each boundary) via least-squares
        regression fora polynomial of degree N.
        """
        x_vals = np.arange(len(data))
        polyb = np.poly1d(np.polyfit(x_vals, data, N))
        polyt = np.poly1d(np.polyfit(x_vals, data, N))
        return [polyb(x_vals), polyt(x_vals)]

    def trend_plot(self, N=3, savefig=None):
        """
        Plot the top and bottom boundaries with overlayed trends
        """
        plt.figure()
        plt.errorbar(self.t_vals, self.yb, yerr=50, color='black')
        plt.errorbar(self.t_vals, self.yt, yerr=50, color='black')
        plt.plot(self.t_vals_sub, self.trend(self.yb_sub, N=N)[0], color='red')
        plt.plot(self.t_vals_sub, self.trend(self.yt_sub, N=N)[1], color='red')
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (km)")
        if savefig is not None:
            plt.savefig(savefig)

    # Calculate average with of the fibril through the given time values.
    def width(self, N=3):
        return np.average(self.trend(self.yt_sub, N)[1] -
                          self.trend(self.yb_sub, N)[0])

    # Subtract the trend from the data
    def detrend(self, N=3, trend_type="poly"):
        """
        Return data minus the trend.

        'poly' uses polynomial trend of degree N,
        'diff' uses difference trend.
        """
        if trend_type == "poly":
            detrend_b = self.yb_sub - self.trend(self.yb_sub, N)[0]
            detrend_t = self.yt_sub - self.trend(self.yt_sub, N)[1]
        elif trend_type == "diff":
            detrend_b = self.yb_sub[1:] - self.yb_sub[:-1]
            detrend_t = self.yt_sub[1:] - self.yt_sub[:-1]
        else:
            raise ValueError('trend_type must be either "poly" or "diff"')
        return [detrend_b, detrend_t]

    def detrend_plot(self, N=3, trend_type="poly", savefig=None):
        # Plot the detrended data
        plt.figure()
        plt.errorbar(self.t_vals, self.detrend(N, trend_type)[0], yerr=50,
                     color='black')
        plt.errorbar(self.t_vals, self.detrend(N, trend_type)[1], yerr=50,
                     color='black')
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (km)")
        if savefig is not None:
            plt.savefig(savefig)

    # Fit a sinusoidal curve
    def sin_fitting(self, p0, N=3, trend_type="poly"):
        """
        Fit sinusoid to each boundary oscillation

        p0 = initial parameters for fitting = [amplitude, frequency, offset]
        trend_type = "poly" or "diff"
        """
        def sin_func(x, a, b, c):  # a=amplitude, b=angular frequency
            return a * np.sin(b * x + c)
        yb_params, yb_params_covariance = curve_fit(sin_func, self.t_vals_sub,
                                                    self.detrend(N, trend_type)[0],
                                                    p0=p0)
        yt_params, yt_params_covariance = curve_fit(sin_func, self.t_vals_sub,
                                                    self.detrend(N, trend_type)[1],
                                                    p0=p0)
        print("\nTOP: Amp = " "%.4g" % yt_params[0] + " km, Freq = "
              "%.4g" % yt_params[1] + " s-1 \n\nBOTTOM: Amp = "
              "%.4g" % yb_params[0] + " km, Freq = "
              "%.4g" % yb_params[1] + " s-1\n")
        return [sin_func(self.t_vals_cont_sub, yb_params[0], yb_params[1],
                         yb_params[2]), yb_params,
                sin_func(self.t_vals_cont_sub, yt_params[0], yt_params[1],
                         yt_params[2]), yt_params]

    def sin_fitting_plot(self, p0, N=3, trend_type="poly",
                         savefig=None):
        sin_fit = self.sin_fitting(p0, N)

        # Bottom boundary oscillation
        plt.figure()
        plt.errorbar(self.t_vals_sub, self.detrend(N, trend_type)[0], yerr=50,
                     color='black')
        plt.plot(self.t_vals_cont_sub, sin_fit[0],
                 color='red')
        plt.ylim([-200, 200])
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (km)")
        if savefig is not None:
            plt.savefig(savefig[0])

        # Top boundary oscillation
        plt.figure()
        plt.errorbar(self.t_vals_sub, self.detrend(N, trend_type)[1], yerr=50,
                     color='black')
        plt.plot(self.t_vals_cont_sub, sin_fit[2],
                 color='red')
        plt.ylim([-200, 200])
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (km)")
        if savefig is not None:
            plt.savefig(savefig[1])

    def AR_inversion(self, p0, N, vA_guess, c_phase, c0, R1, R2, mode):
        """
        Use amplitude ratio technique to acheive an inversion for the Alfven
        speed inside the fibril
        """
        sin_fit = self.sin_fitting(p0, N)
        w = (sin_fit[1][1] + sin_fit[3][1]) / 2.  # freq
        k = w / c_phase  # wavenumber
        x0 = self.width(N) / 2  # half-width of structure
        RA = sin_fit[3][0] / sin_fit[1][0]  # amplitude ratio

        if mode == "saus" and RA >= 0:
            warnings.warn("Warning: You have identified this mode as 'saus' "
                          "but it has a ...........Message")

        vA_sol = asi.alfven_AR_inversion(w, k, vA_guess, c0, R1, R2, x0, RA,
                                         mode)
        print("Amplitude ratio ~ " + "%.4g" % RA)
        print("AR inversion: vA ~ " + "%.4g" % vA_sol + " km s-1")
