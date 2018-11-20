# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield


Creates a class for a given fibril data.

Functions to fit trends and fit sunisoids to detrended data.
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class fibril:
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
        if trend_range is None:
            self.yb_subsection = self.yb
            self.yt_subsection = self.yt
        elif (type(trend_range) != list) or (len(trend_range) != 2):
            raise ValueError("trend_range must be a list of length 2")
        else:
            self.yb_sub = self.yb[trend_range[0]:trend_range[1]]
            self.yt_sub = self.yt[trend_range[0]:trend_range[1]]
            self.t_vals_sub = self.t_vals[trend_range[0]:trend_range[1]]

    # convert units from pixels to km (1px = 50km with ROSA)
    def pix_to_km(self):
        self.yb = self.yb * 50
        self.yt = self.yt * 50

    # smooth over anomalous data points
    def smooth(self, indices_of_anomalies):
        for i in indices_of_anomalies:
            self.yb[i] = 0.5*(self.yb[i+1] + self.yb[i-1])
            self.yt[i] = 0.5*(self.yt[i+1] + self.yt[i-1])

    # fit data with an Nth degree polynomial with least squares regression.
    def trend(data, N=3):
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
        return np.average(self.trend(N)[1] - self.trend(N)[0])

    # Subtract the trend from the data
    def detrend(self, N=3, trend_type="poly"):
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
    def sin_fitting(self, N=3, p0=[100, 0.1, 0], trend_type="poly"):
        """
        p0 = initial parameters for fitting = [amplitude, frequency, offset]

        trend_type = "poly" or "diff"
        """
        def sin_func(x, a, b, c):  # a=amplitude, b=angular frequency
            return a * np.sin(b * x + c)
        yb_params, yb_params_covariance = curve_fit(sin_func, self.t_vals,
                                                    self.detrend(N, trend_type)[0],
                                                    p0=p0)
        yt_params, yt_params_covariance = curve_fit(sin_func, self.t_vals,
                                                    self.detrend(N, trend_type)[1],
                                                    p0=p0)
        return [sin_func(self.t_vals_cont, yb_params[0], yb_params[1], yb_params[2]),
                yb_params,
                sin_func(self.t_vals_cont, yt_params[0], yt_params[1], yt_params[2]),
                yt_params]
        
    def sin_fitting_plot(self, N=3, p0=[100, 0.1, 0], trend_type="poly",
                         savefig=None):
        plt.figure()
        plt.errorbar(self.t_vals, self.detrend(N, trend_type)[0], yerr=50,
                     color='black')
        plt.plot(self.t_vals_cont, self.sin_fitting(N, p0=[100, 0.1, 0.])[0],
                 color='red')
        plt.ylim([-200, 200])
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (km)")
        if savefig is not None:
            plt.savefig(savefig[0])
        ##
        plt.figure()
        plt.errorbar(self.t_vals, self.detrend(N, trend_type)[1], yerr=50,
                     color='black')
        plt.plot(self.t_vals_cont, self.sin_fitting(N, p0=[100, 0.1, 0.])[0],
                 color='red')
        plt.ylim([-200, 200])
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (km)")
        if savefig is not None:
            plt.savefig(savefig[1])
