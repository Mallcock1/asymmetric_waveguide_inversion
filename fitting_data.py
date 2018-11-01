# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:32:57 2018

@author: Matt
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#width of fibril in pixels
width_pix = np.array([16.25,14,11,7.5,8,7,7,5.75,7.5,8.5,9.25,20.25,11.5,11.5,
                      11.75,9.5,8.5,7.5,7.25,7.25,6.25,5.5,7,6.5,7.75,8,6.75,6,
                      7.75,9,12,12.5,11.5,9])

#y-distance from bottom of domain in pixels
yb_pix = np.array([24.75,27.25,29.5,31.5,31,31.75,32.25,33.75,34,34.5,34.75,30,
                   35.5,36,36.25,37,36.75,36.5,37,37,37.75,38.25,38,39.5,41,
                   42.5,43.5,44,43,41.25,39,38,37,36.75])

##smooth over anomalous data point
#width_pix[11] = 0.5*(width_pix[12] + width_pix[10])
#yb_pix[11] = 0.5*(yb_pix[12] + yb_pix[10])

class fibril:
    def __init__(self, yb, yt):
        self.yb = yb
        self.yt = yt      
    
    
    #convert units from pixels to km
    def pix_to_km(self):
        self.yb = self.yb * 50
        self.yt = self.yt * 50
    
    
    #smooth over anomalous data points
    def smooth(self, indices_of_anomalies):
        for i in indices_of_anomalies:
            self.yb[i] = 0.5*(self.yb[i+1] + self.yb[i-1])
            self.yt[i] = 0.5*(self.yt[i+1] + self.yt[i-1])
    
    
    #trend of the data - fit with an Nth degree polynomial with least squares regression.
    def trend(self, N=3):
        x_vals = np.arange(len(self.yb))
        polyb = np.poly1d(np.polyfit(x_vals, self.yb, N))
        polyt = np.poly1d(np.polyfit(x_vals, self.yt, N))
        return [polyb(x_vals), polyt(x_vals)]
            
    
    #subtract the trend from the data
    def detrend(self, N=3, trend_type="poly"):
        if trend_type == "poly":
            detrend_b = self.yb - self.trend(N)[0]
            detrend_t = self.yt - self.trend(N)[1]
        elif trend_type == "diff":
            detrend_b = self.yb[1:] - self.yb[:-1]
            detrend_t = self.yt[1:] - self.yt[:-1]
        else:
            raise ValueError('trend_type must be either "poly" or "diff"')
        return [detrend_b, detrend_t]
    
    
    #fit a sinusoidal curve
    def sin_fitting(self, N=3, trend_type="poly"):
        x_vals = np.arange(len(self.yb))
        def sin_func(x, a, b, c):
            return a * np.sin(b * x + c)
        yb_params, yb_params_covariance = curve_fit(sin_func, x_vals, self.detrend(N, trend_type)[0], p0=[100, 0.1, 0])
        yt_params, yt_params_covariance = curve_fit(sin_func, x_vals, self.detrend(N, trend_type)[1], p0=[100, 0.1, 0])
        x_vals_cont = np.linspace(0, len(self.yb), 1000)
        return [sin_func(x_vals_cont, yb_params[0], yb_params[1], yb_params[2]),
                sin_func(x_vals_cont, yt_params[0], yt_params[1], yt_params[2])]
            
            
yt_pix = yb_pix + width_pix

#fibril1 = fibril(yb_pix[:22], yt_pix[:22] #works well with 3rd degree polynomial
fibril1 = fibril(yb_pix[:22], yt_pix[:22])
fibril1.smooth([11])
fibril1.pix_to_km()

x_vals = np.arange(len(fibril1.yb))
x_vals_cont = np.linspace(0, len(fibril1.yb), 1000)

N=2

plt.figure()
plt.plot(fibril1.yb)
plt.plot(fibril1.yt)
plt.plot(fibril1.trend(N=N)[0])
plt.plot(fibril1.trend(N=N)[1])

plt.figure()
plt.plot(fibril1.detrend(N=N)[0])
plt.plot(fibril1.detrend(N=N)[1])

plt.figure()
plt.plot(fibril1.detrend(N=N)[0])
plt.plot(x_vals_cont, fibril1.sin_fitting(N=N)[0])

plt.figure()
plt.plot(fibril1.detrend(N=N)[1])
plt.plot(x_vals_cont, fibril1.sin_fitting(N=N)[1])
