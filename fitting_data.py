# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:32:57 2018

@author: Matt
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import alfven_speed_inversion as asi

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

cad = 7.68 #temporal cadence of the ROSA instrument in s
t_start = 1000 + 19*cad

class fibril:
    def __init__(self, yb, yt):
        self.yb = yb
        self.yt = yt      
        self.t_vals = t_start + np.arange(len(self.yb)) * cad * 2 #because there is one yellow bar per 2 pixels
        self.t_vals_cont = np.linspace(t_start, self.t_vals[-1], 1000) 
    
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
        
    
    def width(self, N=3):
        return np.average(self.trend(N)[1] - self.trend(N)[0])
            
    
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
        def sin_func(x, a, b, c): #a = amplitude, b = angular frequency
            return a * np.sin(b * x + c)
        yb_params, yb_params_covariance = curve_fit(sin_func, self.t_vals, self.detrend(N, trend_type)[0], p0=[100, 0.01, 0])
        yt_params, yt_params_covariance = curve_fit(sin_func, self.t_vals, self.detrend(N, trend_type)[1], p0=[100, 0.01, 0])
        return [sin_func(self.t_vals_cont, yb_params[0], yb_params[1], yb_params[2]), yb_params,
                sin_func(self.t_vals_cont, yt_params[0], yt_params[1], yt_params[2]), yt_params]
            
            
yt_pix = yb_pix + width_pix

#fibril1 = fibril(yb_pix[:22], yt_pix[:22] #works well with 2nd degree polynomial
fibril1 = fibril(yb_pix[:22], yt_pix[:22])
fibril1.smooth([11])
fibril1.pix_to_km()

fibril1_full = fibril(yb_pix, yt_pix)
fibril1_full.smooth([11])
fibril1_full.pix_to_km()

N=2

sin_fit = fibril1.sin_fitting(N=N)

#plt.figure()
#plt.errorbar(fibril1_full.t_vals, fibril1_full.yb, yerr=50, color='black')
#plt.errorbar(fibril1_full.t_vals, fibril1_full.yt, yerr=50, color='black')
#plt.plot(fibril1.t_vals, fibril1.trend(N=N)[0], color='red')
#plt.plot(fibril1.t_vals, fibril1.trend(N=N)[1], color='red')
#plt.xlabel("Time (s)")
#plt.ylabel("Distance (km)")
#plt.savefig("fibril1_full.png")

#plt.figure()
#plt.plot(fibril1.t_vals, fibril1.detrend(N=N)[0])
#plt.plot(fibril1.t_vals, fibril1.detrend(N=N)[1])
#

#
plt.figure()
plt.errorbar(fibril1.t_vals, fibril1.detrend(N=N)[0], yerr=50, color='black')
#plt.plot(fibril1.t_vals, fibril1.detrend(N=N)[0], color='black')
plt.plot(fibril1.t_vals_cont, sin_fit[0], color='red')
#plt.hlines(0, 1100, 1500, colors='red', linestyles='dashed')
#plt.hlines(abs(sin_fit[1][0]), 1100, 1500, colors='red', linestyles='dashed')
plt.ylim([-200,200])
plt.xlabel("Time (s)")
plt.ylabel("Distance (km)")
plt.savefig("fibril1_detrend_b.png")
##
plt.figure()
plt.errorbar(fibril1.t_vals, fibril1.detrend(N=N)[1], yerr=50, color='black')
#plt.plot(fibril1.t_vals, fibril1.detrend(N=N)[1], color='black')
plt.plot(fibril1.t_vals_cont, sin_fit[2], color='red')
#plt.hlines(0, 1100, 1500, colors='red', linestyles='dashed')
#plt.hlines(abs(sin_fit[3][0]), 1100, 1500, colors='red', linestyles='dashed')
plt.ylim([-200,200])
plt.xlabel("Time (s)")
plt.ylabel("Distance (km)")
plt.savefig("fibril1_detrend_t.png")

print("\n" + "TOP: Amp = " + str("%.4g" % abs(sin_fit[3][0])) + "km, Freq = " + str("%.4g" % sin_fit[3][1]) + "s-1"
+ "\n\n" + "BOTTOM: Amp = " + str("%.4g" % abs(sin_fit[1][0])) + "km,  Freq = " + str("%.4g" % sin_fit[1][1]) + "s-1" + "\n")


##########################################################

w = (sin_fit[1][1] + sin_fit[3][1]) / 2 #freq
k = w / 71. #87. #phase speed estimate from morton 12 supp fig S5
vA_guess = 100. #estimate from morton 12
c0 = 10. # estimate given in Morton 12
R1 = 0.2 #R1 := rho_1 / rho_0
R2 = 0.1 #R2 := rho_2 / rho_0
x0 = fibril1.width(N) / 2
RA = sin_fit[3][0] / sin_fit[1][0]
mode = "saus"

vA_sol = asi.alfven_AR_inversion(w, k, vA_guess, c0, R1, R2, x0, RA, mode)
print("vA ~ " + str(vA_sol))

print("RA guess ~ " + str(np.real(asi.amp_ratio(w, k, vA_guess, c0, R1, R2, x0, mode))))