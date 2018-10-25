# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:32:57 2018

@author: Matt
"""

import numpy as np
import matplotlib.pyplot as plt

#width of fibril in pixels
width_pix = np.array([16.25,14,11,7.5,8,7,7,5.75,7.5,8.5,9.25,20.25,11.5,11.5,
                      11.75,9.5,8.5,7.5,7.25,7.25,6.25,5.5,7,6.5,7.75,8,6.75,6,
                      7.75,9,12,12.5,11.5,9])

#y-distance from bottom of domain in pixels
yb_pix = np.array([24.75,27.25,29.5,31.5,31,31.75,32.25,33.75,34,34.5,34.75,30,
                   35.5,36,36.25,37,36.75,36.5,37,37,37.75,38.25,38,39.5,41,
                   42.5,43.5,44,43,41.25,39,38,37,36.75])

#smooth over anomalous data point
width_pix[11] = 0.5*(width_pix[12] + width_pix[10])
yb_pix[11] = 0.5*(yb_pix[12] + yb_pix[10])


class fibril:
    def __init__(self, width, yb):
        self.width = width
        self.yb = yb
    
    #smooth over anomalous data points
    def smooth(self, indices_of_anomalies):
        for i in indices_of_anomalies:
            self.width[i] = 0.5*(self.width[i+1] + self.width[i-1])
            self.yb[i] = 0.5*(self.yb[i+1] + self.yb[i-1])




fibril1 = fibril(width_pix, yb_pix)

fibril1.smooth([11])

