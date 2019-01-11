# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield
"""

import numpy as np
from scipy.signal import correlate
import time_distance as td
import fibril_inversion as fi
from boundary_data import *

show_time_distance_data = False
do_fibril_inversion = False

show_time_distance_data = True
#do_fibril_inversion = True

############################################################



if show_time_distance_data is True:
    file_path = ('D:/my_work/projects/Asymmetric_slab/'
                 'Application_to_observation/Morton_2012_data/rosa_data')
    slit_coords = [[420, 424, 615.6, 565.6], [591, 580, 265, 371],
                   [775, 775, 430, 490], [712, 687, 51.2, 110]]
    
    morton12 = td.Full_map(file_path, time_range=time_range)
    morton12.crop([0, 0], [720, 780])
    morton12.animate(slit_coords=slit_coords, interval=100)#,
#                     savefig="plots/ani_all.mp4")  # "plots/animation.mp4")
    
#    morton12.animate(slit_coords=slit_coords, interval=100,
#                     savefig="plots/" + fibril_number + "_ani.mp4")  # "plots/animation.mp4")
#    morton12.distancetime(slit_coords=slit_coords, plot=True)
#                          savefig="plots/" + fibril_number + "_dt.png")
#    morton12.intensity_slice(slit_coords=slit_coords, time_slice=[130],
##                             p0=[0.1, 10., 10., -1.0], gauss_fit=True, savefig=None)
#    boundaries = morton12.find_boundaries(slit_coords=slit_coords,
#                                          moving_average=True, num_wtd_av=3,
#                                          p0=p0_gauss, stabilise=stabilise,
#                                          plot=True, savefig="plots/" + fibril_number + "_dt.png")
##    multi_boundaries = morton12.find_multi_slit_boundaries(slit_coords=slit_coords,
##                                                           num_slits=5,
##                                                           slit_distance=10., moving_average=False,
##                                                           p0=[0.1, 10., 10., -1.0], plot=True,
##                                                           savefig=None)
#    
###    [0.5, 45., 10., -1.1]
#    widths = morton12.find_multi_slit_widths(slit_coords=slit_coords,
#                                             moving_average=True, num_wtd_av=3,
#                                             num_slits=3, slit_distance=5.,
#                                             p0=p0_gauss, stabilise=stabilise,
#                                             plot=True)#, savefig="plots/fibril1_widths.png")
#    
#    multi_axes = morton12.find_multi_slit_axis_min_intens(slit_coords=slit_coords,
#                                             moving_average=True, num_wtd_av=3,
#                                             num_slits=3, slit_distance=5.,
#                                             p0=p0_gauss, stabilise=stabilise,
#                                             plot=True)
    
    # Testing cross-correlation time lag between slit width signals
#    lag = np.argmax(correlate(widths[0][0], widths[-1][0]))
    
#    from scipy.interpolate import CubicSpline
#    a = 1305
#    b = 1720
#    N = 100
#    t_vals = np.linspace(a, b, N)
#    
#    cs0 = CubicSpline(widths[0][1], widths[0][0])(t_vals)
#    cs1 = CubicSpline(widths[-2][1], widths[-2][0])(t_vals)
#    
#    cor = correlate(cs0, cs1) #correlate(cs0(t_vals), cs1(t_vals))
#    t_new = np.arange(1-N, N)
#    lag = t_new[np.argmax(cor)]*(b - a)/N
#    plt.figure()
#    plt.plot(t_vals, cs0)
#    plt.plot(t_vals, cs1)
#    plt.figure()
#    plt.plot(cor)
#    phase_speed = 4*50 / lag
#    print("lag is ", lag)
#    print("phase speed is ", phase_speed)

if do_fibril_inversion is True:
    fibril = fi.Fibril(yb, yt, t_vals, trend_range=trend_range)
    fibril.smooth(smooth_indices)
    if unit == "pix":
        fibril.pix_to_km()
    elif unit != "km":
        raise ValueError('unit must be "pix" or "km" in boundary_data.py')

    sin_fit = fibril.sin_fitting(N=N, p0=p0, plot=True)#, savefig=["plots/fibril3_detrend_b.png", "plots/fibril3_detrend_t.png"])
    fibril.trend_plot(N)#, savefig="plots/fibril2_trend.png")
    
    for vA_guess in range(1, 100):
        print(fibril.AR_inversion(p0, N, vA_guess, c_phase, c0, R1, R2, mode)[0])
