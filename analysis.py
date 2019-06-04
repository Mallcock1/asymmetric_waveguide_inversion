# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield

Master code that calls other scripts to produce distance-distance plots,
time-distance plots, and complete seismological inversion.
"""

import numpy as np
from scipy import stats
from scipy.signal import correlate
import time_distance as td
import fibril_inversion as fi
from boundary_data import *
import matplotlib.pyplot as plt


show_time_distance_data = False
do_fibril_inversion = False
include_density_range = False
plot_density_range = False

# Uncomment the things you want
#show_time_distance_data = True
do_fibril_inversion = True
#include_density_range = True
#plot_density_range = True

############################################################
# Specify slit coordinates - place slit perpendicular to dark fibril
#
# Strategy for slit choice:
# - Perpendicular to an isolated fibril that does not radically change
#   direction during time range.
# - At least the width of the fibril on each side

#slit_coords = [[420, 424, 615.6, 565.6], [591, 580, 265, 371],
#               [775, 775, 430, 490], [712, 687, 51.2, 110]]

############################################################

if show_time_distance_data is True:
    file_path = ('D:/my_work/projects/Asymmetric_slab/'
                 'Application_to_observation/Morton_2012_data/rosa_data')

    # Specify viewing box boundary - these must be called as kwargs in Full_map
#   bottom_left=[250, 600], top_right=[400, 700]

    # Retrieve the data
    morton12 = td.Full_map(file_path, time_range=time_range)

    # Crop data to specified viewing box
    morton12.crop()

    # Animate data through time
    morton12.animate(slit_coords=slit_coords, interval=100)#,
#                     savefig="plots/ani_all.mp4")  # "plots/animation.mp4")
#    morton12.animate(slit_coords=slit_coords, interval=100,
#                     savefig="plots/fibril" + fibril_number + "_ani.mp4")  # "plots/animation.mp4")

    # Take slice of data at specified time frame
#    morton12.image(time=52, slit_coords=None, savefig=None)
#    morton12.image(52, slit_coords=slit_coords, multi_slit=True, savefig=None)

    # Create distance-time dataset along slit for the given time_range
#    morton12.distancetime(slit_coords=slit_coords, plot=True)
#                          savefig="plots/fibril" + fibril_number + "_dt.png")

    # Take a slice of the intensity across the slit at a given time
#    morton12.intensity_slice(slit_coords=slit_coords, time_slice=[130],
#                             p0=[0.1, 10., 10., -1.0], gauss_fit=True,
#                             savefig=None)

    # Find boundaries of low intensity structure using gaussian fitting
#    boundaries = morton12.find_boundaries(slit_coords=slit_coords,
#                                          moving_average=True, num_wtd_av=3,
#                                          p0=p0_gauss, stabilise=stabilise,
#                                          plot=True, savefig="plots/fibril" + fibril_number + "_dt.png")
#    multi_boundaries = morton12.find_multi_slit_boundaries(slit_coords=slit_coords,
#                                                           num_slits=5,
#                                                           slit_distance=10., moving_average=False,
#                                                           p0=[0.1, 10., 10., -1.0], plot=True,
#                                                           savefig=None)
#    
###    [0.5, 45., 10., -1.1]

    # Find widths of the fibril through time
#    widths = morton12.find_multi_slit_widths(slit_coords=slit_coords,
#                                             moving_average=True, num_wtd_av=3,
#                                             num_slits=5, slit_distance=5.,
#                                             p0=p0_gauss, stabilise=stabilise,
#                                             plot=True)#, savefig="plots/fibril" + fibril_number + "_widths.png")

#    multi_axes = morton12.find_multi_slit_axis_min_intens(slit_coords=slit_coords,
#                                             moving_average=True, num_wtd_av=3,
#                                             num_slits=3, slit_distance=5.,
#                                             p0=p0_gauss, stabilise=stabilise,
#                                             plot=True)
    
    # Testing cross-correlation time lag between slit width signals to
    # estimate phase speeds
#    lag = np.argmax(correlate(widths[0][0], widths[-1][0]))
#
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
    # Create object of fibril boundary data
    fibril = fi.Fibril(yb, yt, t_vals, trend_range=trend_range)

    # Smooth over anomolous datapoints
    fibril.smooth(smooth_indices)

    # Convert units
    if unit == "pix":
        fibril.pix_to_km()
    elif unit != "km":
        raise ValueError('unit must be "pix" or "km" in boundary_data.py')

    # Fit sinusoids to top and bottom boundary data
    sin_fit = fibril.sin_fitting(N=N, p0=p0)#, plot=True)#,
#                                 savefig=["plots/fibril" + fibril_number + "_detrend_b.png", "plots/fibril" + fibril_number + "_detrend_t.png"])

    # Plot the trend
#    fibril.trend_plot(N, savefig="plots/fibril" + fibril_number + "_trend.png")

    # Alfven speed inversion for 100 initial values to check for consistency
    # Number of initial vA values tried
    N_vA_init = 100

    # Range of initial vA values
    vA_init_min = 1
    vA_init_max = 100
    vA_inits = np.linspace(vA_init_min, vA_init_max, N_vA_init)

    if include_density_range is False:
        vA_inversion = fibril.AR_inversion_multiple_init(p0, N, vA_inits,
                                                         c_phase, c0, R1, R2,
                                                         mode)
        print('Estimated vA in fibril ' + fibril_number + ' is: ' +
              str(vA_inversion[0]) + '.\nThis value was output for ' +
              str(vA_inversion[1]) + ' of the ' + str(N_vA_init) + ' initial' +
                  ' values.')

    if include_density_range is True:
        print('A range of densities will be considered.')

        R_list = [R1, R2]
        # Find the min and max ratio
        R_min = np.min(R_list)
        R_max = np.max(R_list)
        R_min_index = np.argmin(R_list)
        R_max_index = np.argmax(R_list)
        R_diff = R_max - R_min

        print('R' + str(R_max_index + 1) + ' is the largest density ratio ' +
              'so will be varied, holding R' + str(R_min_index + 1) +
              ' constant.')

        # Number of density ratio values tried
        N_R = 10
        # Range of denisity ratio values
        R_range_min = R_max - 0.9*R_diff
        R_range_max = R_max + 5*R_diff
        R_range = np.linspace(R_range_min, R_range_max, N_R)

        # Initialise empty array
        vA_inversion = np.empty((N_R, 2))
        
        if R_max_index == 0:
            for i, R in enumerate(R_range):
                # Calculate inversion ratios R, R2
                vA_inv_i = fibril.AR_inversion_multiple_init(p0, N,vA_inits, c_phase, c0, R, R2, mode)
                vA_inversion[i, :] = vA_inv_i
                print('\nEstimated vA in fibril ' + fibril_number + ' at [R1, R2]=[' +
                      str(round(R, 2)) + ', ' + str(R2) + '] is: ' +
                      str(vA_inversion[i, 0]) + '.\nThis value was output for ' +
                      str(vA_inversion[i, 1]) + ' of the ' + str(N_vA_init) + ' initial' +
                      ' values.')
        
        if R_max_index == 1:
            for i, R in enumerate(R_range):
                # Calculate inversion ratios R1, R
                vA_inv_i = fibril.AR_inversion_multiple_init(p0, N,vA_inits, c_phase, c0, R1, R, mode)
                vA_inversion[i, :] = vA_inv_i
                print('\nEstimated vA in fibril ' + fibril_number + ' at [R1, R2]=[' +
                      str(R1) + ', ' + str(round(R, 2)) + '] is: ' +
                      str(vA_inversion[i, 0]) + '.\nThis value was output for ' +
                      str(vA_inversion[i, 1]) + ' of the ' + str(N_vA_init) + ' initial' +
                      ' values.')
        
        if plot_density_range is True:
            plt.plot(R_range, vA_inversion[:, 0])
