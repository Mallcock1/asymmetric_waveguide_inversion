# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield
"""

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

    morton12 = td.Full_map(file_path, time_range=time_range)
    morton12.crop([0, 0], [720, 780])
#    morton12.animate(slit_coords=slit_coords, interval=100,
#                     savefig="plots/" + fibril_number + "_ani.mp4")  # "plots/animation.mp4")
#    morton12.distancetime(slit_coords=slit_coords, plot=True)
#                          savefig="plots/" + fibril_number + "_dt.png")
#    morton12.intensity_slice(slit_coords=slit_coords, time_slice=list(range(119, 129)),
#                             p0=[0.1, 10., 10., -1.0], gauss_fit=True, savefig=None)
    boundaries = morton12.find_boundaries(slit_coords=slit_coords,
                                          moving_average=True, number_wtd_av=0,
                                          p0=[0.1, 10., 10., -1.0],
                                          plot=True, savefig=None)

if do_fibril_inversion is True:
    fibril = fi.Fibril(yb_pix, yt_pix, t_vals, trend_range=trend_range)
    fibril.smooth(smooth_indices)
    if unit == "pix":
        fibril.pix_to_km()
    elif unit != "km":
        raise ValueError('unit must be "pix" or "km" in boundary_data.py')

    sin_fit = fibril.sin_fitting(N=N, p0=p0, plot=True)
    fibril.trend_plot(N)

#    fibril.AR_inversion(p0, N, vA_guess, c_phase, c0, R1, R2, mode)
