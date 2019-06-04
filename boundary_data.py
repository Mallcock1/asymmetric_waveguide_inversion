# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield

Fibril boundary data from Morton et al 2012 ROSA data

BEWARE: time data might not match exactly with boundary data. This was caused
by migrating data to csv files and is yet to be fixed.
"""

import numpy as np
import pandas as pd

fibril_1 = False
fibril_2 = False
fibril_3 = False
fibril_4 = False
fibril_5 = False

fibril_1 = True
#fibril_2 = True
#fibril_3 = True
#fibril_4 = True
#fibril_5 = True

# Parameters for all fibrils:
# Temporal resolution of the ROSA instrument in s
cad = 7.68
# Degree of the trend polynomial
N = 3
# Initial parameters for Gaussian fitting
p0_gauss = [0.1, None, 10., -1.0]


def stretch_slit(slit_coords, factor):
    """
    Stretch a given slit by a factor
    """

    x0, x1, y0, y1 = slit_coords

    xmid = (x1 + x0)/2
    ymid = (y1 + y0)/2

    x0new = (x0 - xmid)*factor + xmid
    x1new = (x1 - xmid)*factor + xmid
    y0new = (y0 - ymid)*factor + ymid
    y1new = (y1 - ymid)*factor + ymid
    return [x0new, x1new, y0new, y1new]


def shift_slit(slit_coords, shift):
    """
    Shift a given slit
    """

    x0, x1, y0, y1 = slit_coords

    # gradient
    m = (y1 - y0)/(x1 - x0)
    # x and y shift
    sx = shift/np.sqrt(1 + m**2)
    sy = m*shift/np.sqrt(1 + m**2)

    x0new = x0 + sx
    x1new = x1 + sx
    y0new = y0 + sy
    y1new = y1 + sy
    return [x0new, x1new, y0new, y1new]


############################################################
# Fibril 1


if fibril_1 is True:
    print("You have chosen fibril_1")
    fibril_number = "1"

    # Read boundary data from csv file
    data = pd.read_csv('fibril_boundary_data_1.csv')

    # Distances of top and bottom boundary from bottom of slit in km
    yt = np.array(data['Top boundary'])
    yb = np.array(data['Bottom boundary'])
    t_vals = np.array(data['Time'])

    unit = "km"

    slit_coords = [420, 424, 615.6, 565.6]  # This is what the boundary data above is from: [425, 429, 616, 566] (actually it's not)

    smooth_indices = []

    # time range in frames
    t_start = 151
    t_end = 225
    time_range = slice(t_start, t_end)

    trend_range = None

    stabilise = 10

    p0 = [100., 0.01, 0.]

    vA_guess = 100.  # estimate from morton 12
    c0 = 10.  # estimate given in Morton 12
    R1 = 0.1  # R1 := rho_1 / rho_0
    R2 = 0.2  # R2 := rho_2 / rho_0
    c_phase = 63.  # Morton found this: 71.
    mode = "kink"


###########################################################
# Fibril 2

if fibril_2 is True:
    print("You have chosen fibril_2")
    fibril_number = "2"

    # Read boundary data from csv file
    data = pd.read_csv('fibril_boundary_data_2.csv')

    # Distances of top and bottom boundary from bottom of slit in km
    yt = np.array(data['Top boundary'])
    yb = np.array(data['Bottom boundary'])
    t_vals = np.array(data['Time'])

    unit = "km"

    smooth_indices = [24]

    slit_coords = [591, 580, 265, 371]

    t_start = 0
    t_end = 71
    time_range = slice(t_start, t_end)
    trend_range = None

    stabilise = False

    # initial values for sin
    p0 = [100., 0.03, 0.]

    vA_guess = 100.  # estimate from morton 12
    c0 = 10.  # estimate given in Morton 12
    R1 = 0.2  # R1 := rho_1 / rho_0
    R2 = 0.1  # R2 := rho_2 / rho_0
    c_phase = 63.
    mode = "saus"


###########################################################
# Fibril 3

if fibril_3 is True:
    print("You have chosen fibril_3")
    fibril_number = "3"

    # Read boundary data from csv file
    data = pd.read_csv('fibril_boundary_data_3.csv')

    # Distances of top and bottom boundary from bottom of slit in km
    yt = np.array(data['Top boundary'])
    yb = np.array(data['Bottom boundary'])
    t_vals = np.array(data['Time'])

    unit = "km"

    smooth_indices = []

    slit_coords = [775, 775, 430, 490]

    t_start = 132
    t_end = 175
    time_range = slice(t_start, t_end)
    trend_range = None

    p0 = [100., 0.08, 0.]

    vA_guess = 100.  # estimate from morton 12
    c0 = 10.  # estimate given in Morton 12
    R1 = 0.1  # R1 := rho_1 / rho_0
    R2 = 0.2  # R2 := rho_2 / rho_0
    c_phase = 63.
    mode = "saus"


###########################################################
# Fibril 4

if fibril_4 is True:
    print("You have chosen fibril_4")
    fibril_number = "4"

    # Read boundary data from csv file
    data = pd.read_csv('fibril_boundary_data_4.csv')

    # Distances of top and bottom boundary from bottom of slit in km
    yt = np.array(data['Top boundary'])
    yb = np.array(data['Bottom boundary'])
    t_vals = np.array(data['Time'])

    unit = "km"

    smooth_indices = [[22, 23, 24, 25, 26, 27]]

    slit_coords = [712, 687, 51.2, 110]
    slit_coords = stretch_slit(slit_coords, 0.6)

    t_start = 138
    t_end = 238
    time_range = slice(t_start, t_end)
    trend_range = None  # slice(29, -1)

    stabilise = False

    p0 = [[100., 0.01, 0.], [100., 0.01, 0.]]

    vA_guess = 100.  # estimate from morton 12
    c0 = 10.  # estimate given in Morton 12
    R1 = 0.1  # R1 := rho_1 / rho_0
    R2 = 0.2  # R2 := rho_2 / rho_0
    c_phase = 31.
    mode = "saus"


###########################################################
# Fibril 5

if fibril_5 is True:
    print("You have chosen fibril_5")
    fibril_number = "5"

    # Read boundary data from csv file
    data = pd.read_csv('fibril_boundary_data_5.csv')

    # Distances of top and bottom boundary from bottom of slit in km
    yt = np.array(data['Top boundary'])
    yb = np.array(data['Bottom boundary'])
    t_vals = np.array(data['Time'])

    unit = "km"

    smooth_indices = []

    slit_coords = [743, 709, 143, 339]
    slit_coords = stretch_slit(slit_coords, 0.2)
    slit_coords = shift_slit(slit_coords, -40)

    t_start = 16
    t_end = 55
    time_range = slice(t_start, t_end)
    trend_range = None

    stabilise = 5

    p0 = [100., 0.01, 0.]

    vA_guess = 100.  # estimate from morton 12
    c0 = 10.  # estimate given in Morton 12
    R1 = 0.1  # R1 := rho_1 / rho_0
    R2 = 0.2  # R2 := rho_2 / rho_0
    c_phase = 129.
    mode = "kink"
