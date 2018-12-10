# -*- coding: utf-8 -*-
"""
Matthew Allcock
SP2RC, University of Sheffield

fit gaussian function to data
"""

import numpy as np
from scipy.optimize import curve_fit


def gauss_fit(data, p0=[0.5, 45., 10., -1.1], retrn="params"):
    """
    Fits a Gaussian curve to data

    Inputs:
        data = numpy array of very approximately Gaussian data,
        p0 = initial [amplitude, mean, standard deviation, offset],
        retrn = 'params'       ---> return list of least squares params,
              = 'func'         ---> return fitted Gaussian function,
              = 'pos_half_max' ---> return x positions where Gaussian reaches
                                    half maximum.
    """
    # s is the coordinate along the slit
    s_vals = np.arange(len(data))

    # Gaussian curve
    def gauss_func(x, A, x0, sigma, offset):
        return np.abs(A)*np.exp(-(x - x0)**2 / (2*sigma**2)) + offset

    params, params_covariance = curve_fit(gauss_func, s_vals, data, p0=p0)

    def gauss_func_new(x):
        return gauss_func(x, params[0], params[1], params[2], params[3])

    if retrn == "params":
        return params
    elif retrn == "pos_half_max":
        # width at half maximum is np.sqrt(2*np.log(2)) * sigma
        return [[params[1] - np.sqrt(2*np.log(2))*params[2],
                params[1] + np.sqrt(2*np.log(2))*params[2]],
                [params[3] + 0.5*abs(params[0]),
                 params[3] + 0.5*abs(params[0])]]
    elif retrn == "func":
        return gauss_func_new
    else:
        raise ValueError("retrn must be 'params', 'func', or 'pos_half_max'")
