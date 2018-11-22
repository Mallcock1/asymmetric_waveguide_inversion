# -*- coding: utf-8 -*-
"""
Created on Fri Nov 02 13:07:39 2018

@author: Matt

Inversion functions for the Alfven speed in a magnetic slab.
"""

import numpy as np
import scipy as sc
from scipy.optimize import fsolve


def cT(vA, c0):
    return sc.sqrt(c0**2*vA**2 / (c0**2 + vA**2))


def c1(vA, c0, R1):
    return np.sqrt(1/R1 * (c0**2 + 5./6 * vA**2))


def c2(vA, c0, R2):
    return np.sqrt(1/R2 * (c0**2 + 5./6 * vA**2))


def m0(w, k, vA, c0):
    m0function = sc.sqrt((k**2*c0**2 - w**2)*(k**2*vA**2 - w**2) /
                         ((c0**2 + vA**2)*(k**2*cT(vA, c0)**2 - w**2)))
    return m0function


def m1(w, k, vA, c0, R1):
    m1function = sc.sqrt(k**2 - w**2 / c1(vA, c0, R1)**2)
    return m1function


def m2(w, k, vA, c0, R2):
    m2function = sc.sqrt(k**2 - w**2 / c2(vA, c0, R2)**2)
    return m2function


def lamb0(w, k, vA, c0):
    return -(k**2*vA**2 - w**2)*1.j / (m0(w, k, vA, c0)*w)


def lamb1(w, k, vA, c0, R1):
    return R1*w*1.j / m1(w, k, vA, c0, R1)


def lamb2(w, k, vA, c0, R2):
    return R2*w*1.j / m2(w, k, vA, c0, R2)


error_string_kink_saus = "mode argument must be 'kink' or 'saus'"
error_string_subscript = "subscript argument must be 1 or 2"


def disp_rel_sym(w, k, vA, c0, R1, R2, x0, mode, subscript):
    if mode != "kink" and mode != "saus":
        raise ValueError(error_string_kink_saus)
    elif subscript != 1 and subscript != 2:
        raise ValueError(error_string_subscript)
    else:
        if subscript == 1:
            if mode == "kink":
                dispfunction = lamb0(w, k, vA, c0)*sc.tanh(m0(w, k, vA, c0)*x0) + lamb1(w, k, vA, c0, R1)
            elif mode == "saus":
                dispfunction = lamb0(w, k, vA, c0) + lamb1(w, k, vA, c0, R1)*sc.tanh(m0(w, k, vA, c0)*x0)
        elif subscript == 2:
            if mode == "kink":
                dispfunction = lamb0(w, k, vA, c0)*sc.tanh(m0(w, k, vA, c0)*x0) + lamb2(w, k, vA, c0, R2)
            elif mode == "saus":
                dispfunction = lamb0(w, k, vA, c0) + lamb2(w, k, vA, c0, R2)*sc.tanh(m0(w, k, vA, c0)*x0)
        return dispfunction


def amp_ratio(w, k, vA, c0, R1, R2, x0, mode):
    if mode == "kink":
        ampfunction = disp_rel_sym(w, k, vA, c0, R1, R2, x0, 'saus', 1) / disp_rel_sym(w, k, vA, c0, R1, R2, x0, 'saus', 2)
    elif mode == "saus":
        ampfunction = - disp_rel_sym(w, k, vA, c0, R1, R2, x0, 'kink', 1) / disp_rel_sym(w, k, vA, c0, R1, R2, x0, 'kink', 2)
    else:
        raise ValueError(error_string_kink_saus)
    return ampfunction


def alfven_AR_inversion(w, k, vA_guess, c0, R1, R2, x0, RA, mode):
    def inversion_function(vA):
        return np.real(amp_ratio(w, k, vA, c0, R1, R2, x0, mode) - RA)
    vA_sol = fsolve(inversion_function, vA_guess, xtol=1e-08)
    return vA_sol


def alfven_AR_inversion_2var(w, k, vA_guess, c0, R1_guess, R2, x0, RA, mode):
    def inversion_function(vA_R1):
        return np.real(amp_ratio(w, k, vA_R1[0], c0, vA_R1[1], R2, x0, mode) - RA)
    vA_sol, R1_sol = fsolve(inversion_function, [vA_guess, R1_guess], xtol=1e-08)
    return [vA_sol, R1_sol]
