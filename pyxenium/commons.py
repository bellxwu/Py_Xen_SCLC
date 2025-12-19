#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 13:33:11 2025

Description: Potentially useful functions for preprocessing

@author: bellwu
"""
# %% setup environment 
import numpy as np
import pandas as pd
# %% MP edge function
def mp_dist(T, N, pts, var=1.0):
    '''
    Marchenko–Pastur bulk PDF for eigenvalues of S = (1/N) X X^T
    where X is T x N and q = T/N.

    Args:
        T (int): number of rows of X
        N (int): number of columns of X
        pts (int): number of grid points for the PDF curve
        var (float): noise variance sigma^2

    Returns:
        eMin (float): lambda- bulk edge
        eMax (float): lambda+ bulk edge
        pdf (pd.Series): MP density values indexed by eigenvalue grid
    '''
    q = T / N

    eMin = var * (1 - np.sqrt(q))**2
    eMax = var * (1 + np.sqrt(q))**2

    eVal = np.linspace(eMin, eMax, pts)

    rad = (eMax - eVal) * (eVal - eMin)
    rad = np.maximum(rad, 0.0)  # numerical safety at edges

    pdf = (1.0 / (2 * np.pi * var * q * eVal)) * np.sqrt(rad)
    pdf = pd.Series(pdf, index=eVal)
    return eMin, eMax, pdf
# %%