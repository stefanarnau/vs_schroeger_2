#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 16:44:23 2022

@author: plkn
"""

# Imports
import numpy as np
import os
import scipy.io
import sklearn.preprocessing
import sklearn.linear_model
import matplotlib.pyplot as plt

# Paths
path_in = '/mnt/data_dump/schroeger2/results/'

# Regression model and scaler
linreg = sklearn.linear_model.Ridge()
scaler = sklearn.preprocessing.StandardScaler()

# Load erp data
fn = os.path.join(path_in, "erp_data.mat")
erp_data = scipy.io.loadmat(fn)["erp_data"]
erp_times = np.squeeze(scipy.io.loadmat(fn)["erp_times"])
X = scipy.io.loadmat(fn)["ages"]

# Conditions:
# 1: std short
# 2: std long
# 3: dev short
# 4: dev long

# Calculate deviant - standard
erp_std = erp_data[:, (0, 1), :, :].mean(axis=1)
erp_dev = erp_data[:, (2, 3), :, :].mean(axis=1)
erp_diff = erp_dev - erp_std

# Get dimensions
n_subjects, n_channels, n_times = erp_diff.shape

# Reshape
y = np.abs(erp_diff).reshape((n_subjects, -1))

# Scale
X = scaler.fit_transform(X)
y = scaler.fit_transform(y)

# Fit regression model and get coefs
coef_true = linreg.fit(X, y).coef_.reshape((n_channels, n_times))
coef_fake = linreg.fit(np.random.permutation(X), y).coef_.reshape((n_channels, n_times))

# Plotting
maxval = 0.3
plt.contourf(erp_times, np.arange(n_channels) + 1, coef_true, vmin=-maxval, vmax=maxval, cmap="jet")
