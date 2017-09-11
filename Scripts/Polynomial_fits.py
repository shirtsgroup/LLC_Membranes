#!/usr/bin/python

import numpy
import matplotlib.pyplot as plt
import math

# fit a line to the data

def find_A(Z, y):  # Find the matrix of coefficient which give desired fit
    # if A is a matrix of coefficients ([a0, a1, ..., am]), then
    # A = ([Z]^T*[Z])^-1*[Z]^T*Y
    ZT = np.transpose(Z)  # transposed z matrix
    ZTZ = np.dot(ZT, Z)  # multiply ZT and Z
    ZTZinv = np.linalg.inv(ZTZ)  # inverse of ZTZ
    ZTZinvZT = np.dot(ZTZinv, ZT)  # multiply the inverse of Z transposed * Z by Z transposed
    A = np.array(np.dot(ZTZinvZT, y).tolist())
    return A

x = [x for x in frames]  # looking at a chloride concentration of zero
f0 = np.ones((len(frames), 1))
xt = np.transpose(np.matrix(x))
x2 = [x**2 for x in frames]  # adding a quadratic term to the Z matrix
x2t = np.transpose(np.matrix(x2))
x3 = [x**3 for x in frames]  # adding a quadratic term to the Z matrix
x3t = np.transpose(np.matrix(x3))

Z_cub = np.concatenate((f0, xt, x2t, x3t), axis=1)  # used f0 and f1 for a linear fit

A_cub = find_A(Z_cub, out_of_bounds)

def cubic_fit(x, A_cub, y):
    y_val = np.zeros((len(x)))
    Sr = 0
    St = 0
    mean = np.mean(y)
    for i in range(0, len(x)):
        y_val[i] = A_cub[0, 3]*x[i]**3 + A_cub[0, 2]*x[i]**2 + A_cub[0, 1]*x[i] + A_cub[0, 0]
        Sr += (y[i] - (A_cub[0, 3]*x[i]**3 + A_cub[0, 2]*x[i]**2 + A_cub[0, 1]*x[i] + A_cub[0, 0]))**2
        St += (y[i] - mean)**2
    s = math.sqrt(Sr / (len(y) - 4))
    R_squared = 1 - (Sr/St)
    return y_val, R_squared, s

y_cub, R_squared, s = cubic_fit(frames, A_cub, out_of_bounds)