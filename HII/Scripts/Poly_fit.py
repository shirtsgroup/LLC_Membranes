"""
An Adaption of Homework 10, Problem 1, Part 1 from Numerical Analysis, Spring 2016

The general least squares method:

(1) Create a 'Z' matrix -->  a matrix which can be seen as transposed and concatenated lists of each value in a list of
                             independent variables such that the first column is x_i^0, the second column is x_i^1 etc.
                             all the way up to x_i^n, where n is the degree of the polynomial being fit
(2) Calculate 'A', a vector containing the coefficients used for the polynomial fit of the form:
                            a0 + a1*x + a2*x**2 + ... + an*x**n
                            A = [a0, a1, a2, ..., an ]
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sci
import math

# x = [1.01, 1.27, 1.85, 2.38, 2.83, 3.13, 3.96, 4.91]
# y = [0.00, 0.19, 0.58, 0.96, 1.26, 1.47, 2.07, 2.75]


def z_matrix(degree, t):
    Z = np.zeros((len(t), degree + 1))
    for i in range(0, degree + 1):
        for j in range(0, len(t)):
            Z[j, i] = t[j]**i
    return Z


def find_A(x, Y, degree, W = 'none'):  # Find the matrix of coefficients which give desired fit
    # if A is a matrix of coefficients ([a0, a1, ..., an]), then
    # A = ([Z]^T*[Z])^-1*[Z]^T*Y
    # With a weight matrix:
    # A = ([Z]^T*W*[Z])^-1*[Z]^T*W*Y
    # The weight matrix is a square n x n matrix (where n = len(x)) where the diagonals are the inverse of the standard
    # deviation of the data point squared (1/(sigma_i ** 2))
    Z = z_matrix(degree, x)
    ZT = np.transpose(Z)  # transposed z matrix
    if isinstance(W, np.ndarray):
        ZTW = np.dot(ZT, W)
        ZTWZ = np.dot(ZTW, Z)
        ZTWZinv = np.linalg.inv(ZTWZ)
        ZTWZinvZT = np.dot(ZTWZinv, ZT)
        ZTWZinvZTW = np.dot(ZTWZinvZT, W)
        A = np.array(np.dot(ZTWZinvZTW, Y).tolist())
    else:
        ZTZ = np.dot(ZT, Z)  # multiply ZT and Z
        ZTZinv = np.linalg.inv(ZTZ)  # inverse of ZTZ
        ZTZinvZT = np.dot(ZTZinv, ZT)  # multiply the inverse of Z transposed * Z by Z transposed
        A = np.array(np.dot(ZTZinvZT, Y).tolist())
    return A


def poly_eval(degree, x, A):
    y = A[0]
    for i in range(1, degree + 1):
        y += A[i]*x**i
    Sr = (y)
    return y


def poly_fit(x, y, degree, weight='none'):

    if isinstance(weight, np.ndarray):
        A = find_A(x, y, degree, weight)
    else:
        A = find_A(x, y, degree)

    y_val = np.zeros((len(x)))
    St = 0
    Sr = 0
    mean = np.mean(y)
    n = len(y)
    for i in range(0, len(x)):
        y_val[i] = poly_eval(degree, x[i], A)
        Sr += (y[i] - (y_val[i]))**2
        St += (y[i] - mean)**2
    s = math.sqrt(Sr/(len(y) - 2))
    x_square = [i**2 for i in x]
    slope_error = s * np.sqrt(n/(n*sum(x_square) - (sum(x))**2))
    intercept_error = s * np.sqrt(sum(x_square)/(n*sum(x_square) - (sum(x))**2))
    R_squared = 1 - (Sr/St)

    return y_val, R_squared, slope_error, intercept_error, A