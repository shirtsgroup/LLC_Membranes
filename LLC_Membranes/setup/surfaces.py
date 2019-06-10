#!/usr/bin/env python

import numpy as np


def SchwarzD(x, period):
    """
    :param x: a vector of coordinates (x1, x2, x3)
    :param period: length of one period

    :return: An approximation of the Schwarz D "Diamond" infinite periodic minimal surface
    """

    n = 2*np.pi / period

    a = np.sin(n*x[0])*np.sin(n*x[1])*np.sin(n*x[2])
    b = np.sin(n*x[0])*np.cos(n*x[1])*np.cos(n*x[2])
    c = np.cos(n*x[0])*np.sin(n*x[1])*np.cos(n*x[2])
    d = np.cos(n*x[0])*np.cos(n*x[1])*np.sin(n*x[2])

    return a + b + c + d


def gyroid(x, period):

    n = 2*np.pi / period

    a = np.sin(n*x[0])*np.cos(n*x[1])
    b = np.sin(n*x[1])*np.cos(n*x[2])
    c = np.sin(n*x[2])*np.cos(n*x[0])

    return a + b + c


def gridgen(surf, low, high, n, tol=0.05):
    """ Generate an n x n x n grid and reduce it to points that lie close to an implicit surface defined by `surf`

    :param surf: name of implicit surface to approximate
    :param low: lowest coordinate on each axis of grid
    :param high: highest coordinate on each axis of grid (if low is 0, then high is the length of the box vector)
    :param n: number of grid points in each dimension
    :param tol: surface--point distance tolerance. Anything within this cutoff can be used to approximate the surface.

    :type surf: str
    :type low: float
    :type high: float
    :type n: int
    :type tol: float
    """

    # make a cubic grid
    bin_size = (high - low) / n
    x = np.linspace(low, high - bin_size, n)
    y = np.linspace(low, high - bin_size, n)
    z = np.linspace(low, high - bin_size, n)

    if surf.lower() == 'ia3d' or surf.lower() == 'gyroid':

        gyro = gyroid([x[:, None, None], y[None, :, None], z[None, None, :]], high - low)

        gyro_eval = np.zeros([n**3, 3])

        count_gyro = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if -tol < gyro[i, j, k] < tol:
                        gyro_eval[count_gyro, :] = [x[i], y[j], z[k]]
                        count_gyro += 1

        grid = gyro_eval[:count_gyro, :]

    elif surf.lower() == 'pn3m' or surf.lower() == 'schwarzd':

        schwarz = SchwarzD([x[:, None, None], y[None, :, None], z[None, None, :]], high - low)
        schwarz_eval = np.zeros([n**3, 3])
        count_schwarz = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if -tol < schwarz[i, j, k] < tol:
                        schwarz_eval[count_schwarz, :] = [x[i], y[j], z[k]]
                        count_schwarz += 1

        grid = schwarz_eval[:count_schwarz, :]
    else:
        print('The phase you selected is not defined (yet)')
        exit()

    return grid


def gradient(v, surf, period):
    """ Calclate gradient vector, which is normal to the surface at x

    :param v: vector of x, y, z coordinates
    :param surf: which implicit surface is being used to approximate the structure of this phase
    :param period: inverse frequency by which unit cell repeats itself (equal to the box length for 1 period boxes)

    :type v: list or np.ndarray
    :type surf: str
    :type period: float
    """

    x = v[0]
    y = v[1]
    z = v[2]

    n = 2*np.pi / period

    if surf == 'Ia3d' or surf == 'gyroid' or surf == 'ia3d' or surf == 'SchwarzD':

        a = n*np.cos(n*x)*np.cos(n*y) - n*np.sin(n*x)*np.sin(n*z)
        b = -n*np.sin(n*y)*np.sin(n*x) + n*np.cos(n*y)*np.cos(n*z)
        c = -n*np.sin(n*y)*np.sin(n*z) + n*np.cos(n*z)*np.cos(n*x)

    elif surf == 'Pn3m' or surf == 'pn3m':

        a = n*np.cos(n*x)*np.sin(n*y)*np.sin(n*z) + n*np.cos(n*x)*np.cos(n*y)*np.cos(n*z) - n*np.sin(n*x)*np.sin(n*y)*np.cos(n*z) - n*np.sin(n*x)*np.cos(n*y)*np.sin(n*z)
        b = n*np.sin(n*x)*np.cos(n*y)*np.sin(n*z) - n*np.sin(n*x)*np.sin(n*y)*np.cos(n*z) + n*np.cos(n*x)*np.cos(n*y)*np.cos(n*z) - n*np.cos(n*x)*np.sin(n*y)*np.sin(n*z)
        c = n*np.sin(n*x)*np.sin(n*y)*np.cos(n*z) - n*np.sin(n*x)*np.cos(n*y)*np.sin(n*z) - n*np.cos(n*x)*np.sin(n*y)*np.sin(n*z) + n*np.cos(n*x)*np.cos(n*y)*np.cos(n*z)

    return np.array([a, b, c])
