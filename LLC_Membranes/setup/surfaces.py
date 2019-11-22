#!/usr/bin/env python

import numpy as np


class SurfaceError(Exception):
    """ Raised if invalid phase specified """

    def __init__(self, message):

        super().__init__(message)


def SchwarzD(x, period):
    """
    :param x: a vector of coordinates (x1, x2, x3)
    :param period: length of one period

    :return: An approximation of the Schwarz D "Diamond" infinite periodic minimal surface
    """

    n = 2*np.pi / period  # might be just pi / period

    a = np.sin(n*x[0])*np.sin(n*x[1])*np.sin(n*x[2])
    b = np.sin(n*x[0])*np.cos(n*x[1])*np.cos(n*x[2])
    c = np.cos(n*x[0])*np.sin(n*x[1])*np.cos(n*x[2])
    d = np.cos(n*x[0])*np.cos(n*x[1])*np.sin(n*x[2])

    return a + b + c + d


def gyroid(x, period):

    n = 2 * np.pi / period

    a = np.sin(n*x[0])*np.cos(n*x[1])
    b = np.sin(n*x[1])*np.cos(n*x[2])
    c = np.sin(n*x[2])*np.cos(n*x[0])

    return a + b + c


def Sphere(x):

    return np.square(x).sum()


def gridgen(surf, low, high, n, c=0, tol=0.01):
    """ Generate an n x n x n grid and reduce it to points that lie close to an implicit surface defined by `surf`

    :param surf: name of implicit surface to approximate
    :param low: lowest coordinate on each axis of grid
    :param high: highest coordinate on each axis of grid (if low is 0, then high is the length of the box vector)
    :param n: number of grid points in each dimension
    :param c: value of c in F(x, y, z) = c. c = 0 corresponds to zero mean curvature. c < 0 is negative mean curvature
    and c > 0 is positive mean curvature
    :param tol: surface--point distance tolerance. Anything within this cutoff can be used to approximate the surface.

    :type surf: str
    :type low: float
    :type high: float
    :type n: int
    :type c: float
    :type tol: float
    """

    # make a cubic grid
    bin_size = (high - low) / n
    x = np.linspace(low, high - bin_size, n)
    y = np.linspace(low, high - bin_size, n)
    z = np.linspace(low, high - bin_size, n)

    if surf.lower() in ['ia3d', 'gyroid']:

        gyro = gyroid([x[:, None, None], y[None, :, None], z[None, None, :]], high - low)

        gyro_eval = np.zeros([n**3, 3])

        count_gyro = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if abs(gyro[i, j, k] - c) < tol:
                        gyro_eval[count_gyro, :] = [x[i], y[j], z[k]]
                        count_gyro += 1

        grid = gyro_eval[:count_gyro, :]

    elif surf.lower() in ['pn3m', 'schwarzd', 'diamond']:

        schwarz = SchwarzD([x[:, None, None], y[None, :, None], z[None, None, :]], high - low)
        schwarz_eval = np.zeros([n**3, 3])
        count_schwarz = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if -tol < schwarz[i, j, k] - c < tol:
                        schwarz_eval[count_schwarz, :] = [x[i], y[j], z[k]]
                        count_schwarz += 1

        grid = schwarz_eval[:count_schwarz, :]

    elif surf.lower() == 'sphere':

        sphere = Sphere([x[:, None, None], y[None, :, None], z[None, None, :]])#, high - low)
        sphere_eval = np.zeros([n**3, 3])
        count_sphere = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if -tol < sphere[i, j, k] - c**2 < tol:
                        sphere_eval[count_sphere, :] = [x[i], y[j], z[k]]
                        count_sphere += 1

        grid = sphere_eval[:count_sphere, :]

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

    if surf.lower() in ['ia3d', 'gyroid']:

        a = n*np.cos(n*x)*np.cos(n*y) - n*np.sin(n*x)*np.sin(n*z)
        b = -n*np.sin(n*y)*np.sin(n*x) + n*np.cos(n*y)*np.cos(n*z)
        c = -n*np.sin(n*y)*np.sin(n*z) + n*np.cos(n*z)*np.cos(n*x)

    elif surf.lower() in ['pn3m', 'schwarzd', 'diamond']:

        a = n*np.cos(n*x)*np.sin(n*y)*np.sin(n*z) + n*np.cos(n*x)*np.cos(n*y)*np.cos(n*z) - n*np.sin(n*x)*np.sin(n*y)*np.cos(n*z) - n*np.sin(n*x)*np.cos(n*y)*np.sin(n*z)
        b = n*np.sin(n*x)*np.cos(n*y)*np.sin(n*z) - n*np.sin(n*x)*np.sin(n*y)*np.cos(n*z) + n*np.cos(n*x)*np.cos(n*y)*np.cos(n*z) - n*np.cos(n*x)*np.sin(n*y)*np.sin(n*z)
        c = n*np.sin(n*x)*np.sin(n*y)*np.cos(n*z) - n*np.sin(n*x)*np.cos(n*y)*np.sin(n*z) - n*np.cos(n*x)*np.sin(n*y)*np.sin(n*z) + n*np.cos(n*x)*np.cos(n*y)*np.cos(n*z)

    elif surf.lower() == 'sphere':

        a = 2 * x
        b = 2 * y
        c = 2 * z

    else:

        raise SurfaceError('The surface %s is named incorrectly or is not implemented' % surf)

    n = np.array([a, b, c])

    return n / np.linalg.norm(n)
