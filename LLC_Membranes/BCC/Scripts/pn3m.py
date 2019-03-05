#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def SchwarzD(x):
    """
    :param x: a vector of coordinates (x1, x2, x3)
    :return: An approximation of the Schwarz D "Diamond" infinite periodic minimal surface
    """

    a = np.sin(x[0])*np.sin(x[1])*np.sin(x[2])
    b = np.sin(x[0])*np.cos(x[1])*np.cos(x[2])
    c = np.cos(x[0])*np.sin(x[1])*np.cos(x[2])
    d = np.cos(x[0])*np.cos(x[1])*np.sin(x[2])

    return a + b + c + d


def gyroid(x):

    a = np.sin(x[0])*np.cos(x[1])
    b = np.sin(x[1])*np.cos(x[2])
    c = np.sin(x[2])*np.cos(x[0])

    return a + b + c


def gyroidz(x):

    a = np.sin(x[0])*np.cos(x[1])

low = 0
high = 2*np.pi
n = 50

x = np.linspace(low, high, n)
y = np.linspace(low, high, n)
z = np.linspace(low, high, n)

schwarz = SchwarzD([x[:, None, None], y[None, :, None], z[None, None, :]])
gyro = gyroid([x[:, None, None], y[None, :, None], z[None, None, :]])

schwarz_eval = np.zeros([n**3, 3])
gyro_eval = np.zeros([n**3, 3])

count_schwarz = 0
for i in range(n):
    for j in range(n):
        for k in range(n):
            if -0.05 < schwarz[i, j, k] < 0.05:
                schwarz_eval[count_schwarz, :] = [x[i], y[j], z[k]]
                count_schwarz += 1

count_gyro = 0
for i in range(n):
    for j in range(n):
        for k in range(n):
            if -0.05 < gyro[i, j, k] < 0.05:
                gyro_eval[count_gyro, :] = [x[i], y[j], z[k]]
                count_gyro += 1

pn3m = schwarz_eval[:count_schwarz, :]
Ia3d = gyro_eval[:count_gyro, :]


fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

# X, Y = np.meshgrid(pn3m[:, 0], pn3m[:, 1])

# ax.plot_surface(X, Y, Z)

ax.scatter(pn3m[:3, 0], pn3m[:3, 1], pn3m[:3, 2], marker='.')

# ax.plot_trisurf(pn3m[:, 0], pn3m[:, 1], pn3m[:, 2])

# fig = plt.figure(2)
# ax2 = fig.add_subplot(111, projection='3d')
# ax2.scatter(Ia3d[:, 0], Ia3d[:, 1], Ia3d[:, 2], marker='.')

plt.show()