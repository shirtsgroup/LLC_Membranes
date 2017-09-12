from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np


def Pn3m(x, y, z):
    """
    :param x: a vector of coordinates (x1, x2, x3)
    :return: An approximation of the Schwarz D "Diamond" infinite periodic minimal surface
    """

    a = np.sin(x)*np.sin(y)*np.sin(z)
    b = np.sin(x)*np.cos(y)*np.cos(z)
    c = np.cos(x)*np.sin(y)*np.cos(z)
    d = np.cos(x)*np.cos(y)*np.sin(z)

    return a + b + c + d


def gyroid(x, y, z):

    a = np.sin(x)*np.cos(y)
    b = np.sin(y)*np.cos(z)
    c = np.sin(z)*np.cos(x)

    return a + b + c


def sphere(x, y, z):
    # sphere of radius 1
    return x ** 2 + y ** 2 + z ** 2 - 1


def gradient(v, surf):
    """
    :param v: vector of x, y, z coordinates
    :param phase: which implicit surface is being used to approximate the structure of this phase
    :return: The gradient vector (which is normal to the surface at x)
    """

    x = v[0]
    y = v[1]
    z = v[2]

    if surf == 'Ia3d' or surf == 'gyroid' or surf == 'ia3d':

        a = np.cos(x)*np.cos(y) - np.sin(x)*np.sin(z)
        b = -np.sin(y)*np.sin(x) + np.cos(y)*np.cos(z)
        c = -np.sin(y)*np.sin(z) + np.cos(z)*np.cos(x)

    elif surf == 'Pn3m' or surf == 'pn3m':

        a = np.cos(x)*np.sin(y)*np.sin(z) + np.cos(x)*np.cos(y)*np.cos(z) - np.sin(x)*np.sin(y)*np.cos(z) - np.sin(x)*np.cos(y)*np.sin(z)
        b = np.sin(x)*np.cos(y)*np.sin(z) - np.sin(x)*np.sin(y)*np.cos(z) + np.cos(x)*np.cos(y)*np.cos(z) - np.cos(x)*np.sin(y)*np.sin(z)
        c = np.sin(x)*np.sin(y)*np.cos(z) - np.sin(x)*np.cos(y)*np.sin(z) - np.cos(x)*np.sin(y)*np.sin(z) + np.cos(x)*np.cos(y)*np.cos(z)

    elif surf == 'sphere':

        a = 2*x
        b = 2*y
        c = 2*z

    return np.array([a, b, c])


def gridgen(surf, low, high, n):

    # make a cubic grid
    x = np.linspace(low, high, n)
    y = np.linspace(low, high, n)
    z = np.linspace(low, high, n)

    if surf == 'Ia3d' or surf == 'gyroid' or surf == 'ia3d':

        gyro = gyroid(x[:, None, None], y[None, :, None], z[None, None, :])

        gyro_eval = np.zeros([n**3, 3])

        count_gyro = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if -0.05 < gyro[i, j, k] < 0.05:
                        gyro_eval[count_gyro, :] = [x[i], y[j], z[k]]
                        count_gyro += 1

        grid = gyro_eval[:count_gyro, :]

    elif surf == 'Pn3m' or surf == 'pn3m':

        schwarz = Pn3m(x[:, None, None], y[None, :, None], z[None, None, :])
        schwarz_eval = np.zeros([n**3, 3])
        count_schwarz = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if -0.05 < schwarz[i, j, k] < 0.05:
                        schwarz_eval[count_schwarz, :] = [x[i], y[j], z[k]]
                        count_schwarz += 1

        grid = schwarz_eval[:count_schwarz, :]

    elif surf == 'sphere':

        pts = sphere(x[:, None, None], y[None, :, None], z[None, None, :])

        pts_eval = np.zeros([n**3, 3])
        count = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if -0.05 < pts[i, j, k] < 0.05:
                        pts_eval[count, :] = [x[i], y[j], z[k]]
                        count += 1

        grid = pts_eval[:count, :]

    else:
        print('The phase you selected is not defined (yet)')
        exit()

    return grid


def plot_implicit(fn, bbox=(-2.5, 2.5)):
    """
    :param fn: implicit function (plot where fn==0)
    :param bbox: the x,y,and z limits of plotted interval
    :return:
    """
    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    #MRS: trying to get aspect ratio equal
    #fig = plt.figure()
    fig = plt.figure(figsize=plt.figaspect(1.0))
    ax = fig.add_subplot(111, projection='3d')
    A = np.linspace(xmin, xmax, 100)  # resolution of the contour
    B = np.linspace(xmin, xmax, 50)  # number of slices
    A1, A2 = np.meshgrid(A, A)  # grid on which the contour is plotted

    for z in B:  # plot contours in the XY plane
        X, Y = A1, A2
        Z = fn(X, Y, z)
        cset = ax.contour(X, Y, Z+z, [z], zdir='z')
        # [z] defines the only level to plot for this contour for this value of z

    for y in B:  # plot contours in the XZ plane
        X, Z = A1, A2
        Y = fn(X, y, Z)
        cset = ax.contour(X, Y+y, Z, [y], zdir='y')

    for x in B:  # plot contours in the YZ plane
        Y, Z = A1, A2
        X = fn(x, Y, Z)
        cset = ax.contour(X+x, Y, Z, [x], zdir='x')

    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    ax.set_zlim3d(zmin, zmax)
    ax.set_xlim3d(xmin, xmax)
    ax.set_ylim3d(ymin, ymax)

    return ax
    # plt.show()

ax = plot_implicit(Pn3m)

surface = 'Pn3m'

grid = gridgen(surface, -2.5, 2.5, 40)
#
# ax.scatter(grid[:, 0], grid[:, 1], grid[:, 2])

gradv = np.zeros([grid.shape[0], 6])
for i in range(grid.shape[0]):
    gradv[i, :3] = grid[i, :]
    gradv[i, 3:] = gradient(grid[i, :], surface)

# fig = plt.figure(1)
# ax = fig.add_subplot(111, projection='3d')

X, Y, Z, U, V, W = zip(*gradv)
ax.quiver(X, Y, Z, U, V, W, length=0.25)
plt.show()
