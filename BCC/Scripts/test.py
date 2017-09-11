from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

# x^2 + y^2 + z = 0 -- > z = -y**2 - x**2
# gradient vector : grad f = < 2x, 2y, 1 >


def sphere(x):

    return x[0] ** 2 + x[1] ** 2 + x[2] ** 2


def grad(v):

    g = np.array([2*v[0], 2*v[1], 2*v[2]])
    return g / np.linalg.norm(g)

n = 50
x = np.linspace(-1.5, 1.5, n)
y = np.linspace(-1.5, 1.5, n)
z = np.linspace(-1.5, 1.5, n)

pts = sphere([x[:, None, None], y[None, :, None], z[None, None, :]])

pts_eval = np.zeros([n**3, 3])
count = 0
for i in range(n):
    for j in range(n):
        for k in range(n):
            if .95 < pts[i, j, k] < 1.05:
                pts_eval[count, :] = [x[i], y[j], z[k]]
                count += 1

pts = pts_eval[:count, :]

gradv = np.zeros([pts.shape[0], 6])

for i in range(pts.shape[0]):
    gradv[i, :3] = pts[i, :]
    gradv[i, 3:] = grad(pts[i, :])

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

X, Y, Z, U, V, W = zip(*gradv)

ax.quiver(X, Y, Z, U, V, W, length=0.1)
# ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], marker='.')

plt.show()
