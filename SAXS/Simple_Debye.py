import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# A script to try applying the Debye equation for a simple lattice
# I'll choose a simple cubic lattice that is already in fractional coordinates
# Body Centered Cubic:

x_vals = [0, 1, 1, 0, 0, 1, 1, 0, 0.5, 1, 2, 2, 1, 1, 2, 2, 1, 1.5]
y_vals = [0, 0, 0, 0, 1, 1, 1, 1, 0.5, 0, 0, 0, 0, 1, 1, 1, 1, 0.5]
z_vals = [0, 0, 1, 1, 0, 0, 1, 1, 0.5, 0, 0, 0, 0, 1, 1, 1, 1, 0.5]

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# x_mesh, y_mesh, z_mesh = np.meshgrid(x_vals, y_vals, z_vals)
# ax.scatter(x_mesh, y_mesh, z_mesh)
# plt.show()

rij = np.zeros((len(x_vals), len(x_vals)))
for i in range(0, len(x_vals)):  # loops through each atom
    for j in range(0, len(x_vals)):  # compares each atom to all other atoms including itself
        rij[i, j] = math.sqrt((x_vals[i] - x_vals[j])**2 + (y_vals[i] - y_vals[j])**2 + (z_vals[i] - z_vals[j])**2)

plt.hist(rij)
plt.show()

# fitted parameters for Sodium (International Crystallographic Tables, Volume C, Table 6.1.1.3)
a_p = [3.25650, 3.93620, 1.39980, 1.00320]  # [a1, a2, a3, a4]  _p stands for 'parameters'
b_p = [2.66710, 6.11530, 0.200100, 14.0390]  # [b1, b2, b3, b4]
c_p = 0.404000


def scattering_factor(a, b, c, q):
    # d: 1/2d where d is the d-spacing
    sol = q/(4*math.pi)
    f = 0  # intialize summation
    for i in range(0, 4):  # for a1, b1, a2, b2, a3, b3, a4, b4
        f += a[i]*math.exp(-b[i]*sol**2)  # this function was used to fit the data in the Crystallographic tables
    f += c  # the constant c is contained outside the summation
    return f

q = np.linspace(0, 6.5, 1000)
I = np.zeros((len(q), 1))
for k in range(0, len(q)):
    for i in range(0, len(x_vals)):
        for j in range(0, len(y_vals)):
            f = scattering_factor(a_p, b_p, c_p, q[k])  # only one f since there is only one type of atom
            if rij[i, j] != 0:
                I[k] += f*f*(math.sin(q[k]*rij[i, j])/q[k])

plt.plot(q, I)
plt.show()