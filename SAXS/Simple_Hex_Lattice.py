# A simple hexagonal lattice with one point at each pore location

import numpy as np
import math
import matplotlib.pyplot as plt

a = 4.0
rows = 2
columns = 2
angle = 2*math.pi/3

# Make a parallelogram of points
x = np.zeros((rows, columns))
y = np.zeros((rows, columns))
for i in range(0, rows):
    for j in range(0, columns):
        y[i, j] = -i*(a/(rows - 1))*math.sin(angle)
        x[i, j] = j*(a/(rows - 1)) - i*(a/(rows - 1))*math.cos(angle)


# Create matrix of x, y pairs
# generate a bunch of points in the corners:
# x_range = np.linspace(0, 0.5, 11)
# x_uleft = x_range
# x_lleft = x_range +
new_x = np.array([.5, 3.5, 2.5, 5.5, .25, 4.25, 1.75, 5.75, .5, 4, 2, 5.5])
new_y = np.array([0, 0, -4*math.sin(math.pi/3), -4*math.sin(math.pi/3), -.5*math.sin(math.pi/3), -.5*math.sin(math.pi/3), -3.5*math.sin(math.pi/3), -3.5*math.sin(math.pi/3), -.5*math.sin(math.pi/3), -.5*math.sin(math.pi/3), -3.5*math.sin(math.pi/3), -3.5*math.sin(math.pi/3)])
x_ = np.concatenate((x[0, :], x[1, :], new_x), axis=0)
y_ = np.concatenate((y[0, :], y[1, :], new_y), axis=0)
xy = np.zeros((2, len(x_)))
xy[0, :] = x_
xy[1, :] = y_

# plt.scatter(x_, y_)
# plt.show()

# transformation matrix to convert to fractional coordinates
transform = np.matrix(([1/float(a), -1/(float(a)*math.tan(angle))], [0, 1/(float(a)*math.sin(angle))]))

frac = np.dot(transform, xy)
plt.scatter(frac[0, :], frac[1, :])
plt.show()

# Calculate distance of ion from all other ions
dist = []
dist_xindex = []
for i in range(0, len(x_)):
    dist.append([])
    dist_xindex.append(i)
    x_ref = xy[0, i]
    y_ref = xy[1, i]
    for j in range(0, len(x_)):
        dist[i].append(math.sqrt((x_ref - xy[0, j])**2 + (y_ref - xy[1, j])**2))


def Ihk(x_, bins_index):
    # Real part
    real = 0
    for i in range(0, len(x_)):
        real += math.cos(2*math.pi*frac[0, bins_index[i]])

    # Imaginary Part
    imaginary = 0
    for i in range(0, len(x_)):
        imaginary += math.sin(2*math.pi*frac[1, bins_index[i]])

    I = real**2 + imaginary**2
    return I

interval = 0.001  # Distances in angstroms to group together
groups = np.linspace(0, 2*a, 2*a/interval + 1)

bins = []
bins_index = []
for k in range(0, len(groups)):  # looks at each interval which defines the bins
    bins.append([])
    bins_index.append([])
    for j in range(0, len(dist)):  # look at each ion in dist
        for i in range(0, len(x_)):  # look at distance of all other ions from that ion
            if groups[k] <= dist[j][i] < groups[k + 1]:  # if it meets the distance criteria for this bin, add it
                bins[k].append(dist[j][i])
                bins_index[k].append(j)  # append the associated index of this ion so that it has an identity

# print bins
# print bins_index

q_peaks = []  # location where q peaks should be located
Intensity = []
for i in range(0, len(bins)):
    if bins[i]:  # if the list is empty
        a = bins.index(bins[i])
        d = (groups[a] + groups[a + 1])/2
        q = 2*math.pi/d
        if q <= 5:
            q_peaks.append(q)
            Intensity.append(Ihk(bins[i], bins_index[i]))

plt.plot(q_peaks, Intensity)
plt.title('SAXS')
plt.ylabel('Intensity')
plt.xlabel('q')
plt.show()