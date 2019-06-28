#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def check_ellipsoid(pt, abc):
	""" Return whether a point is contained in an ellipsoid"""

	return np.sum(np.square(pt) / np.square(abc)) < 1

abc = np.array([1, 2, 3])

npts = 10000
coords = np.zeros([npts, 3])

count = 0
while count < npts:

	pt = np.array([np.random.uniform(-abc[0], abc[0]), np.random.uniform(-abc[1], abc[1]), np.random.uniform(-abc[2], abc[2])])
	if check_ellipsoid(pt, abc):
		coords[count, :] = pt
		count += 1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2])

print('Emperical Results:')
print("Standard deviation in x: %.2f" % np.std(coords[:, 0]))
print("Standard deviation in y: %.2f" % np.std(coords[:, 1]))
print("Standard deviation in z: %.2f" % np.std(coords[:, 2]))

# Un-normalized principal component analysis
cov = np.cov(coords.T)
eig_vals, eig_vecs = np.linalg.eig(cov)
eig_vals = eig_vals[np.argsort(eig_vals)][::-1]

print('PCA results')
print('Eigenvectors:\n %s' % eig_vecs)
print("Standard deviation of highest variance component: %.2f" % eig_vals[0]**0.5)
print("Standard deviation of second highest variance component: %.2f" % eig_vals[1]**0.5)
print("Standard deviation of lowest variance component: %.2f" % eig_vals[2]**0.5)

plt.show()
