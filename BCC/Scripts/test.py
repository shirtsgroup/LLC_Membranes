#!/usr/bin/env python

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import bcc_class
from llclib import transform
from llclib import file_rw

LC = bcc_class.LC('Dibrpyr14.gro')

n = np.array([-0.5, -0.6, 0.3])
# n = LC.linevector

pts = np.zeros([10, 3])
for k in range(10):
    pts[k, :] = (k/10)*n

normal = transform.translate(pts, pts[0, :], np.array([0, 0, 0]))

# xyz = np.concatenate((LC.xyz, normal))

R = transform.Rvect2vect(LC.linevector, n)  # rotation matrix to rotate monomer in same direction as n

# translate to origin
xyz_origin = transform.translate(LC.xyz, LC.reference, np.array([0, 0, 0]))

xyz_origin = transform.rotate_coords(xyz_origin, R)  # rotate all points in bcc monomer with rotation matrix

v1 = np.array([np.mean(xyz_origin[23:27, 0]), np.mean(xyz_origin[23:27, 1]), np.mean(xyz_origin[23:27, 2])])
v2 = xyz_origin[37, :]
print((v1 - v2) / np.linalg.norm(v1 - v2))
print(n / np.linalg.norm(n))
# file_rw.write_gro_pos(xyz_origin, 'test.gro')
file_rw.write_gro_pos(np.concatenate((xyz_origin, normal)), 'test.gro')
exit()

# find avg location of reference atoms after rotation (since it will change)
ref = np.zeros([3])
for j in range(len(LC.ref_index)):
    ref += xyz[LC.ref_index[j], :]
ref /= len(LC.ref_index)

xyz = transform.translate(xyz, ref, grid[i, :])  # move monomer to grid point w.r.t. reference point on monomer
pt1 = xyz[37, :]
pt2 = np.array([np.mean(xyz[23:27, 0]), np.mean(xyz[23:27, 1]), np.mean(xyz[23:27, 2])])
print(n / np.linalg.norm(n), pt1 - pt2 / (np.linalg.norm(pt1 - pt2)))
xyz = np.concatenate((xyz, normal))
if count == 0:
    x = np.concatenate((grid, xyz))
else:
    x = np.concatenate((x, xyz))

count += 1

