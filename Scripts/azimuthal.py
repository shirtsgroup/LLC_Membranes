#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

from builtins import range
from past.utils import old_div
import numpy as np
import argparse
import matplotlib.pyplot as plt

r_low = 1.2
r_high = 1.6
bins = 90
bounds = 45  # assumed to be symmetric. So bounds will be between +/- bounds

# def radial_integrate(D, Nbins, outputname):
#
#     SF = D[:, :, :, 3]
#
#     R = (D[:, :, :, 0]**2).astype(np.float16) + (D[:, :, :, 1]**2).astype(np.float16) + (D[:, :, :, 2]**2).astype(np.float16)
#     H, E = np.histogram(R, bins=Nbins, weights=SF)
#     Hc, E = np.histogram(R, bins=Nbins)
#     Hc = np.where(Hc != 0, Hc, 1.0)
#     H /= Hc
#     H[:1] = 0.0
#     H /= np.amax(H)
#     plt.plot(E[:-1], H)
#     plt.ylim(0, 0.5)
#     plt.xlim(0.1, 0.5)
#     plt.savefig(outputname, dpi=DPI)

sf = np.load('out_wiggle_sf.npz')

D = sf['kgridplt']

# R = (D[:, :, :, 0]**2).astype(np.float16) + (D[:, :, :, 1]**2).astype(np.float16) + (D[:, :, :, 2]**2).astype(np.float16)

x = D.shape[0]
y = D.shape[1]
z = D.shape[2]

angles = np.zeros(bins)
norm = np.zeros_like(angles)

normal = np.array([0, 0, 1])

for i in range(x):
    for j in range(y):
        for k in range(z):
            v = D[i, j, k, :3]
            r = np.linalg.norm(v)
            if r_high > r > r_low:
                # find angle between vector and xy plane
                vn = np.dot(v, normal)
                nn = np.linalg.norm(normal)
                vv = np.linalg.norm(v)
                angle = np.arcsin(old_div(vn, (nn * vv))) * (old_div(180, np.pi))
                bin = int((bins/2) + ((angle/90)*(bins/2)))
                if bin == bins:
                    bin -= 1
                angles[bin] += D[i, j, k, 3]
                norm[bin] += 1

avg = angles / norm
angles = np.linspace(-90, 90, bins + 1)
bin_angles = [(angles[i] + angles[i + 1])/2 for i in range(bins)]
width = angles[1] - angles[0]

bound1 = int((bins/2)-bounds)
bound2 = int((bins/2) + bounds)

print(np.mean(avg[int((bins/2)-bounds):int((bins/2) + bounds)]))
print(np.max(avg[int((bins/2)-bounds):int((bins/2) + bounds)]))
print(np.max(avg[int((bins/2)-bounds):int((bins/2) + bounds)])/np.mean(avg[int((bins/2)-bounds):int((bins/2) + bounds)]))

plt.bar(bin_angles[bound1:bound2], avg[bound1:bound2], width=width)
# plt.hist(np.array(angles))
plt.show()


# radial_integrate(grid, 750, dir+"radial.png")