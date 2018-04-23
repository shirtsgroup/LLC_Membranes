#!/usr/bin/env python

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

def bounds(pixels):
    # need a 2x2 array of pixel values
    mins = []
    for i in range(pixels.shape[1]):
        mins.append(min(pixels[:, i]))

    Imin = min(mins)

    maxes = []
    for i in range(pixels.shape[1]):
        maxes.append(max(pixels[:, i]))

    Imax = max(maxes)

    return Imin, Imax


im = Image.open("2D-SAXS_cropped.png")

box = (0, 0, 730, 720)
region = im.crop(box).convert('L')

pix = region.load()
pixels = np.zeros(region.size)
pixels = pixels.T

for i in range(pixels.shape[0]):
	for j in range(pixels.shape[1]):
		pixels[i, j] = pix[j, i]

Imin, Imax = bounds(pixels)
cmap = 'jet'
interp = 'gaussian'

q = np.linspace(0, pixels.shape[0], pixels.shape[0])

qbin = .179 / (pixels.shape[0] // 2 - 213)
qmax = qbin *float(pixels.shape[0] // 2)

qx = np.linspace(-qmax, qmax, pixels.shape[0])
qy = np.linspace(-qmax, qmax, pixels.shape[1])

plt.figure()

levels = np.linspace(100, np.amax(pixels), 200)
plt.contourf(qy, qx, pixels, levels=levels, cmap='jet', extend='both')
plt.colorbar()
plt.xlabel('$q_r (\AA^{-1})$')
plt.ylabel('$q_z (\AA^{-1})$')
plt.title('2D small angle X-ray scattering pattern')
#plt.plot([0, 700], [pixels.shape[0] //2, pixels.shape[0]//2])
#plt.plot([pixels.shape[1] //2, pixels.shape[1]//2], [0, 700])
plt.show()
