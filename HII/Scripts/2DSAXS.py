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


im = Image.open("2D-SAXS.png")

# For PIL, (0, 0) is in the upper left hand corner
leftshift = 0
uppershift = 0
rightshift = 108
lowershift = 10
box = (leftshift, uppershift, im.size[0] - rightshift, im.size[1] - lowershift)
region = im.crop(box)#.convert('L')

pix = region.load()
pixels = np.zeros(region.size)
pixels = pixels.T

xsize = region.size[0]
ysize = region.size[1]

for i in range(pixels.shape[0]):
	for j in range(pixels.shape[1]):
		pixels[i, j] = pix[j, i]

Imin, Imax = bounds(pixels)

cmap = 'jet'
interp = 'gaussian'

xsection = pixels[ysize // 2, :]

left_peak = np.argmax(xsection[: xsize //2])
right_peak = np.argmax(xsection[xsize //2:]) + xsize // 2

q = np.linspace(0, pixels.shape[0], pixels.shape[0])

qbin = .179 / (xsize // 2 - left_peak)
qmax = qbin *float(pixels.shape[0] // 2)

qx = np.linspace(-qmax, qmax, pixels.shape[0])
qy = np.linspace(-qmax, qmax, pixels.shape[1])

plt.imshow(pixels, extent=[-qmax, qmax, -qmax, qmax], cmap='jet')
plt.colorbar()
plt.xlabel('$q_r\ (\AA^{-1})$', fontsize=16)
plt.ylabel('$q_z\ (\AA^{-1})$', fontsize=16)
plt.tight_layout()
#plt.savefig('2DSAXS_axes.png')

third_peak_left = 0
while qy[third_peak_left] < -0.38:
	third_peak_left += 1

third_peak_right = 0
while qy[third_peak_right] < 0.38:
	third_peak_right += 1

third_peak_left_I = xsection[np.argmax(xsection[:third_peak_left])]
third_peak_right_I = xsection[np.argmax(xsection[third_peak_right:]) + third_peak_right]
#print('Third peak intensity (left) = %.2f' % third_peak_left_I)
#print('Third peak intensity (right) = %.2f' % third_peak_right_I)

avg_third_peak = (third_peak_left_I + third_peak_right_I) / 2
#print('Average Intensity of 3rd peak : %.2f' % avg_third_peak)

waxs_third_peak = 2.7033  # taken from 2D WAXS pattern

norm = avg_third_peak / waxs_third_peak

print('Average R-pores Intensity = %.2f' % (xsection[np.argmax(xsection)]/norm))

plt.figure()

plt.plot(qy, xsection / norm, linewidth=2)
plt.ylabel('Intensity', fontsize=14)
plt.xlabel('$q_r (\AA^{-1})$', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('saxs_xsection.pdf')
plt.show()
