#!/usr/bin/env python

from PIL import Image
import matplotlib.pyplot as plt
import numpy as np

im = Image.open("2D-SAXS.png")  # open the image

box = (0, 0, 1024, 1024)  # crop out the 2D soft confinement WAXS - this needs to be square
#box = (0, 0, 1042, 1042)  # crop out the 2D soft confinement WAXS - this needs to be square

region = im.crop(box).convert('L')  # convert to grayscale
pix = region.load()  # get the pixel values

region2 = im.crop(box)  # keep it in RGBA format
pix2 = region2.load()  # get color pixel values

pixels = np.zeros(region.size)  # convert pix to a numpy array
pixels = np.zeros([pixels.shape[1], pixels.shape[0]])
pixels2 = np.zeros(pixels.shape)

for i in range(pixels.shape[0]):
    for j in range(pixels.shape[1]):
        pixels[i, j] = pix[j, i]
        # pixels2[i, j] = rgb2gray(pix2[j, i])

# Imin, Imax = bounds(pixels)
cmap = 'jet'
interp = 'gaussian'
plt.figure()
im1 = plt.imshow(pixels, cmap=cmap, interpolation=interp) #, vmin=Imin, vmax=Imax)
plt.show()
exit()
im = Image.open("2D-SAXS.png")
box=(510, 90, 1000, 580)
region = im.crop(box)
pixels = region.load()
print(pixels)
exit()
plt.imshow(pixels, cmap='viridis', interpolation='gaussian')
plt.show()

