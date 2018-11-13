#! /usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from past.utils import old_div
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
import Radial_int_pixels


def rgb2gray(pixels):
    # pixels must be in rgba format
    # brightness  =  sqrt( .241 R2 + .691 G2 + .068 B2 )
    print(pixels)
    I = [x**2 for x in pixels[:3]]
    w = [.241, .691, .068]  # http://www.nbdtech.com/Blog/archive/2008/04/27/Calculating-the-Perceived-Brightness-of-a-Color.aspx

    return np.sqrt(np.dot(I, w))


def mean(pixels):

    means = []
    maxes = []
    for i in range(pixels.shape[0]):
        means.append(np.mean(pixels[i, :]))
        maxes.append(max(pixels[i, :]))

    return np.mean(means), max(maxes)


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


def dist_mat(pixels, dim):

    dim = float(dim)
    pixels = int(pixels)

    A = np.zeros([2, pixels, pixels])

    m = np.shape(A)[1]
    n = np.shape(A)[2]
    shift = m - 1

    for i in range(m):
        for j in range(n):
            A[0, i, j] = (2 * j - shift) * dim / 2
            A[1, i, j] = (- 2 * i + shift) * dim / 2

    angle_matrix = np.zeros([pixels, pixels])
    for i in range(m):
        for j in range(n):
            x, y = A[:, i, j]
            angle_matrix[i, j] = np.arctan(old_div(y, x)) * (old_div(180, np.pi))

    distance_matrix = np.zeros([pixels, pixels])
    for i in range(m):
        for j in range(n):
            x = A[:, i, j]
            dist = np.sqrt(np.dot(x, x))
            distance_matrix[i, j] = dist

    maxes = []
    for i in range(m):
        maxes.append(max(distance_matrix[:, i]))

    max_dist = max(maxes)

    return distance_matrix, max_dist, angle_matrix


def radial_int(pixels, dr, dim, radius_bounds='all', angle_bounds='all'):

    dr = float(dr)
    npixels = pixels.shape[0]
    dim = float(dim)

    if radius_bounds == 'all':
        upper_rbound = npixels * np.sqrt(2) / 2
        lower_rbound = 0
    else:
        upper_rbound = radius_bounds[1]
        lower_rbound = radius_bounds[0]

    if angle_bounds == 'all':
        upper_abound = 180
        lower_abound = 0
    else:
        upper_abound = angle_bounds[1]
        lower_abound = angle_bounds[0]

    distance_mat, max_dist, angle_matrix = dist_mat(npixels, dim)

    intensities = np.zeros([int(old_div(max_dist,dr))])

    bins = len(intensities)

    x = np.linspace(0, max_dist, bins)

    pixels_parsed = np.zeros([npixels, npixels])

    pcount = 0
    for i in range(npixels):
        for j in range(npixels):
            intensity = pixels[i, j]
            distance = distance_mat[j, i]  # indices flipped in this matrix ...
            angle = angle_matrix[i, j]
            if lower_rbound <= distance <= upper_rbound:
                if lower_abound <= abs(angle) <= upper_abound:
                    pixels_parsed[i, j] = intensity
                    pcount += 1
                    bin_no = int(np.floor(distance/max_dist * bins))
                    if bin_no == bins:
                        intensities[bin_no - 1] += intensity
                    else:
                        intensities[bin_no] += intensity

    area_sum = sum(intensities)
    average_pixel_intensity = old_div(area_sum, pcount)

    return intensities, area_sum, average_pixel_intensity, pixels_parsed

#im = Image.open("soft confine WAXS raw.tif")  # open the image
im = Image.open("WAXS.png")  # open the image

box = (510, 90, 1000, 580)  # crop out the 2D soft confinement WAXS - this needs to be square
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
cmap = 'viridis'
interp = 'gaussian'
plt.figure()
im1 = plt.imshow(pixels, cmap=cmap, interpolation=interp) #, vmin=Imin, vmax=Imax)
plt.show()
exit()
rb = [130, 170]  # radius to integrate between. Length is in pixels. Length of 1 pixel = 1
ab = [20, 65]
rb2 = [0, 200]

integration, area_sum, average_pixel_intensity, pixels_parsed = radial_int(pixels, 1, 1, radius_bounds=rb)
integration2, area_sum2, average_pixel_intensity2, pixels_parsed2 = radial_int(pixels, 1, 1, radius_bounds=rb2, angle_bounds=ab)

print('Average intensity in whole ring: %s' % average_pixel_intensity)
print('Average intensity in spots: %s' % average_pixel_intensity2)

Imin, Imax = bounds(pixels)
cmap = 'viridis'
interp='gaussian'
plt.figure()
im1 = plt.imshow(pixels_parsed, cmap=cmap, interpolation=interp, vmin=Imin, vmax=Imax)
plt.figure()
im2 = plt.imshow(pixels_parsed2, cmap=cmap, interpolation=interp, vmin=Imin, vmax=Imax)
plt.figure()
im3 = plt.imshow(pixels, cmap=cmap, interpolation=interp, vmin=Imin, vmax = Imax)

x_axis = np.linspace(0, len(integration), len(integration))

plt.figure()
plt.plot(x_axis, integration)

plt.show()
exit()

Imin, Imax = bounds(pixels)
Imin2, Imax2 = bounds(pixels2)


# plt.figure(1)
# im = plt.imshow(pixels, cmap='Dark2', interpolation='none', vmin=Imin, vmax=Imax)

plt.figure(2)
im = plt.imshow(pixels2, cmap='Dark2', interpolation='none', vmin=Imin2, vmax=Imax2)


# dark = pixels2[160:180, 120:140]
# dark_avg, dark_max = mean(dark)
#
# plt.figure(2)
# im = plt.imshow(dark, cmap='Dark2', interpolation='none', vmin=Imin, vmax=Imax)
#
# light = pixels2[230:270, 90:110]
# light_avg, light_max = mean(light)
# plt.figure(3)
# im = plt.imshow(light, cmap='Dark2', interpolation='none', vmin=Imin, vmax=Imax)


# print dark_avg
# print light_avg
#
# print dark_avg/light_avg
# plt.figure()
# pi_stacking = pixels2[78:90, 200:300]
# im = plt.imshow(pi_stacking, cmap='Dark2', interpolation='none', vmin=Imin, vmax=Imax)
#
# middle_stack = pixels2[177:183, 227:261]

# print mean(middle_stack)
# print mean(dark)
# print mean(light)
# print mean(pi_stacking)
plt.show()

