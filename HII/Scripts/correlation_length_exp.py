#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def lorentz(points, a, b):
    """
    :param p: lorentzian parameters : [full width half max (FWHM), position of maximum]
    :param p: position
    :return:
    """

    w = np.pi / a

    x = (b - points) / (w/2)

    return 1 / (1 + x**2)


def errorfunc(p, points, z):
    return lorentz(p, points) - z


waxs = np.load('WAXS.npy')

center = np.array([530, 473])
hw = 400  # height and width (pixels)
r_high = hw - 155
r_low = hw - 215
bins = 720
db = 180 / float(bins)

waxs = waxs[center[0]-hw:center[0]+hw, center[1]-hw:center[1]+hw]

qpi = 1.7  # q value of pi-stacking reflection
axial = np.copy(waxs[:, int(waxs.shape[0]/2)])  # hold x constant at the center, and get all y values. The array is of the shape [y, x]
mid = axial.shape[0] / 2
axial[int(mid - 150):int(mid + 150)] = 0  # zero out middle values so we get the right maximum
Imax = np.amax(axial)  # max value of Intensity. It will correspond to location of pi-stacking reflection
Imax_pixel = np.where(axial == Imax)[0][0]
pixel_to_q = 1.7 / abs(mid - Imax_pixel)
qmax = mid*pixel_to_q
r_high_pix = 155 * pixel_to_q
r_low_pix = 215 * pixel_to_q

start = 640
end = 800
middle = 401

noise = np.mean(waxs[775:800, middle])  # background noise

x = np.linspace(-400, 400, 800)[start:end]*pixel_to_q
y = (waxs[start:end, middle] - noise)/max(waxs[start:end, middle] - noise)
p = np.array([30, 1.7])

solp, cov_x = curve_fit(lorentz, x, y, p)

print('FWHM: %1.2f' % (np.pi/solp[0]))
print('Position of Maximum: %1.2f' % solp[1])
print('Correlation length: %1.2f +/- %1.2f' % (solp[0], np.sqrt(cov_x[0, 0])))

plt.plot(x, y, label='Raw data')
plt.plot(x, lorentz(x, solp[0], solp[1]), '--', color='black', label='Lorentzian fit')
plt.xlabel('q ($\AA^{-1}$)', fontsize=14)
plt.ylabel('Normalized intensity', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend()
# plt.title('FWHM: %1.3f $\AA^{-1}$' % solp[0])
plt.tight_layout()
plt.savefig('Correlation_length_exp.png')
plt.show()


# plt.figure()
# plt.imshow(waxs, cmap='jet', vmax=0.05, extent=[-qmax, qmax, -qmax, qmax])
# plt.show()