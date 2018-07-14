#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def inverse(x, a):
        return a / x

perfect_z_noise = np.array([0, 0.05, 0.1, .15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])
perfect_z_noise_intensity = np.array([200, 181.06, 134.65, 83.18, 41.51, 7.1, 1.03, 1.02, 1, 1, 1, 1, 1])
perfect_z_noise_FWHM = np.array([.382, .382, .382, .382, .382, .382])  # doesn't change 

plt.figure()
plt.plot(perfect_z_noise, perfect_z_noise_intensity)
plt.ylabel('R-$\pi$ Intensity')
plt.xlabel('Gaussian noise ($\sigma$)')
plt.tight_layout()
plt.savefig('z_noise_perfect.png')

perfect_xy_noise = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
perfect_xy_noise_intensity = np.array([200]) # doesn't change
perfect_xy_noise_FWHM = np.array([0.370, 0.367, 0.360, 0.348, .331, .318, 0.301, .284])

plt.figure()
plt.plot(perfect_xy_noise, perfect_xy_noise_FWHM)
plt.ylabel('R-$\pi$ Intensity')
plt.xlabel('FWHM ($\AA^{-1}$)')
plt.tight_layout()
plt.savefig('xy_noise_perfect.png')
plt.show()
exit()

noise = np.array([0, 0.2, 0.4, 0.6, 0.8, 1])
#random_columns = []

random_columns_z_noise = np.array([0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6, 0.8, 1])
random_columns_z_noise_intensity = np.array([10.08, 9.18, 8.56, 7.67, 6.67, 4.43, 2.96, 1.69, 1.25, 1.07, 1.04, 0.99, 1.01])
random_columns_z_noise_FWHM = [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]

plt.figure()
plt.plot(random_columns_z_noise, random_columns_z_noise_intensity, linewidth=2)
plt.ylabel('R-$\pi$ Intensity')
plt.xlabel('Gaussian noise ($\sigma$)')
plt.savefig('random_columns_z_noise.png')

random_columns_xy_noise = np.array([0, 0.1, 0.2, 0.4, 0.6, 0.8, 1])
random_columns_xy_noise_intensity = np.array([10.08, 9.71, 10.41, 10.50, 10.00, 10.20])
random_columns_xy_noise_FWHM = [np.inf, 2.640, 1.419, 0.715, 0.475, 0.365, 0.285]

x_fine = np.linspace(random_columns_xy_noise[0], random_columns_xy_noise[-1], 1000)

p = [1]

solp_inverse, cov_x = curve_fit(inverse, random_columns_xy_noise[1:], random_columns_xy_noise_FWHM[1:], p)

plt.figure()
plt.plot(random_columns_xy_noise, random_columns_xy_noise_FWHM)
start = 75
plt.plot(x_fine[start:], inverse(x_fine[start:], solp_inverse[0]), '--', color='green', linewidth=2, label='Fit to $a/x$')
plt.xlabel('Gaussian noise ($\sigma$)')
plt.ylabel('FWHM')
plt.legend()
plt.tight_layout()
plt.savefig('random_columns_xy_noise.png')

random_layers_xy_noise = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0])
random_layers_xy_noise_intensity = np.array([200, 200, 200, 200, 200, 200])
random_layers_xy_noise_FWHM = np.array([0.372, 0.361, 0.335, 0.302, 0.267, 0.234])
err_FWHM = np.array([0.017, 0.048, 0.029, 0.013, 0.003, 0.000])

plt.figure()
plt.plot(random_layers_xy_noise, random_layers_xy_noise_FWHM)
plt.xlabel('Gaussian Noise ($\sigma$)')
plt.ylabel('FWHM')
plt.tight_layout()
plt.savefig('random_layers_xy_noise.png')

random_layers_z_noise = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0])
random_layers_z_noise_intensity = np.array([200, 130.68, 40.74, 6.47, 1.38, 1.01, 1.0, 1.0])
# FWHM becomes infinite (but not very intense) above 0.2
random_layers_z_noise_FWHM = np.array([0.372, 0.375, 0.382, 0.400])
err_FWHM = np.array([0.017, 0.057, 0.052])

plt.figure()
plt.plot(random_layers_z_noise, random_layers_z_noise_intensity)
plt.xlabel('Gaussian Noise ($\sigma$)')
plt.ylabel('Intensity')
plt.tight_layout()

plt.savefig('random_layers_z_noise.png')
plt.show()
