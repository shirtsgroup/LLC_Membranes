#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def inverse(x, a):
	return a / x

def squared(x, a):
	return a / x**2

x = np.array([1, 2, 3, 4, 5, 7, 10, 15])
x_fine = np.linspace(x[0], x[-1], 1000)
y = np.array([1.77, 0.999, .642, 0.465, .372, .264, .186, .123])

p = [1]

#solp_squared, cov_x = curve_fit(squared, x, y, p)
solp_inverse, cov_x = curve_fit(inverse, x, y, p)
plt.plot(x, y, color='black', linewidth=2, label='Raw Data')
#plt.plot(x_fine, squared(x_fine, solp_squared[0]), '--')
plt.plot(x_fine, inverse(x_fine, solp_inverse[0]), '--', linewidth=2, label='Fit to $a/x$')

plt.xlabel('Pore Radius ($\AA$)')
plt.ylabel('FWHM ($\AA^{-1}$)')
plt.legend()
plt.tight_layout()
plt.savefig('pore_radius_FWHM.png')
plt.show()
