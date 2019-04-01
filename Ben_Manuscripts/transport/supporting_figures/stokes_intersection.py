#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql

water_size = 0.151391485766 / 2 # from ACH cubic box. std: 0.000407801323403

r_theoretical = np.linspace(0.125, 0.34, 100)
#r_theoretical = np.linspace(min(r), 10, 1000)

wt = 10
# values for methanol
if wt == 10:
	msd_max = 2.836 
	msd_lower = 2.378
	msd_upper = 3.353
elif wt == 5:
	msd_max = 0.111 
	msd_lower = 0.081
	msd_upper = 0.139
rmax = 0.283 / 2

# Corrected version first
alpha = water_size / r_theoretical # radius of solute (residue) to solvent (water) 
f = ((3 * alpha / 2) + (1 / (1 + alpha)))**-1  # f as function of alpha
alpha_max = water_size / rmax  # alpha for solute with max msd
fmax = ((3 * alpha_max / 2) + (1 / (1 + alpha_max)))**-1  # f for solute with max msd
a = fmax * msd_max * rmax  # make the curve pass through the max msd point
a_error_upper = fmax * msd_upper * rmax
a_error_lower = fmax * msd_lower * rmax

y_gierer_wirtz = a / (r_theoretical * f)
#y_gw_lower = a_error_lower / (r_theoretical * f)
#y_gw_upper = a_error_upper / (r_theoretical * f)

plt.plot(r_theoretical, y_gierer_wirtz, color='blue', label='Gierer and Wirtz correction')
#plt.fill_between(r_theoretical, y_gw_lower, y_gw_upper, alpha=opacity, color='blue')


# Make Stokes-Einstein intersect with correction at r_intersect

r_intersects = [0.5, 1.5, 2.5, 3.5, 100] # where we want stokes-einstein and corrected stoke-einstein to intersect

for r_intersect in r_intersects:
	alpha_intersect = water_size / r_intersect
	f_intersect = ((3 * alpha_intersect / 2) + (1 / (1 + alpha_intersect)))**-1
	y_intersect = a / (f_intersect * r_intersect)
	#y_intersect_upper =  a_error_upper / (f_intersect * r_intersect)
	#y_intersect_lower =  a_error_lower / (f_intersect * r_intersect)

	y_stokes = y_intersect * r_intersect / r_theoretical
	#y_stokes_upper = y_intersect_upper * r_intersect / r_theoretical
	#y_stokes_lower = y_intersect_lower * r_intersect / r_theoretical

	plt.plot(r_theoretical, y_stokes, '--', label='Critical Radius: %.1f nm' % r_intersect)
	#plt.fill_between(r_theoretical, y_stokes_lower, y_stokes_upper, color='black', alpha=opacity)

opacity = 0.1

# Plot formatting
plt.xlabel('Radius (nm)', fontsize=14)
plt.ylabel('MSD ($nm^2$)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.legend(fontsize=14)
plt.savefig('stokes_intersection_sensitivity.pdf')
plt.show()
