#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql

path = "/home/bcoscia/PycharmProjects/LLC_Membranes/LLC_Membranes/analysis/"
connection = sql.connect('%s/size.db' % path)
crsr = connection.cursor()
wt = 10  

colors = ['blue', 'red', 'xkcd:green', 'xkcd:orange', 'xkcd:goldenrod']
color_dict = {'MET':colors[0], 'ETH':colors[0], 'PR':colors[0], 'BUT':colors[0], 'GCL':colors[1], 'PG':colors[1], 'GLY':colors[1], 'TET':colors[1], 'RIB':colors[1], 'ACH':colors[2], 'URE':colors[2], 'ACN':colors[2], 'ATO':colors[2], 'SOH':colors[3], 'DMP':colors[3], 'DMS':colors[3], 'THF':colors[4], 'PCB':colors[4], 'EAC':colors[4], 'DMF':colors[4]}

simple_alcohols = False 
diols = False 
ketones = True 
sulfur = False 
nondonors = False
if simple_alcohols:
	restrict = ['MET', 'ETH', 'PR', 'BUT']
	shiftx = {'MET':0.005, 'ETH':0.005, 'PR':0.005, 'BUT':0.005}
	shifty = {'MET':0.03, 'ETH':0.03, 'PR':0.03, 'BUT':0.03}
elif diols:
	restrict = ['GCL', 'PG', 'GLY', 'TET', 'RIB']
	shiftx = {'GCL':0.005, 'PG':0.005, 'GLY':0.005, 'TET':0.005, 'RIB':0.005}
	shifty = {'GCL':0.03, 'PG':0.03, 'GLY':0.03, 'TET':0.03, 'RIB':0.03}
elif ketones:
	restrict = ['ACH', 'URE', 'ACN', 'ATO']
	shiftx = {'ACH':-0.03, 'URE':0.005, 'ACN':0.005, 'ATO':0.005}
	shifty = {'ACH':0.05, 'URE':0.03, 'ACN':0.03, 'ATO':0.03}
elif sulfur:
	restrict = ['SOH', 'DMP', 'DMS']
	shiftx = {'SOH':0.005, 'DMP':0.005, 'DMS':0.005}
	shifty = {'SOH':0.03, 'DMP':0.03, 'DMS':0.03}
elif nondonors:
	restrict = ['THF', 'PCB', 'EAC', 'DMF']
	shiftx = {'THF':-0.02, 'PCB':0.005, 'EAC':0.005, 'DMF':0.005}
	shifty = {'THF':0.075, 'PCB':0.03, 'EAC':0.03, 'DMF':-0.15}
else:
	restrict = []

command = "SELECT name, end_to_end, end_to_end_std FROM size WHERE nsolute = 6"

if restrict:
	command += ' and ('
	for i, n in enumerate(restrict):
		command += " name = '%s'" % n 
		if i < len(restrict) - 1:
			command += " or"
	command += ')'

radii = crsr.execute(command).fetchall()

rnames = [i[0] for i in radii]
r = np.array([i[1] for i in radii])
r_std = np.array([i[2] for i in radii])
connection = sql.connect('%s/../timeseries/msd.db' % path)
crsr = connection.cursor()

command = "SELECT name, MD_TAMSD, MD_TAMSD_CI_lower, MD_TAMSD_CI_upper, sim_length FROM msd WHERE wt_water = %d and name != 'HOH'" % wt

if wt == 10:
	command += "and penalty = 0.25"

if restrict:
	command += ' and ('
	for i, n in enumerate(restrict):
		command += " name = '%s'" % n 
		if i < len(restrict) - 1:
			command += " or"
	command += ')'

msds = crsr.execute(command).fetchall()

names = np.array([i[0] for i in msds])
msd = np.array([i[1] for i in msds])
msd_lower = -np.array([i[2] for i in msds]) + msd
msd_upper = np.array([i[3] for i in msds]) - msd
sim_length = np.array([i[4] for i in msds])

msd *= (400 / (sim_length * 0.4))

ordered = np.argsort(msd)

names = names[ordered].tolist()
msd = msd[ordered]
msd_lower = msd_lower[ordered]
msd_upper = msd_upper[ordered]
sim_length = sim_length[ordered]

scatter_colors = [color_dict[i] for i in names]

r_order = np.array([rnames.index(i) for i in names])

#plt.scatter(msd, r[r_order])
plt.scatter(r[r_order], msd, c=scatter_colors, s=50, zorder=3, edgecolors='black')
plt.errorbar(r[r_order], msd, yerr=(msd_lower, msd_upper), xerr=r_std[r_order], fmt='none', elinewidth=0.5, ecolor='black', marker="none")
if restrict:
	for i, x in enumerate(msd):
		plt.text(r[r_order][i] + shiftx[names[i]], x + shifty[names[i]], names[i])


# Theoretical Calcuations
water_size = 0.151391485766 # from ACH cubic box. std: 0.000407801323403

r_theoretical = np.linspace(0.25, 0.68, 100)
#r_theoretical = np.linspace(min(r), 10, 1000)

# values for methanol
if wt == 10:
	msd_max = 2.836 
	msd_lower = 2.378
	msd_upper = 3.353
elif wt == 5:
	msd_max = 0.111 
	msd_lower = 0.081
	msd_upper = 0.139
rmax = 0.283

# Corrected version first
alpha = water_size / r_theoretical # radius of solute (residue) to solvent (water) 
f = ((3 * alpha / 2) + (1 / (1 + alpha)))**-1  # f as function of alpha
alpha_max = water_size / rmax  # alpha for solute with max msd
fmax = ((3 * alpha_max / 2) + (1 / (1 + alpha_max)))**-1  # f for solute with max msd
a = fmax * msd_max * rmax  # make the curve pass through the max msd point
a_error_upper = fmax * msd_upper * rmax
a_error_lower = fmax * msd_lower * rmax

# Make Stokes-Einstein intersect with correction at r_intersect
r_intersect = 5 # where we want stokes-einstein and corrected stoke-einstein to intersect

alpha_intersect = water_size / r_intersect
f_intersect = ((3 * alpha_intersect / 2) + (1 / (1 + alpha_intersect)))**-1
y_intersect = a / (f_intersect * r_intersect)
y_intersect_upper =  a_error_upper / (f_intersect * r_intersect)
y_intersect_lower =  a_error_lower / (f_intersect * r_intersect)

y_stokes = y_intersect * r_intersect / r_theoretical
y_stokes_upper = y_intersect_upper * r_intersect / r_theoretical
y_stokes_lower = y_intersect_lower * r_intersect / r_theoretical

y_gierer_wirtz = a / (r_theoretical * f)
y_gw_lower = a_error_lower / (r_theoretical * f)
y_gw_upper = a_error_upper / (r_theoretical * f)

opacity = 0.1

plt.plot(r_theoretical, y_stokes, '--', color='black', label='Stokes-Einstein')
plt.fill_between(r_theoretical, y_stokes_lower, y_stokes_upper, color='black', alpha=opacity)
plt.plot(r_theoretical, y_gierer_wirtz, '--', color='blue', label='Gierer and Wirtz correction')
plt.fill_between(r_theoretical, y_gw_lower, y_gw_upper, alpha=opacity, color='blue')

# Plot formatting
plt.xlabel('Radius (nm)', fontsize=14)
plt.ylabel('MSD ($nm^2$)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.legend(fontsize=14)
if restrict:
	if simple_alcohols:
		plt.savefig('msd_radius_simple_alcohols_%dwt.pdf' % wt)
	if diols:
		plt.savefig('msd_radius_diols_%dwt.pdf' % wt)
	if ketones:
		plt.savefig('msd_radius_ketones_%dwt.pdf' % wt)
	if sulfur:
		plt.savefig('msd_radius_sulfur_%dwt.pdf' % wt)
	if nondonors:
		plt.savefig('msd_radius_nondonors_%dwt.pdf' % wt)

else:
	plt.savefig('msd_radius_%dwt.pdf' % wt)
plt.show()
