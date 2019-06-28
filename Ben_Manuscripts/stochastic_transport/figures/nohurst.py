#!/usr/bin/env python

# comparison of MSDs if fractional brownian motion is ignored

import numpy as np
import sqlite3 as sql
import names
import matplotlib.pyplot as plt

connection = sql.connect("../../../timeseries/msd.db")
crsr = connection.cursor()
query = "SELECT DISTINCT name, MD_MSD, MD_MSD_CI_lower, MD_MSD_CI_upper from msd ORDER BY MD_MSD DESC"
output = crsr.execute(query).fetchall()

bar_width = 0.4
index = np.arange(len(output))

ordered_names = [i[0] for i in output]#, dtype=object)
md_msd = np.array([i[1] for i in output])
md_msd_lower = -np.array([i[2] for i in output]) + md_msd
md_msd_upper = np.array([i[3] for i in output]) - md_msd

n = ['ACH', 'ACN', 'ATO', 'BUT', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE', 'DMF']
indices = np.array([n.index(i) for i in ordered_names])

msd = np.zeros([len(ordered_names), 3])
tamsd = np.zeros([len(ordered_names), 3])

#ACH
msd[0, :] = [3.23, 2.89, 2.59]
tamsd[0, :] = [0.98, 0.85, 1.09]

#ACN
msd[1, :] = [3.09, 2.74, 3.44]
tamsd[1, :] = [1, .9, 1.12]

# ATO
msd[2, :] = [1.34, 1.20, 1.51]
tamsd[2, :] = [0.43, 0.38, 0.49]

# BUT
msd[3, :] = [0.85, 0.75, 0.97]
tamsd[3, :] = [0.19, 0.17, 0.21]

# DMP
msd[4, :] = [0.50, 0.53, 0.68]
tamsd[4, :] = [0.17, 0.15, 0.20]

# DMS
msd[5, :] = [2.37, 2.11, 2.63]
tamsd[5, :] = [0.77, 0.68, 0.86]

# EAC
msd[6, :] = [1.01, 0.90, 1.15]
tamsd[6, :] = [0.28, 0.24, 0.31]

# ETH
msd[7, :] = [5.10, 4.57, 5.59]
tamsd[7, :] = [1.90, 1.69, 2.18]

# GCL
msd[8, :] = [11.26, 10.02, 12.67]
tamsd[8, :] = [3.85, 3.56, 4.22]

# GLY
msd[9, :] = [1.71, 1.52, 1.91]
tamsd[9, :] = [.54, .48, .6]

# MET
msd[10, :] = [54.68, 49.48, 60.08]
tamsd[10, :] = [20.60, 18.95, 22.40]

# PCB
msd[11, :] = [0.59, 0.53, 0.65]
tamsd[11, :] = [0.16, 0.15, 0.19]

# PG -- redo
msd[12, :] = [2.04, 1.78, 2.30]
tamsd[12, :] = [0.71, 0.63, 0.81]

# PR
msd[13, :] = [1.54, 1.36, 1.75]
tamsd[13, :] = [.42, .37, .48]

# RIB
msd[14, :] = [0.68, 0.60, 0.75]
tamsd[14, :] = [0.28, 0.24, 0.32]

# SOH
msd[15, :] = [4.48, 3.94, 5.02]
tamsd[15, :] = [1.35, 1.23, 1.51]

# TET
msd[16, :] = [1.25, 1.10, 1.40]
tamsd[16, :] = [0.39, 0.34, 0.43]

# THF
msd[17, :] = [1.15, 1.01, 1.30]
tamsd[17, :] = [0.30, 0.26, 0.34]

# URE
msd[18, :] = [9.62, 8.62, 10.68]
tamsd[18, :] = [3.53, 3.23, 3.83]

# DMF
msd[19, :] = [0.64, 0.55, 0.72]
tamsd[19, :] = [0.19, 0.17, 0.21]

msd[:, 1] = -msd[:, 1] + msd[:, 0]
msd[:, 2] = msd[:, 2] - msd[:, 0]

fig, ax = plt.subplots(figsize=(12, 7))
opacity = 0.8

ax.bar(index, md_msd, bar_width, alpha=opacity, color='red', yerr=(md_msd_lower, md_msd_upper), label='MD')
ax.bar(index + bar_width, msd[indices, 0], bar_width, alpha=opacity, color='blue', yerr=(msd[indices, 1], msd[indices, 2]), label='CTRW')

plt.show()

