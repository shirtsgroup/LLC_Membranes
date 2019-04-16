#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import names

residues = np.array(['ACH', 'ACN', 'ATO', 'BUT', 'DMF', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE'])

length = 1000 # ns

# within 0.65 nm of pore center (random number)
# nhops = np.array([976, 979, 387, 416, 706, 650, 948, 672, 969, 1526, 978, 1729, 308, 1157, 515, 696, 808, 1001, 407, 1082])  # use hop_location.py
# frac_time_spent = np.array([0.54, 0.39, 0.31, 0.36, 0.43, 0.41, 0.37, 0.45, 0.37, 0.50, 0.65, 0.36, 0.40, 0.59, 0.40, 0.67, 0.39, 0.66, 0.42, 0.39])

nhops_total = np.array([1290, 1629, 864, 1125, 1310, 1111, 1716, 1085, 1761, 2166, 1275, 2659, 875, 1568, 925, 840, 1263, 1208, 771, 1671])

# within 0.512 nm of pore center (O3, O4 COM)
nhops = np.array([648, 721, 216, 252, 419, 468, 486, 398, 649, 1189,706, 1143, 126, 778, 291, 840, 605, 785, 286, 820])
#                           ACH   ACN   ATO   BUT   DMF   DMP   DMS   EAC   ETH   GCL   GLY   MET   PCB   PG    PR    RIB   SOH   TET   THF   URE
frac_time_spent = np.array([0.34, 0.25, 0.18, 0.20, 0.24, 0.28, 0.21, 0.24, 0.20, 0.30, 0.45, 0.22, 0.23, 0.40, 0.20, 0.50, 0.24, 0.49, 0.26, 0.26])

# within 0.73 nm of pore center (phenyl ring COM)
#                  ACH   ACN   ATO  BUT  DMF  DMP  DMS   EAC  ETH   GCL   GLY   MET   PCB  PG    PR   RIB  SOH  TET   THF  URE
#nhops = np.array([1065, 1226, 461, 481, 840, 700, 1070, 746, 1196, 1653, 1045, 1929, 372, 1235, 608, 741, 945, 1054, 478, 1185])
#                            ACH   ACN   ATO   BUT   DMF   DMP   DMS   EAC   ETH   GCL   GLY   MET   PCB   PG    PR    RIB   SOH   TET   THF   URE
#frac_time_spent = np.array([0.69, 0.55, 0.44, 0.49, 0.60, 0.51, 0.53, 0.62, 0.50, 0.63, 0.76, 0.49, 0.55, 0.68, 0.56, 0.79, 0.49, 0.76, 0.55, 0.50])

hopfreq = nhops / length
hopfreq_total = nhops_total / length

hopfreq /= 24
hopfreq_total /= 24

freq_total = False 
freq = False 
frac = True 

if freq:
	ordered = np.argsort(hopfreq)[::-1]
elif frac:
	ordered = np.argsort(frac_time_spent)[::-1]
elif freq_total:
	ordered = np.argsort(hopfreq_total)[::-1]

labels = np.array([names.abbreviation[r] for r in residues])[ordered]
colors = np.array([names.color_dict[r] for r in residues])[ordered]

index = np.arange(len(residues))
fig, ax = plt.subplots()
ax.tick_params(labelsize=14)
ax.set_xticks(index)
ax.set_xticklabels(labels, fontsize=14)
[x.set_color(colors[i]) for i, x in enumerate(plt.gca().get_xticklabels())]
plt.xticks(rotation=90)
if freq:
	ax.bar(index, hopfreq[ordered])
	ax.set_ylabel('Hop frequency (ns$^{-1}$)', fontsize=14)
	savename = 'hopfreq.pdf'
elif frac:
	ax.bar(index, frac_time_spent[ordered])
	ax.set_ylabel('Fraction of time spent in pore region', fontsize=14)
	savename = 'frac_time_spent.pdf'
elif freq_total:
	ax.bar(index, hopfreq_total[ordered])
	ax.set_ylabel('Hop frequency (ns$^{-1}$)', fontsize=14)
	savename = 'hopfreq_total.pdf'

plt.tight_layout()
plt.savefig(savename)
plt.show()
