#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#labels = ['S\nxlink', 'PD\nxlink', 'S', 'PD']

# Sandwich xlink, PD xlink, Sandwich, PD
#IC = [12.5 * 10 ** -5, 8.03 * 10 **-5, 16*10**-5, 8.1 * 10 **-5]
#IC_error = [.11 * 10 ** -5, .06 * 10 ** -5, 0.2 * 10 ** -5, 0.1 * 10 ** -5]
#p2p = [ 4.1, 4.1, 4.18, 4.17]
#p2p_error = [0.08, 0.10, 0.13, 0.12]
#dbwl = [4.378, 4.22, 4.33, 4.27]
#dbwl_error = [0.1, 0.07, 0.03, 0.04]
#correlation_length = [4, 8.77, 4.5, 15.5]
#correlation_length_error = [0.62, 2.02, 0.4, 1.9]

# different order
labels = ['S', 'S\nxlink', 'PD', 'PD\nxlink']
print(labels)

# Sandwich, Sandwich xlink, PD, PD xlink
IC = [16*10**-5, 12.5*10**-5, 8.1*10**-5, 8.03 * 10 **-5]
IC_error = [0.2*10**-5, .11 * 10 ** -5, 0.1*10**-5, .06*10**-5]
p2p = [4.18, 4.1, 4.17, 4.1]
p2p_error = [0.13, 0.08, 0.12, 0.10]
dbwl = [4.27, 4.378, 4.33, 4.22]
dbwl_error = [0.03, 0.1, 0.04, 0.07]
correlation_length = [4.5, 4, 15.5, 8.77]
correlation_length_error = [0.4, 0.62, 1.9, 2.02]



centers = [1, 2, 3, 4]
width = 0.9 
opacity = 0.9 
colors = ['xkcd:blue', 'xkcd:green', 'xkcd:red', 'xkcd:orange']

# plot data
f, axarr = plt.subplots(2, 2)
barlist1 = axarr[0, 0].bar(centers, IC, width, yerr=IC_error, alpha=opacity)
barlist2 = axarr[0, 1].bar(centers, p2p, width, yerr=p2p_error, alpha=opacity)
barlist3 = axarr[1, 0].bar(centers, dbwl, width, yerr=dbwl_error, alpha=opacity)
barlist4 = axarr[1, 1].bar(centers, correlation_length, width, yerr=correlation_length_error, alpha=opacity)
for i in range(4):
	barlist1[i].set_color(colors[i])
	barlist2[i].set_color(colors[i])
	barlist3[i].set_color(colors[i])
	barlist4[i].set_color(colors[i])

# format axes
# IC plot
axarr[0, 0].set_yscale('log')
axarr[0, 0].set_ylim(1*10**-5, 5*10**-4)
axarr[0, 0].set_ylabel('Ionic Conductivity (S/m)')

# p2p plot
axarr[0, 1].set_ylim(3.5, 4.5)
axarr[0, 1].set_ylabel('Pore spacing (nm)')

# dbwl plot
axarr[1, 0].set_ylim(4.0, 4.6)
axarr[1, 0].set_ylabel('Monomer stacking distance ($\AA$)')
axarr[1, 0].set_xticks(centers)
axarr[1, 0].set_xticklabels(labels)

# correlation length plot
axarr[1, 1].set_ylim(0, 20)
axarr[1, 1].set_ylabel('Correlation Length ($\AA$)')
axarr[1, 1].set_xticks(centers)
axarr[1, 1].set_xticklabels(labels)


# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
plt.tight_layout()
plt.savefig('xlink_charts.pdf')
plt.show()
