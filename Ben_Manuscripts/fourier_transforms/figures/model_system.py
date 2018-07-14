#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt



# ld = layer displacement

ld = [0, 0.2, 0.4, 0.6, 0.8, 1]
ld_intensity = [25.41, 22.11, 15.56, 10.10, 6.36, 3.73]

plt.plot(ld, ld_intensity)
plt.xlabel('Allowed layer displacement / distance between layers')
plt.ylabel('R-$\pi$ Intensity')
plt.tight_layout()
plt.savefig('ld.png')

# td = thermal disorder
plt.figure()
td_percent = [-100, -50, 0, 50, 100]
td_intensity = [31.34, 30.90, 25.41, 20.07, 13.20]
plt.plot(td_percent, td_intensity)
plt.xlabel('Percent change in thermal disorder')
plt.ylabel('R-$\pi$ Intensity')
plt.tight_layout()
plt.savefig('td.png')

# cl = correlation length
cl = [1, 5, 10, 20, 30, 40]
cl_intensity = [24.22, 24.35, 25.00, 25.78, 26.85, 27.61]

plt.figure()
plt.plot(cl, cl_intensity)
plt.xlabel('Correlation Length ($\AA$)')
plt.ylabel('R-$\pi$ Intensity')
plt.tight_layout()
plt.savefig('cl.png')

# ic = initial configuration

plt.figure()
IC = [-100, -50, 0, 50, 100]
IC_intensity = [266.1, 83.2, 25.78, 10.62, 5.06]
plt.plot(IC, IC_intensity)
plt.xlabel('Percent change in initial configuration noise')
plt.ylabel('R-$\pi$ Intensity')
plt.savefig('IC.png')
plt.tight_layout()

plt.show()
