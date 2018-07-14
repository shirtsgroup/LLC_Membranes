#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

displacement_fraction = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
# same rotation - r
#intensity = [10.46, 12.04, 20.55, 33.98, 60.20, 81.05, 122.53, 140.58, 181.90, 181.92, 193.12] 
intensity = [10.25, 11.75, 20.43, 34.36, 58.13, 87.40, 118.13, 149.90, 176.22, 193.72, 200] 

plt.plot(displacement_fraction, intensity)
plt.xlabel('Allowed column displacement / distance between layers')
plt.ylabel('R-$\pi$ Intensity')
plt.savefig('column_displacement.png')
plt.show()
