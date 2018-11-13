#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


for i in [0, 75, 1]:
	plt.figure()
	data = np.load('columns_sr%s.npz' %i)
	x = data['freq_z']
	y = data['slice']
	if i == 0:
		y *= 1.24
	plt.plot(x, y, linewidth=4)
	plt.xlim(-2.25, 2.25)
	plt.xticks([-2, -1, 0, 1, 2])
	plt.xlabel('$q_r\ (\AA^{-1})$', fontsize=22)
	plt.ylabel('Intensity', fontsize=22)
	plt.gcf().get_axes()[0].tick_params(labelsize=18)
	plt.tight_layout()
	if i == 1:
		plt.savefig('../../structure_paper/figures/sf_qy_sr100.pdf')
	else:
		plt.savefig('../../structure_paper/figures/sf_qy_sr%s.pdf' % i)

plt.show()
exit()

displacement_fraction = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
# same rotation - r
#intensity = [10.46, 12.04, 20.55, 33.98, 60.20, 81.05, 122.53, 140.58, 181.90, 181.92, 193.12] 
intensity = [10.25, 11.75, 20.43, 34.36, 58.13, 87.40, 118.13, 149.90, 176.22, 193.72, 200] 

plt.plot(displacement_fraction, intensity)
plt.xlabel('Allowed column displacement / distance between layers')
plt.ylabel('R-$\pi$ Intensity')
plt.savefig('column_displacement.png')
plt.show()
