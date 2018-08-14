#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

a = []
with open('raw_SAXS.txt', 'r') as f:
	for line in f:
		a.append(line)

q = np.zeros([len(a)])
I = np.zeros([len(a)])

for i in range(len(a)):
	q[i] = float(a[i].split()[0])
	I[i] = float(a[i].split()[1])

start = 0
while q[start] < 0.1:
	start += 1

end = 0
while q[end] < 0.5:
	end += 1

ralkanes = 0.0150233972936
#I /= np.max(I[start:end])
I /= ralkanes
imax = np.where(I == np.max(I[start:end]))[0][0]
print(q[imax])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.annotate('$q*$ = 0.179 $\AA^{-1}$', (0.17, 1.025), fontsize=14)
ax.annotate('$\sqrt{3}q*$', (0.29, 0.5), fontsize=14)
ax.annotate('$\sqrt{4}q*$', (0.34, 0.3), fontsize=14)
plt.plot(q[start:end], I[start:end], linewidth=2.0)
plt.show()
plt.xlabel('$q$ ($\AA^{-1}$)', fontsize=14)
plt.ylabel('Normalized Intensity', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.ylim(0.1, 1.1)
ax.set_aspect('0.5')
plt.tight_layout()
#plt.savefig('SAXS.png')
plt.show()
