#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.set_aspect(1)
x = np.zeros([7, 2])
lxy = 4.5  # distance between hexagon points
lz = 3.7
theta = 60 * (np.pi / 180)
x[1, :] = [lxy, 0]
x[2, :] = [-lxy*np.cos(theta), lz*np.sin(theta)]
x[3, :] = [lxy*np.cos(theta), lz*np.sin(theta)]
x[4, :] = [-lxy*np.cos(theta), -lz*np.sin(theta)]
x[5, :] = [lxy*np.cos(theta), -lz*np.sin(theta)]
x[6, :] = [-lxy, 0]
x[4, :] = [lxy + lxy*np.cos(theta), lz*np.sin(theta)]
x[5, :] = [0, 2*lz*np.sin(theta)]
x[6, :] = [lxy, 2*lz*np.sin(theta)]
x -= x[3, :]
x_real = np.linspace(-10, 10, 100)

m_real_1 = 2*np.sin(theta) 
m_real_2 = -2*np.sin(theta)

#real = np.zeros([2, 6, 100])  # [(x,y), line, points making up line]
#real[0, 0, :] = x_real
#real[1, 0, :] = m_real_1*x_real
#real[0, 1, :] = x_real + l
#real[1, 1, :] = m_real_1*x_real
#real[0, 2, :] = x_real - l 
#real[1, 2, :] = m_real_1*x_real
#real[0, 3, :] = x_real
#real[1, 3, :] = m_real_2*x_real
#real[0, 4, :] = x_real + l
#real[1, 4, :] = m_real_2*x_real
#real[0, 5, :] = x_real - l
#real[1, 5, :] = m_real_2*x_real

#for i in range(6):
#	plt.plot(real[0, i, :], real[1, i, :], '--', color='black')

plt.scatter(x[:, 0], x[:, 1], linewidth=5, color='black')

plt.gcf().get_axes()[0].set_aspect('equal')
plt.ylim(-5, 5)
plt.xlim(-5, 5)
plt.xlabel('x', fontsize=18)
plt.ylabel('y', fontsize=18)
plt.tight_layout()
plt.savefig('hexagonal_packing_realspace.png')

m_reciprocal_1 = - (1/m_real_1)  # lines perpendicular to real space lines
m_reciprocal_2 = - (1/m_real_2)
#l_reciprocal = 2*np.pi / l  # reciprocal lattice spacing

x_reciprocal = np.linspace(-5, 5, 100)
recip = np.zeros([2, 6, 100])  # [(x,y), line, points making up line]
#recip[0, 0, :] = x_reciprocal
#recip[1, 0, :] = m_reciprocal_1*x_reciprocal
#recip[0, 1, :] = x_reciprocal + l_reciprocal/np.sin(theta)
#recip[1, 1, :] = m_reciprocal_1*x_reciprocal
#recip[0, 2, :] = x_reciprocal - l_reciprocal/np.sin(theta)
#recip[1, 2, :] = m_reciprocal_1*x_reciprocal
#recip[0, 3, :] = x_reciprocal
#recip[1, 3, :] = m_reciprocal_2*x_reciprocal
#recip[0, 4, :] = x_reciprocal + l_reciprocal/np.sin(theta)
#recip[1, 4, :] = m_reciprocal_2*x_reciprocal
#recip[0, 5, :] = x_reciprocal + 2* l_reciprocal/np.sin(theta)
#recip[1, 5, :] = m_reciprocal_2*x_reciprocal

#for i in [0, 3]:
#	plt.plot(recip[0, i, :], recip[1, i, :], '--', color='red')

#plt.xlim(-5, 5)
#plt.ylim(-5, 5)

#plt.figure()
hcp = np.zeros([8, 2])
hcp[1::2, 0] = lxy
for i in range(4):
	hcp[2*i : (2*(i + 1)), 1] = i*3.7
hcp[2:4, 0] += lxy/2
hcp[6:, 0] += lxy/2
hist_range = [[0, 2*lxy], [0, 3.7*4]] 
H, xedges, yedges = np.histogram2d(hcp[:, 0], hcp[:, 1], bins=100, range=hist_range)
X = [(xedges[i - 1] + xedges[i])/2 for i in range(1, len(xedges))]
Y = [(xedges[i - 1] + xedges[i])/2 for i in range(1, len(xedges))]
#plt.contourf(X, Y, H)

plt.figure()
xbin = xedges[1] - xedges[0]
ybin = yedges[1] - yedges[0]
ft = np.abs(np.fft.fftn(H)) ** 2
freq_x = np.fft.fftfreq(len(X), d=xbin)
freq_y = np.fft.fftfreq(len(Y), d=ybin)
ndx_x = np.argsort(freq_x)
ndx_y = np.argsort(freq_y)
freq_x = freq_x[ndx_x]
freq_y = freq_y[ndx_y]

ft = ft[ndx_x, :]
ft = ft[:, ndx_y]

plt.contourf(freq_x*2*np.pi, freq_y*2*np.pi, ft, cmap='Greys')
plt.xlim(-2.2, 2.2)
plt.ylim(-2.2, 2.2)
plt.xlabel('$q_x$', fontsize=18)
plt.ylabel('$q_z$', fontsize=18)
plt.gcf().get_axes()[0].set_aspect('equal')
plt.gcf().get_axes()[0].tick_params(labelsize=18)
plt.tight_layout()
plt.savefig('hexagonal_ft.pdf')
plt.show()

#fig, ax = plt.subplots()
#plt.scatter(0.698, -.403, linewidth=5)
#from matplotlib.patches import Circle
#circle = Circle([0.8061, -.4654], radius=l_reciprocal)
#ax.add_patch(circle)

print((180/np.pi)*np.arctan(m_reciprocal_1))
print((180/np.pi)*np.arctan(m_reciprocal_2))

#plt.xlim(-1.5, 1.5)
#plt.ylim(-1.5, 1.5)

#plt.show()
