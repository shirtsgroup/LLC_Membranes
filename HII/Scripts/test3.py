#!/usr/bin/python

import numpy as np
import time
import matplotlib.pyplot as plt
import Radial_int_pixels


def rij(coords, identity):

    atoms = np.shape(coords)[1]
    rij = np.zeros([atoms**2 - atoms, 3], dtype=object)

    count = 0
    for i in range(atoms):
        for k in range(atoms):
            if i != k:
                rij[count, 0] = np.linalg.norm(coords[:, i] - coords[:, k])
                rij[count, 1] = identity[i]
                rij[count, 2] = identity[k]
                count += 1

    return rij

identity = np.zeros([4], dtype=object)
identity[0] = 'C'
identity[1] = 'O'
identity[2] = 'H'
identity[3] = 'NA'

coords = np.zeros([3, 4])
coords[:, 0] = [0, 0, 0]
coords[:, 1] = [1, 0, 0]
coords[:, 2] = [1, 1, 0]
coords[:, 3] = [0, 1, 0]


start = time.time()
rij = rij(coords, identity)
stop = time.time()
print 'rij calculated in %s seconds' % (stop - start)
print rij[0][0]
print type(rij[0][0])

print rij
exit()
so = np.array([1, 0])
print so
s = np.array([np.sqrt(3)/2, 0.5])
print s
r = np.array([-1/np.sqrt(2), 1/np.sqrt(2)])
print r

diff = s - so
print diff

print np.dot(-r, diff)
exit()

def gauss(x, mu, sig):
    return 1/(np.sqrt(2*np.pi) * sig) * np.exp(-(x - mu)**2 / (2*sig**2))


dist_mat, max = Radial_int_pixels.dist_mat(1024)

circle = np.zeros([1024, 1024])

for i in range(1024):
    for j in range(1024):
        if 300 <= dist_mat[i, j] <= 750:
            if 400 <= i <= 600:
                circle[i, j] = 100*gauss((dist_mat[i, j] - 525)/50, 0, 0.5)
        if 800 <= dist_mat[i, j] <= 1200:
            circle[i, j] = 50*gauss((dist_mat[i, j] - 1000)/50, 0, 0.5)

plt.figure(1)
im = plt.imshow(circle, cmap='Greys', interpolation='none', vmin=0, vmax=100)

plt.figure(2)
x, intensities = Radial_int_pixels.radial_int(circle, 1024, 1)
plt.plot(x, intensities)

plt.show()

exit()
f = open('wiggle_traj.gro', 'r')

a = []
for line in f:
    a.append(line)

nT = 751
no_ion = 480
dim = 3

# while i < len(a):
#     if 'trjconv' in a[i]:
#         frames += 1
#     i += 1

x = np.zeros([dim, no_ion, nT])

start = time.time()
pos = []
for i in range(len(a)):
    if 'NA' in a[i]:
        pos.append([float(a[i][20:28]), float(a[i][28:36]), float(a[i][36:44])])
end = time.time()
print 'Positions got in %s seconds' % (end - start)
start = time.time()
for i in range(nT):
    for j in range(no_ion):
        x[:, j, i] = pos[i*no_ion + j]
end = time.time()
print 'Positions arranged into "x" vector in %s seconds' % (end - start)

import Diffusivity
import matplotlib.pyplot as plt
# start = time.time()
# MSDs = np.zeros([nT, no_ion])
# ions = np.shape(x)[1]
# nT -= 1
# for j in range(ions):
#     itau = 0  # counts up to the number of trajectory frames
#     while itau < nT:
#         ncount = (nT-itau)
#         for t in range(nT-itau):  # run for total number of time points from itau to nT
#             xo = x[:, j, t + itau] - x[:, j, t]
#             MSDs[itau, j] += np.dot(xo, xo)  # Should there be a division by 't' in here?
#         MSDs[itau, j] /= ncount
#         itau += 1
# end = time.time()
# print 'MSDs calculated in %s seconds' % (end - start)
# avg = np.zeros([nT])
# for i in range(nT):
#     avg[i] = np.mean(MSDs[i, :])
#
# print len(avg)
# print avg
# import Poly_fit
# t = np.linspace(0, 300000, nT)
# # kb = 1.38 * 10 ** -23
# # y_fit, r_squared, std2, coeff = Poly_fit.poly_fit(t[0:150], avg[0:150], 1)
# # conv = 157899298.191  #
# # ic_cd = (coeff[1]*(1*10**12)*conv)/(2*kb*float(300))
# # print 'Ionic Conductivity: %s' % ic_cd
# plt.plot(t[0:len(avg)], avg)
# # plt.plot(t[0:150], y_fit)
# # plt.suptitle('Collective Diffusion Model')
# # plt.title('Ionic Conductivity = %s' % ic_cd)
# # plt.ylabel('MSD (e$^2$)')
# # plt.xlabel('Time (ps)')
# plt.show()

Nbootstraps = 2000
frontfrac = 0.14
fracshow = 0.4
d = 3
dt = 400
D, MSD, endMSD, limits = Diffusivity.dconst(x, nT, Nbootstraps, frontfrac, fracshow, d, dt)
errorevery = int(np.ceil(fracshow*nT/100.0))  # plot only 100 bars total
plt.errorbar(dt*np.array(range(0,endMSD)),MSD[:endMSD],yerr=[limits[0,:endMSD],limits[1,:endMSD]],errorevery=errorevery)
plt.ylabel('MSD (nm^2)')
plt.xlabel('time (ps)')
plt.suptitle('Stackoverflow/Shirts algorithm MSD results')
plt.title('Diffusivity = %s m$^2$/s' %D)
plt.show()