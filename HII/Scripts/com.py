#!/usr/bin/python

import numpy as np
from Get_Positions import get_positions
import Atom_props
import matplotlib.pyplot as plt

# # p, v, trj_times, box = get_positions('wiggle_traj_temp.gro', 'sys', 'HII', 'no')
#
#
# pos = np.load('pos')
# vel = np.load('vel')
# id = np.load('id')
#
# nT = np.shape(pos)[2]
# nmol = np.shape(pos)[1]
# kb = 1.38 * 10 ** -23
# Nf = 3*nmol - 3
# NA = 6.022 * 10 ** 23
#
# T_array = np.zeros([nT])
#
# for i in range(1, nT):
#     T_sum = 0
#     for j in range(nmol):
#         vs = vel[:, j, i]
#         vi = np.sqrt(np.dot(vs, vs)) * 1000
#         T_sum += ((Atom_props.mass[id[0][i]] / (1000*NA)) * vi ** 2)
#     T_array[i] = T_sum / (kb * 5*Nf)
#
# x = np.linspace(0, 10, len(T_array))
# print np.mean(T_array)
# plt.plot(x, T_array)
# plt.show()
#
# # f1 = open('pos', 'w')
# # f2 = open('vel', 'w')
# # np.save(f1, p)
# # np.save(f2, v)
# # f1.close()
# # f2.close()
#
# exit()
f = open('veloc.xvg', 'r')

a = []
for line in f:
    a.append(line)

velocities = np.zeros([3, 1888])

for i in range(0, len(a)):
    field1_start = 0
    while a[i][field1_start] == ' ':
            field1_start += 1
    field1_end = field1_start
    while a[i][field1_end] != ' ':
            field1_end += 1
    field2_start = field1_end
    while a[i][field2_start] == ' ':
            field2_start += 1
    field2_end = field2_start
    while a[i][field2_end] != ' ':
        field2_end += 1
    field3_start = field2_end
    while a[i][field3_start] == ' ':
        field3_start += 1
    field3_end = field3_start
    while a[i][field3_end] != ' ':
        field3_end += 1
    field4_start = field3_end
    while a[i][field4_start] == ' ':
        field4_start += 1

    velocities[0, i] = float(a[i][field2_start:field2_end])
    velocities[1, i] = float(a[i][field3_start:field3_end])
    velocities[2, i] = float(a[i][field4_start:len(a[0])])

mag_velocities = np.zeros([1888])
for i in range(0, len(a)):
    x = velocities[:, i]
    mag_velocities[i] = np.sqrt(np.dot(x, x))

x = np.linspace(0, 750, 1888)
plt.plot(x, mag_velocities*1000)
plt.show()
avg = np.mean(mag_velocities) * 1000
print avg
# print np.std(mag_velocities)
kb = 1.38 * 10 ** -23  # J/K
NA = 6.022 * 10 ** 23
MW = .86465  # kg/mol
m = (480 / NA) * MW

T = (1/(3*kb))*m*avg**2
print 1/T
# for i in range(0, len(a)):
#     velocities[0, i] = float(a[i][8:20])
#
# print velocities[0, 0]
