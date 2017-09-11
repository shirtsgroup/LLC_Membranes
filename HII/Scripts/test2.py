#!/usr/bin/bash
import numpy as np
import matplotlib.pyplot as plt
import math
# import Poly_fit
import subprocess
import Get_Positions
import Periodic_Images
import Radial_int_pixels
import Atom_props
from pymbar import timeseries
import numpy as np
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import path
from scipy.spatial import ConvexHull
import numdensity

factor = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]
H = []

plt.figure(1)
for j in factor:

    x = np.linspace(0, 50, 51)
    d = np.zeros(51)

    for i in range(51):
        if i % 2 == 0:
            d[i] = 2 + j
        else:
            d[i] = 2 - j

    H.append(numdensity.entropy(x, d))

    plt.plot(x, d)


plt.figure(2)
plt.plot(factor, H)
plt.ylabel('Entropy')
plt.xlabel('"disorder"')
plt.show()

exit()
list = {}

list[1] = 2
print list[1]
exit()
# Get box vertices. Order is important!!! It needs to make a path that when traced from vertex to vertex is a box
box1 = np.zeros([4, 2])
box1[0, :] = [0, 0]
box1[1, :] = [0, 1]
box1[2, :] = [1, 1]
box1[3, :] = [1, 0]

box2 = np.zeros([4, 2])
box2[0, :] = [0.5, 0.5]
box2[1, :] = [0.5, 1.5]
box2[2, :] = [1.5, 1.5]
box2[3, :] = [1.5, 0.5]

import warnings
warnings.filterwarnings("ignore")  # Couldn't get this to work: http://docs.python.org/2/library/warnings.html#temporarily-suppressing-warnings


def ccw(A,B,C):
    """
    Taken from here: http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    """
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])


def intersect(A,B,C,D):
    """
    Return true if line segments AB and CD intersect
    Taken from here: http://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    """
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


def slope(pt1, pt2):

    m = (pt1[1] - pt2[1]) / (pt1[0] - pt2[0])
    b = pt1[1] - m * pt1[0]

    return m, b


def intersection(A,B,C,D):
    """
    :param A: point 1 in vector 1
    :param B: point 2 in vector 1
    :param C: point 1 in vector 2
    :param D: point 2 in vector 2
    :return: The point of intersection (if there is one) or else return False
    """

    if ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D):  # if this is true, the lines intersect

        m1, b1 = slope(A, B)
        m2, b2 = slope(C, D)

        if abs(m2) == inf:  # slope is infinite meaning pt1[x] = pt2[x]
            x_intersect = C[0]
            y_intersect = m1 * x_intersect + b1
        elif abs(m1) == inf:  # same but for other line
            x_intersect = A[0]
            y_intersect = m2 * x_intersect + b2
        else:  # a normal case
            x_intersect = (b1 - b2) / (m2 - m1)
            y_intersect = m1 * x_intersect + b1

        return [x_intersect, y_intersect]

    else:
        return False  # case where they do not intersect


def overlap_pts(shape1, shape2, sides=4):

    pts = np.zeros([sides*2, 2])  # the max number of pts making up the intersecting polygon is somewhere around this. I'll fix it when it fails
    count = 0
    for i in range(sides):
        for j in range(sides):
            pt = intersection(box1[i - 1, :], box1[i, :], box2[j - 1, :], box2[j, :])
            if pt:
                pts[count, :] = pt
                count += 1

    # xyverts1 = np.array(box1)
    # xyverts2 = np.array(box2)
    # create matplotlib path objects that defines the two benzene rings as polygon l in 2D
    p1 = path.Path(shape1)
    p2 = path.Path(shape2)

    for i in range(shape2.shape[0]):
        if p1.contains_point(shape2[i, :]):
            pts[count, :] = shape2[i, :]
            count += 1
        if p2.contains_point(shape1[i, :]):
            pts[count, :] = shape1[i, :]
            count += 1

    return pts[~np.all(pts == 0, axis=1)]  # get rid of excess zeros


def order_pts(pts):

    hull = ConvexHull(pts)

    pts_ordered = np.zeros(pts.shape)
    count = 0
    for vertex in hull.vertices:
        pts_ordered[count, :] = [pts[vertex, 0], pts[vertex, 1]]
        count += 1

    return pts_ordered


def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

pts = overlap_pts(box1, box2)

pts_ordered = order_pts(pts)

# now we have all the points making up the region which will

print PolyArea(pts_ordered[:, 0], pts_ordered[:, 1])

# print pts[:, 0]
# print pts[:, 1]
# print PolyArea(pts[:, 0], pts[:, 1])
exit()
# import random
#
# n = 100010
# xy = np.zeros([n, 2])
# y = np.zeros([n])
# for i in range(n):
#     xy[i, 0] = random.uniform(-1, 2)
#     xy[i, 1] = random.uniform(-1, 2)
#
# x = []
# y = []
#
# for i in range(n):
#     if p.contains_point(xy[i, :]):
#         x.append(xy[i, 0])
#         y.append(xy[i, 1])
#
# plt.scatter(x, y)
# plt.show()
exit()
n_sub = 5
nT_sub = 1000
sigma = 10
q = np.zeros([n_sub, nT_sub])
for i in range(n_sub):
    for j in range(nT_sub):
        m = np.random.randint(-1, 2)
        q[i, j] = m*sigma*np.random.randn()

for i in range(n_sub):
    plt.plot(np.arange(0,nT_sub), q[i, :])
plt.show()

msd = np.zeros([n_sub, nT_sub])
itau = 0  # counts up to the number of trajectory frames
while itau < nT_sub:
    ncount = (nT_sub-itau)
    for t in range(nT_sub-itau):  # run for total number of time points from itau to nT
        for i in range(n_sub):
            xo = q[i, t + itau] - q[i, t]
            msd[i, itau] += np.dot(xo, xo)  # Should there be a division by 't' in here?
    msd[:, itau] /= ncount
    itau += 1


for i in range(n_sub):
    plt.plot(np.arange(0,nT_sub), msd[i, :])
plt.title('Scaling of PLC$\gamma$1 domain')
plt.show()
exit()

# f = open('wiggle.gro', 'r')
# a = []
# for line in f:
#     a.append(line)
# f.close()
#
# pos = np.zeros([3, 69120])
# for i in range(2, len(a) - 1):
#     pos[0, i - 2] = float(a[i][20:28])
#     pos[1, i - 2] = float(a[i][28:36])
#     pos[2, i - 2] = float(a[i][36:44])
#
# f = open('pos_1_frame_0_images', 'w')
# np.save(f, pos)
# f.close()
#
# exit()
# # correlation time stuff for pores
# p2ps = np.load('p2ps')
# exclude = 4
# frames = np.shape(p2ps)[1]
#
# p2p_new = np.zeros([5, frames])
# count = 0
# for i in range(6):
#     if i != exclude:
#         p2p_new[count, :] = p2ps[i, :]
#         count += 1
#
# p2ps = p2p_new
#
# frames = np.shape(p2ps)[1]
# traj_pts = 400 * np.linspace(0, frames - 1, frames)
# nboot = 200
#
# t = 957
#
# g = timeseries.statisticalInefficiency(p2ps[0, t:])
# taus = []
# for i in range(5):
#     tau = timeseries.integratedAutocorrelationTime(p2ps[i, t:])
#     taus.append(tau)
#
# tau = int(max(taus))
#
# ind_trajectories = (frames - t) / tau
# total_trajectories = ind_trajectories * 5
# trajectories = np.zeros([tau, total_trajectories])
#
# for i in range(5):
#     for j in range(ind_trajectories):
#         trajectories[:, i * ind_trajectories + j] = p2ps[i, (t + j*tau):(t + (j + 1)*tau)]
#
# p2p_boot = np.zeros([nboot])
# for i in range(nboot):
#     p2p = 0
#     for j in range(tau):
#         T = ran.randrange(0, total_trajectories)  # pick a random trajectory from all the independent trajectories
#         P = ran.randrange(0, tau)  # choose a random point in that trajectory
#         p2p += trajectories[P, T]
#     p2p_boot[i] = p2p / tau
#
# p2p_avg = np.mean(p2p_boot)
# p2p_std = np.std(p2p_boot)
#
#
# print p2p_avg
# print p2p_std
# # t, g, Neff = timeseries.detectEquilibration(p2ps[0, :])
# # print t
# # print g
# exit()
# # atoms = 4
# # pos = np.zeros([3, 4])
# # pos[:, 0] = [1, 2, 3]
# # pos[:, 1] = [2, 3, 4]
# # pos[:, 2] = [3, 4, 5]
# # pos[:, 3] = [4, 5, 6]
# # print pos
# # for i in range(atoms):
# #     temp = pos[1, i]
# #     pos[1, i] = pos[2, i]
# #     pos[2, i] = temp
# #
# # print pos
# #
# # exit()
# Intensities = np.load('Intensities')
#
# mins = []
# for i in range(1024):
#     mins.append(min(Intensities[:, i]))
#
# Imin = min(mins)
# Imax = 10889207.9094
# q, theta, intensity = Radial_int_pixels.radial_int(Intensities, 1024, .0002, .0001,
#                                  1.18, 1.54)
#
# plt.figure(1)
# im = plt.imshow(Intensities, cmap='Greys', interpolation='none', vmin=Imin, vmax=Imax)
#
# plt.figure(2)
# plt.plot(q, intensity)
#
# plt.show()
#
# exit()
# f = open('/home/bcoscia/Documents/Gromacs/test.gro', 'r')
#
# a = []
# for line in f:
#     a.append(line)
#
# f.close()
#
# pos = np.zeros([3, 480, 1])
# count = 0
# i = 0
# while i < len(a):
#     if str.strip(a[i][5:10]) == 'NA':
#         pos[0, count, 0] = a[i][20:28]
#         pos[1, count, 0] = a[i][28:36]
#         pos[2, count, 0] = a[i][36:44]
#         count += 1
#     i += 1
#
# import Structure_char
#
# pcenters = Structure_char.avg_pore_loc(4, 1, pos, (480/4))
#
# distances = 6  # number of p2p distances to calculate. My algorithm isn't good enough for anything but six yet
# p2ps = Structure_char.p2p(pcenters, distances, 1)
# print p2ps
#
# exit()


pos, _, _, box = Get_Positions.get_positions('wiggle_traj2.gro', 'sys', 'HII', 'no')
print 'box and positions got'
f = open('box_array_thin', 'w')
np.save(f, box)
f.close()
f = open('pos_array_thin', 'w')
np.save(f, pos)
f.close()
# box = np.load('box_array')
# pos = np.load('pos_array612ns')

no_comp = np.shape(pos)[1]

f = open('wiggle_traj2.gro', 'r')
a = []
for line in f:
    a.append(line)
f.close()

id = np.zeros([1, no_comp], dtype=object)

count = 2

for j in range(no_comp):
    atom = str.strip(a[count][10:15])
    id[0, j] = atom
    count += 1

print 'id array made'
f = open('id_thin', 'w')
np.save(f, id)
f.close()

# id = np.load('identity_array612ns')
# id = np.load('id_plain')

for i in range(np.shape(id)[1]):
    if id[0, i].count('C') != 0:
        id[0, i] = 'C'
    elif id[0, i].count('H') != 0:
        id[0, i] = 'H'
    elif id[0, i].count('O') != 0:
        id[0, i] = 'O'

f = open('id_thin_plain', 'w')
np.save(f, id)
f.close()

f = open('Thin_periodic_trajectory.txt', 'w')
frames = 100
atoms = np.shape(id)[1]

for i in range(frames):
    f.write('Frame %s' % (i + 1) + '\n')
    f.write('{:^6}{:^8}{:^8}{:^8}'.format('Atom', 'x', 'y', 'z') + '\n')
    for j in range(atoms):
        f.write('{:6s}{:8.3f}{:8.3f}{:8.3f}'.format(id[0, j], pos[0, j, i], pos[1, j, i], pos[2, j, i]) + '\n')
    f.write('{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}'.format(box[0, i], box[1, i], box[2, i], box[3, i]
                                                                            , box[4, i], box[5, i], box[6, i]
                                                                            , box[7, i], box[8, i]) + '\n')

f.close()
exit()

# pos = Get_Positions.get_positions('wiggle_traj.gro', 'NA', 'HII', 'no')[0]
# f = open('pos_NA_612', 'w')
# np.save(f, pos)
# f.close()
pos = np.load('pos_NA_612')
print 'positions got'
nT_tot = np.shape(pos)[2] - (np.shape(pos)[2] % 100)
print nT_tot
nT = np.linspace(0, nT_tot, 51, dtype=int)
box = np.load('box_array')

images = 10
n = 0
for frame in nT:

    pt_periodic = Periodic_Images.pbcs(pos, images, 60, box[0, frame], box[0, frame], frame)

    pts = np.shape(pt_periodic)[2]
    duplicates = np.shape(pt_periodic)[1]

    f = open('saxs_frames_far/frame_%s.txt' % n, 'w')
    for j in range(duplicates):
        for i in range(pts):
            row = str(pt_periodic[:, j, i])
            row = row.replace("[", "")
            row = row.replace("]", "")
            f.write(row + "\n")

    f.close()

    n += 1

    # f = open('saxs_frames/frame_%s.gro' % frame, 'w')
    #
    # f.write("This is a .gro file\n")
    # f.write("%s\n" % (pts*duplicates))
    # count = 1
    # for j in range(duplicates):
    #     for i in range(pts):
    #         row = str(pt_periodic[:, j, i])
    #         f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:9.3f}'.format(count, 'NA', 'NA', count, pt_periodic[0, j, i],
    #                                                                     pt_periodic[1, j, i], pt_periodic[2, j, i]) + "\n")
    #         count += 1
    #
    # f.write('   0.00000   0.00000  0.00000')
    # f.close()

    print 'frame %s' % n

exit()
NA_list = Get_Positions.get_positions('wiggle.gro', 'NA', 'HII', 'no')[0]
print np.shape(NA_list)

exit()
# pos = np.load('NApos_array612ns')
D = 9.122 * 10 ** -14
D_std = 9.77 * 10 ** -16
q = 1.609 * 10 ** - 19
kb = 1.38 * 10 ** -23
sigma = 1.3 * 10 ** - 5

Conc_props = np.load('C_props_array612ns')
C = Conc_props[0]
C_std = Conc_props[1]
Lc = Conc_props[2]
conv = Conc_props[3]
z_2 = Conc_props[4]
z_1 = Conc_props[5]

print sigma * np.sqrt((C_std / C) ** 2 + ( D_std / D ) ** 2)

exit()
x = subprocess.check_output(["python", "Concentration.py", "-i", "wiggle_traj7.gro", "-c", "NA", "-b", "0.1", "-l", "HII", "-s" "no"]).splitlines()
x = [float(i) for i in x]
print x
exit()
dt = np.load('dt')
print dt
top = 1.5
pt = 1.25
bot = 1
sigma = 0.25
x = np.array([1, 2, 3, 4])
shift = np.array([1, 2, 3, 4])
randos = np.random.randn(1, 4)
print shift
print randos, x
print randos*x + shift
exit()
dq_full = np.load('dq2')
# dq_fifth = np.load('dq_fifth')
nT = len(dq_full)
seg = 30
start = 0
pts = 500
dt = 1


def ic(dt, start, pts):
    dq_fifth = dq_full[start:pts]
    # nT = len(dq_full)
    # dq_cum_full = np.zeros([nT])
    # dq_cum_full[0] = dq_full[0]
    # for i in range(1, nT):
    #     dq_cum_full[i] = dq_full[i] + dq_cum_full[i - 1]
    #
    # msd_full = np.zeros([nT])
    # itau = 0  # counts up to the number of trajectory frames
    # while itau < nT:
    #     ncount = (nT-itau)
    #     for t in range(nT-itau):  # run for total number of time points from itau to nT
    #         xo = dq_cum_full[t + itau] - dq_cum_full[t]
    #         msd_full[itau] += np.dot(xo, xo)  # Should there be a division by 't' in here?
    #     msd_full[itau] /= ncount
    #     itau += 1
    #
    times = np.linspace(dt, (pts - start)*dt, (pts - start))
    # plt.plot(times, msd_full[0:(pts - start)])

    nT = len(dq_fifth)
    dq_cum = np.zeros([nT])
    dq_cum[0] = dq_fifth[0]
    for i in range(1, nT):
        dq_cum[i] = dq_fifth[i] + dq_cum[i - 1]

    msd = np.zeros([nT])
    itau = 0  # counts up to the number of trajectory frames
    while itau < nT:
        ncount = (nT-itau)
        for t in range(nT-itau):  # run for total number of time points from itau to nT
            xo = dq_cum[t + itau] - dq_cum[t]
            msd[itau] += np.dot(xo, xo)  # Should there be a division by 't' in here?
        msd[itau] /= ncount
        itau += 1

    fracshow = .6
    T = 300
    kb = 1.38*10**-23
    y_fit, r_squared, std2, coeff = Poly_fit.poly_fit(times[1:fracshow*(pts - start)], msd[1:fracshow*(pts - start)], 1)
    conv = 157899298.191  #
    ic_cd = (coeff[1]*(1*10**12)*conv)/(2*kb*T)
    return ic_cd

# segs = np.linspace(100, 700, 13)
# for seg in segs:
#     count = 1
#     x = np.arange(nT - seg)
#     ics = []
#     for i in x:
#         ics.append(ic(dt, i, seg + i))
#     plt.figure(count)
#     plt.plot(x, ics)
#     plt.title('Seg = %s' % seg)
#     plt.show(block=False)
#     count += 1

x = np.arange(nT - seg)
ics = []
for i in x:
    ics.append(ic(dt, i, seg + i))

print np.mean(ics)
print np.std(ics)
plt.title('Seg = %s' % seg)
plt.plot(x, ics)
plt.show()

# print 'Ionic Conductivity: %s' % ic_cd
# plt.plot(times, msd)
# plt.plot(times[1:fracshow*(pts - start)], y_fit)
# plt.show()
