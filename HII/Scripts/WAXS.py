#!/usr/bin/env python

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math


waxs = np.load('WAXS.npy')

center = np.array([530, 473])
hw = 400  # height and width (pixels)
r_high = hw - 155
r_low = hw - 215
bins = 720
db = 180 / float(bins)

waxs = waxs[center[0]-hw:center[0]+hw, center[1]-hw:center[1]+hw]

qpi = 1.7  # q value of pi-stacking reflection
axial = np.copy(waxs[:, int(waxs.shape[0]/2)])  # hold x constant at the center, and get all y values. The array is of the shape [y, x]
mid = axial.shape[0] / 2
axial[int(mid - 150):int(mid + 150)] = 0  # zero out middle values so we get the right maximum
Imax = np.amax(axial)  # max value of Intensity. It will correspond to location of pi-stacking reflection
Imax_pixel = np.where(axial == Imax)[0][0]
pixel_to_q = 1.7 / abs(mid - Imax_pixel)
qmax = mid*pixel_to_q
r_high_pix = 155 * pixel_to_q
r_low_pix = 215 * pixel_to_q

# The raw data contains super high intensity specs and some intensity near the center that is meant to be blocked out
# I zero out all of those so that the highest intensity is in the pi-stacking reflection.
for i in range(23):
    m = np.amax(waxs)
    x = np.where(waxs == m)[0][0]
    y = np.where(waxs == m)[1][0]
    waxs[x, y] = 0

# fig, ax = plt.subplots()
# plt.imshow(waxs, cmap='jet')
# plt.show()

####### Isolate to a ring bounded by outer and inner #######
# waxs2 = np.zeros_like(waxs)
#
# outer = 235
# inner = 180
#
# nbins = 45
# bins = np.linspace(-90, 90, nbins)
#
# angles = []
# intensity = 0
# count = 0
# center = waxs2.shape[0] / 2
# for i in range(waxs.shape[0]):
#     for j in range(waxs.shape[1]):
#         if inner < np.linalg.norm([abs(i - center), abs(j - center)]) < outer:
#             intensity += waxs[i, j]
#             waxs2[i, j] = waxs[i, j]
#             count += 1
#             # if (j - center) == 0:
#             #     waxs2[i, j] = waxs[i, j]
#             #     # angles.append((180/np.pi)*np.arctan((i - center)/(j - center)))
#             #     # angles.append(90)
#             #     # intensity.append(waxs[i, j])
#             # elif -60 < (180/np.pi)*np.arctan((i - center)/(j - center)) < 60:
#             #     waxs2[i, j] = waxs[i, j]
#             #     angles.append((180/np.pi)*np.arctan((i - center)/(j - center)))
#             #     intensity.append(waxs[i, j])

#

X = np.linspace(-qmax, qmax, waxs.shape[0])
Y = np.linspace(-qmax, qmax, waxs.shape[1])

inner = 1.1
outer = 1.6

angle = 120
nbins = 45
bins = np.linspace(-90, 90, nbins)

bw = 180 / (nbins - 1)

angles = []
intensity = []
for i in range(waxs.shape[0]):
    for j in range(waxs.shape[1]):
        if inner < np.linalg.norm([X[i], Y[j]]) < outer:
            angles.append((180/np.pi)*np.arctan(Y[j]/X[i]))
            intensity.append(waxs[i, j])

inds = np.digitize(angles, bins)

I = np.zeros([nbins])
counts = np.zeros([nbins])
for i in range(len(inds)):
    I[inds[i] - 1] += intensity[i]
    counts[inds[i] - 1] += 1

#Get average intensity in ring excluding 60 degree slice around top and bottom #######

bin_range = 180 / nbins  # degrees which a single bin covers

start = int((angle/2) / bin_range)  # start at the bin which covers -60 degrees and above
end = nbins - start  # because of symmetry

total_intensity = np.sum(I[start:end])

avg_intensity = total_intensity / np.sum(counts[start:end])

print('New Average Intensity in alkane chain region : %s' % avg_intensity)

# print(np.min(angles))
# inds = np.digitize(angles, bins)
# I = np.zeros([nbins])
# counts = np.zeros([nbins])
# for i in range(len(inds)):
#     I[inds[i]] += intensity[i]
#     counts[inds[i]] += 1

# plt.imshow(waxs2, cmap='jet')
# #################################
# plt.show()

#waxs /= np.amax(waxs)  # normalize with respect to highest intensity in pi-stacking reflection
# waxs /= (intensity/count)

fig = plt.figure()
ax = fig.add_subplot(111)
X = np.linspace(-qmax, qmax, waxs.shape[0])
Y = np.linspace(-qmax, qmax, waxs.shape[1])

factor = 3.0
cbar = 'seismic'
levels = np.linspace(0, factor, 200)
waxs /= avg_intensity

background = 0.5

# restricted = np.zeros_like(waxs)
# for i in range(waxs.shape[0]):
#     for j in range(waxs.shape[1]):
#         if 1.9 < X[i] < 2 and -0.4 < Y[j] < 0.4:
#             # angle = (180/np.pi)*np.arctan(Y[j]/X[i])
#             # if -30 < angle < 30:
#             restricted[i, j] = waxs[i, j]
#
# print(np.mean(restricted))

binarea = (X[1] - X[0]) * (Y[1] - Y[0])
# print('Bin area: %s' % binarea)
# print(np.amax(restricted))
# print(np.count_nonzero(restricted))
# print(np.sum(restricted))
# plt.imshow(restricted)
# plt.show()
# exit()

####################### Move pi-stacking reflections downward ##########################
## Top reflection ##
# qz_lower = 1.5
# qz_upper = 1.85
# qr_lower = -0.75
# qr_upper = 0.75
#
# X_step = X[1] - X[0]
# Y_step = Y[1] - Y[0]
# qz_lower_ndx = len(X) // 2 + int(qz_lower / X_step)
# qz_upper_ndx = len(X) // 2 + int(qz_upper / X_step)
# qr_lower_ndx = len(Y) // 2 - int(qr_upper / Y_step)
# qr_upper_ndx = len(Y) // 2 + int(qr_upper / Y_step)
#
# diff = waxs[qz_lower_ndx:qz_upper_ndx, qr_lower_ndx:qr_upper_ndx] - background
# waxs[qz_lower_ndx:qz_upper_ndx, qr_lower_ndx:qr_upper_ndx] = background
#
# shift = 0.4
# qz_upper -= shift
# qz_lower -= shift
# qz_lower_ndx = len(X) // 2 + int(qz_lower / X_step)
# qz_upper_ndx = len(X) // 2 + int(qz_upper / X_step)
#
# waxs[qz_lower_ndx:qz_upper_ndx, qr_lower_ndx:qr_upper_ndx] += diff
# ## Bottom section ##
# qz_lower = -1.85
# qz_upper = -1.5
#
# X_step = X[1] - X[0]
# Y_step = Y[1] - Y[0]
# qz_lower_ndx = len(X) // 2 + int(qz_lower / X_step)
# qz_upper_ndx = len(X) // 2 + int(qz_upper / X_step)
#
# diff = waxs[qz_lower_ndx:qz_upper_ndx, qr_lower_ndx:qr_upper_ndx] - background
# waxs[qz_lower_ndx:qz_upper_ndx, qr_lower_ndx:qr_upper_ndx] = background
#
# shift = 0.4
# qz_upper += shift
# qz_lower += shift
# qz_lower_ndx = len(X) // 2 + int(qz_lower / X_step)
# qz_upper_ndx = len(X) // 2 + int(qz_upper / X_step)
#
# waxs[qz_lower_ndx:qz_upper_ndx, qr_lower_ndx:qr_upper_ndx] += diff

heatmap = plt.contourf(X, Y, waxs, cmap=cbar, levels=levels, extend='max')
cbar = plt.colorbar(format='%.2f')
plt.savefig('waxs_%s_%.1f.png' % (cbar, factor))
plt.show()
exit()
heatmap = plt.imshow(waxs, cmap='jet', vmax=2.5*(intensity/count), extent=[-qmax, qmax, -qmax, qmax])

# plt.plot(np.zeros([100]), np.linspace(-qmax, qmax, 100), '--', color='black')

# plt.gcf().get_axes()[0].tick_params(labelsize=14)
# plt.xlabel('$q_r\ (\AA^{-1})$', fontsize=14)
# plt.ylabel('$q_z\ (\AA^{-1})$', fontsize=14)
# plt.tight_layout()
# plt.savefig('waxs_dashed.png')
# plt.show()
# exit()
size = 12
txt1 = ax.annotate('R-$\pi$', xy=(0, 1.9), color='w', size=size, weight='bold', ha='center', va='center')
txt2 = ax.annotate('R-helix', xy=(0, -0.6), color='w', size=size, weight='bold', ha='center', va='center')
txt3 = ax.annotate('R-alkanes', xy=(0, -1.3), color='w', size=size, weight='bold', ha='center', va='center')
txt4 = ax.annotate('R-spots', xy=(-1, 1), color='w', size=size, weight='bold', ha='center', va='center')
txt5 = ax.annotate('R-pores', xy=(0, 0.25), color='w', size=size, weight='bold', ha='center', va='center')
txt6 = ax.annotate("", xy=(1, 0.25), xytext=(0.5, 0.25), arrowprops=dict(arrowstyle="simple", connectionstyle="arc3", fc='w'))
txt7 = ax.annotate("", xy=(-1, 0.25), xytext=(-0.5, 0.25), arrowprops=dict(arrowstyle="simple", connectionstyle="arc3", fc='w'))
txt6 = ax.annotate("", xy=(1.3, -0.4), xytext=(0.7, -1.15), arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=0.2", fc='w'))
txt7 = ax.annotate("",  xy=(-1.3, -0.4), xytext=(-0.7, -1.15), arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=-0.2", fc='w'))
import matplotlib.patheffects as PathEffects

lw = 1.5
txt1.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])
txt2.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])
txt3.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])
txt4.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])
txt5.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])

cbar = plt.colorbar(heatmap)
from matplotlib import ticker
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.set_yticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])

plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=14)
plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
plt.show()
np.save('waxs.npy', waxs)
exit()
# fig.savefig('raw_WAXS.png')
# fig.clf()
# #plt.show()
# exit()
x = np.linspace(-qmax, qmax, 2*hw)
y = np.linspace(qmax, -qmax, 2*hw)
xx, yy = np.meshgrid(x, y)

plt.pcolormesh(xx, yy, waxs, cmap='jet')
plt.gcf().get_axes()[0].set_ylim(-2.5, 2.5)
plt.gcf().get_axes()[0].set_xlim(-2.5, 2.5)
plt.gcf().get_axes()[0].set_xlabel('$q_r$', fontsize=14)
plt.gcf().get_axes()[0].set_ylabel('$q_z$', fontsize=14)
plt.gcf().get_axes()[0].set_aspect('equal')
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
# plt.show()
# fig.savefig('raw_WAXS.png')
# fig.clf()
# plt.show()
# exit()
angles = np.zeros([bins])
norm = np.zeros([bins])
ring = np.zeros_like(waxs)

for i in range(waxs.shape[0]):
    for j in range(waxs.shape[1]):
        if r_high > np.linalg.norm([i - hw, j - hw]) > r_low:
            ring[i, j] = waxs[i, j]
            v = np.array([hw, hw]) - np.array([i, j])
            if v[1] != 0:
                angle = np.arctan(float(v[0])/float(v[1])) * (180/np.pi)
                # This appears to be arctan(x/y) but the way I'm using the numpy matrix here, the x coordinate actually
                # represents the y component and vice versa. Here is an example: Consider the 3x3 matrix:
                # [ 1 2 3 ]     If we write the indices in    [(0, 0), (0, 1), (0, 2)]
                # [ 4 5 6 ]     place of the entries it       [(1, 0), (1, 1), (1, 2)]
                # [ 7 8 9 ]     looks like:                   [(2, 0), (2, 1), (2, 2)]
                # If I use the coordinates as unit distances, as I do above, I can measure the angle between points as
                # arctan(y/x). If I am interested in specifically the angle between (0,0) and (2,1) w.r.t. the
                # horizontal, then the y displacement is 2 and the x displacement is 1. The angle is arctan(2/1) = 63.
                # notice that 2 corresponds to the first entry of (2,1) which is typically thought of as x and 1
                # corresponds to the second entry which is typically though of as y.
            else:
                angle = 90

            bin = int((bins/2) + ((angle/90)*(bins/2)))
            if bin == bins:
                bin -= 1

            angles[bin] += waxs[i, j]
            norm[bin] += 1


avg = angles / norm  # normalize intensities so it is on a per count basis
avg /= sum(avg)
angles = np.linspace(-90, 90, bins + 1)  # We will only see angles in the range of -90 to 90 since we use np.arctan
bin_angles = [(angles[i] + angles[i + 1])/2 for i in range(bins)]  # bars will be placed in the middle of the bins
width = angles[1] - angles[0]  # width of bins

max1 = np.max(avg[int((90 - 50)/db):int((90 - 25)/db)])
max2 = np.max(avg[int((90 + 25)/db):int((90 + 50)/db)])
print(90 - np.where(avg == max1)[0][0]*db)
print(-90 + np.where(avg == max2)[0][0]*db)
print(max1, max2)

exit()
print(np.max(avg) / np.mean(avg))

plt.figure()
plt.bar(bin_angles, avg, width=width)
plt.xlabel('Angle (degrees)', fontsize=14)
plt.ylabel('Normalized Intensity', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('integrated_WAXS_ring.png')
plt.show()

# plt.imshow(ring, vmax=0.05)
# plt.show()
