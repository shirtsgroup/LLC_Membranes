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

X = np.linspace(-qmax, qmax, waxs.shape[0])
Y = np.linspace(-qmax, qmax, waxs.shape[1])

# background_intensity = waxs[350, 300]
# #
# for i in range(X.size):
#     for j in range(Y.size):
#         if 1.1 < np.linalg.norm([X[i], Y[j]]) < 1.6:
#             waxs[i, j] = background_intensity

# plt.imshow(waxs)
# plt.show()
# exit()

# ######################## Inverse FT ###############################
# X = X[:-1]
# Y = Y[:-1]
# fbin = X[1] - X[0]  # size of bins in fourier space
# system_size = 2*np.pi / fbin  # fourier_bin = 2*pi / system size therefore system size = 2*pi / fourier_bin
# rbin = system_size / X.shape[0]  # real space bin
#
# # reorder lists so they conform to numpy (https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.ifft2.html#numpy.fft.ifft2)
# start = list(X).index(0)
# X_reordered = np.concatenate((X[start:], X[:start]))
# ndx_x = [list(X).index(i) for i in X_reordered]
#
# start = list(Y).index(0)
# Y_reordered = np.concatenate((Y[start:], Y[:start]))
# ndx_y = [list(Y).index(i) for i in Y_reordered]
#
# waxs_reordered = waxs[ndx_x, :]
# waxs_reordered = waxs_reordered[:, ndx_y]
#
# # inverse fourier transform
# inverse_fft = np.fft.ifft2(waxs_reordered)
#
# inverse_fft = inverse_fft[ndx_x, :]
# inverse_fft = inverse_fft[:, ndx_y]
#
# # fourier transform of inversion as a test
# # ft = np.abs(np.fft.fftn(inverse_fft))**2
# # ft = ft[ndx_x, :]
# # ft = ft[:, ndx_y]
# # plt.imshow(ft)
# # plt.show()
#
# inverse_fft = inverse_fft.real / np.amax(inverse_fft.real)
#
# r = np.linspace(-system_size/2, system_size/2, inverse_fft.shape[0])
#
# plt.plot(r[start:], inverse_fft[start:, start])  # z-distribution function
#
#
# def exponential_decay(x, A, L):
#
#     # return A*np.exp(-(x -rbin) / L)
#     return A*np.exp(-x / L)
#
# # find peaks to fit exponential to
# # import detect_peaks
# # mpd = 100
# # peaks = detect_peaks.detect_peaks(r[start:start + 30], mpd=mpd, show=False)  # adjust mpd if all peaks aren't found
#
# peaks = np.array([0, 6, 12, 19, 25, 41]) + start
#
# # if centers[peaks[0]] < spacing / 2:  # sometimes a peak is found where it shouldn't
# #     peaks = peaks[1:]
#
# # plot locations of peaks
# plt.scatter(r[peaks], inverse_fft[peaks, start], marker='+', c='r', s=200, label='Peak locations')
#
#
# L = 100
# A = 1
# # fit decaying exponential to peaks
# p = np.array([A, L])  # initial guess at parameters
# from scipy.optimize import curve_fit
# solp, cov_x = curve_fit(exponential_decay, r[peaks], inverse_fft[peaks, start], p)
# print(solp)
# # plot decaying exponential fit and decaying exponential with L that we were trying to match
# plt.plot(r[start:], exponential_decay(r[start:], solp[0], solp[1]), '--', c='black', label='Least squares fit')
# # plt.plot(centers, exponential_decay(centers, L), '--', c='blue', label='Theoretical')
# plt.xlim(0, 100)
# plt.ylim(-0.1, 0.2)
# print(r[start:])
# ft = np.abs(np.fft.fft(inverse_fft[start:, start]))
# freq_x = np.fft.fftfreq(r[start:].size, d=rbin)
# ndx = np.argsort(freq_x)
# ft = ft[ndx]
# plt.figure()
# plt.plot(freq_x, ft)
# plt.show()
# exit()
#
# bound1 = 0
# bound2 = 0
# while r[bound1] < -15:
#     bound1 += 1
# while r[bound2] < 15:
#     bound2 += 1
#
# # inverse_fft /= np.mean(inverse_fft)
# levels = np.linspace(np.amin(inverse_fft), 0.08*np.amax(inverse_fft), 200)
# plt.contourf(r[bound1:bound2], r[bound1:bound2], inverse_fft[bound1+1:bound2+1, bound1+1:bound2+1], levels=levels, cmap='seismic', extend='max')
# plt.colorbar()
# plt.xlabel('r ($\AA$)')
# plt.ylabel('z ($\AA$)')
# plt.show()
# exit()
# # plt.plot(z, inverse_fft.real[start:, start])
# # plt.show()
# # plt.imshow(inverse_fft.real/np.amax(inverse_fft.real), vmin=np.amin(inverse_fft.real), vmax=0.5*np.amax(inverse_fft.real))
# # plt.colorbar()
# # plt.show()
# # exit()
########################### END Inverse FT ###################################

# inner = 1.1
# outer = 1.6
#
# angle = 120
# nbins = 45
# bins = np.linspace(-90, 90, nbins)
#
# bw = 180 / (nbins - 1)
#
# angles = []
# intensity = []
# for i in range(waxs.shape[0]):
#     for j in range(waxs.shape[1]):
#         if inner < np.linalg.norm([X[i], Y[j]]) < outer:
#             angles.append((180/np.pi)*np.arctan(Y[j]/X[i]))
#             intensity.append(waxs[i, j])
#
# inds = np.digitize(angles, bins)
#
# I = np.zeros([nbins])
# counts = np.zeros([nbins])
# for i in range(len(inds)):
#     I[inds[i] - 1] += intensity[i]
#     counts[inds[i] - 1] += 1
#
# #Get average intensity in ring excluding 60 degree slice around top and bottom #######
#
# bin_range = 180 / nbins  # degrees which a single bin covers
#
# start = int((angle/2) / bin_range)  # start at the bin which covers -60 degrees and above
# end = nbins - start  # because of symmetry
#
# total_intensity = np.sum(I[start:end])
#
# avg_intensity = total_intensity / np.sum(counts[start:end])
# avg_intensity = waxs[-1, -1]
# print(avg_intensity)


def normalize_alkanes(R, Z, Raw_Intensity, inner, outer, angle):
    """
    Plot angular integration of 2D WAXS data bounded by a circle defined by radii 'inner' and 'outer'
    :param R: points in r direction
    :param Z: points in z direction
    :param Raw_Intensity: values at all (R, Z) points on grid
    :param inner: inside radius of region bounding alkane reflections
    :param outer: outside radius of region bounding alkane reflections
    :return: Intensity values normalized by average intensity inside alkane region
    """

    # temporary override
    # inner = 1.15
    # outer = 1.6

    nbins = 45
    bins = np.linspace(-90, 90, nbins)

    bw = 180 / (nbins - 1)

    angles = []
    intensity = []

    count = 0
    test = np.zeros_like(Raw_Intensity)
    for i in range(R.shape[0]):
        for j in range(Z.shape[0]):
            if inner < np.linalg.norm([R[i], Z[j]]) < outer:
                angles.append((180/np.pi)*np.arctan(Z[j]/R[i]))
                count += 1
                intensity.append(Raw_Intensity[i, j])
                if -60 < angles[count - 1] < 60:
                    test[i, j] = Raw_Intensity[i, j]

    # fancy demo of r-alkanes normalization
    # test /= np.amax(Raw_Intensity)
    # test *= 3.1
    #
    # factor = 3.1
    # colorbar = 'jet'
    # levels = np.linspace(0, factor, 200)
    # fig, ax = plt.subplots()
    # heatmap = plt.contourf(X, Y, test.T, cmap=colorbar, levels=levels, extend='max')
    # cbar = plt.colorbar(format='%.2f')
    # plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=18)
    # plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=18)
    # plt.gcf().get_axes()[0].tick_params(labelsize=14)
    # plt.gcf().get_axes()[0].set_aspect('equal')
    # plt.tight_layout()
    # plt.savefig('ralkanes.pdf')
    # plt.show()
    # exit()

    inds = np.digitize(angles, bins)

    I = np.zeros([nbins])
    counts = np.zeros([nbins])
    for i in range(len(inds)):
        I[inds[i] - 1] += intensity[i]
        counts[inds[i] - 1] += 1

    # Get average intensity in ring excluding 60 degree slice around top and bottom #######

    bin_range = 180 / nbins  # degrees which a single bin covers

    start = int((angle/2) / bin_range)  # start at the bin which covers -60 degrees and above
    end = nbins - start  # because of symmetry

    total_intensity = np.sum(I[start:end])
    avg_intensity = total_intensity / np.sum(counts[start:end])

    print('Average Intensity in alkane chain region : %s' % avg_intensity)

    return avg_intensity

avg_intensity = normalize_alkanes(X, Y, waxs, 1.256, 1.57, 120)

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

# # fig = plt.figure()
# ax = fig.add_subplot(111)
X = np.linspace(-qmax, qmax, waxs.shape[0])
Y = np.linspace(-qmax, qmax, waxs.shape[1])
# np.save('qx', X)
# np.save('qy', Y)
# exit()

factor = 3.1
colorbar = 'jet'
levels = np.linspace(0, factor, 200)
waxs /= avg_intensity


levels_log = np.linspace(np.log10(waxs[-1, -1]), np.log10(np.amax(waxs)), 1000)
#levels_log = np.linspace(np.log10(waxs[-1, -1]), np.log10(np.amax(waxs)) - np.log10(avg_intensity), 1000)

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

# plt.figure()
# plt.contourf(X, Y, np.log10(np.ma.masked_where(waxs==0, waxs)), cmap=colorbar, levels=levels_log)
# plt.colorbar()
# plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=14)
# plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=14)
# plt.gcf().get_axes()[0].tick_params(labelsize=14)
# plt.tight_layout()
# plt.gcf().get_axes()[0].set_aspect('equal')
# plt.savefig('WAXS_log10exp.png')
# plt.show()
# exit()


def onclick(event):

    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))

    # find which index in waxs
    xd = np.argmin(np.abs(X - event.xdata))
    yd = np.argmin(np.abs(Y - event.ydata))

    # print(X[xd], event.xdata)
    # print(Y[yd], event.ydata)
    print(waxs[yd, xd])  # it seems that the transpose is plotted


def Rspots(R, Z, waxs, theta=37, theta_sigma=(7, 5), bounds=(1.256, 1.57)):

    spots = np.copy(waxs)
    inner = bounds[0]
    outer = bounds[1]
    I = []

    for i in range(R.shape[0]):
        for j in range(Z.shape[0]):
            if inner < np.linalg.norm([R[i], Z[j]]) < outer:
                angle = (180 / np.pi) * np.arctan(Z[j] / R[i])
                if (theta - theta_sigma[0]) < angle < (theta + theta_sigma[1]) or \
                        (theta - theta_sigma[0]) < (angle - 2*angle) < (theta + theta_sigma[1]):
                    spots[i, j] = 100
                    I.append(waxs[i, j])

    average_intensity = np.mean(I)

    levels = np.linspace(0, 3.1, 200)
    plt.contourf(R, Z, spots, cmap=colorbar, levels=levels, extend='max')
    cbar = plt.colorbar(format='%.2f')
    plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=18)
    plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=18)
    plt.gcf().get_axes()[0].tick_params(labelsize=14)
    plt.gcf().get_axes()[0].set_aspect('equal')
    plt.tight_layout()
    plt.savefig('rspots.pdf', dpi=300)
    plt.show()

    plt.figure()
    plt.hist(I, bins=25)
    plt.title('Average intensity of R-spots: %.2f' % average_intensity)
    plt.show()
    exit()

    return average_intensity

R_double_top = 0
R_double_bottom = 0
while Y[R_double_bottom] < -1.0:
    R_double_bottom += 1

while Y[R_double_top] < 0:  # The bottom part of the plot is noisier (change 0 to 1.0 to see what I mean)
    R_double_top += 1

# # plot z-slices
# for i in range(-10, 10):
#     plt.plot(Y[R_double_bottom:R_double_top], waxs[R_double_bottom:R_double_top, waxs.shape[0]//2 + i])
#
# plt.show()
# exit()

# plot maximum intensity of each z-slice
# plt.plot(np.linspace(-100, 99, 200), np.amax(waxs[:, waxs.shape[0]//2 - 100:waxs.shape[0]//2 + 100], axis=0))
# plt.plot(np.linspace(-10, 9, 20), np.amax(waxs[R_double_bottom:R_double_top,
#                                           waxs.shape[0]//2 - 10:waxs.shape[0]//2 + 10], axis=0))
# plt.show()

# R-pi / R-double intensity measurement illustration
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(X, waxs[:, waxs.shape[0]//2 + 1], linewidth=2)
# plt.gcf().get_axes()[0].tick_params(labelsize=14)
#
# ax.annotate('R-$\pi$', xy=(-1.7, 2.5), xytext=(-1.25, 2.25), arrowprops=dict(facecolor='black', width=1, headwidth=6,
#                                                                              shrink=0.2, headlength=6), fontsize=14)
# ax.annotate('R-double', xy=(-.81, .906), xytext=(-.3, 0.72), arrowprops=dict(facecolor='black', width=1, headwidth=6,
#                                                                              shrink=0.15, headlength=6), fontsize=14)
#
# plt.xlabel('$q_z~(\AA)$', fontsize=14)
# plt.ylabel('Normalized Intensity', fontsize=14)
# plt.tight_layout()
# plt.savefig('rpi_rdouble.pdf')
# plt.show()


# print('Average R-pi intensity: %.2f' % np.mean(np.amax(waxs[:, waxs.shape[0]//2 - 10:waxs.shape[0]//2 + 10], axis=0)))
# print('Average R-double intensity: %.2f' % np.mean(np.amax(waxs[R_double_bottom:R_double_top,
#                                                            waxs.shape[0]//2 - 10:waxs.shape[0]//2 + 10], axis=0)))
# print('Average R-spots intensity : %.2f' % Rspots(X, Y, waxs, theta=50, theta_sigma=(3, 2), bounds=(1.3, 1.45)))
# exit()
fig, ax = plt.subplots()
heatmap = plt.contourf(X, Y, waxs, cmap=colorbar, levels=levels)#, extend='max')
# cid = fig.canvas.mpl_connect('button_press_event', onclick)
# cbar = plt.colorbar(format='%.2f')
# plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=18)
# plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=18)
# plt.gcf().get_axes()[0].tick_params(labelsize=14)
# plt.gcf().get_axes()[0].set_aspect('equal')
# plt.tight_layout()
# plt.savefig('WAXS_raw.pdf')
# plt.show()
# exit()
# heatmap = plt.imshow(waxs, cmap='jet', vmax=2.5*(intensity/count), extent=[-qmax, qmax, -qmax, qmax])

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
txt2 = ax.annotate('R-double', xy=(0, -0.6), color='w', size=size, weight='bold', ha='center', va='center')
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

cbar = plt.colorbar(heatmap, format='%.2f')
from matplotlib import ticker
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.locator = tick_locator
# cbar.update_ticks()
# cbar.ax.set_yticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])

plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=14)
plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.gcf().get_axes()[0].set_aspect('equal')
plt.tight_layout()
plt.savefig('WAXS_annotated_words.png')
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
