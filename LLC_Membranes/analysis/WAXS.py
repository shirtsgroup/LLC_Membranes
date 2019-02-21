#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.patheffects as PathEffects
from matplotlib import ticker


def normalize_alkanes(R, Z, Raw_Intensity, inner, outer, angle, nbins=45):
    """
    Plot angular integration of 2D WAXS data bounded by a circle defined by radii 'inner' and 'outer'
    :param R: points in r direction
    :param Z: points in z direction
    :param Raw_Intensity: values at all (R, Z) points on grid
    :param inner: inside radius of region bounding alkane reflections
    :param outer: outside radius of region bounding alkane reflections
    :return: Intensity values normalized by average intensity inside alkane region
    """

    bins = np.linspace(-90, 90, nbins)

    angles = []
    intensity = []

    count = 0
    test = np.zeros_like(Raw_Intensity)
    for i in range(R.shape[0]):
        for j in range(Z.shape[0]):
            if inner < np.linalg.norm([R[i], Z[j]]) < outer:
                angles.append((180 / np.pi) * np.arctan(Z[j] / R[i]))
                count += 1
                intensity.append(Raw_Intensity[i, j])
                if -60 < angles[count - 1] < 60:
                    test[i, j] = Raw_Intensity[i, j]

    # fancy demo of r-alkanes normalization
    plt.figure()
    test /= np.amax(Raw_Intensity)
    factor = 3.1
    test *= factor
    colorbar = 'jet'
    levels = np.linspace(0, factor, 200)
    heatmap = plt.contourf(X, Y, test.T, cmap=colorbar, levels=levels, extend='max')
    cbar = plt.colorbar(heatmap, format='%.1f')
    cbar.ax.tick_params(labelsize=14)
    cbar.set_ticks([0, 0.5, 1, 1.5, 2, 2.5, 3.1])
    plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=18)
    plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=18)
    plt.gcf().get_axes()[0].tick_params(labelsize=14)
    plt.gcf().get_axes()[0].set_aspect('equal')
    plt.tight_layout()
    #plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/ralkanes.png')

    inds = np.digitize(angles, bins)

    I = np.zeros([nbins])
    counts = np.zeros([nbins])
    for i in range(len(inds)):
        I[inds[i] - 1] += intensity[i]
        counts[inds[i] - 1] += 1

    # Get average intensity in ring excluding 60 degree slice around top and bottom #######

    bin_range = 180 / nbins  # degrees which a single bin covers

    start = int((angle / 2) / bin_range)  # start at the bin which covers -60 degrees and above
    end = nbins - start  # because of symmetry

    total_intensity = np.sum(I[start:end])
    avg_intensity = total_intensity / np.sum(counts[start:end])

    print('Average Intensity in alkane chain region : %s' % avg_intensity)

    # Plot R-spots integrated intensity curve
    plt.figure()
    I /= (counts*avg_intensity)
    bins -= 90
    plt.plot(np.abs(bins[::-1]), I[::-1], linewidth=2)
    plt.xlabel('Angle with respect to $q_r=0$', fontsize=18)
    plt.ylabel('Normalized integrated intensity', fontsize=18)
    plt.gcf().get_axes()[0].tick_params(labelsize=18)
    plt.tight_layout()
    #plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/angular_integration.pdf')

    return avg_intensity


def onclick(event):

    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))

    # find which index in waxs
    xd = np.argmin(np.abs(X - event.xdata))
    yd = np.argmin(np.abs(Y - event.ydata))

    print(waxs[yd, xd])  # the transpose is plotted


def gaussian(points, mean, sigma, amplitude, yshift):

    return yshift + (amplitude / np.sqrt(2*np.pi*sigma**2)) * np.exp(-(points - mean)**2/(2*sigma**2))


def lorentz(points, a, b, c, yshift):
    """
    :param p: lorentzian parameters : [full width half max (FWHM), position of maximum, maximum heigth]
    :param p: position
    :return:
    """

    w = a / 2

    x = (b - points) / w

    return yshift + (c / (np.pi*w)) / (1 + x ** 2)


waxs = np.load('WAXS.npy')

center = np.array([530, 473])
hw = 400  # height and width (pixels)
# center the plot
r_high = hw - 155
r_low = hw - 215
bins = 720
db = 180 / float(bins)

waxs = waxs[center[0]-hw:center[0]+hw, center[1]-hw:center[1]+hw]  # crop image

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

X = np.linspace(-qmax, qmax, waxs.shape[0])
Y = np.linspace(-qmax, qmax, waxs.shape[1])


avg_intensity = normalize_alkanes(X, Y, waxs, 1.256, 1.57, 120)

print('New Average Intensity in alkane chain region : %s' % avg_intensity)

factor = 3.1
colorbar = 'jet'
levels = np.linspace(0, factor, 200)
waxs /= avg_intensity

################# R-pi / R-double intensity measurement illustration ##################
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(X, waxs[:, waxs.shape[0]//2 + 1], linewidth=2)
plt.gcf().get_axes()[0].tick_params(labelsize=18)

ax.annotate('R-$\pi$', xy=(-1.7, 2.5), xytext=(-1.25, 2.25), arrowprops=dict(facecolor='black', width=1, headwidth=6,
                                                                             shrink=0.2, headlength=6), fontsize=14)
ax.annotate('R-double', xy=(-.81, .906), xytext=(-.3, 0.72), arrowprops=dict(facecolor='black', width=1, headwidth=6,
                                                                             shrink=0.15, headlength=6), fontsize=14)

plt.xlabel('$q_z\ (\AA^{-1})$', fontsize=18)
plt.ylabel('Normalized Intensity', fontsize=18)
plt.tight_layout()
#plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/rpi_rdouble.pdf')

################ Plot and fit lorentzian and gaussian functions to qr and qz experimental Cross-sections ##############

plt.figure()
rsection = waxs[np.argmax(waxs[:, waxs.shape[0]//2])]

p = np.array([0, 0.3, 2.5, 0])  # initial guess at fit parameters
solp_gaussian, cov_x_gaussian = curve_fit(gaussian, X, rsection, p)
p = np.array([0.1, 0, 1, .1])
solp_lorentz, cov_x_lorentz = curve_fit(lorentz, X, rsection, p)

plt.plot(X, rsection, linewidth=2, color='xkcd:blue')
#plt.plot(X, gaussian(X, solp_gaussian[0], solp_gaussian[1], solp_gaussian[2], solp_gaussian[3]), '--', label='Gaussian Fit', linewidth=2)
#plt.plot(X, lorentz(X, solp_lorentz[0], solp_lorentz[1], solp_lorentz[2], solp_lorentz[3]), '--', color='xkcd:orange', label='Lorentzian Fit', linewidth=2)
plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=18)
plt.ylabel('Intensity', fontsize=18)

print("Lorentzian FWHM = %.3f +/- %.3f A^-1" % (solp_lorentz[0], cov_x_lorentz[0, 0] ** 0.5))
print("Gaussian FWHM = %.3f +/- %.3f A^-1" % (2*np.sqrt(2*np.log(2))*solp_gaussian[1],
                                       2 * np.sqrt(2 * np.log(2)) * cov_x_gaussian[1, 1] ** 0.5))
plt.gcf().get_axes()[0].tick_params(labelsize=18)
plt.legend(fontsize=17)
plt.tight_layout()

plt.figure()
start = np.argmin(np.abs(Y - 1.6))
plot_start = Y.size // 2
zsection = waxs[start:, waxs.shape[0]//2]

p = [1.7, 0.3, 2.5, 0]
solp_gaussian, cov_x_gaussian = curve_fit(gaussian, Y[start:], zsection, p)
p = np.array([0.1, 1.7, 1, 0.1])
solp_lorentz, cov_x_lorentz = curve_fit(lorentz, Y[start:], zsection, p)
plt.plot(Y[plot_start:], waxs[plot_start:, waxs.shape[0]//2], linewidth=2, color='xkcd:blue')
plt.plot(Y[start:], gaussian(Y[start:], solp_gaussian[0], solp_gaussian[1], solp_gaussian[2], solp_gaussian[3]), '--', label='Gaussian Fit', linewidth=2)
plt.plot(Y[start:], lorentz(Y[start:], solp_lorentz[0], solp_lorentz[1], solp_lorentz[2], solp_lorentz[3]), '--', color='xkcd:orange', label='Lorentzian Fit', linewidth=2)


print("Lorentzian FWHM = %.3f +/- %.3f A^-1" % (solp_lorentz[0], cov_x_lorentz[0, 0] ** 0.5))
print("Gaussian FWHM = %.3f +/- %.3f A^-1" % (2*np.sqrt(2*np.log(2))*solp_gaussian[1],
                                       2 * np.sqrt(2 * np.log(2)) * cov_x_gaussian[1, 1] ** 0.5))

plt.xlabel('$q_z\ (\AA^{-1}$)', fontsize=18)
plt.ylabel('Intensity', fontsize=18)
plt.gcf().get_axes()[0].tick_params(labelsize=18)
plt.legend(fontsize=18)
plt.tight_layout()

################ Calculate R-pi and R-double from cross-section ########################

print('Average R-pi intensity: %.2f' % np.mean(np.amax(waxs[:, waxs.shape[0]//2 - 10:waxs.shape[0]//2 + 10], axis=0)))

R_double_top = 0
R_double_bottom = 0
while Y[R_double_bottom] < -1.0:
    R_double_bottom += 1

while Y[R_double_top] < 0:  # The bottom part of the plot is noisier (change 0 to 1.0 to see what I mean)
    R_double_top += 1

print('Average R-double intensity: %.2f' % np.mean(np.amax(waxs[R_double_bottom:R_double_top,
                                                           waxs.shape[0]//2 - 10:waxs.shape[0]//2 + 10], axis=0)))


################# Plot normalized WAXS pattern with colorbar #################

fig, ax = plt.subplots()
heatmap = plt.contourf(X, Y, waxs, cmap=colorbar, levels=levels)
cid = fig.canvas.mpl_connect('button_press_event', onclick)
#cbar = plt.colorbar(format='%.2f')
plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=18)
plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=18)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.gcf().get_axes()[0].set_aspect('equal')
#plt.text(-2.2, 2.2, "(c)", fontsize=24, color='white', verticalalignment='center', horizontalalignment='center')
plt.tight_layout()

###################### Annotated Plot ##########################

fig, ax = plt.subplots()
heatmap = plt.contourf(X, Y, waxs, cmap=colorbar, levels=levels)
plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=18)
plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=18)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.gcf().get_axes()[0].set_aspect('equal')

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


lw = 1.5
txt1.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])
txt2.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])
txt3.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])
txt4.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])
txt5.set_path_effects([PathEffects.withStroke(linewidth=lw, foreground='black')])

cbar = plt.colorbar(heatmap, format='%.1f')
tick_locator = ticker.MaxNLocator(nbins=6)
cbar.ax.tick_params(labelsize=14)
cbar.set_ticks([0, 0.5, 1, 1.5, 2, 2.5, 3.1])

plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=18)
plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=18)
plt.gcf().get_axes()[0].tick_params(labelsize=18)
plt.gcf().get_axes()[0].set_aspect('equal')
plt.xticks([-2, -1, 0, 1, 2])
plt.tight_layout()
#plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/WAXS_annotated_words.png')
plt.show()