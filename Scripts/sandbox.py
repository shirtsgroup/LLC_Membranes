#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# First set up the figure, the axis, and the plot element we want to animate



exit()
fig = plt.figure()
ax = plt.axes(xlim=(0, 2), ylim=(-2, 2))
line, = ax.plot([], [], lw=2)
line2, = ax.plot([], [], lw=2)
line3, = ax.plot([], [], lw=5)
line4, = ax.plot([], [], lw=2)
line5, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    line5.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(0, 2, 1000)
    y = np.sin(2 * np.pi / 1.54 * (x - 0.01 * i))
    y2 = np.sin(2 * np.pi / 1.54 * (x - 0.01 * i - 0.5))
    y3 = np.sin(2 * np.pi / 1.54 * (x - 0.01 * i - 1))
    y_steady = np.sin(2*np.pi / 1.54 * (x - 0.01 * i - 1.25))
    ysum = y + y_steady
    line.set_data(x, y)
    line2.set_data(x, y_steady)
    line3.set_data(x, ysum)
    line4.set_data(x, y2)
    line5.set_data(x, y3)
    return line, line2, line3, line4, line5

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
           frames=200, interval=20, blit=True)

plt.show()
exit()

bins = np.load('angle_bins.npy')
angles = np.linspace(-90, 90, len(bins))

# FT = np.fft.fft(n)

# FFT = scipy.fftpack.rfft(n, n=1000)
# print FFT[:10]
#
# a0 = FFT[0]
# an = [FFT[a] for a in range(len(FFT)) if a % 2 != 0]
# bn = [FFT[b] for b in range(len(FFT)) if b % 2 == 0 and b != 0]

# a0 = FT.real[0]
# an = FT.real[1:90]
# bn = FT.imag[1:90]
#
#
# # def fourier(a0, a, b, N, nterms):
# def fourier(A, N, nterms):
#
#     """
#     :param x: transformed values
#     :param a0: first fourier coefficient
#     :param an: real fourier coefficient
#     :param bn: imaginary part
#     :param N: number of data points
#     :return: f(x) as approximated by a fourier series
#     """
#
#     # f = 0.5*a0*np.ones([N - 1])
#     f = np.zeros([N - 1])
#     x = np.linspace(-90, 90, N - 1)
#
#     for i in range(len(x)):
#         for n in range(1, nterms):
#             # f[i] += a[j]*np.cos(j*np.pi*x[i]/90) + b[j]*np.sin(j*np.pi*x[i]/90)
#             f[i] += A[n]*np.exp(1j*2*np.pi*n*x[i]/90)
#
#     return f, x
#
# #f, x = fourier(a0, an, bn, len(n), 8)
# f, x = fourier(FT, len(n), 8)
# # f /= sum(f)

period = 180

def cn(n):
    c = bins*np.exp(-1j*2*n*np.pi*angles/period)
    return c.sum()/c.size


def f(x, Nh):
    f = np.array([2*cn(i)*np.exp(1j*2*i*np.pi*x/period) for i in range(1, Nh + 1)])
    return f.sum()

fapprox = np.array([f(t, 80).real for t in angles]) + np.mean(bins)

plt.plot(angles, bins)
plt.plot(angles, fapprox)
plt.show()