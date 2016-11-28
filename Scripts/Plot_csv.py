#!/usr/bin/python

# Plot a csv file (created to plot azimuthal distribution data from Yale initially)

import argparse
import numpy as np
import pylab as plb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import math
import scipy.optimize as sci

parser = argparse.ArgumentParser(description = 'Run Cylindricity script')

parser.add_argument('-f', '--file', default='1T.csv', help = 'Path to input file')
parser.add_argument('-t', '--toplines', default=2, help='Number of lines before actual data begins')

args = parser.parse_args()

f = open('%s' % args.file, 'r')

a = []
for line in f:
    a.append(line)
f.close()

toplines = int(args.toplines)

# find number of fields

if str(args.file).endswith('.csv'):
    commas = 0
    comma_index = []
    extension = '.csv'
    for i in range(len(a[toplines])):
        if a[toplines][i] == ',':
            comma_index.append(i)
            commas += 1


    fields = commas + 1

    pts = len(a) - toplines
    data = np.zeros([fields, pts])

    for i in range(toplines, len(a)):
        for j in range(fields):
            data_pt = a[i][j*15:(j + 1)*14 + j]
            if data_pt.count("+") != 0:
                data_pt = data_pt.replace("+", "")
            # elif data_pt.count('-') != 0:
            #     data_pt = data_pt.replace('-', '')
            data[j, i - toplines] = float(data_pt)

elif str(args.file).endswith('.txt'):

    fields = 3
    pts = len(a) - toplines
    data = np.zeros([fields, pts])
    extension = '.txt'

    for i in range(toplines, len(a)):
        for j in range(fields):
            data[j, i - toplines] = a[i][11*j:(j + 1)*10 + j]


# Fit the data to a gaussian distribution
xrs = 0  # x_range_shift
y_shift = min(data[1, xrs:pts/2 - xrs]) / max(data[1, :pts/2])
x = ar(data[0, (pts/4 - xrs):(pts/4 + xrs)])# / 100
y = ar(data[1, (pts/4 - xrs):(pts/4 + xrs)]) / max(data[1, :pts/2])

n = len(x)
mean = sum(x*y)/n
sigma = sum(y*(x-mean)**2)/n


def gaus(x, a, x0, sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2)) + y_shift

popt, pcov = curve_fit(gaus, x, y, p0=[1, mean, sigma])

plt.plot(x, y, 'b+:', label='data')
plt.plot(x, gaus(x, *popt), 'ro:', label='fit')
plt.legend()
plt.suptitle('Azimuthal Distribution, %s Magnetic Field' % args.file.replace(extension, ''))
plt.title('a = %.6f $\mu$ = %.2f, $\sigma$ = %.2f' % (popt[0], popt[1], abs(popt[2])))
plt.ylabel('Intensity')
plt.xlabel('Angle (degrees)')
plt.show()
