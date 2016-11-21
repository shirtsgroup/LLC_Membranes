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

parser.add_argument('-f', '--csv', default='1T.csv', help = 'Path to input file')
parser.add_argument('-t', '--toplines', default=2, help='Number of lines before actual data begins')

args = parser.parse_args()

f = open('%s' % args.csv, 'r')

a = []
for line in f:
    a.append(line)
f.close()

toplines = int(args.toplines)

# find number of fields
commas = 0
comma_index = []
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

# Fit the data to a gaussian distribution

x = ar(data[0, :pts/2])
y = ar(data[1, :pts/2])/.101512

n = len(x)
mean = sum(x*y)/n
sigma = sum(y*(x-mean)**2)/n


def gaus(x, a, x0, sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

popt, pcov = curve_fit(gaus, x, y, p0=[1, mean, sigma])

plt.plot(x, y, 'b+:', label='data')
plt.plot(x, gaus(x, *popt), 'ro:', label='fit')
plt.legend()
plt.suptitle('Azimuthal Distribution, %s Magnetic Field' % args.csv.replace('.csv', ''))
plt.title('a = %.6f $\mu$ = %.2f, $\sigma$ = %.2f' % (popt[0], popt[1], abs(popt[2])))
plt.ylabel('Intensity')
plt.xlabel('Angle (degrees)')
plt.show()
