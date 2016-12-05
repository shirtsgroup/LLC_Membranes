#!/usr/bin/python

# Look at output from gmx energy and plots the box vector over the course of the trajectory

import matplotlib.pyplot as plt
import math
import argparse
import Poly_fit

parser = argparse.ArgumentParser(description = 'Crosslink LLC structure')  # allow input from user

# Flags
parser.add_argument('-x', '--xvg', default='energy.xvg', help = '.xvg file to be read')
parser.add_argument('-f', '--fit', default='off', help = 'Whether to do polynomial fit to the data. (Only linear atm')

args = parser.parse_args()


def extract_data(xvg, start_fit, end_fit):
    f = open(xvg, 'r')
    a = []
    for line in f:
        a.append(line)

    x = []  # Independent Variable
    y = []  # Dependent Variable

    # Find where data starts:
    data_start = 0
    while a[data_start].count('#') != 0:
        data_start += 1

    title = a[data_start][12:(len(a[data_start]) - 2)]
    xlabel = a[data_start + 1][19:(len(a[data_start + 1]) - 2)]
    ylabel = a[data_start + 2][19:(len(a[data_start + 2]) - 2)]

    while a[data_start].count('@') != 0:
        data_start += 1

    while a[data_start].count('#') != 0:  # turn off for efield.xvg
        data_start += 1

    x_end = 0
    while a[data_start][x_end] == ' ':
        x_end += 1

    while a[data_start][x_end] != ' ':
        x_end += 1

    for i in range(data_start, len(a)):
        x.append(float(a[i][0:x_end]))
        if xvg == 'rdf.xvg':
            y.append(float(a[i][x_end + 1:x_end + 10]))
        elif xvg == 'efield.xvg':
            y.append(float(a[i][34:46]))
        else:
            y.append(float(a[i][x_end + 1:len(a[i])]))
    print len(y)
    print y.index(max(y))

    if args.fit == 'on':
        dt = x[1] - x[0]
        print dt
        endMSD = int(np.floor(len(y)*end_fit))
        startMSD = int(np.floor(len(y)*start_fit))
        print endMSD, startMSD
        y_fit, r_squared, std, coeff = Poly_fit.poly_fit(dt*np.array(range(startMSD, endMSD)), y[startMSD:endMSD], 1)
        return y_fit, coeff

    return x, y, title, xlabel, ylabel

if __name__ == '__main__':
    import numpy as np
    start_fit = 0.25
    end_fit = 0.1
    if args.fit == 'on':
        x, y, title, xlabel, ylabel, y_fit, coeff = extract_data(args.xvg, start_fit, end_fit)
        plt.plot(x[int(np.floor(start_fit*len(x))):int(np.floor(end_fit*len(x)))], y_fit)
        print 'Diffusion Coefficient: %s' % (coeff[1]/(2*0.1))
    else:
        x, y, title, xlabel, ylabel = extract_data(args.xvg, start_fit, end_fit)
    # print 'Mean Value: %s' % np.mean(y[50:len(y) - 1])
    print 'Mean Value: %s' % np.mean(y[12:50])
    print 'Standard Deviation: %s' % np.std(y[50:len(y) - 1])
    maximum = max(y)
    minimum = min(y[int(len(y)/5.0):len(y)])
    print 'V_drop = %s' % (maximum - minimum)
    print len(x)
    print len(y)
    plt.plot(x, y)
    plt.ylabel('%s' % ylabel, fontsize=22)
    plt.xlabel('%s' % xlabel, fontsize=22)
    plt.xticks(size=18)
    plt.yticks(size=18)
    plt.title('%s' % title, fontsize=22)
    plt.show()


