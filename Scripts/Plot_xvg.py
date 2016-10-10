#!/usr/bin/python

# Look at output from gmx energy and plots the box vector over the course of the trajectory

import matplotlib.pyplot as plt
import math
import argparse

parser = argparse.ArgumentParser(description = 'Crosslink LLC structure')  # allow input from user

# Flags
parser.add_argument('-x', '--xvg', default='energy.xvg', help = '.xvg file to be read')

args = parser.parse_args()

def extract_data(xvg):
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
        y.append(float(a[i][x_end + 1:len(a[i])]))
        # y.append(float(a[i][34:46])) # for efield.xvg

    print y.index(max(y))
    return x, y, title, xlabel, ylabel

x, y, title, xlabel, ylabel = extract_data(args.xvg)

if __name__ == '__main__':
    import numpy as np
    print 'Mean Value: %s' % np.mean(y[50:len(y) - 1])
    print 'Standard Deviation: %s' % np.std(y[50:len(y) - 1])
    maximum = max(y)
    minimum = min(y[int(len(y)/5.0):len(y)])
    print 'V_drop = %s' %(maximum - minimum)
    plt.plot(x, y)
    plt.ylabel('%s' % ylabel)
    plt.xlabel('%s' % xlabel)
    plt.title('%s' % title)
    plt.show()


