#!/usr/bin/python

# Look at output from gmx energy and plots the box vector over the course of the trajectory

import matplotlib.pyplot as plt
import math
import argparse

parser = argparse.ArgumentParser(description = 'Crosslink LLC structure')  # allow input from user

# Flags
parser.add_argument('-x', '--xvg', default='energy.xvg', help = '.xvg file to be read')

args = parser.parse_args()

f = open(args.xvg, 'r')
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

for i in range(data_start, len(a)):
    x.append(float(a[i][0:12]))
    y.append(float(a[i][14:28]))

plt.plot(x, y)
plt.ylabel('%s' % ylabel)
plt.xlabel('%s' % xlabel)
plt.title('%s' % title)
plt.show()


