#!/usr/bin/env python

import matplotlib.pyplot as plt

with open('raw_SWAXS.txt', 'r') as f:

    a = []
    for line in f:
        a.append(line)

it = 0
while a[it].count('WAXS') == 0:
    it += 1

qW = []
IW = []
it += 1
while a[it] != '\n':
    info = a[it].split()
    qW.append(float(info[0]))
    IW.append(float(info[1]))
    it += 1

while a[it].count('SAXS') == 0:
    it += 1

qS = []
IS = []
it += 1
while a[it] != '\n':
    info = a[it].split()
    qS.append(float(info[0]))
    IS.append(float(info[1]))
    it += 1




plt.figure(1)
plt.plot(qS[:-150], IS[:-150])

plt.figure(2)
plt.plot(qW[:-150], IW[:-150])

plt.show()


