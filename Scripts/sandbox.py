#!/home/ben/anaconda2/bin/python2

import time
import numpy as np
import matplotlib.pyplot as plt

a = np.array([6065, 5573, 5091, 4622, 4220, 3810, 3447, 3162, 2900, 2698, 2522, 2386, 2268, 2166, 2082, 2003, 1939,
              1864, 1787, 1719, 1659, 1590, 1523, 1453, 1393, 1345, 1289, 1232, 1184, 1125, 1069, 1022, 955, 895, 851,
              791, 742, 677, 632, 571, 535, 475, 426, 370, 322, 251, 199, 147, 99, 54, 0])

# slope = np.zeros([len(a) - 1])
# for i in range(slope.shape[0]):
#     slope[i] = (a[i + 1] - a[i])/ 0.01


x = np.linspace(0, 0.5, 51)
# xp = np.linspace(0, 0.49, 50)

plt.plot(x, a)
plt.title('Number of leftover waters vs. Buffer')
plt.xlabel('Buffer')
plt.ylabel('Leftover Waters')
# plt.plot(xp, slope)
plt.show()


