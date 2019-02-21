#!/usr/bin/env python

import sqlite3 as sql
import matplotlib.pyplot as plt

x = [0.840116, 0.747754, 0.727045, 0.730229, 0.65957, 0.674671, 0.686392, 0.694758, 0.700637, 0.704781, 0.708275, 0.710127]

y = [1281, 1050, 947, 964, 494, 568, 602, 612, 622, 620, 643, 641]

plt.scatter(sorted(x), sorted(y))
plt.show()
