#!/usr/bin/env python

import numpy as np

with open('bcc.gro') as f:
	a = []
	for line in f:
		a.append(line)

b = []
b.append(a[0])
b.append(a[1])
count = 2
mon = 1
for i in range(2, 79155):
	b.append('{:>5}'.format(mon) + a[i][5:])
	if (i - 1) % 136 == 0 and (i - 2) != 0:
		mon += 1

mon += 1
for i in range(79155, len(a) - 1):
	b.append('{:>5}'.format(mon) + a[i][5:])
	#a[i] = a[i].replace(a[i][:5], '{:>5}'.format(mon))
	mon += 1

with open('test.gro', 'w') as f:
	for line in b:
		f.write(line)
