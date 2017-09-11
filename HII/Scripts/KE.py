#!/usr/bin/python

import numpy as np
import math
import Atom_props as atom


f = open('wiggle.gro', 'r')

a = []
for line in f:
	a.append(line)

start = 0
while a[start].count('HII') == 0:
	start += 1

v = []
for i in range(start, len(a) - 1):
	vx = float(a[i][45:52])
	vy = float(a[i][52:60])
	vz = float(a[i][60:len(a)-1])
	v.append(math.sqrt(vx**2 + vy**2 + vz**2))

m = []
for i in range(start, len(a) - 1):
	m.append(atom.mass[str.strip(a[i][9:15])])

KE_list=[]
for i in range(0, len(m)):
	KE_list.append(0.5*m[i]*v[i]**2)

KE = sum(KE_list)
print 'Total Kinetic Energy: %s KJ/mol' %KE

