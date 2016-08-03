# Look at output from gmx energy and plots the box vector over the course of the trajectory

import matplotlib.pyplot as plt
import math

fx = open("4.4nm_box_x.xvg")
a = []
for line in fx:
    a.append(line)
fy = open("4.4nm_box_y.xvg")
b = []
for line in fy:
    b.append(line
             )
time_x = []  # time (ps)
box_x = []  # energy (kJ/mol)

for i in range(23, len(a)):
    time_x.append(float(a[i][0:12]))
    box_x.append(float(a[i][14:28]))

time_y = []
box_y = []

for i in range(23, len(b)):
    time_y.append(float(b[i][0:12]))
    box_y.append(float(b[i][14:28])/math.sin(math.pi/3))

plt.plot(time_x, box_x)
plt.plot(time_y, box_y)
plt.ylabel('Box Vector (nm)')
plt.xlabel('Time (ps)')
plt.title('GROMACS Energies')
plt.show()


