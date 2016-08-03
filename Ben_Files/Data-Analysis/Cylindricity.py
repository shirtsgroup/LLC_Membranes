# A script to measure the cylindricity of membrane pores based on their standard deviation from some central axis
# Also measures distance between pores
import numpy as np
import math

# Pore Numbering : Each corner is the location of the center of the pore
#  Pore 2 ---------- Pore 1
#         \         \
#          \         \
#           \         \
#     Pore 4 ---------- Pore 3

f = open("Monomer1.gro", "r")  # .gro file whose positions of Na ions will be read
a = []  # list to hold lines of file
for line in f:
    a.append(line)

no_atoms = 137  # number of atoms in a single monomer
no_layers = 20  # number of layers in the membrane structure
mon_per_layer = 6  # number of monomers in each layer
no_pores = 4  # number of pores
# no_monomers =  # in case there are random distributions of monomers in each layer
no_ion = no_layers*mon_per_layer*no_pores
no_monomers = no_ion  # number of monomers (= number of ions since there is one monomer per ion in the case of sodium)

x = []  # list to hold x positions of ions
y = []  # list to hold y positions of ions
z = []  # list to hold z positions of ions. Z axis runs parallel to pore
start = no_atoms*no_layers*mon_per_layer*no_pores + 2  # This is the line number where the sodium ions will be listed
for i in range(start, start + no_ion):  # adds all sodium coordinates to respective x, y and z lists
    x.append(float(a[i][20:28]))
    y.append(float(a[i][28:36]))
    z.append(float(a[i][36:44]))


# find the central axis in each pore:
x_axis = []
y_axis = []
for k in range(0, no_pores):  # calculates average x and y values in each pore. Taken to be the center of the pore
    sum_x = 0
    sum_y = 0
    for i in range(k * no_ion/4, (k + 1) * no_ion/4):
        sum_x += x[i]
        sum_y += y[i]
    x_avg = sum_x / (no_layers*mon_per_layer)
    y_avg = sum_y / (no_layers*mon_per_layer)
    x_axis.append(x_avg)  # average x coordinates [x_coord_pore1, x_coord_pore2, x_coord_pore3, x_coord_pore4]
    y_axis.append(y_avg)  # average y coordinates [y_coord_pore1, y_coord_pore2, y_coord_pore3, y_coord_pore4]

# find distance between pores
pore12 = math.sqrt((x_axis[0] - x_axis[1])**2 + (y_axis[0] - y_axis[1])**2)
pore13 = math.sqrt((x_axis[0] - x_axis[2])**2 + (y_axis[0] - y_axis[2])**2)
pore34 = math.sqrt((x_axis[2] - x_axis[3])**2 + (y_axis[2] - y_axis[3])**2)
pore42 = math.sqrt((x_axis[1] - x_axis[3])**2 + (y_axis[1] - y_axis[3])**2)
avg_dist_bw_pores = 0.25*(pore12 + pore13 + pore34 + pore42)
pore14 = math.sqrt((x_axis[0] - x_axis[3])**2 + (y_axis[0] - y_axis[3])**2)
pore23 = math.sqrt((x_axis[1] - x_axis[2])**2 + (y_axis[1] - y_axis[2])**2)


# Find the deviation in distance from central axis
deviations = []  # list to hold the standard deviation in distance from center of pore
avg = []  # hold the average distance of ions from the central axis in each pore
for k in range(0, no_pores):
    for j in range(0, no_layers):
        dist_list = []  # list to hold distances of ions from center of pore. Reset on each loop
        for i in range(0, mon_per_layer):  # loops through each pore individually
            dist = math.sqrt((x[i + mon_per_layer*j + mon_per_layer*no_layers*k] - x_axis[k])**2 +
                             (y[i + mon_per_layer*j + mon_per_layer*no_layers*k] - y_axis[k])**2)
            #print k, x_axis[k], y_axis[k], x[i + mon_per_layer*j + mon_per_layer*no_layers*k], y[i + mon_per_layer*j + mon_per_layer*no_layers*k], dist
            # magnitude of vector pointing from center to point where ion is stationed
            dist_list.append(dist)  # add the distance from the center to a list
    avg_dist = np.mean(dist_list)
    avg.append(avg_dist)
    deviations.append(np.std(dist_list))  # add standard deviation for pore to list [pore1, pore2, pore3, pore4]

Average_Deviation = np.mean(deviations)

# Lots of output!

print 'The standard deviation in distance of ions from the center of the pore are as follows:'
print 'Pore 1: %s angstroms' %deviations[0]
print 'Pore 2: %s angstroms' %deviations[1]
print 'Pore 3: %s angstroms' %deviations[2]
print 'Pore 4: %s angstroms' %deviations[3]
print ''
print 'The average deviation from the pore center is %s angstroms' %Average_Deviation
print ''
print 'The distance between each pores is as follows:'
print 'Distance from pore 1 to pore 2: %s angstroms' %pore12
print 'Distance from pore 1 to pore 3: %s angstroms' %pore13
print 'Distance from pore 3 to pore 4: %s angstroms' %pore34
print 'Distance from pore 4 to pore 2: %s angstroms' %pore42
print 'Distance from pore 4 to pore 1: %s angstroms' %pore14
print 'Distance from pore 2 to pore 3: %s angstroms' %pore23
print ''
print 'The average pore to pore distance (not including diagonals) is: %s angstroms' %avg_dist_bw_pores
print 'Pore Orientation for reference: '
print ' Pore 2 ------------- Pore 1'
print '        \            \  '
print '         \            \  '
print '          \            \  '
print '           \            \  '
print '     Pore 4 ------------- Pore 3'





