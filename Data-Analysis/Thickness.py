# Find the membrane thickness based on a maximum and minimum z position
# Also gives the amount of space needed for water molecules given the amount of space wanted between periodic boundaries
# in the z direction

water_layer = 6  # nm of water wanted between layers

f = open("new_box.gro", "r")  # .gro file whose positions of Na ions will be read
a = []  # list to hold lines of file
for line in f:
    a.append(line)

z = []  # list to hold z positions of all atoms

for i in range(2, len(a) - 1):
    z.append(float(a[i][36:44]))

z_max = max(z)
z_min = min(z)

thickness = z_max - z_min

print 'The membrane is %s nm thick' %thickness
print ''
print 'The z box vector should be %s nm' %(thickness + water_layer)