#!/usr/bin/python

# Find the membrane thickness based on a maximum and minimum z position
# Also gives the amount of space needed for water molecules given the amount of space wanted between periodic boundaries
# in the z direction

import argparse

# Input arguments to python file by choosing which monomer to assemble with and number of layers. Can be used in a shell
# script

parser = argparse.ArgumentParser(description = 'Measure Membrane Thickness')
parser.add_argument('-i', '--input', help = 'Path to input file')
parser.add_argument('-w', '--water', default=6, help = 'nm of water between layers')
parser.add_argument('-m', '--no_monomers', default=6, help= 'number of monomers per layer')
parser.add_argument('-l', '--layers', default=20, help= 'Number of stacked layers in unit cell')
parser.add_argument('-p', '--pores', default=4, help= 'Number of pores')
parser.add_argument('-a', '--atoms', default=137, help='Number of atoms excluding ions')
args = parser.parse_args()

water_layer = int(args.water)  # nm of water wanted between layers

f = open(args.input, "r")  # .gro file whose positions of Na ions will be read
a = []  # list to hold lines of file
for line in f:
    a.append(line)

line = 0
while a[line].count('HII') == 0:
    line += 1

z = []  # list to hold z positions of all atoms

benz_carbs = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']

for i in range(line, len(a) -1):
    if str.strip(a[i][11:15]) in benz_carbs:
        z.append(float(a[i][36:44]))
# for i in range(0, tot_monomers):
#     for k in range(0, 5):
#         z.append(float(a[i*args.atoms + k + line][36:44]))


z_max = max(z)
z_min = min(z)

thickness = z_max - z_min
tot_thickness = thickness + water_layer

print 'Membrane Thickness: %s nm' %tot_thickness