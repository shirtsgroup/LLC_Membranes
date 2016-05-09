# Find the membrane thickness based on a maximum and minimum z position
# Also gives the amount of space needed for water molecules given the amount of space wanted between periodic boundaries
# in the z direction

import argparse

# Input arguments to python file by choosing which monomer to assemble with and number of layers. Can be used in a shell
# script

parser = argparse.ArgumentParser(description = 'Measure Membrane Thickness')
parser.add_argument('-i', '--input', help = 'Path to input file')
parser.add_argument('-w', '--water', help = 'nm of water between layers')
args = parser.parse_args()

water_layer = int(args.water)  # nm of water wanted between layers

f = open(args.input, "r")  # .gro file whose positions of Na ions will be read
a = []  # list to hold lines of file
for line in f:
    a.append(line)

z = []  # list to hold z positions of all atoms

for i in range(2, len(a) - 1):
    z.append(float(a[i][36:44]))

z_max = max(z)
z_min = min(z)

thickness = z_max - z_min
tot_thickness = thickness + water_layer

print tot_thickness