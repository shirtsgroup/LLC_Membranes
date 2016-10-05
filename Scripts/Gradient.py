#!/usr/bin/python

import argparse
import Thickness
import random
import os.path
import subprocess

# Randomly replaces solvent molecules with solute like gmx solvate, but is capable of creating a concentration gradient

parser = argparse.ArgumentParser(description = 'Run Cylindricity script')

parser.add_argument('-i', '--input', default='wiggle.gro', help = 'Path to input file')
parser.add_argument('-p', '--percent', default=25, help = 'Percent solvent to replace with solute')
parser.add_argument('-s', '--solvent', default='Ar', help = 'Solvent which will be swapped with solute')
parser.add_argument('-S', '--Solute', default='NA', help = 'Name of atom replacing solvent molecules')
parser.add_argument('-t', '--top', default='NaPore.top', help = 'Name of topology file to be modified')
parser.add_argument('-d', '--double', default='on', help='Is this a double layer system? (on/off)')

args = parser.parse_args()

f = open(args.input, 'r')  # The input should be a single layered system regardless

a = []

for line in f:
    a.append(line)
f.close()

# Find planes defining top and bottom of membrane
thick, z_max_1, z_min_1 = Thickness.thickness(a)

# Find the z box dimension
z_box = float(a[len(a) - 1][20:30])

# Make the system into a membrane double layer
p = subprocess.Popen(["gmx", "genconf", "-f", "%s" % args.input, "-nbox", "1", "1", "2", "-o",
                      "%s_double.gro" % args.input.replace('.gro', '')])
p.wait()  # wait for this process to finish

f = open('%s_double.gro' % args.input.replace('.gro', ''))  # Now we need to read in coordinates for the double layer

a = []  # redefine a

for line in f:
    a.append(line)
f.close()

# Find the bounds on the other membrane by translating the other bounds by z_box. This works since gmx genconf simply
# stacks replicas of the system in whatever direction you specify translated by the corresponding box coordinate

z_max_2 = z_max_1 + z_box
z_min_2 = z_min_1 + z_box

z_sol = []
Pore_Ions = []  # We will simultaneously create this index group
top_lines = 0
while a[top_lines].count('HII') == 0:
    top_lines += 1

count = -top_lines  # Things get tricky when there are greater than 100000 atoms so I count the number of lines as a way
# to keep track of atom numbers. The top line numbers are subtracted
for i in range(0, len(a)):
    count += 1
    if a[i].count(args.solvent) != 0:
        z = float(a[i][36:44])
        atom_no = count
        if z_min_2 > z > z_max_1:  # check to see if the solvent molecule is between the two membrane layers
            z_sol.append(count)  # if it is, record it
    if a[i].count(args.Solute) != 0:
        z = float(a[i][36:44])
        if z_max_1 >= z >= z_min_1:  # if the ions are between the two planes, i.e. they are a part of the pore
            Pore_Ions.append(count)

replacement = []
i = 0
while i < int(len(z_sol)*float(args.percent)/100):
    replacement.append(random.choice(z_sol))  # randomly choose an atom number from z_sol
    index = z_sol.index(replacement[i])  # get the index of the chosen number
    del z_sol[index]  # delete the number from z_sol so it can't be chosen again
    i += 1

for i in range(0, len(replacement)):
    a[replacement[i] + 1] = a[replacement[i] + 1][0:5] + a[replacement[i] + 1][5:10].replace('%s' % args.solvent, '%s' % args.Solute) + \
        a[replacement[i] + 1][10:15].replace('%s' % args.solvent, '%s' % args.Solute) + a[replacement[i] + 1][15:len(a[replacement[i] + 1])]

target = open('double_layer.gro', 'w')

for line in a:
    target.write(line)

exit()
# Now modify the topology accordingly

if os.path.isfile('%s' % args.top):  # check that the file exists first
    f = open(args.top, 'r')
    b = []
    for line in f:
        b.append(line)

    # find the [ molecules ] section
    mol_ind = 0
    while b[mol_ind].count('[ molecules ]') == 0:
        mol_ind += 1
    solv_ind = mol_ind + 1
    while b[solv_ind].count(args.solvent) == 0:
        solv_ind += 1
    sol_ind = mol_ind + 1
    while b[sol_ind].count(args.Solute) == 0:
        sol_ind += 1

    # Extract the number of molecules of solvent and solute present
    no_solvent = int(b[solv_ind][5:len(b[solv_ind])])
    no_solute = int(b[sol_ind][5:len(b[sol_ind])])

    # Now edit them to reflect the replacement of solvent with solute
    no_solvent -= len(replacement)
    no_solute += len(replacement)
    b[solv_ind] = b[solv_ind][0:5] + '               ' + str(no_solvent) + '\n' # spaces in between purely for visual formatting purposes
    b[sol_ind] = b[sol_ind][0:5] + '               ' + str(no_solute) + '\n'

    # Now rewrite the topology
    target = open('NaPore2.top', 'w')

    for line in b:
        target.write(line)

    target.close()
else:
    print 'The topology file, %s, does not exist' % args.top

# Create an index file
p1 = subprocess.Popen(["echo", "q"], stdout=subprocess.PIPE)
p2 = subprocess.Popen(["gmx", "make_ndx", "-f", "grad.gro", "-o", "indices.ndx"], stdin=p1.stdout)
p2.wait()

f = open('indices.ndx', 'r')
c = []
for line in f:
    c.append(line)

c.append('\n')
c.append('[ Pore_Ions ]' + '\n')
for i in range(0, len(Pore_Ions)):
    if i != 0 and (i + 1) % 13 == 0:
        c.append(str(Pore_Ions[i]) + '\n')
    else:
        c.append(str(Pore_Ions[i]) + ' ')

target = open('indices.ndx', 'w')

for line in c:
    target.write(line)

target.close()