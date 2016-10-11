#!/usr/bin/python

import argparse
import Thickness
import random
import os.path
import subprocess
import fix_gro

# Randomly replaces solvent molecules with solute like gmx solvate, but is capable of creating a concentration gradient

parser = argparse.ArgumentParser(description = 'Run Cylindricity script')

parser.add_argument('-i', '--input', default='wiggle.gro', help = 'Path to input file')
parser.add_argument('-p', '--percent', default=25, help = 'Percent solvent to replace with solute')
parser.add_argument('-l', '--LC', default='HII', help='Type of liquid crystal')
parser.add_argument('-s', '--solvent', default='Ar', help = 'Solvent which will be swapped with solute')
parser.add_argument('-S', '--Solute', default='NA', help = 'Name of atom replacing solvent molecules')
parser.add_argument('-t', '--top', default='NaPore.top', help = 'Name of topology file to be modified')
parser.add_argument('-d', '--double', default='on', help='Is this a double layer system? (on/off)')
parser.add_argument('-m', '--mdp', default='wiggle.mdp', help='Name of .mdp file which needs to be modified')

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

# Make lists of atoms numbers of interest

atoms_of_interest = ['%s' % args.LC, '%s' % args.solvent, '%s' % args.Solute]
top_lines = 0

while str.strip(a[top_lines][5:10]) not in atoms_of_interest:
    top_lines += 1

# generate index groups based on the single layer membrane
LC = 0  # list to hold atom number of atoms in membrane -- not including solvent
solvent = 0
solute = 0
no_atoms = 0
for i in range(top_lines, len(a) - 1):
    no_atoms += 1
    if str.strip(a[i][5:10]) == args.solvent:
        solvent += 1
    if str.strip(a[i][5:10]) == args.Solute:
        solute += 1
    if str.strip(a[i][5:10]) == args.LC:
        LC += 1

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
atom_no = []  # atom numbers of all
top_lines = 0
while a[top_lines].count('HII') == 0:
    top_lines += 1

count = -top_lines  # Things get tricky when there are greater than 100000 atoms so I count the number of lines as a way
# to keep track of atom numbers. The top line numbers are subtracted
for i in range(0, len(a)):
    count += 1
    if a[i].count(args.solvent) != 0:
        z = float(a[i][36:44])
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

fix_gro.reorder(a, 'double_layer.gro')  # reorders and then rewrites the file with

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
    LC_ind = mol_ind + 1
    while b[LC_ind].count(args.LC) == 0:
        LC_ind += 1
    solv_ind = mol_ind + 1
    while b[solv_ind].count(args.solvent) == 0:
        solv_ind += 1
    sol_ind = mol_ind + 1
    while b[sol_ind].count(args.Solute) == 0:
        sol_ind += 1

    # Extract the number of molecules of solvent and solute present
    if args.double == 'on':
        multiplier = 2
    else:
        multiplier = 1

    no_LC = int(b[LC_ind][5:len(b[LC_ind])])*multiplier
    no_solvent = int(b[solv_ind][5:len(b[solv_ind])])*multiplier
    no_solute = int(b[sol_ind][5:len(b[sol_ind])])*multiplier

    # Now edit them to reflect the replacement of solvent with solute
    no_solvent -= len(replacement)
    no_solute += len(replacement)
    b[LC_ind] = b[LC_ind][0:5] + '                ' + str(no_LC) + '\n'  # spaces in between purely for aesthetics
    b[solv_ind] = b[solv_ind][0:5] + '               ' + str(no_solvent) + '\n'
    b[sol_ind] = b[sol_ind][0:5] + '               ' + str(no_solute) + '\n'

    # Now rewrite the topology
    target = open('NaPore2.top', 'w')

    for line in b:
        target.write(line)
    target.close()

else:
    print 'The topology file, %s, does not exist' % args.top

# generate index groups for the soon to be made top membrane layer. They are just a translation of the bottom layer.
# Using gmx make_ndx with the -twin flag would do something similar but I don't like all the default groups that are
# made so I am making an index file from scratch

# In the following, I take advantage of the fact that among atoms of the same type, fix_gro.py writes them in order of
# their appearance in the disorganized .gro file

# Create an index file with default groups
# p1 = subprocess.Popen(["echo", "q"], stdout=subprocess.PIPE)
# p2 = subprocess.Popen(["gmx", "make_ndx", "-f", "double_layer.gro", "-o", "indices.ndx"], stdin=p1.stdout)
# p2.wait()

f = open('double_layer.gro', 'r')

a = []
for line in f:
    a.append(line)
f.close()

top_lines = 0
while str.strip(a[top_lines][5:10]) not in atoms_of_interest:
    top_lines += 1

LC_start = -top_lines
Solvent_start = -top_lines
Solute_start = -top_lines

while a[LC_start + top_lines].count('%s' % args.LC) == 0:
    LC_start += 1
while a[Solvent_start + top_lines].count('%s' % args.solvent) == 0:
    Solvent_start += 1
while a[Solute_start + top_lines].count('%s' % args.Solute) == 0:
    Solute_start += 1

bot_membrane = []
top_membrane = []
for i in range(LC_start, LC + LC_start):
    bot_membrane.append(i + 1)
    top_membrane.append(i + 1 + LC)

# Now add the ions associated with the bottom membrane
for i in range(Solute_start, solute + Solute_start):
    bot_membrane.append(i + 1)

# Now add the ions associated with the top membrane
start = Solute_start + solute + 1 + top_lines
while str.strip(a[start][5:10]) == args.Solute and z_max_1 < float(a[start][36:44]) < z_min_2:
    start += 1

start -= top_lines

for i in range(start, start + solute):
    top_membrane.append(i + 1)

solute_ind = []  # This one is more complicated since the location are a little less organized. Instead, we index them
# based on the z positions of the solute
count = -top_lines
for i in range(0, len(a)):
    count += 1
    if str.strip(a[i][5:10]) == args.Solute and z_max_1 < float(a[i][36:44]) < z_min_2:
        solute_ind.append(count)

solvent_ind = []
for i in range(Solvent_start, multiplier*solvent + Solvent_start - len(solute_ind)):  # need to subtract out the solute
    # that was turned to solvent and multiply by the multiplier if the system was stacked into a duplicate layer
    if z_max_1 < float(a[i][36:44]) < z_min_2:
        solvent_ind.append(i + 1)

System = []
for i in range(0, no_atoms*multiplier):
    System.append(i + 1)

ndx = []


def add_index(ndx, list, name):
    if len(ndx) != 0:
        ndx.append('\n' + '\n')
    ndx.append('[ %s ]' % name + '\n')
    for i in range(0, len(list)):
        if i != 0 and (i + 1) % 13 == 0:  # looks complicated just so the formatting looks nicer
            ndx.append(str(list[i]) + '\n')
        else:
            ndx.append(str(list[i]) + ' ')
    return ndx

ndx = add_index(ndx, System, 'System')
ndx = add_index(ndx, top_membrane, 'Top_Membrane')
ndx = add_index(ndx, bot_membrane, 'Bottom_Membrane')
ndx = add_index(ndx, solute_ind, 'Solute')
ndx = add_index(ndx, solvent_ind, 'Solvent')
ndx.append('\n')  # index file should end with a newline

target = open('indices.ndx', 'w')

for line in ndx:
    target.write(line)

target.close()

# Now set up the .mdp file for Computational Electrophysiology

f = open('%s' % args.mdp, 'r')

a = []
for line in f:
    a.append(line)

f.close()
# find the temperature coupling groups
tc_grp_index = 0
while a[tc_grp_index].count('tc-grps') == 0:
    tc_grp_index += 1

tau_t_index = 0
while a[tau_t_index].count('tau-t') == 0:
    tau_t_index += 1

ref_t_index = 0
while a[ref_t_index].count('ref-t') == 0:
    ref_t_index += 1

# The following has definite room for improvement
# a[tc_grp_index] = 'tc-grps = Top_Membrane Bottom_Membrane Solute Solvent' + '\n'
# a[tau_t_index] = 'tau-t = 0.1 0.1 0.1 0.1' + '\n'
# a[ref_t_index] = 'ref-t = 300 300 300 300' + '\n'
a[tc_grp_index] = 'tc-grps = System' + '\n'
a[tau_t_index] = 'tau-t = 0.1' + '\n'
a[ref_t_index] = 'ref-t = 300' + '\n'

a.append(';Computation Electrophysiology Parameters' + '\n')
a.append('swapcoords     = Z' + '\n')
a.append('swap-frequency = 50' + '\n')
a.append('split-group0   = Bottom_Membrane' + '\n')
a.append('split-group1   = Top_Membrane' + '\n')
a.append('swap-group     = Solute' + '\n')
a.append('solvent-group  = Solvent' + '\n')

f = open('compel.mdp', 'w')

for line in a:
    f.write(line)
f.close()