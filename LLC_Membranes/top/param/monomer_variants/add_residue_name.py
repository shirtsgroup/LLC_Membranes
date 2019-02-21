#!/usr/bin/env python

input_pdb = "Molecule_X1.pdb"
residue_name = 'HX1'
output_pdb = "Molecule_X1_renamed.pdb"

with open(input_pdb, 'r') as f:
	a = []
	for line in f:
		a.append(line)

for i in range(len(a)):
	if a[i].count('ATOM') > 0:
		a[i] = a[i][:18] + residue_name + a[i][21:]

with open(output_pdb, 'w') as f:
	for line in a:
		f.write(line)
