#!/usr/bin/env python
import argparse
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Add missing fields to .pdb for fullrmc run')

    parser.add_argument('-i', '--pdb', default='wiggle.pdb', help='Name of .pdb file to modify')

    args = parser.parse_args()

    return args


def get_element(atom):

    if atom[0] == 'C':
        element = 'C'
    elif atom[0] == 'O':
        element = 'O'
    elif atom[0] == 'H':
        element = 'H'
    elif atom == 'NA':
        element = 'na'
    else:
        print("Element %s is not defined" % atom)
        exit()

    return element

if __name__ == "__main__":

    args = initialize()

    a = []
    with open(args.pdb, 'r') as f:
        for line in f:
            a.append(line)

    # lazy way to get box vectors
    t = md.load(args.pdb)
    box = t.unitcell_vectors[0, :, :]*10

    atom_start = 0
    while a[atom_start].count('ATOM') == 0:
        atom_start += 1

    for i in range(atom_start, atom_start + t.n_atoms):
        atom = str.strip(a[i][12:15])
        sequence_number = a[i][22:26]
        element = get_element(atom)
        a[i] = a[i][:72] + sequence_number + element + '\n'

    end = -1
    while a[end].count('ATOM') == 0:
        end -= 1

    with open("initial_elements.pdb", 'w') as f:

        f.write(a[0])
        f.write(a[1])
        f.write("REMARK    Boundary Conditions: {:2.4f} {:2.4f} {:2.4f} {:2.4f} {:2.4f} {:2.4f} {:2.4f} {:2.4f} "
                "{:2.4f}\n".format(box[0, 0], box[0, 1], box[0, 2], box[1, 0], box[1, 1], box[1, 2], box[2, 0],
                                   box[2, 1], box[2, 2]))
        if end != -1:
            for line in a[2:end+1]:
                f.write(line)
        else:
            for line in a[2:]:
                f.write(line)
