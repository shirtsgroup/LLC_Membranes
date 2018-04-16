#!/usr/bin/env python

import argparse
import mdtraj as md
import os


def initialize():

    parser = argparse.ArgumentParser(description='Generate topology file from coordinate file. This only works if all'
                                                 'residues have a topology in HII/top')

    parser.add_argument('-g', '--gro', help='Name of coordinate file to write topology file for')
    parser.add_argument('-o', '--output', default='topol.top', help='Name of topology to output')
    parser.add_argument('-d', '--description', default='Simulation box', help='Description of system to put under '
                                                                              '[ system ] directive')
    parser.add_argument('-ff', '--forcefield', default='gaff', help='Name of forcefield to use')

    args = parser.parse_args()

    return args


class SystemTopology(object):

    def __init__(self, gro, ff='gaff'):
        """
        :param gro: (str) coordinate file for which to create topology
        :param ff: (str) forcefield to use (default=gaff)
        """

        t = md.load(gro)  # load coordinate file

        self.script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        self.top_location = "%s/../top/topologies" % self.script_location  # topology location
        self.ff_location = "%s/../top/Forcefields" % self.script_location  # forcefield location
        self.forcefield = ff  # which forcefield to use

        with open('%s/ions.txt' % self.top_location) as f:
            self.ions = []
            for line in f:
                if line[0] != '#':
                    self.ions.append(str.strip(line))

        residues = [a.residue.name for a in t.topology.atoms]

        for i, r in enumerate(residues):
            if r == 'HOH':
                residues[i] = 'SOL'

        unique_residues = list(set(residues))

        residue_count = {}
        for i in unique_residues:
            if i in self.ions:
                natoms = 1
            else:
                natoms = md.load('%s/%s.gro' % (self.top_location, i)).xyz.shape[1]  # number of atoms that make up a single residue of type i
            residue_count[i] = residues.count(i) // natoms

        self.residues = unique_residues
        self.residue_count = residue_count

    def write_top(self, name='topol.top', description='Simulation Box'):

        top = []
        top.append(';Forcefield\n')
        top.append('#include "%s/%s/%s.itp"\n' %(self.ff_location, self.forcefield, self.forcefield))
        top.append('\n')

        ion_count = 0
        for r in self.residues:
            if r in self.ions:
                if ion_count == 0:
                    top.append(';Ion Topology\n')
                    top.append('#include "%s/ions.itp"\n' % self.top_location)
                    top.append('\n')
            else:
                top.append(';%s Topology\n' % r)
                top.append('#include "%s/%s.itp"\n' % (self.top_location, r))
                top.append('\n')

        top.append('[ system ]\n')
        top.append('%s\n' % description)
        top.append('\n')

        top.append('[ molecules ]\n')
        top.append(';Compounds     nmols\n')

        for r in self.residues:
            top.append('{:10s}{:>10d}\n'.format(r, self.residue_count[r]))

        with open(name, 'w') as f:
            for line in top:
                f.write(line)

    def add_residue(self, residue, n=1, write=False, topname='topol.top', top_description='Simulation Box'):
        """
        Add molecule(s) of a single residue to the topology
        :param residue : name of residue object
        :param n : number of molecules to add
        :param write : write new topology file
        :param topname : name of output topology if written
        :param top_description : system description to be written into top if desired
        """

        resname = residue.resname
        if residue.is_ion:
            self.ions.append(resname)

        if resname not in self.residues:
            self.residues.append(resname)
            self.residue_count[resname] = n
        else:
            self.residue_count[resname] += n

        if write:
            self.write_top(name=topname, description=top_description)


if __name__ == "__main__":

    args = initialize()
    t = SystemTopology(args.gro, ff=args.forcefield)
    t.write_top(name=args.output, description=args.description)