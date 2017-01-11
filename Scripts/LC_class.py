#!/usr/bin/python


class LC(object):
    """A Liquid Crystal monomer has the following attributes which are relevant to crosslinking:

    Attributes:
        name: A string representing the monomer's name.
        atoms: An integer accounting for the number of atoms in a single monomer.
        build_mon: Monomer used to build the unit cell
        images: Number of periodic images to be used in calculations
        c1_atoms: A list of atoms which will be involved in crosslinking as 'c1' -- See xlink.py
        c2_atoms: A list of atoms which will be involved in crosslinking as 'c2' -- See xlink.py
        tails: Number of tails on each monomer
        residues: A list of the minimum residue names present in a typical structure
        no_vsites: A string indicating whether there are dummy atoms associated with this monomer.
    """

    def __init__(self, name, atoms, build_mon, images, c1_atoms, c2_atoms, tails, residues, no_vsites=0):
        self.name = name
        self.atoms = atoms
        self.build_mon = build_mon
        self.no_vsites = no_vsites
        self.images = images
        self.c1_atoms = c1_atoms
        self.c2_atoms = c2_atoms
        self.tails = tails
        self.residues = residues
        self.tot_atoms = self.atoms + self.no_vsites
        self.topology = '%s.itp' % self.build_mon


    def vsites(self, no_vsites):
        """Returns the total number of atoms per monomer including virtual sites"""
        self.atoms += self.no_vsites
        return self.atoms

HII = LC('Hexagonal', 137, 'monomer2', 9, ['C20', 'C34', 'C48'], ['C19', 'C33', 'C47'], 3, ['HII', 'NA'])
SOL = LC('Water', 3, 'spc216.gro', 9, ['na'], ['na'], 1, ['OW', 'HW1', 'HW2'])