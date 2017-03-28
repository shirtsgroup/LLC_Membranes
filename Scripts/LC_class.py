#!/usr/bin/env python


class LC(object):
    """A Liquid Crystal monomer has the following attributes which are relevant to building and crosslinking:

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

    def __init__(self, name, counterion, atoms, build_mon, images, c1_atoms, c2_atoms, tails, residues, valence, planeatoms,
                 lineatoms, ref_atom_index, no_vsites=0):
        self.name = name  # residue name
        self.counterion = counterion  # name of counter ion in monomer salt
        self.atoms = atoms  # number of atoms in the monomer residue (no counterion)
        self.build_mon = build_mon  # the name of the monomer coordinate file used for building
        self.images = images
        self.c1_atoms = c1_atoms  # c1 atoms (closest to tail ends) used in crosslinking vinyl groups
        self.c2_atoms = c2_atoms  # c2 atoms (adjacent to c1) used in crosslinking vinyl groups
        self.tails = tails  # number of tails
        self.residues = residues  # residues present in the full structure including counter ions
        self.valence = valence  # valence of counter ions
        self.planeatoms = planeatoms  # atoms defining a plane. Used for orientation during build
        self.lineatoms = lineatoms  # indices of atoms used to create a straight line while building. (see build.py)
        self.ref_atom_index = ref_atom_index  # atom index number used as a reference for rotation during build
        self.no_vsites = no_vsites
        self.tot_atoms = self.atoms + self.no_vsites
        self.topology = '%s.itp' % self.build_mon

    def vsites(self, no_vsites):
        """Returns the total number of atoms per monomer including virtual sites"""
        self.atoms += self.no_vsites
        return self.atoms

NAcarb8V = LC('HII', 'NA', 110, 'NAcarb8V.pdb', 9, ['C14', 'C26', 'C38'], ['C13', 'C25', 'C37'], 3, ['HII', 'NA'],
                1, ['C2', 'C15', 'C39'], [1, 2], 2)
NAcarb9V = LC('HII', 'NA', 119, 'NAcarb9V.pdb', 9, ['C15', 'C28', 'C41'], ['C14', 'C27', 'C40'], 3, ['HII', 'NA'],
                1, ['C2', 'C16', 'C42'], [1, 2], 2)
NAcarb10V = LC('HII', 'NA', 128, 'NAcarb10V.pdb', 9, ['C16', 'C30', 'C44'], ['C15', 'C29', 'C43'], 3, ['HII', 'NA'],
                1, ['C2', 'C17', 'C45'], [1, 2], 2)
NAcarb11V = LC('HII', 'NA', 137, 'monomer6.gro', 9, ['C20', 'C34', 'C48'], ['C19', 'C33', 'C47'], 3, ['HII', 'NA'], 1,
               ['C', 'C2', 'C4'], [0, 3], 9)
NAcarb11Vd = LC('HII', 'NA', 143, 'NAcarb11V_1_dummy.gro', 9, ['C20', 'C34', 'C48'], ['C19', 'C33', 'C47'], 3, ['HII', 'NA'],
                1, ['C', 'C2', 'C4'], [0, 3, 9], 9)
NAcarb11Vflat = LC('HII', 'NA', 137, 'flat.gro', 9, ['C20', 'C34', 'C48'], ['C19', 'C33', 'C47'], 3, ['HII', 'NA'], 1,
               ['C', 'C2', 'C4'], [0, 3], 9)

# SOL should be deleted. Need to check for usage first
SOL = LC('Water', 'H', 3, 'spc216.gro', 9, ['na'], ['na'], 1, ['OW', 'HW1', 'HW2'], 0, ['OW', 'HW1', 'HW2'], ['C', 'C3', 'C6'], 9)