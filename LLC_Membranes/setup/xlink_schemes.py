#!/usr/bin/env python

import sys
import numpy as np
import importlib

""" Specify reactants and products of cross-linking reactions
"""


class UndefinedReactionError(Exception):
    """ Raised if an undefined reaction is attempted """

    def __init__(self, message):

        super().__init__(message)


class XlinkReaction:

    def __init__(self, monomer):
        """ Define the cross-linking reaction scheme based on the monomer.

        :param monomer: name of monomer being cross-linked. It should have its own class defined with attributes.

        :type monomer: str

        Attributes:
            scheme: a class that defines the reaction scheme followed by the monomer
        """

        if monomer == 'MOL':  # maybe use a dict here instead
            monomer = 'Dibrpyr14'

        if monomer == 'Dibrpyr14':
            self.scheme = DieneScheme(monomer)
        else:
            self.scheme = None


class DieneScheme:

    def __init__(self, monomer):
        """ Define the cross-linking reaction of a monomer with terminal diene tails

        :param monomer: name of monomer

        :type monomer: str

        Attributes:
            **monomer**: a class defining the constituent atoms. The name of the class is the same as the monomer

            **weights**: dictionary with the relative weights of one type of reaction. For example, a C1 and C2 atom \
            may be chosen to bond but there may be multiple types of addition (e.g. head-to-tail and 1,4 addition)

            **reaction_weights**: the relative weights of all types of reactions. For example, a C1 and C2 may be \
            eligible to react and perhaps another C2 is eligible to react with a C3. These weights control the number \
            of each type of distinct reaction so that head-to-tail dominates over 13 addition for example.

            **radical_reaction_weights**: the same as reaction weights, but applied to the radical reactions.
        """

        self.monomer = getattr(importlib.import_module("LLC_Membranes.setup.xlink_schemes"), '%s' % monomer)()

        # keys: type of carbon carbon bond; values: dict with (keys: name of reaction, values: likelihood of reaction)
        self.weights = {'C1-C2': {'head2tail': 0.5, 'head2head': 0.5},
                        'C1-C2_radical': {'radical_c2': 1}}  # probabilities must sum to 1
        self.reaction_weights = {'head2tail': 0.5, 'head2head': 0.5, '13addition': 0}
        self.radical_reaction_weights = {'radical_c2': 1.}

    def determine_chains(self, c):
        """ Determine to which chain carbon atoms come from

        :param c: list of names or single name of carbon atoms

        :type c: list, str
        """

        if isinstance(c, str):
            c = [c]

        chains = [None for _ in c]
        for k in self.monomer.chains.keys():
            for i, x in enumerate(c):
                if x in self.monomer.chains[k].keys():
                    chains[i] = k

        return chains

    def determine_reaction_type(self, c1, c2, radical=False):
        """ Determine the type of reaction based on which carbons are bonding

        :param c1: name of carbon in first chain
        :param c2: name of carbon in second chain
        :param radical: True if this reaction involves a radical

        :type c1: str
        :type c2: str
        :type radical: bool
        """

        chain1, chain2 = self.determine_chains([c1, c2])
        bond = '%s-%s' % (self.monomer.chains[chain1][c1], self.monomer.chains[chain2][c2])

        if radical:
            bond += '_radical'

        try:
            return np.random.choice(list(self.weights[bond].keys()), p=list(self.weights[bond].values()))
        except KeyError:
            return False  # bond between chain1 and chain2 is not defined

    def react(self, reaction_type, atoms):

        if reaction_type == 'head2tail':

            return self.head2tail(atoms)

        elif reaction_type == 'radical_c2':

            return self.radical_c2(atoms)

        elif reaction_type == 'terminate':

            return self.terminate(atoms)

        elif reaction_type == 'head2head':

            return self.head2head(atoms)

        else:
            raise UndefinedReactionError('Reaction %s undefined' % reaction_type)  # replace with custom exception

    def head2tail(self, atoms):
        """ Define the head-to-tail addition of the terminal double bonds (c1 -- c2)

        :param atoms: dictionary of atom names (keys) and their index (values) in the context of an entire unit cell

        :type atoms: dict

        :return: atom indices and their corresponding types after reaction
        """

        c1, c2 = atoms.keys()
        c1_ndx, c2_ndx = atoms.values()

        chain1, chain2 = self.determine_chains([c1, c2])

        # to get indexing right
        c1_ndx -= self.monomer.indices[chain1]['C1']
        c2_ndx -= self.monomer.indices[chain2]['C2']

        types = {'chain1': {'C1': 'c3', 'C2': 'c2', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'ha',
                            'H4': 'ha', 'H5': 'ha'},
                 'chain2': {'C1': 'c3', 'C2': 'c3', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'hc',
                            'H4': 'ha', 'H5': 'ha', 'D1': 'hc'}}

        # update types
        reacted_types = {'chain1': {c1_ndx + self.monomer.indices[chain1][a]: types['chain1'][a] for a in
                                    types['chain1'].keys()},
                         'chain2': {c2_ndx + self.monomer.indices[chain2][a]: types['chain2'][a] for a in
                                    types['chain2'].keys()}}

        # bond between carbons. Format [c1, c2, type]
        bonds = [[c1_ndx + self.monomer.indices[chain1]['C1'], c2_ndx + self.monomer.indices[chain2]['C2'], 'carbon']]

        # dummy bonds - 1 new bond between dummy atoms and carbon
        bonds += [[c2_ndx + self.monomer.indices[chain2]['C1'], c2_ndx + self.monomer.indices[chain2]['D1'], 'dummy']]

        # define indices of left-over radicals
        radicals = [c1_ndx + self.monomer.indices[chain1]['C2']]

        chain1_impropers = ['C1']
        chain2_impropers = ['C1', 'C2']
        rm_improper = []
        for c in chain1_impropers:
            rm_improper.append([c1_ndx + self.monomer.indices[chain1][x] for x in self.monomer.impropers[chain1][c]])
        for c in chain2_impropers:
            rm_improper.append([c2_ndx + self.monomer.indices[chain2][x] for x in self.monomer.impropers[chain2][c]])

        # define terminated atoms
        terminated = [c1_ndx + self.monomer.indices[chain1]['C1'], c2_ndx + self.monomer.indices[chain2]['C2'], c2_ndx +
                      self.monomer.indices[chain2]['C1']]

        return reacted_types, bonds, radicals, rm_improper, terminated

    def head2head(self, atoms):
        """ Define the head-to-head addition of the terminal double bonds (c1 -- c1). Note that this reaction actually
        occurs based on proximity of c1 and c2 atoms

        :param atoms: dictionary of atom names (keys) and their index (values) in the context of an entire unit cell

        :type atoms: dict

        :return: atom indices and their corresponding types after reaction
        """

        c1, c2 = atoms.keys()
        c1_ndx, c2_ndx = atoms.values()

        chain1, chain2 = self.determine_chains([c1, c2])

        # to get indexing right
        c1_ndx -= self.monomer.indices[chain1]['C1']
        c2_ndx -= self.monomer.indices[chain2]['C2']

        types = {'chain1': {'C1': 'c3', 'C2': 'c2', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'ha',
                            'H4': 'ha', 'H5': 'ha'},
                 'chain2': {'C1': 'c3', 'C2': 'c2', 'C3': 'c2', 'C4': 'c3', 'H1': 'hc', 'H2': 'hc', 'H3': 'ha',
                            'H4': 'ha', 'H5': 'hc', 'D4': 'hc'}}

        # update types
        reacted_types = {'chain1': {c1_ndx + self.monomer.indices[chain1][a]: types['chain1'][a] for a in
                                    types['chain1'].keys()},
                         'chain2': {c2_ndx + self.monomer.indices[chain2][a]: types['chain2'][a] for a in
                                    types['chain2'].keys()}}

        # bond between carbons
        bonds = [[c1_ndx + self.monomer.indices[chain1]['C1'], c2_ndx + self.monomer.indices[chain2]['C1'], 'carbon']]

        # dummy bonds - 1 new bond between dummy atoms and carbon
        bonds += [[c2_ndx + self.monomer.indices[chain2]['C4'], c2_ndx + self.monomer.indices[chain2]['D4'], 'dummy']]

        # define indices of left-over radicals
        radicals = [c1_ndx + self.monomer.indices[chain1]['C2']]

        chain1_impropers = ['C1']  # [1]
        chain2_impropers = ['C1', 'C4']  # [1, 2]
        rm_improper = []
        for c in chain1_impropers:
            rm_improper.append([c1_ndx + self.monomer.indices[chain1][x] for x in self.monomer.impropers[chain1][c]])
        for c in chain2_impropers:
            rm_improper.append([c2_ndx + self.monomer.indices[chain2][x] for x in self.monomer.impropers[chain2][c]])

        # define terminated atoms
        terminated = [c1_ndx + self.monomer.indices[chain1]['C1'], c2_ndx + self.monomer.indices[chain2]['C1'],
                      c2_ndx + self.monomer.indices[chain2]['C2']]  # C2 terminated for now even though still alkene

        return reacted_types, bonds, radicals, rm_improper, terminated

    def radical_c2(self, atoms):
        """ Define the reaction of a radical c2 with unreacted c1

        :param atoms: dictionary of atom names (keys) and their index (values) in the context of an entire unit cell

        :type atoms: dict

        :return: reacted_types, dummy_bonds, radicals, rm_improper, terminated
        """

        c1, c2 = atoms.keys()
        c1_ndx, c2_ndx = atoms.values()

        chain1, chain2 = self.determine_chains([c1, c2])

        # to get indexing right
        c1_ndx -= self.monomer.indices[chain1]['C1']
        c2_ndx -= self.monomer.indices[chain2]['C2']

        # types after reaction
        types = {'chain1': {'C1': 'c3', 'C2': 'c2', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'ha',
                            'H4': 'ha', 'H5': 'ha'},  # chain1 contains c1
                 'chain2': {'C1': 'c3', 'C2': 'c3', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'hc',
                            'H4': 'ha', 'H5': 'ha'}}  # chain2 contains c2 radical

        # update types
        reacted_types = {'chain1': {c1_ndx + self.monomer.indices[chain1][a]: types['chain1'][a]
                                    for a in types['chain1'].keys()},
                         'chain2': {c2_ndx + self.monomer.indices[chain2][a]: types['chain2'][a]
                                    for a in types['chain2'].keys()}}

        # new bonds
        bonds = [[c1_ndx + self.monomer.indices[chain1]['C1'], c2_ndx + self.monomer.indices[chain2]['C2'], 'carbon']]

        # no dummy bonds to add

        # define indices of left-over radicals
        radicals = [c1_ndx + self.monomer.indices[chain1]['C2']]

        chain1_impropers = ['C1']  # [1]
        chain2_impropers = ['C2']  # [2]
        rm_improper = []
        for c in chain1_impropers:
            rm_improper.append([c1_ndx + self.monomer.indices[chain1][x] for x in self.monomer.impropers[chain1][c]])
        for c in chain2_impropers:
            rm_improper.append([c2_ndx + self.monomer.indices[chain2][x] for x in self.monomer.impropers[chain2][c]])

        # define terminated atoms
        terminated = [c1_ndx + self.monomer.indices[chain1]['C1'], c2_ndx + self.monomer.indices[chain2]['C2']]

        return reacted_types, bonds, radicals, rm_improper, terminated

    def terminate(self, atoms):
        """ Define the termination reaction. i.e. a dummy atom attaches to a radical carbon atoms. The hybridization of
        the carbon changes to sp3

        :param atoms: dictionary with atom name (key) and its index (value) in the context of an entire unit cell

        :type atoms: dict

        :return: reacted_types, dummy_bonds, radicals, rm_improper, terminated
        """

        c = list(atoms.keys())[0]  # name of carbon atom being terminated
        c_ndx = list(atoms.values())[0]  # serial index of carbon begin terminated

        chain = self.determine_chains(c)[0]  # which chain carbon atom is on
        c_name = self.monomer.chains[chain][c]

        # to get indexing right
        c_ndx -= self.monomer.indices[chain][c_name]

        # types after reaction. Keeping this dictionary format so it integrates easily with xlinking algorithm
        types = {'chain': {self.monomer.chains[chain][c]: 'c3', self.monomer.dummy_connectivity[chain][c]: 'hc'}}

        for i in self.monomer.hydrogen_connectivity[c]:  # turn already attached carbon(s) to c3
            types['chain'][i] = 'hc'

        # update types
        reacted_types = {'chain': {c_ndx + self.monomer.indices[chain][a]: types['chain'][a]
                                   for a in types['chain'].keys()}}

        # add dummy atom bond
        bonds = [[c_ndx + self.monomer.indices[chain]['C2'], c_ndx + self.monomer.indices[chain]['D2'], 'dummy']]

        radicals = []

        rm_improper = [[c_ndx + self.monomer.indices[chain][x] for x in self.monomer.impropers[chain][c_name]]]

        # define terminated atoms
        terminated = [c_ndx + self.monomer.indices[chain][c_name]]

        return reacted_types, bonds, radicals, rm_improper, terminated


class Dibrpyr14:

    def __init__(self):
        """ Define monomer attributes necessary for cross-linking
        TODO: clean this up! Some attributes can be combined or written more intuitively
        Attributes:
            **chains**: a dictionary of dictionaries of the form {'a': {'C': 'C1', 'C1': 'C2 ...}, 'b': {'C45': 'C1, \
            'C44': 'C2' ...}}. This monomer contains two alkane chain tails. 'a' and 'b' are labels referring to each \
            chain respectively. Each sub-dictionary is used to map the names of relevant monomer atoms to generalized \
            atom names that are recognized by the reaction schemes. For a diene tail, C1 refers to the terminal carbon \
            atom, C2 refers to the second-to-last carbon atom and so on. H1 and H2 refer to hydrogen atoms attached to \
            C1. H3 is attached to C2, H4 to C3 and H5 to C4. These generalized names are the values in each \
            sub-dictionary. The keys are the names of the carbon atoms in the context of the monomer, so the names as \
            they appear in the .gro file. Take a look at Dibrpyr14.gro if this is still unclear.

            **nchains**: the number of chains that a monomer has. This makes code more readable elsewhere

            **chain_numbers**: Each chain should have its own number starting from 0. This is used to determine \
            whether a given chain has already been involved in a reaction.

            **indices**: this dictionary of dictionaries stores the indices of relevant monomer atoms in the context \
            of a single monomer. It is constructed in a similar way to self.chains, but the generalized carbon atoms \
            names are the keys (instead of the values) and the index is the index of generalized carbon atom. For \
            example, C1 in chain 'b' is index 49, meaning it is the 50th atom in the topology file for a single \
            Dibrpyr14 residue. (See Dibrpyr14.gro or Dibrpyr14.itp to confirm this and all other entries). The indices \
            include those corresponding to the dummy atoms (D1, D2, D3, D4).

            **dummy_connectivity**: this dictionary of dictionaries defines the connectivity between carbon atoms and \
            the dummy atoms. For each chain, the sub-dictionary contains the name of each carbon atom (in the context \
            of the monomer) as keys and the generalized name of the attached dummy atom as the value. In general, for \
            the diene reaction, D1 is attached to C1, D2 to C2 and so on.

            **hydrogen_connectivity**: similar to dummy_connectivity, this defines which hydrogen atoms are bonded to \
            which carbon. This dictionary does not distinguish between chains. Each key is the name of a carbon atom, \
            as you'd find in the corresponding Dibrpyr14 .gro or .itp file, with a list of generalized hydrogen atoms \
            that are attached to that carbon.

            **dummy_mass**: the mass of the dummy atom (hydrogen)

            **carbons**: a dictionary with keys that are the generalized carbon atom names and values that are a list \
            of the actual carbon atom names as you'd find them in the .gro or .itp file.

            **bonds_with**: a list of carbon pairs. Each sub-list contains two lists: a list of carbon atoms in the \
            1st entry that can bond with carbon atoms in the second list entry.

            **impropers**: a dictionary of dictionaries. Each sub-dictionary refers to a chain. The sub-directories \
            contain keys corresponding to each carbon atom. Each value is a list of atoms that make up the improper \
            dihedral which holds the originally sp2 carbon in a planar configuration. When an sp2 carbon is converted \
            to sp3, the improper must be removed. The order of the improper is written in the same way that you'd find \
            it in Dibrpyr14.itp. Required ordering was a choice that I made to speed up some calculations. Be very \
            careful defining this attribute.
        """

        # names of atoms that make up relevant segements of each chain
        self.chains = {'a': {'C': 'C1', 'C1': 'C2', 'C2': 'C3', 'C3': 'C4', 'C4': 'C5', 'H': 'H1', 'H1': 'H2',
                             'H2': 'H3', 'H3': 'H4', 'H4': 'H5'},
                       'b': {'C45': 'C1', 'C44': 'C2', 'C43': 'C3', 'C42': 'C4', 'C41': 'C5', 'H81': 'H1', 'H80': 'H2',
                             'H79': 'H3', 'H78': 'H4', 'H77': 'H5'}
                       }

        self.nchains = len(list(self.chains.keys()))

        self.chain_numbers = {'a': 0, 'b': 1}  # used to number chains

        # self.initial_types = {'C1': 'c2', 'C2': 'ce', 'C3': 'ce', 'C4': 'c2', 'H1': 'ha', 'H2': 'ha', 'H3': 'ha',
        #                       'H4': 'ha', 'H5': 'ha'}

        # all indices numbered from 0. D1, D2, ... correspond to dummies attached to C1, C2, ... respectively
        self.indices = {'a': {'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'C5': 4, 'H1': 52, 'H2': 53, 'H3': 54, 'H4': 55,
                              'H5': 56, 'D1': 136, 'D2': 137, 'D3': 138, 'D4': 139},
                        'b': {'C1': 49, 'C2': 48, 'C3': 47, 'C4': 46, 'C5': 45, 'H1': 133, 'H2': 132, 'H3': 131,
                              'H4': 130, 'H5': 129, 'D1': 140, 'D2': 141, 'D3': 142, 'D4': 143}
                        }

        self.dummy_connectivity = {'a': {'C': 'D1', 'C1': 'D2', 'C2': 'D3', 'C3': 'D4'},
                                   'b': {'C45': 'D1', 'C44': 'D2', 'C43': 'D3', 'C42': 'D4'}}

        self.hydrogen_connectivity = {'C': ['H1', 'H2'], 'C1': ['H3'], 'C2': ['H4'], 'C3': ['H5'],
                                      'C45': ['H1', 'H2'], 'C44': ['H3'], 'C43': ['H4'], 'C42': ['H5']}

        self.dummy_mass = 1.008  # mass of hydrogen

        # write these in order of priority
        # for efficiency, don't repeat things. For example self.carbons['C1']: self.carbons['C2'] is the same as
        # self.carbons['C2']: self.carbons['C1']. Otherwise, computational expense goes up and a new reaction has
        # to be defined below.
        self.carbons = {'C1': ['C', 'C45'], 'C2': ['C1', 'C44'], 'C3': ['C2', 'C43'], 'C4': ['C3', 'C42']}
        self.bonds_with = [[self.carbons['C1'], self.carbons['C2']]]

        # define which improper dihedrals to remove -- written in same order as .itp file!!!
        # note that the order of the atoms may be different for each chain
        # NOTE: C3 not tested
        self.impropers = {'a': {'C1': ['H2', 'C1', 'H1', 'C2'], 'C2': ['C1', 'C3', 'C2', 'H3'],
                          'C3': ['C4', 'C2', 'C3', 'H4'], 'C4': ['C5', 'C3', 'C4', 'H5']},
                          'b': {'C1': ['C2', 'H2', 'C1', 'H1'], 'C2': ['C1', 'C3', 'C2', 'H3'],
                          'C3': ['C4', 'C2', 'C3', 'H4'], 'C4': ['C5', 'C3', 'C4', 'H5']}}
