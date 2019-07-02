#!/usr/bin/env python

import sys
import numpy as np

""" Specify reactants and products of cross-linking reactions
"""


class UndefinedReactionError(Exception):
    """ Raised if an undefined reaction is attempted """

    def __init__(self, message):

        super().__init__(message)


class XlinkReaction:

    def __init__(self, monomer):

        if monomer in ['Dibrpyr14', 'MOL']:
            self.scheme = DieneScheme(monomer)
        else:
            self.scheme = None


class DieneScheme:

    def __init__(self, monomer):
        """ Define the cross-linking reaction of a monomer with terminal diene tails

        :param monomer: name of monomer

        :type monomer: str
        """

        # keys: type of carbon carbon bond; values: dict with (keys: name of reaction, values: likelihood of reaction)
        self.weights = {'C1-C2': {'head2tail': 1., '14addition': 0},
                        'C1-C2_radical': {'radical_c2': 1}}  # probabilities must sum to 1
        self.reaction_weights = {'head2tail': 1., '14addition': 0, '13addition': 0}
        self.radical_reaction_weights = {'radical_c2': 1.}

        if monomer in ['Dibrpyr14', 'MOL']:  # could build this into annotations. Might be better like this

            # names of atoms that make up relevant segements of each chain
            self.chains = {'a': {'C': 'C1', 'C1': 'C2', 'C2': 'C3', 'C3': 'C4', 'H': 'H1', 'H1': 'H2', 'H2': 'H3',
                                 'H3': 'H4', 'H4': 'H5'},
                           'b': {'C45': 'C1', 'C44': 'C2', 'C43': 'C3', 'C42': 'C4', 'H81': 'H1', 'H80': 'H2',
                                 'H79': 'H3', 'H78': 'H4', 'H77': 'H5'}
                           }

            self.nchains = len(list(self.chains.keys()))

            self.chain_numbers = {'a': 0, 'b': 1}  # used to number chains

            self.carbons = {'C1': ['C', 'C45'], 'C2': ['C1', 'C44'], 'C3': ['C2', 'C43'], 'C4': ['C3', 'C42']}

            self.initial_types = {'C1': 'c2', 'C2': 'ce', 'C3': 'ce', 'C4': 'c2', 'H1': 'ha', 'H2': 'ha', 'H3': 'ha',
                                  'H4': 'ha', 'H5': 'ha'}

            # all indices numbered from 0. D1, D2, ... correspond to dummies attached to C1, C2, ... respectively
            self.indices = {'a': {'C1': 0, 'C2': 1, 'C3': 2, 'C4': 3, 'H1': 52, 'H2': 53, 'H3': 54, 'H4': 55, 'H5': 56,
                                  'D1': 136, 'D2': 137, 'D3': 138, 'D4': 139},
                            'b': {'C1': 49, 'C2': 48, 'C3': 47, 'C4': 46, 'H1': 133, 'H2': 132, 'H3': 131, 'H4': 130,
                                  'H5': 129, 'D1': 140, 'D2': 141, 'D3': 142, 'D4': 143}
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
            self.bonds_with = [[self.carbons['C1'], self.carbons['C2']]]

    def determine_chains(self, c):
        """ Determine to which chain carbon atoms come from

        :param c: list of names or single name of carbon atoms

        :type c: list, str
        """

        if isinstance(c, str):
            c = [c]

        chains = [None for _ in c]
        for k in self.chains.keys():
            for i, x in enumerate(c):
                if x in self.chains[k].keys():
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
        bond = '%s-%s' % (self.chains[chain1][c1], self.chains[chain2][c2])

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

        else:
            raise UndefinedReactionError('Reaction %s Undefined' % reaction_type)  # replace with custom exception

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
        c1_ndx -= self.indices[chain1]['C1']
        c2_ndx -= self.indices[chain2]['C2']

        types = {'chain1': {'C1': 'c3', 'C2': 'c2', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'ha',
                            'H4': 'ha', 'H5': 'ha'},
                 'chain2': {'C1': 'c3', 'C2': 'c3', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'hc',
                            'H4': 'ha', 'H5': 'ha', 'D1': 'hc', 'D2': 'hc'}}

        # update types
        reacted_types = {'chain1': {c1_ndx + self.indices[chain1][a]: types['chain1'][a] for a in types['chain1'].keys()},
                         'chain2': {c2_ndx + self.indices[chain2][a]: types['chain2'][a] for a in types['chain2'].keys()}}

        # update bonds - 2 new bonds between dummy atoms and carbon
        dummy_bonds = [[c2_ndx + self.indices[chain2]['C2'], c2_ndx + self.indices[chain2]['D2']],
                       [c2_ndx + self.indices[chain2]['C1'], c2_ndx + self.indices[chain2]['D1']]]

        # define indices of left-over radicals
        radicals = [c1_ndx + self.indices[chain1]['C2']]

        # define which improper dihedrals to remove -- written in same order as .itp file!!!
        # note that the order of the atoms may be different for each chain
        impropers = {'a': {1: ['H2', 'C1', 'H1', 'C2'], 2: ['C1', 'C3', 'C2', 'H3']},
                     'b': {1: ['C2', 'H2', 'C1', 'H1'], 2: ['C1', 'C3', 'C2', 'H3']}}

        chain1_impropers = [1]
        chain2_impropers = [1, 2]
        rm_improper = []
        for c in chain1_impropers:
            rm_improper.append([c1_ndx + self.indices[chain1][x] for x in impropers[chain1][c]])
        for c in chain2_impropers:
            rm_improper.append([c2_ndx + self.indices[chain2][x] for x in impropers[chain2][c]])

        # define terminated atoms
        terminated = [c1_ndx + self.indices[chain1]['C1'], c2_ndx + self.indices[chain2]['C2'], c2_ndx +
                      self.indices[chain2]['C1']]

        return reacted_types, dummy_bonds, radicals, rm_improper, terminated

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
        c1_ndx -= self.indices[chain1]['C1']
        c2_ndx -= self.indices[chain2]['C2']

        # types after reaction
        types = {'chain1': {'C1': 'c3', 'C2': 'c2', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'ha',
                            'H4': 'ha', 'H5': 'ha'},  # chain1 contains c1
                 'chain2': {'C1': 'c3', 'C2': 'c3', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'hc',
                            'H4': 'ha', 'H5': 'ha'}}  # chain2 contains c2 radical

        # update types
        reacted_types = {'chain1': {c1_ndx + self.indices[chain1][a]: types['chain1'][a] for a in types['chain1'].keys()},
                         'chain2': {c2_ndx + self.indices[chain2][a]: types['chain2'][a] for a in types['chain2'].keys()}}

        # update bonds - no new bonds between dummy atoms and carbon
        dummy_bonds = []

        # define indices of left-over radicals
        radicals = [c1_ndx + self.indices[chain1]['C2']]

        # define which improper dihedrals to remove -- written in same order as .itp file!!!
        # note that the order of the atoms may be different for each chain
        impropers = {'a': {1: ['H2', 'C1', 'H1', 'C2'], 2: ['C1', 'C3', 'C2', 'H3']},
                     'b': {1: ['C2', 'H2', 'C1', 'H1'], 2: ['C1', 'C3', 'C2', 'H3']}}

        chain1_impropers = [1]
        chain2_impropers = [2]
        rm_improper = []
        for c in chain1_impropers:
            rm_improper.append([c1_ndx + self.indices[chain1][x] for x in impropers[chain1][c]])
        for c in chain2_impropers:
            rm_improper.append([c2_ndx + self.indices[chain2][x] for x in impropers[chain2][c]])

        # define terminated atoms
        terminated = [c1_ndx + self.indices[chain1]['C1'], c2_ndx + self.indices[chain2]['C2']]

        return reacted_types, dummy_bonds, radicals, rm_improper, terminated

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
        c_name = self.chains[chain][c]

        # to get indexing right
        c_ndx -= self.indices[chain][c_name]

        # types after reaction. Keeping this dictionary format so it integrates easily with xlinking algorithm
        types = {'chain': {self.chains[chain][c]: 'c3', self.dummy_connectivity[chain][c]: 'hc'}}

        for i in self.hydrogen_connectivity[c]:  # turn already attached carbon(s) to c3
            types['chain'][i] = 'hc'

        # update types
        reacted_types = {'chain': {c_ndx + self.indices[chain][a]: types['chain'][a] for a in types['chain'].keys()}}

        # add dummy atom bond
        dummy_bonds = [list(reacted_types['chain'].keys())]

        # no radicals are produced (only eliminated) -- TODO: make sure the elimination is accounted for
        radicals = []
        # stopped here

        # define which improper dihedrals to remove -- written in same order as .itp file!!!
        # note that the order of the atoms may be different for each chain
        # NOTE: C3 not tested
        impropers = {'a': {'C1': ['H2', 'C1', 'H1', 'C2'], 'C2': ['C1', 'C3', 'C2', 'H3'],
                           'C3': ['C4', 'C2', 'C3', 'H4']},
                     'b': {'C1': ['C2', 'H2', 'C1', 'H1'], 'C2': ['C1', 'C3', 'C2', 'H3'],
                           'C3': ['C4', 'C2', 'C3', 'H4']}}

        rm_improper = [[c_ndx + self.indices[chain][x] for x in impropers[chain][c_name]]]

        # define terminated atoms
        terminated = [c_ndx + self.indices[chain][c_name]]

        return reacted_types, dummy_bonds, radicals, rm_improper, terminated

