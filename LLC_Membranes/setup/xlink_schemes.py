#!/usr/bin/env python

import sys

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

        if monomer in ['Dibrpyr14', 'MOL']:  # could build this into annotations. Might be better like this

            # names of atoms that make up relevant segements of each chain
            self.chains = {'a': {'C': 'C1', 'C1': 'C2', 'C2': 'C3', 'C3': 'C4', 'H': 'H1', 'H1': 'H2', 'H2': 'H3',
                                 'H3': 'H4', 'H4': 'H5'},
                           'b': {'C45': 'C1', 'C44': 'C2', 'C43': 'C3', 'C42': 'C4', 'H81': 'H1', 'H80': 'H2',
                                 'H79': 'H3', 'H78': 'H4', 'H77': 'H5'}
                           }

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

            self.dummy_mass = 1.008  # mass of hydrogen

    # def determine_chains(self, c1, c2):
    # TODO: delete after testing below
    #
    #     chain1 = None
    #     chain2 = None
    #     for k in self.chains.keys():
    #         if c1 in self.chains[k].values():
    #             chain1 = k
    #         if c2 in self.chains[k].values():
    #             chain2 = k
    #
    #     return chain1, chain2

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
                            'H4': 'ha', 'H5': 'ha', 'D1': 'hc'},
                 'chain2': {'C1': 'c3', 'C2': 'c3', 'C3': 'c2', 'C4': 'c2', 'H1': 'hc', 'H2': 'hc', 'H3': 'hc',
                            'H4': 'ha', 'H5': 'ha', 'D1': 'hc', 'D2': 'hc'}}

        # update types
        reacted_types = {'chain1': {c1_ndx + self.indices[chain1][a]: types['chain1'][a] for a in types['chain1'].keys()},
                         'chain2': {c2_ndx + self.indices[chain2][a]: types['chain2'][a] for a in types['chain2'].keys()}}

        # update bonds - 3 new bonds between dummy atoms and carbon
        dummy_bonds = [[c1_ndx + self.indices[chain1]['C1'], c1_ndx + self.indices[chain1]['D1']],
                       [c2_ndx + self.indices[chain2]['C2'], c2_ndx + self.indices[chain2]['D2']],
                       [c2_ndx + self.indices[chain2]['C1'], c2_ndx + self.indices[chain2]['D1']]]

        # define indices of left-over radicals
        radicals = [c1_ndx + self.indices[chain1]['C2']]

        # define which improper dihedrals to remove -- written in same order as .itp file!!!
        # note that the order of the atoms may be different for each chain
        impropers = {'a': {1: ['H2', 'C1', 'H1', 'C2'], 2: ['C1', 'C3', 'C2', 'H3']},
                     'b': {1: ['C2', 'H2', 'C1', 'H1'], 2: ['C1', 'C3', 'C2', 'H3']}}

        # step through this for chain2's second iteration (1st is right. all after are wrong)
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

        # step through this for chain2's second iteration (1st is right. all after are wrong)
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

        c = list(atoms.keys())[0]
        c_ndx = list(atoms.values())[0]

        chain = self.determine_chains(c)[0]

        # to get indexing right
        c_ndx -= self.indices[chain][self.chains[chain][c]]

        # types after reaction
        types = {'chain': {self.chains[chain][c]: 'c3', {self.dummy_connectivity[c]: 'hc'}}}

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

        # step through this for chain2's second iteration (1st is right. all after are wrong)
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

