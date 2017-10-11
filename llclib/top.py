#! /usr/bin/env python

import numpy as np


def water_indices(pos, max, min, buffer, t):
    """
    :param pos: xyz coordinates of all atoms in system
    :param max: maximum z coordinate before water will be removed
    :param min: minimum z coordinate for water to not be removed
    :param buffer: Percent of membrane thickness to adjust max and min by
    :param t: mdtraj trajectory object
    :return: indices of water molecules to be removed
    """

    indices = [a.index for a in t.topology.atoms if a.name == 'O' and 'HOH' in str(a.residue)]

    thickness = max - min
    b = thickness * buffer
    max -= b
    min += b

    remove = []

    for i in range(len(indices)):
        if min > pos[0, indices[i], 2] or pos[0, indices[i], 2] > max:
            remove.append(indices[i])

    for i in range(len(remove)):
        remove.append(remove[i] + 1)
        remove.append(remove[i] + 2)

    return remove


def keep(list, t, exclude=False):

    """
    :param list: list of atom indices to keep (or not keep if exclude is True)
    :param t: mdtraj trajectory object
    :param exclude: If this is set to true, then all atoms will be kept except for those in list
    :return: new trajectory object
    """

    if exclude:
        keep = [a.index for a in t.topology.atoms if a.index not in list]
    else:
        keep = [a.index for a in t.topology.atoms if a.index in list]


    return t.atom_slice(keep), keep