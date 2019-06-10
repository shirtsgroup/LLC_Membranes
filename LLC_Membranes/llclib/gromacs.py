#!/usr/bin/env python
"""
Run GROMACS commands with subprocess
"""

import subprocess
import os

script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def simulate(mdp, top, gro, out, verbose=False, em_energy=False, mpi=False, nprocesses=4):

    gmx = "gmx"
    if mpi:
        gmx = "gmx_mpi mpirun -np %d" % nprocesses

    grompp = '%s grompp -f %s -c %s -p %s -o %s' % (gmx, mdp, gro, top, out)

    if verbose:
        p1 = subprocess.Popen(grompp.split())
    else:
        p1 = subprocess.Popen(grompp.split(), stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)

    p1.wait()

    mdrun = 'gmx mdrun -deffnm %s' % out

    if verbose:
        mdrun += ' -v'
        p2 = subprocess.Popen(mdrun.split())
    else:
        p2 = subprocess.Popen(mdrun.split(), stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)

    p2.wait()

    if em_energy:

        nrg = subprocess.check_output(
            ["awk", "/Potential Energy/ {print $4}", "%s.log" % out])  # get Potential energy from em.log

        try:
            return float(nrg.decode("utf-8"))
        except ValueError:
            return 0  # If the system did not energy minimize, the above statement will not work because nrg will be an
            # empty string. Make nrg=0 so placement gets attempted again


def insert_molecules(gro, solute, n, out, scale=0.4, mpi=False, nprocesses=4):
    """ Insert n solutes into a .gro file

    :param gro: name of coordinate file where solutes will be placed
    :param solute: name of solute to add to gro
    :param n: number of solutes to add to gro
    :param out: name of output configuration
    :param scale: Scale factor to multiply Van der Waals radii from the database in share/gromacs/top/vdwradii.dat. The
    default value of 0.57 yields density close to 1000 g/l for proteins in water

    :type gro: str
    :type solute: sol
    :type n: int
    :type out: str
    :type scale: float
    """

    gmx = "gmx"
    if mpi:
        gmx = "gmx_mpi mpirun -np %d" % nprocesses

    insert = "%s insert-molecules -f %s -ci %s/../top/topologies/%s -nmol %d -o %s -scale %s" % (gmx, gro,
                                                                                                 script_location,
                                                                                                 solute, n, out, scale)

    p = subprocess.Popen(insert.split(), stdout=open('inserted.txt', 'w'), stderr=subprocess.STDOUT)
    p.wait()

    nadded = 0
    with open('inserted.txt', 'r') as f:
        for line in f:
            if line.count('Added') > 0:
                nadded = int(line.split()[1])

    return nadded
