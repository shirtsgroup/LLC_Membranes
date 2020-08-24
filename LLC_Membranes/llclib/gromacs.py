#!/usr/bin/env python
"""
Run GROMACS commands with subprocess
"""

import subprocess
import os

script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def simulate(mdp, top, gro, out, verbose=True, em_energy=False, mpi=False, nprocesses=4, dd=None, restraints=False):
    """ A wrapper for running GROMACS molecular dynamics simulations

    :param mdp: name of GROMACS Molecular Dynamics Paramters (mdp) file
    :param top: name of GROMACS topology (.top)
    :param gro: name of GROMACS coordinate file (.gro)
    :param out: name of output simulation files
    :param verbose: if True, prints simulation output to the screen
    :param em_energy: if this is an energy minimzation and this argument is True, return the final total energy
    :param mpi: if True, run simulation in parallel using MPI
    :param nprocesses: number of MPI process for an MPI simulation
    :param dd: domain decomposition grid for parallelization. If this is not specified, GROMACS decides (which usually
    works)
    :param restraints: True if position restraints are applied

    :type mdp: str
    :type top: str
    :type gro: str
    :type out: str
    :type verbose: bool
    :type em_energy: bool
    :type mpi: bool
    :type nprocesses: int
    :type dd: list
    :type restraints: bool
    """

    gmx = "gmx"
    if mpi:
        gmx = "mpirun -np %d gmx_mpi" % nprocesses

    grompp = '%s grompp -f %s -c %s -p %s -o %s' % (gmx, mdp, gro, top, out)

    if restraints:
        grompp += ' -r %s' % gro

    if verbose:
        p1 = subprocess.Popen(grompp.split())
    else:
        p1 = subprocess.Popen(grompp.split(), stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)

    p1.wait()

    mdrun = '%s mdrun -deffnm %s' % (gmx, out)

    if dd and mpi:
        mdrun += ' -dd %s %s %s' % tuple(dd)

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
            return 1  # If the system did not energy minimize, the above statement will not work because nrg will be an
            # empty string. Make nrg=1 so placement gets attempted again


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
        gmx = "mpirun -np %d gmx_mpi" % nprocesses

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


def editconf(gro, out, d=None, center=True, box_type='cubic'):
    """ Run gmx editconf. See their documentation: http://manual.gromacs.org/documentation/2019/onlinehelp/gmx-editconf.html

    :param gro: name of input coordinate file to put a box around
    :param out: name of output file
    :param d: if not None, distance between solute and box
    :param center: center solute in box
    :param box_type: type of box (only cubic is implemented)

    :type gro: str
    :type out: str
    :type d: float
    :type center: bool
    :type box_type: str
    """

    editconf = 'gmx editconf -f %s -o %s -bt %s' % (gro, out, box_type)
    if center:
        editconf += ' -c'
    if d is not None:
        editconf += ' -d %f' % d

    p = subprocess.Popen(editconf.split())
    p.wait()
